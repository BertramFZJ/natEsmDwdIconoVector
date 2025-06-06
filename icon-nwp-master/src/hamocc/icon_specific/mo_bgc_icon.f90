! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! @brief Main biogeochemical subroutine, called at each time step
!
! This subroutine calls all routines that calculate changes of pelagic biogeochemical
! tracers due to local processes (like photosythesis, heterotrophic
! processes, N-fixation, and denitrification), the air-sea gas
! exchange of carbon dioxide, oxygen, dinitrogen, and the
! benthic processes. It further calls the computation of the vertical displacement of
! particles
!
! called by mo_hydro_ocean_run:perform_ho_stepping

#include "icon_definitions.inc"
#include "omp_definitions.inc"

MODULE mo_bgc_icon
#ifdef _OPENMP
  USE omp_lib
#endif

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish,message_to_own_unit

  USE mo_model_domain,        ONLY: t_patch,t_patch_3D

  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range

  USE mo_bgc_icon_comm,       ONLY: update_icon, update_bgc, hamocc_state, &
       &                            set_bgc_tendencies_output, set_bgc_tendencies_output_sedon
  USE mo_dynamics_config,     ONLY: nold
  USE mo_hamocc_nml,          ONLY: i_settling, l_cyadyn,l_bgc_check,io_stdo_bgc,l_implsed, &
       &                            l_dynamic_pi, l_pdm_settling
  USE mo_ocean_nml,           ONLY: lsediment_only
  USE mo_control_bgc,         ONLY: ndtdaybgc,  &
       &                        ldtrunbgc, bgc_nproma, bgc_zlevs

  USE mo_cyano,               ONLY: cyano, cyadyn
  USE mo_bgc_surface,         ONLY: gasex, update_weathering, dust_deposition, &
&                                   nitrogen_deposition, update_linage
  USE mo_bgc_bcond,           ONLY: ext_data_bgc
  USE mo_hamocc_diagnostics,  ONLY: get_inventories, get_omz
  USE mo_carchm,              ONLY: calc_dissol
  USE mo_powach,              ONLY: powach, powach_impl
  USE mo_sedmnt, ONLY         : ini_bottom
  USE mo_timer, ONLY          : timer_bgc_up_bgc, timer_bgc_swr, timer_bgc_wea,timer_bgc_depo, &
    &                           timer_bgc_chemcon, timer_bgc_ocprod, timer_bgc_sett,timer_bgc_cya,&
    &                           timer_bgc_gx, timer_bgc_calc, timer_bgc_powach, timer_bgc_up_ic, &
    &                           timer_bgc_tend, timer_start, timer_stop, timers_level,timer_bgc_agg
  USE mo_settling,   ONLY    : settling, settling_pdm
  USE mo_aggregates, ONLY     : mean_aggregate_sinking_speed
  USE mo_ocean_hamocc_couple_state, ONLY: t_ocean_to_hamocc_state, t_hamocc_to_ocean_state, &
    & t_hamocc_ocean_state
!   USE mo_util_dbg_prnt,          ONLY: dbg_print
  USE mo_chemcon, ONLY         : chemcon
  USE mo_ocprod, ONLY: ocprod
  USE mo_sedshi, ONLY: sedshi
  USE mo_hamocc_swr_absorption, ONLY: swr_absorption

  USE mo_bgc_memory_types, ONLY: t_bgc_memory, t_sediment_memory, t_aggregates_memory, &
    & bgc_local_memory, sediment_local_memory, aggregates_memory, bgc_memory_copies
  USE mo_fortran_tools, ONLY: set_acc_host_or_device

  USE mpi

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: bgc_icon

  INTEGER :: loadBalanceWeight = 0
  INTEGER :: loadBalanceNum = 0
  LOGICAL :: printLoadBalance = .TRUE.
  INTEGER :: numSteps = 0
  DOUBLE PRECISION :: exeTime = 0.0D0

CONTAINS

SUBROUTINE BGC_ICON(p_patch_3D, hamocc_ocean_state, ssh, pddpo, ptiestu, lacc)

  IMPLICIT NONE

  TYPE(t_patch_3d ),TARGET, INTENT(in)   :: p_patch_3d
  TYPE(t_hamocc_ocean_state), TARGET     :: hamocc_ocean_state
  REAL(wp), INTENT(IN) :: pddpo(bgc_nproma, bgc_zlevs, p_patch_3d%p_patch_2d(1)%nblks_c)
  REAL(wp), INTENT(IN) :: ptiestu(bgc_nproma, bgc_zlevs, p_patch_3d%p_patch_2d(1)%nblks_c)
  REAL(wp), INTENT(IN) :: ssh(bgc_nproma, p_patch_3d%p_patch_2d(1)%nblks_c)
  LOGICAL, INTENT(IN), OPTIONAL          :: lacc

  ! Local variables
  INTEGER ::  jb
  INTEGER :: start_index, end_index
  INTEGER :: itrig_chemcon
  !INTEGER, POINTER :: levels(:)
  INTEGER :: levels(bgc_nproma)

  TYPE(t_bgc_memory), POINTER :: local_bgc_memory
  TYPE(t_sediment_memory), POINTER :: local_sediment_memory
  TYPE(t_aggregates_memory), POINTER :: local_aggregate_memory
  INTEGER :: local_memory_idx, test_memory_copies

  TYPE(t_subset_range), POINTER :: all_cells
  INTEGER :: alloc_cell_blocks
  TYPE(t_patch),POINTER    :: p_patch

  TYPE(t_ocean_to_hamocc_state), POINTER           :: ocean_to_hamocc_state
  TYPE(t_hamocc_to_ocean_state), POINTER           :: hamocc_to_ocean_state

  INTEGER :: i,jc,k,kpke
  LOGICAL :: lzacc

  DOUBLE PRECISION :: startTimer, finishTimer

  CHARACTER(LEN=*), PARAMETER  :: str_module = 'BGC_ICON'  ! Output of module for 1 line debug
!   INTEGER :: idt_src
!   TYPE(t_subset_range), POINTER :: owned_cells

  CALL set_acc_host_or_device(lzacc, lacc)

  ocean_to_hamocc_state => hamocc_ocean_state%ocean_to_hamocc_state
  hamocc_to_ocean_state => hamocc_ocean_state%hamocc_to_ocean_state
  !-----------------------------------------------------------------------
  p_patch   => p_patch_3D%p_patch_2d(1)
  alloc_cell_blocks =  p_patch%alloc_cell_blocks
  !--------------------------------------------
  !
  !----------------------------------------------------------------------
  !
  all_cells => p_patch%cells%ALL
!   owned_cells => p_patch%cells%owned



  !----------------------------------------------------------------------




  IF (lsediment_only) THEN
    ! trigger chemcon at depth only once per run cycle
    itrig_chemcon=merge(1,0,ldtrunbgc<1)
  ELSE
    ! trigger chemcon at depth only once per day
    itrig_chemcon=mod(ldtrunbgc,ndtdaybgc)+1
  ENDIF

  !
  !
  !----------------------------------------------------------------------

  !---------DEBUG DIAGNOSTICS-------------------------------------------
!   idt_src=1  ! output print level (1-5, fix)
!   CALL dbg_print('h_old'           ,ocean_to_hamocc_state%h_old, str_module,idt_src, owned_cells)
!   CALL dbg_print('co2_mixing'      ,ocean_to_hamocc_state%co2_mixing_ratio, str_module,idt_src, owned_cells)
!   CALL dbg_print('short_wave_flux' ,ocean_to_hamocc_state%short_wave_flux, str_module,idt_src, owned_cells)
!   CALL dbg_print('concSum'         ,ocean_to_hamocc_state%ice_concentration_sum, str_module,idt_src, owned_cells)
!   CALL dbg_print('salinity'        ,ocean_to_hamocc_state%salinity, str_module,idt_src, owned_cells)
!   CALL dbg_print('temperature'     ,ocean_to_hamocc_state%temperature, str_module,idt_src, owned_cells)
!   CALL dbg_print('wind10m'         ,ocean_to_hamocc_state%wind10m, str_module,idt_src, owned_cells)
 !----------------------------------------------------------------------

IF(l_bgc_check)THEN
 CALL message_to_own_unit('1. before bgc','inventories',io_stdo_bgc)
 CALL get_inventories(hamocc_state, ssh, pddpo, hamocc_state%p_prog(nold(1))%tracer, &
   &                  p_patch_3d, 0._wp, 0._wp, lacc=lzacc)
ENDIF

kpke = MAXVAL(p_patch_3D%p_patch_1d(1)%dolic_c(:,:))

IF (.not. lsediment_only) THEN

local_memory_idx = 0
test_memory_copies = 1
local_bgc_memory => bgc_local_memory(local_memory_idx)
local_sediment_memory => sediment_local_memory(local_memory_idx)
local_aggregate_memory => aggregates_memory(local_memory_idx)

!$ACC DATA CREATE(levels) IF(lzacc)

!DIR$ INLINE
#ifdef _OPENMP
!ICON_OMP_PARALLEL PRIVATE(local_memory_idx, local_bgc_memory, local_sediment_memory, local_aggregate_memory)
local_memory_idx = omp_get_thread_num()
! write(0,*) "local_memory_idx=", local_memory_idx
local_bgc_memory => bgc_local_memory(local_memory_idx)
local_sediment_memory => sediment_local_memory(local_memory_idx)
local_aggregate_memory => aggregates_memory(local_memory_idx)

!ICON_OMP_SINGLE
test_memory_copies = OMP_GET_NUM_THREADS()
IF (test_memory_copies /= bgc_memory_copies) &
  & CALL finish(str_module, "test_memory_copies /= bgc_memory_copies")
!ICON_OMP_END_SINGLE

!ICON_OMP_DO PRIVATE(levels, start_index, end_index)
#endif

  DO jb = all_cells%start_block, all_cells%end_block

        CALL get_index_range(all_cells, jb, start_index, end_index)
        !  tracer 1: potential temperature
        !  tracer 2: salinity
        !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        levels(start_index:end_index) = p_patch_3D%p_patch_1d(1)%dolic_c(start_index:end_index,jb)
        !$ACC END KERNELS

        start_detail_timer(timer_bgc_up_bgc,5)
        CALL update_bgc(local_bgc_memory, local_sediment_memory, start_index,end_index,levels,&
             & pddpo(:,:,jb),&  ! cell thickness
             &jb, hamocc_state%p_prog(nold(1))%tracer, &
             & ocean_to_hamocc_state%co2_mixing_ratio(:,jb)                        & ! co2mixing ratio
             & ,hamocc_state%p_diag,hamocc_state%p_sed, hamocc_state%p_tend, kpke, lacc=lzacc)
        stop_detail_timer(timer_bgc_up_bgc,5)

        CALL ini_bottom(local_bgc_memory, start_index, end_index, levels, pddpo(:,:,jb), lacc=lzacc)

        start_detail_timer(timer_bgc_swr,5)
       ! Net solar radiation update and swr_frac
        CALL swr_absorption(local_bgc_memory, start_index,end_index,levels,                   &
 &                          ocean_to_hamocc_state%short_wave_flux(:,jb),                                 & ! SW radiation
 &                          ocean_to_hamocc_state%ice_concentration_sum(:,jb),                           & ! sea ice concentration
 &                          pddpo(:,:,jb), lacc=lzacc)                                                  ! level thickness

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=1,kpke
          DO jc = start_index,end_index
            hamocc_to_ocean_state%swr_fraction(jc,k,jb) = local_bgc_memory%swr_frac(jc,k)
          END DO
        END DO
        !$ACC END PARALLEL

        stop_detail_timer(timer_bgc_swr,5)

       ! Linear age
        CALL update_linage(local_bgc_memory, levels, start_index, end_index,  & ! index range, levels, salinity
   &                 pddpo(:,:,jb), lacc=lzacc)! cell thickness (check for z0)

       ! Biogeochemistry

        start_detail_timer(timer_bgc_wea,5)
       ! Weathering fluxes
        CALL update_weathering(local_bgc_memory, start_index, end_index,  & ! index range, levels, salinity
   &                 pddpo(:,:,jb),&! cell thickness (check for z0)
   &                 ssh(:,jb), lacc=lzacc) ! surface_height

        stop_detail_timer(timer_bgc_wea,5)

        start_detail_timer(timer_bgc_depo,5)
      ! Dust deposition
        CALL dust_deposition(local_bgc_memory, start_index, end_index,  & ! index range, levels,
   &                 pddpo(:,:,jb), &! cell thickness (check for z0)
   &                 ssh(:,jb),& ! surface_height
   &                 ext_data_bgc%dusty(:,jb), lacc=lzacc)       ! dust input
      ! Nitrogen deposition
        CALL nitrogen_deposition(local_bgc_memory, start_index, end_index,  & ! index range, levels
   &                 pddpo(:,:,jb), &! cell thickness (check for z0)
   &                 ssh(:,jb),& ! surface_height
   &                 ext_data_bgc%nitro(:,jb), lacc=lzacc)      ! nitrogen input

        stop_detail_timer(timer_bgc_depo,5)
       !----------------------------------------------------------------------
       ! Calculate chemical properties

        start_detail_timer(timer_bgc_chemcon,5)
        CALL chemcon(local_bgc_memory, start_index, end_index,levels,  ocean_to_hamocc_state%salinity(:,:,jb), & ! index range, levels, salinity
   &                 ocean_to_hamocc_state%temperature(:,:,jb),                              & ! pot. temperature
   &                 pddpo(:,:,jb),                    & ! cell thickness
   &                 ptiestu(:,:,jb) ,  &           ! depths at interface
   &                 itrig_chemcon, lacc=lzacc)
        stop_detail_timer(timer_bgc_chemcon,5)
       !----------------------------------------------------------------------
       ! Calculate plankton dynamics and particle settling

       IF(i_settling==2)then
         ! sinking speeds from MAGO aggregation scheme (MARMA)

         start_detail_timer(timer_bgc_agg,5)
         CALL mean_aggregate_sinking_speed (local_bgc_memory, local_aggregate_memory, levels, start_index, end_index, &
   &                 pddpo(:,:,jb), & ! cell thickness
   &                 ocean_to_hamocc_state%press_hyd(:,:,jb),                 & ! hydrostatic pressure
   &                 ocean_to_hamocc_state%temperature(:,:,jb),               & ! pot. temperature
   &                 ocean_to_hamocc_state%salinity(:,:,jb))                    ! salinity

        stop_detail_timer(timer_bgc_agg,5)
       ENDIF

        start_detail_timer(timer_bgc_ocprod,5)
         ! plankton dynamics and remineralization
         CALL ocprod(local_bgc_memory, levels, start_index,end_index, ocean_to_hamocc_state%temperature(:,:,jb),&
   &               pddpo(:,:,jb), & ! cell thickness
   &               ssh(:,jb),& ! surface height
   &               ptiestu(:,:,jb),& ! depths at interface
   &               l_dynamic_pi, kpke, lacc=lzacc) ! depths at interface
        stop_detail_timer(timer_bgc_ocprod,5)

        start_detail_timer(timer_bgc_sett,5)
         ! particle settling
        IF (l_pdm_settling)then
         CALL settling_pdm(local_bgc_memory, local_sediment_memory, levels,start_index, end_index, &
   &                   pddpo(:,:,jb)) ! cell thickness

        ELSE
         CALL settling(local_bgc_memory, local_sediment_memory, levels,start_index, end_index, &
   &                   pddpo(:,:,jb),& ! cell thickness
   &                   ssh(:,jb), lacc=lzacc)                   ! surface height

        ENDIF
        stop_detail_timer(timer_bgc_sett,5)





       !----------------------------------------------------------------------
       ! Calculate N2 fixation

       start_detail_timer(timer_bgc_cya,5)
       IF (l_cyadyn) THEN
        ! dynamic cyanobacteria
        CALL cyadyn(local_bgc_memory, levels, start_index,end_index, &  ! vertical range, cell range,
     &               pddpo(:,:,jb),& ! cell thickness
     &               ssh(:,jb), &                 ! surface height
     &               ocean_to_hamocc_state%temperature(:,:,jb), &        ! pot. temperature
     &               ptiestu(:,:,jb),& ! depths at interface
     &               l_dynamic_pi, kpke, lacc=lzacc) ! depths at interface
       ELSE
        ! diagnostic N2 fixation
        CALL cyano (local_bgc_memory, start_index, end_index,pddpo(:,:,jb),&
     &               ssh(:,jb))                 ! surface height
       endif
       stop_detail_timer(timer_bgc_cya,5)


       !----------------------------------------------------------------------
       ! Calculate gas exchange

        start_detail_timer(timer_bgc_gx,5)
        CALL gasex(local_bgc_memory, start_index, end_index,    &
  &               pddpo(:,:,jb),&  ! cell thickness
  &               ssh(:,jb), &                   ! surface height
  &               ocean_to_hamocc_state%temperature(:,:,jb), &          ! pot. temperature
  &               ocean_to_hamocc_state%salinity(:,:,jb), &          ! salinity
  &               ocean_to_hamocc_state%wind10m(:,jb)            , &          ! 10m wind speed
  &               ocean_to_hamocc_state%ice_concentration_sum(:,jb), lacc=lzacc)                              ! sea ice concentration

       stop_detail_timer(timer_bgc_gx,5)
        !----------------------------------------------------------------------
        ! Calculate carbonate dissolution

        IF(printLoadBalance) THEN
            loadBalanceWeight = loadBalanceWeight + SUM(levels(start_index : end_index))
            loadBalanceNum = loadBalanceNum + end_index - start_index + 1
        END IF

        start_detail_timer(timer_bgc_calc,5)
        startTimer = MPI_Wtime()
        CALL calc_dissol(local_bgc_memory, start_index, end_index, levels,   &
   &               pddpo(:,:,jb),& ! cell thickness
   &               ocean_to_hamocc_state%salinity(:,:,jb),         &  ! salinity
   &               ptiestu(:,:,jb), lacc=lzacc) !depths at interface
        finishTimer = MPI_Wtime()
        exeTime = exeTime + finishTimer - startTimer
        stop_detail_timer(timer_bgc_calc,5)
       !----------------------------------------------------------------------
        ! Calculate sediment dynamics
        start_detail_timer(timer_bgc_powach,5)
        if(l_implsed)then
         CALL powach_impl(local_bgc_memory, local_sediment_memory,  start_index, end_index,    &
   &               ocean_to_hamocc_state%salinity(:,:,jb))          ! salinity
         else

        CALL powach(local_bgc_memory, local_sediment_memory, start_index, end_index, &
   &               ocean_to_hamocc_state%salinity(:,:,jb),          &! salinity
   &               pddpo(:,:,jb), lacc=lzacc)  ! cell thickness
         endif
        stop_detail_timer(timer_bgc_powach,5)

        if(mod(ldtrunbgc,ndtdaybgc).eq.0) CALL sedshi(local_bgc_memory, local_sediment_memory, &
                                                      start_index, end_index, lacc=lzacc)

        start_detail_timer(timer_bgc_up_ic,5)
        CALL update_icon(local_bgc_memory, start_index,end_index,levels,&
  &               pddpo(:,:,jb),&  ! cell thickness
  &               jb, hamocc_state%p_prog(nold(1))%tracer,            &
  &               hamocc_to_ocean_state%co2_flux(:,jb), lacc=lzacc)          ! co2flux for coupling
        stop_detail_timer(timer_bgc_up_ic,5)

        start_detail_timer(timer_bgc_tend,5)
        CALL set_bgc_tendencies_output(local_bgc_memory, local_sediment_memory, local_aggregate_memory, &
          & start_index,end_index,levels, &
  &          p_patch_3D%p_patch_1d(1)%prism_thick_c(:,:,jb),&  ! cell thickness
  &                                   jb, &
  &                                   hamocc_state%p_tend,            &
  &                                   hamocc_state%p_diag,            &
  &                                   hamocc_state%p_sed,             &
  &                                   hamocc_state%p_agg, kpke, lacc=lzacc)

        stop_detail_timer(timer_bgc_tend,5)
 ENDDO
#ifdef _OPENMP
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL
#endif

IF(printLoadBalance) THEN
    WRITE(*,*) "DISSOL_LB ", loadBalanceWeight, " ", loadBalanceNum
    printLoadBalance = .FALSE.
END IF

numSteps = numSteps + 1
IF(numSteps == 200) THEN
    WRITE(*,*) "DISSOL_REPORT ", exeTime, " ", loadBalanceNum, " ", loadBalanceWeight
END IF

! O2 min depth & value diagnostics
CALL get_omz(hamocc_state, p_patch_3d, pddpo, ssh, lacc=lzacc)

!$ACC WAIT
!$ACC END DATA

ELSE
! offline sediment
!DIR$ INLINE
  DO jb = all_cells%start_block, all_cells%end_block

        local_memory_idx = 0
        local_bgc_memory => bgc_local_memory(local_memory_idx)
        local_sediment_memory => sediment_local_memory(local_memory_idx)


        CALL get_index_range(all_cells, jb, start_index, end_index)
        !  tracer 1: potential temperature
        !  tracer 2: salinity
        levels(start_index:end_index) = p_patch_3D%p_patch_1d(1)%dolic_c(start_index:end_index,jb)

        start_detail_timer(timer_bgc_up_bgc,5)

        CALL update_bgc(local_bgc_memory, local_sediment_memory, start_index,end_index,levels,&
             & pddpo(:,:,jb),&  ! cell thickness
             &jb, hamocc_state%p_prog(nold(1))%tracer, &
             & ocean_to_hamocc_state%co2_mixing_ratio(:,jb)                        & ! co2mixing ratio
             & ,hamocc_state%p_diag,hamocc_state%p_sed, hamocc_state%p_tend, kpke)

        stop_detail_timer(timer_bgc_up_bgc,5)

        CALL ini_bottom(local_bgc_memory, start_index,end_index,levels,pddpo(:,:,jb))

       !----------------------------------------------------------------------
       ! Calculate chemical properties

        start_detail_timer(timer_bgc_chemcon,5)
        CALL chemcon(local_bgc_memory, start_index, end_index,levels,  ocean_to_hamocc_state%salinity(:,:,jb), & ! index range, levels, salinity
   &                 ocean_to_hamocc_state%temperature(:,:,jb),                              & ! pot. temperature
   &                 pddpo(:,:,jb),                    & ! cell thickness
   &                 ptiestu(:,:,jb) ,  &           ! depths at interface
   &                 itrig_chemcon)
        stop_detail_timer(timer_bgc_chemcon,5)

       !----------------------------------------------------------------------
        ! Calculate sediment dynamics
        start_detail_timer(timer_bgc_powach,5)
        if(l_implsed)then
         CALL powach_impl(local_bgc_memory, local_sediment_memory, start_index, end_index,    &
   &               ocean_to_hamocc_state%salinity(:,:,jb))          ! salinity
         else

        CALL powach(local_bgc_memory, local_sediment_memory, start_index, end_index, &
   &               ocean_to_hamocc_state%salinity(:,:,jb),          &! salinity
   &               pddpo(:,:,jb))  ! cell thickness
         endif
        stop_detail_timer(timer_bgc_powach,5)

        if(mod(ldtrunbgc,ndtdaybgc).eq.0) CALL sedshi(local_bgc_memory, local_sediment_memory, start_index,end_index)


        start_detail_timer(timer_bgc_tend,5)
        CALL set_bgc_tendencies_output_sedon(local_bgc_memory, local_sediment_memory, start_index,end_index, &
  &               p_patch_3D%p_patch_1d(1)%prism_thick_c(:,:,jb),&  ! cell thickness
  &                                   jb, &
  &                                   hamocc_state%p_tend,            &
  &                                   hamocc_state%p_sed)

        stop_detail_timer(timer_bgc_tend,5)

 ENDDO

ENDIF  ! lsediment_only

  ! Increment bgc time step counter of run (initialized in INI_BGC).
  !
  ldtrunbgc = ldtrunbgc + 1

  IF(l_bgc_check)THEN
   CALL message_to_own_unit('2. after bgc','inventories',io_stdo_bgc)
   CALL get_inventories(hamocc_state, ssh, pddpo, hamocc_state%p_prog(nold(1))%tracer, &
 &                      p_patch_3d, 1._wp, 1._wp, lacc=lzacc)
  ENDIF


END SUBROUTINE


END MODULE
