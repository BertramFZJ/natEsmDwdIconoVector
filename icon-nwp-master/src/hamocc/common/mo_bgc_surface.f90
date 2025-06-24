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

! @brief module contains gas exchange, weathering fluxes,
!        dust & nitrogen deposition

MODULE mo_bgc_surface

  USE mo_kind, ONLY           : wp
  USE mo_control_bgc, ONLY    : dtbgc, bgc_nproma, bgc_zlevs
  USE mo_bgc_memory_types, ONLY  : t_bgc_memory
  USE mo_fortran_tools, ONLY  : set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gasex, update_weathering, dust_deposition, nitrogen_deposition,&
&           update_linage


contains

SUBROUTINE update_linage (local_bgc_mem, klev,start_idx,end_idx, pddpo, lacc)

! update linear age tracer
  USE mo_param1_bgc, ONLY     : iagesc
  USE mo_control_bgc, ONLY    : dtbgc

  ! Arguments
  TYPE(t_bgc_memory), POINTER    :: local_bgc_mem

  INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)
  INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir)
  INTEGER, INTENT(in), TARGET    :: klev(bgc_nproma)       !<  vertical levels

  REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]
  LOGICAL, INTENT(IN), OPTIONAL  :: lacc

  INTEGER :: jc,k, kpke
  REAL(wp) :: fac001
  LOGICAL :: lzacc

  CALL set_acc_host_or_device(lzacc, lacc)

  fac001 = dtbgc/(86400._wp*365._wp)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO jc = start_idx, end_idx
     kpke=klev(jc)
     !$ACC LOOP SEQ
     DO k = 2, kpke
      if(pddpo(jc,k) > EPSILON(0.5_wp)) then
         local_bgc_mem%bgctra(jc,k,iagesc) = local_bgc_mem%bgctra(jc,k,iagesc) + fac001
      endif
     ENDDO
     if(pddpo(jc,1) > EPSILON(0.5_wp)) local_bgc_mem%bgctra(jc,1,iagesc) = 0._wp
  ENDDO
  !$ACC END PARALLEL


END SUBROUTINE update_linage

SUBROUTINE update_weathering (local_bgc_mem, start_idx,end_idx, pddpo, za, lacc)
! apply weathering rates

  USE mo_memory_bgc, ONLY : calcinp, orginp, silinp
  USE mo_param1_bgc, ONLY     : isco212, ialkali, idoc, isilica,  &
 &                              korginp, ksilinp, kcalinp

  ! Arguments
  TYPE(t_bgc_memory), POINTER    :: local_bgc_mem

  INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)
  INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir)

  REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]
  REAL(wp), INTENT(in), TARGET   :: za(bgc_nproma)      !< surface height
  LOGICAL, INTENT(IN), OPTIONAL  :: lacc

  ! Local variables

  INTEGER :: jc
  LOGICAL :: lzacc

  CALL set_acc_host_or_device(lzacc, lacc)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO jc = start_idx, end_idx

  if(pddpo(jc,1) > EPSILON(0.5_wp)) then

    local_bgc_mem%bgctra(jc,1,idoc) = local_bgc_mem%bgctra(jc,1,idoc) + orginp / (pddpo(jc,1) + za(jc))
    local_bgc_mem%bgctra(jc,1,isco212) = local_bgc_mem%bgctra(jc,1,isco212) + calcinp / (pddpo(jc,1) + za(jc))
    local_bgc_mem%bgctra(jc,1,ialkali) = local_bgc_mem%bgctra(jc,1,ialkali) + 2._wp * calcinp / (pddpo(jc,1) + za(jc))
    local_bgc_mem%bgctra(jc,1,isilica) = local_bgc_mem%bgctra(jc,1,isilica) + silinp / (pddpo(jc,1) + za(jc))

    local_bgc_mem%bgcflux(jc,korginp) = orginp / (pddpo(jc,1) + za(jc))
    local_bgc_mem%bgcflux(jc,ksilinp) = silinp / (pddpo(jc,1) + za(jc))
    local_bgc_mem%bgcflux(jc,kcalinp) = calcinp / (pddpo(jc,1) + za(jc))

  endif

  ENDDO
  !$ACC END PARALLEL

END SUBROUTINE

SUBROUTINE nitrogen_deposition (local_bgc_mem, start_idx,end_idx, pddpo, za, nitinput, lacc)
! apply nitrogen deposition
  USE mo_param1_bgc, ONLY     : iano3, ialkali, kn2b,knitinp
  USE mo_bgc_constants, ONLY  : rmnit

  !Arguments
  TYPE(t_bgc_memory), POINTER    :: local_bgc_mem

  INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)
  INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir)

  REAL(wp),INTENT(in) :: nitinput(bgc_nproma )                         !< nitrogen input
  REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]
  REAL(wp), INTENT(in), TARGET   :: za(bgc_nproma)                   !< surface height
  LOGICAL, INTENT(IN), OPTIONAL :: lacc

  ! Local variables

  INTEGER :: jc
  REAL(wp) :: ninp
  LOGICAL :: lzacc

  CALL set_acc_host_or_device(lzacc, lacc)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO jc = start_idx, end_idx

  if(pddpo(jc,1) > EPSILON(0.5_wp)) then

      ! ndepo : CCMI wet+dry dep of NHx and NOy in kg (N) m-2 s-1
       ninp = nitinput(jc) / rmnit* dtbgc/(pddpo(jc,1)+za(jc)) ! kmol N m-3 time_step-1

       local_bgc_mem%bgctra(jc,1,iano3) = local_bgc_mem%bgctra(jc,1,iano3) + ninp
       local_bgc_mem%bgctra(jc,1,ialkali) = local_bgc_mem%bgctra(jc,1,ialkali) - ninp
       local_bgc_mem%bgctend(jc,1,kn2b)   = local_bgc_mem%bgctend(jc,1,kn2b) - ninp * (pddpo(jc,1) + za(jc))
       local_bgc_mem%bgcflux(jc,knitinp) = ninp

  endif

  ENDDO
  !$ACC END PARALLEL


END SUBROUTINE
SUBROUTINE dust_deposition (local_bgc_mem, start_idx,end_idx, pddpo, za, dustinp, lacc)
! apply dust deposition
  USE mo_memory_bgc, ONLY      : perc_diron
  USE mo_param1_bgc, ONLY     : iiron, idust
  USE mo_control_bgc, ONLY    : dtb

  !Arguments
  TYPE(t_bgc_memory), POINTER    :: local_bgc_mem


  INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)
  INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir)

  REAL(wp),INTENT(in) :: dustinp(bgc_nproma )                        !< dust input
  REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]
  REAL(wp), INTENT(in), TARGET   :: za(bgc_nproma)                   !< surface height
  LOGICAL, INTENT(IN), OPTIONAL  :: lacc

  ! Local variables

  INTEGER :: jc
  LOGICAL :: lzacc

  CALL set_acc_host_or_device(lzacc, lacc)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO jc = start_idx, end_idx

  if(pddpo(jc,1) > EPSILON(0.5_wp)) then

   local_bgc_mem%bgctra(jc,1,iiron) = local_bgc_mem%bgctra(jc,1,iiron) + dustinp(jc)*dtb/365._wp/(pddpo(jc,1)+za(jc)) *perc_diron

   local_bgc_mem%bgctra(jc,1,idust) = local_bgc_mem%bgctra(jc,1,idust) + dustinp(jc)*dtb/365._wp/(pddpo(jc,1)+za(jc))

  endif

 ENDDO
 !$ACC END PARALLEL

END SUBROUTINE


SUBROUTINE gasex (local_bgc_mem, start_idx,end_idx, pddpo, za, ptho, psao,  &
     &              pfu10, psicomo, lacc)
!! @brief Computes sea-air gass exchange
!!         for oxygen, O2, N2, N2O, DMS, and CO2.
!!


  USE mo_param1_bgc, ONLY     : igasnit, ian2o,  iatmco2, iphosph,         &
       &                        ioxygen, isco212, isilica,       &
       &                        ialkali, kcflux, koflux, knflux,          &
       &                        kn2oflux, idms, kdmsflux,kpco2, &
       &                        iammo, knh3flux, kcflux_cpl

  USE mo_memory_bgc, ONLY         : kg_denom

  USE mo_hamocc_nml, ONLY     : l_cpl_co2, atm_co2, atm_o2, atm_n2, &
       &                        l_N_cycle

  USE mo_bgc_constants, ONLY : cmh2ms

  USE mo_carchm,    ONLY: update_hi, update_hi_sub

  IMPLICIT NONE

  !! Arguments
  TYPE(t_bgc_memory), POINTER    :: local_bgc_mem

  INTEGER, INTENT(in)  :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)
  INTEGER, INTENT(in)  :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir)

  REAL(wp),INTENT(in) :: pddpo(bgc_nproma,bgc_zlevs) !< size of scalar grid cell (3rd REAL) [m]
  REAL(wp),INTENT(in) :: psao(bgc_nproma,bgc_zlevs)  !< salinity
  REAL(wp),INTENT(in) :: ptho(bgc_nproma,bgc_zlevs)  !< potential temperature
  REAL(wp),INTENT(in) :: pfu10(bgc_nproma)           !< forcing field wind speed
  REAL(wp),INTENT(in) :: psicomo(bgc_nproma)         !< sea ice concentration
  REAL(wp),INTENT(in) :: za(bgc_nproma)              !< sea surface height
  LOGICAL, INTENT(IN), OPTIONAL :: lacc

  !! Local variables

  INTEGER :: j

  REAL(wp) :: fluxd, fluxu
  REAL(wp) :: kwco2(start_idx:end_idx), kwo2, kwdms, kwn2o
  REAL(wp) :: scco2, sco2, scdms, scn2o
  REAL(wp) :: oxflux, niflux, nlaughflux, dmsflux
  REAL(wp) :: ato2, atn2, atco2(start_idx:end_idx), pco2
  REAL(wp) :: thickness(start_idx:end_idx)
  LOGICAL :: lzacc

  ! for extended N-cycle
  REAL (wp):: kgammo, kh_nh3i, kh_nh3, pka_nh3, ka_nh3, nh3sw, ammoflux
  REAL (wp):: ecoef, tabs

  !! Vectorization local variables
  LOGICAL :: vmask(start_idx:end_idx)

  CALL set_acc_host_or_device(lzacc, lacc)

  !
  !---------------------------------------------------------------------
  !

  !NEC$ nomove
  DO j = start_idx, end_idx
    IF( pddpo(j, 1) .GT. EPSILON(0.5_wp) ) THEN
        vmask(j) = .TRUE.
    ELSE
        vmask(j) = .FALSE.
    END IF
  END DO

  !NEC$ nomove
  DO j = start_idx, end_idx

    IF( vmask(j) ) THEN

        !*********************************************************************
        !
        !  Compute the Schmidt number of CO2 in seawater and the transfer
        !  (piston) velocity using the formulation presented in
        !   Wanninkhof 2014
        !*********************************************************************

        scco2 = 2116.8_wp - 136.25_wp*ptho(j,1) + 4.7353_wp*ptho(j,1)**2 &
        &   - 0.092307_wp*ptho(j,1)**3 + 0.0007555_wp*ptho(j,1)**4

        sco2 =  1920.4_wp - 135.6_wp*ptho(j,1)  + 5.2122_wp*ptho(j,1)**2 &
        & - 0.10939_wp*ptho(j,1)**3 + 0.00093777_wp*ptho(j,1)**4

        scn2o =  2356.2_wp - 166.38_wp*ptho(j,1)  + 6.3952_wp*ptho(j,1)**2 &
        & - 0.13422_wp*ptho(j,1)**3 + 0.0011506_wp*ptho(j,1)**4

        scdms =  2855.7_wp - 177.63_wp*ptho(j,1)  + 6.0438_wp*ptho(j,1)**2 &
        & - 0.11645_wp*ptho(j,1)**3 + 0.00094743_wp*ptho(j,1)**4

        !
        !  Compute the transfer (piston) velocity in m/s
        !  660 = Schmidt number of CO2 @ 20 degC in seawater

        ! Ignore ice for co2 here because coupling needs unscaled flux
        kwco2(j) =                      cmh2ms * pfu10( j)**2        &
        &           * (660._wp / scco2)**0.5_wp

        kwo2  = (1._wp - psicomo( j)) * cmh2ms * pfu10( j)**2        &
        &           * (660._wp / sco2)**0.5_wp

        kwdms = (1._wp - psicomo(j)) * cmh2ms * pfu10( j)**2        &
        &           * (660._wp / scdms)**0.5_wp

        kwn2o = (1._wp - psicomo(j)) * cmh2ms * pfu10( j)**2        &
        &           * (660._wp / scn2o)**0.5_wp

        IF(l_cpl_co2) THEN
            atco2(j) = local_bgc_mem%atm(j,iatmco2)
        END IF

        IF(.NOT. l_cpl_co2) THEN
            atco2(j) = atm_co2
        END IF

        ato2  = atm_o2
        atn2  = atm_n2

        !*********************************************************************
        !
        ! Calculate air-sea exchange for O2, N2, N2O
        !
        !*********************************************************************

        ! Surface flux of oxygen
        ! (Meiner-Reimer et. al, 2005, Eq. 74)

        oxflux = kwo2 * dtbgc * (local_bgc_mem%bgctra(j,1,ioxygen)                &
        &  -local_bgc_mem%satoxy(j,1) * (ato2 / 196800._wp)) ! *ppao(i,j)/101300. ! sea level pressure normalization

        local_bgc_mem%bgcflux(j,koflux) = oxflux/dtbgc

        local_bgc_mem%bgctra(j,1,ioxygen) = local_bgc_mem%bgctra(j,1,ioxygen)                 &
        &                - oxflux/(pddpo(j,1)+za(j))

        ! Surface flux of gaseous nitrogen (same piston velocity as for O2)
        ! (Meiner-Reimer et. al, 2005, Eq. 75)

        niflux = kwo2 * dtbgc * (local_bgc_mem%bgctra(j,1,igasnit)                &
        & -local_bgc_mem%satn2(j)*(atn2/802000._wp)) ! *ppao(i,j)/101300.

        local_bgc_mem%bgcflux(j,knflux) = niflux/dtbgc

        local_bgc_mem%bgctra(j,1,igasnit) = local_bgc_mem%bgctra(j,1,igasnit)                   &
        &                - niflux/(pddpo(j,1)+za(j))

        ! Surface flux of laughing gas (same piston velocity as for O2 and N2)
        ! (Meiner-Reimer et. al, 2005, Eq. 76)

        nlaughflux = kwn2o * dtbgc * (local_bgc_mem%bgctra(j,1,ian2o)              &
        &     - local_bgc_mem%satn2o(j))

        local_bgc_mem%bgctra(j,1,ian2o) = local_bgc_mem%bgctra(j,1,ian2o)                     &
        &              - nlaughflux/(pddpo(j,1)+za(j))
        local_bgc_mem%bgcflux(j,kn2oflux) = nlaughflux/dtbgc

        ! Surface flux of dms

        ! (Meiner-Reimer et. al, 2005, Eq. 77)

        dmsflux = kwdms*dtbgc*local_bgc_mem%bgctra(j,1,idms)

        local_bgc_mem%bgctra(j,1,idms) = local_bgc_mem%bgctra(j,1,idms) - dmsflux/pddpo(j,1)

        local_bgc_mem%bgcflux(j,kdmsflux) = dmsflux/dtbgc

    END IF

  END DO

#if 0

  !NEC$ nomove
  DO j = start_idx, end_idx

    IF( vmask(j) ) THEN

        !*********************************************************************
        !
        ! Calculate air sea exchange for CO2
        !
        !*********************************************************************

        ! Update local_bgc_mem%hi

        local_bgc_mem%hi(j,1) = update_hi(local_bgc_mem%hi(j,1), local_bgc_mem%bgctra(j,1,isco212), local_bgc_mem%aksurf(j,1) , &
        &  local_bgc_mem%aksurf(j,2),  local_bgc_mem%aksurf(j,4), local_bgc_mem%aksurf(j,7), local_bgc_mem%aksurf(j,6), local_bgc_mem%aksurf(j,5),&
        &  local_bgc_mem%aksurf(j,8), local_bgc_mem%aksurf(j,9),local_bgc_mem%aksurf(j,10), psao(j,1) , local_bgc_mem%aksurf(j,3), &
        &  local_bgc_mem%bgctra(j,1,isilica), local_bgc_mem%bgctra(j,1,iphosph),local_bgc_mem%bgctra(j,1,ialkali) )

    END IF

  END DO

#else

    CALL update_hi_sub(local_bgc_mem%hi(start_idx:end_idx,1), local_bgc_mem%bgctra(start_idx:end_idx,1,isco212), &
        &  local_bgc_mem%aksurf(start_idx:end_idx,1), local_bgc_mem%aksurf(start_idx:end_idx,2), local_bgc_mem%aksurf(start_idx:end_idx,4), &
        &  local_bgc_mem%aksurf(start_idx:end_idx,7), local_bgc_mem%aksurf(start_idx:end_idx,6), local_bgc_mem%aksurf(start_idx:end_idx,5), &
        &  local_bgc_mem%aksurf(start_idx:end_idx,8), local_bgc_mem%aksurf(start_idx:end_idx,9), local_bgc_mem%aksurf(start_idx:end_idx,10), &
        &  psao(start_idx:end_idx,1) , local_bgc_mem%aksurf(start_idx:end_idx,3), local_bgc_mem%bgctra(start_idx:end_idx,1,isilica), &
        &  local_bgc_mem%bgctra(start_idx:end_idx,1,iphosph), local_bgc_mem%bgctra(start_idx:end_idx,1,ialkali), &
        &  local_bgc_mem%hi(start_idx:end_idx,1), start_idx, end_idx, vmask(start_idx:end_idx))

#endif

  !NEC$ nomove
  DO j = start_idx, end_idx

    IF( vmask(j) ) THEN

        !
        ! Calculate pCO2 [ppmv] from total dissolved inorganic carbon (DIC: SCO212)
        ! the calculation also includes solubility
        !
        pco2=  local_bgc_mem%bgctra(j,1,isco212)  /((1._wp + local_bgc_mem%aksurf(j,1) * (1._wp + &
        & local_bgc_mem%aksurf(j,2)/local_bgc_mem%hi(j,1))/local_bgc_mem%hi(j,1)) * local_bgc_mem%solco2(j))

        fluxd=atco2(j)*kwco2(j)*dtbgc*local_bgc_mem%solco2(j) !
        fluxu=pco2 *kwco2(j)*dtbgc*local_bgc_mem%solco2(j) !

        ! new concentrations ocean (kmol/m3 -->ppm)
        thickness(j) = pddpo(j,1) + za(j)
        local_bgc_mem%bgctra(j,1,isco212) = local_bgc_mem%bgctra(j,1,isco212)+ (1._wp - psicomo(j)) * (fluxd-fluxu)/thickness(j)
        local_bgc_mem%bgcflux(j,kcflux) = (1._wp - psicomo(j)) * (fluxu-fluxd)/dtbgc
        local_bgc_mem%bgcflux(j,kcflux_cpl) = (fluxu-fluxd)/dtbgc
        local_bgc_mem%bgcflux(j,kpco2) = pco2

    END IF

  END DO

  IF (l_N_cycle) THEN

    ! RSE: INACTIVE CODE BLOCK
    ! WRITE(0,*) "ERROR: The program has entered an untested block of code."
    ! STOP

    !NEC$ nomove
    DO j = start_idx, end_idx

        IF( vmask(j) ) THEN

            ! Surface flux of ammonia  ! taken from Johnson et al, GBC,2008
            ! with atm. NH3 set to zero F = kgammo*KH_nh3 * NH3_seawater
            ! with NH3_seatwater = NH4* Ka/(Ka+hi) (Ka dissociation coef.)

            ! gas phase tranfer velocity
            kgammo = (1._wp - psicomo(j))*pfu10(j)/kg_denom

            ! Henry's law coefficient
            tabs = ptho(j,1) + 273.15_wp
            ecoef = 4092._wp/tabs - 9.70_wp
            kh_nh3i = 17.93_wp*(tabs/273.15_wp)*exp(ecoef)
            kh_nh3 = 1./kh_nh3i

            pka_nh3 = -0.467_wp + 0.00113_wp*psao(j,1) + 2887.9_wp/tabs
            ka_nh3 = 10._wp**(-pka_nh3)

            ! NH3 in seawater
            nh3sw = local_bgc_mem%bgctra(j,1,iammo)*ka_nh3/(ka_nh3+local_bgc_mem%hi(j,1))

            ammoflux = max(0._wp, dtbgc*kgammo*kh_nh3*nh3sw)
            local_bgc_mem%bgctra(j,1,iammo) = local_bgc_mem%bgctra(j,1,iammo) - ammoflux/thickness(j)

            ! LR: from mpiom, do not know what this is for ?!
            ! atm(i,j,iatmn2) = atm(i,j,iatmn2) + ammoflux*contppm/2._wp   !closing mass balance

            ! LR: from mpiom, don't think this is needed
            ! nh3flux(j) = nh3flux(i,j) + ammoflux  ! LR: from mpiom

            local_bgc_mem%bgcflux(j,knh3flux) = ammoflux/dtbgc
        END IF

    END DO

  END IF ! l_N_cycle

END SUBROUTINE
END MODULE mo_bgc_surface
