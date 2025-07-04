! ********************************************************************************************************
! ********************************************************************************************************
! ********************************************************************************************************
!
! Module updated during the natESM sprint #19
! 1. Moved the declaration of the iteration counter niter_atgen inside the solve_at_general subroutine
!    to prevent race conditions.
! 2. Unifying the processing of the surface layer and subsurface layers into a single loop within
!    the subroutine calc_dissol.
! 3. Changed the nesting order of the loops over layers and cells in the subroutine calc_dissol.
!    OpenACC directives have been added, but the correctness of GPU offloading has not been verified.
!
! ********************************************************************************************************
! ********************************************************************************************************
! ********************************************************************************************************

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

! @brief Inorganic carbon cycle.
!        Calc dissolution, update of hydrogen ions

MODULE mo_carchm

USE mo_bgc_memory_types, ONLY  : t_bgc_memory

USE mo_memory_bgc, ONLY     : rrrcl

USE mo_kind, ONLY           : wp
USE mo_control_bgc, ONLY    : bgc_nproma, bgc_zlevs
USE mo_param1_bgc, ONLY     : icalc, ialkali, isco212, isilica, iphosph, klysocl
USE mo_hamocc_nml, ONLY     : hion_solver, dremcalc
USE mo_fortran_tools, ONLY  : set_acc_host_or_device

IMPLICIT NONE

PRIVATE :: anw_infsup, equation_at, ahini_for_at, solve_at_general
PRIVATE :: solve_at_general_sub

PUBLIC  :: calc_dissol, update_hi
PUBLIC  :: update_hi_sub

!===============================================================================
! The Parameters needed for the mocsy Funtions/Subroutines
!===============================================================================
! General parameters
REAL(wp), PARAMETER :: pp_rdel_ah_target = 1.E-8_wp
REAL(wp), PARAMETER :: pp_ln10 = 2.302585092994045684018_wp

! Maximum number of iterations for each method
INTEGER, PARAMETER :: jp_maxniter_atgen    = 50

!===============================================================================


CONTAINS

SUBROUTINE calc_dissol (local_bgc_mem, start_idx, end_idx, klevs, pddpo, psao, ptiestu, lacc)

!! Computes calcium carbonate dissolution

  IMPLICIT NONE

  !! Arguments
  TYPE(t_bgc_memory), POINTER    :: local_bgc_mem

  INTEGER, INTENT(in) :: start_idx             !< start index for j loop (ICON cells, MPIOM lat dir)
  INTEGER, INTENT(in) :: end_idx               !< end index  for j loop  (ICON cells, MPIOM lat dir)
  INTEGER, INTENT(in) :: klevs(bgc_nproma)     !<  vertical levels


  REAL(wp),INTENT(in) :: pddpo(bgc_nproma,bgc_zlevs) !< size of scalar grid cell (3rd REAL) [m]
  REAL(wp),INTENT(in) :: psao(bgc_nproma,bgc_zlevs)  !< salinity
  REAL(wp),INTENT(in) :: ptiestu(bgc_nproma,bgc_zlevs)  !< depth of scalar grid cell [m]
  LOGICAL, INTENT(IN), OPTIONAL :: lacc

  !! Local variables

  INTEGER :: k, j, kpke

  REAL(wp) :: supsat, undsa, dissol
  REAL(wp) :: supsatup,satdiff,depthdiff   ! needed to calculate depth of lysocline
  INTEGER  :: iflag(start_idx:end_idx), kpkeMAX
  LOGICAL  :: lzacc

  !! Vectorization local variables
  LOGICAL :: vmask(start_idx:end_idx)

  CALL set_acc_host_or_device(lzacc, lacc)

  !
  !*********************************************************************
  !
  ! Dissolution of calcite, whole water column
  !
  !*********************************************************************

  iflag(:) = 0  
  kpkeMAX = MAXVAL(klevs(start_idx:end_idx))

  !  Dissolution in surface layer and subsurface layers
  !NEC$ nomove
  DO k = 1, kpkeMAX

    DO j = start_idx, end_idx
        IF( (k <= klevs(j)) .AND. (pddpo(j,k) > EPSILON(0.5_wp)) ) THEN
            vmask(j) = .TRUE.
        ELSE
            vmask(j) = .FALSE.
        END IF
    END DO

#if 0
    DO j = start_idx, end_idx

        IF( vmask(j) ) THEN

            CALL update_hi_sub(local_bgc_mem%hi(j,k), local_bgc_mem%bgctra(j,k,isco212), local_bgc_mem%ak13(j,k) , &
            &    local_bgc_mem%ak23(j,k), local_bgc_mem%akw3(j,k),local_bgc_mem%aks3(j,k),local_bgc_mem%akf3(j,k), local_bgc_mem%aksi3(j,k),&
            &    local_bgc_mem%ak1p3(j,k),local_bgc_mem%ak2p3(j,k),local_bgc_mem%ak3p3(j,k),psao(j,k) , local_bgc_mem%akb3(j,k), &
            &    local_bgc_mem%bgctra(j,k,isilica),local_bgc_mem%bgctra(j,k,iphosph),local_bgc_mem%bgctra(j,k,ialkali), local_bgc_mem%hi(j,k) )

        END IF ! vmask(j) == .TRUE.

    END DO ! j
#else

    CALL update_hi_sub(local_bgc_mem%hi(start_idx:end_idx,k), local_bgc_mem%bgctra(start_idx:end_idx,k,isco212), &
    & local_bgc_mem%ak13(start_idx:end_idx,k), local_bgc_mem%ak23(start_idx:end_idx,k), &
    & local_bgc_mem%akw3(start_idx:end_idx,k), local_bgc_mem%aks3(start_idx:end_idx,k), &
    & local_bgc_mem%akf3(start_idx:end_idx,k), local_bgc_mem%aksi3(start_idx:end_idx,k), &
    & local_bgc_mem%ak1p3(start_idx:end_idx,k), local_bgc_mem%ak2p3(start_idx:end_idx,k), local_bgc_mem%ak3p3(start_idx:end_idx,k), &
    & psao(start_idx:end_idx,k), local_bgc_mem%akb3(start_idx:end_idx,k), local_bgc_mem%bgctra(start_idx:end_idx,k,isilica), &
    & local_bgc_mem%bgctra(start_idx:end_idx,k,iphosph), local_bgc_mem%bgctra(start_idx:end_idx,k,ialkali), &
    & local_bgc_mem%hi(start_idx:end_idx,k), start_idx, end_idx, vmask(start_idx:end_idx))

#endif

    !NEC$ nomove
    DO j = start_idx, end_idx

        IF( vmask(j) ) THEN

            local_bgc_mem%co3(j,k) = local_bgc_mem%bgctra(j,k,isco212)/(1._wp+local_bgc_mem%hi(j,k) * &
            & (1._wp+local_bgc_mem%hi(j,k)/local_bgc_mem%ak13(j,k))/local_bgc_mem%ak23(j,k))

            supsat = local_bgc_mem%co3(j,k)-97._wp*local_bgc_mem%aksp(j,k)   ! 97. = 1./1.03e-2 (MEAN TOTAL [CA++] IN SEAWATER [kmol/m3])
            undsa  = MAX(0._wp, -supsat)

            dissol = MIN(undsa,dremcalc*local_bgc_mem%bgctra(j,k,icalc))
            local_bgc_mem%bgctra(j,k,icalc)   = local_bgc_mem%bgctra(j,k,icalc)-dissol
            local_bgc_mem%bgctra(j,k,ialkali) = local_bgc_mem%bgctra(j,k,ialkali)+2._wp*dissol

            local_bgc_mem%bgctra(j,k,isco212) = local_bgc_mem%bgctra(j,k,isco212)+dissol

            IF (k == 1 .AND. supsat < 0._wp .AND. iflag(j) == 0) THEN
                iflag(j) = 1
                local_bgc_mem%bgcflux(j,klysocl) = ptiestu(j,1)
            END IF

            IF (k >  1 .AND. supsat < 0._wp .AND. iflag(j) == 0) THEN
                iflag(j) = 1
                supsatup  = local_bgc_mem%co3(j,k-1)-97._wp*local_bgc_mem%aksp(j,k-1)
                depthdiff = 0.5_wp * (pddpo(j,k)+pddpo(j,k-1))
                satdiff   = supsatup-supsat
                local_bgc_mem%bgcflux(j,klysocl) = ptiestu(j,k-1)+depthdiff*(supsatup/satdiff)  ! depth of lysokline
            END IF

        END IF ! vmask(j) == .TRUE.

    END DO ! j

  END DO ! k

END SUBROUTINE


FUNCTION update_hi(hi,c,ak1,ak2,akw,aks,akf,aksi,ak1p,ak2p,ak3p,s,akb,sit,pt,alk) RESULT(h)
 !$ACC ROUTINE SEQ
 ! update hydrogen ion concentration

 REAL(wp), INTENT(IN)  :: hi, c, ak1, ak2, akw, aks, akf, aksi, ak1p, ak2p, ak3p, s, akb, sit, pt, alk
 REAL(wp)              :: h
 ! REAL(wp), INTENT(OUT) :: h

 ! LOCAL
 REAL(wp) :: ah1, bt, sti, ft, hso4, hf, hsi, hpo4, ab, aw, ac, ah2o, ah2, erel
 INTEGER  :: iter, jit

 bt    = rrrcl*s
! sulfate Morris & Riley (1966)
 sti   = 0.14_wp *  s*1.025_wp/1.80655_wp  / 96.062_wp
! fluoride Riley (1965)
 ft    = 0.000067_wp * s*1.025_wp/1.80655_wp / 18.9984_wp

IF (hion_solver == 0) THEN

  ah1 = hi + epsilon(1._wp)
  iter  = 0

  DO jit = 1,20

     hso4 = sti / ( 1._wp + aks / ( ah1 / ( 1._wp + sti / aks ) ) )
     hf   = 1._wp / ( 1._wp + akf / ah1 )
     hsi  = 1._wp/ ( 1._wp + ah1 / aksi )
     hpo4 = ( ak1p * ak2p * ( ah1 + 2._wp * ak3p ) - ah1**3 ) /    &
           &        ( ah1**3 + ak1p * ah1**2 + ak1p * ak2p * ah1 + ak1p * ak2p*ak3p )
     ab   = bt / ( 1._wp + ah1 / akb )
     aw   = akw / ah1 - ah1 / ( 1._wp + sti / aks )
     ac   = alk + hso4 - sit * hsi - ab - aw + ft * hf - pt * hpo4
     ah2o = SQRT( ( c - ac )**2 + 4._wp * ( ac * ak2 / ak1 ) * ( 2._wp*c - ac ) )
     ah2  = 0.5_wp * ak1 / ac *( ( c - ac ) + ah2o )
     erel = ( ah2 - ah1 ) / ah2

     IF (abs( erel ).ge.5.e-5_wp) THEN
        ah1 = ah2
        iter = iter + 1
     ELSE
        ah1 = ah2
        EXIT
     END IF

   ENDDO

   IF(ah1.gt.0._wp) THEN
    h = max(1.e-20_wp,ah1)
   ELSE
    ! WARNING: Value is undefined
   END IF

END IF

IF (hion_solver == 1) THEN

  h = solve_at_general(alk,c,bt,pt,sit,sti,ft,ak1,ak2,akb,akw,aks,akf,ak1p,ak2p,ak3p,aksi,hi)

END IF

END FUNCTION

SUBROUTINE update_hi_sub(hi,c,ak1,ak2,akw,aks,akf,aksi,ak1p,ak2p,ak3p,s,akb,sit,pt,alk,h,start_idx,end_idx,vmask)
!$ACC ROUTINE SEQ
! update hydrogen ion concentration

    INTEGER, INTENT(IN)   :: start_idx, end_idx
    REAL(wp), INTENT(IN)  :: hi(start_idx:end_idx), c(start_idx:end_idx), ak1(start_idx:end_idx), ak2(start_idx:end_idx)
    REAL(wp), INTENT(IN)  :: akw(start_idx:end_idx), aks(start_idx:end_idx), akf(start_idx:end_idx), aksi(start_idx:end_idx)
    REAL(wp), INTENT(IN)  :: ak1p(start_idx:end_idx), ak2p(start_idx:end_idx), ak3p(start_idx:end_idx), s(start_idx:end_idx)
    REAL(wp), INTENT(IN)  :: akb(start_idx:end_idx), sit(start_idx:end_idx), pt(start_idx:end_idx), alk(start_idx:end_idx)
    LOGICAL, INTENT(IN)   :: vmask(start_idx:end_idx)
    REAL(wp), INTENT(OUT) :: h(start_idx:end_idx)

    ! LOCAL
    REAL(wp) :: bt(start_idx:end_idx), sti(start_idx:end_idx), ft(start_idx:end_idx)
    INTEGER  :: j

    DO j = start_idx, end_idx

        IF( vmask(j) ) THEN
            bt(j)    = rrrcl*s(j)
            ! sulfate Morris & Riley (1966)
            sti(j)   = 0.14_wp *  s(j)*1.025_wp/1.80655_wp  / 96.062_wp
            ! fluoride Riley (1965)
            ft(j)    = 0.000067_wp * s(j)*1.025_wp/1.80655_wp / 18.9984_wp
        END IF

    END DO

#if 0

    bt    = rrrcl*s
    ! sulfate Morris & Riley (1966)
    sti   = 0.14_wp *  s*1.025_wp/1.80655_wp  / 96.062_wp
    ! fluoride Riley (1965)
    ft    = 0.000067_wp * s*1.025_wp/1.80655_wp / 18.9984_wp

    IF (hion_solver == 0) THEN

        ah1 = hi + epsilon(1._wp)
        iter  = 0

        DO jit = 1,20

            hso4 = sti / ( 1._wp + aks / ( ah1 / ( 1._wp + sti / aks ) ) )
            hf   = 1._wp / ( 1._wp + akf / ah1 )
            hsi  = 1._wp/ ( 1._wp + ah1 / aksi )
            hpo4 = ( ak1p * ak2p * ( ah1 + 2._wp * ak3p ) - ah1**3 ) /    &
            &      ( ah1**3 + ak1p * ah1**2 + ak1p * ak2p * ah1 + ak1p * ak2p*ak3p )
            ab   = bt / ( 1._wp + ah1 / akb )
            aw   = akw / ah1 - ah1 / ( 1._wp + sti / aks )
            ac   = alk + hso4 - sit * hsi - ab - aw + ft * hf - pt * hpo4
            ah2o = SQRT( ( c - ac )**2 + 4._wp * ( ac * ak2 / ak1 ) * ( 2._wp*c - ac ) )
            ah2  = 0.5_wp * ak1 / ac *( ( c - ac ) + ah2o )
            erel = ( ah2 - ah1 ) / ah2

            IF (abs( erel ).ge.5.e-5_wp) THEN
                ah1 = ah2
                iter = iter + 1
            ELSE
                ah1 = ah2
                EXIT
            END IF

        ENDDO

        IF(ah1.gt.0._wp) h = max(1.e-20_wp,ah1)

    END IF

#else

    IF (hion_solver == 0) THEN
        WRITE(0,*) "ERROR: The program entered a code path that was not vectorized"
        STOP
    END IF

#endif

    IF (hion_solver == 1) THEN

#if 0
        DO j = start_idx, end_idx

            IF( vmask(j) ) THEN
                CALL solve_at_general_sub(h(j), alk(j),c(j),bt(j),pt(j),sit(j),sti(j),ft(j),ak1(j),ak2(j), &
                &                         akb(j),akw(j),aks(j),akf(j),ak1p(j),ak2p(j),ak3p(j),aksi(j),hi(j))
            END IF

        END DO
#else
        CALL solve_at_general_sub(start_idx, end_idx, vmask(start_idx:end_idx), h(start_idx:end_idx), &
        & alk(start_idx:end_idx), c(start_idx:end_idx), bt(start_idx:end_idx), pt(start_idx:end_idx), &
        & sit(start_idx:end_idx), sti(start_idx:end_idx), ft(start_idx:end_idx), ak1(start_idx:end_idx), ak2(start_idx:end_idx), &
        & akb(start_idx:end_idx), akw(start_idx:end_idx), aks(start_idx:end_idx), akf(start_idx:end_idx), &
        & ak1p(start_idx:end_idx), ak2p(start_idx:end_idx), ak3p(start_idx:end_idx), aksi(start_idx:end_idx), hi(start_idx:end_idx))
#endif

    END IF

END SUBROUTINE

!===============================================================================
! Routines from the mocsy package
!===============================================================================

SUBROUTINE anw_infsup(p_dictot, p_bortot,                                     &
                      p_po4tot, p_siltot,                                     &
                      p_so4tot, p_flutot,                                     &
                      p_alknw_inf, p_alknw_sup)
!$ACC ROUTINE SEQ

! Subroutine returns the lower and upper bounds of "non-water-selfionization"
! contributions to total alkalinity (the infimum and the supremum), i.e
! inf(TA - [OH-] + [H+]) and sup(TA - [OH-] + [H+])

IMPLICIT NONE

! Argument variables
REAL(wp), INTENT(IN)  :: p_dictot
REAL(wp), INTENT(IN)  :: p_bortot
REAL(wp), INTENT(IN)  :: p_po4tot
REAL(wp), INTENT(IN)  :: p_siltot
REAL(wp), INTENT(IN)  :: p_so4tot
REAL(wp), INTENT(IN)  :: p_flutot
REAL(wp), INTENT(OUT) :: p_alknw_inf
REAL(wp), INTENT(OUT) :: p_alknw_sup

p_alknw_inf =  -p_po4tot - p_so4tot - p_flutot
p_alknw_sup =   p_dictot + p_dictot + p_bortot &
              + p_po4tot + p_po4tot + p_siltot !&
!             + p_nh4tot + p_h2stot

RETURN
END SUBROUTINE anw_infsup

!===============================================================================

FUNCTION equation_at(p_alktot, p_h,       p_dictot, p_bortot,                 &
                     p_po4tot, p_siltot,                                      &
                     p_so4tot, p_flutot,                                      &
                     K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,              &
                     p_deriveqn)
!$ACC ROUTINE SEQ

  ! Purpose: Compute total alkalinity from ion concentrations and equilibrium constants

IMPLICIT NONE
REAL(wp) :: equation_at

! Argument variables
REAL(wp), INTENT(IN)            :: p_alktot
REAL(wp), INTENT(IN)            :: p_h
REAL(wp), INTENT(IN)            :: p_dictot
REAL(wp), INTENT(IN)            :: p_bortot
REAL(wp), INTENT(IN)            :: p_po4tot
REAL(wp), INTENT(IN)            :: p_siltot
REAL(wp), INTENT(IN)            :: p_so4tot
REAL(wp), INTENT(IN)            :: p_flutot
REAL(wp), INTENT(IN)            :: K1, K2, Kb, Kw, Ks, Kf
REAL(wp), INTENT(IN)            :: K1p, K2p, K3p, Ksi
REAL(wp), INTENT(OUT), OPTIONAL :: p_deriveqn

! Local variables
!-----------------
REAL(wp) :: znumer_dic, zdnumer_dic, zdenom_dic, zalk_dic, zdalk_dic
REAL(wp) :: znumer_bor, zdnumer_bor, zdenom_bor, zalk_bor, zdalk_bor
REAL(wp) :: znumer_po4, zdnumer_po4, zdenom_po4, zalk_po4, zdalk_po4
REAL(wp) :: znumer_sil, zdnumer_sil, zdenom_sil, zalk_sil, zdalk_sil
REAL(wp) :: znumer_nh4, zdnumer_nh4, zdenom_nh4, zalk_nh4, zdalk_nh4
REAL(wp) :: znumer_h2s, zdnumer_h2s, zdenom_h2s, zalk_h2s, zdalk_h2s
REAL(wp) :: znumer_so4, zdnumer_so4, zdenom_so4, zalk_so4, zdalk_so4
REAL(wp) :: znumer_flu, zdnumer_flu, zdenom_flu, zalk_flu, zdalk_flu
REAL(wp) ::                                      zalk_wat, zdalk_wat
REAL(wp) :: aphscale

! TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
aphscale = 1._wp + p_so4tot/Ks

! H2CO3 - HCO3 - CO3 : n=2, m=0
znumer_dic = 2._wp*K1*K2 + p_h*       K1
zdenom_dic =       K1*K2 + p_h*(      K1 + p_h)
zalk_dic   = p_dictot * (znumer_dic/zdenom_dic)

! B(OH)3 - B(OH)4 : n=1, m=0
znumer_bor =       Kb
zdenom_bor =       Kb + p_h
zalk_bor   = p_bortot * (znumer_bor/zdenom_bor)

! H3PO4 - H2PO4 - HPO4 - PO4 : n=3, m=1
znumer_po4 = 3._wp*K1p*K2p*K3p + p_h*(2._wp*K1p*K2p + p_h* K1p)
zdenom_po4 =       K1p*K2p*K3p + p_h*(      K1p*K2p + p_h*(K1p + p_h))
zalk_po4   = p_po4tot * (znumer_po4/zdenom_po4 - 1._wp) ! Zero level of H3PO4 = 1

! H4SiO4 - H3SiO4 : n=1, m=0
znumer_sil =       Ksi
zdenom_sil =       Ksi + p_h
zalk_sil   = p_siltot * (znumer_sil/zdenom_sil)

! HSO4 - SO4 : n=1, m=1
znumer_so4 =       Ks
zdenom_so4 =       Ks + p_h
zalk_so4   = p_so4tot * (znumer_so4/zdenom_so4 - 1._wp)

! HF - F : n=1, m=1
znumer_flu =       Kf
zdenom_flu =       Kf + p_h
zalk_flu   = p_flutot * (znumer_flu/zdenom_flu - 1._wp)

! H2O - OH
zalk_wat   = Kw/p_h - p_h/aphscale

equation_at =    zalk_dic + zalk_bor + zalk_po4 + zalk_sil &
               + zalk_so4 + zalk_flu                       &
               + zalk_wat - p_alktot

IF(PRESENT(p_deriveqn)) THEN
   ! H2CO3 - HCO3 - CO3 : n=2
   zdnumer_dic = K1*K1*K2 + p_h*(4._wp*K1*K2               &
                          + p_h*       K1    )
   zdalk_dic   = -p_dictot*(zdnumer_dic/zdenom_dic**2)

   ! B(OH)3 - B(OH)4 : n=1
   zdnumer_bor = Kb
   zdalk_bor   = -p_bortot*(zdnumer_bor/zdenom_bor**2)

   ! H3PO4 - H2PO4 - HPO4 - PO4 : n=3
   zdnumer_po4 = K1p*K2p*K1p*K2p*K3p + p_h*(4._wp*K1p*K1p*K2p*K3p                &
                                     + p_h*(9._wp*K1p*K2p*K3p + K1p*K1p*K2p      &
                                     + p_h*(4._wp*K1p*K2p                        &
                                     + p_h*       K1p)))
   zdalk_po4   = -p_po4tot * (zdnumer_po4/zdenom_po4**2)

   ! H4SiO4 - H3SiO4 : n=1
   zdnumer_sil = Ksi
   zdalk_sil   = -p_siltot * (zdnumer_sil/zdenom_sil**2)

!  ! NH4 - NH3 : n=1
!  zdnumer_nh4 = Knh4
!  zdalk_nh4   = -p_nh4tot * (zdnumer_nh4/zdenom_nh4**2)

!  ! H2S - HS : n=1
!  zdnumer_h2s = api1_h2s
!  zdalk_h2s   = -p_h2stot * (zdnumer_h2s/zdenom_h2s**2)

   ! HSO4 - SO4 : n=1
   zdnumer_so4 = Ks
   zdalk_so4   = -p_so4tot * (zdnumer_so4/zdenom_so4**2)

   ! HF - F : n=1
   zdnumer_flu = Kf
   zdalk_flu   = -p_flutot * (zdnumer_flu/zdenom_flu**2)

!  p_deriveqn =   zdalk_dic + zdalk_bor + zdalk_po4 + zdalk_sil &
!               + zdalk_nh4 + zdalk_h2s + zdalk_so4 + zdalk_flu &
!               - Kw/p_h**2 - 1._wp/aphscale
   p_deriveqn =   zdalk_dic + zdalk_bor + zdalk_po4 + zdalk_sil &
                                        + zdalk_so4 + zdalk_flu &
                - Kw/p_h**2 - 1._wp/aphscale
ENDIF
RETURN
END FUNCTION equation_at


!===============================================================================

SUBROUTINE ahini_for_at(p_alkcb, p_dictot, p_bortot, K1, K2, Kb, p_hini)
!$ACC ROUTINE SEQ

! Subroutine returns the root for the 2nd order approximation of the
! DIC -- B_T -- A_CB equation for [H+] (reformulated as a cubic polynomial)
! around the local minimum, if it exists.

! Returns * 1E-03_wp if p_alkcb <= 0
!         * 1E-10_wp if p_alkcb >= 2*p_dictot + p_bortot
!         * 1E-07_wp if 0 < p_alkcb < 2*p_dictot + p_bortot
!                    and the 2nd order approximation does not have a solution

!USE MOD_CHEMCONST, ONLY : api1_dic, api2_dic, api1_bor

IMPLICIT NONE

! Argument variables
!--------------------
REAL(wp), INTENT(IN)   ::  p_alkcb, p_dictot, p_bortot
REAL(wp), INTENT(IN)   ::  K1, K2, Kb
REAL(wp), INTENT(OUT)  ::  p_hini

! Local variables
!-----------------
REAL(wp)  ::  zca, zba
REAL(wp)  ::  zd, zsqrtd, zhmin
REAL(wp)  ::  za2, za1, za0

IF (p_alkcb <= 0._wp) THEN
  p_hini = 1.e-3_wp
ELSEIF (p_alkcb >= (2._wp*p_dictot + p_bortot)) THEN
  p_hini = 1.e-10_wp
ELSE
  zca = p_dictot/p_alkcb
  zba = p_bortot/p_alkcb

  ! Coefficients of the cubic polynomial
  za2 = Kb*(1._wp - zba) + K1*(1._wp-zca)
  za1 = K1*Kb*(1._wp - zba - zca) + K1*K2*(1._wp - (zca+zca))
  za0 = K1*K2*Kb*(1._wp - zba - (zca+zca))
                                        ! Taylor expansion around the minimum
  zd = za2*za2 - 3._wp*za1              ! Discriminant of the quadratic equation
                                        ! for the minimum close to the root

  IF(zd > 0._wp) THEN                   ! If the discriminant is positive
    zsqrtd = SQRT(zd)
    IF(za2 < 0) THEN
      zhmin = (-za2 + zsqrtd)/3._wp
    ELSE
      zhmin = -za1/(za2 + zsqrtd)
    ENDIF
    p_hini = zhmin + SQRT(-(za0 + zhmin*(za1 + zhmin*(za2 + zhmin)))/zsqrtd)
  ELSE
    p_hini = 1.e-7_wp
  ENDIF

ENDIF
RETURN
END SUBROUTINE ahini_for_at

!===============================================================================

FUNCTION solve_at_general(p_alktot, p_dictot, p_bortot,                       &
    p_po4tot, p_siltot,                                 &
    p_so4tot, p_flutot,                                 &
    K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,         &
    p_hini,   p_val)
!$ACC ROUTINE SEQ

! Purpose: Compute [H degree] ion concentration from sea-water ion concentrations,
!          alkalinity, DIC, and equilibrium constants
! Universal pH solver that converges from any given initial value,
! determines upper an lower bounds for the solution if required

IMPLICIT NONE
REAL(wp) :: SOLVE_AT_GENERAL

! Argument variables
!--------------------
REAL(wp), INTENT(IN)            :: p_alktot
REAL(wp), INTENT(IN)            :: p_dictot
REAL(wp), INTENT(IN)            :: p_bortot
REAL(wp), INTENT(IN)            :: p_po4tot
REAL(wp), INTENT(IN)            :: p_siltot
!REAL(wp), INTENT(IN)            :: p_nh4tot
!REAL(wp), INTENT(IN)            :: p_h2stot
REAL(wp), INTENT(IN)            :: p_so4tot
REAL(wp), INTENT(IN)            :: p_flutot
REAL(wp), INTENT(IN)            :: K1, K2, Kb, Kw, Ks, Kf
REAL(wp), INTENT(IN)            :: K1p, K2p, K3p, Ksi
REAL(wp), INTENT(IN), OPTIONAL  :: p_hini
REAL(wp), INTENT(OUT), OPTIONAL :: p_val

! Local variables
!-----------------
REAL(wp)  ::  zh_ini, zh, zh_prev, zh_lnfactor
REAL(wp)  ::  zalknw_inf, zalknw_sup
REAL(wp)  ::  zh_min, zh_max
REAL(wp)  ::  zdelta, zh_delta
REAL(wp)  ::  zeqn, zdeqndh, zeqn_absmin
REAL(wp)  ::  aphscale
LOGICAL   ::  l_exitnow
REAL(wp), PARAMETER :: pz_exp_threshold = 1.0_wp
INTEGER   ::  niter_atgen

! TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
aphscale = 1._wp + p_so4tot/Ks

IF(PRESENT(p_hini)) THEN
zh_ini = p_hini
ELSE
CALL ahini_for_at(p_alktot, p_dictot, p_bortot, K1, K2, Kb, zh_ini)
ENDIF

CALL anw_infsup(p_dictot, p_bortot,                                           &
p_po4tot, p_siltot,                                           &
p_so4tot, p_flutot,                                           &
zalknw_inf, zalknw_sup)

zdelta = (p_alktot-zalknw_inf)**2 + 4._wp*Kw/aphscale

IF(p_alktot >= zalknw_inf) THEN
zh_min = 2._wp*Kw /( p_alktot-zalknw_inf + SQRT(zdelta) )
ELSE
zh_min = aphscale*(-(p_alktot-zalknw_inf) + SQRT(zdelta) ) / 2._wp
ENDIF

zdelta = (p_alktot-zalknw_sup)**2 + 4._wp*Kw/aphscale

IF(p_alktot <= zalknw_sup) THEN
zh_max = aphscale*(-(p_alktot-zalknw_sup) + SQRT(zdelta) ) / 2._wp
ELSE
zh_max = 2._wp*Kw /( p_alktot-zalknw_sup + SQRT(zdelta) )
ENDIF

zh = MAX(MIN(zh_max, zh_ini), zh_min)
niter_atgen        = 0                 ! Reset counters of iterations
zeqn_absmin        = HUGE(1._wp)

!NEC$ nomove
DO
IF(niter_atgen >= jp_maxniter_atgen) THEN
zh = -1._wp
EXIT
ENDIF

zh_prev = zh
zeqn = equation_at(p_alktot, zh,       p_dictot, p_bortot,                  &
p_po4tot, p_siltot,                                      &
p_so4tot, p_flutot,                                      &
K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,              &
P_DERIVEQN = zdeqndh)

! Adapt bracketing interval
IF(zeqn > 0._wp) THEN
zh_min = zh_prev
ELSEIF(zeqn < 0._wp) THEN
zh_max = zh_prev
ELSE
! zh is the root; unlikely but, one never knows
EXIT
ENDIF

! Now determine the next iterate zh
niter_atgen = niter_atgen + 1

IF(ABS(zeqn) >= 0.5_wp*zeqn_absmin) THEN
! if the function evaluation at the current point is
! not decreasing faster than with a bisection step (at least linearly)
! in absolute value take one bisection step on [ph_min, ph_max]
! ph_new = (ph_min + ph_max)/2d0
!
! In terms of [H]_new:
! [H]_new = 10**(-ph_new)
!         = 10**(-(ph_min + ph_max)/2d0)
!         = SQRT(10**(-(ph_min + phmax)))
!         = SQRT(zh_max * zh_min)
zh = SQRT(zh_max * zh_min)
zh_lnfactor = (zh - zh_prev)/zh_prev ! Required to test convergence below
ELSE
! dzeqn/dpH = dzeqn/d[H] * d[H]/dpH
!           = -zdeqndh * LOG(10) * [H]
! \Delta pH = -zeqn/(zdeqndh*d[H]/dpH) = zeqn/(zdeqndh*[H]*LOG(10))
!
! pH_new = pH_old + \deltapH
!
! [H]_new = 10**(-pH_new)
!         = 10**(-pH_old - \Delta pH)
!         = [H]_old * 10**(-zeqn/(zdeqndh*[H]_old*LOG(10)))
!         = [H]_old * EXP(-LOG(10)*zeqn/(zdeqndh*[H]_old*LOG(10)))
!         = [H]_old * EXP(-zeqn/(zdeqndh*[H]_old))

zh_lnfactor = -zeqn/(zdeqndh*zh_prev)

IF(ABS(zh_lnfactor) > pz_exp_threshold) THEN
zh          = zh_prev*EXP(zh_lnfactor)
ELSE
zh_delta    = zh_lnfactor*zh_prev
zh          = zh_prev + zh_delta
ENDIF

IF( zh < zh_min ) THEN
! if [H]_new < [H]_min
! i.e., if ph_new > ph_max then
! take one bisection step on [ph_prev, ph_max]
! ph_new = (ph_prev + ph_max)/2d0
! In terms of [H]_new:
! [H]_new = 10**(-ph_new)
!         = 10**(-(ph_prev + ph_max)/2d0)
!         = SQRT(10**(-(ph_prev + phmax)))
!         = SQRT([H]_old*10**(-ph_max))
!         = SQRT([H]_old * zh_min)
zh                = SQRT(zh_prev * zh_min)
zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below
ENDIF

IF( zh > zh_max ) THEN
! if [H]_new > [H]_max
! i.e., if ph_new < ph_min, then
! take one bisection step on [ph_min, ph_prev]
! ph_new = (ph_prev + ph_min)/2d0
! In terms of [H]_new:
! [H]_new = 10**(-ph_new)
!         = 10**(-(ph_prev + ph_min)/2d0)
!         = SQRT(10**(-(ph_prev + ph_min)))
!         = SQRT([H]_old*10**(-ph_min))
!         = SQRT([H]_old * zhmax)
zh                = SQRT(zh_prev * zh_max)
zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below
ENDIF
ENDIF

zeqn_absmin = MIN( ABS(zeqn), zeqn_absmin)

! Stop iterations once |\delta{[H]}/[H]| < rdel
! <=> |(zh - zh_prev)/zh_prev| = |EXP(-zeqn/(zdeqndh*zh_prev)) -1| < rdel
! |EXP(-zeqn/(zdeqndh*zh_prev)) -1| ~ |zeqn/(zdeqndh*zh_prev)|

! Alternatively:
! |\Delta pH| = |zeqn/(zdeqndh*zh_prev*LOG(10))|
!             ~ 1/LOG(10) * |\Delta [H]|/[H]
!             < 1/LOG(10) * rdel

! Hence |zeqn/(zdeqndh*zh)| < rdel

! rdel <-- pp_rdel_ah_target

l_exitnow = (ABS(zh_lnfactor) < pp_rdel_ah_target)

IF(l_exitnow) EXIT
ENDDO

solve_at_general = zh

IF(PRESENT(p_val)) THEN
IF(zh > 0._wp) THEN
p_val = equation_at(p_alktot, zh,       p_dictot, p_bortot,              &
    p_po4tot, p_siltot,                                  &
    p_so4tot, p_flutot,                                  &
    K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi)
ELSE
p_val = HUGE(1._wp)
ENDIF
ENDIF
RETURN
END FUNCTION solve_at_general

SUBROUTINE solve_at_general_sub(start_idx, end_idx, vmask, resSAG,                  &
                                p_alktot, p_dictot, p_bortot,                       &
                                p_po4tot, p_siltot,                                 &
                                p_so4tot, p_flutot,                                 &
                                K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,         &
                                p_hini,   p_val)
!$ACC ROUTINE SEQ

IMPLICIT NONE

! Argument variables
!--------------------
INTEGER, INTENT(IN)             :: start_idx, end_idx
LOGICAL, INTENT(IN)             :: vmask(start_idx:end_idx)
REAL(wp), INTENT(OUT)           :: resSAG(start_idx:end_idx)
REAL(wp), INTENT(IN)            :: p_alktot(start_idx:end_idx)
REAL(wp), INTENT(IN)            :: p_dictot(start_idx:end_idx)
REAL(wp), INTENT(IN)            :: p_bortot(start_idx:end_idx)
REAL(wp), INTENT(IN)            :: p_po4tot(start_idx:end_idx)
REAL(wp), INTENT(IN)            :: p_siltot(start_idx:end_idx)
REAL(wp), INTENT(IN)            :: p_so4tot(start_idx:end_idx)
REAL(wp), INTENT(IN)            :: p_flutot(start_idx:end_idx)
REAL(wp), INTENT(IN)            :: K1(start_idx:end_idx), K2(start_idx:end_idx), Kb(start_idx:end_idx)
REAL(wp), INTENT(IN)            :: Kw(start_idx:end_idx), Ks(start_idx:end_idx), Kf(start_idx:end_idx)
REAL(wp), INTENT(IN)            :: K1p(start_idx:end_idx), K2p(start_idx:end_idx), K3p(start_idx:end_idx), Ksi(start_idx:end_idx)

REAL(wp), INTENT(IN), OPTIONAL  :: p_hini(start_idx:end_idx)
REAL(wp), INTENT(OUT), OPTIONAL :: p_val(start_idx:end_idx)

! Local variables
!-----------------
REAL(wp)  ::  zh_ini(start_idx:end_idx), zh(start_idx:end_idx), zh_prev(start_idx:end_idx), zh_lnfactor(start_idx:end_idx)
REAL(wp)  ::  zalknw_inf(start_idx:end_idx), zalknw_sup(start_idx:end_idx)
REAL(wp)  ::  zh_min(start_idx:end_idx), zh_max(start_idx:end_idx)
REAL(wp)  ::  zdelta(start_idx:end_idx), zh_delta(start_idx:end_idx)
REAL(wp)  ::  zeqn(start_idx:end_idx), zdeqndh(start_idx:end_idx), zeqn_absmin(start_idx:end_idx)
REAL(wp)  ::  aphscale(start_idx:end_idx)
LOGICAL   ::  l_exitnow(start_idx:end_idx)
REAL(wp), PARAMETER :: pz_exp_threshold = 1.0_wp

INTEGER   ::  niter_atgen
INTEGER   ::  j
LOGICAL   ::  imask(start_idx:end_idx)
INTEGER   ::  numActiveCells


DO j = start_idx, end_idx
    IF(vmask(j)) THEN

        aphscale(j) = 1._wp + p_so4tot(j)/Ks(j)

        IF(PRESENT(p_hini)) THEN
            zh_ini(j) = p_hini(j)
        ELSE
            CALL ahini_for_at(p_alktot(j), p_dictot(j), p_bortot(j), K1(j), K2(j), Kb(j), zh_ini(j))
        ENDIF

        CALL anw_infsup(p_dictot(j), p_bortot(j),                                           &
                        p_po4tot(j), p_siltot(j),                                           &
                        p_so4tot(j), p_flutot(j),                                           &
                        zalknw_inf(j), zalknw_sup(j))

        zdelta(j) = (p_alktot(j)-zalknw_inf(j))**2 + 4._wp*Kw(j)/aphscale(j)

        IF(p_alktot(j) >= zalknw_inf(j)) THEN
            zh_min(j) = 2._wp*Kw(j) /( p_alktot(j)-zalknw_inf(j) + SQRT(zdelta(j)) )
        ELSE
            zh_min(j) = aphscale(j)*(-(p_alktot(j)-zalknw_inf(j)) + SQRT(zdelta(j)) ) / 2._wp
        ENDIF

        zdelta(j) = (p_alktot(j)-zalknw_sup(j))**2 + 4._wp*Kw(j)/aphscale(j)

        IF(p_alktot(j) <= zalknw_sup(j)) THEN
            zh_max(j) = aphscale(j)*(-(p_alktot(j)-zalknw_sup(j)) + SQRT(zdelta(j)) ) / 2._wp
        ELSE
            zh_max(j) = 2._wp*Kw(j) /( p_alktot(j)-zalknw_sup(j) + SQRT(zdelta(j)) )
        ENDIF

        zh(j) = MAX(MIN(zh_max(j), zh_ini(j)), zh_min(j))

        zeqn_absmin(j)        = HUGE(1._wp)

    END IF ! vmask
END DO ! j

numActiveCells     = 0
niter_atgen        = 0                 ! Reset counters of iterations
DO j = start_idx, end_idx
    imask(j) = vmask(j)
END DO

DO ! Loop over iterations

    IF(niter_atgen >= jp_maxniter_atgen) THEN
        DO j = start_idx, end_idx
            IF(imask(j)) THEN
                zh(j) = -1._wp
            END IF
        END DO
        EXIT
    ENDIF

    ! Now determine the next iterate zh
    niter_atgen = niter_atgen + 1

    DO j = start_idx, end_idx
    IF(imask(j)) THEN

        zh_prev(j) = zh(j)
        zeqn(j) = equation_at(p_alktot(j), zh(j), p_dictot(j), p_bortot(j),                              &
                              p_po4tot(j), p_siltot(j),                                                  &
                              p_so4tot(j), p_flutot(j),                                                  &
                              K1(j), K2(j), Kb(j), Kw(j), Ks(j), Kf(j), K1p(j), K2p(j), K3p(j), Ksi(j),  &
                              P_DERIVEQN = zdeqndh(j))

        ! Adapt bracketing interval
        IF(zeqn(j) > 0._wp) THEN
            zh_min(j) = zh_prev(j)
        ELSEIF(zeqn(j) < 0._wp) THEN
            zh_max(j) = zh_prev(j)
        ELSE
            ! zh is the root; unlikely but, one never knows
            imask(j) = .FALSE.
            CYCLE
        ENDIF

        IF(ABS(zeqn(j)) >= 0.5_wp*zeqn_absmin(j)) THEN
            zh(j) = SQRT(zh_max(j) * zh_min(j))
            zh_lnfactor(j) = (zh(j) - zh_prev(j))/zh_prev(j) ! Required to test convergence below
        ELSE
            zh_lnfactor(j) = -zeqn(j)/(zdeqndh(j)*zh_prev(j))

            IF(ABS(zh_lnfactor(j)) > pz_exp_threshold) THEN
                zh(j)          = zh_prev(j)*EXP(zh_lnfactor(j))
            ELSE
                zh_delta(j)    = zh_lnfactor(j)*zh_prev(j)
                zh(j)          = zh_prev(j) + zh_delta(j)
            ENDIF

            IF( zh(j) < zh_min(j) ) THEN
                zh(j)                = SQRT(zh_prev(j) * zh_min(j))
                zh_lnfactor(j)       = (zh(j) - zh_prev(j))/zh_prev(j) ! Required to test convergence below
            ENDIF

            IF( zh(j) > zh_max(j) ) THEN
                zh(j)                = SQRT(zh_prev(j) * zh_max(j))
                zh_lnfactor(j)       = (zh(j) - zh_prev(j))/zh_prev(j) ! Required to test convergence below
            ENDIF
        ENDIF

        zeqn_absmin(j) = MIN( ABS(zeqn(j)), zeqn_absmin(j))

        l_exitnow(j) = (ABS(zh_lnfactor(j)) < pp_rdel_ah_target)

        IF(l_exitnow(j)) THEN
            imask(j) = .FALSE.
        END IF
    END IF
    END DO ! j

    numActiveCells = 0
    DO j = start_idx, end_idx
        IF(imask(j)) THEN
            numActiveCells = numActiveCells + 1
        END IF
    END DO

    IF(numActiveCells == 0) THEN
        EXIT
    END IF
END DO ! iterative process

DO j = start_idx, end_idx
    IF(vmask(j)) THEN

        resSAG(j) = zh(j)

        IF(PRESENT(p_val)) THEN
            IF(zh(j) > 0._wp) THEN
                p_val(j) = equation_at(p_alktot(j), zh(j), p_dictot(j), p_bortot(j),                           &
                                       p_po4tot(j), p_siltot(j),                                               &
                                       p_so4tot(j), p_flutot(j),                                               &
                                       K1(j), K2(j), Kb(j), Kw(j), Ks(j), Kf(j), K1p(j), K2p(j), K3p(j), Ksi(j))
            ELSE
                p_val(j) = HUGE(1._wp)
            ENDIF
        ENDIF
    END IF
END DO

END SUBROUTINE solve_at_general_sub

END MODULE
