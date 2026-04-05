!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Set disperson types for DFT calculations
!> \author JGH (04.2014)
! **************************************************************************************************
MODULE qs_dispersion_utils

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE input_constants,                 ONLY: &
        vdw_nl_DRSLL, vdw_nl_LMKLL, vdw_nl_RVV10, vdw_pairpot_dftd2, vdw_pairpot_dftd3, &
        vdw_pairpot_dftd3bj, vdw_pairpot_dftd4, xc_vdw_fun_nonloc, xc_vdw_fun_pairpot
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE physcon,                         ONLY: bohr,&
                                              kjmol
   USE qs_dispersion_types,             ONLY: qs_atom_dispersion_type,&
                                              qs_dispersion_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE string_utilities,                ONLY: uppercase
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dispersion_utils'

   PUBLIC :: qs_dispersion_env_set, qs_write_dispersion
   PUBLIC :: qs_scaling_dftd4
   PUBLIC :: cellhash

! **************************************************************************************************
CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param dispersion_env ...
!> \param xc_section ...
! **************************************************************************************************
   SUBROUTINE qs_dispersion_env_set(dispersion_env, xc_section)
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      TYPE(section_vals_type), POINTER                   :: xc_section

      LOGICAL                                            :: exfun, explicit
      REAL(dp), POINTER                                  :: params(:), scal(:)
      TYPE(section_vals_type), POINTER                   :: nl_section, pp_section, vdw_section, &
                                                            xc_fun_section

      CPASSERT(ASSOCIATED(dispersion_env))

      ! set general defaults
      dispersion_env%doabc = .FALSE.
      dispersion_env%c9cnst = .FALSE.
      dispersion_env%lrc = .FALSE.
      dispersion_env%srb = .FALSE.
      dispersion_env%verbose = .FALSE.
      dispersion_env%nd3_exclude_pair = 0
      NULLIFY (dispersion_env%c6ab, dispersion_env%maxci, dispersion_env%r0ab, dispersion_env%rcov, &
               dispersion_env%r2r4, dispersion_env%cn, dispersion_env%cnkind, dispersion_env%cnlist, &
               dispersion_env%d3_exclude_pair)
      NULLIFY (dispersion_env%q_mesh, dispersion_env%kernel, dispersion_env%d2phi_dk2, &
               dispersion_env%d2y_dx2)
      NULLIFY (dispersion_env%sab_vdw, dispersion_env%sab_cn)
      NULLIFY (dispersion_env%dftd_section)
      NULLIFY (vdw_section, xc_fun_section)
      vdw_section => section_vals_get_subs_vals(xc_section, "vdw_potential")
      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      CALL section_vals_val_get(vdw_section, "POTENTIAL_TYPE", i_val=dispersion_env%type)
      IF (dispersion_env%type == xc_vdw_fun_pairpot) THEN
         NULLIFY (pp_section)
         pp_section => section_vals_get_subs_vals(vdw_section, "PAIR_POTENTIAL")
         CALL section_vals_val_get(pp_section, "VERBOSE_OUTPUT", l_val=dispersion_env%verbose)
         CALL section_vals_val_get(pp_section, "TYPE", i_val=dispersion_env%pp_type)
         IF (dispersion_env%pp_type == vdw_pairpot_dftd2) THEN
            ! functional parameters for Grimme D2 type
            CALL section_vals_val_get(pp_section, "EXP_PRE", r_val=dispersion_env%exp_pre)
            CALL section_vals_val_get(pp_section, "SCALING", explicit=explicit)
            IF (.NOT. explicit) THEN
               CALL section_vals_val_get(pp_section, "REFERENCE_FUNCTIONAL", explicit=exfun)
               IF (.NOT. exfun) THEN
                  CALL cp_abort(__LOCATION__, "D2 vdW REFERENCE_FUNCTIONAL expected a second parameter "// &
                                "(example REFERENCE_FUNCTIONAL PBE). "// &
                                "Go to https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3/ "// &
                                "for a full list of supported functionals")
               END IF
               CALL qs_scaling_dftd2(dispersion_env%scaling, vdw_section)
            ELSE
               CALL section_vals_val_get(pp_section, "SCALING", r_val=dispersion_env%scaling)
            END IF
         ELSE
            dispersion_env%exp_pre = 0._dp
            dispersion_env%scaling = 0._dp
         END IF
         IF (dispersion_env%pp_type == vdw_pairpot_dftd3 .OR. &
             dispersion_env%pp_type == vdw_pairpot_dftd3bj) THEN
            ! functional parameters for Grimme DFT-D3 type
            CALL section_vals_val_get(pp_section, "EPS_CN", r_val=dispersion_env%eps_cn)
            CALL section_vals_val_get(pp_section, "CALCULATE_C9_TERM", l_val=dispersion_env%doabc)
            CALL section_vals_val_get(pp_section, "REFERENCE_C9_TERM", l_val=dispersion_env%c9cnst)
            CALL section_vals_val_get(pp_section, "LONG_RANGE_CORRECTION", l_val=dispersion_env%lrc)
            CALL section_vals_val_get(pp_section, "SHORT_RANGE_CORRECTION", l_val=dispersion_env%srb)
            CALL section_vals_val_get(pp_section, "SHORT_RANGE_CORRECTION_PARAMETERS", r_vals=params)
            dispersion_env%srb_params(1:4) = params(1:4)
            ! KG corrections
            CALL section_vals_val_get(pp_section, "MOLECULE_CORRECTION", l_val=dispersion_env%domol)
            CALL section_vals_val_get(pp_section, "MOLECULE_CORRECTION_C8", r_val=dispersion_env%kgc8)
            IF (dispersion_env%pp_type == vdw_pairpot_dftd3) THEN
               CALL section_vals_val_get(pp_section, "D3_SCALING", explicit=explicit)
            ELSE IF (dispersion_env%pp_type == vdw_pairpot_dftd3bj) THEN
               CALL section_vals_val_get(pp_section, "D3BJ_SCALING", explicit=explicit)
            END IF
            IF (.NOT. explicit) THEN
               CALL section_vals_val_get(pp_section, "REFERENCE_FUNCTIONAL", explicit=exfun)
               IF (.NOT. exfun) THEN
                  CALL cp_abort(__LOCATION__, "D3 vdW REFERENCE_FUNCTIONAL expected a second parameter "// &
                                "(example REFERENCE_FUNCTIONAL PBE). "// &
                                "Go to https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3/ for a full "// &
                                "list of supported functionals")
               ELSE
                  CALL section_vals_val_get(vdw_section, &
                                            "PAIR_POTENTIAL%REFERENCE_FUNCTIONAL", &
                                            c_val=dispersion_env%ref_functional)
               END IF
               IF (dispersion_env%pp_type == vdw_pairpot_dftd3) THEN
                  CALL qs_scaling_dftd3(dispersion_env%s6, dispersion_env%sr6, &
                                        dispersion_env%s8, vdw_section)
               ELSE IF (dispersion_env%pp_type == vdw_pairpot_dftd3bj) THEN
                  CALL qs_scaling_dftd3bj(dispersion_env%s6, dispersion_env%a1, dispersion_env%s8, &
                                          dispersion_env%a2, vdw_section)
               END IF
            ELSE
               IF (dispersion_env%pp_type == vdw_pairpot_dftd3) THEN
                  ! zero damping
                  CALL section_vals_val_get(pp_section, "D3_SCALING", r_vals=scal)
                  dispersion_env%s6 = scal(1)
                  dispersion_env%sr6 = scal(2)
                  dispersion_env%s8 = scal(3)
                  dispersion_env%a1 = 0.0_dp
                  dispersion_env%a2 = 0.0_dp
               ELSE IF (dispersion_env%pp_type == vdw_pairpot_dftd3bj) THEN
                  ! BJ damping
                  CALL section_vals_val_get(pp_section, "D3BJ_SCALING", r_vals=scal)
                  dispersion_env%s6 = scal(1)
                  dispersion_env%a1 = scal(2)
                  dispersion_env%s8 = scal(3)
                  dispersion_env%a2 = scal(4)
                  dispersion_env%sr6 = 0.0_dp
               END IF
            END IF
         ELSE
            dispersion_env%s6 = 0._dp
            dispersion_env%sr6 = 0._dp
            dispersion_env%s8 = 0._dp
            dispersion_env%s9 = 0._dp
            dispersion_env%a1 = 0._dp
            dispersion_env%a2 = 0._dp
            dispersion_env%eps_cn = 0._dp
         END IF
         IF (dispersion_env%pp_type == vdw_pairpot_dftd4) THEN
            CALL section_vals_val_get(pp_section, "D4_SCALING", explicit=explicit)
            IF (.NOT. explicit) THEN
               CALL section_vals_val_get(pp_section, "REFERENCE_FUNCTIONAL", explicit=exfun)
               IF (.NOT. exfun) THEN
                  CPABORT("Missing REFERENCE_FUNCTIONAL or D4_SCALING for D4")
               ELSE
                  CALL qs_scaling_dftd4(dispersion_env%s6, dispersion_env%a1, &
                                        dispersion_env%s8, dispersion_env%a2, vdw_section)
                  dispersion_env%sr6 = 0.0_dp
                  dispersion_env%ref_functional = "none"
               END IF
            ELSE
               CALL section_vals_val_get(pp_section, "D4_SCALING", r_vals=scal)
               dispersion_env%s6 = scal(1)
               dispersion_env%a1 = scal(2)
               dispersion_env%s8 = scal(3)
               dispersion_env%a2 = scal(4)
               dispersion_env%sr6 = 0.0_dp
               dispersion_env%ref_functional = "none"
            END IF
            IF (TRIM(ADJUSTL(dispersion_env%ref_functional)) == "none") THEN
               CPABORT("Missing REFERENCE_FUNCTIONAL or D4_SCALING for D4")
            END IF
            CALL section_vals_val_get(pp_section, "EPS_CN", r_val=dispersion_env%eps_cn)
            CALL section_vals_val_get(pp_section, "D4_REFERENCE_CODE", &
                                      l_val=dispersion_env%d4_reference_code)
            CALL section_vals_val_get(pp_section, "D4_DEBUG", l_val=dispersion_env%d4_debug)
            CALL section_vals_val_get(pp_section, "D4_CUTOFF", r_val=dispersion_env%rc_d4)
            CALL section_vals_val_get(pp_section, "D4_CN_CUTOFF", r_val=dispersion_env%rc_cn)
            CALL section_vals_val_get(pp_section, "FACTOR_S9_TERM", r_val=dispersion_env%s9)
            !C9 term default=T for D4
            CALL section_vals_val_get(pp_section, "CALCULATE_C9_TERM", explicit=exfun)
            IF (exfun) THEN
               CALL section_vals_val_get(pp_section, "CALCULATE_C9_TERM", l_val=dispersion_env%doabc)
            ELSE
               dispersion_env%doabc = .TRUE.
            END IF
         END IF
         CALL section_vals_val_get(pp_section, "R_CUTOFF", r_val=dispersion_env%rc_disp)
         CALL section_vals_val_get(pp_section, "PARAMETER_FILE_NAME", &
                                   c_val=dispersion_env%parameter_file_name)
         ! set DFTD section for output handling
         dispersion_env%dftd_section => pp_section
      ELSE IF (dispersion_env%type == xc_vdw_fun_nonloc) THEN
         NULLIFY (nl_section)
         nl_section => section_vals_get_subs_vals(vdw_section, "NON_LOCAL")
         CALL section_vals_val_get(nl_section, "VERBOSE_OUTPUT", l_val=dispersion_env%verbose)
         CALL section_vals_val_get(nl_section, "KERNEL_FILE_NAME", &
                                   c_val=dispersion_env%kernel_file_name)
         CALL section_vals_val_get(nl_section, "TYPE", i_val=dispersion_env%nl_type)
         CALL section_vals_val_get(nl_section, "CUTOFF", r_val=dispersion_env%pw_cutoff)
         CALL section_vals_val_get(nl_section, "PARAMETERS", r_vals=params)
         CALL section_vals_val_get(nl_section, "SCALE", r_val=dispersion_env%scale_rvv10)
         dispersion_env%b_value = params(1)
         dispersion_env%c_value = params(2)
      END IF
   END SUBROUTINE qs_dispersion_env_set

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dispersion_env ...
!> \param ounit ...
! **************************************************************************************************
   SUBROUTINE qs_write_dispersion(qs_env, dispersion_env, ounit)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      INTEGER, INTENT(in), OPTIONAL                      :: ounit

      CHARACTER(LEN=2)                                   :: symbol
      INTEGER                                            :: i, ikind, nkind, output_unit
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_atom_dispersion_type), POINTER             :: disp
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: dft_section

      IF (PRESENT(ounit)) THEN
         output_unit = ounit
      ELSE
         NULLIFY (logger)
         logger => cp_get_default_logger()

         dft_section => section_vals_get_subs_vals(qs_env%input, "DFT")
         output_unit = cp_print_key_unit_nr(logger, dft_section, &
                                            "PRINT%DFT_CONTROL_PARAMETERS", extension=".Log")
      END IF

      IF (output_unit > 0) THEN
         ! vdW type specific output
         IF (dispersion_env%type == xc_vdw_fun_pairpot) THEN
            WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T67,'Pair Potential')")
            ! Pair potentials
            IF (dispersion_env%pp_type == vdw_pairpot_dftd2) THEN
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T35,'DFT-D2')")
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T35,'Potential Form: S. Grimme, JCC 27: 1787 (2006)')")
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T35,'Cutoff Radius [Bohr]:',T73,F8.2)") dispersion_env%rc_disp
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T35,'Scaling Factor:',T73,F8.4)") dispersion_env%scaling
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T35,'Exp Prefactor for Damping:',T73,F8.1)") dispersion_env%exp_pre
               CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)
               nkind = SIZE(atomic_kind_set)
               DO ikind = 1, nkind
                  CALL get_atomic_kind(atomic_kind_set(ikind), element_symbol=symbol)
                  CALL get_qs_kind(qs_kind_set(ikind), dispersion=disp)
                  IF (disp%defined) THEN
                     WRITE (output_unit, fmt="(' vdW PARAMETER| ',T18,'Atom=',A2, "// &
                            "T28,'C6[J*nm^6*mol^-1]=',F8.4,T63,'r(vdW)[A]=',F8.4)") &
                        symbol, disp%c6/(1000._dp*bohr**6/kjmol), disp%vdw_radii/bohr
                  ELSE
                     WRITE (output_unit, fmt="(' vdW PARAMETER| ',T20,'Atom=',A2,T70,'not defined')")
                  END IF
               END DO
            ELSE IF (dispersion_env%pp_type == vdw_pairpot_dftd3) THEN
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'DFT-D3 (Version 3.1)')")
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'Potential Form: S. Grimme et al, JCP 132: 154104 (2010)')")
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'Zero Damping')")
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'Cutoff Radius [Bohr]:',T73,F8.2)") dispersion_env%rc_disp
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'s6 Scaling Factor:',T73,F8.4)") dispersion_env%s6
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'sr6 Scaling Factor:',T73,F8.4)") dispersion_env%sr6
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'s8 Scaling Factor:',T73,F8.4)") dispersion_env%s8
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'Cutoff for CN calculation:',T69,E12.4)") dispersion_env%eps_cn
               IF (dispersion_env%nd3_exclude_pair > 0) THEN
                  DO i = 1, dispersion_env%nd3_exclude_pair
                     WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'Excluded Pairs: ',T76,I2,' ',I2)") &
                        dispersion_env%d3_exclude_pair(i, :)
                  END DO
               END IF
            ELSE IF (dispersion_env%pp_type == vdw_pairpot_dftd3bj) THEN
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'DFT-D3 (Version 3.1)')")
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'Potential Form: S. Grimme et al, JCP 132: 154104 (2010)')")
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'BJ Damping: S. Grimme et al, JCC 32: 1456 (2011)')")
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'Cutoff Radius [Bohr]:',T73,F8.2)") dispersion_env%rc_disp
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'s6 Scaling Factor:',T73,F8.4)") dispersion_env%s6
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'a1 Damping Factor:',T73,F8.4)") dispersion_env%a1
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'s8 Scaling Factor:',T73,F8.4)") dispersion_env%s8
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'a2 Damping Factor:',T73,F8.4)") dispersion_env%a2
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'Cutoff for CN calculation:',T69,E12.4)") dispersion_env%eps_cn
               IF (dispersion_env%nd3_exclude_pair > 0) THEN
                  DO i = 1, dispersion_env%nd3_exclude_pair
                     WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'Excluded Kind Pairs: ',T76,I2,' ',I2)") &
                        dispersion_env%d3_exclude_pair(i, :)
                  END DO
               END IF
            ELSE IF (dispersion_env%pp_type == vdw_pairpot_dftd4) THEN
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'DFT-D4(Version 4.1.0)')")
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'see https://github.com/dftd4/dftd4')")
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'E. Caldeweyher et al, PCCP 22: 8499 (2020)')")
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'E. Caldeweyher et al, JCP 150: 154122 (2019)')")
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'E. Caldeweyher et al, JCP 147: 034112 (2017)')")
            END IF
         ELSE IF (dispersion_env%type == xc_vdw_fun_nonloc) THEN
            WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T61,'Non-local Functional')")
            WRITE (output_unit, &
                   fmt="(' vdW POTENTIAL| ','Implementation: G. Roman-Perez, J. Soler, PRL 103: 096102 (2009)')")
            WRITE (output_unit, &
                   fmt="(' vdW POTENTIAL| ',T38,' T. Thonhauser et al, PRB 76: 125112 (2007)')")
            WRITE (output_unit, &
                   fmt="(' vdW POTENTIAL| ',T22,' R. Sabatini et al, J.Phys:Condens Matter 24: 424209 (2012)')")
            WRITE (output_unit, &
                   fmt="(' vdW POTENTIAL| ',T16,' Based on QE implementation by Brian Kolb, Timo Thonhauser (2009)')")
            SELECT CASE (dispersion_env%nl_type)
            CASE DEFAULT
               ! unknown functional
               CPABORT("")
            CASE (vdw_nl_DRSLL)
               WRITE (output_unit, &
                      fmt="(' vdW POTENTIAL| ','DRSLL Functional:           M. Dion et al, PRL 92: 246401 (2004)')")
            CASE (vdw_nl_LMKLL)
               WRITE (output_unit, &
                      fmt="(' vdW POTENTIAL| ','LMKLL Functional:            K. Lee et al, PRB 82: 081101 (2010)')")
            CASE (vdw_nl_RVV10)
               WRITE (output_unit, &
                      fmt="(' vdW POTENTIAL| ','RVV10 Functional:    R. Sabatini et al, PRB 87: 041108(R) (2013)')")
            END SELECT
            IF (dispersion_env%verbose) THEN
               WRITE (output_unit, &
                      fmt="(' vdW POTENTIAL| ','         Carrying out vdW-DF run using the following parameters:')")
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ','Nqs =',I8,'        Nr_points =',I8,'       r_max =',F10.3)") &
                  dispersion_env%nqs, dispersion_env%nr_points, dispersion_env%r_max
               WRITE (output_unit, fmt="(' vdW POTENTIAL| ','q_mesh =')")
               WRITE (output_unit, fmt="(8X,4F18.8)") (dispersion_env%q_mesh(i), i=1, dispersion_env%nqs)
               WRITE (output_unit, &
                      fmt="(' vdW POTENTIAL| ','Density cutoff for convolution [a.u.]:',T71,F10.1)") &
                  dispersion_env%pw_cutoff
            END IF
         END IF
      END IF
      IF (.NOT. PRESENT(ounit)) THEN
         CALL cp_print_key_finished_output(output_unit, logger, dft_section, &
                                           "PRINT%DFT_CONTROL_PARAMETERS")
      END IF

   END SUBROUTINE qs_write_dispersion

! **************************************************************************************************
!> \brief ...
!> \param scaling ...
!> \param vdw_section ...
! **************************************************************************************************
   SUBROUTINE qs_scaling_dftd2(scaling, vdw_section)
      REAL(KIND=dp), INTENT(inout)                       :: scaling
      TYPE(section_vals_type), POINTER                   :: vdw_section

      CHARACTER(LEN=default_string_length)               :: functional

      CALL section_vals_val_get(vdw_section, "PAIR_POTENTIAL%REFERENCE_FUNCTIONAL", c_val=functional)

      SELECT CASE (TRIM(functional))
      CASE DEFAULT
         ! unknown functional
         CPABORT("No DFT-D2 s6 value available for this functional:"//TRIM(functional))
      CASE ("BLYP")
         scaling = 1.20_dp
      CASE ("B3LYP")
         scaling = 1.05_dp
      CASE ("TPSS")
         scaling = 1.00_dp
      CASE ("PBE")
         scaling = 0.75_dp
      CASE ("PBE0")
         scaling = 0.6_dp
      CASE ("B2PLYP")
         scaling = 0.55_dp
      CASE ("BP86")
         scaling = 1.05_dp
      CASE ("B97")
         scaling = 1.25_dp
      END SELECT

   END SUBROUTINE qs_scaling_dftd2

! **************************************************************************************************
!> \brief ...
!> \param s6 ...
!> \param sr6 ...
!> \param s8 ...
!> \param vdw_section ...
! **************************************************************************************************
   SUBROUTINE qs_scaling_dftd3(s6, sr6, s8, vdw_section)

      REAL(KIND=dp), INTENT(inout)                       :: s6, sr6, s8
      TYPE(section_vals_type), POINTER                   :: vdw_section

      CHARACTER(LEN=default_string_length)               :: functional

      CALL section_vals_val_get(vdw_section, "PAIR_POTENTIAL%REFERENCE_FUNCTIONAL", c_val=functional)
      CALL uppercase(functional)
      ! values for different functionals from:
      ! https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3
      ! L. Goerigk et al. PCCP 2017, 32147-32744, SI
      SELECT CASE (TRIM(functional))
      CASE DEFAULT
         ! unknown functional
         CPABORT("No DFT-D3 values available for this functional:"//TRIM(functional))
      CASE ("B1B95")
         s6 = 1.000_dp
         sr6 = 1.613_dp
         s8 = 1.868_dp
      CASE ("B2GPPLYP")
         ! L. Goerigk and S. Grimme
         ! J. Chem. Theory Comput. 2011, 7, 291-309; doi:10.1021/ct100466k
         s6 = 0.56_dp
         sr6 = 1.586_dp
         s8 = 0.760_dp
      CASE ("B2PLYP")
         ! L. Goerigk and S. Grimme
         ! J. Chem. Theory Comput. 2011, 7, 291-309; doi:10.1021/ct100466k
         s6 = 0.64_dp
         sr6 = 1.427_dp
         s8 = 1.022_dp
      CASE ("DSD-BLYP")
         ! L. Goerigk and S. Grimme
         ! J. Chem. Theory Comput. 2011, 7, 291-309; doi:10.1021/ct100466k
         s6 = 0.50_dp
         sr6 = 1.569_dp
         s8 = 0.705_dp
      CASE ("B3LYP")
         ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
         ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
         s6 = 1.000_dp
         sr6 = 1.261_dp
         s8 = 1.703_dp
      CASE ("B97-D")
         ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
         ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
         s6 = 1.000_dp
         sr6 = 0.892_dp
         s8 = 0.909_dp
      CASE ("BHLYP")
         s6 = 1.000_dp
         sr6 = 1.370_dp
         s8 = 1.442_dp
      CASE ("BLYP")
         ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
         ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
         s6 = 1.000_dp
         sr6 = 1.094_dp
         s8 = 1.682_dp
      CASE ("BP86")
         ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
         ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
         s6 = 1.000_dp
         sr6 = 1.139_dp
         s8 = 1.683_dp
      CASE ("BPBE")
         s6 = 1.000_dp
         sr6 = 1.087_dp
         s8 = 2.033_dp
      CASE ("MPWLYP")
         s6 = 1.000_dp
         sr6 = 1.239_dp
         s8 = 1.098_dp
      CASE ("PBE")
         ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
         ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
         s6 = 1.000_dp
         sr6 = 1.217_dp
         s8 = 0.722_dp
      CASE ("PBEHPBE")
         s6 = 1.000_dp
         sr6 = 1.5703_dp
         s8 = 1.4010_dp
      CASE ("PBE0")
         ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
         ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
         s6 = 1.000_dp
         sr6 = 1.287_dp
         s8 = 0.928_dp
      CASE ("PW6B95")
         ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
         ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
         s6 = 1.000_dp
         sr6 = 1.532_dp
         s8 = 0.862_dp
      CASE ("PWB6K")
         s6 = 1.000_dp
         sr6 = 1.660_dp
         s8 = 0.550_dp
      CASE ("REVPBE")
         ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
         ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
         s6 = 1.000_dp
         sr6 = 0.923_dp
         s8 = 1.010_dp
      CASE ("RPBE")
         s6 = 1.000_dp
         sr6 = 0.872_dp
         s8 = 0.514_dp
      CASE ("TPSS")
         ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
         ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
         s6 = 1.000_dp
         sr6 = 1.166_dp
         s8 = 1.105_dp
      CASE ("TPSS0")
         ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
         ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
         s6 = 1.000_dp
         sr6 = 1.252_dp
         s8 = 1.242_dp
      CASE ("TPSSH")
         s6 = 1.000_dp
         sr6 = 1.223_dp
         s8 = 1.219_dp
      CASE ("B1LYP")
         s6 = 1.000_dp
         sr6 = 1.3725_dp
         s8 = 1.9467_dp
      CASE ("B1P86")
         s6 = 1.000_dp
         sr6 = 1.1815_dp
         s8 = 1.1209_dp
      CASE ("B3P86")
         s6 = 1.000_dp
         sr6 = 1.1897_dp
         s8 = 1.1961_dp
      CASE ("B3PW91")
         s6 = 1.000_dp
         sr6 = 1.176_dp
         s8 = 1.775_dp
      CASE ("BMK")
         s6 = 1.000_dp
         sr6 = 1.931_dp
         s8 = 2.168_dp
      CASE ("CAMB3LYP")
         s6 = 1.000_dp
         sr6 = 1.378_dp
         s8 = 1.217_dp
      CASE ("LCWPBE")
         s6 = 1.000_dp
         sr6 = 1.355_dp
         s8 = 1.279_dp
      CASE ("M052X")
         s6 = 1.000_dp
         sr6 = 1.417_dp
         s8 = 0.000_dp
      CASE ("M05")
         s6 = 1.000_dp
         sr6 = 1.373_dp
         s8 = 0.595_dp
      CASE ("M062X")
         s6 = 1.000_dp
         sr6 = 1.619_dp
         s8 = 0.000_dp
      CASE ("M06HF")
         s6 = 1.000_dp
         sr6 = 1.446_dp
         s8 = 0.000_dp
      CASE ("M06L")
         s6 = 1.000_dp
         sr6 = 1.581_dp
         s8 = 0.000_dp
      CASE ("M06N")
         s6 = 1.000_dp
         sr6 = 1.325_dp
         s8 = 0.000_dp
      CASE ("HCTH120")
         s6 = 1.000_dp
         sr6 = 1.221_dp
         s8 = 1.206_dp
      CASE ("HCTH407")
         s6 = 1.000_dp
         sr6 = 4.0426_dp
         s8 = 2.7694_dp
      CASE ("MPW2PLYP")
         s6 = 1.000_dp
         sr6 = 1.5527_dp
         s8 = 0.7529_dp
      CASE ("PKZB")
         s6 = 1.000_dp
         sr6 = 0.6327_dp
         s8 = 0.000_dp
      CASE ("PTPSS")
         s6 = 0.750_dp
         sr6 = 1.541_dp
         s8 = 0.879_dp
      CASE ("PWPB95")
         s6 = 0.820_dp
         sr6 = 1.557_dp
         s8 = 0.705_dp
      CASE ("OLYP")
         s6 = 1.000_dp
         sr6 = 0.806_dp
         s8 = 1.764_dp
      CASE ("OPBE")
         s6 = 1.000_dp
         sr6 = 0.837_dp
         s8 = 2.055_dp
      CASE ("OTPSS")
         s6 = 1.000_dp
         sr6 = 1.128_dp
         s8 = 1.494_dp
      CASE ("PBE1KCIS")
         s6 = 1.000_dp
         sr6 = 3.6355_dp
         s8 = 1.7934_dp
      CASE ("PBE38")
         s6 = 1.000_dp
         sr6 = 1.333_dp
         s8 = 0.998_dp
      CASE ("PBEH1PBE")
         s6 = 1.000_dp
         sr6 = 1.3719_dp
         s8 = 1.0430_dp
      CASE ("PBESOL")
         s6 = 1.000_dp
         sr6 = 1.345_dp
         s8 = 0.612_dp
      CASE ("REVSSB")
         s6 = 1.000_dp
         sr6 = 1.221_dp
         s8 = 0.560_dp
      CASE ("REVTPSS")
         s6 = 1.000_dp
         sr6 = 1.3491_dp
         s8 = 1.3666_dp
      CASE ("SSB")
         s6 = 1.000_dp
         sr6 = 1.215_dp
         s8 = 0.663_dp
      CASE ("B97-1")
         s6 = 1.000_dp
         sr6 = 3.7924_dp
         s8 = 1.6418_dp
      CASE ("B97-2")
         s6 = 1.000_dp
         sr6 = 1.7066_dp
         s8 = 2.4661_dp
      CASE ("B98")
         s6 = 1.000_dp
         sr6 = 2.6895_dp
         s8 = 1.9078_dp
      CASE ("BOP")
         s6 = 1.000_dp
         sr6 = 0.929_dp
         s8 = 1.975_dp
      CASE ("HISS")
         s6 = 1.000_dp
         sr6 = 1.3338_dp
         s8 = 0.7615_dp
      CASE ("HSE03")
         s6 = 1.000_dp
         sr6 = 1.3944_dp
         s8 = 1.0156_dp
      CASE ("HSE06")
         s6 = 1.000_dp
         sr6 = 1.129_dp
         s8 = 0.109_dp
      CASE ("M08HX")
         s6 = 1.000_dp
         sr6 = 1.6247_dp
         s8 = 0.000_dp
      CASE ("MN15L")
         s6 = 1.000_dp
         sr6 = 3.3388_dp
         s8 = 0.000_dp
      CASE ("MPWPW91")
         s6 = 1.0000_dp
         sr6 = 1.3725_dp
         s8 = 1.9467_dp
      CASE ("MPW1B95")
         s6 = 1.000_dp
         sr6 = 1.605_dp
         s8 = 1.118_dp
      CASE ("MPW1KCIS")
         s6 = 1.000_dp
         sr6 = 1.7231_dp
         s8 = 2.2917_dp
      CASE ("MPW1LYP")
         s6 = 1.000_dp
         sr6 = 2.0512_dp
         s8 = 1.9529_dp
      CASE ("MPW1PW91")
         s6 = 1.000_dp
         sr6 = 1.2892_dp
         s8 = 1.4758_dp
      CASE ("MPWB1K")
         s6 = 1.000_dp
         sr6 = 1.671_dp
         s8 = 1.061_dp
      CASE ("MPWKCIS1K")
         s6 = 1.000_dp
         sr6 = 1.4853_dp
         s8 = 1.7553_dp
      CASE ("O3LYP")
         s6 = 1.000_dp
         sr6 = 1.4060_dp
         s8 = 1.8058_dp
      CASE ("PW1PW")
         s6 = 1.000_dp
         sr6 = 1.4968_dp
         s8 = 1.1786_dp
      CASE ("PW91P86")
         s6 = 1.0000_dp
         sr6 = 2.1040_dp
         s8 = 0.8747_dp
      CASE ("REVPBE0")
         s6 = 1.000_dp
         sr6 = 0.949_dp
         s8 = 0.792_dp
      CASE ("REVPBE38")
         s6 = 1.000_dp
         sr6 = 1.021_dp
         s8 = 0.862_dp
      CASE ("REVTPSSh")
         s6 = 1.000_dp
         sr6 = 1.3224_dp
         s8 = 1.2504_dp
      CASE ("REVTPSS0")
         s6 = 1.000_dp
         sr6 = 1.2881_dp
         s8 = 1.0649_dp
      CASE ("TPSS1KCIS")
         s6 = 1.000_dp
         sr6 = 1.7729_dp
         s8 = 2.0902_dp
      CASE ("THCTHHYB")
         s6 = 1.000_dp
         sr6 = 1.5001_dp
         s8 = 1.6302_dp
      CASE ("RPW86PBE")
         s6 = 1.000_dp
         sr6 = 1.224_dp
         s8 = 0.901_dp
      CASE ("SCAN")
         s6 = 1.000_dp
         sr6 = 1.324_dp
         s8 = 0.000_dp
      CASE ("THCTH")
         s6 = 1.000_dp
         sr6 = 0.932_dp
         s8 = 0.5662_dp
      CASE ("XLYP")
         s6 = 1.0000_dp
         sr6 = 0.9384_dp
         s8 = 0.7447_dp
      CASE ("X3LYP")
         s6 = 1.000_dp
         sr6 = 1.0000_dp
         s8 = 0.2990_dp
      END SELECT

   END SUBROUTINE qs_scaling_dftd3

! **************************************************************************************************
!> \brief ...
!> \param s6 ...
!> \param a1 ...
!> \param s8 ...
!> \param a2 ...
!> \param vdw_section ...
! **************************************************************************************************
   SUBROUTINE qs_scaling_dftd3bj(s6, a1, s8, a2, vdw_section)
      REAL(KIND=dp), INTENT(inout)                       :: s6, a1, s8, a2
      TYPE(section_vals_type), POINTER                   :: vdw_section

      CHARACTER(LEN=default_string_length)               :: functional

      CALL section_vals_val_get(vdw_section, "PAIR_POTENTIAL%REFERENCE_FUNCTIONAL", c_val=functional)

      ! values for different functionals from:
      ! http://www.thch.uni-bonn.de/tc/downloads/DFT-D3/functionalsbj.html
      ! L. Goerigk et al. PCCP 2017, 32147-32744, SI
      SELECT CASE (TRIM(functional))
      CASE DEFAULT
         ! unknown functional
         CPABORT("No DFT-D3(BJ) values available for this functional:"//TRIM(functional))
      CASE ("B1B95")
         s6 = 1.0000_dp
         a1 = 0.2092_dp
         s8 = 1.4507_dp
         a2 = 5.5545_dp
      CASE ("B2GPPLYP")
         s6 = 0.5600_dp
         a1 = 0.0000_dp
         s8 = 0.2597_dp
         a2 = 6.3332_dp
      CASE ("B3PW91")
         s6 = 1.0000_dp
         a1 = 0.4312_dp
         s8 = 2.8524_dp
         a2 = 4.4693_dp
      CASE ("BHLYP")
         s6 = 1.0000_dp
         a1 = 0.2793_dp
         s8 = 1.0354_dp
         a2 = 4.9615_dp
      CASE ("BMK")
         s6 = 1.0000_dp
         a1 = 0.1940_dp
         s8 = 2.0860_dp
         a2 = 5.9197_dp
      CASE ("BOP")
         s6 = 1.0000_dp
         a1 = 0.4870_dp
         s8 = 3.2950_dp
         a2 = 3.5043_dp
      CASE ("BPBE")
         s6 = 1.0000_dp
         a1 = 0.4567_dp
         s8 = 4.0728_dp
         a2 = 4.3908_dp
      CASE ("B97-3c")
         s6 = 1.0000_dp
         a1 = 0.3700_dp
         s8 = 1.5000_dp
         a2 = 4.1000_dp
      CASE ("CAMB3LYP")
         s6 = 1.0000_dp
         a1 = 0.3708_dp
         s8 = 2.0674_dp
         a2 = 5.4743_dp
      CASE ("DSDBLYP")
         s6 = 0.5000_dp
         a1 = 0.0000_dp
         s8 = 0.2130_dp
         a2 = 6.0519_dp
      CASE ("DSDPBEP86")
         s6 = 0.4180_dp
         a1 = 0.0000_dp
         s8 = 0.0000_dp
         a2 = 5.6500_dp
      CASE ("DSDPBEB95")
         s6 = 0.6100_dp
         a1 = 0.0000_dp
         s8 = 0.0000_dp
         a2 = 6.2000_dp
      CASE ("LCWPBE")
         s6 = 1.0000_dp
         a1 = 0.3919_dp
         s8 = 1.8541_dp
         a2 = 5.0897_dp
      CASE ("LCWhPBE")
         s6 = 1.0000_dp
         a1 = 0.2746_dp
         s8 = 1.1908_dp
         a2 = 5.3157_dp
      CASE ("MPW1B95")
         s6 = 1.0000_dp
         a1 = 0.1955_dp
         s8 = 1.0508_dp
         a2 = 6.4177_dp
      CASE ("MPW2PLYP")
         s6 = 0.6600_dp
         a1 = 0.4105_dp
         s8 = 0.6223_dp
         a2 = 5.0136_dp
      CASE ("MPWB1K")
         s6 = 1.0000_dp
         a1 = 0.1474_dp
         s8 = 0.9499_dp
         a2 = 6.6223_dp
      CASE ("MPWLYP")
         s6 = 1.0000_dp
         a1 = 0.4831_dp
         s8 = 2.0077_dp
         a2 = 4.5323_dp
      CASE ("OLYP")
         s6 = 1.0000_dp
         a1 = 0.5299_dp
         s8 = 2.6205_dp
         a2 = 2.8065_dp
      CASE ("OPBE")
         s6 = 1.0000_dp
         a1 = 0.5512_dp
         s8 = 3.3816_dp
         a2 = 2.9444_dp
      CASE ("OTPSS")
         s6 = 1.0000_dp
         a1 = 0.4634_dp
         s8 = 2.7495_dp
         a2 = 4.3153_dp
      CASE ("PBE38")
         s6 = 1.0000_dp
         a1 = 0.3995_dp
         s8 = 1.4623_dp
         a2 = 5.1405_dp
      CASE ("PBEsol")
         s6 = 1.0000_dp
         a1 = 0.4466_dp
         s8 = 2.9491_dp
         a2 = 6.1742_dp
      CASE ("PTPSS")
         s6 = 0.7500_dp
         a1 = 0.0000_dp
         s8 = 0.2804_dp
         a2 = 6.5745_dp
      CASE ("PWB6K")
         s6 = 1.0000_dp
         a1 = 0.1805_dp
         s8 = 0.9383_dp
         a2 = 7.7627_dp
      CASE ("revSSB")
         s6 = 1.0000_dp
         a1 = 0.4720_dp
         s8 = 0.4389_dp
         a2 = 4.0986_dp
      CASE ("SSB")
         s6 = 1.0000_dp
         a1 = -0.0952_dp
         s8 = -0.1744_dp
         a2 = 5.2170_dp
      CASE ("TPSSh")
         s6 = 1.0000_dp
         a1 = 0.4529_dp
         s8 = 2.2382_dp
         a2 = 4.6550_dp
      CASE ("HCTH120")
         s6 = 1.0000_dp
         a1 = 0.3563_dp
         s8 = 1.0821_dp
         a2 = 4.3359_dp
      CASE ("B2PLYP")
         s6 = 0.6400_dp
         a1 = 0.3065_dp
         s8 = 0.9147_dp
         a2 = 5.0570_dp
      CASE ("B1LYP")
         s6 = 1.0000_dp
         a1 = 0.1986_dp
         s8 = 2.1167_dp
         a2 = 5.3875_dp
      CASE ("B1P86")
         s6 = 1.0000_dp
         a1 = 0.4724_dp
         s8 = 3.5681_dp
         a2 = 4.9858_dp
      CASE ("B3LYP")
         s6 = 1.0000_dp
         a1 = 0.3981_dp
         s8 = 1.9889_dp
         a2 = 4.4211_dp
      CASE ("B3P86")
         s6 = 1.0000_dp
         a1 = 0.4601_dp
         s8 = 3.3211_dp
         a2 = 4.9294_dp
      CASE ("B97-1")
         s6 = 1.0000_dp
         a1 = 0.0000_dp
         s8 = 0.4814_dp
         a2 = 6.2279_dp
      CASE ("B97-2")
         s6 = 1.0000_dp
         a1 = 0.0000_dp
         s8 = 0.9448_dp
         a2 = 5.4603_dp
      CASE ("B97-D")
         s6 = 1.0000_dp
         a1 = 0.5545_dp
         s8 = 2.2609_dp
         a2 = 3.2297_dp
      CASE ("B98")
         s6 = 1.0000_dp
         a1 = 0.0000_dp
         s8 = 0.7086_dp
         a2 = 6.0672_dp
      CASE ("BLYP")
         s6 = 1.0000_dp
         a1 = 0.4298_dp
         s8 = 2.6996_dp
         a2 = 4.2359_dp
      CASE ("BP86")
         s6 = 1.0000_dp
         a1 = 0.3946_dp
         s8 = 3.2822_dp
         a2 = 4.8516_dp
      CASE ("DSD-BLYP")
         s6 = 0.5000_dp
         a1 = 0.0000_dp
         s8 = 0.2130_dp
         a2 = 6.0519_dp
      CASE ("HCTH407")
         s6 = 1.0000_dp
         a1 = 0.0000_dp
         s8 = 0.6490_dp
         a2 = 4.8162_dp
      CASE ("HISS")
         s6 = 1.0000_dp
         a1 = 0.0000_dp
         s8 = 1.6112_dp
         a2 = 7.3539_dp
      CASE ("HSE03")
         s6 = 1.0000_dp
         a1 = 0.0000_dp
         s8 = 1.1243_dp
         a2 = 6.8889_dp
      CASE ("HSE06")
         s6 = 1.0000_dp
         a1 = 0.3830_dp
         s8 = 2.3100_dp
         a2 = 5.6850_dp
      CASE ("M11")
         s6 = 1.0000_dp
         a1 = 0.0000_dp
         s8 = 2.8112_dp
         a2 = 10.1389_dp
      CASE ("MN12SX")
         s6 = 1.0000_dp
         a1 = 0.0983_dp
         s8 = 1.1674_dp
         a2 = 8.0259_dp
      CASE ("MN15")
         s6 = 1.0000_dp
         a1 = 2.0971_dp
         s8 = 0.7862_dp
         a2 = 7.5923_dp
      CASE ("mPWPW91")
         s6 = 1.0000_dp
         a1 = 0.3168_dp
         s8 = 1.7974_dp
         a2 = 4.7732_dp
      CASE ("MPW1PW91")
         s6 = 1.0000_dp
         a1 = 0.3342_dp
         s8 = 1.8744_dp
         a2 = 4.9819_dp
      CASE ("MPW1KCIS")
         s6 = 1.0000_dp
         a1 = 0.0576_dp
         s8 = 1.0893_dp
         a2 = 5.5314_dp
      CASE ("MPWKCIS1K")
         s6 = 1.0000_dp
         a1 = 0.0855_dp
         s8 = 1.2875_dp
         a2 = 5.8961_dp
      CASE ("N12SX")
         s6 = 1.0000_dp
         a1 = 0.3283_dp
         s8 = 2.4900_dp
         a2 = 5.7898_dp
      CASE ("O3LYP")
         s6 = 1.0000_dp
         a1 = 0.0963_dp
         s8 = 1.8171_dp
         a2 = 5.9940_dp
      CASE ("PBE0")
         s6 = 1.0000_dp
         a1 = 0.4145_dp
         s8 = 1.2177_dp
         a2 = 4.8593_dp
      CASE ("PBE")
         s6 = 1.0000_dp
         a1 = 0.4289_dp
         s8 = 0.7875_dp
         a2 = 4.4407_dp
      CASE ("PBEhPBE")
         s6 = 1.0000_dp
         a1 = 0.0000_dp
         s8 = 1.1152_dp
         a2 = 6.7184_dp
      CASE ("PBEh1PBE")
         s6 = 1.0000_dp
         a1 = 0.0000_dp
         s8 = 1.4877_dp
         a2 = 7.0385_dp
      CASE ("PBE1KCIS")
         s6 = 1.0000_dp
         a1 = 0.0000_dp
         s8 = 0.7688_dp
         a2 = 6.2794_dp
      CASE ("PW6B95")
         s6 = 1.0000_dp
         a1 = 0.2076_dp
         s8 = 0.7257_dp
         a2 = 6.3750_dp
      CASE ("PWPB95")
         s6 = 0.8200_dp
         a1 = 0.0000_dp
         s8 = 0.2904_dp
         a2 = 7.3141_dp
      CASE ("revPBE0")
         s6 = 1.0000_dp
         a1 = 0.4679_dp
         s8 = 1.7588_dp
         a2 = 3.7619_dp
      CASE ("revPBE38")
         s6 = 1.0000_dp
         a1 = 0.4309_dp
         s8 = 1.4760_dp
         a2 = 3.9446_dp
      CASE ("revPBE")
         s6 = 1.0000_dp
         a1 = 0.5238_dp
         s8 = 2.3550_dp
         a2 = 3.5016_dp
      CASE ("revTPSS")
         s6 = 1.0000_dp
         a1 = 0.4426_dp
         s8 = 1.4023_dp
         a2 = 4.4723_dp
      CASE ("revTPSS0")
         s6 = 1.0000_dp
         a1 = 0.2218_dp
         s8 = 1.6151_dp
         a2 = 5.7985_dp
      CASE ("revTPSSh")
         s6 = 1.0000_dp
         a1 = 0.2660_dp
         s8 = 1.4076_dp
         a2 = 5.3761_dp
      CASE ("RPBE")
         s6 = 1.0000_dp
         a1 = 0.1820_dp
         s8 = 0.8318_dp
         a2 = 4.0094_dp
      CASE ("RPW86PBE")
         s6 = 1.0000_dp
         a1 = 0.4613_dp
         s8 = 1.3845_dp
         a2 = 4.5062_dp
      CASE ("SCAN")
         s6 = 1.0000_dp
         a1 = 0.538_dp
         s8 = 0.0000_dp
         a2 = 5.420_dp
      CASE ("SOGGA11X")
         s6 = 1.0000_dp
         a1 = 0.1330_dp
         s8 = 1.1426_dp
         a2 = 5.7381_dp
      CASE ("TPSS0")
         s6 = 1.0000_dp
         a1 = 0.3768_dp
         s8 = 1.2576_dp
         a2 = 4.5865_dp
      CASE ("TPSS1KCIS")
         s6 = 1.0000_dp
         a1 = 0.0000_dp
         s8 = 1.0542_dp
         a2 = 6.0201_dp
      CASE ("TPSS")
         s6 = 1.0000_dp
         a1 = 0.4535_dp
         s8 = 1.9435_dp
         a2 = 4.4752_dp
      CASE ("tHCTH")
         s6 = 1.0000_dp
         a1 = 0.0000_dp
         s8 = 1.2626_dp
         a2 = 5.6162_dp
      CASE ("tHCTHhyb")
         s6 = 1.0000_dp
         a1 = 0.0000_dp
         s8 = 0.9585_dp
         a2 = 6.2303_dp
      CASE ("XLYP")
         s6 = 1.0000_dp
         a1 = 0.0809_dp
         s8 = 1.5669_dp
         a2 = 5.3166_dp
      CASE ("X3LYP")
         s6 = 1.0000_dp
         a1 = 0.2022_dp
         s8 = 1.5744_dp
         a2 = 5.4184_dp
      END SELECT

   END SUBROUTINE qs_scaling_dftd3bj

! **************************************************************************************************
!> \brief ...
!> \param s6 ...
!> \param a1 ...
!> \param s8 ...
!> \param a2 ...
!> \param vdw_section ...
! **************************************************************************************************
   SUBROUTINE qs_scaling_dftd4(s6, a1, s8, a2, vdw_section)
      REAL(KIND=dp), INTENT(INOUT)                       :: s6, a1, s8, a2
      TYPE(section_vals_type), POINTER                   :: vdw_section

      CHARACTER(LEN=default_string_length)               :: functional

      CALL section_vals_val_get(vdw_section, "PAIR_POTENTIAL%REFERENCE_FUNCTIONAL", c_val=functional)
      CALL uppercase(functional)

      ! DFT-D4 BJ-ATM rational damping parameters (for use with C9 three-body term).
      ! References: SAW190103/SAW190107/SAW191202/SAW211021/SAW211107 = S.A. Wildman fits;
      !             DOI references given inline for more recent parametrizations.
      SELECT CASE (TRIM(functional))
      CASE DEFAULT
         CPABORT("No DFT-D4 parameters available for this functional: "//TRIM(functional))
      CASE ("HF")
         s6 = 1.0000_dp
         s8 = 1.61679827_dp
         a1 = 0.44959224_dp
         a2 = 3.35743605_dp
      CASE ("AM05")
         s6 = 1.0000_dp
         s8 = 1.71885838_dp
         a1 = 0.47901431_dp
         a2 = 5.96771581_dp
      CASE ("BLYP")
         s6 = 1.0000_dp
         s8 = 2.34076671_dp
         a1 = 0.44488865_dp
         a2 = 4.09330090_dp
      CASE ("BPBE")
         s6 = 1.0000_dp
         s8 = 3.64405246_dp
         a1 = 0.52905620_dp
         a2 = 4.11311891_dp
      CASE ("BP86", "BP")
         s6 = 1.0000_dp
         s8 = 3.35497927_dp
         a1 = 0.43645861_dp
         a2 = 4.92406854_dp
      CASE ("BPW", "BPW91")
         s6 = 1.0000_dp
         s8 = 3.24571506_dp
         a1 = 0.50050454_dp
         a2 = 4.12346483_dp
      CASE ("GLYP")
         s6 = 1.0000_dp
         s8 = 4.23798924_dp
         a1 = 0.38426465_dp
         a2 = 4.38412863_dp
      CASE ("LB94")
         s6 = 1.0000_dp
         s8 = 2.59538499_dp
         a1 = 0.42088944_dp
         a2 = 3.28193223_dp
      CASE ("MPWLYP")
         s6 = 1.0000_dp
         s8 = 1.25842942_dp
         a1 = 0.25773894_dp
         a2 = 5.02319542_dp
      CASE ("MPWPW91", "MPWPW")
         s6 = 1.0000_dp
         s8 = 1.82596836_dp
         a1 = 0.34526745_dp
         a2 = 4.84620734_dp
      CASE ("OLYP")
         s6 = 1.0000_dp
         s8 = 2.74836820_dp
         a1 = 0.60184498_dp
         a2 = 2.53292167_dp
      CASE ("OPBE")
         s6 = 1.0000_dp
         s8 = 3.06917417_dp
         a1 = 0.68267534_dp
         a2 = 2.22849018_dp
      CASE ("PBE")
         s6 = 1.0000_dp
         s8 = 0.95948085_dp
         a1 = 0.38574991_dp
         a2 = 4.80688534_dp
      CASE ("PBESOL")
         s6 = 1.0000_dp
         s8 = 1.71885698_dp
         a1 = 0.47901421_dp
         a2 = 5.96771589_dp
      CASE ("PW86PBE")
         s6 = 1.0000_dp
         s8 = 1.21362856_dp
         a1 = 0.40510366_dp
         a2 = 4.66737724_dp
      CASE ("PW91")
         s6 = 1.0000_dp
         s8 = 0.77283111_dp
         a1 = 0.39581542_dp
         a2 = 4.93405761_dp
      CASE ("PWP", "PW91P86")
         s6 = 1.0000_dp
         s8 = 0.32801227_dp
         a1 = 0.35874687_dp
         a2 = 6.05861168_dp
      CASE ("REVPBE")
         s6 = 1.0000_dp
         s8 = 1.74676530_dp
         a1 = 0.53634900_dp
         a2 = 3.07261485_dp
      CASE ("RPBE")
         s6 = 1.0000_dp
         s8 = 1.31183787_dp
         a1 = 0.46169493_dp
         a2 = 3.15711757_dp
      CASE ("RPW86PBE")
         s6 = 1.0000_dp
         s8 = 1.12624034_dp
         a1 = 0.38151218_dp
         a2 = 4.75480472_dp
      CASE ("XLYP")
         s6 = 1.0000_dp
         s8 = 1.62972054_dp
         a1 = 0.11268673_dp
         a2 = 5.40786417_dp
      CASE ("B97D", "B97-D")
         s6 = 1.0000_dp
         s8 = 1.69460052_dp
         a1 = 0.28904684_dp
         a2 = 4.13407323_dp
      CASE ("B97")
         s6 = 1.0000_dp
         s8 = 0.87854260_dp
         a1 = 0.29319126_dp
         a2 = 4.51647719_dp
      CASE ("B97M")
         s6 = 1.0000_dp
         s8 = 0.6633_dp
         a1 = 0.4288_dp
         a2 = 3.9935_dp
      CASE ("M06L")
         s6 = 1.0000_dp
         s8 = 0.59493760_dp
         a1 = 0.71422359_dp
         a2 = 6.35314182_dp
      CASE ("REVTPSS")
         s6 = 1.0000_dp
         s8 = 1.53089454_dp
         a1 = 0.44880597_dp
         a2 = 4.64042317_dp
      CASE ("RSCAN")
         s6 = 1.0000_dp
         s8 = 0.87728975_dp
         a1 = 0.49116966_dp
         a2 = 5.75859346_dp
      CASE ("R2SCAN")
         s6 = 1.0000_dp
         s8 = 0.60187490_dp
         a1 = 0.51559235_dp
         a2 = 5.77342911_dp
      CASE ("SCAN")
         s6 = 1.0000_dp
         s8 = 1.46126056_dp
         a1 = 0.62930855_dp
         a2 = 6.31284039_dp
      CASE ("TPSS")
         s6 = 1.0000_dp
         s8 = 1.76596355_dp
         a1 = 0.42822303_dp
         a2 = 4.54257102_dp
      CASE ("B1B95")
         s6 = 1.0000_dp
         s8 = 1.27701162_dp
         a1 = 0.40554715_dp
         a2 = 4.63323074_dp
      CASE ("B1LYP")
         s6 = 1.0000_dp
         s8 = 1.98553711_dp
         a1 = 0.39309040_dp
         a2 = 4.55465145_dp
      CASE ("B1P86", "B1P")
         s6 = 1.0000_dp
         s8 = 3.36115015_dp
         a1 = 0.48665293_dp
         a2 = 5.05219572_dp
      CASE ("B1PW91", "B1PW")
         s6 = 1.0000_dp
         s8 = 3.02227550_dp
         a1 = 0.47396846_dp
         a2 = 4.49845309_dp
      CASE ("B3LYP")
         s6 = 1.0000_dp
         s8 = 2.02929367_dp
         a1 = 0.40868035_dp
         a2 = 4.53807137_dp
      CASE ("B3P86", "B3P")
         s6 = 1.0000_dp
         s8 = 3.08822155_dp
         a1 = 0.47324238_dp
         a2 = 4.98682134_dp
      CASE ("B3PW91", "B3PW")
         s6 = 1.0000_dp
         s8 = 2.88364295_dp
         a1 = 0.46990860_dp
         a2 = 4.51641422_dp
      CASE ("BHLYP")
         s6 = 1.0000_dp
         s8 = 1.65281646_dp
         a1 = 0.27263660_dp
         a2 = 5.48634586_dp
      CASE ("CAMB3LYP", "CAM-B3LYP")
         s6 = 1.0000_dp
         s8 = 1.66041301_dp
         a1 = 0.40267156_dp
         a2 = 5.17432195_dp
      CASE ("CAMQTP01", "CAM-QTP01")
         s6 = 1.0000_dp
         s8 = 1.156_dp
         a1 = 0.461_dp
         a2 = 6.375_dp
      CASE ("HSE03")
         s6 = 1.0000_dp
         s8 = 1.19812280_dp
         a1 = 0.38662939_dp
         a2 = 5.22925796_dp
      CASE ("HSE06")
         s6 = 1.0000_dp
         s8 = 1.19528249_dp
         a1 = 0.38663183_dp
         a2 = 5.19133469_dp
      CASE ("HSE12")
         s6 = 1.0000_dp
         s8 = 1.23500792_dp
         a1 = 0.39226921_dp
         a2 = 5.22036266_dp
      CASE ("HSE12S")
         s6 = 1.0000_dp
         s8 = 1.23767762_dp
         a1 = 0.39989137_dp
         a2 = 5.34809245_dp
      CASE ("HSESOL")
         s6 = 1.0000_dp
         s8 = 1.82207807_dp
         a1 = 0.45646268_dp
         a2 = 5.59662251_dp
      CASE ("LCBLYP", "LC-BLYP")
         s6 = 1.0000_dp
         s8 = 1.60344180_dp
         a1 = 0.45769839_dp
         a2 = 7.86924893_dp
      CASE ("LCWPBE", "LC-WPBE")
         s6 = 1.0000_dp
         s8 = 1.170_dp
         a1 = 0.378_dp
         a2 = 4.816_dp
      CASE ("LCWPBEH", "LC-WPBEH")
         s6 = 1.0000_dp
         s8 = 1.318_dp
         a1 = 0.386_dp
         a2 = 5.010_dp
      CASE ("M06")
         s6 = 1.0000_dp
         s8 = 0.16366729_dp
         a1 = 0.53456413_dp
         a2 = 6.06192174_dp
      CASE ("MN12SX", "MN12-SX")
         s6 = 1.0000_dp
         s8 = 0.85964873_dp
         a1 = 0.62662681_dp
         a2 = 5.62088906_dp
      CASE ("MPW1B95")
         s6 = 1.0000_dp
         s8 = 0.50093024_dp
         a1 = 0.41585097_dp
         a2 = 4.99154869_dp
      CASE ("MPW1LYP")
         s6 = 1.0000_dp
         s8 = 1.15591153_dp
         a1 = 0.25603493_dp
         a2 = 5.32083895_dp
      CASE ("MPW1PW91", "MPW1PW")
         s6 = 1.0000_dp
         s8 = 1.80841716_dp
         a1 = 0.42961819_dp
         a2 = 4.68892341_dp
      CASE ("MPWB1K")
         s6 = 1.0000_dp
         s8 = 0.57338313_dp
         a1 = 0.44687975_dp
         a2 = 5.21266777_dp
      CASE ("O3LYP")
         s6 = 1.0000_dp
         s8 = 1.75762508_dp
         a1 = 0.10348980_dp
         a2 = 6.16233282_dp
      CASE ("PBE0")
         s6 = 1.0000_dp
         s8 = 1.20065498_dp
         a1 = 0.40085597_dp
         a2 = 5.02928789_dp
      CASE ("PBE0-2", "PBE02")
         s6 = 0.5000_dp
         s8 = 0.64299082_dp
         a1 = 0.76542115_dp
         a2 = 5.78578675_dp
      CASE ("PBE0-DH", "PBE0DH")
         s6 = 0.8750_dp
         s8 = 0.96811578_dp
         a1 = 0.47592488_dp
         a2 = 5.08622873_dp
      CASE ("PW1PW")
         s6 = 1.0000_dp
         s8 = 0.96850170_dp
         a1 = 0.42427511_dp
         a2 = 5.02060636_dp
      CASE ("PW6B95")
         s6 = 1.0000_dp
         s8 = -0.31926054_dp
         a1 = 0.04142919_dp
         a2 = 5.84655608_dp
      CASE ("PWP1")
         s6 = 1.0000_dp
         s8 = 0.60492565_dp
         a1 = 0.46855837_dp
         a2 = 5.76921413_dp
      CASE ("REVPBE0")
         s6 = 1.0000_dp
         s8 = 1.57185414_dp
         a1 = 0.38705966_dp
         a2 = 4.11028876_dp
      CASE ("REVPBE0-DH", "REVPBE0DH")
         s6 = 0.8750_dp
         s8 = 1.24456037_dp
         a1 = 0.36730560_dp
         a2 = 4.71126482_dp
      CASE ("REVPBE38")
         s6 = 1.0000_dp
         s8 = 1.66597472_dp
         a1 = 0.39476833_dp
         a2 = 4.39026628_dp
      CASE ("WB97", "OMEGAB97")
         s6 = 1.0000_dp
         s8 = 6.55792598_dp
         a1 = 0.76666802_dp
         a2 = 8.36027334_dp
      CASE ("WB97X", "OMEGAB97X")
         s6 = 1.0000_dp
         s8 = 0.5093_dp
         a1 = 0.0662_dp
         a2 = 5.4487_dp
      CASE ("WB97X-2008", "WB97X_2008")
         s6 = 1.0000_dp
         s8 = -0.07519516_dp
         a1 = 0.45094893_dp
         a2 = 6.78425255_dp
      CASE ("WB97X-REV", "WB97X_REV")
         s6 = 1.0000_dp
         s8 = 0.4485_dp
         a1 = 0.3306_dp
         a2 = 4.279_dp
      CASE ("WB97X-3C", "WB97X_3C")
         s6 = 1.0000_dp
         s8 = 0.0_dp
         a1 = 0.2464_dp
         a2 = 4.737_dp
      CASE ("WB97M", "OMEGAB97M")
         s6 = 1.0000_dp
         s8 = 0.7761_dp
         a1 = 0.7514_dp
         a2 = 2.7099_dp
      CASE ("WB97M-REV", "WB97M_REV")
         s6 = 1.0000_dp
         s8 = 0.842_dp
         a1 = 0.359_dp
         a2 = 4.668_dp
      CASE ("X3LYP")
         s6 = 1.0000_dp
         s8 = 1.54701429_dp
         a1 = 0.20318443_dp
         a2 = 5.61852648_dp
      CASE ("REVTPSS0")
         s6 = 1.0000_dp
         s8 = 1.54664499_dp
         a1 = 0.45890964_dp
         a2 = 4.78426405_dp
      CASE ("REVTPSSH")
         s6 = 1.0000_dp
         s8 = 1.52740307_dp
         a1 = 0.45161957_dp
         a2 = 4.70779483_dp
      CASE ("TPSS0")
         s6 = 1.0000_dp
         s8 = 1.62438102_dp
         a1 = 0.40329022_dp
         a2 = 4.80537871_dp
      CASE ("TPSSH")
         s6 = 1.0000_dp
         s8 = 1.85897750_dp
         a1 = 0.44286966_dp
         a2 = 4.60230534_dp
      CASE ("PWPB95")
         s6 = 0.8200_dp
         s8 = -0.34639127_dp
         a1 = 0.41080636_dp
         a2 = 3.83878274_dp
      CASE ("R2SCANH")
         s6 = 1.0000_dp
         s8 = 0.8324_dp
         a1 = 0.4944_dp
         a2 = 5.9019_dp
      CASE ("R2SCAN0")
         s6 = 1.0000_dp
         s8 = 0.8992_dp
         a1 = 0.4778_dp
         a2 = 5.8779_dp
      CASE ("R2SCAN50")
         s6 = 1.0000_dp
         s8 = 1.0471_dp
         a1 = 0.4574_dp
         a2 = 5.8969_dp
      CASE ("R2SCAN-3C", "R2SCAN_3C")
         s6 = 1.00_dp
         s8 = 0.00_dp
         a1 = 0.42_dp
         a2 = 5.65_dp
      CASE ("WR2SCAN")
         s6 = 1.0000_dp
         s8 = 1.0000_dp
         a1 = 0.3834_dp
         a2 = 5.7889_dp
      CASE ("R2SCAN0-DH", "R2SCAN0DH")
         s6 = 0.9424_dp
         s8 = 0.3856_dp
         a1 = 0.4271_dp
         a2 = 5.8565_dp
      CASE ("R2SCAN-CIDH", "R2SCANCIDH")
         s6 = 0.8666_dp
         s8 = 0.5336_dp
         a1 = 0.4171_dp
         a2 = 5.9125_dp
      CASE ("R2SCAN-QIDH", "R2SCANQIDH")
         s6 = 0.7867_dp
         s8 = 0.2955_dp
         a1 = 0.4001_dp
         a2 = 5.8300_dp
      CASE ("R2SCAN0-2", "R2SCAN02")
         s6 = 0.7386_dp
         s8 = 0.0000_dp
         a1 = 0.4030_dp
         a2 = 5.5142_dp
      CASE ("PR2SCAN50")
         s6 = 0.7964_dp
         s8 = 0.3421_dp
         a1 = 0.4663_dp
         a2 = 5.7916_dp
      CASE ("PR2SCAN69")
         s6 = 0.7167_dp
         s8 = 0.0000_dp
         a1 = 0.4644_dp
         a2 = 5.2563_dp
      CASE ("KPR2SCAN50")
         s6 = 0.8402_dp
         s8 = 0.1212_dp
         a1 = 0.4382_dp
         a2 = 5.8232_dp
      CASE ("WPR2SCAN50")
         s6 = 0.8143_dp
         s8 = 0.3842_dp
         a1 = 0.4135_dp
         a2 = 5.8773_dp
      CASE ("LH07TSVWN", "LH07T-SVWN")
         s6 = 1.0000_dp
         s8 = 2.09333001_dp
         a1 = 0.35025189_dp
         a2 = 4.34166515_dp
      CASE ("LH07SSVWN", "LH07S-SVWN")
         s6 = 1.0000_dp
         s8 = 3.16675531_dp
         a1 = 0.35965552_dp
         a2 = 4.31947614_dp
      CASE ("LH12CTSSIRPW92", "LH12CT-SSIRPW92")
         s6 = 1.0000_dp
         s8 = 2.48973402_dp
         a1 = 0.34026075_dp
         a2 = 3.96948081_dp
      CASE ("LH12CTSSIFPW92", "LH12CT-SSIFPW92")
         s6 = 1.0000_dp
         s8 = 2.68467610_dp
         a1 = 0.34190416_dp
         a2 = 3.91039666_dp
      CASE ("LH14TCALPBE", "LH14T-CALPBE")
         s6 = 1.0000_dp
         s8 = 1.28130770_dp
         a1 = 0.38822021_dp
         a2 = 4.92501211_dp
      CASE ("LH20T")
         s6 = 1.000_dp
         s8 = 0.113_dp
         a1 = 0.479_dp
         a2 = 4.635_dp
      CASE ("B2PLYP")
         s6 = 0.6400_dp
         s8 = 1.16888646_dp
         a1 = 0.44154604_dp
         a2 = 4.73114642_dp
      CASE ("B2GPPLYP", "B2GP-PLYP")
         s6 = 0.5600_dp
         s8 = 0.94633372_dp
         a1 = 0.42907301_dp
         a2 = 5.18802602_dp
      CASE ("MPW2PLYP")
         s6 = 0.7500_dp
         s8 = 0.45788846_dp
         a1 = 0.42997704_dp
         a2 = 5.07650682_dp
      CASE ("DSDBLYP", "DSD-BLYP")
         s6 = 0.5400_dp
         s8 = 0.63018237_dp
         a1 = 0.47591835_dp
         a2 = 4.73713781_dp
      CASE ("DSDPBE", "DSD-PBE")
         s6 = 0.4500_dp
         s8 = 0.70584116_dp
         a1 = 0.45787085_dp
         a2 = 4.44566742_dp
      CASE ("DSDPBEB95", "DSD-PBEB95")
         s6 = 0.5400_dp
         s8 = -0.14668670_dp
         a1 = 0.46394587_dp
         a2 = 3.64913860_dp
      CASE ("DSDPBEP86", "DSD-PBEP86")
         s6 = 0.4700_dp
         s8 = 0.37586675_dp
         a1 = 0.53698768_dp
         a2 = 5.13022435_dp
      CASE ("DSDSVWN", "DSD-SVWN")
         s6 = 0.4100_dp
         s8 = 0.72914436_dp
         a1 = 0.51347412_dp
         a2 = 5.11858541_dp
      CASE ("DODBLYP", "DOD-BLYP")
         s6 = 0.4700_dp
         s8 = 1.31146043_dp
         a1 = 0.43407294_dp
         a2 = 4.27914360_dp
      CASE ("DODPBE", "DOD-PBE")
         s6 = 0.4800_dp
         s8 = 0.92051454_dp
         a1 = 0.43037052_dp
         a2 = 4.38067238_dp
      CASE ("DODPBEB95", "DOD-PBEB95")
         s6 = 0.5600_dp
         s8 = 0.01574635_dp
         a1 = 0.43745720_dp
         a2 = 3.69180763_dp
      CASE ("DODPBEP86", "DOD-PBEP86")
         s6 = 0.4600_dp
         s8 = 0.71405681_dp
         a1 = 0.42408665_dp
         a2 = 4.52884439_dp
      CASE ("DODSVWN", "DOD-SVWN")
         s6 = 0.4200_dp
         s8 = 0.94500207_dp
         a1 = 0.47449026_dp
         a2 = 5.05316093_dp
      CASE ("REVDSDBLYP", "REVDSD-BLYP")
         s6 = 0.6141_dp
         s8 = 0.0000_dp
         a1 = 0.3800_dp
         a2 = 3.5200_dp
      CASE ("REVDSDPBE", "REVDSD-PBE", "REVDSD-PBEPBE")
         s6 = 0.6706_dp
         s8 = 0.0000_dp
         a1 = 0.4000_dp
         a2 = 3.6000_dp
      CASE ("REVDSDPBEP86", "REVDSD-PBEP86")
         s6 = 0.5132_dp
         s8 = 0.0000_dp
         a1 = 0.4400_dp
         a2 = 3.6000_dp
      CASE ("REVDODPBEP86", "REVDOD-PBEP86")
         s6 = 0.5552_dp
         s8 = 0.0000_dp
         a1 = 0.4400_dp
         a2 = 3.6000_dp
      END SELECT

   END SUBROUTINE qs_scaling_dftd4

! **************************************************************************************************
!> \brief ...
!> \param cell ...
!> \param ncell ...
!> \return ...
! **************************************************************************************************
   FUNCTION cellhash(cell, ncell) RESULT(hash)
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell, ncell
      INTEGER                                            :: hash

      INTEGER                                            :: ix, iy, iz, nx, ny, nz

      CPASSERT(ALL(ABS(cell) <= ncell))

      ix = cell(1)
      IF (ix /= 0) THEN
         ix = 2*ABS(ix) - (1 + SIGN(1, ix))/2
      END IF
      iy = cell(2)
      IF (iy /= 0) THEN
         iy = 2*ABS(iy) - (1 + SIGN(1, iy))/2
      END IF
      iz = cell(3)
      IF (iz /= 0) THEN
         iz = 2*ABS(iz) - (1 + SIGN(1, iz))/2
      END IF

      nx = 2*ncell(1) + 1
      ny = 2*ncell(2) + 1
      nz = 2*ncell(3) + 1

      hash = ix*ny*nz + iy*nz + iz + 1

   END FUNCTION cellhash
! **************************************************************************************************

END MODULE qs_dispersion_utils

