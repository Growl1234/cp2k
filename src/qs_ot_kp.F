!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Orbital transformation support routines for k-point calculations
! **************************************************************************************************
MODULE qs_ot_kp
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_log_handling,                 ONLY: cp_to_string
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_mo_types,                     ONLY: get_mo_set
   USE qs_scf_diagonalization,          ONLY: diag_kp_smat
   USE qs_scf_types,                    ONLY: qs_scf_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ot_kp'
   REAL(KIND=dp), PARAMETER, PRIVATE    :: eps_ortho_check = 1.0E-8_dp

   PUBLIC :: ot_scf_mini_kp

CONTAINS

! **************************************************************************************************
!> \brief Prepare the first k-point OT state. The actual OT update is intentionally not implemented
!>        in this preparatory step.
!> \param matrix_s   overlap matrices in real-space image representation
!> \param kpoints    k-point environment
!> \param scf_env    SCF environment
! **************************************************************************************************
   SUBROUTINE ot_scf_mini_kp(matrix_s, kpoints, scf_env)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ot_scf_mini_kp'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: max_deviation
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: fmwork

      CALL timeset(routineN, handle)

      fmwork => scf_env%scf_work1

      ! Build S(k)^(1/2) for all local k-points and use it to verify that the
      ! diagonalization-based initial occupied subspaces satisfy C(k)^H S(k) C(k) = I.
      CALL diag_kp_smat(matrix_s, kpoints, fmwork)
      CALL ot_kp_check_orthonormality(kpoints, max_deviation)
      IF (max_deviation > eps_ortho_check) THEN
         CALL cp_abort(__LOCATION__, &
                       "KP code: initial OT-k subspace is not orthonormal; max |C^H S C - I| = "// &
                       cp_to_string(max_deviation))
      END IF

      CALL timestop(handle)

      CALL cp_abort(__LOCATION__, &
                    "KP code: OT step for k-points not implemented yet")

   END SUBROUTINE ot_scf_mini_kp

! **************************************************************************************************
!> \brief Pack the real/imaginary k-point MO coefficients into a complex full matrix.
!> \param kp            local k-point environment
!> \param ispin         spin index
!> \param use_real_wfn  whether the k-point calculation uses real wave functions
!> \param cmos          complex MO coefficient matrix
! **************************************************************************************************
   SUBROUTINE ot_kp_pack_mos(kp, ispin, use_real_wfn, cmos)

      TYPE(kpoint_env_type), POINTER                     :: kp
      INTEGER, INTENT(IN)                                :: ispin
      LOGICAL, INTENT(IN)                                :: use_real_wfn
      TYPE(cp_cfm_type), INTENT(IN)                      :: cmos

      TYPE(cp_fm_type), POINTER                          :: imos, rmos

      CALL get_mo_set(kp%mos(1, ispin), mo_coeff=rmos)
      IF (use_real_wfn) THEN
         CALL cp_fm_to_cfm(msourcer=rmos, mtarget=cmos)
      ELSE
         CALL get_mo_set(kp%mos(2, ispin), mo_coeff=imos)
         CALL cp_fm_to_cfm(msourcer=rmos, msourcei=imos, mtarget=cmos)
      END IF

   END SUBROUTINE ot_kp_pack_mos

! **************************************************************************************************
!> \brief Unpack a complex full matrix into the real/imaginary k-point MO coefficient storage.
!> \param cmos          complex MO coefficient matrix
!> \param kp            local k-point environment
!> \param ispin         spin index
!> \param use_real_wfn  whether the k-point calculation uses real wave functions
! **************************************************************************************************
   SUBROUTINE ot_kp_unpack_mos(cmos, kp, ispin, use_real_wfn)

      TYPE(cp_cfm_type), INTENT(IN)                      :: cmos
      TYPE(kpoint_env_type), POINTER                     :: kp
      INTEGER, INTENT(IN)                                :: ispin
      LOGICAL, INTENT(IN)                                :: use_real_wfn

      TYPE(cp_fm_type), POINTER                          :: imos, rmos

      CALL get_mo_set(kp%mos(1, ispin), mo_coeff=rmos)
      IF (use_real_wfn) THEN
         CALL cp_cfm_to_fm(cmos, mtargetr=rmos)
      ELSE
         CALL get_mo_set(kp%mos(2, ispin), mo_coeff=imos)
         CALL cp_cfm_to_fm(cmos, mtargetr=rmos, mtargeti=imos)
      END IF

   END SUBROUTINE ot_kp_unpack_mos

! **************************************************************************************************
!> \brief Verify C(k)^H S(k) C(k) = I using the S(k)^(1/2) matrices stored in kpoint_env.
!> \param kpoints       k-point environment
!> \param max_deviation maximum absolute deviation from the unit overlap matrix
! **************************************************************************************************
   SUBROUTINE ot_kp_check_orthonormality(kpoints, max_deviation)

      TYPE(kpoint_type), POINTER                         :: kpoints
      REAL(KIND=dp), INTENT(OUT)                         :: max_deviation

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ot_kp_check_orthonormality'
      COMPLEX(KIND=dp), PARAMETER                        :: cone = (1.0_dp, 0.0_dp), &
                                                            czero = (0.0_dp, 0.0_dp)

      COMPLEX(KIND=dp)                                   :: unit_element
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER, &
         CONTIGUOUS                                      :: p_local
      INTEGER                                            :: handle, icol, ikp, ispin, irow, &
                                                            kplocal, nao, ncol_local, nmo, &
                                                            nrow_local, nspin
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: use_real_wfn
      REAL(KIND=dp)                                      :: local_deviation
      TYPE(cp_cfm_type)                                  :: cmos, coverlap, cs_half, scmos
      TYPE(cp_fm_struct_type), POINTER                   :: mo_struct, s_struct
      TYPE(cp_fm_type), POINTER                          :: rmos
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (col_indices, mo_struct, para_env, row_indices, rmos, s_struct)
      max_deviation = 0.0_dp
      CALL get_kpoint_info(kpoints, kp_range=kp_range, use_real_wfn=use_real_wfn)
      kplocal = kp_range(2) - kp_range(1) + 1
      kp => kpoints%kp_env(1)%kpoint_env
      nspin = SIZE(kp%mos, 2)

      DO ikp = 1, kplocal
         kp => kpoints%kp_env(ikp)%kpoint_env
         IF (use_real_wfn) THEN
            CALL cp_fm_get_info(kp%shalf, matrix_struct=s_struct, nrow_global=nao)
            CALL cp_cfm_create(cs_half, s_struct)
            CALL cp_fm_to_cfm(msourcer=kp%shalf, mtarget=cs_half)
         ELSE
            CALL cp_cfm_get_info(kp%cshalf, matrix_struct=s_struct, nrow_global=nao)
            CALL cp_cfm_create(cs_half, s_struct)
            CALL cp_cfm_to_cfm(kp%cshalf, cs_half)
         END IF

         DO ispin = 1, nspin
            CALL get_mo_set(kp%mos(1, ispin), mo_coeff=rmos)
            CALL cp_fm_get_info(rmos, matrix_struct=mo_struct, ncol_global=nmo)

            CALL cp_cfm_create(cmos, mo_struct)
            CALL cp_cfm_create(scmos, mo_struct)
            CALL cp_cfm_create(coverlap, mo_struct, nrow=nmo, ncol=nmo, set_zero=.TRUE.)

            CALL ot_kp_pack_mos(kp, ispin, use_real_wfn, cmos)
            CALL cp_cfm_gemm('N', 'N', nao, nmo, nao, cone, cs_half, cmos, czero, scmos)
            CALL cp_cfm_gemm('C', 'N', nmo, nmo, nao, cone, scmos, scmos, czero, coverlap)

            CALL cp_cfm_get_info(coverlap, nrow_local=nrow_local, ncol_local=ncol_local, &
                                 row_indices=row_indices, col_indices=col_indices, &
                                 local_data=p_local, para_env=para_env)
            local_deviation = 0.0_dp
            DO icol = 1, ncol_local
               DO irow = 1, nrow_local
                  IF (row_indices(irow) == col_indices(icol)) THEN
                     unit_element = cone
                  ELSE
                     unit_element = czero
                  END IF
                  local_deviation = MAX(local_deviation, &
                                        ABS(p_local(irow, icol) - unit_element))
               END DO
            END DO
            CALL para_env%max(local_deviation)
            max_deviation = MAX(max_deviation, local_deviation)

            CALL cp_cfm_release(coverlap)
            CALL cp_cfm_release(scmos)
            CALL cp_cfm_release(cmos)
         END DO

         CALL cp_cfm_release(cs_half)
      END DO

      IF (ASSOCIATED(kpoints%para_env_inter_kp)) CALL kpoints%para_env_inter_kp%max(max_deviation)

      CALL timestop(handle)

   END SUBROUTINE ot_kp_check_orthonormality

END MODULE qs_ot_kp
