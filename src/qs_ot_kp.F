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
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm,&
                                              cp_cfm_scale_and_add,&
                                              cp_cfm_scale_and_add_fm
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_desymmetrize,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_type_no_symmetry,&
                                              dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_log_handling,                 ONLY: cp_to_string
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              fm_pool_create_fm,&
                                              fm_pool_give_back_fm
   USE cp_fm_types,                     ONLY: copy_info_type,&
                                              cp_fm_cleanup_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_finish_copy_general,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_start_copy_general,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE mathconstants,                   ONLY: gaussi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_matrix_pools,                 ONLY: mpools_get,&
                                              qs_matrix_pools_type
   USE qs_mo_types,                     ONLY: get_mo_set
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_scf_diagonalization,          ONLY: diag_kp_smat
   USE qs_scf_types,                    ONLY: qs_scf_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ot_kp'
   REAL(KIND=dp), PARAMETER, PRIVATE    :: eps_ortho_check = 1.0E-8_dp, &
                                             eps_hamiltonian_check = 1.0E-6_dp, &
                                             eps_tangent_check = 1.0E-8_dp

   PUBLIC :: ot_scf_mini_kp

CONTAINS

! **************************************************************************************************
!> \brief Prepare the first k-point OT state. The actual OT update is intentionally not implemented
!>        in this preparatory step.
!> \param matrix_ks  Kohn-Sham matrices in real-space image representation
!> \param matrix_s   overlap matrices in real-space image representation
!> \param kpoints    k-point environment
!> \param scf_env    SCF environment
! **************************************************************************************************
   SUBROUTINE ot_scf_mini_kp(matrix_ks, matrix_s, kpoints, scf_env)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ot_scf_mini_kp'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: max_deviation, max_descent_deviation, &
                                                            max_eval_deviation, max_gradient_norm, &
                                                            max_offdiag, max_residual, &
                                                            max_tangent_deviation
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

      ! Rebuild H(k), apply it to the initialized MOs, and verify that C(k)^H H(k) C(k)
      ! is consistent with the eigenvalues produced by the initialization diagonalization.
      CALL ot_kp_check_hamiltonian(matrix_ks, matrix_s, kpoints, fmwork, max_offdiag, &
                                      max_eval_deviation, max_residual, max_tangent_deviation, &
                                      max_gradient_norm, max_descent_deviation)
      IF (MAX(max_offdiag, max_eval_deviation, max_residual) > eps_hamiltonian_check .OR. &
          max_tangent_deviation > eps_tangent_check .OR. &
          max_descent_deviation > eps_tangent_check) THEN
         CALL cp_abort(__LOCATION__, &
                       "KP code: initial OT-k Hamiltonian/gradient check failed; max offdiag = "// &
                       cp_to_string(max_offdiag)//", max eval deviation = "// &
                       cp_to_string(max_eval_deviation)//", max residual = "// &
                       cp_to_string(max_residual)//", max |C^H G| = "// &
                       cp_to_string(max_tangent_deviation)//", max gradient norm = "// &
                       cp_to_string(max_gradient_norm)//", max descent deviation = "// &
                       cp_to_string(max_descent_deviation))
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


! **************************************************************************************************
!> \brief Verify H(k)C(k) by checking that initialized MOs diagonalize the current H(k).
!> \param matrix_ks          Kohn-Sham matrices in real-space image representation
!> \param matrix_s           overlap matrices in real-space image representation
!> \param kpoints            k-point environment
!> \param max_offdiag        maximum absolute off-diagonal element of C^H H C
!> \param max_eval_deviation maximum absolute deviation between diag(C^H H C) and eigenvalues
!> \param max_residual          maximum absolute residual of H C - S C (C^H H C)
!> \param max_tangent_deviation maximum absolute deviation of C^H G from zero
!> \param max_gradient_norm     maximum Frobenius norm of the tentative tangent gradient
!> \param max_descent_deviation maximum deviation of <G,-G> from -||G||^2
! **************************************************************************************************
   SUBROUTINE ot_kp_check_hamiltonian(matrix_ks, matrix_s, kpoints, fmwork, max_offdiag, &
                                      max_eval_deviation, max_residual, max_tangent_deviation, &
                                      max_gradient_norm, max_descent_deviation)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: fmwork
      REAL(KIND=dp), INTENT(OUT)                         :: max_offdiag, max_eval_deviation, &
                                                            max_residual, max_tangent_deviation, &
                                                            max_gradient_norm, max_descent_deviation

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ot_kp_check_hamiltonian'
      COMPLEX(KIND=dp), PARAMETER                        :: cone = (1.0_dp, 0.0_dp), &
                                                            czero = (0.0_dp, 0.0_dp)

      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER, CONTIGUOUS :: d_local, g_local, h_local, r_local
      INTEGER                                            :: handle, icol, igroup, ik, ikp, &
                                                            indx, ispin, irow, kplocal, &
                                                            nao, ncol_local, nmo, &
                                                            nrow_local, nspin, nkp, &
                                                            nkp_groups
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :), POINTER                  :: kp_dist
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: my_kpgrp, use_real_wfn
      REAL(KIND=dp)                                      :: descent_dev_local, descent_local, eval_dev_local, &
                                                            gradient_norm_local, gradient_norm_sq_local, &
                                                            offdiag_local, residual_local, tangent_local
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(copy_info_type), ALLOCATABLE, DIMENSION(:, :) :: info
      TYPE(cp_cfm_type)                                  :: cksmat, cmos, csmat, dxmos, gsmall, hcmos, hsmall, scmos
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_ao_fm_pools
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct, mo_struct
      TYPE(cp_fm_type)                                   :: fmdummy, fmlocal, rksmat, rsmat
      TYPE(cp_fm_type), POINTER                          :: rmos
      TYPE(dbcsr_type), POINTER                          :: cmatrix, rmatrix, tmpmat
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_matrix_pools_type), POINTER                :: mpools

      CALL timeset(routineN, handle)

      NULLIFY (cell_to_index, col_indices, eigenvalues, kp_dist, mpools, &
               para_env, row_indices, rmos, sab_nl, xkp)
      max_offdiag = 0.0_dp
      max_eval_deviation = 0.0_dp
      max_residual = 0.0_dp
      max_tangent_deviation = 0.0_dp
      max_gradient_norm = 0.0_dp
      max_descent_deviation = 0.0_dp

      CALL get_kpoint_info(kpoints, nkp=nkp, xkp=xkp, use_real_wfn=use_real_wfn, kp_range=kp_range, &
                           nkp_groups=nkp_groups, kp_dist=kp_dist, sab_nl=sab_nl, &
                           cell_to_index=cell_to_index, mpools=mpools)
      CPASSERT(ASSOCIATED(sab_nl))
      kplocal = kp_range(2) - kp_range(1) + 1
      nspin = SIZE(matrix_ks, 1)

      ALLOCATE (rmatrix, cmatrix, tmpmat)
      CALL dbcsr_create(rmatrix, template=matrix_ks(1, 1)%matrix, matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(cmatrix, template=matrix_ks(1, 1)%matrix, matrix_type=dbcsr_type_antisymmetric)
      CALL dbcsr_create(tmpmat, template=matrix_ks(1, 1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL cp_dbcsr_alloc_block_from_nbl(rmatrix, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(cmatrix, sab_nl)

      CALL mpools_get(mpools, ao_ao_fm_pools=ao_ao_fm_pools)
      CALL fm_pool_create_fm(ao_ao_fm_pools(1)%pool, fmlocal)
      CALL cp_fm_get_info(fmlocal, matrix_struct=matrix_struct)
      CALL cp_cfm_create(cksmat, matrix_struct)
      CALL cp_cfm_create(csmat, matrix_struct)
      IF (use_real_wfn) THEN
         CALL cp_fm_create(rksmat, matrix_struct)
         CALL cp_fm_create(rsmat, matrix_struct)
      END IF

      para_env => kpoints%blacs_env_all%para_env
      ALLOCATE (info(kplocal*nspin*nkp_groups, 4))

      indx = 0
      DO ikp = 1, kplocal
         DO ispin = 1, nspin
            DO igroup = 1, nkp_groups
               ik = kp_dist(1, igroup) + ikp - 1
               my_kpgrp = (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2))
               indx = indx + 1
               IF (use_real_wfn) THEN
                  CALL dbcsr_set(rmatrix, 0.0_dp)
                  CALL rskp_transform(rmatrix=rmatrix, rsmat=matrix_ks, ispin=ispin, &
                                      xkp=xkp(1:3, ik), cell_to_index=cell_to_index, sab_nl=sab_nl)
                  CALL dbcsr_desymmetrize(rmatrix, tmpmat)
                  CALL copy_dbcsr_to_fm(tmpmat, fmwork(1))
                  CALL dbcsr_set(rmatrix, 0.0_dp)
                  CALL rskp_transform(rmatrix=rmatrix, rsmat=matrix_s, ispin=1, &
                                      xkp=xkp(1:3, ik), cell_to_index=cell_to_index, sab_nl=sab_nl)
                  CALL dbcsr_desymmetrize(rmatrix, tmpmat)
                  CALL copy_dbcsr_to_fm(tmpmat, fmwork(3))
                  IF (my_kpgrp) THEN
                     CALL cp_fm_start_copy_general(fmwork(1), rksmat, para_env, info(indx, 1))
                     CALL cp_fm_start_copy_general(fmwork(3), rsmat, para_env, info(indx, 3))
                  ELSE
                     CALL cp_fm_start_copy_general(fmwork(1), fmdummy, para_env, info(indx, 1))
                     CALL cp_fm_start_copy_general(fmwork(3), fmdummy, para_env, info(indx, 3))
                  END IF
               ELSE
                  CALL dbcsr_set(rmatrix, 0.0_dp)
                  CALL dbcsr_set(cmatrix, 0.0_dp)
                  CALL rskp_transform(rmatrix=rmatrix, cmatrix=cmatrix, rsmat=matrix_ks, ispin=ispin, &
                                      xkp=xkp(1:3, ik), cell_to_index=cell_to_index, sab_nl=sab_nl)
                  CALL dbcsr_desymmetrize(rmatrix, tmpmat)
                  CALL copy_dbcsr_to_fm(tmpmat, fmwork(1))
                  CALL dbcsr_desymmetrize(cmatrix, tmpmat)
                  CALL copy_dbcsr_to_fm(tmpmat, fmwork(2))
                  CALL dbcsr_set(rmatrix, 0.0_dp)
                  CALL dbcsr_set(cmatrix, 0.0_dp)
                  CALL rskp_transform(rmatrix=rmatrix, cmatrix=cmatrix, rsmat=matrix_s, ispin=1, &
                                      xkp=xkp(1:3, ik), cell_to_index=cell_to_index, sab_nl=sab_nl)
                  CALL dbcsr_desymmetrize(rmatrix, tmpmat)
                  CALL copy_dbcsr_to_fm(tmpmat, fmwork(3))
                  CALL dbcsr_desymmetrize(cmatrix, tmpmat)
                  CALL copy_dbcsr_to_fm(tmpmat, fmwork(4))
                  IF (my_kpgrp) THEN
                     CALL cp_fm_start_copy_general(fmwork(1), fmlocal, para_env, info(indx, 1))
                     CALL cp_fm_start_copy_general(fmwork(2), fmlocal, para_env, info(indx, 2))
                     CALL cp_fm_start_copy_general(fmwork(3), fmlocal, para_env, info(indx, 3))
                     CALL cp_fm_start_copy_general(fmwork(4), fmlocal, para_env, info(indx, 4))
                  ELSE
                     CALL cp_fm_start_copy_general(fmwork(1), fmdummy, para_env, info(indx, 1))
                     CALL cp_fm_start_copy_general(fmwork(2), fmdummy, para_env, info(indx, 2))
                     CALL cp_fm_start_copy_general(fmwork(3), fmdummy, para_env, info(indx, 3))
                     CALL cp_fm_start_copy_general(fmwork(4), fmdummy, para_env, info(indx, 4))
                  END IF
               END IF
            END DO
         END DO
      END DO

      indx = 0
      DO ikp = 1, kplocal
         DO ispin = 1, nspin
            DO igroup = 1, nkp_groups
               ik = kp_dist(1, igroup) + ikp - 1
               my_kpgrp = (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2))
               indx = indx + 1
               IF (my_kpgrp) THEN
                  IF (use_real_wfn) THEN
                     CALL cp_fm_finish_copy_general(rksmat, info(indx, 1))
                     CALL cp_fm_to_cfm(msourcer=rksmat, mtarget=cksmat)
                     CALL cp_fm_finish_copy_general(rsmat, info(indx, 3))
                     CALL cp_fm_to_cfm(msourcer=rsmat, mtarget=csmat)
                  ELSE
                     CALL cp_fm_finish_copy_general(fmlocal, info(indx, 1))
                     CALL cp_cfm_scale_and_add_fm(z_zero, cksmat, z_one, fmlocal)
                     CALL cp_fm_finish_copy_general(fmlocal, info(indx, 2))
                     CALL cp_cfm_scale_and_add_fm(z_one, cksmat, gaussi, fmlocal)
                     CALL cp_fm_finish_copy_general(fmlocal, info(indx, 3))
                     CALL cp_cfm_scale_and_add_fm(z_zero, csmat, z_one, fmlocal)
                     CALL cp_fm_finish_copy_general(fmlocal, info(indx, 4))
                     CALL cp_cfm_scale_and_add_fm(z_one, csmat, gaussi, fmlocal)
                  END IF
               END IF
            END DO

            kp => kpoints%kp_env(ikp)%kpoint_env
            CALL get_mo_set(kp%mos(1, ispin), mo_coeff=rmos, eigenvalues=eigenvalues)
            CALL cp_fm_get_info(rmos, matrix_struct=mo_struct, nrow_global=nao, ncol_global=nmo)
            CALL cp_cfm_create(cmos, mo_struct)
            CALL cp_cfm_create(hcmos, mo_struct)
            CALL cp_cfm_create(scmos, mo_struct)
            CALL cp_cfm_create(hsmall, mo_struct, nrow=nmo, ncol=nmo, set_zero=.TRUE.)
            CALL cp_cfm_create(gsmall, mo_struct, nrow=nmo, ncol=nmo, set_zero=.TRUE.)
            CALL cp_cfm_create(dxmos, mo_struct)
            CALL ot_kp_pack_mos(kp, ispin, use_real_wfn, cmos)
            CALL cp_cfm_gemm('N', 'N', nao, nmo, nao, cone, cksmat, cmos, czero, hcmos)
            CALL cp_cfm_gemm('N', 'N', nao, nmo, nao, cone, csmat, cmos, czero, scmos)
            CALL cp_cfm_gemm('C', 'N', nmo, nmo, nao, cone, cmos, hcmos, czero, hsmall)

            CALL cp_cfm_get_info(hsmall, nrow_local=nrow_local, ncol_local=ncol_local, &
                                 row_indices=row_indices, col_indices=col_indices, &
                                 local_data=h_local, para_env=para_env)
            offdiag_local = 0.0_dp
            eval_dev_local = 0.0_dp
            DO icol = 1, ncol_local
               DO irow = 1, nrow_local
                  IF (row_indices(irow) == col_indices(icol)) THEN
                     eval_dev_local = MAX(eval_dev_local, &
                                          ABS(h_local(irow, icol) - eigenvalues(row_indices(irow))))
                  ELSE
                     offdiag_local = MAX(offdiag_local, ABS(h_local(irow, icol)))
                  END IF
               END DO
            END DO
            ! Use the small projected Hamiltonian to form a Galerkin residual,
            ! R = H C - S C (C^H H C).  This is the quantity that will become
            ! the input for the k-point OT gradient in the next implementation step.
            CALL cp_cfm_gemm('N', 'N', nao, nmo, nmo, -cone, scmos, hsmall, cone, hcmos)
            CALL cp_cfm_get_info(hcmos, nrow_local=nrow_local, ncol_local=ncol_local, &
                                 local_data=r_local, para_env=para_env)
            residual_local = 0.0_dp
            gradient_norm_sq_local = 0.0_dp
            DO icol = 1, ncol_local
               DO irow = 1, nrow_local
                  residual_local = MAX(residual_local, ABS(r_local(irow, icol)))
                  gradient_norm_sq_local = gradient_norm_sq_local + ABS(r_local(irow, icol))**2
               END DO
            END DO

            ! Treat the residual as the first unpreconditioned tangent gradient candidate.
            ! It should satisfy C^H G = 0 for the diagonalization-based initial state.
            CALL cp_cfm_gemm('C', 'N', nmo, nmo, nao, cone, cmos, hcmos, czero, gsmall)
            CALL cp_cfm_get_info(gsmall, nrow_local=nrow_local, ncol_local=ncol_local, &
                                 local_data=g_local, para_env=para_env)
            tangent_local = 0.0_dp
            DO icol = 1, ncol_local
               DO irow = 1, nrow_local
                  tangent_local = MAX(tangent_local, ABS(g_local(irow, icol)))
               END DO
            END DO

            ! The first steepest-descent direction would be D = -G.  Check the complex
            ! inner-product convention that will later be used by CG and line search.
            CALL cp_cfm_to_cfm(hcmos, dxmos)
            CALL cp_cfm_scale_and_add(-cone, dxmos)
            CALL cp_cfm_get_info(dxmos, nrow_local=nrow_local, ncol_local=ncol_local, &
                                 local_data=d_local, para_env=para_env)
            descent_local = 0.0_dp
            DO icol = 1, ncol_local
               DO irow = 1, nrow_local
                  descent_local = descent_local + REAL(CONJG(r_local(irow, icol))*d_local(irow, icol), KIND=dp)
               END DO
            END DO

            CALL para_env%max(offdiag_local)
            CALL para_env%max(eval_dev_local)
            CALL para_env%max(residual_local)
            CALL para_env%max(tangent_local)
            CALL para_env%sum(gradient_norm_sq_local)
            CALL para_env%sum(descent_local)
            gradient_norm_local = SQRT(MAX(gradient_norm_sq_local, 0.0_dp))
            descent_dev_local = ABS(descent_local + gradient_norm_sq_local)
            max_offdiag = MAX(max_offdiag, offdiag_local)
            max_eval_deviation = MAX(max_eval_deviation, eval_dev_local)
            max_residual = MAX(max_residual, residual_local)
            max_tangent_deviation = MAX(max_tangent_deviation, tangent_local)
            max_gradient_norm = MAX(max_gradient_norm, gradient_norm_local)
            max_descent_deviation = MAX(max_descent_deviation, descent_dev_local)

            CALL cp_cfm_release(dxmos)
            CALL cp_cfm_release(gsmall)
            CALL cp_cfm_release(hsmall)
            CALL cp_cfm_release(scmos)
            CALL cp_cfm_release(hcmos)
            CALL cp_cfm_release(cmos)
         END DO
      END DO

      indx = 0
      DO ikp = 1, kplocal
         DO ispin = 1, nspin
            DO igroup = 1, nkp_groups
               indx = indx + 1
               IF (use_real_wfn) THEN
                  CALL cp_fm_cleanup_copy_general(info(indx, 1))
                  CALL cp_fm_cleanup_copy_general(info(indx, 3))
               ELSE
                  CALL cp_fm_cleanup_copy_general(info(indx, 1))
                  CALL cp_fm_cleanup_copy_general(info(indx, 2))
                  CALL cp_fm_cleanup_copy_general(info(indx, 3))
                  CALL cp_fm_cleanup_copy_general(info(indx, 4))
               END IF
            END DO
         END DO
      END DO
      DEALLOCATE (info)

      IF (ASSOCIATED(kpoints%para_env_inter_kp)) THEN
         CALL kpoints%para_env_inter_kp%max(max_offdiag)
         CALL kpoints%para_env_inter_kp%max(max_eval_deviation)
         CALL kpoints%para_env_inter_kp%max(max_residual)
         CALL kpoints%para_env_inter_kp%max(max_tangent_deviation)
         CALL kpoints%para_env_inter_kp%max(max_gradient_norm)
         CALL kpoints%para_env_inter_kp%max(max_descent_deviation)
      END IF

      CALL dbcsr_deallocate_matrix(rmatrix)
      CALL dbcsr_deallocate_matrix(cmatrix)
      CALL dbcsr_deallocate_matrix(tmpmat)
      IF (use_real_wfn) THEN
         CALL cp_fm_release(rksmat)
         CALL cp_fm_release(rsmat)
      END IF
      CALL cp_cfm_release(csmat)
      CALL cp_cfm_release(cksmat)
      CALL fm_pool_give_back_fm(ao_ao_fm_pools(1)%pool, fmlocal)

      CALL timestop(handle)

   END SUBROUTINE ot_kp_check_hamiltonian


END MODULE qs_ot_kp
