!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Pair-local LRI backend for k-point HFX.
!> \author SY Wang (05.2026)
! **************************************************************************************************
MODULE hfx_lri_kp
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE hfx_lri_types,                   ONLY: hfx_lri_kp_release
   USE hfx_types,                       ONLY: hfx_ri_type
   USE input_constants,                 ONLY: do_potential_short
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_neighbor_lists,               ONLY: setup_neighbor_list
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hfx_lri_update_ks_kp

CONTAINS

! **************************************************************************************************
!> \brief Updates the k-point KS matrices with the pair-local LRI HFX contribution.
!> \param qs_env ...
!> \param ri_data ...
!> \param ks_matrix ...
!> \param ehfx ...
!> \param rho_ao ...
!> \param geometry_did_change ...
!> \param nspins ...
!> \param hf_fraction ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_update_ks_kp(qs_env, ri_data, ks_matrix, ehfx, rho_ao, &
                                   geometry_did_change, nspins, hf_fraction)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      REAL(KIND=dp), INTENT(OUT)                         :: ehfx
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      LOGICAL, INTENT(IN)                                :: geometry_did_change
      INTEGER, INTENT(IN)                                :: nspins
      REAL(KIND=dp), INTENT(IN)                          :: hf_fraction

      ehfx = 0.0_dp

      IF (ri_data%hfx_pot%potential_type /= do_potential_short) THEN
         CPABORT("LRI-HFXk is currently restricted to short-range/HSE exchange")
      END IF

      ! Keep dummy arguments referenced in this scaffold.
      IF (ASSOCIATED(qs_env) .AND. ASSOCIATED(ks_matrix) .AND. ASSOCIATED(rho_ao) .AND. &
          geometry_did_change .AND. nspins > 0 .AND. hf_fraction >= 0.0_dp) THEN
         CONTINUE
      END IF

      IF (geometry_did_change .OR. .NOT. ASSOCIATED(ri_data%lri_kp)) THEN
         CALL hfx_lri_pre_scf_kp(qs_env, ri_data)
      END IF

      CPABORT("LRI-HFXk pair-local contraction is not implemented yet")

   END SUBROUTINE hfx_lri_update_ks_kp

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_pre_scf_kp(qs_env, ri_data)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_lri_pre_scf_kp'

      INTEGER                                            :: handle, iatom, ikind, ilist, jatom, &
                                                            jkind, jneighbor, mepos, nao_i, nao_j, &
                                                            naux_i, naux_j, naux_pair, nkind, &
                                                            nlist, nneighbor
      INTEGER(int_8)                                     :: c_nreal
      REAL(dp), DIMENSION(3)                             :: rab
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_AO, basis_set_RI
      TYPE(gto_basis_set_type), POINTER                  :: obas_i, obas_j, rbas_i, rbas_j
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: soo_list
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (qs_kind_set, soo_list, nl_iterator)
      NULLIFY (obas_i, obas_j, rbas_i, rbas_j)
      NULLIFY (basis_set_AO, basis_set_RI)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, nkind=nkind)

      CALL hfx_lri_kp_release(ri_data%lri_kp)
      ALLOCATE (ri_data%lri_kp)

      ALLOCATE (basis_set_AO(nkind), basis_set_RI(nkind))
      CALL basis_set_list_setup(basis_set_AO, ri_data%orb_basis_type, qs_kind_set)
      CALL basis_set_list_setup(basis_set_RI, ri_data%ri_basis_type, qs_kind_set)

      ! First prototype: AO-AO neighbor list only.
      ! Later this cutoff should be tightened for HSE short-range LRI.
      CALL setup_neighbor_list(soo_list, basis_set_AO, qs_env=qs_env, mic=.FALSE., symmetric=.TRUE.)

      CALL neighbor_list_iterator_create(nl_iterator, soo_list)

      mepos = 0
      DO WHILE (neighbor_list_iterate(nl_iterator, mepos) == 0)

         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, &
                                nlist=nlist, ilist=ilist, nnode=nneighbor, inode=jneighbor, &
                                iatom=iatom, jatom=jatom, r=rab)

         obas_i => basis_set_AO(ikind)%gto_basis_set
         obas_j => basis_set_AO(jkind)%gto_basis_set
         rbas_i => basis_set_RI(ikind)%gto_basis_set
         rbas_j => basis_set_RI(jkind)%gto_basis_set

         IF (.NOT. ASSOCIATED(obas_i)) CYCLE
         IF (.NOT. ASSOCIATED(obas_j)) CYCLE
         IF (.NOT. ASSOCIATED(rbas_i)) CYCLE
         IF (.NOT. ASSOCIATED(rbas_j)) CYCLE

         nao_i = obas_i%nsgf
         nao_j = obas_j%nsgf
         naux_i = rbas_i%nsgf
         naux_j = rbas_j%nsgf

         ! Pair-local / two-center auxiliary space:
         ! Aux(ij) = RI(i) union RI(j).
         IF (iatom == jatom .AND. SUM(rab*rab) < 1.0E-20_dp) THEN
            naux_pair = naux_i
         ELSE
            naux_pair = naux_i + naux_j
         END IF

         c_nreal = INT(nao_i, int_8)*INT(nao_j, int_8)*INT(naux_pair, int_8)

         ri_data%lri_kp%npairs = ri_data%lri_kp%npairs + 1
         ri_data%lri_kp%sum_nao_pair = ri_data%lri_kp%sum_nao_pair + INT(nao_i, int_8)*INT(nao_j, int_8)
         ri_data%lri_kp%sum_naux_pair = ri_data%lri_kp%sum_naux_pair + INT(naux_pair, int_8)
         ri_data%lri_kp%c_nreal_est = ri_data%lri_kp%c_nreal_est + c_nreal

      END DO

      CALL neighbor_list_iterator_release(nl_iterator)
      CALL release_neighbor_list_sets(soo_list)

      IF (ASSOCIATED(basis_set_AO)) DEALLOCATE (basis_set_AO)
      IF (ASSOCIATED(basis_set_RI)) DEALLOCATE (basis_set_RI)

      ri_data%lri_kp%c_mem_mib_est = &
         REAL(ri_data%lri_kp%c_nreal_est, KIND=dp)*REAL(dp, KIND=dp)/1024.0_dp**2

      ri_data%lri_kp%initialized = .TRUE.

      IF (ri_data%unit_nr > 0) THEN
         WRITE (ri_data%unit_nr, "(/,T2,A)") "LRI-HFXk pair-local setup"
         WRITE (ri_data%unit_nr, "(T4,A,T46,I12)") "number of AO atom pairs", ri_data%lri_kp%npairs
         WRITE (ri_data%unit_nr, "(T4,A,T46,I12)") "sum AO-pair block size", ri_data%lri_kp%sum_nao_pair
         WRITE (ri_data%unit_nr, "(T4,A,T46,I12)") "sum pair-local aux size", ri_data%lri_kp%sum_naux_pair
         WRITE (ri_data%unit_nr, "(T4,A,T46,ES12.4)") "estimated C storage / MiB", ri_data%lri_kp%c_mem_mib_est
      END IF

      CALL timestop(handle)

   END SUBROUTINE hfx_lri_pre_scf_kp

END MODULE hfx_lri_kp
