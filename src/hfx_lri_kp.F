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
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_block_p,&
                                              dbcsr_p_type
   USE gamma,                           ONLY: init_md_ftable
   USE hfx_lri_types,                   ONLY: hfx_lri_kp_release,&
                                              hfx_lri_pair_block_type
   USE hfx_types,                       ONLY: hfx_ri_type
   USE input_constants,                 ONLY: do_potential_short
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE libint_2c_3c,                    ONLY: eri_2center,&
                                              eri_3center
   USE libint_wrapper,                  ONLY: cp_libint_cleanup_2eri,&
                                              cp_libint_cleanup_3eri,&
                                              cp_libint_init_2eri,&
                                              cp_libint_init_3eri,&
                                              cp_libint_set_contrdepth,&
                                              cp_libint_t
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: ncoset
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
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

   INTEGER, PARAMETER :: lri_debug_max_pair_blocks = 1000
   INTEGER, PARAMETER :: lri_debug_max_k_blocks = 500
   REAL(dp), PARAMETER :: lri_debug_rho_threshold = 1.0E-8_dp

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

      CALL hfx_lri_debug_density_block(qs_env, ri_data, rho_ao)
      CALL hfx_lri_debug_ks_writeback_probe(qs_env, ri_data, ks_matrix, hf_fraction)

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

      INTEGER :: dbg_iatom, dbg_ikind, dbg_j_cols, dbg_j_rows, dbg_jatom, dbg_ji_ncart, &
         dbg_jj_ncart, dbg_jkind, dbg_m_dim, dbg_m_eig_small, dbg_m_ncart, dbg_nao_i, dbg_nao_j, &
         dbg_naux_i, dbg_naux_j, dbg_naux_pair, dbg_pair_pair_index, dbg_sph_j_cols, &
         dbg_sph_j_rows, dbg_sph_m_dim, dbg_stored_index, handle, iatom, ibasis, ikind, ilist, &
         ipair_block, jatom, jkind, jneighbor, mepos, nao_i, nao_j, naux_i, naux_j, naux_pair, &
         nkind, nlist, nneighbor, npair_blocks_local, npair_blocks_max, npair_blocks_sum, &
         npairs_local, npairs_max, npairs_sum
      INTEGER(int_8) :: c_nreal, c_nreal_local, c_nreal_max, c_nreal_stored_local, &
         c_nreal_stored_max, c_nreal_stored_sum, c_nreal_sum, sum_nao_pair_local, &
         sum_nao_pair_max, sum_nao_pair_sum, sum_naux_pair_local, sum_naux_pair_max, &
         sum_naux_pair_sum
      INTEGER, DIMENSION(3)                              :: cell_j, dbg_cell
      LOGICAL                                            :: dbg_integrals_done, dbg_pair_found
      REAL(dp) :: c_mem_mib_local, c_mem_mib_max, c_mem_mib_stored_local, c_mem_mib_stored_max, &
         c_mem_mib_stored_sum, c_mem_mib_sum, dbg_c_norm, dbg_j_norm, dbg_ji_norm, dbg_jj_norm, &
         dbg_m_cond_diag, dbg_m_dense_norm, dbg_m_eig_cond, dbg_m_eig_max, dbg_m_eig_min, &
         dbg_m_norm, dbg_pair_pair_cvc12_norm, dbg_pair_pair_r2, dbg_pair_pair_v12_norm, dbg_r2, &
         dbg_sph_c_norm, dbg_sph_cvc_norm, dbg_sph_j_norm, dbg_sph_m_norm, dbg_stored_c_norm
      REAL(dp), DIMENSION(3)                             :: dbg_rab, rab
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_AO, basis_set_RI
      TYPE(gto_basis_set_type), POINTER                  :: obas_i, obas_j, rbas_i, rbas_j
      TYPE(hfx_lri_pair_block_type)                      :: dbg_pair_block
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: soo_list
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (qs_kind_set, soo_list, nl_iterator, para_env, dft_control)
      NULLIFY (obas_i, obas_j, rbas_i, rbas_j)
      NULLIFY (basis_set_AO, basis_set_RI)

      IF (ALLOCATED(dbg_pair_block%c_fit)) DEALLOCATE (dbg_pair_block%c_fit)
      dbg_pair_found = .FALSE.
      dbg_iatom = 0
      dbg_jatom = 0
      dbg_ikind = 0
      dbg_jkind = 0
      dbg_nao_i = 0
      dbg_nao_j = 0
      dbg_naux_i = 0
      dbg_naux_j = 0
      dbg_naux_pair = 0
      dbg_cell(:) = 0
      dbg_r2 = HUGE(1.0_dp)
      dbg_rab(:) = 0.0_dp
      dbg_integrals_done = .FALSE.
      dbg_m_ncart = 0
      dbg_ji_ncart = 0
      dbg_jj_ncart = 0
      dbg_m_dim = 0
      dbg_j_rows = 0
      dbg_j_cols = 0
      dbg_sph_m_dim = 0
      dbg_sph_j_rows = 0
      dbg_sph_j_cols = 0
      dbg_m_norm = 0.0_dp
      dbg_ji_norm = 0.0_dp
      dbg_jj_norm = 0.0_dp
      dbg_m_dense_norm = 0.0_dp
      dbg_j_norm = 0.0_dp
      dbg_m_cond_diag = 0.0_dp
      dbg_m_eig_min = 0.0_dp
      dbg_m_eig_max = 0.0_dp
      dbg_m_eig_cond = 0.0_dp
      dbg_m_eig_small = 0
      dbg_c_norm = 0.0_dp
      dbg_sph_m_norm = 0.0_dp
      dbg_sph_j_norm = 0.0_dp
      dbg_sph_c_norm = 0.0_dp
      dbg_sph_cvc_norm = 0.0_dp
      dbg_stored_index = 0
      dbg_stored_c_norm = 0.0_dp
      dbg_pair_pair_index = 0
      dbg_pair_pair_r2 = HUGE(1.0_dp)
      dbg_pair_pair_v12_norm = 0.0_dp
      dbg_pair_pair_cvc12_norm = 0.0_dp

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, nkind=nkind, para_env=para_env, &
                      dft_control=dft_control)

      CALL hfx_lri_kp_release(ri_data%lri_kp)
      ALLOCATE (ri_data%lri_kp)
      ALLOCATE (ri_data%lri_kp%pair_blocks(lri_debug_max_pair_blocks))

      ALLOCATE (basis_set_AO(nkind), basis_set_RI(nkind))
      CALL basis_set_list_setup(basis_set_AO, ri_data%orb_basis_type, qs_kind_set)
      CALL basis_set_list_setup(basis_set_RI, ri_data%ri_basis_type, qs_kind_set)

      ! Use RI-HFX screening radii while building the LRI-HFXk pair list.
      ! They are restored to the standard QS radii before leaving this routine.
      DO ibasis = 1, SIZE(basis_set_AO)
         obas_i => basis_set_AO(ibasis)%gto_basis_set
         rbas_i => basis_set_RI(ibasis)%gto_basis_set
         IF (ASSOCIATED(obas_i)) THEN
            CALL init_interaction_radii_orb_basis(obas_i, ri_data%eps_pgf_orb)
         END IF
         IF (ASSOCIATED(rbas_i)) THEN
            CALL init_interaction_radii_orb_basis(rbas_i, ri_data%eps_pgf_orb)
         END IF
      END DO

      ! First prototype: AO-AO neighbor list only.
      ! Later this cutoff should be tightened for HSE short-range LRI.
      CALL setup_neighbor_list(soo_list, basis_set_AO, qs_env=qs_env, mic=.FALSE., symmetric=.TRUE.)

      CALL neighbor_list_iterator_create(nl_iterator, soo_list)

      mepos = 0
      DO WHILE (neighbor_list_iterate(nl_iterator, mepos) == 0)

         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, &
                                nlist=nlist, ilist=ilist, nnode=nneighbor, inode=jneighbor, &
                                iatom=iatom, jatom=jatom, r=rab, cell=cell_j)

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

         IF ((.NOT. dbg_pair_found) .OR. SUM(rab*rab) < dbg_r2) THEN
            dbg_pair_found = .TRUE.
            dbg_r2 = SUM(rab*rab)
            dbg_iatom = iatom
            dbg_jatom = jatom
            dbg_ikind = ikind
            dbg_jkind = jkind
            dbg_nao_i = nao_i
            dbg_nao_j = nao_j
            dbg_naux_i = naux_i
            dbg_naux_j = naux_j
            dbg_naux_pair = naux_pair
            dbg_cell(:) = cell_j(:)
            dbg_rab(:) = rab(:)
         END IF

         c_nreal = INT(nao_i, int_8)*INT(nao_j, int_8)*INT(naux_pair, int_8)

         IF (ri_data%lri_kp%npair_blocks < lri_debug_max_pair_blocks) THEN
            ipair_block = ri_data%lri_kp%npair_blocks + 1
            CALL hfx_lri_init_pair_block(ri_data%lri_kp%pair_blocks(ipair_block), iatom, jatom, &
                                         ikind, jkind, cell_j, rab, nao_i, nao_j, naux_pair)
            CALL hfx_lri_debug_pair_integrals(basis_set_AO, basis_set_RI, ikind, jkind, rab, &
                                              ri_data, ri_data%lri_kp%pair_blocks(ipair_block), &
                                              dbg_m_ncart, dbg_ji_ncart, dbg_jj_ncart, &
                                              dbg_m_norm, dbg_ji_norm, dbg_jj_norm, dbg_m_dim, &
                                              dbg_j_rows, dbg_j_cols, dbg_m_dense_norm, &
                                              dbg_j_norm, dbg_m_cond_diag, dbg_m_eig_min, &
                                              dbg_m_eig_max, dbg_m_eig_cond, dbg_m_eig_small, &
                                              dbg_c_norm, dbg_sph_m_dim, dbg_sph_j_rows, &
                                              dbg_sph_j_cols, dbg_sph_m_norm, dbg_sph_j_norm, &
                                              dbg_sph_c_norm, dbg_sph_cvc_norm, debug_probe=.FALSE.)
            IF (ALLOCATED(ri_data%lri_kp%pair_blocks(ipair_block)%c_fit)) THEN
               ri_data%lri_kp%npair_blocks = ipair_block
               ri_data%lri_kp%c_nreal_stored = ri_data%lri_kp%c_nreal_stored + &
                                               INT(SIZE(ri_data%lri_kp%pair_blocks(ipair_block)%c_fit), int_8)
            END IF
         END IF

         ri_data%lri_kp%npairs = ri_data%lri_kp%npairs + 1
         ri_data%lri_kp%sum_nao_pair = ri_data%lri_kp%sum_nao_pair + INT(nao_i, int_8)*INT(nao_j, int_8)
         ri_data%lri_kp%sum_naux_pair = ri_data%lri_kp%sum_naux_pair + INT(naux_pair, int_8)
         ri_data%lri_kp%c_nreal_est = ri_data%lri_kp%c_nreal_est + c_nreal

      END DO

      CALL neighbor_list_iterator_release(nl_iterator)
      CALL release_neighbor_list_sets(soo_list)

      IF (dbg_pair_found) THEN
         CALL hfx_lri_init_pair_block(dbg_pair_block, dbg_iatom, dbg_jatom, dbg_ikind, dbg_jkind, &
                                      dbg_cell, dbg_rab, dbg_nao_i, dbg_nao_j, dbg_naux_pair)
         CALL hfx_lri_debug_pair_integrals(basis_set_AO, basis_set_RI, dbg_ikind, dbg_jkind, &
                                           dbg_rab, ri_data, dbg_pair_block, dbg_m_ncart, &
                                           dbg_ji_ncart, dbg_jj_ncart, dbg_m_norm, dbg_ji_norm, &
                                           dbg_jj_norm, dbg_m_dim, dbg_j_rows, dbg_j_cols, &
                                           dbg_m_dense_norm, dbg_j_norm, dbg_m_cond_diag, &
                                           dbg_m_eig_min, dbg_m_eig_max, dbg_m_eig_cond, &
                                           dbg_m_eig_small, dbg_c_norm, dbg_sph_m_dim, &
                                           dbg_sph_j_rows, dbg_sph_j_cols, dbg_sph_m_norm, &
                                           dbg_sph_j_norm, dbg_sph_c_norm, dbg_sph_cvc_norm)
         dbg_integrals_done = .TRUE.
      END IF

      IF (dbg_pair_found .AND. ri_data%lri_kp%npair_blocks > 0) THEN
         dbg_stored_index = hfx_lri_find_pair_block(ri_data%lri_kp%pair_blocks, &
                                                    ri_data%lri_kp%npair_blocks, &
                                                    dbg_iatom, dbg_jatom, dbg_cell)
      END IF

      IF (dbg_stored_index > 0) THEN
         DO ipair_block = 1, ri_data%lri_kp%npair_blocks
            IF (ipair_block == dbg_stored_index) CYCLE
            IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks(ipair_block)%c_fit)) CYCLE
            ! The temporary real V12 probe is currently limited to pair blocks
            ! sharing the same first center, so their auxiliary-center positions
            ! can be expressed with the stored rab vectors alone. Pick the nearest
            ! such pair rather than the first one in the neighbor-list order.
            IF (ri_data%lri_kp%pair_blocks(ipair_block)%iatom /= &
                ri_data%lri_kp%pair_blocks(dbg_stored_index)%iatom) CYCLE
            IF (SUM((ri_data%lri_kp%pair_blocks(ipair_block)%rab - &
                     ri_data%lri_kp%pair_blocks(dbg_stored_index)%rab)**2) >= dbg_pair_pair_r2) CYCLE
            dbg_pair_pair_index = ipair_block
            dbg_pair_pair_r2 = SUM((ri_data%lri_kp%pair_blocks(ipair_block)%rab - &
                                    ri_data%lri_kp%pair_blocks(dbg_stored_index)%rab)**2)
         END DO
         IF (dbg_pair_pair_index > 0) THEN
            CALL hfx_lri_contract_pair_pair_cvc_realv(basis_set_RI, &
                                                      ri_data%lri_kp%pair_blocks(dbg_stored_index), &
                                                      ri_data%lri_kp%pair_blocks(dbg_pair_pair_index), &
                                                      ri_data, dbg_pair_pair_v12_norm, &
                                                      dbg_pair_pair_cvc12_norm)
         END IF
      END IF

      ! Restore the standard QS screening radii after the LRI-HFXk setup.
      IF (ASSOCIATED(dft_control)) THEN
         DO ibasis = 1, SIZE(basis_set_AO)
            obas_i => basis_set_AO(ibasis)%gto_basis_set
            rbas_i => basis_set_RI(ibasis)%gto_basis_set
            IF (ASSOCIATED(obas_i)) THEN
               CALL init_interaction_radii_orb_basis(obas_i, dft_control%qs_control%eps_pgf_orb)
            END IF
            IF (ASSOCIATED(rbas_i)) THEN
               CALL init_interaction_radii_orb_basis(rbas_i, dft_control%qs_control%eps_pgf_orb)
            END IF
         END DO
      END IF

      IF (ASSOCIATED(basis_set_AO)) DEALLOCATE (basis_set_AO)
      IF (ASSOCIATED(basis_set_RI)) DEALLOCATE (basis_set_RI)

      npairs_local = ri_data%lri_kp%npairs
      sum_nao_pair_local = ri_data%lri_kp%sum_nao_pair
      sum_naux_pair_local = ri_data%lri_kp%sum_naux_pair
      c_nreal_local = ri_data%lri_kp%c_nreal_est
      c_nreal_stored_local = ri_data%lri_kp%c_nreal_stored
      npair_blocks_local = ri_data%lri_kp%npair_blocks
      c_mem_mib_local = REAL(c_nreal_local, KIND=dp)*8.0_dp/1024.0_dp**2
      c_mem_mib_stored_local = REAL(c_nreal_stored_local, KIND=dp)*8.0_dp/1024.0_dp**2

      npairs_sum = npairs_local
      npairs_max = npairs_local
      sum_nao_pair_sum = sum_nao_pair_local
      sum_nao_pair_max = sum_nao_pair_local
      sum_naux_pair_sum = sum_naux_pair_local
      sum_naux_pair_max = sum_naux_pair_local
      c_nreal_sum = c_nreal_local
      c_nreal_max = c_nreal_local
      c_nreal_stored_sum = c_nreal_stored_local
      c_nreal_stored_max = c_nreal_stored_local
      npair_blocks_sum = npair_blocks_local
      npair_blocks_max = npair_blocks_local

      IF (ASSOCIATED(para_env)) THEN
         CALL para_env%sum(npairs_sum)
         CALL para_env%max(npairs_max)
         CALL para_env%sum(sum_nao_pair_sum)
         CALL para_env%max(sum_nao_pair_max)
         CALL para_env%sum(sum_naux_pair_sum)
         CALL para_env%max(sum_naux_pair_max)
         CALL para_env%sum(c_nreal_sum)
         CALL para_env%max(c_nreal_max)
         CALL para_env%sum(c_nreal_stored_sum)
         CALL para_env%max(c_nreal_stored_max)
         CALL para_env%sum(npair_blocks_sum)
         CALL para_env%max(npair_blocks_max)
      END IF

      c_mem_mib_sum = REAL(c_nreal_sum, KIND=dp)*8.0_dp/1024.0_dp**2
      c_mem_mib_max = REAL(c_nreal_max, KIND=dp)*8.0_dp/1024.0_dp**2
      c_mem_mib_stored_sum = REAL(c_nreal_stored_sum, KIND=dp)*8.0_dp/1024.0_dp**2
      c_mem_mib_stored_max = REAL(c_nreal_stored_max, KIND=dp)*8.0_dp/1024.0_dp**2

      ri_data%lri_kp%c_mem_mib_est = c_mem_mib_local
      ri_data%lri_kp%c_mem_mib_stored = c_mem_mib_stored_local
      ri_data%lri_kp%initialized = .TRUE.

      IF (ri_data%unit_nr > 0) THEN
         WRITE (ri_data%unit_nr, "(/,T2,A)") "LRI-HFXk pair-local setup"
         WRITE (ri_data%unit_nr, "(T4,A,T28,A)") "orbital basis", TRIM(ri_data%orb_basis_type)
         WRITE (ri_data%unit_nr, "(T4,A,T28,A)") "RI basis", TRIM(ri_data%ri_basis_type)
         WRITE (ri_data%unit_nr, "(/,T4,A,T30,A,T43,A,T56,A)") &
            "quantity", "local", "global", "max/rank"
         WRITE (ri_data%unit_nr, "(T4,A,T28,I10,T41,I10,T54,I10)") &
            "AO pairs", npairs_local, npairs_sum, npairs_max
         WRITE (ri_data%unit_nr, "(T4,A,T28,I10,T41,I10,T54,I10)") &
            "AO-pair size", sum_nao_pair_local, sum_nao_pair_sum, sum_nao_pair_max
         WRITE (ri_data%unit_nr, "(T4,A,T28,I10,T41,I10,T54,I10)") &
            "pair aux size", sum_naux_pair_local, sum_naux_pair_sum, sum_naux_pair_max
         WRITE (ri_data%unit_nr, "(T4,A,T28,F10.2,T41,F10.2,T54,F10.2)") &
            "C storage/MiB", c_mem_mib_local, c_mem_mib_sum, c_mem_mib_max
         WRITE (ri_data%unit_nr, "(T4,A,T28,I10,T41,I10,T54,I10)") &
            "stored blocks", npair_blocks_local, npair_blocks_sum, npair_blocks_max
         WRITE (ri_data%unit_nr, "(T4,A,T28,F10.4,T41,F10.4,T54,F10.4)") &
            "stored C/MiB", c_mem_mib_stored_local, c_mem_mib_stored_sum, c_mem_mib_stored_max
         IF (npairs_sum > 0) THEN
            WRITE (ri_data%unit_nr, "(T4,A,T28,F10.2,T41,F10.2,T54,F10.2)") &
               "avg AO-pair", &
               REAL(sum_nao_pair_local, KIND=dp)/MAX(1.0_dp, REAL(npairs_local, KIND=dp)), &
               REAL(sum_nao_pair_sum, KIND=dp)/REAL(npairs_sum, KIND=dp), &
               REAL(sum_nao_pair_max, KIND=dp)/MAX(1.0_dp, REAL(npairs_max, KIND=dp))
            WRITE (ri_data%unit_nr, "(T4,A,T28,F10.2,T41,F10.2,T54,F10.2)") &
               "avg pair aux", &
               REAL(sum_naux_pair_local, KIND=dp)/MAX(1.0_dp, REAL(npairs_local, KIND=dp)), &
               REAL(sum_naux_pair_sum, KIND=dp)/REAL(npairs_sum, KIND=dp), &
               REAL(sum_naux_pair_max, KIND=dp)/MAX(1.0_dp, REAL(npairs_max, KIND=dp))
         END IF
         IF (dbg_pair_found) THEN
            WRITE (ri_data%unit_nr, "(/,T4,A)") "debug local AO pair"
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
               "atoms", dbg_iatom, dbg_jatom
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
               "kinds", dbg_ikind, dbg_jkind
            WRITE (ri_data%unit_nr, "(T6,A,T28,3I6)") &
               "cell", dbg_cell
            WRITE (ri_data%unit_nr, "(T6,A,T28,F10.4)") &
               "distance", SQRT(dbg_r2)
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
               "AO functions", dbg_nao_i, dbg_nao_j
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10,T54,I10)") &
               "RI functions", dbg_naux_i, dbg_naux_j, dbg_naux_pair
            IF (dbg_integrals_done) THEN
               WRITE (ri_data%unit_nr, "(/,T4,A)") "debug spherical contracted blocks"
               WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,ES10.3)") &
                  "metric dim/norm", dbg_sph_m_dim, dbg_sph_m_norm
               WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
                  "J rows/cols", dbg_sph_j_rows, dbg_sph_j_cols
               WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3)") &
                  "J norm", dbg_sph_j_norm
               WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3)") &
                  "C=J Minv norm", dbg_sph_c_norm
               IF (ALLOCATED(dbg_pair_block%c_fit)) THEN
                  WRITE (ri_data%unit_nr, "(/,T4,A)") "debug pair coefficient block"
                  WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10,T54,I10)") &
                     "C dims", SIZE(dbg_pair_block%c_fit, 1), SIZE(dbg_pair_block%c_fit, 2), &
                     SIZE(dbg_pair_block%c_fit, 3)
                  WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3,T41,F10.4)") &
                     "C norm/MiB", SQRT(SUM(dbg_pair_block%c_fit*dbg_pair_block%c_fit)), &
                     REAL(SIZE(dbg_pair_block%c_fit), KIND=dp)*8.0_dp/1024.0_dp**2
               END IF
               WRITE (ri_data%unit_nr, "(/,T4,A)") "debug real pair-pair CVC probe"
               WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
                  "pair indices", dbg_stored_index, dbg_pair_pair_index
               IF (dbg_pair_pair_index > 0) THEN
                  WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3,T41,ES10.3)") &
                     "V12/CVC norm", dbg_pair_pair_v12_norm, dbg_pair_pair_cvc12_norm
               END IF
            END IF
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE hfx_lri_pre_scf_kp

! **************************************************************************************************
!> \brief ...
!> \param basis_set_AO ...
!> \param basis_set_RI ...
!> \param ikind ...
!> \param jkind ...
!> \param rab ...
!> \param ri_data ...
!> \param pair_block ...
!> \param metric_ncart ...
!> \param threec_i_ncart ...
!> \param threec_j_ncart ...
!> \param metric_norm ...
!> \param threec_i_norm ...
!> \param threec_j_norm ...
!> \param metric_dim ...
!> \param threec_rows ...
!> \param threec_cols ...
!> \param metric_dense_norm ...
!> \param threec_norm ...
!> \param metric_diag_cond ...
!> \param metric_eig_min ...
!> \param metric_eig_max ...
!> \param metric_eig_cond ...
!> \param metric_eig_small ...
!> \param c_norm ...
!> \param sph_metric_dim ...
!> \param sph_j_rows ...
!> \param sph_j_cols ...
!> \param sph_metric_norm ...
!> \param sph_j_norm ...
!> \param sph_c_norm ...
!> \param sph_cvc_norm ...
!> \param debug_probe ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_debug_pair_integrals(basis_set_AO, basis_set_RI, ikind, jkind, rab, ri_data, &
                                           pair_block, metric_ncart, threec_i_ncart, threec_j_ncart, &
                                           metric_norm, threec_i_norm, threec_j_norm, metric_dim, &
                                           threec_rows, threec_cols, metric_dense_norm, threec_norm, &
                                           metric_diag_cond, metric_eig_min, metric_eig_max, &
                                           metric_eig_cond, metric_eig_small, c_norm, sph_metric_dim, &
                                           sph_j_rows, sph_j_cols, sph_metric_norm, sph_j_norm, &
                                           sph_c_norm, sph_cvc_norm, debug_probe)

      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_AO, basis_set_RI
      INTEGER, INTENT(IN)                                :: ikind, jkind
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rab
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      TYPE(hfx_lri_pair_block_type), INTENT(INOUT)       :: pair_block
      INTEGER, INTENT(OUT)                               :: metric_ncart, threec_i_ncart, &
                                                            threec_j_ncart
      REAL(dp), INTENT(OUT)                              :: metric_norm, threec_i_norm, threec_j_norm
      INTEGER, INTENT(OUT)                               :: metric_dim, threec_rows, threec_cols
      REAL(dp), INTENT(OUT)                              :: metric_dense_norm, threec_norm, &
                                                            metric_diag_cond, metric_eig_min, &
                                                            metric_eig_max, metric_eig_cond
      INTEGER, INTENT(OUT)                               :: metric_eig_small
      REAL(dp), INTENT(OUT)                              :: c_norm
      INTEGER, INTENT(OUT)                               :: sph_metric_dim, sph_j_rows, sph_j_cols
      REAL(dp), INTENT(OUT)                              :: sph_metric_norm, sph_j_norm, sph_c_norm, &
                                                            sph_cvc_norm
      LOGICAL, INTENT(IN), OPTIONAL                      :: debug_probe

      INTEGER                                            :: idim, maxl_ai, maxl_aj, maxl_ri, &
                                                            maxl_rj, mmax, ncart_ai, ncart_aj, &
                                                            ncart_ri, ncart_rj, offset_j, &
                                                            sph_eig_small
      LOGICAL                                            :: do_debug_probe, onsite_aux
      REAL(dp) :: diag_abs, diag_max, diag_min, metric_norm2, sph_eig_cond, sph_eig_max, &
         sph_eig_min, threec_i_norm2, threec_j_norm2
      REAL(dp), ALLOCATABLE, DIMENSION(:, :) :: c_sph, cvc_block, j_block, j_sph, metric_block, &
         metric_sph, sphi_aux, sphi_i, sphi_j, sphi_pair, sphi_ri_i, sphi_ri_j, tmp_cvc, tmp_sph
      REAL(dp), DIMENSION(3)                             :: r0
      TYPE(cp_libint_t)                                  :: lib2, lib3
      TYPE(gto_basis_set_type), POINTER                  :: obas_i, obas_j, rbas_i, rbas_j

      NULLIFY (obas_i, obas_j, rbas_i, rbas_j)

      obas_i => basis_set_AO(ikind)%gto_basis_set
      obas_j => basis_set_AO(jkind)%gto_basis_set
      rbas_i => basis_set_RI(ikind)%gto_basis_set
      rbas_j => basis_set_RI(jkind)%gto_basis_set

      metric_ncart = 0
      threec_i_ncart = 0
      threec_j_ncart = 0
      metric_dim = 0
      threec_rows = 0
      threec_cols = 0
      metric_norm = 0.0_dp
      threec_i_norm = 0.0_dp
      threec_j_norm = 0.0_dp
      metric_dense_norm = 0.0_dp
      threec_norm = 0.0_dp
      metric_diag_cond = 0.0_dp
      metric_eig_min = 0.0_dp
      metric_eig_max = 0.0_dp
      metric_eig_cond = 0.0_dp
      metric_eig_small = 0
      c_norm = 0.0_dp
      sph_metric_dim = 0
      sph_j_rows = 0
      sph_j_cols = 0
      sph_metric_norm = 0.0_dp
      sph_j_norm = 0.0_dp
      sph_c_norm = 0.0_dp

      do_debug_probe = .TRUE.
      IF (PRESENT(debug_probe)) do_debug_probe = debug_probe

      IF (.NOT. ASSOCIATED(obas_i)) RETURN
      IF (.NOT. ASSOCIATED(obas_j)) RETURN
      IF (.NOT. ASSOCIATED(rbas_i)) RETURN
      IF (.NOT. ASSOCIATED(rbas_j)) RETURN

      CALL get_gto_basis_set(obas_i, maxl=maxl_ai)
      CALL get_gto_basis_set(obas_j, maxl=maxl_aj)
      CALL get_gto_basis_set(rbas_i, maxl=maxl_ri)
      CALL get_gto_basis_set(rbas_j, maxl=maxl_rj)

      mmax = MAX(maxl_ri + maxl_rj, maxl_ai + maxl_aj + MAX(maxl_ri, maxl_rj))
      CALL init_md_ftable(nmax=mmax)

      CALL cp_libint_init_2eri(lib2, MAX(maxl_ri, maxl_rj))
      CALL cp_libint_set_contrdepth(lib2, 1)
      CALL cp_libint_init_3eri(lib3, MAX(maxl_ai, maxl_aj, maxl_ri, maxl_rj))
      CALL cp_libint_set_contrdepth(lib3, 1)

      r0(:) = 0.0_dp
      metric_norm2 = 0.0_dp
      threec_i_norm2 = 0.0_dp
      threec_j_norm2 = 0.0_dp
      sph_eig_min = 0.0_dp
      sph_eig_max = 0.0_dp
      sph_eig_cond = 0.0_dp
      sph_eig_small = 0
      onsite_aux = (NORM2(rab) <= 1.0E-12_dp .AND. ikind == jkind)

      IF (do_debug_probe) THEN
         CALL hfx_lri_accumulate_2c_norm(rbas_i, rbas_i, r0, r0, lib2, ri_data, metric_ncart, &
                                         metric_norm2)
         IF (.NOT. onsite_aux) THEN
            CALL hfx_lri_accumulate_2c_norm(rbas_i, rbas_j, r0, rab, lib2, ri_data, metric_ncart, &
                                            metric_norm2, prefactor=2.0_dp)
            CALL hfx_lri_accumulate_2c_norm(rbas_j, rbas_j, rab, rab, lib2, ri_data, metric_ncart, &
                                            metric_norm2)
         END IF

         CALL hfx_lri_accumulate_3c_norm(obas_i, obas_j, rbas_i, r0, rab, r0, lib3, ri_data, &
                                         threec_i_ncart, threec_i_norm2)
         IF (.NOT. onsite_aux) THEN
            CALL hfx_lri_accumulate_3c_norm(obas_i, obas_j, rbas_j, r0, rab, rab, lib3, ri_data, &
                                            threec_j_ncart, threec_j_norm2)
         END IF
      END IF

      ncart_ai = hfx_lri_count_cart_functions(obas_i)
      ncart_aj = hfx_lri_count_cart_functions(obas_j)
      ncart_ri = hfx_lri_count_cart_functions(rbas_i)
      ncart_rj = hfx_lri_count_cart_functions(rbas_j)
      metric_dim = ncart_ri
      IF (.NOT. onsite_aux) metric_dim = metric_dim + ncart_rj
      threec_rows = ncart_ai*ncart_aj
      threec_cols = metric_dim

      IF (metric_dim > 0 .AND. threec_rows > 0) THEN
         ALLOCATE (metric_block(metric_dim, metric_dim), SOURCE=0.0_dp)
         ALLOCATE (j_block(threec_rows, threec_cols), SOURCE=0.0_dp)

         CALL hfx_lri_fill_2c_block(rbas_i, rbas_i, r0, r0, lib2, ri_data, metric_block, 0, 0)
         CALL hfx_lri_fill_3c_block(obas_i, obas_j, rbas_i, r0, rab, r0, lib3, ri_data, &
                                    j_block, 0, ncart_aj)
         IF (.NOT. onsite_aux) THEN
            offset_j = ncart_ri
            CALL hfx_lri_fill_2c_block(rbas_i, rbas_j, r0, rab, lib2, ri_data, metric_block, 0, &
                                       offset_j, symmetric=.TRUE.)
            CALL hfx_lri_fill_2c_block(rbas_j, rbas_j, rab, rab, lib2, ri_data, metric_block, &
                                       offset_j, offset_j)
            CALL hfx_lri_fill_3c_block(obas_i, obas_j, rbas_j, r0, rab, rab, lib3, ri_data, &
                                       j_block, offset_j, ncart_aj)
         END IF

         IF (do_debug_probe) THEN
            metric_dense_norm = SQRT(SUM(metric_block*metric_block))
            threec_norm = SQRT(SUM(j_block*j_block))
            diag_min = HUGE(1.0_dp)
            diag_max = 0.0_dp
            DO idim = 1, metric_dim
               diag_abs = ABS(metric_block(idim, idim))
               diag_max = MAX(diag_max, diag_abs)
               IF (diag_abs > 1.0E-300_dp) diag_min = MIN(diag_min, diag_abs)
            END DO
            IF (diag_min < HUGE(1.0_dp)) metric_diag_cond = diag_max/diag_min

            CALL hfx_lri_probe_metric_inverse(metric_block, j_block, metric_eig_min, metric_eig_max, &
                                              metric_eig_cond, metric_eig_small, c_norm)
         END IF

         CALL hfx_lri_build_sphi_matrix(obas_i, sphi_i)
         CALL hfx_lri_build_sphi_matrix(obas_j, sphi_j)
         CALL hfx_lri_build_sphi_matrix(rbas_i, sphi_ri_i)
         IF (.NOT. onsite_aux) THEN
            CALL hfx_lri_build_sphi_matrix(rbas_j, sphi_ri_j)
         ELSE
            ALLOCATE (sphi_ri_j(0, 0))
         END IF

         sph_metric_dim = SIZE(sphi_ri_i, 2)
         IF (.NOT. onsite_aux) sph_metric_dim = sph_metric_dim + SIZE(sphi_ri_j, 2)
         sph_j_rows = SIZE(sphi_i, 2)*SIZE(sphi_j, 2)
         sph_j_cols = sph_metric_dim

         IF (sph_metric_dim > 0 .AND. sph_j_rows > 0) THEN
            ALLOCATE (sphi_aux(metric_dim, sph_metric_dim), SOURCE=0.0_dp)
            sphi_aux(1:ncart_ri, 1:SIZE(sphi_ri_i, 2)) = sphi_ri_i
            IF (.NOT. onsite_aux) THEN
               sphi_aux(ncart_ri + 1:metric_dim, SIZE(sphi_ri_i, 2) + 1:sph_metric_dim) = sphi_ri_j
            END IF
            CALL hfx_lri_build_pair_sphi(sphi_i, sphi_j, sphi_pair)

            ALLOCATE (tmp_sph(metric_dim, sph_metric_dim), SOURCE=0.0_dp)
            tmp_sph = MATMUL(metric_block, sphi_aux)
            ALLOCATE (metric_sph(sph_metric_dim, sph_metric_dim), SOURCE=0.0_dp)
            metric_sph = MATMUL(TRANSPOSE(sphi_aux), tmp_sph)
            DEALLOCATE (tmp_sph)

            ALLOCATE (tmp_sph(threec_rows, sph_metric_dim), SOURCE=0.0_dp)
            tmp_sph = MATMUL(j_block, sphi_aux)
            ALLOCATE (j_sph(sph_j_rows, sph_j_cols), SOURCE=0.0_dp)
            j_sph = MATMUL(TRANSPOSE(sphi_pair), tmp_sph)
            DEALLOCATE (tmp_sph)

            IF (do_debug_probe) THEN
               sph_metric_norm = SQRT(SUM(metric_sph*metric_sph))
               sph_j_norm = SQRT(SUM(j_sph*j_sph))
            END IF
            CALL hfx_lri_probe_metric_inverse(metric_sph, j_sph, sph_eig_min, sph_eig_max, &
                                              sph_eig_cond, sph_eig_small, sph_c_norm, c_out=c_sph)
            IF (do_debug_probe) THEN
               ALLOCATE (tmp_cvc(SIZE(c_sph, 1), sph_metric_dim), SOURCE=0.0_dp)
               tmp_cvc = MATMUL(c_sph, metric_sph)
               ALLOCATE (cvc_block(SIZE(c_sph, 1), SIZE(c_sph, 1)), SOURCE=0.0_dp)
               cvc_block = MATMUL(tmp_cvc, TRANSPOSE(c_sph))
               sph_cvc_norm = SQRT(SUM(cvc_block*cvc_block))
               DEALLOCATE (cvc_block, tmp_cvc)
            END IF
            CALL hfx_lri_pack_pair_coefficients(pair_block, c_sph, SIZE(sphi_i, 2), &
                                                SIZE(sphi_j, 2), sph_metric_dim)
            CALL hfx_lri_store_pair_metric(pair_block, metric_sph)
            DEALLOCATE (c_sph)

            DEALLOCATE (j_sph, metric_sph, sphi_aux, sphi_pair)
         END IF
         DEALLOCATE (sphi_ri_j, sphi_ri_i, sphi_j, sphi_i)

         DEALLOCATE (j_block, metric_block)
      END IF

      CALL cp_libint_cleanup_2eri(lib2)
      CALL cp_libint_cleanup_3eri(lib3)

      metric_norm = SQRT(metric_norm2)
      threec_i_norm = SQRT(threec_i_norm2)
      threec_j_norm = SQRT(threec_j_norm2)

   END SUBROUTINE hfx_lri_debug_pair_integrals

! **************************************************************************************************
!> \brief ...
!> \param basis_a ...
!> \param basis_b ...
!> \param ra ...
!> \param rb ...
!> \param lib ...
!> \param ri_data ...
!> \param ncart_total ...
!> \param norm2_total ...
!> \param prefactor ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_accumulate_2c_norm(basis_a, basis_b, ra, rb, lib, ri_data, ncart_total, &
                                         norm2_total, prefactor)

      TYPE(gto_basis_set_type), POINTER                  :: basis_a, basis_b
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rb
      TYPE(cp_libint_t), INTENT(INOUT)                   :: lib
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      INTEGER, INTENT(INOUT)                             :: ncart_total
      REAL(dp), INTENT(INOUT)                            :: norm2_total
      REAL(dp), INTENT(IN), OPTIONAL                     :: prefactor

      INTEGER                                            :: iset, jset, ncoa, ncob, nseta, nsetb
      INTEGER, DIMENSION(:), POINTER                     :: lmax_a, lmax_b, lmin_a, lmin_b, npgf_a, &
                                                            npgf_b
      REAL(dp)                                           :: dab, pfac
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: sij
      REAL(dp), DIMENSION(:, :), POINTER                 :: rpgf_a, rpgf_b, zeta_a, zeta_b

      CALL get_gto_basis_set(basis_a, lmax=lmax_a, lmin=lmin_a, npgf=npgf_a, nset=nseta, &
                             pgf_radius=rpgf_a, zet=zeta_a)
      CALL get_gto_basis_set(basis_b, lmax=lmax_b, lmin=lmin_b, npgf=npgf_b, nset=nsetb, &
                             pgf_radius=rpgf_b, zet=zeta_b)

      pfac = 1.0_dp
      IF (PRESENT(prefactor)) pfac = prefactor
      dab = NORM2(rb - ra)

      DO iset = 1, nseta
         ncoa = npgf_a(iset)*ncoset(lmax_a(iset))
         DO jset = 1, nsetb
            ncob = npgf_b(jset)*ncoset(lmax_b(jset))
            IF (ncoa*ncob == 0) CYCLE
            ALLOCATE (sij(ncoa, ncob), SOURCE=0.0_dp)
            CALL eri_2center(sij, lmin_a(iset), lmax_a(iset), npgf_a(iset), zeta_a(:, iset), &
                             rpgf_a(:, iset), ra, lmin_b(jset), lmax_b(jset), npgf_b(jset), &
                             zeta_b(:, jset), rpgf_b(:, jset), rb, dab, lib, ri_data%hfx_pot)
            ncart_total = ncart_total + SIZE(sij)
            norm2_total = norm2_total + pfac*SUM(sij*sij)
            DEALLOCATE (sij)
         END DO
      END DO

   END SUBROUTINE hfx_lri_accumulate_2c_norm

! **************************************************************************************************
!> \brief ...
!> \param basis_a ...
!> \param basis_b ...
!> \param basis_c ...
!> \param ra ...
!> \param rb ...
!> \param rc ...
!> \param lib ...
!> \param ri_data ...
!> \param ncart_total ...
!> \param norm2_total ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_accumulate_3c_norm(basis_a, basis_b, basis_c, ra, rb, rc, lib, ri_data, &
                                         ncart_total, norm2_total)

      TYPE(gto_basis_set_type), POINTER                  :: basis_a, basis_b, basis_c
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rb, rc
      TYPE(cp_libint_t), INTENT(INOUT)                   :: lib
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      INTEGER, INTENT(INOUT)                             :: ncart_total
      REAL(dp), INTENT(INOUT)                            :: norm2_total

      INTEGER                                            :: iset, jset, kset, ncoa, ncob, ncoc, &
                                                            nseta, nsetb, nsetc
      INTEGER, DIMENSION(:), POINTER                     :: lmax_a, lmax_b, lmax_c, lmin_a, lmin_b, &
                                                            lmin_c, npgf_a, npgf_b, npgf_c
      REAL(dp)                                           :: dab, dac, dbc
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: sijk
      REAL(dp), DIMENSION(:, :), POINTER                 :: rpgf_a, rpgf_b, rpgf_c, zeta_a, zeta_b, &
                                                            zeta_c

      CALL get_gto_basis_set(basis_a, lmax=lmax_a, lmin=lmin_a, npgf=npgf_a, nset=nseta, &
                             pgf_radius=rpgf_a, zet=zeta_a)
      CALL get_gto_basis_set(basis_b, lmax=lmax_b, lmin=lmin_b, npgf=npgf_b, nset=nsetb, &
                             pgf_radius=rpgf_b, zet=zeta_b)
      CALL get_gto_basis_set(basis_c, lmax=lmax_c, lmin=lmin_c, npgf=npgf_c, nset=nsetc, &
                             pgf_radius=rpgf_c, zet=zeta_c)

      dab = NORM2(rb - ra)
      dac = NORM2(rc - ra)
      dbc = NORM2(rc - rb)

      DO iset = 1, nseta
         ncoa = npgf_a(iset)*ncoset(lmax_a(iset))
         DO jset = 1, nsetb
            ncob = npgf_b(jset)*ncoset(lmax_b(jset))
            DO kset = 1, nsetc
               ncoc = npgf_c(kset)*ncoset(lmax_c(kset))
               IF (ncoa*ncob*ncoc == 0) CYCLE
               ALLOCATE (sijk(ncoa, ncob, ncoc), SOURCE=0.0_dp)
               CALL eri_3center(sijk, lmin_a(iset), lmax_a(iset), npgf_a(iset), zeta_a(:, iset), &
                                rpgf_a(:, iset), ra, lmin_b(jset), lmax_b(jset), npgf_b(jset), &
                                zeta_b(:, jset), rpgf_b(:, jset), rb, lmin_c(kset), lmax_c(kset), &
                                npgf_c(kset), zeta_c(:, kset), rpgf_c(:, kset), rc, dab, dac, dbc, &
                                lib, ri_data%hfx_pot)
               ncart_total = ncart_total + SIZE(sijk)
               norm2_total = norm2_total + SUM(sijk*sijk)
               DEALLOCATE (sijk)
            END DO
         END DO
      END DO

   END SUBROUTINE hfx_lri_accumulate_3c_norm

! **************************************************************************************************
!> \brief Counts raw Cartesian primitive-function slots used by the debug probe.
!> \param basis ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION hfx_lri_count_cart_functions(basis) RESULT(ncart)

      TYPE(gto_basis_set_type), POINTER                  :: basis

      INTEGER                                            :: iset, nset
      INTEGER, DIMENSION(:), POINTER                     :: lmax, npgf

      ncart = 0
      CALL get_gto_basis_set(basis, lmax=lmax, npgf=npgf, nset=nset)
      DO iset = 1, nset
         ncart = ncart + npgf(iset)*ncoset(lmax(iset))
      END DO

   END FUNCTION hfx_lri_count_cart_functions

! **************************************************************************************************
!> \brief Fills a raw Cartesian two-center integral block in a dense debug metric.
!> \param basis_a ...
!> \param basis_b ...
!> \param ra ...
!> \param rb ...
!> \param lib ...
!> \param ri_data ...
!> \param mat ...
!> \param offset_a ...
!> \param offset_b ...
!> \param symmetric ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_fill_2c_block(basis_a, basis_b, ra, rb, lib, ri_data, mat, offset_a, &
                                    offset_b, symmetric)

      TYPE(gto_basis_set_type), POINTER                  :: basis_a, basis_b
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rb
      TYPE(cp_libint_t), INTENT(INOUT)                   :: lib
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: mat
      INTEGER, INTENT(IN)                                :: offset_a, offset_b
      LOGICAL, INTENT(IN), OPTIONAL                      :: symmetric

      INTEGER                                            :: ia0, ib0, iset, jset, ncoa, ncob, nseta, &
                                                            nsetb
      INTEGER, DIMENSION(:), POINTER                     :: lmax_a, lmax_b, lmin_a, lmin_b, npgf_a, &
                                                            npgf_b
      LOGICAL                                            :: do_sym
      REAL(dp)                                           :: dab
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: sij
      REAL(dp), DIMENSION(:, :), POINTER                 :: rpgf_a, rpgf_b, zeta_a, zeta_b

      CALL get_gto_basis_set(basis_a, lmax=lmax_a, lmin=lmin_a, npgf=npgf_a, nset=nseta, &
                             pgf_radius=rpgf_a, zet=zeta_a)
      CALL get_gto_basis_set(basis_b, lmax=lmax_b, lmin=lmin_b, npgf=npgf_b, nset=nsetb, &
                             pgf_radius=rpgf_b, zet=zeta_b)

      do_sym = .FALSE.
      IF (PRESENT(symmetric)) do_sym = symmetric
      dab = NORM2(rb - ra)
      ia0 = offset_a
      DO iset = 1, nseta
         ncoa = npgf_a(iset)*ncoset(lmax_a(iset))
         ib0 = offset_b
         DO jset = 1, nsetb
            ncob = npgf_b(jset)*ncoset(lmax_b(jset))
            IF (ncoa*ncob > 0) THEN
               ALLOCATE (sij(ncoa, ncob), SOURCE=0.0_dp)
               CALL eri_2center(sij, lmin_a(iset), lmax_a(iset), npgf_a(iset), zeta_a(:, iset), &
                                rpgf_a(:, iset), ra, lmin_b(jset), lmax_b(jset), npgf_b(jset), &
                                zeta_b(:, jset), rpgf_b(:, jset), rb, dab, lib, ri_data%hfx_pot)
               mat(ia0 + 1:ia0 + ncoa, ib0 + 1:ib0 + ncob) = sij
               IF (do_sym) THEN
                  mat(ib0 + 1:ib0 + ncob, ia0 + 1:ia0 + ncoa) = TRANSPOSE(sij)
               END IF
               DEALLOCATE (sij)
            END IF
            ib0 = ib0 + ncob
         END DO
         ia0 = ia0 + ncoa
      END DO

   END SUBROUTINE hfx_lri_fill_2c_block

! **************************************************************************************************
!> \brief Fills a raw Cartesian three-center integral block in a dense debug J matrix.
!> \param basis_a ...
!> \param basis_b ...
!> \param basis_c ...
!> \param ra ...
!> \param rb ...
!> \param rc ...
!> \param lib ...
!> \param ri_data ...
!> \param mat ...
!> \param offset_c ...
!> \param ncart_b ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_fill_3c_block(basis_a, basis_b, basis_c, ra, rb, rc, lib, ri_data, mat, &
                                    offset_c, ncart_b)

      TYPE(gto_basis_set_type), POINTER                  :: basis_a, basis_b, basis_c
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rb, rc
      TYPE(cp_libint_t), INTENT(INOUT)                   :: lib
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: mat
      INTEGER, INTENT(IN)                                :: offset_c, ncart_b

      INTEGER                                            :: ia, ia0, ib, ib0, ic0, iset, jset, kset, &
                                                            ncoa, ncob, ncoc, nseta, nsetb, nsetc
      INTEGER, DIMENSION(:), POINTER                     :: lmax_a, lmax_b, lmax_c, lmin_a, lmin_b, &
                                                            lmin_c, npgf_a, npgf_b, npgf_c
      REAL(dp)                                           :: dab, dac, dbc
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: sijk
      REAL(dp), DIMENSION(:, :), POINTER                 :: rpgf_a, rpgf_b, rpgf_c, zeta_a, zeta_b, &
                                                            zeta_c

      CALL get_gto_basis_set(basis_a, lmax=lmax_a, lmin=lmin_a, npgf=npgf_a, nset=nseta, &
                             pgf_radius=rpgf_a, zet=zeta_a)
      CALL get_gto_basis_set(basis_b, lmax=lmax_b, lmin=lmin_b, npgf=npgf_b, nset=nsetb, &
                             pgf_radius=rpgf_b, zet=zeta_b)
      CALL get_gto_basis_set(basis_c, lmax=lmax_c, lmin=lmin_c, npgf=npgf_c, nset=nsetc, &
                             pgf_radius=rpgf_c, zet=zeta_c)

      dab = NORM2(rb - ra)
      dac = NORM2(rc - ra)
      dbc = NORM2(rc - rb)
      ia0 = 0
      DO iset = 1, nseta
         ncoa = npgf_a(iset)*ncoset(lmax_a(iset))
         ib0 = 0
         DO jset = 1, nsetb
            ncob = npgf_b(jset)*ncoset(lmax_b(jset))
            ic0 = offset_c
            DO kset = 1, nsetc
               ncoc = npgf_c(kset)*ncoset(lmax_c(kset))
               IF (ncoa*ncob*ncoc > 0) THEN
                  ALLOCATE (sijk(ncoa, ncob, ncoc), SOURCE=0.0_dp)
                  CALL eri_3center(sijk, lmin_a(iset), lmax_a(iset), npgf_a(iset), zeta_a(:, iset), &
                                   rpgf_a(:, iset), ra, lmin_b(jset), lmax_b(jset), npgf_b(jset), &
                                   zeta_b(:, jset), rpgf_b(:, jset), rb, lmin_c(kset), lmax_c(kset), &
                                   npgf_c(kset), zeta_c(:, kset), rpgf_c(:, kset), rc, dab, dac, &
                                   dbc, lib, ri_data%hfx_pot)
                  DO ia = 1, ncoa
                     DO ib = 1, ncob
                        mat((ia0 + ia - 1)*ncart_b + ib0 + ib, &
                            ic0 + 1:ic0 + ncoc) = sijk(ia, ib, :)
                     END DO
                  END DO
                  DEALLOCATE (sijk)
               END IF
               ic0 = ic0 + ncoc
            END DO
            ib0 = ib0 + ncob
         END DO
         ia0 = ia0 + ncoa
      END DO

   END SUBROUTINE hfx_lri_fill_3c_block

! **************************************************************************************************
!> \brief Builds the raw-Cartesian to spherical/contracted transformation for one basis.
!> \param basis ...
!> \param transform ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_build_sphi_matrix(basis, transform)

      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: transform

      INTEGER                                            :: cart0, iset, ncart_set, nset, nsgf, sgf0
      INTEGER, DIMENSION(:), POINTER                     :: lmax, npgf, nsgf_set
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf
      REAL(dp), DIMENSION(:, :), POINTER                 :: sphi

      CALL get_gto_basis_set(basis, first_sgf=first_sgf, lmax=lmax, npgf=npgf, nset=nset, &
                             nsgf=nsgf, nsgf_set=nsgf_set, sphi=sphi)
      ALLOCATE (transform(hfx_lri_count_cart_functions(basis), nsgf), SOURCE=0.0_dp)
      cart0 = 0
      DO iset = 1, nset
         ncart_set = npgf(iset)*ncoset(lmax(iset))
         sgf0 = first_sgf(1, iset)
         IF (ncart_set > 0 .AND. nsgf_set(iset) > 0) THEN
            transform(cart0 + 1:cart0 + ncart_set, sgf0:sgf0 + nsgf_set(iset) - 1) = &
               sphi(1:ncart_set, sgf0:sgf0 + nsgf_set(iset) - 1)
         END IF
         cart0 = cart0 + ncart_set
      END DO

   END SUBROUTINE hfx_lri_build_sphi_matrix

! **************************************************************************************************
!> \brief Builds the product transformation for an AO pair.
!> \param sphi_a ...
!> \param sphi_b ...
!> \param pair_sphi ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_build_pair_sphi(sphi_a, sphi_b, pair_sphi)

      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: sphi_a, sphi_b
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: pair_sphi

      INTEGER                                            :: ca, cb, ncart_a, ncart_b, nsgf_a, &
                                                            nsgf_b, sa, sb

      ncart_a = SIZE(sphi_a, 1)
      ncart_b = SIZE(sphi_b, 1)
      nsgf_a = SIZE(sphi_a, 2)
      nsgf_b = SIZE(sphi_b, 2)
      ALLOCATE (pair_sphi(ncart_a*ncart_b, nsgf_a*nsgf_b), SOURCE=0.0_dp)

      DO sa = 1, nsgf_a
         DO sb = 1, nsgf_b
            DO ca = 1, ncart_a
               DO cb = 1, ncart_b
                  pair_sphi((ca - 1)*ncart_b + cb, (sa - 1)*nsgf_b + sb) = &
                     sphi_a(ca, sa)*sphi_b(cb, sb)
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE hfx_lri_build_pair_sphi

! **************************************************************************************************
!> \brief Diagonalizes a debug metric block and probes C = J M^-1.
!> \param metric_block ...
!> \param j_block ...
!> \param eig_min ...
!> \param eig_max ...
!> \param eig_cond ...
!> \param eig_small ...
!> \param c_norm ...
!> \param c_out ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_probe_metric_inverse(metric_block, j_block, eig_min, eig_max, eig_cond, &
                                           eig_small, c_norm, c_out)

      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: metric_block, j_block
      REAL(dp), INTENT(OUT)                              :: eig_min, eig_max, eig_cond
      INTEGER, INTENT(OUT)                               :: eig_small
      REAL(dp), INTENT(OUT)                              :: c_norm
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT), OPTIONAL                           :: c_out

      INTEGER                                            :: i, info, lwork, ncol, nrow
      REAL(dp)                                           :: eig_abs, eig_cutoff
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: eigenvalues, work
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: c_block, eigvec, tmp

      eig_min = 0.0_dp
      eig_max = 0.0_dp
      eig_cond = 0.0_dp
      eig_small = 0
      c_norm = 0.0_dp

      ncol = SIZE(metric_block, 1)
      nrow = SIZE(j_block, 1)
      IF (ncol == 0 .OR. nrow == 0) RETURN

      ALLOCATE (eigvec(ncol, ncol), eigenvalues(ncol))
      eigvec(:, :) = 0.5_dp*(metric_block + TRANSPOSE(metric_block))

      lwork = -1
      ALLOCATE (work(1))
      CALL dsyev('V', 'L', ncol, eigvec, ncol, eigenvalues, work, lwork, info)
      IF (info /= 0) THEN
         CPABORT("DSYEV workspace query failed in LRI-HFXk debug metric")
      END IF
      lwork = MAX(1, INT(work(1)))
      DEALLOCATE (work)

      ALLOCATE (work(lwork))
      CALL dsyev('V', 'L', ncol, eigvec, ncol, eigenvalues, work, lwork, info)
      DEALLOCATE (work)
      IF (info /= 0) THEN
         CPABORT("DSYEV failed in LRI-HFXk debug metric")
      END IF

      eig_max = MAXVAL(ABS(eigenvalues))
      eig_cutoff = MAX(1.0E-12_dp*MAX(1.0_dp, eig_max), 1.0E-14_dp)
      eig_min = HUGE(1.0_dp)
      DO i = 1, ncol
         eig_abs = ABS(eigenvalues(i))
         IF (eig_abs <= eig_cutoff) THEN
            eig_small = eig_small + 1
         ELSE
            eig_min = MIN(eig_min, eig_abs)
         END IF
      END DO
      IF (eig_min < HUGE(1.0_dp) .AND. eig_max > 0.0_dp) eig_cond = eig_max/eig_min

      ALLOCATE (tmp(nrow, ncol), c_block(nrow, ncol))
      tmp(:, :) = MATMUL(j_block, eigvec)
      DO i = 1, ncol
         eig_abs = ABS(eigenvalues(i))
         IF (eig_abs > eig_cutoff) THEN
            tmp(:, i) = tmp(:, i)/eigenvalues(i)
         ELSE
            tmp(:, i) = 0.0_dp
         END IF
      END DO
      c_block(:, :) = MATMUL(tmp, TRANSPOSE(eigvec))
      c_norm = SQRT(SUM(c_block*c_block))
      IF (PRESENT(c_out)) THEN
         IF (ALLOCATED(c_out)) DEALLOCATE (c_out)
         ALLOCATE (c_out(SIZE(c_block, 1), SIZE(c_block, 2)))
         c_out(:, :) = c_block(:, :)
      END IF

      DEALLOCATE (c_block, tmp, eigenvalues, eigvec)

   END SUBROUTINE hfx_lri_probe_metric_inverse

! **************************************************************************************************
!> \brief Initializes a debug pair-local LRI coefficient block.
!> \param pair_block ...
!> \param iatom ...
!> \param jatom ...
!> \param ikind ...
!> \param jkind ...
!> \param cell ...
!> \param rab ...
!> \param nao_i ...
!> \param nao_j ...
!> \param naux_pair ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_init_pair_block(pair_block, iatom, jatom, ikind, jkind, cell, rab, nao_i, &
                                      nao_j, naux_pair)

      TYPE(hfx_lri_pair_block_type), INTENT(INOUT)       :: pair_block
      INTEGER, INTENT(IN)                                :: iatom, jatom, ikind, jkind
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rab
      INTEGER, INTENT(IN)                                :: nao_i, nao_j, naux_pair

      IF (ALLOCATED(pair_block%c_fit)) DEALLOCATE (pair_block%c_fit)
      IF (ALLOCATED(pair_block%v_fit)) DEALLOCATE (pair_block%v_fit)

      pair_block%iatom = iatom
      pair_block%jatom = jatom
      pair_block%ikind = ikind
      pair_block%jkind = jkind
      pair_block%cell(:) = cell(:)
      pair_block%rab(:) = rab(:)
      pair_block%nao_i = nao_i
      pair_block%nao_j = nao_j
      pair_block%naux_pair = naux_pair

   END SUBROUTINE hfx_lri_init_pair_block

! **************************************************************************************************
!> \brief Packs a 2D spherical coefficient matrix into C(mu,nu,P).
!> \param pair_block ...
!> \param c_mat ...
!> \param nao_i ...
!> \param nao_j ...
!> \param naux_pair ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_pack_pair_coefficients(pair_block, c_mat, nao_i, nao_j, naux_pair)

      TYPE(hfx_lri_pair_block_type), INTENT(INOUT)       :: pair_block
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: c_mat
      INTEGER, INTENT(IN)                                :: nao_i, nao_j, naux_pair

      INTEGER                                            :: imu, inu, irow

      IF (SIZE(c_mat, 1) /= nao_i*nao_j) THEN
         CPABORT("Inconsistent row dimension while packing LRI-HFXk pair coefficients")
      END IF
      IF (SIZE(c_mat, 2) /= naux_pair) THEN
         CPABORT("Inconsistent auxiliary dimension while packing LRI-HFXk pair coefficients")
      END IF

      IF (pair_block%nao_i /= nao_i .OR. pair_block%nao_j /= nao_j .OR. &
          pair_block%naux_pair /= naux_pair) THEN
         CPABORT("Inconsistent pair metadata while packing LRI-HFXk coefficients")
      END IF

      IF (ALLOCATED(pair_block%c_fit)) DEALLOCATE (pair_block%c_fit)
      ALLOCATE (pair_block%c_fit(nao_i, nao_j, naux_pair))

      irow = 0
      DO inu = 1, nao_j
         DO imu = 1, nao_i
            irow = irow + 1
            pair_block%c_fit(imu, inu, :) = c_mat(irow, :)
         END DO
      END DO

   END SUBROUTINE hfx_lri_pack_pair_coefficients

! **************************************************************************************************
!> \brief Stores the spherical auxiliary coupling metric for a pair block.
!> \param pair_block ...
!> \param metric_sph ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_store_pair_metric(pair_block, metric_sph)

      TYPE(hfx_lri_pair_block_type), INTENT(INOUT)       :: pair_block
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: metric_sph

      IF (SIZE(metric_sph, 1) /= pair_block%naux_pair .OR. &
          SIZE(metric_sph, 2) /= pair_block%naux_pair) THEN
         CPABORT("Inconsistent metric dimension while storing LRI-HFXk pair metric")
      END IF

      IF (ALLOCATED(pair_block%v_fit)) DEALLOCATE (pair_block%v_fit)
      ALLOCATE (pair_block%v_fit(pair_block%naux_pair, pair_block%naux_pair))
      pair_block%v_fit(:, :) = metric_sph(:, :)

   END SUBROUTINE hfx_lri_store_pair_metric

! **************************************************************************************************
!> \brief Reads a stored-pair density block from the current real-space density matrices.
!> \brief Stores a debug K block together with its target pair metadata.
!> \param k_pair ...
!> \param pair_block ...
!> \param pair_index ...
!> \param k_block ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_store_debug_k_block(k_pair, pair_block, pair_index, k_block)

      USE hfx_lri_types, ONLY: hfx_lri_k_block_type
      TYPE(hfx_lri_k_block_type), INTENT(INOUT)          :: k_pair
      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair_block
      INTEGER, INTENT(IN)                                :: pair_index
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: k_block

      IF (ALLOCATED(k_pair%k_block)) DEALLOCATE (k_pair%k_block)

      k_pair%pair_index = pair_index
      k_pair%iatom = pair_block%iatom
      k_pair%jatom = pair_block%jatom
      k_pair%ikind = pair_block%ikind
      k_pair%jkind = pair_block%jkind
      k_pair%cell(:) = pair_block%cell(:)
      k_pair%nao_i = SIZE(k_block, 1)
      k_pair%nao_j = SIZE(k_block, 2)

      ALLOCATE (k_pair%k_block(SIZE(k_block, 1), SIZE(k_block, 2)))
      k_pair%k_block(:, :) = k_block(:, :)

   END SUBROUTINE hfx_lri_store_debug_k_block

! **************************************************************************************************
!> \brief Accumulates a debug K block for one target pair using nearest-origin density pairs.
!> \param basis_set_RI ...
!> \param rho_ao ...
!> \param cell_to_index ...
!> \param ipair_target ...
!> \param ri_data ...
!> \param k_accum_block ...
!> \param used_pairs ...
!> \param candidate_pairs ...
!> \param rho_accum_norm ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_accumulate_target_k_block(basis_set_RI, rho_ao, cell_to_index, ipair_target, &
                                                ri_data, k_accum_block, used_pairs, &
                                                candidate_pairs, rho_accum_norm)

      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      INTEGER, INTENT(IN)                                :: ipair_target
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: k_accum_block
      INTEGER, INTENT(OUT)                               :: used_pairs, candidate_pairs
      REAL(dp), INTENT(OUT)                              :: rho_accum_norm

      INTEGER                                            :: img_pair, ipair, ipair_best, lb1, lb2, &
                                                            lb3, search_iter, ub1, ub2, ub3
      LOGICAL                                            :: found
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: tried_density_pair
      REAL(dp)                                           :: best_r2, pair_r2, rho_norm2_accum, &
                                                            rho_pair_norm
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: k_pair_block
      REAL(dp), DIMENSION(:, :), POINTER                 :: rho_block

      NULLIFY (rho_block)

      used_pairs = 0
      candidate_pairs = 0
      rho_accum_norm = 0.0_dp

      IF (.NOT. ASSOCIATED(ri_data%lri_kp)) RETURN
      IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks)) RETURN
      IF (ipair_target < 1 .OR. ipair_target > ri_data%lri_kp%npair_blocks) RETURN
      IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks(ipair_target)%c_fit)) RETURN
      IF (.NOT. ASSOCIATED(rho_ao)) RETURN

      ALLOCATE (k_accum_block(ri_data%lri_kp%pair_blocks(ipair_target)%nao_i, &
                              ri_data%lri_kp%pair_blocks(ipair_target)%nao_j), SOURCE=0.0_dp)
      ALLOCATE (tried_density_pair(ri_data%lri_kp%npair_blocks), SOURCE=.FALSE.)
      rho_norm2_accum = 0.0_dp

      DO search_iter = 1, ri_data%lri_kp%npair_blocks
         ipair_best = 0
         best_r2 = HUGE(1.0_dp)
         DO ipair = 1, ri_data%lri_kp%npair_blocks
            IF (tried_density_pair(ipair)) CYCLE
            IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks(ipair)%c_fit)) CYCLE
            IF (ri_data%lri_kp%pair_blocks(ipair)%iatom /= &
                ri_data%lri_kp%pair_blocks(ipair_target)%iatom) CYCLE
            pair_r2 = SUM((ri_data%lri_kp%pair_blocks(ipair)%rab - &
                           ri_data%lri_kp%pair_blocks(ipair_target)%rab)**2)
            IF (pair_r2 >= best_r2) CYCLE
            ipair_best = ipair
            best_r2 = pair_r2
         END DO
         IF (ipair_best == 0) EXIT
         tried_density_pair(ipair_best) = .TRUE.
         candidate_pairs = candidate_pairs + 1

         img_pair = 1
         IF (ASSOCIATED(cell_to_index)) THEN
            lb1 = LBOUND(cell_to_index, 1); ub1 = UBOUND(cell_to_index, 1)
            lb2 = LBOUND(cell_to_index, 2); ub2 = UBOUND(cell_to_index, 2)
            lb3 = LBOUND(cell_to_index, 3); ub3 = UBOUND(cell_to_index, 3)
            IF (.NOT. (ALL(ri_data%lri_kp%pair_blocks(ipair_best)%cell >= [lb1, lb2, lb3]) .AND. &
                       ALL(ri_data%lri_kp%pair_blocks(ipair_best)%cell <= [ub1, ub2, ub3]))) CYCLE
            img_pair = cell_to_index(ri_data%lri_kp%pair_blocks(ipair_best)%cell(1), &
                                     ri_data%lri_kp%pair_blocks(ipair_best)%cell(2), &
                                     ri_data%lri_kp%pair_blocks(ipair_best)%cell(3))
         END IF
         IF (img_pair < 1 .OR. img_pair > SIZE(rho_ao, 2)) CYCLE

         CALL dbcsr_get_block_p(matrix=rho_ao(1, img_pair)%matrix, &
                                row=ri_data%lri_kp%pair_blocks(ipair_best)%iatom, &
                                col=ri_data%lri_kp%pair_blocks(ipair_best)%jatom, &
                                block=rho_block, found=found)
         IF (.NOT. found) CYCLE
         rho_pair_norm = SQRT(SUM(rho_block*rho_block))
         IF (rho_pair_norm < lri_debug_rho_threshold) CYCLE

         CALL hfx_lri_contract_pair_pair_density_realv(basis_set_RI, &
                                                       ri_data%lri_kp%pair_blocks(ipair_target), &
                                                       ri_data%lri_kp%pair_blocks(ipair_best), &
                                                       ri_data, rho_block, k_pair_block)
         IF (.NOT. ALLOCATED(k_pair_block)) CYCLE
         k_accum_block(:, :) = k_accum_block(:, :) + k_pair_block(:, :)
         rho_norm2_accum = rho_norm2_accum + rho_pair_norm*rho_pair_norm
         used_pairs = used_pairs + 1
         DEALLOCATE (k_pair_block)
      END DO

      IF (ALLOCATED(tried_density_pair)) DEALLOCATE (tried_density_pair)
      rho_accum_norm = SQRT(rho_norm2_accum)
      IF (used_pairs <= 0) THEN
         IF (ALLOCATED(k_accum_block)) DEALLOCATE (k_accum_block)
      END IF

   END SUBROUTINE hfx_lri_accumulate_target_k_block

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ri_data ...
!> \param rho_ao ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_debug_density_block(qs_env, ri_data, rho_ao)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao

      INTEGER :: candidate_pairs, candidate_pairs_tmp, img, img_pair, ipair, ipair_best, &
         ipair_sel, ipair_target, lb1, lb2, lb3, max_k_blocks, multi_candidate_pairs, &
         multi_used_pairs, nkind, search_iter, ub1, ub2, ub3, used_pairs, used_pairs_tmp
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: tried_density_pair
      REAL(dp) :: best_r2, c_norm, k_accum_norm, k_blocks_norm, k_norm, pair_r2, rho_accum_norm, &
         rho_blocks_norm, rho_norm, rho_norm2_accum, rho_pair_norm, rho_tmp_norm
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: k_accum_block, k_block, k_pair_block, &
                                                            k_tmp_block
      REAL(dp), DIMENSION(:, :), POINTER                 :: rho_block
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_RI
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (basis_set_RI, cell_to_index, kpoints, qs_kind_set, rho_block)

      IF (.NOT. ASSOCIATED(ri_data%lri_kp)) THEN
         IF (ri_data%unit_nr > 0) WRITE (ri_data%unit_nr, "(/,T4,A)") &
            "debug density block skipped: LRI-HFXk state is not associated"
         RETURN
      END IF
      IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks)) THEN
         IF (ri_data%unit_nr > 0) WRITE (ri_data%unit_nr, "(/,T4,A)") &
            "debug density block skipped: pair blocks are not allocated"
         RETURN
      END IF
      IF (ri_data%lri_kp%npair_blocks <= 0) THEN
         IF (ri_data%unit_nr > 0) WRITE (ri_data%unit_nr, "(/,T4,A)") &
            "debug density block skipped: no pair blocks are available"
         RETURN
      END IF
      IF (.NOT. ASSOCIATED(rho_ao)) THEN
         IF (ri_data%unit_nr > 0) WRITE (ri_data%unit_nr, "(/,T4,A)") &
            "debug density block skipped: rho_ao is not associated"
         RETURN
      END IF

      ipair_sel = 0
      DO ipair = 1, ri_data%lri_kp%npair_blocks
         IF (ri_data%lri_kp%pair_blocks(ipair)%iatom == &
             ri_data%lri_kp%pair_blocks(ipair)%jatom .AND. &
             ALL(ri_data%lri_kp%pair_blocks(ipair)%cell == 0)) THEN
            ipair_sel = ipair
            EXIT
         END IF
      END DO
      IF (ipair_sel == 0) ipair_sel = 1
      IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks(ipair_sel)%c_fit)) THEN
         IF (ri_data%unit_nr > 0) THEN
            WRITE (ri_data%unit_nr, "(/,T4,A)") &
               "debug density block skipped: selected C block is not allocated"
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10)") "pair index", ipair_sel
         END IF
         RETURN
      END IF

      CALL get_qs_env(qs_env=qs_env, kpoints=kpoints, qs_kind_set=qs_kind_set, nkind=nkind)
      IF (ASSOCIATED(kpoints)) THEN
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
      END IF

      ALLOCATE (basis_set_RI(nkind))
      CALL basis_set_list_setup(basis_set_RI, ri_data%ri_basis_type, qs_kind_set)

      img = 1
      IF (ASSOCIATED(cell_to_index)) THEN
         lb1 = LBOUND(cell_to_index, 1); ub1 = UBOUND(cell_to_index, 1)
         lb2 = LBOUND(cell_to_index, 2); ub2 = UBOUND(cell_to_index, 2)
         lb3 = LBOUND(cell_to_index, 3); ub3 = UBOUND(cell_to_index, 3)
         IF (ALL(ri_data%lri_kp%pair_blocks(ipair_sel)%cell >= [lb1, lb2, lb3]) .AND. &
             ALL(ri_data%lri_kp%pair_blocks(ipair_sel)%cell <= [ub1, ub2, ub3])) THEN
            img = cell_to_index(ri_data%lri_kp%pair_blocks(ipair_sel)%cell(1), &
                                ri_data%lri_kp%pair_blocks(ipair_sel)%cell(2), &
                                ri_data%lri_kp%pair_blocks(ipair_sel)%cell(3))
         END IF
      END IF
      IF (img < 1 .OR. img > SIZE(rho_ao, 2)) THEN
         IF (ri_data%unit_nr > 0) THEN
            WRITE (ri_data%unit_nr, "(/,T4,A)") &
               "debug density block skipped: invalid image index"
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") "img/nimages", img, SIZE(rho_ao, 2)
            WRITE (ri_data%unit_nr, "(T6,A,T28,3I6)") "cell", &
               ri_data%lri_kp%pair_blocks(ipair_sel)%cell
         END IF
         IF (ASSOCIATED(basis_set_RI)) DEALLOCATE (basis_set_RI)
         RETURN
      END IF

      CALL dbcsr_get_block_p(matrix=rho_ao(1, img)%matrix, &
                             row=ri_data%lri_kp%pair_blocks(ipair_sel)%iatom, &
                             col=ri_data%lri_kp%pair_blocks(ipair_sel)%jatom, &
                             block=rho_block, found=found)

      rho_norm = 0.0_dp
      k_norm = 0.0_dp
      IF (found) rho_norm = SQRT(SUM(rho_block*rho_block))
      IF (found .AND. ALLOCATED(ri_data%lri_kp%pair_blocks(ipair_sel)%v_fit)) THEN
         CALL hfx_lri_contract_real_density_self(ri_data%lri_kp%pair_blocks(ipair_sel), rho_block, &
                                                 k_block)
         k_norm = SQRT(SUM(k_block*k_block))
         DEALLOCATE (k_block)
      END IF
      c_norm = SQRT(SUM(ri_data%lri_kp%pair_blocks(ipair_sel)%c_fit* &
                        ri_data%lri_kp%pair_blocks(ipair_sel)%c_fit))

      used_pairs = 0
      candidate_pairs = 0
      rho_norm2_accum = 0.0_dp
      IF (ALLOCATED(k_accum_block)) DEALLOCATE (k_accum_block)
      ALLOCATE (k_accum_block(ri_data%lri_kp%pair_blocks(ipair_sel)%nao_i, &
                              ri_data%lri_kp%pair_blocks(ipair_sel)%nao_j), SOURCE=0.0_dp)
      ALLOCATE (tried_density_pair(ri_data%lri_kp%npair_blocks), SOURCE=.FALSE.)
      DO search_iter = 1, ri_data%lri_kp%npair_blocks

         ipair_best = 0
         best_r2 = HUGE(1.0_dp)
         DO ipair = 1, ri_data%lri_kp%npair_blocks
            IF (tried_density_pair(ipair)) CYCLE
            IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks(ipair)%c_fit)) CYCLE
            IF (ri_data%lri_kp%pair_blocks(ipair)%iatom /= &
                ri_data%lri_kp%pair_blocks(ipair_sel)%iatom) CYCLE
            pair_r2 = SUM((ri_data%lri_kp%pair_blocks(ipair)%rab - &
                           ri_data%lri_kp%pair_blocks(ipair_sel)%rab)**2)
            IF (pair_r2 >= best_r2) CYCLE
            ipair_best = ipair
            best_r2 = pair_r2
         END DO
         IF (ipair_best == 0) EXIT
         tried_density_pair(ipair_best) = .TRUE.
         candidate_pairs = candidate_pairs + 1

         img_pair = 1
         IF (ASSOCIATED(cell_to_index)) THEN
            IF (.NOT. (ALL(ri_data%lri_kp%pair_blocks(ipair_best)%cell >= [lb1, lb2, lb3]) .AND. &
                       ALL(ri_data%lri_kp%pair_blocks(ipair_best)%cell <= [ub1, ub2, ub3]))) CYCLE
            img_pair = cell_to_index(ri_data%lri_kp%pair_blocks(ipair_best)%cell(1), &
                                     ri_data%lri_kp%pair_blocks(ipair_best)%cell(2), &
                                     ri_data%lri_kp%pair_blocks(ipair_best)%cell(3))
         END IF
         IF (img_pair < 1 .OR. img_pair > SIZE(rho_ao, 2)) CYCLE

         CALL dbcsr_get_block_p(matrix=rho_ao(1, img_pair)%matrix, &
                                row=ri_data%lri_kp%pair_blocks(ipair_best)%iatom, &
                                col=ri_data%lri_kp%pair_blocks(ipair_best)%jatom, &
                                block=rho_block, found=found)
         IF (.NOT. found) CYCLE
         rho_pair_norm = SQRT(SUM(rho_block*rho_block))
         IF (rho_pair_norm < lri_debug_rho_threshold) CYCLE

         CALL hfx_lri_contract_pair_pair_density_realv(basis_set_RI, &
                                                       ri_data%lri_kp%pair_blocks(ipair_sel), &
                                                       ri_data%lri_kp%pair_blocks(ipair_best), &
                                                       ri_data, rho_block, k_pair_block)
         IF (.NOT. ALLOCATED(k_pair_block)) CYCLE
         k_accum_block(:, :) = k_accum_block(:, :) + k_pair_block(:, :)
         rho_norm2_accum = rho_norm2_accum + rho_pair_norm*rho_pair_norm
         used_pairs = used_pairs + 1
         DEALLOCATE (k_pair_block)
      END DO
      IF (ALLOCATED(tried_density_pair)) DEALLOCATE (tried_density_pair)
      rho_accum_norm = SQRT(rho_norm2_accum)
      k_accum_norm = SQRT(SUM(k_accum_block*k_accum_block))

      CALL hfx_lri_store_debug_k_block(ri_data%lri_kp%debug_k_block, &
                                       ri_data%lri_kp%pair_blocks(ipair_sel), &
                                       ipair_sel, k_accum_block)

      IF (ALLOCATED(ri_data%lri_kp%debug_k_blocks)) THEN
         DO ipair = 1, SIZE(ri_data%lri_kp%debug_k_blocks)
            IF (ALLOCATED(ri_data%lri_kp%debug_k_blocks(ipair)%k_block)) THEN
               DEALLOCATE (ri_data%lri_kp%debug_k_blocks(ipair)%k_block)
            END IF
         END DO
         DEALLOCATE (ri_data%lri_kp%debug_k_blocks)
      END IF
      ri_data%lri_kp%ndebug_k_blocks = 0
      max_k_blocks = MIN(lri_debug_max_k_blocks, ri_data%lri_kp%npair_blocks)
      IF (max_k_blocks > 0) THEN
         ALLOCATE (ri_data%lri_kp%debug_k_blocks(max_k_blocks))
      END IF
      multi_used_pairs = 0
      multi_candidate_pairs = 0
      rho_blocks_norm = 0.0_dp
      k_blocks_norm = 0.0_dp
      DO ipair_target = 1, ri_data%lri_kp%npair_blocks
         IF (ri_data%lri_kp%ndebug_k_blocks >= max_k_blocks) EXIT
         IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks(ipair_target)%c_fit)) CYCLE
         CALL hfx_lri_accumulate_target_k_block(basis_set_RI, rho_ao, cell_to_index, ipair_target, &
                                                ri_data, k_tmp_block, used_pairs_tmp, &
                                                candidate_pairs_tmp, rho_tmp_norm)
         IF (.NOT. ALLOCATED(k_tmp_block)) CYCLE
         IF (used_pairs_tmp <= 0) THEN
            DEALLOCATE (k_tmp_block)
            CYCLE
         END IF
         ri_data%lri_kp%ndebug_k_blocks = ri_data%lri_kp%ndebug_k_blocks + 1
         CALL hfx_lri_store_debug_k_block( &
            ri_data%lri_kp%debug_k_blocks(ri_data%lri_kp%ndebug_k_blocks), &
            ri_data%lri_kp%pair_blocks(ipair_target), ipair_target, k_tmp_block)
         multi_used_pairs = multi_used_pairs + used_pairs_tmp
         multi_candidate_pairs = multi_candidate_pairs + candidate_pairs_tmp
         rho_blocks_norm = rho_blocks_norm + rho_tmp_norm*rho_tmp_norm
         k_blocks_norm = k_blocks_norm + SUM(k_tmp_block*k_tmp_block)
         DEALLOCATE (k_tmp_block)
      END DO
      rho_blocks_norm = SQRT(rho_blocks_norm)
      k_blocks_norm = SQRT(k_blocks_norm)

      IF (ri_data%unit_nr > 0) THEN
         WRITE (ri_data%unit_nr, "(/,T4,A)") "debug density block"
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
            "pair/image", ipair_sel, img
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
            "atoms", ri_data%lri_kp%pair_blocks(ipair_sel)%iatom, &
            ri_data%lri_kp%pair_blocks(ipair_sel)%jatom
         WRITE (ri_data%unit_nr, "(T6,A,T28,3I6)") &
            "cell", ri_data%lri_kp%pair_blocks(ipair_sel)%cell
         IF (found) THEN
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
               "rho dims", SIZE(rho_block, 1), SIZE(rho_block, 2)
         ELSE
            WRITE (ri_data%unit_nr, "(T6,A,T28,A)") "rho dims", "not found"
         END IF
         WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3,T41,ES10.3)") &
            "rho/C norm", rho_norm, c_norm
         WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3)") &
            "real K norm", k_norm
         WRITE (ri_data%unit_nr, "(/,T4,A)") "debug nearest-origin K accumulation"
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
            "used/candidates", used_pairs, candidate_pairs
         WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3)") &
            "rho threshold", lri_debug_rho_threshold
         WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3,T41,ES10.3)") &
            "rho/K norm", rho_accum_norm, k_accum_norm
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
            "stored K dims", SIZE(ri_data%lri_kp%debug_k_block%k_block, 1), &
            SIZE(ri_data%lri_kp%debug_k_block%k_block, 2)
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10)") &
            "stored K pair", ri_data%lri_kp%debug_k_block%pair_index
         WRITE (ri_data%unit_nr, "(/,T4,A)") "debug multi-K block accumulation"
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
            "stored/max blocks", ri_data%lri_kp%ndebug_k_blocks, max_k_blocks
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
            "used/candidates", multi_used_pairs, multi_candidate_pairs
         WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3,T41,ES10.3)") &
            "rho/K norm", rho_blocks_norm, k_blocks_norm
      END IF

      IF (ALLOCATED(k_accum_block)) DEALLOCATE (k_accum_block)
      IF (ASSOCIATED(basis_set_RI)) DEALLOCATE (basis_set_RI)

   END SUBROUTINE hfx_lri_debug_density_block

! **************************************************************************************************
!> \brief Probes whether the stored debug K block can be located in the k-point KS matrix.
!> \param qs_env ...
!> \param ri_data ...
!> \param ks_matrix ...
!> \param hf_fraction ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_debug_ks_writeback_probe(qs_env, ri_data, ks_matrix, hf_fraction)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      REAL(dp), INTENT(IN)                               :: hf_fraction

      INTEGER                                            :: added_blocks, img, ipair, lb1, lb2, lb3, &
                                                            ub1, ub2, ub3
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found
      REAL(dp)                                           :: after_norm, before_norm, k_norm, &
                                                            multi_add_norm2, multi_after_norm2, &
                                                            multi_before_norm2, scale
      REAL(dp), DIMENSION(:, :), POINTER                 :: ks_block
      TYPE(kpoint_type), POINTER                         :: kpoints

      NULLIFY (cell_to_index, kpoints, ks_block)

      IF (.NOT. ASSOCIATED(ri_data%lri_kp)) THEN
         IF (ri_data%unit_nr > 0) WRITE (ri_data%unit_nr, "(/,T4,A)") &
            "debug KS writeback skipped: LRI-HFXk state is not associated"
         RETURN
      END IF
      IF (.NOT. ALLOCATED(ri_data%lri_kp%debug_k_block%k_block)) THEN
         IF (ri_data%unit_nr > 0) WRITE (ri_data%unit_nr, "(/,T4,A)") &
            "debug KS writeback skipped: debug K block is not allocated"
         RETURN
      END IF
      IF (.NOT. ASSOCIATED(ks_matrix)) THEN
         IF (ri_data%unit_nr > 0) WRITE (ri_data%unit_nr, "(/,T4,A)") &
            "debug KS writeback skipped: KS matrix is not associated"
         RETURN
      END IF

      CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)
      IF (ASSOCIATED(kpoints)) THEN
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
      END IF

      img = 1
      IF (ASSOCIATED(cell_to_index)) THEN
         lb1 = LBOUND(cell_to_index, 1); ub1 = UBOUND(cell_to_index, 1)
         lb2 = LBOUND(cell_to_index, 2); ub2 = UBOUND(cell_to_index, 2)
         lb3 = LBOUND(cell_to_index, 3); ub3 = UBOUND(cell_to_index, 3)
         IF (ALL(ri_data%lri_kp%debug_k_block%cell >= [lb1, lb2, lb3]) .AND. &
             ALL(ri_data%lri_kp%debug_k_block%cell <= [ub1, ub2, ub3])) THEN
            img = cell_to_index(ri_data%lri_kp%debug_k_block%cell(1), &
                                ri_data%lri_kp%debug_k_block%cell(2), &
                                ri_data%lri_kp%debug_k_block%cell(3))
         END IF
      END IF
      IF (img < 1 .OR. img > SIZE(ks_matrix, 2)) THEN
         IF (ri_data%unit_nr > 0) THEN
            WRITE (ri_data%unit_nr, "(/,T4,A)") &
               "debug KS writeback skipped: invalid image index"
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") "img/nimages", img, SIZE(ks_matrix, 2)
            WRITE (ri_data%unit_nr, "(T6,A,T28,3I6)") "cell", ri_data%lri_kp%debug_k_block%cell
         END IF
         RETURN
      END IF

      CALL dbcsr_get_block_p(matrix=ks_matrix(1, img)%matrix, &
                             row=ri_data%lri_kp%debug_k_block%iatom, &
                             col=ri_data%lri_kp%debug_k_block%jatom, &
                             block=ks_block, found=found)

      k_norm = SQRT(SUM(ri_data%lri_kp%debug_k_block%k_block* &
                        ri_data%lri_kp%debug_k_block%k_block))
      scale = hf_fraction
      before_norm = 0.0_dp
      after_norm = 0.0_dp
      IF (found) THEN
         IF (SIZE(ks_block, 1) /= SIZE(ri_data%lri_kp%debug_k_block%k_block, 1) .OR. &
             SIZE(ks_block, 2) /= SIZE(ri_data%lri_kp%debug_k_block%k_block, 2)) THEN
            CPABORT("Inconsistent KS block size in LRI-HFXk debug writeback")
         END IF
         before_norm = SQRT(SUM(ks_block*ks_block))
         ks_block(:, :) = ks_block(:, :) + scale*ri_data%lri_kp%debug_k_block%k_block(:, :)
         after_norm = SQRT(SUM(ks_block*ks_block))
      END IF

      added_blocks = 0
      multi_before_norm2 = 0.0_dp
      multi_add_norm2 = 0.0_dp
      multi_after_norm2 = 0.0_dp
      IF (ALLOCATED(ri_data%lri_kp%debug_k_blocks)) THEN
         DO ipair = 1, ri_data%lri_kp%ndebug_k_blocks
            IF (.NOT. ALLOCATED(ri_data%lri_kp%debug_k_blocks(ipair)%k_block)) CYCLE
            IF (ri_data%lri_kp%debug_k_blocks(ipair)%pair_index == &
                ri_data%lri_kp%debug_k_block%pair_index) CYCLE
            img = 1
            IF (ASSOCIATED(cell_to_index)) THEN
               IF (.NOT. (ALL(ri_data%lri_kp%debug_k_blocks(ipair)%cell >= [lb1, lb2, lb3]) .AND. &
                          ALL(ri_data%lri_kp%debug_k_blocks(ipair)%cell <= [ub1, ub2, ub3]))) CYCLE
               img = cell_to_index(ri_data%lri_kp%debug_k_blocks(ipair)%cell(1), &
                                   ri_data%lri_kp%debug_k_blocks(ipair)%cell(2), &
                                   ri_data%lri_kp%debug_k_blocks(ipair)%cell(3))
            END IF
            IF (img < 1 .OR. img > SIZE(ks_matrix, 2)) CYCLE
            CALL dbcsr_get_block_p(matrix=ks_matrix(1, img)%matrix, &
                                   row=ri_data%lri_kp%debug_k_blocks(ipair)%iatom, &
                                   col=ri_data%lri_kp%debug_k_blocks(ipair)%jatom, &
                                   block=ks_block, found=found)
            IF (.NOT. found) CYCLE
            IF (SIZE(ks_block, 1) /= SIZE(ri_data%lri_kp%debug_k_blocks(ipair)%k_block, 1) .OR. &
                SIZE(ks_block, 2) /= SIZE(ri_data%lri_kp%debug_k_blocks(ipair)%k_block, 2)) THEN
               CPABORT("Inconsistent KS block size in LRI-HFXk multi-block debug writeback")
            END IF
            k_norm = SQRT(SUM(ri_data%lri_kp%debug_k_blocks(ipair)%k_block* &
                              ri_data%lri_kp%debug_k_blocks(ipair)%k_block))
            multi_before_norm2 = multi_before_norm2 + SUM(ks_block*ks_block)
            multi_add_norm2 = multi_add_norm2 + (ABS(scale)*k_norm)**2
            ks_block(:, :) = ks_block(:, :) + scale*ri_data%lri_kp%debug_k_blocks(ipair)%k_block(:, :)
            multi_after_norm2 = multi_after_norm2 + SUM(ks_block*ks_block)
            added_blocks = added_blocks + 1
         END DO
      END IF

      IF (ri_data%unit_nr > 0) THEN
         WRITE (ri_data%unit_nr, "(/,T4,A)") "debug KS writeback prototype"
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
            "pair/image", ri_data%lri_kp%debug_k_block%pair_index, img
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
            "atoms", ri_data%lri_kp%debug_k_block%iatom, &
            ri_data%lri_kp%debug_k_block%jatom
         WRITE (ri_data%unit_nr, "(T6,A,T28,3I6)") &
            "cell", ri_data%lri_kp%debug_k_block%cell
         IF (found) THEN
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
               "KS dims", SIZE(ks_block, 1), SIZE(ks_block, 2)
         ELSE
            WRITE (ri_data%unit_nr, "(T6,A,T28,A)") "KS dims", "not found"
         END IF
         WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3,T41,ES10.3)") &
            "before/K norm", before_norm, k_norm
         WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3,T41,ES10.3)") &
            "scale/add norm", scale, ABS(scale)*k_norm
         WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3)") &
            "after norm", after_norm
         WRITE (ri_data%unit_nr, "(/,T4,A)") "debug multi-K writeback prototype"
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
            "added/stored", added_blocks, ri_data%lri_kp%ndebug_k_blocks
         WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3,T41,ES10.3)") &
            "before/add norm", SQRT(multi_before_norm2), SQRT(multi_add_norm2)
         WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3)") &
            "after norm", SQRT(multi_after_norm2)
      END IF

   END SUBROUTINE hfx_lri_debug_ks_writeback_probe

! **************************************************************************************************
!> \brief Contracts a target pair with a density pair and a real pair-pair auxiliary metric.
!> \param basis_set_RI ...
!> \param pair_target ...
!> \param pair_density ...
!> \param ri_data ...
!> \param rho_block ...
!> \param k_block ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_contract_pair_pair_density_realv(basis_set_RI, pair_target, pair_density, &
                                                       ri_data, rho_block, k_block)

      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair_target, pair_density
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: rho_block
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: k_block

      INTEGER :: iao, irow, jao, maxl_li, maxl_lj, maxl_ri, maxl_rj, naux_left, naux_right, &
         ncart_li, ncart_lj, ncart_ri, ncart_rj, nraw_left, nraw_right, offset_left_j, &
         offset_right_j
      LOGICAL                                            :: onsite_left, onsite_right
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: aux_density, aux_target, d_vec, k_vec
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: c_density, c_target, sphi_left, sphi_li, &
                                                            sphi_lj, sphi_ri, sphi_right, sphi_rj, &
                                                            tmp, v_raw, v_sph
      REAL(dp), DIMENSION(3)                             :: r0
      TYPE(cp_libint_t)                                  :: lib2
      TYPE(gto_basis_set_type), POINTER                  :: rbas_li, rbas_lj, rbas_ri, rbas_rj

      NULLIFY (rbas_li, rbas_lj, rbas_ri, rbas_rj)

      IF (.NOT. ALLOCATED(pair_target%c_fit)) RETURN
      IF (.NOT. ALLOCATED(pair_density%c_fit)) RETURN
      IF (pair_target%iatom /= pair_density%iatom) RETURN
      IF (SIZE(rho_block, 1) < pair_density%nao_i .OR. &
          SIZE(rho_block, 2) < pair_density%nao_j) RETURN

      rbas_li => basis_set_RI(pair_target%ikind)%gto_basis_set
      rbas_lj => basis_set_RI(pair_target%jkind)%gto_basis_set
      rbas_ri => basis_set_RI(pair_density%ikind)%gto_basis_set
      rbas_rj => basis_set_RI(pair_density%jkind)%gto_basis_set
      IF (.NOT. ASSOCIATED(rbas_li)) RETURN
      IF (.NOT. ASSOCIATED(rbas_lj)) RETURN
      IF (.NOT. ASSOCIATED(rbas_ri)) RETURN
      IF (.NOT. ASSOCIATED(rbas_rj)) RETURN

      onsite_left = (pair_target%iatom == pair_target%jatom .AND. NORM2(pair_target%rab) <= 1.0E-12_dp)
      onsite_right = (pair_density%iatom == pair_density%jatom .AND. NORM2(pair_density%rab) <= 1.0E-12_dp)

      CALL get_gto_basis_set(rbas_li, maxl=maxl_li)
      CALL get_gto_basis_set(rbas_lj, maxl=maxl_lj)
      CALL get_gto_basis_set(rbas_ri, maxl=maxl_ri)
      CALL get_gto_basis_set(rbas_rj, maxl=maxl_rj)
      CALL init_md_ftable(nmax=MAX(maxl_li + maxl_ri, maxl_li + maxl_rj, &
                                   maxl_lj + maxl_ri, maxl_lj + maxl_rj))
      CALL cp_libint_init_2eri(lib2, MAX(maxl_li, maxl_lj, maxl_ri, maxl_rj))
      CALL cp_libint_set_contrdepth(lib2, 1)

      ncart_li = hfx_lri_count_cart_functions(rbas_li)
      ncart_lj = hfx_lri_count_cart_functions(rbas_lj)
      ncart_ri = hfx_lri_count_cart_functions(rbas_ri)
      ncart_rj = hfx_lri_count_cart_functions(rbas_rj)
      nraw_left = ncart_li
      IF (.NOT. onsite_left) nraw_left = nraw_left + ncart_lj
      nraw_right = ncart_ri
      IF (.NOT. onsite_right) nraw_right = nraw_right + ncart_rj

      r0(:) = 0.0_dp
      ALLOCATE (v_raw(nraw_left, nraw_right), SOURCE=0.0_dp)
      CALL hfx_lri_fill_2c_block(rbas_li, rbas_ri, r0, r0, lib2, ri_data, v_raw, 0, 0)
      IF (.NOT. onsite_right) THEN
         offset_right_j = ncart_ri
         CALL hfx_lri_fill_2c_block(rbas_li, rbas_rj, r0, pair_density%rab, lib2, ri_data, v_raw, 0, &
                                    offset_right_j)
      END IF
      IF (.NOT. onsite_left) THEN
         offset_left_j = ncart_li
         CALL hfx_lri_fill_2c_block(rbas_lj, rbas_ri, pair_target%rab, r0, lib2, ri_data, v_raw, &
                                    offset_left_j, 0)
         IF (.NOT. onsite_right) THEN
            offset_right_j = ncart_ri
            CALL hfx_lri_fill_2c_block(rbas_lj, rbas_rj, pair_target%rab, pair_density%rab, lib2, &
                                       ri_data, v_raw, offset_left_j, offset_right_j)
         END IF
      END IF

      CALL hfx_lri_build_sphi_matrix(rbas_li, sphi_li)
      CALL hfx_lri_build_sphi_matrix(rbas_ri, sphi_ri)
      naux_left = SIZE(sphi_li, 2)
      naux_right = SIZE(sphi_ri, 2)
      IF (.NOT. onsite_left) THEN
         CALL hfx_lri_build_sphi_matrix(rbas_lj, sphi_lj)
         naux_left = naux_left + SIZE(sphi_lj, 2)
      ELSE
         ALLOCATE (sphi_lj(0, 0))
      END IF
      IF (.NOT. onsite_right) THEN
         CALL hfx_lri_build_sphi_matrix(rbas_rj, sphi_rj)
         naux_right = naux_right + SIZE(sphi_rj, 2)
      ELSE
         ALLOCATE (sphi_rj(0, 0))
      END IF

      IF (naux_left /= pair_target%naux_pair .OR. naux_right /= pair_density%naux_pair) THEN
         CPABORT("Inconsistent auxiliary dimension in LRI-HFXk pair-pair density contraction")
      END IF

      ALLOCATE (sphi_left(nraw_left, naux_left), SOURCE=0.0_dp)
      sphi_left(1:ncart_li, 1:SIZE(sphi_li, 2)) = sphi_li
      IF (.NOT. onsite_left) THEN
         sphi_left(ncart_li + 1:nraw_left, SIZE(sphi_li, 2) + 1:naux_left) = sphi_lj
      END IF
      ALLOCATE (sphi_right(nraw_right, naux_right), SOURCE=0.0_dp)
      sphi_right(1:ncart_ri, 1:SIZE(sphi_ri, 2)) = sphi_ri
      IF (.NOT. onsite_right) THEN
         sphi_right(ncart_ri + 1:nraw_right, SIZE(sphi_ri, 2) + 1:naux_right) = sphi_rj
      END IF

      ALLOCATE (tmp(nraw_left, naux_right), SOURCE=0.0_dp)
      tmp = MATMUL(v_raw, sphi_right)
      ALLOCATE (v_sph(naux_left, naux_right), SOURCE=0.0_dp)
      v_sph = MATMUL(TRANSPOSE(sphi_left), tmp)
      DEALLOCATE (tmp, v_raw)

      CALL hfx_lri_flatten_pair_coefficients(pair_target, c_target)
      CALL hfx_lri_flatten_pair_coefficients(pair_density, c_density)

      ALLOCATE (d_vec(pair_density%nao_i*pair_density%nao_j), SOURCE=0.0_dp)
      irow = 0
      DO jao = 1, pair_density%nao_j
         DO iao = 1, pair_density%nao_i
            irow = irow + 1
            d_vec(irow) = rho_block(iao, jao)
         END DO
      END DO

      ALLOCATE (aux_density(naux_right), aux_target(naux_left), k_vec(SIZE(c_target, 1)), &
                SOURCE=0.0_dp)
      aux_density = MATMUL(TRANSPOSE(c_density), d_vec)
      aux_target = MATMUL(v_sph, aux_density)
      k_vec = MATMUL(c_target, aux_target)

      ALLOCATE (k_block(pair_target%nao_i, pair_target%nao_j), SOURCE=0.0_dp)
      irow = 0
      DO jao = 1, pair_target%nao_j
         DO iao = 1, pair_target%nao_i
            irow = irow + 1
            k_block(iao, jao) = k_vec(irow)
         END DO
      END DO

      DEALLOCATE (k_vec, aux_target, aux_density, d_vec, c_density, c_target, v_sph, sphi_right, &
                  sphi_left)
      DEALLOCATE (sphi_rj, sphi_ri, sphi_lj, sphi_li)
      CALL cp_libint_cleanup_2eri(lib2)

   END SUBROUTINE hfx_lri_contract_pair_pair_density_realv

! **************************************************************************************************
!> \brief Probes a two-pair C V C^T contraction using a real pair-pair auxiliary metric.
!> \param basis_set_RI ...
!> \param pair_left ...
!> \param pair_right ...
!> \param ri_data ...
!> \param v_norm ...
!> \param cvc_norm ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_contract_pair_pair_cvc_realv(basis_set_RI, pair_left, pair_right, ri_data, &
                                                   v_norm, cvc_norm)

      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair_left, pair_right
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      REAL(dp), INTENT(OUT)                              :: v_norm, cvc_norm

      INTEGER :: maxl_li, maxl_lj, maxl_ri, maxl_rj, naux_left, naux_right, ncart_li, ncart_lj, &
         ncart_ri, ncart_rj, nraw_left, nraw_right, offset_left_j, offset_right_j
      LOGICAL                                            :: onsite_left, onsite_right
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: c_left, c_right, cvc_block, sphi_left, &
                                                            sphi_li, sphi_lj, sphi_ri, sphi_right, &
                                                            sphi_rj, tmp, v_raw, v_sph
      REAL(dp), DIMENSION(3)                             :: r0
      TYPE(cp_libint_t)                                  :: lib2
      TYPE(gto_basis_set_type), POINTER                  :: rbas_li, rbas_lj, rbas_ri, rbas_rj

      NULLIFY (rbas_li, rbas_lj, rbas_ri, rbas_rj)
      v_norm = 0.0_dp
      cvc_norm = 0.0_dp

      IF (.NOT. ALLOCATED(pair_left%c_fit)) RETURN
      IF (.NOT. ALLOCATED(pair_right%c_fit)) RETURN
      IF (pair_left%iatom /= pair_right%iatom) RETURN

      rbas_li => basis_set_RI(pair_left%ikind)%gto_basis_set
      rbas_lj => basis_set_RI(pair_left%jkind)%gto_basis_set
      rbas_ri => basis_set_RI(pair_right%ikind)%gto_basis_set
      rbas_rj => basis_set_RI(pair_right%jkind)%gto_basis_set
      IF (.NOT. ASSOCIATED(rbas_li)) RETURN
      IF (.NOT. ASSOCIATED(rbas_lj)) RETURN
      IF (.NOT. ASSOCIATED(rbas_ri)) RETURN
      IF (.NOT. ASSOCIATED(rbas_rj)) RETURN

      onsite_left = (pair_left%iatom == pair_left%jatom .AND. NORM2(pair_left%rab) <= 1.0E-12_dp)
      onsite_right = (pair_right%iatom == pair_right%jatom .AND. NORM2(pair_right%rab) <= 1.0E-12_dp)

      CALL get_gto_basis_set(rbas_li, maxl=maxl_li)
      CALL get_gto_basis_set(rbas_lj, maxl=maxl_lj)
      CALL get_gto_basis_set(rbas_ri, maxl=maxl_ri)
      CALL get_gto_basis_set(rbas_rj, maxl=maxl_rj)
      CALL init_md_ftable(nmax=MAX(maxl_li + maxl_ri, maxl_li + maxl_rj, &
                                   maxl_lj + maxl_ri, maxl_lj + maxl_rj))
      CALL cp_libint_init_2eri(lib2, MAX(maxl_li, maxl_lj, maxl_ri, maxl_rj))
      CALL cp_libint_set_contrdepth(lib2, 1)

      ncart_li = hfx_lri_count_cart_functions(rbas_li)
      ncart_lj = hfx_lri_count_cart_functions(rbas_lj)
      ncart_ri = hfx_lri_count_cart_functions(rbas_ri)
      ncart_rj = hfx_lri_count_cart_functions(rbas_rj)
      nraw_left = ncart_li
      IF (.NOT. onsite_left) nraw_left = nraw_left + ncart_lj
      nraw_right = ncart_ri
      IF (.NOT. onsite_right) nraw_right = nraw_right + ncart_rj

      r0(:) = 0.0_dp
      ALLOCATE (v_raw(nraw_left, nraw_right), SOURCE=0.0_dp)
      CALL hfx_lri_fill_2c_block(rbas_li, rbas_ri, r0, r0, lib2, ri_data, v_raw, 0, 0)
      IF (.NOT. onsite_right) THEN
         offset_right_j = ncart_ri
         CALL hfx_lri_fill_2c_block(rbas_li, rbas_rj, r0, pair_right%rab, lib2, ri_data, v_raw, 0, &
                                    offset_right_j)
      END IF
      IF (.NOT. onsite_left) THEN
         offset_left_j = ncart_li
         CALL hfx_lri_fill_2c_block(rbas_lj, rbas_ri, pair_left%rab, r0, lib2, ri_data, v_raw, &
                                    offset_left_j, 0)
         IF (.NOT. onsite_right) THEN
            offset_right_j = ncart_ri
            CALL hfx_lri_fill_2c_block(rbas_lj, rbas_rj, pair_left%rab, pair_right%rab, lib2, &
                                       ri_data, v_raw, offset_left_j, offset_right_j)
         END IF
      END IF

      CALL hfx_lri_build_sphi_matrix(rbas_li, sphi_li)
      CALL hfx_lri_build_sphi_matrix(rbas_ri, sphi_ri)
      naux_left = SIZE(sphi_li, 2)
      naux_right = SIZE(sphi_ri, 2)
      IF (.NOT. onsite_left) THEN
         CALL hfx_lri_build_sphi_matrix(rbas_lj, sphi_lj)
         naux_left = naux_left + SIZE(sphi_lj, 2)
      ELSE
         ALLOCATE (sphi_lj(0, 0))
      END IF
      IF (.NOT. onsite_right) THEN
         CALL hfx_lri_build_sphi_matrix(rbas_rj, sphi_rj)
         naux_right = naux_right + SIZE(sphi_rj, 2)
      ELSE
         ALLOCATE (sphi_rj(0, 0))
      END IF

      IF (naux_left /= pair_left%naux_pair .OR. naux_right /= pair_right%naux_pair) THEN
         CPABORT("Inconsistent auxiliary dimension in LRI-HFXk pair-pair V probe")
      END IF

      ALLOCATE (sphi_left(nraw_left, naux_left), SOURCE=0.0_dp)
      sphi_left(1:ncart_li, 1:SIZE(sphi_li, 2)) = sphi_li
      IF (.NOT. onsite_left) THEN
         sphi_left(ncart_li + 1:nraw_left, SIZE(sphi_li, 2) + 1:naux_left) = sphi_lj
      END IF
      ALLOCATE (sphi_right(nraw_right, naux_right), SOURCE=0.0_dp)
      sphi_right(1:ncart_ri, 1:SIZE(sphi_ri, 2)) = sphi_ri
      IF (.NOT. onsite_right) THEN
         sphi_right(ncart_ri + 1:nraw_right, SIZE(sphi_ri, 2) + 1:naux_right) = sphi_rj
      END IF

      ALLOCATE (tmp(nraw_left, naux_right), SOURCE=0.0_dp)
      tmp = MATMUL(v_raw, sphi_right)
      ALLOCATE (v_sph(naux_left, naux_right), SOURCE=0.0_dp)
      v_sph = MATMUL(TRANSPOSE(sphi_left), tmp)
      DEALLOCATE (tmp, v_raw)

      CALL hfx_lri_flatten_pair_coefficients(pair_left, c_left)
      CALL hfx_lri_flatten_pair_coefficients(pair_right, c_right)
      ALLOCATE (tmp(SIZE(c_left, 1), naux_right), SOURCE=0.0_dp)
      tmp = MATMUL(c_left, v_sph)
      ALLOCATE (cvc_block(SIZE(c_left, 1), SIZE(c_right, 1)), SOURCE=0.0_dp)
      cvc_block = MATMUL(tmp, TRANSPOSE(c_right))
      v_norm = SQRT(SUM(v_sph*v_sph))
      cvc_norm = SQRT(SUM(cvc_block*cvc_block))

      DEALLOCATE (cvc_block, tmp, c_right, c_left, v_sph, sphi_right, sphi_left)
      DEALLOCATE (sphi_rj, sphi_ri, sphi_lj, sphi_li)
      CALL cp_libint_cleanup_2eri(lib2)

   END SUBROUTINE hfx_lri_contract_pair_pair_cvc_realv

! **************************************************************************************************
!> \brief Flattens C(mu,nu,P) into C(mu nu,P).
!> \param pair_block ...
!> \param c_mat ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_flatten_pair_coefficients(pair_block, c_mat)

      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair_block
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: c_mat

      INTEGER                                            :: iao, irow, jao

      IF (.NOT. ALLOCATED(pair_block%c_fit)) THEN
         ALLOCATE (c_mat(0, 0))
         RETURN
      END IF

      ALLOCATE (c_mat(pair_block%nao_i*pair_block%nao_j, pair_block%naux_pair), SOURCE=0.0_dp)
      irow = 0
      DO jao = 1, pair_block%nao_j
         DO iao = 1, pair_block%nao_i
            irow = irow + 1
            c_mat(irow, :) = pair_block%c_fit(iao, jao, :)
         END DO
      END DO

   END SUBROUTINE hfx_lri_flatten_pair_coefficients

! **************************************************************************************************
!> \brief Contracts a stored self-pair C/V block with a real density block.
!> \param pair_block ...
!> \param rho_block ...
!> \param k_block ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_contract_real_density_self(pair_block, rho_block, k_block)

      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair_block
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: rho_block
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: k_block

      INTEGER                                            :: iao, jao, jaux, nao_i, nao_j, naux
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: aux_density, v_aux

      IF (.NOT. ALLOCATED(pair_block%c_fit)) THEN
         CPABORT("Missing C block in LRI-HFXk density contraction probe")
      END IF
      IF (.NOT. ALLOCATED(pair_block%v_fit)) THEN
         CPABORT("Missing V block in LRI-HFXk density contraction probe")
      END IF

      nao_i = pair_block%nao_i
      nao_j = pair_block%nao_j
      naux = pair_block%naux_pair
      IF (SIZE(rho_block, 1) < nao_i .OR. SIZE(rho_block, 2) < nao_j) THEN
         CPABORT("Inconsistent density block dimension in LRI-HFXk contraction probe")
      END IF

      ALLOCATE (aux_density(naux), v_aux(naux), k_block(nao_i, nao_j), SOURCE=0.0_dp)

      DO jaux = 1, naux
         DO jao = 1, nao_j
            DO iao = 1, nao_i
               aux_density(jaux) = aux_density(jaux) + &
                                   pair_block%c_fit(iao, jao, jaux)*rho_block(iao, jao)
            END DO
         END DO
      END DO

      v_aux = MATMUL(pair_block%v_fit, aux_density)

      DO jao = 1, nao_j
         DO iao = 1, nao_i
            k_block(iao, jao) = SUM(pair_block%c_fit(iao, jao, :)*v_aux(:))
         END DO
      END DO

      DEALLOCATE (v_aux, aux_density)

   END SUBROUTINE hfx_lri_contract_real_density_self

! **************************************************************************************************
!> \brief Finds a stored pair-local LRI coefficient block by atom pair and cell.
!> \param pair_blocks ...
!> \param npair_blocks ...
!> \param iatom ...
!> \param jatom ...
!> \param cell ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION hfx_lri_find_pair_block(pair_blocks, npair_blocks, iatom, jatom, cell) RESULT(ipair_found)

      TYPE(hfx_lri_pair_block_type), DIMENSION(:), &
         INTENT(IN)                                      :: pair_blocks
      INTEGER, INTENT(IN)                                :: npair_blocks, iatom, jatom
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell

      INTEGER                                            :: ipair

      ipair_found = 0
      DO ipair = 1, MIN(npair_blocks, SIZE(pair_blocks))
         IF (pair_blocks(ipair)%iatom == iatom .AND. pair_blocks(ipair)%jatom == jatom .AND. &
             ALL(pair_blocks(ipair)%cell == cell)) THEN
            ipair_found = ipair
            RETURN
         END IF
      END DO

   END FUNCTION hfx_lri_find_pair_block

END MODULE hfx_lri_kp
