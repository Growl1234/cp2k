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
!$ USE OMP_LIB, ONLY: omp_get_max_threads, &
!$                    omp_get_thread_num
   USE basis_set_types, ONLY: get_gto_basis_set, &
                              gto_basis_set_p_type, &
                              gto_basis_set_type
   USE coulomb_integral_interface, ONLY: coulomb_cutoff_screen_factor, &
                                         compute_coulomb_2c, &
                                         compute_coulomb_3c, &
                                         coulomb_integral_cleanup, &
                                         coulomb_integral_context_type, &
                                         coulomb_integral_init
   USE cp_control_types, ONLY: dft_control_type
   USE cp_dbcsr_api, ONLY: dbcsr_get_block_p, &
                           dbcsr_p_type
   USE gamma, ONLY: init_md_ftable
   USE hfx_lri_contract, ONLY: hfx_lri_contract_exchange, &
                               hfx_lri_contract_exchange_add, &
                               hfx_lri_contract_outer_accumulate, &
                               hfx_lri_contract_outer_finalize, &
                               hfx_lri_contract_outer_prepare, &
                               hfx_lri_density_cell, &
                               hfx_lri_product_index, &
                               hfx_lri_right_cell
   USE hfx_lri_types, ONLY: hfx_lri_kp_release, &
                            hfx_lri_exchange_workspace_type, &
                            hfx_lri_metric_cache_entry_type, &
                            hfx_lri_metric_cache_type, &
                            hfx_lri_metric_table_type, &
                            hfx_lri_pair_block_type
   USE hfx_types, ONLY: hfx_ri_type
   USE input_constants, ONLY: do_potential_short
   USE integral_library_types, ONLY: active_integral_library, &
                                     library_libint
   USE kinds, ONLY: dp, &
                    int_8
   USE kpoint_types, ONLY: get_kpoint_info, &
                           kpoint_type
   USE message_passing, ONLY: mp_para_env_type
   USE orbital_pointers, ONLY: ncoset
   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   USE qs_integral_utils, ONLY: basis_set_list_setup
   USE qs_interactions, ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types, ONLY: qs_kind_type
   USE qs_neighbor_list_types, ONLY: get_iterator_info, &
                                     neighbor_list_iterate, &
                                     neighbor_list_iterator_create, &
                                     neighbor_list_iterator_p_type, &
                                     neighbor_list_iterator_release, &
                                     neighbor_list_set_p_type, &
                                     release_neighbor_list_sets
   USE qs_neighbor_lists, ONLY: setup_neighbor_list
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hfx_lri_update_ks_kp

   INTEGER, PARAMETER :: lri_debug_level = 0
   INTEGER, PARAMETER :: lri_metric_cache_capacity = 128
   INTEGER, PARAMETER :: lri_metric_cache_ways = 4
   REAL(dp), PARAMETER :: lri_metric_cache_position_tol = 1.0E-12_dp
   REAL(dp), PARAMETER :: lri_debug_rho_threshold = 1.0E-8_dp
   INTEGER, PARAMETER :: lri_debug_max_general_targets = 20
   INTEGER, PARAMETER :: lri_debug_max_general_tests = 200
   REAL(dp), PARAMETER :: lri_debug_k_threshold_low = 1.0E-8_dp
   REAL(dp), PARAMETER :: lri_debug_k_threshold_mid = 1.0E-6_dp
   REAL(dp), PARAMETER :: lri_debug_k_threshold_high = 1.0E-4_dp
   REAL(dp), PARAMETER :: lri_debug_k_threshold_large = 1.0E-2_dp

   TYPE lri_pair_lookup_type
      INTEGER                                            :: capacity = 0
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_count, atom_offset, atom_pairs, slot
   END TYPE lri_pair_lookup_type

   TYPE lri_density_block_type
      INTEGER                                            :: iatom = 0, ikind = 0, jatom = 0, jkind = 0
      INTEGER, DIMENSION(3)                              :: cell = 0
      REAL(dp)                                           :: norm = 0.0_dp
      REAL(dp), DIMENSION(3)                             :: rab = 0.0_dp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: values
   END TYPE lri_density_block_type

   TYPE lri_density_cache_type
      INTEGER                                            :: capacity = 0, nblocks = 0
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: slot
      TYPE(lri_density_block_type), ALLOCATABLE, DIMENSION(:) :: blocks
   END TYPE lri_density_cache_type

   TYPE lri_density_ref_type
      INTEGER                                            :: iblock = 0, latom = 0
   END TYPE lri_density_ref_type

   TYPE lri_density_index_type
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_count, atom_offset
      TYPE(lri_density_ref_type), ALLOCATABLE, DIMENSION(:) :: refs
   END TYPE lri_density_index_type

   TYPE lri_basis_transform_type
      REAL(dp)                                           :: radius = 0.0_dp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: sphi
   END TYPE lri_basis_transform_type

CONTAINS

! **************************************************************************************************
!> \brief Moves all allocatable buffers between persistent and thread-private LRI workspaces.
!> \param source ...
!> \param destination ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_exchange_workspace_move(source, destination)

      TYPE(hfx_lri_exchange_workspace_type), &
         INTENT(INOUT)                                   :: source, destination

      destination%pair_screened = source%pair_screened
      destination%metric_cache%clock = source%metric_cache%clock
      destination%metric_cache%hits = source%metric_cache%hits
      destination%metric_cache%misses = source%metric_cache%misses
      destination%metric_cache%screened = source%metric_cache%screened
      CALL MOVE_ALLOC(source%metric_cache%entries, destination%metric_cache%entries)
      CALL MOVE_ALLOC(source%v_outer, destination%v_outer)
      CALL MOVE_ALLOC(source%v_sph, destination%v_sph)
      CALL MOVE_ALLOC(source%contract%density_right, destination%contract%density_right)
      CALL MOVE_ALLOC(source%contract%left_coeff, destination%contract%left_coeff)
      CALL MOVE_ALLOC(source%contract%right_coeff, destination%contract%right_coeff)
      CALL MOVE_ALLOC(source%contract%weighted_left, destination%contract%weighted_left)
      CALL MOVE_ALLOC(source%contract%weighted_right_flat, &
                      destination%contract%weighted_right_flat)
      CALL MOVE_ALLOC(source%contract%left_flat, destination%contract%left_flat)
      CALL MOVE_ALLOC(source%contract%left_metric_flat, destination%contract%left_metric_flat)
      CALL MOVE_ALLOC(source%contract%outer_density_sum, destination%contract%outer_density_sum)

      source%pair_screened = 0_int_8
      source%metric_cache%clock = 0_int_8
      source%metric_cache%hits = 0_int_8
      source%metric_cache%misses = 0_int_8
      source%metric_cache%screened = 0_int_8

   END SUBROUTINE hfx_lri_exchange_workspace_move

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

      TYPE(mp_para_env_type), POINTER                    :: para_env

      ehfx = 0.0_dp
      NULLIFY (para_env)

      CALL get_qs_env(qs_env=qs_env, para_env=para_env)

      IF (ri_data%hfx_pot%potential_type /= do_potential_short) THEN
         CPABORT("LRI-HFXk is currently restricted to short-range/HSE exchange")
      END IF
      IF (active_integral_library%coulomb2_library /= library_libint .OR. &
          active_integral_library%coulomb3_library /= library_libint) THEN
         CPABORT("LRI-HFXk currently requires LIBINT for COULOMB_2C and COULOMB_3C")
      END IF
      IF (nspins /= 1) THEN
         CPABORT("LRI-HFXk currently supports only spin-unpolarized calculations")
      END IF
      IF (.NOT. ASSOCIATED(para_env)) THEN
         CPABORT("LRI-HFXk could not obtain the MPI environment")
      END IF
      IF (para_env%num_pe /= 1) THEN
         CPABORT("LRI-HFXk distributed pair exchange is not implemented yet; use one MPI rank")
      END IF

      IF (geometry_did_change .OR. .NOT. ASSOCIATED(ri_data%lri_kp)) THEN
         CALL hfx_lri_pre_scf_kp(qs_env, ri_data)
      ELSE IF (.NOT. ri_data%lri_kp%initialized) THEN
         CALL hfx_lri_pre_scf_kp(qs_env, ri_data)
      END IF

      CALL hfx_lri_build_rank_local_k_blocks(qs_env, ri_data, rho_ao)
      CALL hfx_lri_compute_rank_local_energy(qs_env, ri_data, hf_fraction, ehfx)
      CALL hfx_lri_write_rank_local_k_blocks(qs_env, ri_data, ks_matrix, hf_fraction)

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
         ipair_block, ithread, jatom, jkind, jneighbor, maxl_ao, maxl_ao_kind, maxl_ri, &
         maxl_ri_kind, mepos, mmax, nao_i, nao_j, naux_i, naux_j, naux_pair, nkind, nlist, &
         nneighbor, npair_blocks_local, npair_blocks_max, npair_blocks_sum, npair_candidates, &
         npairs_local
      INTEGER :: npairs_max, npairs_sum, nthreads
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
         dbg_sph_c_norm, dbg_sph_cvc_norm, dbg_sph_j_norm, dbg_sph_m_norm
      REAL(dp), DIMENSION(3)                             :: dbg_rab, rab
      TYPE(coulomb_integral_context_type)                :: integral_context
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_AO, basis_set_RI
      TYPE(gto_basis_set_type), POINTER                  :: obas_i, obas_j, rbas_i, rbas_j
      TYPE(hfx_lri_exchange_workspace_type)              :: exchange_workspace
      TYPE(hfx_lri_pair_block_type)                      :: dbg_pair_block
      TYPE(lri_basis_transform_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: ao_transform_cache, ri_transform_cache
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
      dbg_pair_pair_index = 0
      dbg_pair_pair_r2 = HUGE(1.0_dp)
      dbg_pair_pair_v12_norm = 0.0_dp
      dbg_pair_pair_cvc12_norm = 0.0_dp

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, nkind=nkind, para_env=para_env, &
                      dft_control=dft_control)

      CALL hfx_lri_kp_release(ri_data%lri_kp)
      ALLOCATE (ri_data%lri_kp)
      nthreads = 1
!$    nthreads = omp_get_max_threads()
      ALLOCATE (ri_data%lri_kp%thread_workspaces(nthreads))

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

      CALL hfx_lri_build_ri_transform_cache(basis_set_AO, ao_transform_cache)
      CALL hfx_lri_build_ri_transform_cache(basis_set_RI, ri_transform_cache)
      maxl_ao = 0
      maxl_ri = 0
      DO ibasis = 1, SIZE(basis_set_AO)
         obas_i => basis_set_AO(ibasis)%gto_basis_set
         rbas_i => basis_set_RI(ibasis)%gto_basis_set
         IF (ASSOCIATED(obas_i)) THEN
            CALL get_gto_basis_set(obas_i, maxl=maxl_ao_kind)
            maxl_ao = MAX(maxl_ao, maxl_ao_kind)
         END IF
         IF (ASSOCIATED(rbas_i)) THEN
            CALL get_gto_basis_set(rbas_i, maxl=maxl_ri_kind)
            maxl_ri = MAX(maxl_ri, maxl_ri_kind)
         END IF
      END DO
      mmax = MAX(2*maxl_ri, 2*maxl_ao + maxl_ri)
      CALL init_md_ftable(nmax=mmax)

      ! First prototype: AO-AO neighbor list only.
      ! Later this cutoff should be tightened for HSE short-range LRI.
      CALL setup_neighbor_list(soo_list, basis_set_AO, qs_env=qs_env, mic=.FALSE., symmetric=.TRUE.)

      ! Count first so that every significant pair is retained. The previous fixed-size diagnostic
      ! buffer silently truncated larger systems and made both the energy and the workload depend on
      ! an implementation constant.
      CALL neighbor_list_iterator_create(nl_iterator, soo_list)

      npair_candidates = 0
      mepos = 0
      DO WHILE (neighbor_list_iterate(nl_iterator, mepos) == 0)
         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind)
         obas_i => basis_set_AO(ikind)%gto_basis_set
         obas_j => basis_set_AO(jkind)%gto_basis_set
         rbas_i => basis_set_RI(ikind)%gto_basis_set
         rbas_j => basis_set_RI(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(obas_i)) CYCLE
         IF (.NOT. ASSOCIATED(obas_j)) CYCLE
         IF (.NOT. ASSOCIATED(rbas_i)) CYCLE
         IF (.NOT. ASSOCIATED(rbas_j)) CYCLE
         npair_candidates = npair_candidates + 1
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      ALLOCATE (ri_data%lri_kp%pair_blocks(npair_candidates))
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

         IF (lri_debug_level >= 1) THEN
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
         END IF

         c_nreal = INT(nao_i, int_8)*INT(nao_j, int_8)*INT(naux_pair, int_8)

         ipair_block = ri_data%lri_kp%npair_blocks + 1
         CALL hfx_lri_init_pair_block(ri_data%lri_kp%pair_blocks(ipair_block), iatom, jatom, &
                                      ikind, jkind, cell_j, rab, nao_i, nao_j, naux_i, naux_j, &
                                      naux_pair)
         ri_data%lri_kp%npair_blocks = ipair_block

         ri_data%lri_kp%npairs = ri_data%lri_kp%npairs + 1
         ri_data%lri_kp%sum_nao_pair = ri_data%lri_kp%sum_nao_pair + INT(nao_i, int_8)*INT(nao_j, int_8)
         ri_data%lri_kp%sum_naux_pair = ri_data%lri_kp%sum_naux_pair + INT(naux_pair, int_8)
         ri_data%lri_kp%c_nreal_est = ri_data%lri_kp%c_nreal_est + c_nreal

      END DO

      CALL neighbor_list_iterator_release(nl_iterator)
      CALL release_neighbor_list_sets(soo_list)

      ! The neighbor-list iterator is serial, but every fitted pair is independent. Reuse the
      ! immutable basis transforms and give every worker its own Libint engines.
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP          SHARED(ao_transform_cache, basis_set_AO, basis_set_RI, maxl_ao, maxl_ri, &
      !$OMP                 ri_data, ri_transform_cache) &
      !$OMP          PRIVATE(exchange_workspace, integral_context, ipair_block, ithread)
      ithread = 1
!$    ithread = omp_get_thread_num() + 1
      CALL hfx_lri_exchange_workspace_move(ri_data%lri_kp%thread_workspaces(ithread), &
                                           exchange_workspace)
      CALL coulomb_integral_init(integral_context, max_am_2c=maxl_ri, &
                                 max_am_3c=MAX(maxl_ao, maxl_ri))
      !$OMP DO SCHEDULE(DYNAMIC, 1)
      DO ipair_block = 1, ri_data%lri_kp%npair_blocks
         CALL hfx_lri_build_pair_coefficients( &
            basis_set_AO, basis_set_RI, ri_data%lri_kp%pair_blocks(ipair_block), ri_data, &
            integral_context, ao_transform_cache, ri_transform_cache, &
            exchange_workspace%metric_cache)
      END DO
      !$OMP END DO
      CALL coulomb_integral_cleanup(integral_context)
      CALL hfx_lri_exchange_workspace_move(exchange_workspace, &
                                           ri_data%lri_kp%thread_workspaces(ithread))
      !$OMP END PARALLEL

      ri_data%lri_kp%c_nreal_stored = 0_int_8
      npair_blocks_local = 0
      DO ipair_block = 1, ri_data%lri_kp%npair_blocks
         IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks(ipair_block)%c_fit)) CYCLE
         npair_blocks_local = npair_blocks_local + 1
         ri_data%lri_kp%c_nreal_stored = ri_data%lri_kp%c_nreal_stored + &
                                         INT(SIZE(ri_data%lri_kp%pair_blocks(ipair_block)%c_fit), int_8)
      END DO

      IF (dbg_pair_found .AND. lri_debug_level >= 1) THEN
         CALL hfx_lri_init_pair_block(dbg_pair_block, dbg_iatom, dbg_jatom, dbg_ikind, dbg_jkind, &
                                      dbg_cell, dbg_rab, dbg_nao_i, dbg_nao_j, dbg_naux_i, dbg_naux_j, &
                                      dbg_naux_pair)
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

      IF (lri_debug_level >= 1) THEN
         IF (dbg_pair_found .AND. ri_data%lri_kp%npair_blocks > 0) THEN
            dbg_stored_index = hfx_lri_find_pair_block(ri_data%lri_kp%pair_blocks, &
                                                       ri_data%lri_kp%npair_blocks, &
                                                       dbg_iatom, dbg_jatom, dbg_cell)
         END IF

         IF (dbg_stored_index > 0) THEN
            DO ipair_block = 1, ri_data%lri_kp%npair_blocks
               IF (ipair_block == dbg_stored_index) CYCLE
               IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks(ipair_block)%c_fit)) CYCLE
               ! Pick the nearest stored pair for the real V12 probe. The actual
               ! auxiliary-center positions are now taken from pair-block metadata.
               IF (SUM((ri_data%lri_kp%pair_blocks(ipair_block)%rab - &
                        ri_data%lri_kp%pair_blocks(dbg_stored_index)%rab)**2) >= dbg_pair_pair_r2) CYCLE
               dbg_pair_pair_index = ipair_block
               dbg_pair_pair_r2 = SUM((ri_data%lri_kp%pair_blocks(ipair_block)%rab - &
                                       ri_data%lri_kp%pair_blocks(dbg_stored_index)%rab)**2)
            END DO
            IF (dbg_pair_pair_index > 0) THEN
               CALL hfx_lri_contract_pair_pair_cvc_realv( &
                  basis_set_RI, ri_data%lri_kp%pair_blocks(dbg_stored_index), &
                  ri_data%lri_kp%pair_blocks(dbg_pair_pair_index), ri_data, &
                  dbg_pair_pair_v12_norm, dbg_pair_pair_cvc12_norm)
            END IF
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
         IF (dbg_pair_found .AND. lri_debug_level >= 1) THEN
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
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10)") &
               "aux centers", dbg_pair_block%naux_centers
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
      TYPE(coulomb_integral_context_type)                :: integral_context
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

      CALL coulomb_integral_init(integral_context, max_am_2c=MAX(maxl_ri, maxl_rj), &
                                 max_am_3c=MAX(maxl_ai, maxl_aj, maxl_ri, maxl_rj))

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
         CALL hfx_lri_accumulate_2c_norm(rbas_i, rbas_i, r0, r0, integral_context, ri_data, metric_ncart, &
                                         metric_norm2)
         IF (.NOT. onsite_aux) THEN
            CALL hfx_lri_accumulate_2c_norm(rbas_i, rbas_j, r0, rab, integral_context, ri_data, metric_ncart, &
                                            metric_norm2, prefactor=2.0_dp)
            CALL hfx_lri_accumulate_2c_norm(rbas_j, rbas_j, rab, rab, integral_context, ri_data, metric_ncart, &
                                            metric_norm2)
         END IF

         CALL hfx_lri_accumulate_3c_norm(obas_i, obas_j, rbas_i, r0, rab, r0, integral_context, ri_data, &
                                         threec_i_ncart, threec_i_norm2)
         IF (.NOT. onsite_aux) THEN
            CALL hfx_lri_accumulate_3c_norm(obas_i, obas_j, rbas_j, r0, rab, rab, integral_context, ri_data, &
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

         CALL hfx_lri_fill_2c_block(rbas_i, rbas_i, r0, r0, integral_context, ri_data, metric_block, 0, 0)
         CALL hfx_lri_fill_3c_block(obas_i, obas_j, rbas_i, r0, rab, r0, integral_context, ri_data, &
                                    j_block, 0, ncart_aj)
         IF (.NOT. onsite_aux) THEN
            offset_j = ncart_ri
            CALL hfx_lri_fill_2c_block(rbas_i, rbas_j, r0, rab, integral_context, ri_data, metric_block, 0, &
                                       offset_j, symmetric=.TRUE.)
            CALL hfx_lri_fill_2c_block(rbas_j, rbas_j, rab, rab, integral_context, ri_data, metric_block, &
                                       offset_j, offset_j)
            CALL hfx_lri_fill_3c_block(obas_i, obas_j, rbas_j, r0, rab, rab, integral_context, ri_data, &
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

            CALL hfx_lri_probe_metric_inverse(metric_block, j_block, ri_data%eps_eigval, &
                                              metric_eig_min, metric_eig_max, &
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
            CALL hfx_lri_probe_metric_inverse(metric_sph, j_sph, ri_data%eps_eigval, &
                                              sph_eig_min, sph_eig_max, &
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

      CALL coulomb_integral_cleanup(integral_context)

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
!> \param integral_context ...
!> \param ri_data ...
!> \param ncart_total ...
!> \param norm2_total ...
!> \param prefactor ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_accumulate_2c_norm(basis_a, basis_b, ra, rb, integral_context, ri_data, ncart_total, &
                                         norm2_total, prefactor)

      TYPE(gto_basis_set_type), POINTER                  :: basis_a, basis_b
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rb
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: integral_context
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      INTEGER, INTENT(INOUT)                             :: ncart_total
      REAL(dp), INTENT(INOUT)                            :: norm2_total
      REAL(dp), INTENT(IN), OPTIONAL                     :: prefactor

      INTEGER                                            :: iset, jset, ncoa, ncob, nseta, nsetb
      INTEGER, DIMENSION(:), POINTER                     :: lmax_a, lmax_b, lmin_a, lmin_b, npgf_a, &
                                                            npgf_b
      REAL(dp)                                           :: pfac
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: sij
      REAL(dp), DIMENSION(:, :), POINTER                 :: rpgf_a, rpgf_b, zeta_a, zeta_b

      CALL get_gto_basis_set(basis_a, lmax=lmax_a, lmin=lmin_a, npgf=npgf_a, nset=nseta, &
                             pgf_radius=rpgf_a, zet=zeta_a)
      CALL get_gto_basis_set(basis_b, lmax=lmax_b, lmin=lmin_b, npgf=npgf_b, nset=nsetb, &
                             pgf_radius=rpgf_b, zet=zeta_b)

      pfac = 1.0_dp
      IF (PRESENT(prefactor)) pfac = prefactor

      DO iset = 1, nseta
         ncoa = npgf_a(iset)*ncoset(lmax_a(iset))
         DO jset = 1, nsetb
            ncob = npgf_b(jset)*ncoset(lmax_b(jset))
            IF (ncoa*ncob == 0) CYCLE
            ALLOCATE (sij(ncoa, ncob), SOURCE=0.0_dp)
            CALL compute_coulomb_2c(integral_context, lmin_a(iset), lmax_a(iset), &
                                    lmin_b(jset), lmax_b(jset), npgf_a(iset), npgf_b(jset), &
                                    zeta_a(:, iset), zeta_b(:, jset), rpgf_a(:, iset), &
                                    rpgf_b(:, jset), ra, rb, sij, ri_data%hfx_pot)
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
!> \param integral_context ...
!> \param ri_data ...
!> \param ncart_total ...
!> \param norm2_total ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_accumulate_3c_norm(basis_a, basis_b, basis_c, ra, rb, rc, integral_context, ri_data, &
                                         ncart_total, norm2_total)

      TYPE(gto_basis_set_type), POINTER                  :: basis_a, basis_b, basis_c
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rb, rc
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: integral_context
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      INTEGER, INTENT(INOUT)                             :: ncart_total
      REAL(dp), INTENT(INOUT)                            :: norm2_total

      INTEGER                                            :: iset, jset, kset, ncoa, ncob, ncoc, &
                                                            nseta, nsetb, nsetc
      INTEGER, DIMENSION(:), POINTER                     :: lmax_a, lmax_b, lmax_c, lmin_a, lmin_b, &
                                                            lmin_c, npgf_a, npgf_b, npgf_c
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: sijk
      REAL(dp), DIMENSION(:, :), POINTER                 :: rpgf_a, rpgf_b, rpgf_c, zeta_a, zeta_b, &
                                                            zeta_c

      CALL get_gto_basis_set(basis_a, lmax=lmax_a, lmin=lmin_a, npgf=npgf_a, nset=nseta, &
                             pgf_radius=rpgf_a, zet=zeta_a)
      CALL get_gto_basis_set(basis_b, lmax=lmax_b, lmin=lmin_b, npgf=npgf_b, nset=nsetb, &
                             pgf_radius=rpgf_b, zet=zeta_b)
      CALL get_gto_basis_set(basis_c, lmax=lmax_c, lmin=lmin_c, npgf=npgf_c, nset=nsetc, &
                             pgf_radius=rpgf_c, zet=zeta_c)

      DO iset = 1, nseta
         ncoa = npgf_a(iset)*ncoset(lmax_a(iset))
         DO jset = 1, nsetb
            ncob = npgf_b(jset)*ncoset(lmax_b(jset))
            DO kset = 1, nsetc
               ncoc = npgf_c(kset)*ncoset(lmax_c(kset))
               IF (ncoa*ncob*ncoc == 0) CYCLE
               ALLOCATE (sijk(ncoa, ncob, ncoc), SOURCE=0.0_dp)
               CALL compute_coulomb_3c(integral_context, lmin_a(iset), lmax_a(iset), &
                                       lmin_b(jset), lmax_b(jset), lmin_c(kset), lmax_c(kset), &
                                       npgf_a(iset), npgf_b(jset), npgf_c(kset), zeta_a(:, iset), &
                                       zeta_b(:, jset), zeta_c(:, kset), rpgf_a(:, iset), &
                                       rpgf_b(:, jset), rpgf_c(:, kset), ra, rb, rc, sijk, &
                                       ri_data%hfx_pot)
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
!> \param integral_context ...
!> \param ri_data ...
!> \param mat ...
!> \param offset_a ...
!> \param offset_b ...
!> \param symmetric ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_fill_2c_block(basis_a, basis_b, ra, rb, integral_context, ri_data, mat, offset_a, &
                                    offset_b, symmetric)

      TYPE(gto_basis_set_type), POINTER                  :: basis_a, basis_b
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rb
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: integral_context
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: mat
      INTEGER, INTENT(IN)                                :: offset_a, offset_b
      LOGICAL, INTENT(IN), OPTIONAL                      :: symmetric

      INTEGER                                            :: ia0, ib0, iset, jset, ncoa, ncob, nseta, &
                                                            nsetb
      INTEGER, DIMENSION(:), POINTER                     :: lmax_a, lmax_b, lmin_a, lmin_b, npgf_a, &
                                                            npgf_b
      LOGICAL                                            :: do_sym
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: sij
      REAL(dp), DIMENSION(:, :), POINTER                 :: rpgf_a, rpgf_b, zeta_a, zeta_b

      CALL get_gto_basis_set(basis_a, lmax=lmax_a, lmin=lmin_a, npgf=npgf_a, nset=nseta, &
                             pgf_radius=rpgf_a, zet=zeta_a)
      CALL get_gto_basis_set(basis_b, lmax=lmax_b, lmin=lmin_b, npgf=npgf_b, nset=nsetb, &
                             pgf_radius=rpgf_b, zet=zeta_b)

      do_sym = .FALSE.
      IF (PRESENT(symmetric)) do_sym = symmetric
      ia0 = offset_a
      DO iset = 1, nseta
         ncoa = npgf_a(iset)*ncoset(lmax_a(iset))
         ib0 = offset_b
         DO jset = 1, nsetb
            ncob = npgf_b(jset)*ncoset(lmax_b(jset))
            IF (ncoa*ncob > 0) THEN
               ALLOCATE (sij(ncoa, ncob), SOURCE=0.0_dp)
               CALL compute_coulomb_2c(integral_context, lmin_a(iset), lmax_a(iset), &
                                       lmin_b(jset), lmax_b(jset), npgf_a(iset), npgf_b(jset), &
                                       zeta_a(:, iset), zeta_b(:, jset), rpgf_a(:, iset), &
                                       rpgf_b(:, jset), ra, rb, sij, ri_data%hfx_pot)
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
!> \param integral_context ...
!> \param ri_data ...
!> \param mat ...
!> \param offset_c ...
!> \param ncart_b ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_fill_3c_block(basis_a, basis_b, basis_c, ra, rb, rc, integral_context, ri_data, mat, &
                                    offset_c, ncart_b)

      TYPE(gto_basis_set_type), POINTER                  :: basis_a, basis_b, basis_c
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rb, rc
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: integral_context
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: mat
      INTEGER, INTENT(IN)                                :: offset_c, ncart_b

      INTEGER                                            :: ia, ia0, ib, ib0, ic0, iset, jset, kset, &
                                                            ncoa, ncob, ncoc, nseta, nsetb, nsetc
      INTEGER, DIMENSION(:), POINTER                     :: lmax_a, lmax_b, lmax_c, lmin_a, lmin_b, &
                                                            lmin_c, npgf_a, npgf_b, npgf_c
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: sijk
      REAL(dp), DIMENSION(:, :), POINTER                 :: rpgf_a, rpgf_b, rpgf_c, zeta_a, zeta_b, &
                                                            zeta_c

      CALL get_gto_basis_set(basis_a, lmax=lmax_a, lmin=lmin_a, npgf=npgf_a, nset=nseta, &
                             pgf_radius=rpgf_a, zet=zeta_a)
      CALL get_gto_basis_set(basis_b, lmax=lmax_b, lmin=lmin_b, npgf=npgf_b, nset=nsetb, &
                             pgf_radius=rpgf_b, zet=zeta_b)
      CALL get_gto_basis_set(basis_c, lmax=lmax_c, lmin=lmin_c, npgf=npgf_c, nset=nsetc, &
                             pgf_radius=rpgf_c, zet=zeta_c)

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
                  CALL compute_coulomb_3c(integral_context, lmin_a(iset), lmax_a(iset), &
                                          lmin_b(jset), lmax_b(jset), lmin_c(kset), lmax_c(kset), &
                                          npgf_a(iset), npgf_b(jset), npgf_c(kset), zeta_a(:, iset), &
                                          zeta_b(:, jset), zeta_c(:, kset), rpgf_a(:, iset), &
                                          rpgf_b(:, jset), rpgf_c(:, kset), ra, rb, rc, sijk, &
                                          ri_data%hfx_pot)
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
!> \brief Builds each basis-kind Cartesian-to-spherical transformation once.
!> \param basis_set_RI ...
!> \param transform_cache ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_build_ri_transform_cache(basis_set_RI, transform_cache)

      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(lri_basis_transform_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(OUT)                       :: transform_cache

      INTEGER                                            :: ikind
      TYPE(gto_basis_set_type), POINTER                  :: basis

      NULLIFY (basis)
      ALLOCATE (transform_cache(SIZE(basis_set_RI)))
      DO ikind = 1, SIZE(basis_set_RI)
         basis => basis_set_RI(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis)) CYCLE
         CALL get_gto_basis_set(basis, kind_radius=transform_cache(ikind)%radius)
         CALL hfx_lri_build_sphi_matrix(basis, transform_cache(ikind)%sphi)
      END DO

   END SUBROUTINE hfx_lri_build_ri_transform_cache

! **************************************************************************************************
!> \brief Applies the two AO transformations to a three-center block without forming their
!>        Kronecker product.
!> \param j_cart ...
!> \param sphi_a ...
!> \param sphi_b ...
!> \param j_sph ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_transform_ao_pair(j_cart, sphi_a, sphi_b, j_sph)

      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: j_cart, sphi_a, sphi_b
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: j_sph

      INTEGER                                            :: iaux, naux, ncart_a, ncart_b, nsgf_a, &
                                                            nsgf_b
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: cart_matrix, sph_matrix, tmp_b

      ncart_a = SIZE(sphi_a, 1)
      ncart_b = SIZE(sphi_b, 1)
      nsgf_a = SIZE(sphi_a, 2)
      nsgf_b = SIZE(sphi_b, 2)
      naux = SIZE(j_cart, 2)
      IF (SIZE(j_cart, 1) /= ncart_a*ncart_b) THEN
         CPABORT("Inconsistent AO-product dimension in LRI-HFXk three-center transformation")
      END IF

      ALLOCATE (cart_matrix(ncart_b, ncart_a))
      ALLOCATE (tmp_b(nsgf_b, ncart_a))
      ALLOCATE (sph_matrix(nsgf_b, nsgf_a))
      ALLOCATE (j_sph(nsgf_a*nsgf_b, naux))
      DO iaux = 1, naux
         ! hfx_lri_fill_3c_block stores the second AO index contiguously.
         cart_matrix(:, :) = RESHAPE(j_cart(:, iaux), [ncart_b, ncart_a])
         tmp_b(:, :) = MATMUL(TRANSPOSE(sphi_b), cart_matrix)
         sph_matrix(:, :) = MATMUL(tmp_b, sphi_a)
         ! The resulting row order is (a-1)*nsgf_b+b, matching hfx_lri_product_index.
         j_sph(:, iaux) = RESHAPE(sph_matrix, [nsgf_a*nsgf_b])
      END DO

   END SUBROUTINE hfx_lri_transform_ao_pair

! **************************************************************************************************
!> \brief Builds one pair-local fitting coefficient block using reusable transforms and a
!>        thread-private integral context.
!> \param basis_set_AO ...
!> \param basis_set_RI ...
!> \param pair_block ...
!> \param ri_data ...
!> \param integral_context ...
!> \param ao_transform_cache ...
!> \param ri_transform_cache ...
!> \param metric_cache ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_build_pair_coefficients(basis_set_AO, basis_set_RI, pair_block, ri_data, &
                                              integral_context, ao_transform_cache, &
                                              ri_transform_cache, metric_cache)

      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_AO, basis_set_RI
      TYPE(hfx_lri_pair_block_type), INTENT(INOUT)       :: pair_block
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: integral_context
      TYPE(lri_basis_transform_type), DIMENSION(:), &
         INTENT(IN)                                      :: ao_transform_cache, ri_transform_cache
      TYPE(hfx_lri_metric_cache_type), INTENT(INOUT)     :: metric_cache

      INTEGER                                            :: cache_slot, eig_small, ikind, jkind, &
                                                            metric_dim, ncart_ai, ncart_aj, &
                                                            ncart_ri, ncart_rj, offset_j
      LOGICAL                                            :: cache_transpose, metric_zero
      REAL(dp)                                           :: c_norm, eig_cond, eig_max, eig_min
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: c_sph, j_aux_cart, j_cart, j_sph, &
                                                            metric_sph, sphi_aux
      REAL(dp), DIMENSION(3)                             :: r0
      TYPE(gto_basis_set_type), POINTER                  :: obas_i, obas_j, rbas_i, rbas_j

      NULLIFY (obas_i, obas_j, rbas_i, rbas_j)
      ikind = pair_block%ikind
      jkind = pair_block%jkind
      IF (ikind < 1 .OR. ikind > SIZE(basis_set_AO) .OR. &
          ikind > SIZE(basis_set_RI) .OR. jkind < 1 .OR. &
          jkind > SIZE(basis_set_AO) .OR. jkind > SIZE(basis_set_RI)) THEN
         CPABORT("Invalid basis kind in LRI-HFXk pair coefficient build")
      END IF

      obas_i => basis_set_AO(ikind)%gto_basis_set
      obas_j => basis_set_AO(jkind)%gto_basis_set
      rbas_i => basis_set_RI(ikind)%gto_basis_set
      rbas_j => basis_set_RI(jkind)%gto_basis_set
      IF (.NOT. ASSOCIATED(obas_i) .OR. .NOT. ASSOCIATED(obas_j) .OR. &
          .NOT. ASSOCIATED(rbas_i) .OR. .NOT. ASSOCIATED(rbas_j)) THEN
         CPABORT("Missing basis in LRI-HFXk pair coefficient build")
      END IF
      IF (.NOT. ALLOCATED(ao_transform_cache(ikind)%sphi) .OR. &
          .NOT. ALLOCATED(ao_transform_cache(jkind)%sphi) .OR. &
          .NOT. ALLOCATED(ri_transform_cache(ikind)%sphi) .OR. &
          .NOT. ALLOCATED(ri_transform_cache(jkind)%sphi)) THEN
         CPABORT("Missing basis transformation in LRI-HFXk pair coefficient build")
      END IF

      ncart_ai = SIZE(ao_transform_cache(ikind)%sphi, 1)
      ncart_aj = SIZE(ao_transform_cache(jkind)%sphi, 1)
      ncart_ri = SIZE(ri_transform_cache(ikind)%sphi, 1)
      ncart_rj = SIZE(ri_transform_cache(jkind)%sphi, 1)
      metric_dim = ncart_ri
      IF (.NOT. pair_block%onsite_pair) metric_dim = metric_dim + ncart_rj
      IF (pair_block%naux_pair /= pair_block%naux_i + &
          MERGE(0, pair_block%naux_j, pair_block%onsite_pair)) THEN
         CPABORT("Inconsistent auxiliary dimension in LRI-HFXk pair coefficient build")
      END IF

      r0(:) = 0.0_dp
      ALLOCATE (j_cart(ncart_ai*ncart_aj, metric_dim), SOURCE=0.0_dp)
      CALL hfx_lri_fill_3c_block(obas_i, obas_j, rbas_i, r0, pair_block%rab, r0, &
                                 integral_context, ri_data, j_cart, 0, ncart_aj)
      IF (.NOT. pair_block%onsite_pair) THEN
         offset_j = ncart_ri
         CALL hfx_lri_fill_3c_block(obas_i, obas_j, rbas_j, r0, pair_block%rab, &
                                    pair_block%rab, integral_context, ri_data, j_cart, offset_j, &
                                    ncart_aj)
      END IF

      ALLOCATE (sphi_aux(metric_dim, pair_block%naux_pair), SOURCE=0.0_dp)
      sphi_aux(1:ncart_ri, 1:pair_block%naux_i) = ri_transform_cache(ikind)%sphi
      IF (.NOT. pair_block%onsite_pair) THEN
         sphi_aux(ncart_ri + 1:metric_dim, pair_block%naux_i + 1:pair_block%naux_pair) = &
            ri_transform_cache(jkind)%sphi
      END IF

      ! Assemble S^T M S from cached center-center blocks. This is algebraically identical to
      ! transforming the full Cartesian pair metric, but avoids recomputing the on-site metric for
      ! every AO pair and removes the dense Cartesian pair-metric workspace.
      ALLOCATE (metric_sph(pair_block%naux_pair, pair_block%naux_pair), SOURCE=0.0_dp)
      CALL hfx_lri_metric_cache_get( &
         basis_set_RI, ri_transform_cache, ri_data, integral_context, ikind, ikind, r0, r0, &
         metric_cache, cache_slot, cache_transpose, metric_zero)
      IF (.NOT. metric_zero) THEN
         IF (ANY(SHAPE(metric_cache%entries(cache_slot)%values) /= &
                 [pair_block%naux_i, pair_block%naux_i])) THEN
            CPABORT("Inconsistent on-site metric dimensions in LRI-HFXk coefficient build")
         END IF
         IF (cache_transpose) THEN
            metric_sph(1:pair_block%naux_i, 1:pair_block%naux_i) = &
               TRANSPOSE(metric_cache%entries(cache_slot)%values)
         ELSE
            metric_sph(1:pair_block%naux_i, 1:pair_block%naux_i) = &
               metric_cache%entries(cache_slot)%values
         END IF
      END IF
      IF (.NOT. pair_block%onsite_pair) THEN
         CALL hfx_lri_metric_cache_get( &
            basis_set_RI, ri_transform_cache, ri_data, integral_context, ikind, jkind, &
            r0, pair_block%rab, metric_cache, cache_slot, cache_transpose, metric_zero)
         IF (.NOT. metric_zero) THEN
            IF (cache_transpose) THEN
               IF (ANY(SHAPE(metric_cache%entries(cache_slot)%values) /= &
                       [pair_block%naux_j, pair_block%naux_i])) THEN
                  CPABORT("Inconsistent cross metric dimensions in LRI-HFXk coefficient build")
               END IF
               metric_sph(1:pair_block%naux_i, pair_block%naux_i + 1:pair_block%naux_pair) = &
                  TRANSPOSE(metric_cache%entries(cache_slot)%values)
            ELSE
               IF (ANY(SHAPE(metric_cache%entries(cache_slot)%values) /= &
                       [pair_block%naux_i, pair_block%naux_j])) THEN
                  CPABORT("Inconsistent cross metric dimensions in LRI-HFXk coefficient build")
               END IF
               metric_sph(1:pair_block%naux_i, pair_block%naux_i + 1:pair_block%naux_pair) = &
                  metric_cache%entries(cache_slot)%values
            END IF
            metric_sph(pair_block%naux_i + 1:pair_block%naux_pair, 1:pair_block%naux_i) = &
               TRANSPOSE(metric_sph(1:pair_block%naux_i, &
                                    pair_block%naux_i + 1:pair_block%naux_pair))
         END IF
         CALL hfx_lri_metric_cache_get( &
            basis_set_RI, ri_transform_cache, ri_data, integral_context, jkind, jkind, &
            pair_block%rab, pair_block%rab, metric_cache, cache_slot, cache_transpose, metric_zero)
         IF (.NOT. metric_zero) THEN
            IF (ANY(SHAPE(metric_cache%entries(cache_slot)%values) /= &
                    [pair_block%naux_j, pair_block%naux_j])) THEN
               CPABORT("Inconsistent second on-site metric dimensions in LRI-HFXk coefficient build")
            END IF
            IF (cache_transpose) THEN
               metric_sph(pair_block%naux_i + 1:pair_block%naux_pair, &
                          pair_block%naux_i + 1:pair_block%naux_pair) = &
                  TRANSPOSE(metric_cache%entries(cache_slot)%values)
            ELSE
               metric_sph(pair_block%naux_i + 1:pair_block%naux_pair, &
                          pair_block%naux_i + 1:pair_block%naux_pair) = &
                  metric_cache%entries(cache_slot)%values
            END IF
         END IF
      END IF

      ALLOCATE (j_aux_cart(ncart_ai*ncart_aj, pair_block%naux_pair))
      j_aux_cart(:, :) = MATMUL(j_cart, sphi_aux)
      DEALLOCATE (j_cart, sphi_aux)
      CALL hfx_lri_transform_ao_pair(j_aux_cart, ao_transform_cache(ikind)%sphi, &
                                     ao_transform_cache(jkind)%sphi, j_sph)
      DEALLOCATE (j_aux_cart)

      CALL hfx_lri_probe_metric_inverse(metric_sph, j_sph, ri_data%eps_eigval, &
                                        eig_min, eig_max, eig_cond, &
                                        eig_small, c_norm, c_out=c_sph)
      IF (.NOT. ALLOCATED(c_sph)) THEN
         CPABORT("Failed to build LRI-HFXk pair fitting coefficients")
      END IF
      CALL hfx_lri_pack_pair_coefficients(pair_block, c_sph, pair_block%nao_i, &
                                          pair_block%nao_j, pair_block%naux_pair)
      IF (c_norm < ri_data%filter_eps) DEALLOCATE (pair_block%c_fit)

      ! metric_sph is deliberately not retained. It is only a setup workspace, and retaining one
      ! dense pair metric per AO pair costs O(Npair*Naux_pair**2) memory without helping the SCF path.

   END SUBROUTINE hfx_lri_build_pair_coefficients

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
!> \brief Diagonalizes a pair metric and forms C = J M^-1.
!> \param metric_block ...
!> \param j_block ...
!> \param eig_threshold ...
!> \param eig_min ...
!> \param eig_max ...
!> \param eig_cond ...
!> \param eig_small ...
!> \param c_norm ...
!> \param c_out ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_probe_metric_inverse(metric_block, j_block, eig_threshold, &
                                           eig_min, eig_max, eig_cond, &
                                           eig_small, c_norm, c_out)

      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: metric_block, j_block
      REAL(dp), INTENT(IN)                               :: eig_threshold
      REAL(dp), INTENT(OUT)                              :: eig_min, eig_max, eig_cond
      INTEGER, INTENT(OUT)                               :: eig_small
      REAL(dp), INTENT(OUT)                              :: c_norm
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT), OPTIONAL                           :: c_out

      INTEGER                                            :: i, info, lwork, ncol, nrow
      REAL(dp)                                           :: eig_abs
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: eigenvalues, work
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: c_block, eigvec, tmp

      eig_min = 0.0_dp
      eig_max = 0.0_dp
      eig_cond = 0.0_dp
      eig_small = 0
      c_norm = 0.0_dp
      CPASSERT(eig_threshold >= 0.0_dp)

      ncol = SIZE(metric_block, 1)
      nrow = SIZE(j_block, 1)
      IF (ncol == 0 .OR. nrow == 0) RETURN

      ALLOCATE (eigvec(ncol, ncol), eigenvalues(ncol))
      eigvec(:, :) = 0.5_dp*(metric_block + TRANSPOSE(metric_block))

      lwork = -1
      ALLOCATE (work(1))
      CALL dsyev('V', 'L', ncol, eigvec, ncol, eigenvalues, work, lwork, info)
      IF (info /= 0) THEN
         CPABORT("DSYEV workspace query failed in LRI-HFXk pair metric")
      END IF
      lwork = MAX(1, INT(work(1)))
      DEALLOCATE (work)

      ALLOCATE (work(lwork))
      CALL dsyev('V', 'L', ncol, eigvec, ncol, eigenvalues, work, lwork, info)
      DEALLOCATE (work)
      IF (info /= 0) THEN
         CPABORT("DSYEV failed in LRI-HFXk pair metric")
      END IF

      eig_max = MAXVAL(ABS(eigenvalues))
      eig_min = HUGE(1.0_dp)
      DO i = 1, ncol
         eig_abs = ABS(eigenvalues(i))
         IF (eigenvalues(i) < eig_threshold) THEN
            eig_small = eig_small + 1
         ELSE
            eig_min = MIN(eig_min, eig_abs)
         END IF
      END DO
      IF (eig_min < HUGE(1.0_dp) .AND. eig_max > 0.0_dp) eig_cond = eig_max/eig_min

      ALLOCATE (tmp(nrow, ncol), c_block(nrow, ncol))
      tmp(:, :) = MATMUL(j_block, eigvec)
      DO i = 1, ncol
         IF (eigenvalues(i) >= eig_threshold) THEN
            tmp(:, i) = tmp(:, i)/eigenvalues(i)
         ELSE
            tmp(:, i) = 0.0_dp
         END IF
      END DO
      c_block(:, :) = MATMUL(tmp, TRANSPOSE(eigvec))
      c_norm = SQRT(SUM(c_block*c_block))
      IF (PRESENT(c_out)) THEN
         IF (ALLOCATED(c_out)) DEALLOCATE (c_out)
         CALL MOVE_ALLOC(c_block, c_out)
      END IF

      IF (ALLOCATED(c_block)) DEALLOCATE (c_block)
      DEALLOCATE (tmp, eigenvalues, eigvec)

   END SUBROUTINE hfx_lri_probe_metric_inverse

! **************************************************************************************************
!> \brief Initializes a pair-local LRI coefficient block.
!> \param pair_block ...
!> \param iatom ...
!> \param jatom ...
!> \param ikind ...
!> \param jkind ...
!> \param cell ...
!> \param rab ...
!> \param nao_i ...
!> \param nao_j ...
!> \param naux_i ...
!> \param naux_j ...
!> \param naux_pair ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_init_pair_block(pair_block, iatom, jatom, ikind, jkind, cell, rab, nao_i, &
                                      nao_j, naux_i, naux_j, naux_pair)

      TYPE(hfx_lri_pair_block_type), INTENT(INOUT)       :: pair_block
      INTEGER, INTENT(IN)                                :: iatom, jatom, ikind, jkind
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rab
      INTEGER, INTENT(IN)                                :: nao_i, nao_j, naux_i, naux_j, naux_pair

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
      pair_block%naux_i = naux_i
      pair_block%naux_j = naux_j
      pair_block%naux_pair = naux_pair

      pair_block%onsite_pair = (iatom == jatom .AND. SUM(rab*rab) < 1.0E-20_dp)
      pair_block%naux_centers = 1
      pair_block%aux_atom(:) = 0
      pair_block%aux_kind(:) = 0
      pair_block%aux_cell(:, :) = 0
      pair_block%aux_r(:, :) = 0.0_dp
      pair_block%aux_offset(:) = 0
      pair_block%aux_size(:) = 0

      pair_block%aux_atom(1) = iatom
      pair_block%aux_kind(1) = ikind
      pair_block%aux_cell(:, 1) = 0
      pair_block%aux_r(:, 1) = 0.0_dp
      pair_block%aux_offset(1) = 0
      pair_block%aux_size(1) = naux_i

      IF (.NOT. pair_block%onsite_pair) THEN
         pair_block%naux_centers = 2
         pair_block%aux_atom(2) = jatom
         pair_block%aux_kind(2) = jkind
         pair_block%aux_cell(:, 2) = cell(:)
         pair_block%aux_r(:, 2) = rab(:)
         pair_block%aux_offset(2) = naux_i
         pair_block%aux_size(2) = naux_j
      END IF

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

      DO imu = 1, nao_i
         DO inu = 1, nao_j
            irow = hfx_lri_product_index(imu, inu, nao_j)
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
!> \brief Transfers an allocated K block together with its target-pair metadata.
!> \param k_pair ...
!> \param pair_block ...
!> \param pair_index ...
!> \param k_block ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_store_k_block(k_pair, pair_block, pair_index, k_block)

      USE hfx_lri_types, ONLY: hfx_lri_k_block_type
      TYPE(hfx_lri_k_block_type), INTENT(INOUT)          :: k_pair
      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair_block
      INTEGER, INTENT(IN)                                :: pair_index
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: k_block

      IF (ALLOCATED(k_pair%k_block)) DEALLOCATE (k_pair%k_block)
      IF (.NOT. ALLOCATED(k_block)) CPABORT("Cannot store an unallocated LRI-HFXk K block")

      k_pair%pair_index = pair_index
      k_pair%iatom = pair_block%iatom
      k_pair%jatom = pair_block%jatom
      k_pair%ikind = pair_block%ikind
      k_pair%jkind = pair_block%jkind
      k_pair%cell(:) = pair_block%cell(:)
      k_pair%nao_i = SIZE(k_block, 1)
      k_pair%nao_j = SIZE(k_block, 2)

      CALL MOVE_ALLOC(k_block, k_pair%k_block)

   END SUBROUTINE hfx_lri_store_k_block

! **************************************************************************************************
!> \brief Hashes an atom-pair/cell key into a lookup-table slot.
!> \param iatom ...
!> \param jatom ...
!> \param cell ...
!> \param capacity ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION hfx_lri_pair_hash(iatom, jatom, cell, capacity) RESULT(islot)

      INTEGER, INTENT(IN)                                :: iatom, jatom
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell
      INTEGER, INTENT(IN)                                :: capacity

      INTEGER                                            :: idir
      INTEGER(int_8)                                     :: hash_value, table_size

      CPASSERT(capacity > 0)
      table_size = INT(capacity, int_8)
      hash_value = MODULO(INT(iatom, int_8), table_size)
      hash_value = MODULO(104729_int_8*hash_value + INT(jatom, int_8), table_size)
      DO idir = 1, 3
         hash_value = MODULO(104729_int_8*hash_value + INT(cell(idir), int_8), table_size)
      END DO
      islot = INT(hash_value) + 1

   END FUNCTION hfx_lri_pair_hash

! **************************************************************************************************
!> \brief Hashes a translationally invariant two-center metric key.
!> \param kind_left ...
!> \param kind_right ...
!> \param displacement ...
!> \param capacity ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION hfx_lri_metric_hash(kind_left, kind_right, displacement, capacity) RESULT(islot)

      INTEGER, INTENT(IN)                                :: kind_left, kind_right
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: displacement
      INTEGER, INTENT(IN)                                :: capacity

      INTEGER                                            :: idir
      INTEGER(int_8)                                     :: hash_value, position_key, table_size

      CPASSERT(capacity > 0)
      table_size = INT(capacity, int_8)
      hash_value = MODULO(INT(kind_left, int_8), table_size)
      hash_value = MODULO(104729_int_8*hash_value + INT(kind_right, int_8), table_size)
      DO idir = 1, 3
         position_key = NINT(displacement(idir)/lri_metric_cache_position_tol, KIND=int_8)
         hash_value = MODULO(104729_int_8*hash_value + MODULO(position_key, table_size), &
                             table_size)
      END DO
      islot = INT(hash_value) + 1

   END FUNCTION hfx_lri_metric_hash

! **************************************************************************************************
!> \brief Chooses one representative of V_AB(R)=transpose(V_BA(-R)).
!> \param kind_left_in ...
!> \param kind_right_in ...
!> \param displacement_in ...
!> \param kind_left ...
!> \param kind_right ...
!> \param displacement ...
!> \param reversed ...
! **************************************************************************************************
   PURE SUBROUTINE hfx_lri_metric_key_canonical(kind_left_in, kind_right_in, displacement_in, &
                                                kind_left, kind_right, displacement, reversed)

      INTEGER, INTENT(IN)                                :: kind_left_in, kind_right_in
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: displacement_in
      INTEGER, INTENT(OUT)                               :: kind_left, kind_right
      REAL(dp), DIMENSION(3), INTENT(OUT)                :: displacement
      LOGICAL, INTENT(OUT), OPTIONAL                     :: reversed

      INTEGER                                            :: idir
      INTEGER(int_8)                                     :: position_key
      LOGICAL                                            :: reverse_key

      reverse_key = kind_left_in > kind_right_in
      IF (kind_left_in == kind_right_in) THEN
         DO idir = 1, 3
            position_key = NINT(displacement_in(idir)/lri_metric_cache_position_tol, KIND=int_8)
            IF (position_key < 0_int_8) THEN
               reverse_key = .TRUE.
               EXIT
            ELSE IF (position_key > 0_int_8) THEN
               EXIT
            END IF
         END DO
      END IF
      IF (reverse_key) THEN
         kind_left = kind_right_in
         kind_right = kind_left_in
         displacement(:) = -displacement_in(:)
      ELSE
         kind_left = kind_left_in
         kind_right = kind_right_in
         displacement(:) = displacement_in(:)
      END IF
      IF (PRESENT(reversed)) reversed = reverse_key

   END SUBROUTINE hfx_lri_metric_key_canonical

! **************************************************************************************************
!> \brief Rebuilds the open-addressed index of the rank-shared metric table.
!> \param metric_table ...
!> \param requested_capacity ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_metric_table_rehash(metric_table, requested_capacity)

      TYPE(hfx_lri_metric_table_type), INTENT(INOUT)     :: metric_table
      INTEGER, INTENT(IN)                                :: requested_capacity

      INTEGER                                            :: ientry, islot, nprobe
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: new_slot

      metric_table%capacity = MAX(17, requested_capacity)
      ALLOCATE (new_slot(metric_table%capacity), SOURCE=0)
      DO ientry = 1, metric_table%nentries
         islot = hfx_lri_metric_hash( &
                 metric_table%entries(ientry)%kind_left, &
                 metric_table%entries(ientry)%kind_right, &
                 metric_table%entries(ientry)%displacement, metric_table%capacity)
         DO nprobe = 1, metric_table%capacity
            IF (new_slot(islot) == 0) THEN
               new_slot(islot) = ientry
               EXIT
            END IF
            islot = MOD(islot, metric_table%capacity) + 1
         END DO
         IF (nprobe > metric_table%capacity) CPABORT("LRI-HFXk metric table rehash failed")
      END DO
      CALL MOVE_ALLOC(new_slot, metric_table%slot)

   END SUBROUTINE hfx_lri_metric_table_rehash

! **************************************************************************************************
!> \brief Adds one geometry-only metric key. Values are filled in a later parallel pass.
!> \param metric_table ...
!> \param kind_left_in ...
!> \param kind_right_in ...
!> \param displacement_in ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_metric_table_add_key(metric_table, kind_left_in, kind_right_in, &
                                           displacement_in)

      TYPE(hfx_lri_metric_table_type), INTENT(INOUT)     :: metric_table
      INTEGER, INTENT(IN)                                :: kind_left_in, kind_right_in
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: displacement_in

      INTEGER                                            :: ientry, islot, kind_left, kind_right, &
                                                            nprobe, old_size
      REAL(dp), DIMENSION(3)                             :: displacement
      TYPE(hfx_lri_metric_cache_entry_type), &
         ALLOCATABLE, DIMENSION(:)                       :: new_entries

      CALL hfx_lri_metric_key_canonical(kind_left_in, kind_right_in, displacement_in, &
                                        kind_left, kind_right, displacement)
      IF (.NOT. ALLOCATED(metric_table%entries)) THEN
         ALLOCATE (metric_table%entries(256))
         CALL hfx_lri_metric_table_rehash(metric_table, 1031)
      END IF

      islot = hfx_lri_metric_hash(kind_left, kind_right, displacement, metric_table%capacity)
      DO nprobe = 1, metric_table%capacity
         ientry = metric_table%slot(islot)
         IF (ientry == 0) EXIT
         IF (metric_table%entries(ientry)%kind_left == kind_left .AND. &
             metric_table%entries(ientry)%kind_right == kind_right .AND. &
             MAXVAL(ABS(metric_table%entries(ientry)%displacement - displacement)) <= &
             lri_metric_cache_position_tol) RETURN
         islot = MOD(islot, metric_table%capacity) + 1
      END DO
      IF (nprobe > metric_table%capacity) CPABORT("LRI-HFXk metric key table is full")

      IF (10*(metric_table%nentries + 1) > 6*metric_table%capacity) THEN
         CALL hfx_lri_metric_table_rehash(metric_table, 2*metric_table%capacity + 1)
         islot = hfx_lri_metric_hash(kind_left, kind_right, displacement, metric_table%capacity)
         DO WHILE (metric_table%slot(islot) /= 0)
            islot = MOD(islot, metric_table%capacity) + 1
         END DO
      END IF
      IF (metric_table%nentries == SIZE(metric_table%entries)) THEN
         old_size = SIZE(metric_table%entries)
         ALLOCATE (new_entries(2*old_size))
         new_entries(1:old_size) = metric_table%entries(1:old_size)
         CALL MOVE_ALLOC(new_entries, metric_table%entries)
      END IF

      metric_table%nentries = metric_table%nentries + 1
      ientry = metric_table%nentries
      metric_table%entries(ientry)%valid = .TRUE.
      metric_table%entries(ientry)%kind_left = kind_left
      metric_table%entries(ientry)%kind_right = kind_right
      metric_table%entries(ientry)%displacement(:) = displacement(:)
      metric_table%slot(islot) = ientry

   END SUBROUTINE hfx_lri_metric_table_add_key

! **************************************************************************************************
!> \brief Finds one immutable shared metric block, accepting its transpose-equivalent key.
!> \param metric_table ...
!> \param kind_left_in ...
!> \param kind_right_in ...
!> \param displacement_in ...
!> \param ientry ...
!> \param use_transpose ...
!> \param found ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_metric_table_find(metric_table, kind_left_in, kind_right_in, &
                                        displacement_in, ientry, use_transpose, found)

      TYPE(hfx_lri_metric_table_type), INTENT(IN)        :: metric_table
      INTEGER, INTENT(IN)                                :: kind_left_in, kind_right_in
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: displacement_in
      INTEGER, INTENT(OUT)                               :: ientry
      LOGICAL, INTENT(OUT)                               :: use_transpose, found

      INTEGER                                            :: candidate, islot, kind_left, kind_right, &
                                                            nprobe
      REAL(dp), DIMENSION(3)                             :: displacement

      ientry = 0
      found = .FALSE.
      CALL hfx_lri_metric_key_canonical(kind_left_in, kind_right_in, displacement_in, &
                                        kind_left, kind_right, displacement, use_transpose)
      IF (.NOT. metric_table%initialized .OR. metric_table%capacity <= 0 .OR. &
          .NOT. ALLOCATED(metric_table%slot)) RETURN
      islot = hfx_lri_metric_hash(kind_left, kind_right, displacement, metric_table%capacity)
      DO nprobe = 1, metric_table%capacity
         candidate = metric_table%slot(islot)
         IF (candidate == 0) RETURN
         IF (metric_table%entries(candidate)%kind_left == kind_left .AND. &
             metric_table%entries(candidate)%kind_right == kind_right .AND. &
             MAXVAL(ABS(metric_table%entries(candidate)%displacement - displacement)) <= &
             lri_metric_cache_position_tol) THEN
            ientry = candidate
            found = .TRUE.
            RETURN
         END IF
         islot = MOD(islot, metric_table%capacity) + 1
      END DO

   END SUBROUTINE hfx_lri_metric_table_find

! **************************************************************************************************
!> \brief Returns the logical key represented by a signed pair reference. A negative reference is
!>        the zero-copy transposed/translated view (j,i,-R) of the stored pair (i,j,R).
!> \param pair_blocks ...
!> \param pair_ref ...
!> \param iatom ...
!> \param jatom ...
!> \param cell ...
! **************************************************************************************************
   PURE SUBROUTINE hfx_lri_pair_ref_key(pair_blocks, pair_ref, iatom, jatom, cell)

      TYPE(hfx_lri_pair_block_type), DIMENSION(:), &
         INTENT(IN)                                      :: pair_blocks
      INTEGER, INTENT(IN)                                :: pair_ref
      INTEGER, INTENT(OUT)                               :: iatom, jatom
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: cell

      INTEGER                                            :: ipair

      ipair = ABS(pair_ref)
      IF (pair_ref > 0) THEN
         iatom = pair_blocks(ipair)%iatom
         jatom = pair_blocks(ipair)%jatom
         cell(:) = pair_blocks(ipair)%cell(:)
      ELSE
         iatom = pair_blocks(ipair)%jatom
         jatom = pair_blocks(ipair)%iatom
         cell(:) = -pair_blocks(ipair)%cell(:)
      END IF

   END SUBROUTINE hfx_lri_pair_ref_key

! **************************************************************************************************
!> \brief Builds O(1) pair lookup and contiguous per-first-atom pair lists.
!> \param pair_lookup ...
!> \param pair_blocks ...
!> \param npair_blocks ...
!> \param natom ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_pair_lookup_create(pair_lookup, pair_blocks, npair_blocks, natom)

      TYPE(lri_pair_lookup_type), INTENT(OUT)            :: pair_lookup
      TYPE(hfx_lri_pair_block_type), DIMENSION(:), &
         INTENT(IN)                                      :: pair_blocks
      INTEGER, INTENT(IN)                                :: npair_blocks, natom

      INTEGER                                            :: existing_iatom, existing_jatom, iatom, &
                                                            ipair, islot, nprobe, nrefs, pair_ref, &
                                                            ref_iatom, ref_jatom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cursor
      INTEGER, DIMENSION(3)                              :: existing_cell, ref_cell
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: include_reverse

      CPASSERT(npair_blocks >= 0)
      CPASSERT(npair_blocks <= SIZE(pair_blocks))
      CPASSERT(natom > 0)

      ! The stored neighbor list is symmetric and therefore contains only one representative of
      ! each pair. Add signed reverse aliases to the lookup without copying any C tensor.
      pair_lookup%capacity = MAX(17, 4*npair_blocks + 1)
      ALLOCATE (pair_lookup%slot(pair_lookup%capacity), SOURCE=0)
      ALLOCATE (pair_lookup%atom_count(natom), SOURCE=0)
      ALLOCATE (pair_lookup%atom_offset(natom + 1), SOURCE=1)
      ALLOCATE (include_reverse(npair_blocks), SOURCE=.FALSE.)

      ! Insert the stored orientations first, so a directly stored key wins over an alias.
      DO ipair = 1, npair_blocks
         IF (.NOT. ALLOCATED(pair_blocks(ipair)%c_fit)) CYCLE
         pair_ref = ipair
         CALL hfx_lri_pair_ref_key(pair_blocks, pair_ref, ref_iatom, ref_jatom, ref_cell)
         iatom = ref_iatom
         IF (iatom < 1 .OR. iatom > natom) THEN
            CPABORT("Invalid first-atom index in an LRI-HFXk pair block")
         END IF
         pair_lookup%atom_count(iatom) = pair_lookup%atom_count(iatom) + 1
         islot = hfx_lri_pair_hash(ref_iatom, ref_jatom, ref_cell, pair_lookup%capacity)
         DO nprobe = 1, pair_lookup%capacity
            IF (pair_lookup%slot(islot) == 0) THEN
               pair_lookup%slot(islot) = pair_ref
               EXIT
            END IF
            CALL hfx_lri_pair_ref_key(pair_blocks, pair_lookup%slot(islot), &
                                      existing_iatom, existing_jatom, existing_cell)
            IF (existing_iatom == ref_iatom .AND. existing_jatom == ref_jatom .AND. &
                ALL(existing_cell == ref_cell)) THEN
               CPABORT("Duplicate atom-pair/cell key in the LRI-HFXk pair list")
            END IF
            islot = MOD(islot, pair_lookup%capacity) + 1
         END DO
         IF (nprobe > pair_lookup%capacity) THEN
            CPABORT("LRI-HFXk pair lookup table is full")
         END IF
      END DO

      ! Add missing reverse orientations as negative references. Self pairs in the origin cell and
      ! any unexpectedly explicit reverse pair already have a direct entry and are not duplicated.
      DO ipair = 1, npair_blocks
         IF (.NOT. ALLOCATED(pair_blocks(ipair)%c_fit)) CYCLE
         pair_ref = -ipair
         CALL hfx_lri_pair_ref_key(pair_blocks, pair_ref, ref_iatom, ref_jatom, ref_cell)
         islot = hfx_lri_pair_hash(ref_iatom, ref_jatom, ref_cell, pair_lookup%capacity)
         DO nprobe = 1, pair_lookup%capacity
            IF (pair_lookup%slot(islot) == 0) THEN
               pair_lookup%slot(islot) = pair_ref
               include_reverse(ipair) = .TRUE.
               pair_lookup%atom_count(ref_iatom) = pair_lookup%atom_count(ref_iatom) + 1
               EXIT
            END IF
            CALL hfx_lri_pair_ref_key(pair_blocks, pair_lookup%slot(islot), &
                                      existing_iatom, existing_jatom, existing_cell)
            IF (existing_iatom == ref_iatom .AND. existing_jatom == ref_jatom .AND. &
                ALL(existing_cell == ref_cell)) EXIT
            islot = MOD(islot, pair_lookup%capacity) + 1
         END DO
         IF (nprobe > pair_lookup%capacity) THEN
            CPABORT("LRI-HFXk pair lookup table is full while adding reverse views")
         END IF
      END DO

      pair_lookup%atom_offset(1) = 1
      DO iatom = 1, natom
         pair_lookup%atom_offset(iatom + 1) = pair_lookup%atom_offset(iatom) + &
                                              pair_lookup%atom_count(iatom)
      END DO
      nrefs = pair_lookup%atom_offset(natom + 1) - 1
      ALLOCATE (pair_lookup%atom_pairs(nrefs), SOURCE=0)
      ALLOCATE (cursor(natom))
      cursor(:) = pair_lookup%atom_offset(1:natom)
      DO ipair = 1, npair_blocks
         IF (.NOT. ALLOCATED(pair_blocks(ipair)%c_fit)) CYCLE
         pair_ref = ipair
         CALL hfx_lri_pair_ref_key(pair_blocks, pair_ref, ref_iatom, ref_jatom, ref_cell)
         pair_lookup%atom_pairs(cursor(ref_iatom)) = pair_ref
         cursor(ref_iatom) = cursor(ref_iatom) + 1
         IF (include_reverse(ipair)) THEN
            pair_ref = -ipair
            CALL hfx_lri_pair_ref_key(pair_blocks, pair_ref, ref_iatom, ref_jatom, ref_cell)
            pair_lookup%atom_pairs(cursor(ref_iatom)) = pair_ref
            cursor(ref_iatom) = cursor(ref_iatom) + 1
         END IF
      END DO
      DEALLOCATE (cursor)
      DEALLOCATE (include_reverse)

   END SUBROUTINE hfx_lri_pair_lookup_create

! **************************************************************************************************
!> \brief Finds a pair block from its first atom, second atom, and lattice cell.
!> \param pair_lookup ...
!> \param pair_blocks ...
!> \param iatom ...
!> \param jatom ...
!> \param cell ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION hfx_lri_pair_lookup_find(pair_lookup, pair_blocks, iatom, jatom, cell) &
      RESULT(ipair_found)

      TYPE(lri_pair_lookup_type), INTENT(IN)             :: pair_lookup
      TYPE(hfx_lri_pair_block_type), DIMENSION(:), &
         INTENT(IN)                                      :: pair_blocks
      INTEGER, INTENT(IN)                                :: iatom, jatom
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell

      INTEGER                                            :: existing_iatom, existing_jatom, ipair, &
                                                            islot, nprobe
      INTEGER, DIMENSION(3)                              :: existing_cell

      ipair_found = 0
      IF (pair_lookup%capacity <= 0) RETURN
      islot = hfx_lri_pair_hash(iatom, jatom, cell, pair_lookup%capacity)
      DO nprobe = 1, pair_lookup%capacity
         ipair = pair_lookup%slot(islot)
         IF (ipair == 0) RETURN
         CALL hfx_lri_pair_ref_key(pair_blocks, ipair, existing_iatom, existing_jatom, existing_cell)
         IF (existing_iatom == iatom .AND. existing_jatom == jatom .AND. &
             ALL(existing_cell == cell)) THEN
            ipair_found = ipair
            RETURN
         END IF
         islot = MOD(islot, pair_lookup%capacity) + 1
      END DO

   END FUNCTION hfx_lri_pair_lookup_find

! **************************************************************************************************
!> \brief Adds one logical real-space density block to a sparse cache.
!> \param density_cache ...
!> \param iatom ...
!> \param jatom ...
!> \param ikind ...
!> \param jkind ...
!> \param cell ...
!> \param rab ...
!> \param source ...
!> \param source_transposed ...
!> \param factor ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_density_cache_add(density_cache, iatom, jatom, ikind, jkind, cell, rab, &
                                        source, source_transposed, factor)

      TYPE(lri_density_cache_type), INTENT(INOUT)        :: density_cache
      INTEGER, INTENT(IN)                                :: iatom, jatom, ikind, jkind
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rab
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: source
      LOGICAL, INTENT(IN)                                :: source_transposed
      REAL(dp), INTENT(IN)                               :: factor

      INTEGER                                            :: iblock, islot, nprobe

      CPASSERT(ALLOCATED(density_cache%blocks))
      CPASSERT(ALLOCATED(density_cache%slot))
      CPASSERT(density_cache%capacity == SIZE(density_cache%slot))
      islot = hfx_lri_pair_hash(iatom, jatom, cell, density_cache%capacity)
      DO nprobe = 1, density_cache%capacity
         iblock = density_cache%slot(islot)
         IF (iblock == 0) THEN
            density_cache%nblocks = density_cache%nblocks + 1
            iblock = density_cache%nblocks
            IF (iblock > SIZE(density_cache%blocks)) THEN
               CPABORT("LRI-HFXk sparse density cache capacity was underestimated")
            END IF
            density_cache%slot(islot) = iblock
            density_cache%blocks(iblock)%iatom = iatom
            density_cache%blocks(iblock)%jatom = jatom
            density_cache%blocks(iblock)%ikind = ikind
            density_cache%blocks(iblock)%jkind = jkind
            density_cache%blocks(iblock)%cell(:) = cell(:)
            density_cache%blocks(iblock)%rab(:) = rab(:)
            IF (source_transposed) THEN
               ALLOCATE (density_cache%blocks(iblock)%values(SIZE(source, 2), SIZE(source, 1)), &
                         SOURCE=0.0_dp)
            ELSE
               ALLOCATE (density_cache%blocks(iblock)%values(SIZE(source, 1), SIZE(source, 2)), &
                         SOURCE=0.0_dp)
            END IF
            EXIT
         END IF
         IF (density_cache%blocks(iblock)%iatom == iatom .AND. &
             density_cache%blocks(iblock)%jatom == jatom .AND. &
             ALL(density_cache%blocks(iblock)%cell == cell)) EXIT
         islot = MOD(islot, density_cache%capacity) + 1
      END DO
      IF (nprobe > density_cache%capacity) CPABORT("LRI-HFXk sparse density hash table is full")
      IF (density_cache%blocks(iblock)%ikind /= ikind .OR. &
          density_cache%blocks(iblock)%jkind /= jkind .OR. &
          MAXVAL(ABS(density_cache%blocks(iblock)%rab - rab)) > &
          lri_metric_cache_position_tol) THEN
         CPABORT("Inconsistent geometry in LRI-HFXk density cache")
      END IF

      IF (source_transposed) THEN
         IF (ANY(SHAPE(density_cache%blocks(iblock)%values) /= &
                 [SIZE(source, 2), SIZE(source, 1)])) THEN
            CPABORT("Inconsistent transposed block shape in LRI-HFXk density cache")
         END IF
         density_cache%blocks(iblock)%values(:, :) = &
            density_cache%blocks(iblock)%values + factor*TRANSPOSE(source)
      ELSE
         IF (ANY(SHAPE(density_cache%blocks(iblock)%values) /= SHAPE(source))) THEN
            CPABORT("Inconsistent block shape in LRI-HFXk density cache")
         END IF
         density_cache%blocks(iblock)%values(:, :) = &
            density_cache%blocks(iblock)%values + factor*source
      END IF

   END SUBROUTINE hfx_lri_density_cache_add

! **************************************************************************************************
!> \brief Reconstructs and caches only the logical density blocks present in CP2K's symmetric
!>        real-space representation. This is the sparse equivalent of get_pmat_images().
!> \param rho_ao ...
!> \param cell_to_index ...
!> \param sab_nl ...
!> \param filter_eps ...
!> \param density_cache ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_density_cache_create(rho_ao, cell_to_index, sab_nl, filter_eps, density_cache)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      REAL(dp), INTENT(IN)                               :: filter_eps
      TYPE(lri_density_cache_type), INTENT(OUT)          :: density_cache

      INTEGER                                            :: iatom, iblock, icol, ikind, image, &
                                                            image_rev, irow, jatom, jkind, mepos, &
                                                            nentry
      INTEGER, DIMENSION(3)                              :: cell, cell_rev
      LOGICAL                                            :: found, source_transposed
      REAL(dp)                                           :: factor
      REAL(dp), DIMENSION(3)                             :: rab
      REAL(dp), DIMENSION(:, :), POINTER                 :: rho_block
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      NULLIFY (nl_iterator, rho_block)
      CPASSERT(ASSOCIATED(rho_ao))
      CPASSERT(ASSOCIATED(sab_nl))
      CPASSERT(filter_eps >= 0.0_dp)

      ! The symmetric list contains at most one stored representative per physical pair. Each
      ! representative contributes its direct orientation and, where present, its Hermitian image.
      CALL neighbor_list_iterator_create(nl_iterator, sab_nl)
      nentry = 0
      mepos = 0
      DO WHILE (neighbor_list_iterate(nl_iterator, mepos) == 0)
         nentry = nentry + 1
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      ALLOCATE (density_cache%blocks(2*nentry))
      density_cache%capacity = MAX(17, 4*nentry + 1)
      ALLOCATE (density_cache%slot(density_cache%capacity), SOURCE=0)
      density_cache%nblocks = 0

      CALL neighbor_list_iterator_create(nl_iterator, sab_nl)
      mepos = 0
      DO WHILE (neighbor_list_iterate(nl_iterator, mepos) == 0)
         CALL get_iterator_info(nl_iterator, mepos=mepos, iatom=iatom, jatom=jatom, &
                                ikind=ikind, jkind=jkind, cell=cell, r=rab)
         image = 1
         IF (ASSOCIATED(cell_to_index)) THEN
            IF (ANY(cell < LBOUND(cell_to_index)) .OR. ANY(cell > UBOUND(cell_to_index))) CYCLE
            image = cell_to_index(cell(1), cell(2), cell(3))
         END IF
         IF (image < 1 .OR. image > SIZE(rho_ao, 2)) CYCLE

         factor = 1.0_dp
         cell_rev(:) = -cell(:)
         image_rev = 0
         IF (ASSOCIATED(cell_to_index)) THEN
            IF (ALL(cell_rev >= LBOUND(cell_to_index)) .AND. &
                ALL(cell_rev <= UBOUND(cell_to_index))) THEN
               image_rev = cell_to_index(cell_rev(1), cell_rev(2), cell_rev(3))
            END IF
         ELSE
            image_rev = 1
         END IF
         IF (iatom == jatom .AND. image_rev >= 1 .AND. image_rev <= SIZE(rho_ao, 2)) THEN
            factor = 0.5_dp
         END IF

         irow = iatom
         icol = jatom
         source_transposed = .FALSE.
         IF (iatom > jatom) THEN
            irow = jatom
            icol = iatom
            source_transposed = .TRUE.
         END IF
         NULLIFY (rho_block)
         found = .FALSE.
         CALL dbcsr_get_block_p(matrix=rho_ao(1, image)%matrix, row=irow, col=icol, &
                                block=rho_block, found=found)
         IF (.NOT. found) CYCLE

         CALL hfx_lri_density_cache_add(density_cache, iatom, jatom, ikind, jkind, cell, rab, &
                                        rho_block, source_transposed, factor)
         IF (image_rev >= 1 .AND. image_rev <= SIZE(rho_ao, 2)) THEN
            CALL hfx_lri_density_cache_add(density_cache, jatom, iatom, jkind, ikind, cell_rev, -rab, &
                                           rho_block,.NOT. source_transposed, factor)
         END IF
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      ! Match the block-norm filtering used by the extended RI-HFXk density conversion. This is
      ! the user's EPS_FILTER accuracy control, not an additional LRI-specific approximation.
      DO iblock = 1, density_cache%nblocks
         density_cache%blocks(iblock)%norm = &
            SQRT(SUM(density_cache%blocks(iblock)%values**2))
         IF (density_cache%blocks(iblock)%norm < filter_eps) THEN
            DEALLOCATE (density_cache%blocks(iblock)%values)
            density_cache%blocks(iblock)%norm = 0.0_dp
         END IF
      END DO
   END SUBROUTINE hfx_lri_density_cache_create

! **************************************************************************************************
!> \brief Finds a logical block in the sparse real-space density cache.
!> \param density_cache ...
!> \param iatom ...
!> \param jatom ...
!> \param cell ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION hfx_lri_density_cache_find(density_cache, iatom, jatom, cell) RESULT(iblock_found)

      TYPE(lri_density_cache_type), INTENT(IN)           :: density_cache
      INTEGER, INTENT(IN)                                :: iatom, jatom
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell

      INTEGER                                            :: iblock, islot, nprobe

      iblock_found = 0
      IF (.NOT. ALLOCATED(density_cache%slot)) RETURN
      IF (density_cache%capacity <= 0) RETURN
      islot = hfx_lri_pair_hash(iatom, jatom, cell, density_cache%capacity)
      DO nprobe = 1, density_cache%capacity
         iblock = density_cache%slot(islot)
         IF (iblock == 0) RETURN
         IF (density_cache%blocks(iblock)%iatom == iatom .AND. &
             density_cache%blocks(iblock)%jatom == jatom .AND. &
             ALL(density_cache%blocks(iblock)%cell == cell)) THEN
            iblock_found = iblock
            RETURN
         END IF
         islot = MOD(islot, density_cache%capacity) + 1
      END DO

   END FUNCTION hfx_lri_density_cache_find

! **************************************************************************************************
!> \brief Indexes only density blocks that are actually present, grouped by their first atom.
!> \param density_cache ...
!> \param natom ...
!> \param density_index ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_density_index_create(density_cache, natom, density_index)

      TYPE(lri_density_cache_type), INTENT(IN)           :: density_cache
      INTEGER, INTENT(IN)                                :: natom
      TYPE(lri_density_index_type), INTENT(OUT)          :: density_index

      INTEGER                                            :: iatom, iblock, iref
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cursor

      ALLOCATE (density_index%atom_count(natom), SOURCE=0)
      DO iblock = 1, density_cache%nblocks
         IF (.NOT. ALLOCATED(density_cache%blocks(iblock)%values)) CYCLE
         iatom = density_cache%blocks(iblock)%iatom
         IF (iatom < 1 .OR. iatom > natom) CPABORT("Invalid atom in LRI-HFXk density cache")
         density_index%atom_count(iatom) = density_index%atom_count(iatom) + 1
      END DO

      ALLOCATE (density_index%atom_offset(natom + 1))
      density_index%atom_offset(1) = 1
      DO iatom = 1, natom
         density_index%atom_offset(iatom + 1) = density_index%atom_offset(iatom) + &
                                                density_index%atom_count(iatom)
      END DO
      ALLOCATE (density_index%refs(density_index%atom_offset(natom + 1) - 1))
      ALLOCATE (cursor(natom))
      cursor(:) = density_index%atom_offset(1:natom)
      DO iblock = 1, density_cache%nblocks
         IF (.NOT. ALLOCATED(density_cache%blocks(iblock)%values)) CYCLE
         iatom = density_cache%blocks(iblock)%iatom
         iref = cursor(iatom)
         density_index%refs(iref)%iblock = iblock
         density_index%refs(iref)%latom = density_cache%blocks(iblock)%jatom
         cursor(iatom) = iref + 1
      END DO
      DEALLOCATE (cursor)

   END SUBROUTINE hfx_lri_density_index_create

! **************************************************************************************************
!> \brief Builds the exact symmetric set of AO blocks used by the k-point KS matrices. Fitting-pair
!>        locality must not be used to truncate the output K matrix: exchange can couple two local
!>        fitting pairs whose outer AO centers do not overlap.
!> \param sab_nl ...
!> \param cell_to_index ...
!> \param nimages ...
!> \param basis_set_AO ...
!> \param target_pairs ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_target_pairs_create(sab_nl, cell_to_index, nimages, basis_set_AO, target_pairs)

      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      INTEGER, INTENT(IN)                                :: nimages
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_AO
      TYPE(hfx_lri_pair_block_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(OUT)                       :: target_pairs

      INTEGER                                            :: iatom, ikind, image, ipair, jatom, &
                                                            jkind, mepos, ntarget
      INTEGER, DIMENSION(3)                              :: cell
      REAL(dp), DIMENSION(3)                             :: rab
      TYPE(gto_basis_set_type), POINTER                  :: basis_i, basis_j
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      NULLIFY (basis_i, basis_j, nl_iterator)
      CPASSERT(ASSOCIATED(sab_nl))
      CPASSERT(nimages > 0)

      ntarget = 0
      CALL neighbor_list_iterator_create(nl_iterator, sab_nl)
      mepos = 0
      DO WHILE (neighbor_list_iterate(nl_iterator, mepos) == 0)
         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, cell=cell)
         IF (ikind < 1 .OR. ikind > SIZE(basis_set_AO) .OR. &
             jkind < 1 .OR. jkind > SIZE(basis_set_AO)) CYCLE
         basis_i => basis_set_AO(ikind)%gto_basis_set
         basis_j => basis_set_AO(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_i) .OR. .NOT. ASSOCIATED(basis_j)) CYCLE
         image = 1
         IF (ASSOCIATED(cell_to_index)) THEN
            IF (ANY(cell < LBOUND(cell_to_index)) .OR. ANY(cell > UBOUND(cell_to_index))) CYCLE
            image = cell_to_index(cell(1), cell(2), cell(3))
         END IF
         IF (image < 1 .OR. image > nimages) CYCLE
         ntarget = ntarget + 1
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      ALLOCATE (target_pairs(ntarget))
      ipair = 0
      CALL neighbor_list_iterator_create(nl_iterator, sab_nl)
      mepos = 0
      DO WHILE (neighbor_list_iterate(nl_iterator, mepos) == 0)
         CALL get_iterator_info(nl_iterator, mepos=mepos, iatom=iatom, jatom=jatom, &
                                ikind=ikind, jkind=jkind, cell=cell, r=rab)
         IF (ikind < 1 .OR. ikind > SIZE(basis_set_AO) .OR. &
             jkind < 1 .OR. jkind > SIZE(basis_set_AO)) CYCLE
         basis_i => basis_set_AO(ikind)%gto_basis_set
         basis_j => basis_set_AO(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_i) .OR. .NOT. ASSOCIATED(basis_j)) CYCLE
         image = 1
         IF (ASSOCIATED(cell_to_index)) THEN
            IF (ANY(cell < LBOUND(cell_to_index)) .OR. ANY(cell > UBOUND(cell_to_index))) CYCLE
            image = cell_to_index(cell(1), cell(2), cell(3))
         END IF
         IF (image < 1 .OR. image > nimages) CYCLE

         ipair = ipair + 1
         target_pairs(ipair)%iatom = iatom
         target_pairs(ipair)%jatom = jatom
         target_pairs(ipair)%ikind = ikind
         target_pairs(ipair)%jkind = jkind
         target_pairs(ipair)%cell(:) = cell(:)
         target_pairs(ipair)%rab(:) = rab(:)
         target_pairs(ipair)%nao_i = basis_i%nsgf
         target_pairs(ipair)%nao_j = basis_j%nsgf
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)
      CPASSERT(ipair == ntarget)

   END SUBROUTINE hfx_lri_target_pairs_create

! **************************************************************************************************
!> \brief Collects the atom-center metric keys needed by one oriented fitting pair.
!> \param metric_table ...
!> \param pair_left ...
!> \param left_reversed ...
!> \param right_kind ...
!> \param right_position ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_metric_table_add_pair_center(metric_table, pair_left, left_reversed, &
                                                   right_kind, right_position)

      TYPE(hfx_lri_metric_table_type), INTENT(INOUT)     :: metric_table
      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair_left
      LOGICAL, INTENT(IN)                                :: left_reversed
      INTEGER, INTENT(IN)                                :: right_kind
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: right_position

      INTEGER                                            :: ic_left
      REAL(dp), DIMENSION(3)                             :: left_origin, left_position

      left_origin(:) = 0.0_dp
      IF (left_reversed) left_origin(:) = pair_left%rab(:)
      DO ic_left = 1, pair_left%naux_centers
         left_position(:) = pair_left%aux_r(:, ic_left) - left_origin(:)
         CALL hfx_lri_metric_table_add_key( &
            metric_table, pair_left%aux_kind(ic_left), right_kind, &
            right_position - left_position)
      END DO

   END SUBROUTINE hfx_lri_metric_table_add_pair_center

! **************************************************************************************************
!> \brief Enumerates unique two-center metric keys without evaluating any integral.  The two loops
!>        cover respectively the auxiliary center fixed at the K target and the center fixed at L.
!> \param metric_table ...
!> \param pair_blocks ...
!> \param pair_lookup ...
!> \param density_cache ...
!> \param target_pairs ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_metric_table_collect(metric_table, pair_blocks, pair_lookup, &
                                           density_cache, target_pairs)

      TYPE(hfx_lri_metric_table_type), INTENT(INOUT)     :: metric_table
      TYPE(hfx_lri_pair_block_type), DIMENSION(:), &
         INTENT(IN)                                      :: pair_blocks
      TYPE(lri_pair_lookup_type), INTENT(IN)             :: pair_lookup
      TYPE(lri_density_cache_type), INTENT(IN)           :: density_cache
      TYPE(hfx_lri_pair_block_type), DIMENSION(:), &
         INTENT(IN)                                      :: target_pairs

      INTEGER                                            :: iblock, ipair_left, ipair_left_ref, &
                                                            ipair_target, ipos_left, &
                                                            jatom_density, target_iatom
      LOGICAL                                            :: left_reversed
      REAL(dp), DIMENSION(3)                             :: left_rab, right_position

      IF (metric_table%initialized) RETURN

      ! V between either center of C_IJ and the first (K) center of C_KL.  It depends only on the
      ! target I-K and the left pair, not on L.
      DO ipair_target = 1, SIZE(target_pairs)
         target_iatom = target_pairs(ipair_target)%iatom
         DO ipos_left = pair_lookup%atom_offset(target_iatom), &
            pair_lookup%atom_offset(target_iatom + 1) - 1
            ipair_left_ref = pair_lookup%atom_pairs(ipos_left)
            ipair_left = ABS(ipair_left_ref)
            left_reversed = ipair_left_ref < 0
            CALL hfx_lri_metric_table_add_pair_center( &
               metric_table, pair_blocks(ipair_left), left_reversed, &
               target_pairs(ipair_target)%jkind, target_pairs(ipair_target)%rab)
         END DO
      END DO

      ! V to the second (L) center of C_KL.  L is fixed by an I-J fitting pair and a J-L density
      ! block, so this enumeration avoids the target/right-pair Cartesian product entirely.
      DO target_iatom = 1, SIZE(pair_lookup%atom_count)
         DO ipos_left = pair_lookup%atom_offset(target_iatom), &
            pair_lookup%atom_offset(target_iatom + 1) - 1
            ipair_left_ref = pair_lookup%atom_pairs(ipos_left)
            ipair_left = ABS(ipair_left_ref)
            left_reversed = ipair_left_ref < 0
            IF (left_reversed) THEN
               jatom_density = pair_blocks(ipair_left)%iatom
               left_rab(:) = -pair_blocks(ipair_left)%rab(:)
            ELSE
               jatom_density = pair_blocks(ipair_left)%jatom
               left_rab(:) = pair_blocks(ipair_left)%rab(:)
            END IF
            DO iblock = 1, density_cache%nblocks
               IF (density_cache%blocks(iblock)%iatom /= jatom_density) CYCLE
               right_position(:) = left_rab(:) + density_cache%blocks(iblock)%rab(:)
               CALL hfx_lri_metric_table_add_pair_center( &
                  metric_table, pair_blocks(ipair_left), left_reversed, &
                  density_cache%blocks(iblock)%jkind, right_position)
            END DO
         END DO
      END DO

   END SUBROUTINE hfx_lri_metric_table_collect

! **************************************************************************************************
!> \brief Evaluates one canonical table entry with a thread-private integral context.
!> \param metric_entry ...
!> \param basis_set_RI ...
!> \param ri_transform_cache ...
!> \param ri_data ...
!> \param integral_context ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_metric_table_compute_entry(metric_entry, basis_set_RI, ri_transform_cache, &
                                                 ri_data, integral_context)

      TYPE(hfx_lri_metric_cache_entry_type), &
         INTENT(INOUT)                                   :: metric_entry
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(lri_basis_transform_type), DIMENSION(:), &
         INTENT(IN)                                      :: ri_transform_cache
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: integral_context

      INTEGER                                            :: kind_left, kind_right
      REAL(dp)                                           :: screen_distance
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: tmp, v_raw
      REAL(dp), DIMENSION(3)                             :: r0
      TYPE(gto_basis_set_type), POINTER                  :: basis_left, basis_right

      NULLIFY (basis_left, basis_right)
      kind_left = metric_entry%kind_left
      kind_right = metric_entry%kind_right
      r0(:) = 0.0_dp
      metric_entry%zero = .FALSE.
      metric_entry%norm2 = 0.0_dp
      screen_distance = ri_transform_cache(kind_left)%radius + &
                        ri_transform_cache(kind_right)%radius + &
                        coulomb_cutoff_screen_factor*ri_data%hfx_pot%cutoff_radius
      IF (SUM(metric_entry%displacement**2) > screen_distance*screen_distance) THEN
         metric_entry%zero = .TRUE.
         RETURN
      END IF

      basis_left => basis_set_RI(kind_left)%gto_basis_set
      basis_right => basis_set_RI(kind_right)%gto_basis_set
      IF (.NOT. ASSOCIATED(basis_left) .OR. .NOT. ASSOCIATED(basis_right) .OR. &
          .NOT. ALLOCATED(ri_transform_cache(kind_left)%sphi) .OR. &
          .NOT. ALLOCATED(ri_transform_cache(kind_right)%sphi)) THEN
         CPABORT("Missing basis data while building the LRI-HFXk metric table")
      END IF

      ALLOCATE (v_raw(SIZE(ri_transform_cache(kind_left)%sphi, 1), &
                      SIZE(ri_transform_cache(kind_right)%sphi, 1)), SOURCE=0.0_dp)
      CALL hfx_lri_fill_2c_block(basis_left, basis_right, r0, metric_entry%displacement, &
                                 integral_context, ri_data, v_raw, 0, 0)
      ALLOCATE (tmp(SIZE(v_raw, 1), SIZE(ri_transform_cache(kind_right)%sphi, 2)))
      tmp(:, :) = MATMUL(v_raw, ri_transform_cache(kind_right)%sphi)
      ALLOCATE (metric_entry%values(SIZE(ri_transform_cache(kind_left)%sphi, 2), &
                                    SIZE(ri_transform_cache(kind_right)%sphi, 2)))
      metric_entry%values(:, :) = MATMUL(TRANSPOSE(ri_transform_cache(kind_left)%sphi), tmp)
      DEALLOCATE (tmp, v_raw)
      metric_entry%norm2 = SUM(metric_entry%values**2)
      IF (metric_entry%norm2 == 0.0_dp) THEN
         DEALLOCATE (metric_entry%values)
         metric_entry%zero = .TRUE.
      END IF

   END SUBROUTINE hfx_lri_metric_table_compute_entry

! **************************************************************************************************
!> \brief Builds the immutable rank-shared metric table in parallel and reports its actual storage.
!> \param metric_table ...
!> \param pair_blocks ...
!> \param pair_lookup ...
!> \param density_cache ...
!> \param target_pairs ...
!> \param basis_set_RI ...
!> \param ri_transform_cache ...
!> \param ri_data ...
!> \param maxl_ri ...
!> \param stored_entries ...
!> \param screened_entries ...
!> \param storage_mib ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_metric_table_build(metric_table, pair_blocks, pair_lookup, density_cache, &
                                         target_pairs, basis_set_RI, ri_transform_cache, ri_data, &
                                         maxl_ri, stored_entries, screened_entries, storage_mib)

      TYPE(hfx_lri_metric_table_type), INTENT(INOUT)     :: metric_table
      TYPE(hfx_lri_pair_block_type), DIMENSION(:), &
         INTENT(IN)                                      :: pair_blocks
      TYPE(lri_pair_lookup_type), INTENT(IN)             :: pair_lookup
      TYPE(lri_density_cache_type), INTENT(IN)           :: density_cache
      TYPE(hfx_lri_pair_block_type), DIMENSION(:), &
         INTENT(IN)                                      :: target_pairs
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(lri_basis_transform_type), DIMENSION(:), &
         INTENT(IN)                                      :: ri_transform_cache
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      INTEGER, INTENT(IN)                                :: maxl_ri
      INTEGER, INTENT(OUT)                               :: stored_entries, screened_entries
      REAL(dp), INTENT(OUT)                              :: storage_mib

      INTEGER                                            :: ientry
      INTEGER(int_8)                                     :: stored_reals
      TYPE(coulomb_integral_context_type)                :: integral_context
      TYPE(hfx_lri_metric_cache_entry_type), &
         ALLOCATABLE, DIMENSION(:)                       :: compact_entries

      stored_entries = 0
      screened_entries = 0
      stored_reals = 0_int_8
      storage_mib = 0.0_dp
      IF (.NOT. metric_table%initialized) THEN
         CALL hfx_lri_metric_table_collect( &
            metric_table, pair_blocks, pair_lookup, density_cache, target_pairs)
         IF (metric_table%nentries <= 0) CPABORT("LRI-HFXk collected no metric-table entries")
         IF (SIZE(metric_table%entries) /= metric_table%nentries) THEN
            ALLOCATE (compact_entries(metric_table%nentries))
            compact_entries(:) = metric_table%entries(1:metric_table%nentries)
            CALL MOVE_ALLOC(compact_entries, metric_table%entries)
         END IF

         !$OMP PARALLEL DEFAULT(NONE) &
         !$OMP          SHARED(basis_set_RI, maxl_ri, metric_table, ri_data, ri_transform_cache) &
         !$OMP          PRIVATE(ientry, integral_context)
         CALL coulomb_integral_init(integral_context, max_am_2c=maxl_ri)
         !$OMP DO SCHEDULE(DYNAMIC, 1)
         DO ientry = 1, metric_table%nentries
            CALL hfx_lri_metric_table_compute_entry( &
               metric_table%entries(ientry), basis_set_RI, ri_transform_cache, &
               ri_data, integral_context)
         END DO
         !$OMP END DO
         CALL coulomb_integral_cleanup(integral_context)
         !$OMP END PARALLEL
         metric_table%initialized = .TRUE.
      END IF

      DO ientry = 1, metric_table%nentries
         IF (metric_table%entries(ientry)%zero) THEN
            screened_entries = screened_entries + 1
         ELSE IF (ALLOCATED(metric_table%entries(ientry)%values)) THEN
            stored_entries = stored_entries + 1
            stored_reals = stored_reals + INT(SIZE(metric_table%entries(ientry)%values), int_8)
         END IF
      END DO
      storage_mib = REAL(stored_reals, dp)*8.0_dp/1024.0_dp**2

   END SUBROUTINE hfx_lri_metric_table_build

! **************************************************************************************************
!> \brief Accumulates an exchange K block for one output pair.
!> \param basis_set_RI ...
!> \param pair_lookup indexed logical fitting-pair views
!> \param density_cache sparse, desymmetrized real-space density
!> \param density_index nonzero density blocks grouped by their first atom
!> \param ri_transform_cache ...
!> \param integral_context ...
!> \param exchange_workspace ...
!> \param pair_target requested K_IK lattice block
!> \param ri_data ...
!> \param k_accum_block ...
!> \param used_pairs ...
!> \param candidate_pairs ...
!> \param rho_accum_norm ...
!> \param skipped_image ...
!> \param skipped_cell_bounds ...
!> \param skipped_img_index ...
!> \param skipped_missing_rho ...
!> \param skipped_rho_dim ...
!> \param skipped_small_rho ...
!> \param skipped_no_kblock ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_accumulate_target_k_block(basis_set_RI, pair_lookup, density_cache, density_index, &
                                                ri_transform_cache, integral_context, exchange_workspace, &
                                                pair_target, ri_data, &
                                                k_accum_block, used_pairs, &
                                                candidate_pairs, rho_accum_norm, skipped_image, &
                                                skipped_cell_bounds, skipped_img_index, &
                                                skipped_missing_rho, skipped_rho_dim, skipped_small_rho, &
                                                skipped_no_kblock)

      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(lri_pair_lookup_type), INTENT(IN)             :: pair_lookup
      TYPE(lri_density_cache_type), INTENT(IN)           :: density_cache
      TYPE(lri_density_index_type), INTENT(IN)           :: density_index
      TYPE(lri_basis_transform_type), DIMENSION(:), &
         INTENT(IN)                                      :: ri_transform_cache
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: integral_context
      TYPE(hfx_lri_exchange_workspace_type), &
         INTENT(INOUT)                                   :: exchange_workspace
      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair_target
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: k_accum_block
      INTEGER, INTENT(OUT)                               :: used_pairs, candidate_pairs
      REAL(dp), INTENT(OUT)                              :: rho_accum_norm
      INTEGER, INTENT(OUT) :: skipped_image, skipped_cell_bounds, skipped_img_index, &
         skipped_missing_rho, skipped_rho_dim, skipped_small_rho, skipped_no_kblock

      INTEGER :: iblock, indexed_pairs, ipair_left_ref, ipair_right_ref, ipos_left, ipos_right, &
         iref, jatom_density, latom_density, left_iatom, right_iatom, right_jatom, target_iatom, &
         target_jatom
      INTEGER, DIMENSION(3)                              :: cell_density, cell_left, cell_right
      LOGICAL                                            :: contraction_done, outer_initialized, &
                                                            outer_ready, outer_sum_used
      REAL(dp)                                           :: outer_norm2, rho_norm2_accum, &
                                                            rho_pair_norm

      used_pairs = 0
      candidate_pairs = 0
      skipped_image = 0
      skipped_cell_bounds = 0
      skipped_img_index = 0
      skipped_missing_rho = 0
      skipped_rho_dim = 0
      skipped_small_rho = 0
      skipped_no_kblock = 0
      rho_accum_norm = 0.0_dp

      IF (.NOT. ASSOCIATED(ri_data%lri_kp)) RETURN
      IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks)) RETURN
      IF (.NOT. ALLOCATED(pair_lookup%atom_count)) RETURN
      IF (.NOT. ALLOCATED(pair_lookup%atom_offset)) RETURN
      IF (.NOT. ALLOCATED(pair_lookup%atom_pairs)) RETURN

      ALLOCATE (k_accum_block(pair_target%nao_i, pair_target%nao_j), SOURCE=0.0_dp)
      rho_norm2_accum = 0.0_dp
      indexed_pairs = 0
      target_iatom = pair_target%iatom
      target_jatom = pair_target%jatom

      ! For K_IK(R_IK), combine the fitting pairs I-J(R_IJ) and K-L(R_KL)
      ! with D_JL(R_JL), where R_JL = R_IK + R_KL - R_IJ. The second
      ! fitting pair is translated by R_IK when its auxiliary metric is built.
      DO ipos_left = pair_lookup%atom_offset(target_iatom), pair_lookup%atom_offset(target_iatom + 1) - 1
         ipair_left_ref = pair_lookup%atom_pairs(ipos_left)
         CALL hfx_lri_pair_ref_key(ri_data%lri_kp%pair_blocks, ipair_left_ref, &
                                   left_iatom, jatom_density, cell_left)
         CPASSERT(left_iatom == target_iatom)
         outer_initialized = .FALSE.
         outer_ready = .FALSE.
         outer_sum_used = .FALSE.
         outer_norm2 = 0.0_dp

         ! Both traversals enumerate the identical sparse join
         !   C_IJ(a) D_JL(c) C_KL(b), c = target + b - a.
         ! Iterate whichever adjacency is shorter and find the other operand by hash.
         IF (density_index%atom_count(jatom_density) <= pair_lookup%atom_count(target_jatom)) THEN
            candidate_pairs = candidate_pairs + density_index%atom_count(jatom_density)
            DO iref = density_index%atom_offset(jatom_density), &
               density_index%atom_offset(jatom_density + 1) - 1
               iblock = density_index%refs(iref)%iblock
               latom_density = density_index%refs(iref)%latom
               cell_right(:) = hfx_lri_right_cell( &
                               density_cache%blocks(iblock)%cell, &
                               pair_target%cell, cell_left)
               ipair_right_ref = hfx_lri_pair_lookup_find( &
                                 pair_lookup, ri_data%lri_kp%pair_blocks, &
                                 target_jatom, latom_density, cell_right)
               IF (ipair_right_ref == 0) CYCLE
               indexed_pairs = indexed_pairs + 1

               CALL hfx_lri_accumulate_matched_quartet( &
                  basis_set_RI, ri_data%lri_kp%pair_blocks, pair_target, &
                  ipair_left_ref, ipair_right_ref, density_cache%blocks(iblock), ri_data, &
                  integral_context, ri_transform_cache, exchange_workspace, k_accum_block, &
                  outer_initialized, outer_ready, outer_sum_used, outer_norm2, &
                  contraction_done)
               IF (.NOT. contraction_done) THEN
                  skipped_no_kblock = skipped_no_kblock + 1
                  CYCLE
               END IF
               rho_pair_norm = density_cache%blocks(iblock)%norm
               rho_norm2_accum = rho_norm2_accum + rho_pair_norm*rho_pair_norm
               used_pairs = used_pairs + 1
            END DO
         ELSE
            candidate_pairs = candidate_pairs + pair_lookup%atom_count(target_jatom)
            DO ipos_right = pair_lookup%atom_offset(target_jatom), &
               pair_lookup%atom_offset(target_jatom + 1) - 1
               ipair_right_ref = pair_lookup%atom_pairs(ipos_right)
               CALL hfx_lri_pair_ref_key(ri_data%lri_kp%pair_blocks, ipair_right_ref, &
                                         right_iatom, right_jatom, cell_right)
               CPASSERT(right_iatom == target_jatom)
               cell_density(:) = hfx_lri_density_cell(pair_target%cell, cell_left, cell_right)
               iblock = hfx_lri_density_cache_find( &
                        density_cache, jatom_density, right_jatom, cell_density)
               IF (iblock == 0) CYCLE
               indexed_pairs = indexed_pairs + 1

               CALL hfx_lri_accumulate_matched_quartet( &
                  basis_set_RI, ri_data%lri_kp%pair_blocks, pair_target, &
                  ipair_left_ref, ipair_right_ref, density_cache%blocks(iblock), ri_data, &
                  integral_context, ri_transform_cache, exchange_workspace, k_accum_block, &
                  outer_initialized, outer_ready, outer_sum_used, outer_norm2, &
                  contraction_done)
               IF (.NOT. contraction_done) THEN
                  skipped_no_kblock = skipped_no_kblock + 1
                  CYCLE
               END IF
               rho_pair_norm = density_cache%blocks(iblock)%norm
               rho_norm2_accum = rho_norm2_accum + rho_pair_norm*rho_pair_norm
               used_pairs = used_pairs + 1
            END DO
         END IF
         IF (outer_sum_used) THEN
            CALL hfx_lri_contract_outer_finalize(k_accum_block, exchange_workspace%contract)
         END IF
      END DO

      skipped_cell_bounds = MAX(0, candidate_pairs - indexed_pairs)
      skipped_image = skipped_cell_bounds
      rho_accum_norm = SQRT(rho_norm2_accum)
      IF (used_pairs <= 0) THEN
         IF (ALLOCATED(k_accum_block)) DEALLOCATE (k_accum_block)
      END IF

   END SUBROUTINE hfx_lri_accumulate_target_k_block

! **************************************************************************************************
!> \brief Validates and accumulates one matched C-D-C sparse-join entry.
!> \param basis_set_RI ...
!> \param pair_blocks ...
!> \param pair_target ...
!> \param ipair_left_ref ...
!> \param ipair_right_ref ...
!> \param density_block ...
!> \param ri_data ...
!> \param integral_context ...
!> \param ri_transform_cache ...
!> \param exchange_workspace ...
!> \param k_accum_block ...
!> \param outer_initialized ...
!> \param outer_ready ...
!> \param outer_sum_used ...
!> \param outer_norm2 ...
!> \param contraction_done ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_accumulate_matched_quartet( &
      basis_set_RI, pair_blocks, pair_target, ipair_left_ref, ipair_right_ref, density_block, &
      ri_data, integral_context, ri_transform_cache, exchange_workspace, k_accum_block, &
      outer_initialized, outer_ready, outer_sum_used, outer_norm2, &
      contraction_done)

      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(hfx_lri_pair_block_type), DIMENSION(:), &
         INTENT(IN)                                      :: pair_blocks
      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair_target
      INTEGER, INTENT(IN)                                :: ipair_left_ref, ipair_right_ref
      TYPE(lri_density_block_type), INTENT(IN)           :: density_block
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: integral_context
      TYPE(lri_basis_transform_type), DIMENSION(:), &
         INTENT(IN)                                      :: ri_transform_cache
      TYPE(hfx_lri_exchange_workspace_type), &
         INTENT(INOUT)                                   :: exchange_workspace
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: k_accum_block
      LOGICAL, INTENT(INOUT)                             :: outer_initialized, outer_ready, &
                                                            outer_sum_used
      REAL(dp), INTENT(INOUT)                            :: outer_norm2
      LOGICAL, INTENT(OUT)                               :: contraction_done

      INTEGER                                            :: ipair_left, ipair_right, jatom_density, &
                                                            latom_density, left_iatom, left_nao_j, &
                                                            right_iatom, right_nao_j
      INTEGER, DIMENSION(3)                              :: cell_density, cell_left, cell_right
      LOGICAL                                            :: left_reversed, right_reversed

      contraction_done = .FALSE.
      ipair_left = ABS(ipair_left_ref)
      ipair_right = ABS(ipair_right_ref)
      IF (ipair_left < 1 .OR. ipair_left > SIZE(pair_blocks) .OR. &
          ipair_right < 1 .OR. ipair_right > SIZE(pair_blocks)) THEN
         CPABORT("Invalid signed pair reference in LRI-HFXk sparse contraction")
      END IF
      left_reversed = ipair_left_ref < 0
      right_reversed = ipair_right_ref < 0
      CALL hfx_lri_pair_ref_key(pair_blocks, ipair_left_ref, left_iatom, jatom_density, cell_left)
      CALL hfx_lri_pair_ref_key(pair_blocks, ipair_right_ref, right_iatom, latom_density, cell_right)
      cell_density(:) = hfx_lri_density_cell(pair_target%cell, cell_left, cell_right)
      IF (left_iatom /= pair_target%iatom .OR. right_iatom /= pair_target%jatom .OR. &
          density_block%iatom /= jatom_density .OR. density_block%jatom /= latom_density .OR. &
          ANY(density_block%cell /= cell_density)) THEN
         CPABORT("Inconsistent atom/cell topology in LRI-HFXk sparse contraction")
      END IF
      IF (.NOT. ALLOCATED(density_block%values)) THEN
         CPABORT("Missing density values in LRI-HFXk sparse contraction")
      END IF

      left_nao_j = MERGE(pair_blocks(ipair_left)%nao_i, &
                         pair_blocks(ipair_left)%nao_j, left_reversed)
      right_nao_j = MERGE(pair_blocks(ipair_right)%nao_i, &
                          pair_blocks(ipair_right)%nao_j, right_reversed)
      IF (ANY(SHAPE(density_block%values) /= [left_nao_j, right_nao_j])) THEN
         CPABORT("Inconsistent density dimensions in LRI-HFXk sparse contraction")
      END IF
      IF (density_block%norm == 0.0_dp) RETURN

      CALL hfx_lri_accumulate_pair_pair_density_realv( &
         basis_set_RI, pair_target, pair_blocks(ipair_left), left_reversed, &
         pair_blocks(ipair_right), right_reversed, ri_data, density_block%values, &
         k_accum_block, integral_context, ri_transform_cache, exchange_workspace, &
         outer_initialized, outer_ready, outer_sum_used, outer_norm2, &
         contraction_done)

   END SUBROUTINE hfx_lri_accumulate_matched_quartet

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ri_data ...
!> \param rho_ao ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_build_rank_local_k_blocks(qs_env, ri_data, rho_ao)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_lri_build_rank_local_k_blocks'

      INTEGER :: candidate_pairs, candidate_pairs_tmp, filtered_k_blocks, general_candidates, &
         general_density_pairs, general_k_active_high, general_k_active_large, &
         general_k_active_low, general_k_active_mid, general_rho_active, general_targets, &
         general_tests, handle, ikind, img, img_pair, img_rev, ipair, ipair_best, &
         ipair_left_ref_est, ipair_sel, ipair_target, ithread, lb1, lb2, lb3, left_second_atom, &
         max_k_blocks, maxl_ri, maxl_ri_kind, metric_table_screened, metric_table_stored, &
         multi_candidate_pairs, multi_skipped_cell_bounds, multi_skipped_image, &
         multi_skipped_img_index
      INTEGER :: multi_skipped_missing_rho, multi_skipped_no_kblock, multi_skipped_rho_dim, &
         multi_skipped_small_rho, multi_used_pairs, natom, nearest_candidates, &
         nearest_k_active_high, nearest_k_active_large, nearest_k_active_low, &
         nearest_k_active_mid, nearest_rho_active, nearest_targets, nearest_tests, &
         nk_blocks_built, nkind, nthreads, rho_energy_block, search_iter, skipped_cell_bounds, &
         skipped_cell_bounds_tmp, skipped_image, skipped_image_tmp, skipped_img_index, &
         skipped_img_index_tmp, skipped_missing_rho, skipped_missing_rho_tmp, skipped_no_kblock, &
         skipped_no_kblock_tmp, skipped_rho_dim
      INTEGER :: skipped_rho_dim_tmp, skipped_small_rho, skipped_small_rho_tmp, ub1, ub2, ub3, &
         used_pairs, used_pairs_tmp
      INTEGER(int_8) :: metric_cache_hits, metric_cache_misses, metric_cache_screened, &
         pair_screened, scan_candidates_est, target_left_est
      INTEGER, DIMENSION(3)                              :: cell_rev
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found, found_rev
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: tried_density_pair
      REAL(dp) :: best_r2, c_norm, exchange_trace_raw, general_k_max, general_k_rms, k_accum_norm, &
         k_blocks_norm, k_norm, k_tmp_norm, metric_table_mib, nearest_k_max, nearest_k_rms, &
         pair_r2, rho_accum_norm, rho_blocks_norm, rho_norm, rho_norm2_accum, rho_pair_norm, &
         rho_tmp_norm, trace_weight
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: k_accum_block, k_block, k_pair_block, &
                                                            k_tmp_block, rho_work
      REAL(dp), DIMENSION(:, :), POINTER                 :: rho_block, rho_block_rev
      TYPE(coulomb_integral_context_type)                :: integral_context
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_AO, basis_set_RI
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(hfx_lri_exchange_workspace_type)              :: exchange_workspace
      TYPE(hfx_lri_pair_block_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: target_pairs
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(lri_basis_transform_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: ri_transform_cache
      TYPE(lri_density_cache_type)                       :: density_cache
      TYPE(lri_density_index_type)                       :: density_index
      TYPE(lri_pair_lookup_type)                         :: pair_lookup
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (basis_set_AO, basis_set_RI, cell_to_index, dft_control, kpoints, qs_kind_set, rho_block, &
               rho_block_rev, ri_basis, sab_nl)

      IF (.NOT. ASSOCIATED(ri_data%lri_kp)) CPABORT("LRI-HFXk state is not associated")
      IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks)) THEN
         CPABORT("LRI-HFXk fitting-pair blocks are not allocated")
      END IF
      IF (ri_data%lri_kp%npair_blocks <= 0) CPABORT("LRI-HFXk has no fitting-pair blocks")
      IF (.NOT. ASSOCIATED(rho_ao)) CPABORT("LRI-HFXk density matrix is not associated")

      CALL timeset(routineN, handle)

      ipair_sel = 1
      IF (lri_debug_level >= 1) THEN
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
            CALL timestop(handle)
            RETURN
         END IF
      END IF

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, kpoints=kpoints, &
                      qs_kind_set=qs_kind_set, nkind=nkind, natom=natom)
      IF (ASSOCIATED(kpoints)) THEN
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index, sab_nl=sab_nl)
      END IF

      ALLOCATE (basis_set_AO(nkind), basis_set_RI(nkind))
      CALL basis_set_list_setup(basis_set_AO, ri_data%orb_basis_type, qs_kind_set)
      CALL basis_set_list_setup(basis_set_RI, ri_data%ri_basis_type, qs_kind_set)
      ! The basis radii are mutable shared metadata. Use the LRI threshold during every Libint call
      ! and restore the ordinary QS value before returning.
      DO ikind = 1, nkind
         ri_basis => basis_set_RI(ikind)%gto_basis_set
         IF (ASSOCIATED(ri_basis)) THEN
            CALL init_interaction_radii_orb_basis(ri_basis, ri_data%eps_pgf_orb)
         END IF
      END DO

      CALL hfx_lri_pair_lookup_create(pair_lookup, ri_data%lri_kp%pair_blocks, &
                                      ri_data%lri_kp%npair_blocks, natom)
      CALL hfx_lri_density_cache_create(rho_ao, cell_to_index, sab_nl, ri_data%filter_eps, density_cache)
      CALL hfx_lri_density_index_create(density_cache, natom, density_index)
      CALL hfx_lri_target_pairs_create(sab_nl, cell_to_index, SIZE(rho_ao, 2), &
                                       basis_set_AO, target_pairs)
      IF (SIZE(target_pairs) == 0) THEN
         CPABORT("LRI-HFXk found no writable target blocks in the k-point neighbor list")
      END IF
      CALL hfx_lri_build_ri_transform_cache(basis_set_RI, ri_transform_cache)

      nthreads = 1
!$    nthreads = omp_get_max_threads()
      IF (ALLOCATED(ri_data%lri_kp%thread_workspaces)) THEN
         IF (SIZE(ri_data%lri_kp%thread_workspaces) /= nthreads) THEN
            DEALLOCATE (ri_data%lri_kp%thread_workspaces)
         END IF
      END IF
      IF (.NOT. ALLOCATED(ri_data%lri_kp%thread_workspaces)) THEN
         ALLOCATE (ri_data%lri_kp%thread_workspaces(nthreads))
      END IF

      maxl_ri = 0
      DO ikind = 1, nkind
         ri_basis => basis_set_RI(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(ri_basis)) CYCLE
         CALL get_gto_basis_set(ri_basis, maxl=maxl_ri_kind)
         maxl_ri = MAX(maxl_ri, maxl_ri_kind)
      END DO
      CALL init_md_ftable(nmax=2*maxl_ri)

      IF (ri_data%unit_nr > 0 .AND. .NOT. ri_data%lri_kp%metric_table%initialized) THEN
         WRITE (ri_data%unit_nr, "(/,T2,A)") "LRI-HFXk shared metric setup"
         WRITE (ri_data%unit_nr, "(T4,A,T28,I10)") "target blocks", SIZE(target_pairs)
         WRITE (ri_data%unit_nr, "(T4,A,T28,I10)") "density blocks", density_cache%nblocks
         WRITE (ri_data%unit_nr, "(T4,A,T28,I10)") "active density", SIZE(density_index%refs)
         FLUSH (ri_data%unit_nr)
      END IF
      CALL hfx_lri_metric_table_build( &
         ri_data%lri_kp%metric_table, ri_data%lri_kp%pair_blocks, pair_lookup, density_cache, &
         target_pairs, basis_set_RI, ri_transform_cache, ri_data, maxl_ri, &
         metric_table_stored, metric_table_screened, metric_table_mib)
      ! C construction used the small private caches.  Once the complete immutable table exists,
      ! discard those replicas; K construction normally performs no private cache allocation.
      DO ithread = 1, SIZE(ri_data%lri_kp%thread_workspaces)
         IF (ALLOCATED(ri_data%lri_kp%thread_workspaces(ithread)%metric_cache%entries)) THEN
            DEALLOCATE (ri_data%lri_kp%thread_workspaces(ithread)%metric_cache%entries)
         END IF
      END DO
      IF (ri_data%unit_nr > 0) THEN
         WRITE (ri_data%unit_nr, "(T4,A,T28,I10)") &
            "unique metric keys", ri_data%lri_kp%metric_table%nentries
         WRITE (ri_data%unit_nr, "(T4,A,T28,I10)") "stored metric blocks", metric_table_stored
         WRITE (ri_data%unit_nr, "(T4,A,T28,I10)") "distance-screened", metric_table_screened
         WRITE (ri_data%unit_nr, "(T4,A,T28,F10.2)") "metric storage/MiB", metric_table_mib
         FLUSH (ri_data%unit_nr)
      END IF

      IF (ALLOCATED(ri_data%lri_kp%diagnostic_k_block%k_block)) THEN
         DEALLOCATE (ri_data%lri_kp%diagnostic_k_block%k_block)
      END IF
      IF (lri_debug_level >= 1) THEN
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
            IF (ALLOCATED(target_pairs)) DEALLOCATE (target_pairs)
            IF (ASSOCIATED(dft_control)) THEN
               DO ikind = 1, nkind
                  ri_basis => basis_set_RI(ikind)%gto_basis_set
                  IF (ASSOCIATED(ri_basis)) THEN
                     CALL init_interaction_radii_orb_basis(ri_basis, dft_control%qs_control%eps_pgf_orb)
                  END IF
               END DO
            END IF
            IF (ASSOCIATED(basis_set_AO)) DEALLOCATE (basis_set_AO)
            IF (ASSOCIATED(basis_set_RI)) DEALLOCATE (basis_set_RI)
            CALL timestop(handle)
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
         skipped_image = 0
         skipped_cell_bounds = 0
         skipped_img_index = 0
         skipped_missing_rho = 0
         skipped_rho_dim = 0
         skipped_small_rho = 0
         skipped_no_kblock = 0
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
                          ALL(ri_data%lri_kp%pair_blocks(ipair_best)%cell <= [ub1, ub2, ub3]))) THEN
                  skipped_cell_bounds = skipped_cell_bounds + 1
                  skipped_image = skipped_image + 1
                  CYCLE
               END IF
               img_pair = cell_to_index(ri_data%lri_kp%pair_blocks(ipair_best)%cell(1), &
                                        ri_data%lri_kp%pair_blocks(ipair_best)%cell(2), &
                                        ri_data%lri_kp%pair_blocks(ipair_best)%cell(3))
            END IF
            IF (img_pair < 1 .OR. img_pair > SIZE(rho_ao, 2)) THEN
               skipped_img_index = skipped_img_index + 1
               skipped_image = skipped_image + 1
               CYCLE
            END IF

            CALL dbcsr_get_block_p(matrix=rho_ao(1, img_pair)%matrix, &
                                   row=ri_data%lri_kp%pair_blocks(ipair_best)%iatom, &
                                   col=ri_data%lri_kp%pair_blocks(ipair_best)%jatom, &
                                   block=rho_block, found=found)
            IF (ALLOCATED(rho_work)) DEALLOCATE (rho_work)
            IF (found) THEN
               IF (SIZE(rho_block, 1) < ri_data%lri_kp%pair_blocks(ipair_best)%nao_i .OR. &
                   SIZE(rho_block, 2) < ri_data%lri_kp%pair_blocks(ipair_best)%nao_j) THEN
                  skipped_rho_dim = skipped_rho_dim + 1
                  CYCLE
               END IF
               ALLOCATE (rho_work(ri_data%lri_kp%pair_blocks(ipair_best)%nao_i, &
                                  ri_data%lri_kp%pair_blocks(ipair_best)%nao_j))
               rho_work(:, :) = rho_block(1:ri_data%lri_kp%pair_blocks(ipair_best)%nao_i, &
                                          1:ri_data%lri_kp%pair_blocks(ipair_best)%nao_j)
            ELSE
               ! In the image representation the transposed block may be stored as D_ji(-R)^T.
               ! Keep the single-debug accumulation path consistent with the multi-target K builder.
               found_rev = .FALSE.
               cell_rev(:) = -ri_data%lri_kp%pair_blocks(ipair_best)%cell(:)
               img_rev = 0
               IF (ASSOCIATED(cell_to_index)) THEN
                  IF (ALL(cell_rev >= [lb1, lb2, lb3]) .AND. ALL(cell_rev <= [ub1, ub2, ub3])) THEN
                     img_rev = cell_to_index(cell_rev(1), cell_rev(2), cell_rev(3))
                  END IF
               ELSE
                  img_rev = 1
               END IF
               IF (img_rev >= 1 .AND. img_rev <= SIZE(rho_ao, 2)) THEN
                  CALL dbcsr_get_block_p(matrix=rho_ao(1, img_rev)%matrix, &
                                         row=ri_data%lri_kp%pair_blocks(ipair_best)%jatom, &
                                         col=ri_data%lri_kp%pair_blocks(ipair_best)%iatom, &
                                         block=rho_block_rev, found=found_rev)
               END IF
               IF (.NOT. found_rev) THEN
                  skipped_missing_rho = skipped_missing_rho + 1
                  CYCLE
               END IF
               IF (SIZE(rho_block_rev, 1) < ri_data%lri_kp%pair_blocks(ipair_best)%nao_j .OR. &
                   SIZE(rho_block_rev, 2) < ri_data%lri_kp%pair_blocks(ipair_best)%nao_i) THEN
                  skipped_rho_dim = skipped_rho_dim + 1
                  CYCLE
               END IF
               ALLOCATE (rho_work(ri_data%lri_kp%pair_blocks(ipair_best)%nao_i, &
                                  ri_data%lri_kp%pair_blocks(ipair_best)%nao_j))
               rho_work(:, :) = TRANSPOSE(rho_block_rev(1:ri_data%lri_kp%pair_blocks(ipair_best)%nao_j, &
                                                        1:ri_data%lri_kp%pair_blocks(ipair_best)%nao_i))
            END IF
            rho_pair_norm = SQRT(SUM(rho_work*rho_work))
            IF (rho_pair_norm < lri_debug_rho_threshold) THEN
               skipped_small_rho = skipped_small_rho + 1
               CYCLE
            END IF
            IF (ri_data%lri_kp%pair_blocks(ipair_best)%iatom /= &
                ri_data%lri_kp%pair_blocks(ipair_sel)%jatom) CYCLE

            CALL hfx_lri_contract_pair_pair_density_realv(basis_set_RI, &
                                                          ri_data%lri_kp%pair_blocks(ipair_sel), &
                                                          ri_data%lri_kp%pair_blocks(ipair_sel), &
                                                          ri_data%lri_kp%pair_blocks(ipair_best), &
                                                          ri_data, rho_work, k_pair_block)
            IF (.NOT. ALLOCATED(k_pair_block)) THEN
               skipped_no_kblock = skipped_no_kblock + 1
               CYCLE
            END IF
            k_accum_block(:, :) = k_accum_block(:, :) + k_pair_block(:, :)
            rho_norm2_accum = rho_norm2_accum + rho_pair_norm*rho_pair_norm
            used_pairs = used_pairs + 1
            DEALLOCATE (k_pair_block)
         END DO
         IF (ALLOCATED(rho_work)) DEALLOCATE (rho_work)
         IF (ALLOCATED(tried_density_pair)) DEALLOCATE (tried_density_pair)
         rho_accum_norm = SQRT(rho_norm2_accum)
         k_accum_norm = SQRT(SUM(k_accum_block*k_accum_block))

         CALL hfx_lri_store_k_block(ri_data%lri_kp%diagnostic_k_block, &
                                    ri_data%lri_kp%pair_blocks(ipair_sel), &
                                    ipair_sel, k_accum_block)
      END IF

      IF (ALLOCATED(ri_data%lri_kp%k_blocks)) THEN
         DO ipair = 1, SIZE(ri_data%lri_kp%k_blocks)
            IF (ALLOCATED(ri_data%lri_kp%k_blocks(ipair)%k_block)) THEN
               DEALLOCATE (ri_data%lri_kp%k_blocks(ipair)%k_block)
            END IF
         END DO
         DEALLOCATE (ri_data%lri_kp%k_blocks)
      END IF
      ri_data%lri_kp%nk_blocks = 0
      max_k_blocks = SIZE(target_pairs)
      IF (max_k_blocks > 0) THEN
         ALLOCATE (ri_data%lri_kp%k_blocks(max_k_blocks))
      END IF
      multi_used_pairs = 0
      multi_candidate_pairs = 0
      multi_skipped_image = 0
      multi_skipped_cell_bounds = 0
      multi_skipped_img_index = 0
      multi_skipped_missing_rho = 0
      multi_skipped_rho_dim = 0
      multi_skipped_small_rho = 0
      multi_skipped_no_kblock = 0
      nk_blocks_built = 0
      filtered_k_blocks = 0
      metric_cache_hits = 0_int_8
      metric_cache_misses = 0_int_8
      metric_cache_screened = 0_int_8
      pair_screened = 0_int_8
      rho_blocks_norm = 0.0_dp
      k_blocks_norm = 0.0_dp

      scan_candidates_est = 0_int_8
      target_left_est = 0_int_8
      DO ipair_target = 1, SIZE(target_pairs)
         DO ipair = pair_lookup%atom_offset(target_pairs(ipair_target)%iatom), &
            pair_lookup%atom_offset(target_pairs(ipair_target)%iatom + 1) - 1
            ipair_left_ref_est = pair_lookup%atom_pairs(ipair)
            IF (ipair_left_ref_est > 0) THEN
               left_second_atom = ri_data%lri_kp%pair_blocks(ipair_left_ref_est)%jatom
            ELSE
               left_second_atom = ri_data%lri_kp%pair_blocks(-ipair_left_ref_est)%iatom
            END IF
            target_left_est = target_left_est + 1_int_8
            scan_candidates_est = scan_candidates_est + INT( &
                                  MIN(density_index%atom_count(left_second_atom), &
                                      pair_lookup%atom_count(target_pairs(ipair_target)%jatom)), int_8)
         END DO
      END DO
      IF (ri_data%unit_nr > 0) THEN
         WRITE (ri_data%unit_nr, "(/,T2,A)") "LRI-HFXk K contraction"
         WRITE (ri_data%unit_nr, "(T4,A,T28,I14)") "target-left groups", target_left_est
         WRITE (ri_data%unit_nr, "(T4,A,T28,I14)") "sparse-join probes", scan_candidates_est
         FLUSH (ri_data%unit_nr)
      END IF

      ! A Coulomb-integral engine owns mutable scratch state, so every OpenMP thread gets a
      ! private context. Each iteration writes one distinct K-block slot and therefore needs
      ! neither a lock nor a thread-sized replica of the complete result.
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP          SHARED(basis_set_RI, density_cache, density_index, maxl_ri, pair_lookup, &
      !$OMP                 ri_data, ri_transform_cache, target_pairs) &
      !$OMP          PRIVATE(candidate_pairs_tmp, exchange_workspace, integral_context, &
      !$OMP                  ipair_target, ithread, k_tmp_block, k_tmp_norm, &
      !$OMP                  rho_tmp_norm, skipped_cell_bounds_tmp, skipped_image_tmp, &
      !$OMP                  skipped_img_index_tmp, skipped_missing_rho_tmp, &
      !$OMP                  skipped_no_kblock_tmp, skipped_rho_dim_tmp, skipped_small_rho_tmp, &
      !$OMP                  used_pairs_tmp) &
      !$OMP          REDUCTION(+:k_blocks_norm, multi_candidate_pairs, multi_skipped_cell_bounds, &
      !$OMP                      multi_skipped_image, multi_skipped_img_index, &
      !$OMP                      multi_skipped_missing_rho, multi_skipped_no_kblock, &
      !$OMP                      multi_skipped_rho_dim, multi_skipped_small_rho, multi_used_pairs, &
      !$OMP                      filtered_k_blocks, nk_blocks_built, rho_blocks_norm, metric_cache_hits, &
      !$OMP                      metric_cache_misses, metric_cache_screened, pair_screened)
      ithread = 1
!$    ithread = omp_get_thread_num() + 1
      CALL hfx_lri_exchange_workspace_move(ri_data%lri_kp%thread_workspaces(ithread), &
                                           exchange_workspace)
      exchange_workspace%pair_screened = 0_int_8
      exchange_workspace%metric_cache%hits = 0_int_8
      exchange_workspace%metric_cache%misses = 0_int_8
      exchange_workspace%metric_cache%screened = 0_int_8
      CALL coulomb_integral_init(integral_context, max_am_2c=maxl_ri)

      !$OMP DO SCHEDULE(DYNAMIC, 1)
      DO ipair_target = 1, SIZE(target_pairs)
         CALL hfx_lri_accumulate_target_k_block(basis_set_RI, pair_lookup, density_cache, density_index, &
                                                ri_transform_cache, integral_context, exchange_workspace, &
                                                target_pairs(ipair_target), ri_data, &
                                                k_tmp_block, used_pairs_tmp, &
                                                candidate_pairs_tmp, rho_tmp_norm, skipped_image_tmp, &
                                                skipped_cell_bounds_tmp, skipped_img_index_tmp, &
                                                skipped_missing_rho_tmp, skipped_rho_dim_tmp, &
                                                skipped_small_rho_tmp, skipped_no_kblock_tmp)
         multi_used_pairs = multi_used_pairs + used_pairs_tmp
         multi_candidate_pairs = multi_candidate_pairs + candidate_pairs_tmp
         multi_skipped_image = multi_skipped_image + skipped_image_tmp
         multi_skipped_cell_bounds = multi_skipped_cell_bounds + skipped_cell_bounds_tmp
         multi_skipped_img_index = multi_skipped_img_index + skipped_img_index_tmp
         multi_skipped_missing_rho = multi_skipped_missing_rho + skipped_missing_rho_tmp
         multi_skipped_rho_dim = multi_skipped_rho_dim + skipped_rho_dim_tmp
         multi_skipped_small_rho = multi_skipped_small_rho + skipped_small_rho_tmp
         multi_skipped_no_kblock = multi_skipped_no_kblock + skipped_no_kblock_tmp
         IF (.NOT. ALLOCATED(k_tmp_block)) CYCLE
         IF (used_pairs_tmp <= 0) THEN
            DEALLOCATE (k_tmp_block)
            CYCLE
         END IF
         k_tmp_norm = SQRT(SUM(k_tmp_block*k_tmp_block))
         IF (k_tmp_norm < ri_data%filter_eps) THEN
            DEALLOCATE (k_tmp_block)
            filtered_k_blocks = filtered_k_blocks + 1
            CYCLE
         END IF
         CALL hfx_lri_store_k_block( &
            ri_data%lri_kp%k_blocks(ipair_target), &
            target_pairs(ipair_target), ipair_target, k_tmp_block)
         nk_blocks_built = nk_blocks_built + 1
         rho_blocks_norm = rho_blocks_norm + rho_tmp_norm*rho_tmp_norm
         k_blocks_norm = k_blocks_norm + k_tmp_norm*k_tmp_norm
      END DO
      !$OMP END DO

      metric_cache_hits = metric_cache_hits + exchange_workspace%metric_cache%hits
      metric_cache_misses = metric_cache_misses + exchange_workspace%metric_cache%misses
      metric_cache_screened = metric_cache_screened + exchange_workspace%metric_cache%screened
      pair_screened = pair_screened + exchange_workspace%pair_screened
      CALL coulomb_integral_cleanup(integral_context)
      CALL hfx_lri_exchange_workspace_move(exchange_workspace, &
                                           ri_data%lri_kp%thread_workspaces(ithread))
      !$OMP END PARALLEL

      ri_data%lri_kp%nk_blocks = nk_blocks_built
      rho_blocks_norm = SQRT(rho_blocks_norm)
      k_blocks_norm = SQRT(k_blocks_norm)

      ! Form the full real-space trace from the symmetric sab_nl representatives. Off-diagonal
      ! atom blocks have one stored orientation and therefore account for their transpose partner;
      ! diagonal atom blocks contain every lattice image explicitly.
      exchange_trace_raw = 0.0_dp
      DO ipair = 1, SIZE(ri_data%lri_kp%k_blocks)
         IF (.NOT. ALLOCATED(ri_data%lri_kp%k_blocks(ipair)%k_block)) CYCLE
         rho_energy_block = hfx_lri_density_cache_find( &
                            density_cache, ri_data%lri_kp%k_blocks(ipair)%iatom, &
                            ri_data%lri_kp%k_blocks(ipair)%jatom, &
                            ri_data%lri_kp%k_blocks(ipair)%cell)
         IF (rho_energy_block == 0) CYCLE
         IF (.NOT. ALLOCATED(density_cache%blocks(rho_energy_block)%values)) CYCLE
         IF (ANY(SHAPE(density_cache%blocks(rho_energy_block)%values) /= &
                 SHAPE(ri_data%lri_kp%k_blocks(ipair)%k_block))) THEN
            CPABORT("Inconsistent density/K dimensions in LRI-HFXk exchange trace")
         END IF
         trace_weight = MERGE(1.0_dp, 2.0_dp, &
                              ri_data%lri_kp%k_blocks(ipair)%iatom == &
                              ri_data%lri_kp%k_blocks(ipair)%jatom)
         exchange_trace_raw = exchange_trace_raw + trace_weight* &
                              SUM(density_cache%blocks(rho_energy_block)%values* &
                                  ri_data%lri_kp%k_blocks(ipair)%k_block)
      END DO
      ri_data%lri_kp%exchange_trace_raw = exchange_trace_raw

      nearest_targets = 0
      nearest_candidates = 0
      nearest_tests = 0
      nearest_rho_active = 0
      nearest_k_active_low = 0
      nearest_k_active_mid = 0
      nearest_k_active_high = 0
      nearest_k_active_large = 0
      nearest_k_rms = 0.0_dp
      nearest_k_max = 0.0_dp
      general_targets = 0
      general_density_pairs = 0
      general_candidates = 0
      general_tests = 0
      general_rho_active = 0
      general_k_active_low = 0
      general_k_active_mid = 0
      general_k_active_high = 0
      general_k_active_large = 0
      general_k_rms = 0.0_dp
      general_k_max = 0.0_dp

      IF (lri_debug_level >= 1) THEN
         CALL hfx_lri_debug_nearest_origin_contributions(basis_set_RI, rho_ao, cell_to_index, sab_nl, &
                                                         ri_data, nearest_targets, nearest_candidates, &
                                                         nearest_tests, nearest_rho_active, &
                                                         nearest_k_active_low, nearest_k_active_mid, &
                                                         nearest_k_active_high, nearest_k_active_large, &
                                                         nearest_k_rms, nearest_k_max)

         CALL hfx_lri_debug_general_pair_contributions(basis_set_RI, rho_ao, cell_to_index, ri_data, &
                                                       general_targets, general_density_pairs, &
                                                       general_candidates, general_tests, &
                                                       general_rho_active, general_k_active_low, &
                                                       general_k_active_mid, general_k_active_high, &
                                                       general_k_active_large, general_k_rms, &
                                                       general_k_max)
      END IF

      IF (ri_data%unit_nr > 0) THEN
         IF (lri_debug_level >= 1) THEN
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
            WRITE (ri_data%unit_nr, "(/,T4,A)") "LRI-HFXk diagnostic K accumulation"
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
               "used/candidates", used_pairs, candidate_pairs
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
               "skip cell/img", skipped_cell_bounds, skipped_img_index
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
               "skip rho/dim", skipped_missing_rho, skipped_rho_dim
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
               "skip small/noK", skipped_small_rho, skipped_no_kblock
            WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3)") &
               "rho threshold", lri_debug_rho_threshold
            WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3,T41,ES10.3)") &
               "rho/K norm", rho_accum_norm, k_accum_norm
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
               "stored K dims", SIZE(ri_data%lri_kp%diagnostic_k_block%k_block, 1), &
               SIZE(ri_data%lri_kp%diagnostic_k_block%k_block, 2)
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10)") &
               "stored K pair", ri_data%lri_kp%diagnostic_k_block%pair_index
         END IF
         WRITE (ri_data%unit_nr, "(/,T4,A)") "LRI-HFXk rank-local K construction"
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
            "stored/max blocks", ri_data%lri_kp%nk_blocks, max_k_blocks
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10)") &
            "filtered K blocks", filtered_k_blocks
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
            "used/candidates", multi_used_pairs, multi_candidate_pairs
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
            "skip cell/img", multi_skipped_cell_bounds, multi_skipped_img_index
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
            "skip rho/dim", multi_skipped_missing_rho, multi_skipped_rho_dim
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
            "skip small/noK", multi_skipped_small_rho, multi_skipped_no_kblock
         WRITE (ri_data%unit_nr, "(T6,A,T28,I14,T45,I14)") &
            "metric hit/miss", metric_cache_hits, metric_cache_misses
         WRITE (ri_data%unit_nr, "(T6,A,T28,I14)") &
            "metric dist-screen", metric_cache_screened
         WRITE (ri_data%unit_nr, "(T6,A,T28,I14)") &
            "pair dist-screen", pair_screened
         WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3,T41,ES10.3)") &
            "rho/K norm", rho_blocks_norm, k_blocks_norm
         WRITE (ri_data%unit_nr, "(T6,A,T28,ES14.6)") &
            "raw Tr[D K]", exchange_trace_raw
         IF (lri_debug_level >= 1) THEN
            WRITE (ri_data%unit_nr, "(/,T4,A)") "LRI-HFXk nearest-origin contribution probe"
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
               "targets/tests", nearest_targets, nearest_tests
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
               "candidates/rho", nearest_candidates, nearest_rho_active
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,ES10.3)") &
               "K > threshold", nearest_k_active_low, lri_debug_k_threshold_low
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,ES10.3)") &
               "K > threshold", nearest_k_active_mid, lri_debug_k_threshold_mid
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,ES10.3)") &
               "K > threshold", nearest_k_active_high, lri_debug_k_threshold_high
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,ES10.3)") &
               "K > threshold", nearest_k_active_large, lri_debug_k_threshold_large
            WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3,T41,ES10.3)") &
               "K rms/max", nearest_k_rms, nearest_k_max
            WRITE (ri_data%unit_nr, "(/,T4,A)") "LRI-HFXk general pair-pair contribution probe"
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
               "target/density", general_targets, general_density_pairs
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
               "candidates/tests", general_candidates, general_tests
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10)") &
               "rho active", general_rho_active
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,ES10.3)") &
               "K > threshold", general_k_active_low, lri_debug_k_threshold_low
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,ES10.3)") &
               "K > threshold", general_k_active_mid, lri_debug_k_threshold_mid
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,ES10.3)") &
               "K > threshold", general_k_active_high, lri_debug_k_threshold_high
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,ES10.3)") &
               "K > threshold", general_k_active_large, lri_debug_k_threshold_large
            WRITE (ri_data%unit_nr, "(T6,A,T28,ES10.3,T41,ES10.3)") &
               "K rms/max", general_k_rms, general_k_max
         END IF
      END IF

      IF (ALLOCATED(k_accum_block)) DEALLOCATE (k_accum_block)
      IF (ALLOCATED(target_pairs)) DEALLOCATE (target_pairs)
      IF (ASSOCIATED(dft_control)) THEN
         DO ikind = 1, nkind
            ri_basis => basis_set_RI(ikind)%gto_basis_set
            IF (ASSOCIATED(ri_basis)) THEN
               CALL init_interaction_radii_orb_basis(ri_basis, dft_control%qs_control%eps_pgf_orb)
            END IF
         END DO
      END IF
      IF (ASSOCIATED(basis_set_AO)) DEALLOCATE (basis_set_AO)
      IF (ASSOCIATED(basis_set_RI)) DEALLOCATE (basis_set_RI)

      CALL timestop(handle)

   END SUBROUTINE hfx_lri_build_rank_local_k_blocks

! **************************************************************************************************
!> \brief Probes nearest-origin target K norms without changing the current K path.
!> \param basis_set_RI ...
!> \param rho_ao ...
!> \param cell_to_index ...
!> \param sab_nl ...
!> \param ri_data ...
!> \param n_target_pairs ...
!> \param n_candidates ...
!> \param n_tests ...
!> \param n_rho_active ...
!> \param n_k_active_low ...
!> \param n_k_active_mid ...
!> \param n_k_active_high ...
!> \param n_k_active_large ...
!> \param k_norm_rms ...
!> \param k_norm_max ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_debug_nearest_origin_contributions(basis_set_RI, rho_ao, cell_to_index, sab_nl, &
                                                         ri_data, n_target_pairs, n_candidates, &
                                                         n_tests, n_rho_active, n_k_active_low, &
                                                         n_k_active_mid, n_k_active_high, &
                                                         n_k_active_large, k_norm_rms, k_norm_max)

      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      INTEGER, INTENT(OUT) :: n_target_pairs, n_candidates, n_tests, n_rho_active, n_k_active_low, &
         n_k_active_mid, n_k_active_high, n_k_active_large
      REAL(dp), INTENT(OUT)                              :: k_norm_rms, k_norm_max

      INTEGER :: candidate_pairs, ikind, ipair_target, maxl_ri, maxl_ri_kind, natom, &
         skipped_cell_bounds, skipped_image, skipped_img_index, skipped_missing_rho, &
         skipped_no_kblock, skipped_rho_dim, skipped_small_rho, used_pairs
      REAL(dp)                                           :: k_norm, k_norm2, rho_accum_norm
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: k_block
      TYPE(coulomb_integral_context_type)                :: integral_context
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(hfx_lri_exchange_workspace_type)              :: exchange_workspace
      TYPE(lri_basis_transform_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: ri_transform_cache
      TYPE(lri_density_cache_type)                       :: density_cache
      TYPE(lri_density_index_type)                       :: density_index
      TYPE(lri_pair_lookup_type)                         :: pair_lookup

      NULLIFY (ri_basis)

      n_target_pairs = 0
      n_candidates = 0
      n_tests = 0
      n_rho_active = 0
      n_k_active_low = 0
      n_k_active_mid = 0
      n_k_active_high = 0
      n_k_active_large = 0
      k_norm_rms = 0.0_dp
      k_norm_max = 0.0_dp
      k_norm2 = 0.0_dp

      IF (.NOT. ASSOCIATED(ri_data%lri_kp)) RETURN
      IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks)) RETURN
      IF (ri_data%lri_kp%npair_blocks <= 0) RETURN
      IF (.NOT. ASSOCIATED(rho_ao)) RETURN

      natom = 0
      DO ipair_target = 1, ri_data%lri_kp%npair_blocks
         natom = MAX(natom, ri_data%lri_kp%pair_blocks(ipair_target)%iatom, &
                     ri_data%lri_kp%pair_blocks(ipair_target)%jatom)
      END DO
      CALL hfx_lri_pair_lookup_create(pair_lookup, ri_data%lri_kp%pair_blocks, &
                                      ri_data%lri_kp%npair_blocks, natom)
      CALL hfx_lri_density_cache_create(rho_ao, cell_to_index, sab_nl, ri_data%filter_eps, density_cache)
      CALL hfx_lri_density_index_create(density_cache, natom, density_index)
      CALL hfx_lri_build_ri_transform_cache(basis_set_RI, ri_transform_cache)
      maxl_ri = 0
      DO ikind = 1, SIZE(basis_set_RI)
         ri_basis => basis_set_RI(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(ri_basis)) CYCLE
         CALL get_gto_basis_set(ri_basis, maxl=maxl_ri_kind)
         maxl_ri = MAX(maxl_ri, maxl_ri_kind)
      END DO
      CALL init_md_ftable(nmax=2*maxl_ri)
      CALL coulomb_integral_init(integral_context, max_am_2c=maxl_ri)

      DO ipair_target = 1, ri_data%lri_kp%npair_blocks
         IF (n_target_pairs >= lri_debug_max_general_targets) EXIT
         IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks(ipair_target)%c_fit)) CYCLE
         n_target_pairs = n_target_pairs + 1

         CALL hfx_lri_accumulate_target_k_block(basis_set_RI, pair_lookup, density_cache, density_index, &
                                                ri_transform_cache, integral_context, exchange_workspace, &
                                                ri_data%lri_kp%pair_blocks(ipair_target), ri_data, &
                                                k_block, used_pairs, candidate_pairs, &
                                                rho_accum_norm, skipped_image, skipped_cell_bounds, &
                                                skipped_img_index, skipped_missing_rho, &
                                                skipped_rho_dim, skipped_small_rho, skipped_no_kblock)
         n_candidates = n_candidates + candidate_pairs
         n_rho_active = n_rho_active + used_pairs
         IF (.NOT. ALLOCATED(k_block)) CYCLE
         n_tests = n_tests + 1
         k_norm = SQRT(SUM(k_block*k_block))
         k_norm2 = k_norm2 + k_norm*k_norm
         k_norm_max = MAX(k_norm_max, k_norm)
         IF (k_norm > lri_debug_k_threshold_low) n_k_active_low = n_k_active_low + 1
         IF (k_norm > lri_debug_k_threshold_mid) n_k_active_mid = n_k_active_mid + 1
         IF (k_norm > lri_debug_k_threshold_high) n_k_active_high = n_k_active_high + 1
         IF (k_norm > lri_debug_k_threshold_large) n_k_active_large = n_k_active_large + 1
         DEALLOCATE (k_block)
      END DO

      IF (n_tests > 0) THEN
         k_norm_rms = SQRT(k_norm2/REAL(n_tests, KIND=dp))
      END IF
      CALL coulomb_integral_cleanup(integral_context)

   END SUBROUTINE hfx_lri_debug_nearest_origin_contributions

! **************************************************************************************************
!> \brief Probes general pair-pair contribution norms without changing the current K path.
!> \param basis_set_RI ...
!> \param rho_ao ...
!> \param cell_to_index ...
!> \param ri_data ...
!> \param n_target_pairs ...
!> \param n_density_pairs ...
!> \param n_candidates ...
!> \param n_tests ...
!> \param n_rho_active ...
!> \param n_k_active_low ...
!> \param n_k_active_mid ...
!> \param n_k_active_high ...
!> \param n_k_active_large ...
!> \param k_norm_rms ...
!> \param k_norm_max ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_debug_general_pair_contributions(basis_set_RI, rho_ao, cell_to_index, ri_data, &
                                                       n_target_pairs, n_density_pairs, n_candidates, &
                                                       n_tests, n_rho_active, n_k_active_low, &
                                                       n_k_active_mid, n_k_active_high, &
                                                       n_k_active_large, k_norm_rms, k_norm_max)

      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      INTEGER, INTENT(OUT) :: n_target_pairs, n_density_pairs, n_candidates, n_tests, &
         n_rho_active, n_k_active_low, n_k_active_mid, n_k_active_high, n_k_active_large
      REAL(dp), INTENT(OUT)                              :: k_norm_rms, k_norm_max

      INTEGER                                            :: candidate_index, img_pair, img_rev, &
                                                            ipair_density, ipair_target, lb1, lb2, &
                                                            lb3, sample_stride, ub1, ub2, ub3
      INTEGER, DIMENSION(3)                              :: cell_rev
      LOGICAL                                            :: found, found_rev
      REAL(dp)                                           :: k_norm, k_norm2, rho_norm
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: k_block, rho_work
      REAL(dp), DIMENSION(:, :), POINTER                 :: rho_block, rho_block_rev

      NULLIFY (rho_block, rho_block_rev)

      n_target_pairs = 0
      n_density_pairs = 0
      n_candidates = 0
      n_tests = 0
      n_rho_active = 0
      n_k_active_low = 0
      n_k_active_mid = 0
      n_k_active_high = 0
      n_k_active_large = 0
      k_norm_rms = 0.0_dp
      k_norm_max = 0.0_dp
      k_norm2 = 0.0_dp

      IF (.NOT. ASSOCIATED(ri_data%lri_kp)) RETURN
      IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks)) RETURN
      IF (ri_data%lri_kp%npair_blocks <= 0) RETURN
      IF (.NOT. ASSOCIATED(rho_ao)) RETURN

      DO ipair_density = 1, ri_data%lri_kp%npair_blocks
         IF (ALLOCATED(ri_data%lri_kp%pair_blocks(ipair_density)%c_fit)) THEN
            n_density_pairs = n_density_pairs + 1
         END IF
      END DO

      DO ipair_target = 1, ri_data%lri_kp%npair_blocks
         IF (n_target_pairs >= lri_debug_max_general_targets) EXIT
         IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks(ipair_target)%c_fit)) CYCLE
         n_target_pairs = n_target_pairs + 1
      END DO

      n_candidates = n_target_pairs*n_density_pairs
      sample_stride = MAX(1, n_candidates/MAX(1, lri_debug_max_general_tests))
      candidate_index = 0
      n_target_pairs = 0

      DO ipair_target = 1, ri_data%lri_kp%npair_blocks
         IF (n_target_pairs >= lri_debug_max_general_targets) EXIT
         IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks(ipair_target)%c_fit)) CYCLE
         n_target_pairs = n_target_pairs + 1
         DO ipair_density = 1, ri_data%lri_kp%npair_blocks
            IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks(ipair_density)%c_fit)) CYCLE
            candidate_index = candidate_index + 1
            IF (MOD(candidate_index - 1, sample_stride) /= 0) CYCLE
            IF (n_tests >= lri_debug_max_general_tests) CYCLE
            IF (ri_data%lri_kp%pair_blocks(ipair_density)%iatom /= &
                ri_data%lri_kp%pair_blocks(ipair_target)%jatom) CYCLE
            n_tests = n_tests + 1

            img_pair = 1
            IF (ASSOCIATED(cell_to_index)) THEN
               lb1 = LBOUND(cell_to_index, 1); ub1 = UBOUND(cell_to_index, 1)
               lb2 = LBOUND(cell_to_index, 2); ub2 = UBOUND(cell_to_index, 2)
               lb3 = LBOUND(cell_to_index, 3); ub3 = UBOUND(cell_to_index, 3)
               IF (.NOT. (ALL(ri_data%lri_kp%pair_blocks(ipair_density)%cell >= [lb1, lb2, lb3]) .AND. &
                          ALL(ri_data%lri_kp%pair_blocks(ipair_density)%cell <= [ub1, ub2, ub3]))) CYCLE
               img_pair = cell_to_index(ri_data%lri_kp%pair_blocks(ipair_density)%cell(1), &
                                        ri_data%lri_kp%pair_blocks(ipair_density)%cell(2), &
                                        ri_data%lri_kp%pair_blocks(ipair_density)%cell(3))
            END IF
            IF (img_pair < 1 .OR. img_pair > SIZE(rho_ao, 2)) CYCLE

            CALL dbcsr_get_block_p(matrix=rho_ao(1, img_pair)%matrix, &
                                   row=ri_data%lri_kp%pair_blocks(ipair_density)%iatom, &
                                   col=ri_data%lri_kp%pair_blocks(ipair_density)%jatom, &
                                   block=rho_block, found=found)
            IF (ALLOCATED(rho_work)) DEALLOCATE (rho_work)
            IF (found) THEN
               IF (SIZE(rho_block, 1) < ri_data%lri_kp%pair_blocks(ipair_density)%nao_i .OR. &
                   SIZE(rho_block, 2) < ri_data%lri_kp%pair_blocks(ipair_density)%nao_j) CYCLE
               ALLOCATE (rho_work(ri_data%lri_kp%pair_blocks(ipair_density)%nao_i, &
                                  ri_data%lri_kp%pair_blocks(ipair_density)%nao_j))
               rho_work(:, :) = rho_block(1:ri_data%lri_kp%pair_blocks(ipair_density)%nao_i, &
                                          1:ri_data%lri_kp%pair_blocks(ipair_density)%nao_j)
            ELSE
               found_rev = .FALSE.
               cell_rev(:) = -ri_data%lri_kp%pair_blocks(ipair_density)%cell(:)
               img_rev = 0
               IF (ASSOCIATED(cell_to_index)) THEN
                  IF (ALL(cell_rev >= [lb1, lb2, lb3]) .AND. ALL(cell_rev <= [ub1, ub2, ub3])) THEN
                     img_rev = cell_to_index(cell_rev(1), cell_rev(2), cell_rev(3))
                  END IF
               ELSE
                  img_rev = 1
               END IF
               IF (img_rev >= 1 .AND. img_rev <= SIZE(rho_ao, 2)) THEN
                  CALL dbcsr_get_block_p(matrix=rho_ao(1, img_rev)%matrix, &
                                         row=ri_data%lri_kp%pair_blocks(ipair_density)%jatom, &
                                         col=ri_data%lri_kp%pair_blocks(ipair_density)%iatom, &
                                         block=rho_block_rev, found=found_rev)
               END IF
               IF (.NOT. found_rev) CYCLE
               IF (SIZE(rho_block_rev, 1) < ri_data%lri_kp%pair_blocks(ipair_density)%nao_j .OR. &
                   SIZE(rho_block_rev, 2) < ri_data%lri_kp%pair_blocks(ipair_density)%nao_i) CYCLE
               ALLOCATE (rho_work(ri_data%lri_kp%pair_blocks(ipair_density)%nao_i, &
                                  ri_data%lri_kp%pair_blocks(ipair_density)%nao_j))
               rho_work(:, :) = TRANSPOSE(rho_block_rev(1:ri_data%lri_kp%pair_blocks(ipair_density)%nao_j, &
                                                        1:ri_data%lri_kp%pair_blocks(ipair_density)%nao_i))
            END IF

            rho_norm = SQRT(SUM(rho_work*rho_work))
            IF (rho_norm < lri_debug_rho_threshold) CYCLE
            n_rho_active = n_rho_active + 1

            CALL hfx_lri_contract_pair_pair_density_realv(basis_set_RI, &
                                                          ri_data%lri_kp%pair_blocks(ipair_target), &
                                                          ri_data%lri_kp%pair_blocks(ipair_target), &
                                                          ri_data%lri_kp%pair_blocks(ipair_density), &
                                                          ri_data, rho_work, k_block)
            IF (.NOT. ALLOCATED(k_block)) CYCLE
            k_norm = SQRT(SUM(k_block*k_block))
            k_norm2 = k_norm2 + k_norm*k_norm
            k_norm_max = MAX(k_norm_max, k_norm)
            IF (k_norm > lri_debug_k_threshold_low) n_k_active_low = n_k_active_low + 1
            IF (k_norm > lri_debug_k_threshold_mid) n_k_active_mid = n_k_active_mid + 1
            IF (k_norm > lri_debug_k_threshold_high) n_k_active_high = n_k_active_high + 1
            IF (k_norm > lri_debug_k_threshold_large) n_k_active_large = n_k_active_large + 1
            DEALLOCATE (k_block)
         END DO
      END DO

      IF (ALLOCATED(rho_work)) DEALLOCATE (rho_work)
      IF (n_tests > 0) THEN
         k_norm_rms = SQRT(k_norm2/REAL(n_tests, KIND=dp))
      END IF

   END SUBROUTINE hfx_lri_debug_general_pair_contributions

! **************************************************************************************************
!> \brief Adds the completed pair-local exchange blocks to the k-point KS matrices.
!> \param qs_env ...
!> \param ri_data ...
!> \param ks_matrix ...
!> \param hf_fraction ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_write_rank_local_k_blocks(qs_env, ri_data, ks_matrix, hf_fraction)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      REAL(dp), INTENT(IN)                               :: hf_fraction

      INTEGER                                            :: added_blocks, icol, img, ipair, irow, &
                                                            transposed_blocks
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found, use_transpose
      REAL(dp)                                           :: scale
      REAL(dp), DIMENSION(:, :), POINTER                 :: ks_block
      TYPE(kpoint_type), POINTER                         :: kpoints

      NULLIFY (cell_to_index, kpoints, ks_block)

      IF (.NOT. ASSOCIATED(ri_data%lri_kp)) THEN
         CPABORT("LRI-HFXk state is not associated during KS writeback")
      END IF
      IF (.NOT. ASSOCIATED(ks_matrix)) THEN
         CPABORT("LRI-HFXk KS matrix is not associated")
      END IF

      CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)
      IF (.NOT. ASSOCIATED(kpoints)) CPABORT("LRI-HFXk could not obtain k-point metadata")
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
      IF (.NOT. ASSOCIATED(cell_to_index)) CPABORT("LRI-HFXk could not obtain the image map")

      ! The input density is spin-summed RKS, so its same-spin density is rho/2.
      scale = -0.5_dp*hf_fraction

      added_blocks = 0
      transposed_blocks = 0
      IF (ALLOCATED(ri_data%lri_kp%k_blocks)) THEN
         DO ipair = 1, SIZE(ri_data%lri_kp%k_blocks)
            IF (.NOT. ALLOCATED(ri_data%lri_kp%k_blocks(ipair)%k_block)) CYCLE

            IF (ANY(ri_data%lri_kp%k_blocks(ipair)%cell < LBOUND(cell_to_index)) .OR. &
                ANY(ri_data%lri_kp%k_blocks(ipair)%cell > UBOUND(cell_to_index))) THEN
               CPABORT("LRI-HFXk target cell is outside the k-point image map")
            END IF
            img = cell_to_index(ri_data%lri_kp%k_blocks(ipair)%cell(1), &
                                ri_data%lri_kp%k_blocks(ipair)%cell(2), &
                                ri_data%lri_kp%k_blocks(ipair)%cell(3))
            IF (img < 1 .OR. img > SIZE(ks_matrix, 2)) THEN
               CPABORT("LRI-HFXk target has no corresponding KS image")
            END IF

            irow = MIN(ri_data%lri_kp%k_blocks(ipair)%iatom, &
                       ri_data%lri_kp%k_blocks(ipair)%jatom)
            icol = MAX(ri_data%lri_kp%k_blocks(ipair)%iatom, &
                       ri_data%lri_kp%k_blocks(ipair)%jatom)
            use_transpose = ri_data%lri_kp%k_blocks(ipair)%iatom > &
                            ri_data%lri_kp%k_blocks(ipair)%jatom
            NULLIFY (ks_block)
            CALL dbcsr_get_block_p(matrix=ks_matrix(1, img)%matrix, row=irow, col=icol, &
                                   block=ks_block, found=found)
            IF (.NOT. found) CPABORT("LRI-HFXk target block is absent from the KS matrix")

            IF (use_transpose) THEN
               IF (ANY(SHAPE(ks_block) /= &
                       [ri_data%lri_kp%k_blocks(ipair)%nao_j, &
                        ri_data%lri_kp%k_blocks(ipair)%nao_i])) THEN
                  CPABORT("Inconsistent transposed LRI-HFXk KS block dimensions")
               END IF
               ks_block(:, :) = ks_block(:, :) + scale*TRANSPOSE(ri_data%lri_kp%k_blocks(ipair)%k_block(:, :))
               transposed_blocks = transposed_blocks + 1
            ELSE
               IF (ANY(SHAPE(ks_block) /= SHAPE(ri_data%lri_kp%k_blocks(ipair)%k_block))) THEN
                  CPABORT("Inconsistent LRI-HFXk KS block dimensions")
               END IF
               ks_block(:, :) = ks_block(:, :) + scale*ri_data%lri_kp%k_blocks(ipair)%k_block(:, :)
            END IF
            added_blocks = added_blocks + 1
         END DO
      END IF

      IF (added_blocks /= ri_data%lri_kp%nk_blocks) THEN
         CPABORT("Not all LRI-HFXk K blocks were written to the KS matrix")
      END IF

      IF (ri_data%unit_nr > 0) THEN
         WRITE (ri_data%unit_nr, "(/,T4,A)") "LRI-HFXk rank-local KS writeback"
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
            "added/stored", added_blocks, ri_data%lri_kp%nk_blocks
         WRITE (ri_data%unit_nr, "(T6,A,T28,I10)") &
            "transposed", transposed_blocks
      END IF

   END SUBROUTINE hfx_lri_write_rank_local_k_blocks

! **************************************************************************************************
!> \brief Computes the RKS HFX energy from the completed pair-local exchange trace.
!> \param qs_env ...
!> \param ri_data ...
!> \param hf_fraction ...
!> \param ehfx ...
!> \note For a spin-summed RKS density, E_x = -hf_fraction*Tr[D K_raw]/4.
! **************************************************************************************************
   SUBROUTINE hfx_lri_compute_rank_local_energy(qs_env, ri_data, hf_fraction, ehfx)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      REAL(dp), INTENT(IN)                               :: hf_fraction
      REAL(dp), INTENT(OUT)                              :: ehfx

      REAL(dp)                                           :: exchange_trace
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (para_env)

      ehfx = 0.0_dp
      IF (.NOT. ASSOCIATED(ri_data%lri_kp)) RETURN
      CALL get_qs_env(qs_env=qs_env, para_env=para_env)
      exchange_trace = ri_data%lri_kp%exchange_trace_raw
      IF (ASSOCIATED(para_env)) CALL para_env%sum(exchange_trace)
      ! For an RKS spin-summed density, F_x = -hf_fraction*K_raw/2 and
      ! E_x = Tr[D F_x]/2 = -hf_fraction*Tr[D K_raw]/4.
      ehfx = -0.25_dp*hf_fraction*exchange_trace

      IF (ri_data%unit_nr > 0) THEN
         WRITE (ri_data%unit_nr, "(/,T4,A)") "LRI-HFXk exchange energy"
         WRITE (ri_data%unit_nr, "(T6,A,T28,ES14.6)") &
            "Tr[D K_raw]", exchange_trace
         WRITE (ri_data%unit_nr, "(T6,A,T28,ES14.6)") &
            "HF fraction", hf_fraction
         WRITE (ri_data%unit_nr, "(T6,A,T28,ES14.6)") &
            "E scaled", ehfx
      END IF

   END SUBROUTINE hfx_lri_compute_rank_local_energy

! **************************************************************************************************
!> \brief Adds one pair-pair-density contribution directly to an existing K block.
!> \param basis_set_RI ...
!> \param pair_output ...
!> \param pair_left ...
!> \param left_reversed ...
!> \param pair_right ...
!> \param right_reversed ...
!> \param ri_data ...
!> \param rho_block ...
!> \param k_block ...
!> \param integral_context ...
!> \param ri_transform_cache ...
!> \param exchange_workspace ...
!> \param outer_initialized ...
!> \param outer_ready ...
!> \param outer_sum_used ...
!> \param outer_norm2 ...
!> \param contraction_done ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_accumulate_pair_pair_density_realv( &
      basis_set_RI, pair_output, pair_left, left_reversed, pair_right, right_reversed, &
      ri_data, rho_block, k_block, &
      integral_context, ri_transform_cache, exchange_workspace, &
      outer_initialized, outer_ready, outer_sum_used, outer_norm2, contraction_done)

      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair_output, pair_left
      LOGICAL, INTENT(IN)                                :: left_reversed
      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair_right
      LOGICAL, INTENT(IN)                                :: right_reversed
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: rho_block
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: k_block
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: integral_context
      TYPE(lri_basis_transform_type), DIMENSION(:), &
         INTENT(IN)                                      :: ri_transform_cache
      TYPE(hfx_lri_exchange_workspace_type), &
         INTENT(INOUT)                                   :: exchange_workspace
      LOGICAL, INTENT(INOUT)                             :: outer_initialized, outer_ready, &
                                                            outer_sum_used
      REAL(dp), INTENT(INOUT)                            :: outer_norm2
      LOGICAL, INTENT(OUT)                               :: contraction_done

      INTEGER :: inner_center, inner_offset, inner_size, left_iatom, left_nao_i, left_nao_j, &
         outer_center, outer_kind, outer_offset, outer_size, right_iatom, right_nao_i, right_nao_j
      LOGICAL                                            :: inner_done, outer_done
      REAL(dp)                                           :: inner_norm2, v_norm2
      REAL(dp), DIMENSION(3)                             :: right_inner_position

      contraction_done = .FALSE.
      IF (.NOT. ALLOCATED(pair_left%c_fit) .OR. .NOT. ALLOCATED(pair_right%c_fit)) RETURN
      left_iatom = MERGE(pair_left%jatom, pair_left%iatom, left_reversed)
      right_iatom = MERGE(pair_right%jatom, pair_right%iatom, right_reversed)
      left_nao_i = MERGE(pair_left%nao_j, pair_left%nao_i, left_reversed)
      left_nao_j = MERGE(pair_left%nao_i, pair_left%nao_j, left_reversed)
      right_nao_i = MERGE(pair_right%nao_j, pair_right%nao_i, right_reversed)
      right_nao_j = MERGE(pair_right%nao_i, pair_right%nao_j, right_reversed)
      IF (left_iatom /= pair_output%iatom .OR. right_iatom /= pair_output%jatom) THEN
         CPABORT("Inconsistent atom labels in LRI-HFXk exchange contraction")
      END IF
      IF (left_nao_i /= pair_output%nao_i .OR. right_nao_i /= pair_output%nao_j) THEN
         CPABORT("Inconsistent output dimensions in LRI-HFXk exchange contraction")
      END IF
      IF (SIZE(rho_block, 1) < left_nao_j .OR. SIZE(rho_block, 2) < right_nao_j) RETURN

      ! Reject a quartet before any metric-cache traffic or tensor packing when every one of its
      ! auxiliary-center blocks is outside the exact shell-distance cutoff used by eri_2center.
      ! This changes neither EPS_FILTER_2C nor the contracted result: the integral backend would
      ! return four identically zero center blocks for the same geometry.
      IF (.NOT. hfx_lri_pair_pair_may_couple( &
          pair_left, pair_right, left_reversed, right_reversed, pair_output%rab, &
          ri_data, ri_transform_cache)) THEN
         exchange_workspace%pair_screened = exchange_workspace%pair_screened + 1_int_8
         RETURN
      END IF

      ! In the logical K-L orientation, the first auxiliary center is fixed at the output atom K.
      ! Its metric coupling is independent of L and can therefore be precontracted and accumulated
      ! across every right pair belonging to this target/left combination.
      IF (pair_right%onsite_pair) THEN
         outer_center = 1
         inner_center = 0
      ELSE IF (right_reversed) THEN
         outer_center = 2
         inner_center = 1
      ELSE
         outer_center = 1
         inner_center = 2
      END IF
      outer_offset = pair_right%aux_offset(outer_center)
      outer_size = pair_right%aux_size(outer_center)
      outer_kind = pair_right%aux_kind(outer_center)
      IF (outer_size <= 0 .OR. outer_offset < 0 .OR. &
          outer_offset + outer_size > pair_right%naux_pair) THEN
         CPABORT("Invalid outer auxiliary slice in LRI-HFXk exchange contraction")
      END IF
      IF (outer_kind /= pair_output%jkind) THEN
         CPABORT("Outer auxiliary kind does not match the LRI-HFXk target atom")
      END IF

      IF (.NOT. outer_initialized) THEN
         CALL hfx_lri_build_pair_center_v_sph_context( &
            basis_set_RI, pair_left, left_reversed, outer_kind, pair_output%rab, outer_size, &
            ri_data, exchange_workspace%v_outer, integral_context, ri_transform_cache, &
            exchange_workspace%metric_cache, ri_data%lri_kp%metric_table, outer_norm2)
         outer_initialized = .TRUE.
      END IF
      IF (.NOT. ALLOCATED(exchange_workspace%v_outer)) THEN
         CPABORT("Missing outer-center metric in LRI-HFXk exchange contraction")
      END IF
      IF (ANY(SHAPE(exchange_workspace%v_outer) /= [pair_left%naux_pair, outer_size])) THEN
         CPABORT("Inconsistent outer-center metric dimensions in LRI-HFXk exchange contraction")
      END IF

      inner_norm2 = 0.0_dp
      inner_size = 0
      inner_offset = 0
      IF (inner_center > 0) THEN
         inner_size = pair_right%aux_size(inner_center)
         inner_offset = pair_right%aux_offset(inner_center)
         IF (inner_size <= 0 .OR. inner_offset < 0 .OR. &
             inner_offset + inner_size > pair_right%naux_pair) THEN
            CPABORT("Invalid inner auxiliary slice in LRI-HFXk exchange contraction")
         END IF
         right_inner_position(:) = pair_right%aux_r(:, inner_center) + pair_output%rab(:)
         IF (right_reversed) right_inner_position(:) = right_inner_position(:) - pair_right%rab(:)
         CALL hfx_lri_build_pair_center_v_sph_context( &
            basis_set_RI, pair_left, left_reversed, pair_right%aux_kind(inner_center), &
            right_inner_position, inner_size, ri_data, exchange_workspace%v_sph, &
            integral_context, ri_transform_cache, exchange_workspace%metric_cache, &
            ri_data%lri_kp%metric_table, inner_norm2)
      END IF

      ! The outer and inner column slices exactly partition the former full V block. Their squared
      ! norms therefore reproduce the old EPS_FILTER_2C decision without a split-block ambiguity.
      v_norm2 = outer_norm2 + inner_norm2
      IF (v_norm2 < ri_data%filter_eps_2c*ri_data%filter_eps_2c) RETURN

      outer_done = .FALSE.
      IF (outer_norm2 > 0.0_dp) THEN
         IF (.NOT. outer_ready) THEN
            CALL hfx_lri_contract_outer_prepare(pair_left%c_fit, exchange_workspace%v_outer, &
                                                exchange_workspace%contract, left_reversed)
            outer_ready = .TRUE.
         END IF
         CALL hfx_lri_contract_outer_accumulate( &
            pair_right%c_fit(:, :, outer_offset + 1:outer_offset + outer_size), &
            rho_block(1:left_nao_j, 1:right_nao_j), exchange_workspace%contract, &
            first=.NOT. outer_sum_used, transpose_right=right_reversed)
         outer_sum_used = .TRUE.
         outer_done = .TRUE.
      END IF

      inner_done = .FALSE.
      IF (inner_size > 0 .AND. inner_norm2 > 0.0_dp) THEN
         IF (.NOT. ALLOCATED(exchange_workspace%v_sph)) THEN
            CPABORT("Missing inner-center metric in LRI-HFXk exchange contraction")
         END IF
         IF (ANY(SHAPE(exchange_workspace%v_sph) /= [pair_left%naux_pair, inner_size])) THEN
            CPABORT("Inconsistent inner-center metric dimensions in LRI-HFXk exchange contraction")
         END IF
         CALL hfx_lri_contract_exchange_add( &
            pair_left%c_fit, &
            pair_right%c_fit(:, :, inner_offset + 1:inner_offset + inner_size), &
            exchange_workspace%v_sph, rho_block(1:left_nao_j, 1:right_nao_j), &
            k_block, exchange_workspace%contract, left_reversed, right_reversed)
         inner_done = .TRUE.
      END IF
      contraction_done = outer_done .OR. inner_done

   END SUBROUTINE hfx_lri_accumulate_pair_pair_density_realv

! **************************************************************************************************
!> \brief Contracts two fitting pairs and a density block into one exchange K block.
!> \param basis_set_RI ...
!> \param pair_output ...
!> \param pair_left ...
!> \param pair_right ...
!> \param ri_data ...
!> \param rho_block ...
!> \param k_block ...
!> \param integral_context ...
!> \param ri_transform_cache ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_contract_pair_pair_density_realv(basis_set_RI, pair_output, pair_left, &
                                                       pair_right, ri_data, rho_block, k_block, &
                                                       integral_context, ri_transform_cache)

      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair_output, pair_left, pair_right
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: rho_block
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: k_block
      TYPE(coulomb_integral_context_type), &
         INTENT(INOUT), OPTIONAL                         :: integral_context
      TYPE(lri_basis_transform_type), DIMENSION(:), &
         INTENT(IN), OPTIONAL                            :: ri_transform_cache

      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: v_sph

      IF (.NOT. ALLOCATED(pair_left%c_fit)) RETURN
      IF (.NOT. ALLOCATED(pair_right%c_fit)) RETURN
      IF (pair_left%iatom /= pair_output%iatom .OR. &
          pair_right%iatom /= pair_output%jatom) THEN
         CPABORT("Inconsistent atom labels in LRI-HFXk exchange contraction")
      END IF
      IF (pair_left%nao_i /= pair_output%nao_i .OR. &
          pair_right%nao_i /= pair_output%nao_j) THEN
         CPABORT("Inconsistent output dimensions in LRI-HFXk exchange contraction")
      END IF
      IF (SIZE(rho_block, 1) < pair_left%nao_j .OR. &
          SIZE(rho_block, 2) < pair_right%nao_j) RETURN

      IF (PRESENT(integral_context) .AND. PRESENT(ri_transform_cache)) THEN
         CALL hfx_lri_build_pair_pair_v_sph_context( &
            basis_set_RI, pair_left, pair_right, ri_data, v_sph, integral_context, &
            ri_transform_cache, &
            right_shift=pair_output%rab)
      ELSE
         CALL hfx_lri_build_pair_pair_v_sph( &
            basis_set_RI, pair_left, pair_right, ri_data, v_sph, &
            right_shift=pair_output%rab)
      END IF
      IF (.NOT. ALLOCATED(v_sph)) RETURN
      IF (SIZE(v_sph, 1) /= pair_left%naux_pair .OR. &
          SIZE(v_sph, 2) /= pair_right%naux_pair) THEN
         CPABORT("Inconsistent auxiliary dimension in LRI-HFXk pair-pair density contraction")
      END IF

      CALL hfx_lri_contract_exchange( &
         pair_left%c_fit, pair_right%c_fit, v_sph, &
         rho_block(1:pair_left%nao_j, 1:pair_right%nao_j), k_block)

      DEALLOCATE (v_sph)

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

      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: c_left, c_right, cvc_block, tmp, v_sph

      v_norm = 0.0_dp
      cvc_norm = 0.0_dp

      IF (.NOT. ALLOCATED(pair_left%c_fit)) RETURN
      IF (.NOT. ALLOCATED(pair_right%c_fit)) RETURN

      CALL hfx_lri_build_pair_pair_v_sph(basis_set_RI, pair_left, pair_right, ri_data, v_sph)
      IF (.NOT. ALLOCATED(v_sph)) RETURN
      IF (SIZE(v_sph, 1) /= pair_left%naux_pair .OR. SIZE(v_sph, 2) /= pair_right%naux_pair) THEN
         CPABORT("Inconsistent auxiliary dimension in LRI-HFXk pair-pair V probe")
      END IF

      CALL hfx_lri_flatten_pair_coefficients(pair_left, c_left)
      CALL hfx_lri_flatten_pair_coefficients(pair_right, c_right)
      ALLOCATE (tmp(SIZE(c_left, 1), pair_right%naux_pair), SOURCE=0.0_dp)
      tmp = MATMUL(c_left, v_sph)
      ALLOCATE (cvc_block(SIZE(c_left, 1), SIZE(c_right, 1)), SOURCE=0.0_dp)
      cvc_block = MATMUL(tmp, TRANSPOSE(c_right))
      v_norm = SQRT(SUM(v_sph*v_sph))
      cvc_norm = SQRT(SUM(cvc_block*cvc_block))

      DEALLOCATE (cvc_block, tmp, c_right, c_left, v_sph)

   END SUBROUTINE hfx_lri_contract_pair_pair_cvc_realv

! **************************************************************************************************
!> \brief Returns a cached spherical two-center auxiliary block.
!> \param basis_set_RI ...
!> \param ri_transform_cache ...
!> \param ri_data ...
!> \param integral_context ...
!> \param kind_left ...
!> \param kind_right ...
!> \param r_left ...
!> \param r_right ...
!> \param metric_cache ...
!> \param islot ...
!> \param use_transpose ...
!> \param is_zero ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_metric_cache_get(basis_set_RI, ri_transform_cache, ri_data, &
                                       integral_context, kind_left, kind_right, &
                                       r_left, r_right, metric_cache, islot, &
                                       use_transpose, is_zero)

      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(lri_basis_transform_type), DIMENSION(:), &
         INTENT(IN)                                      :: ri_transform_cache
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: integral_context
      INTEGER, INTENT(IN)                                :: kind_left, kind_right
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: r_left, r_right
      TYPE(hfx_lri_metric_cache_type), INTENT(INOUT)     :: metric_cache
      INTEGER, INTENT(OUT)                               :: islot
      LOGICAL, INTENT(OUT)                               :: use_transpose, is_zero

      INTEGER                                            :: direct_set, first_slot, ientry, &
                                                            last_slot, nsets, oldest_slot, &
                                                            reverse_set
      INTEGER(int_8)                                     :: oldest_stamp
      REAL(dp)                                           :: block_norm2, screen_distance
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: tmp, v_raw
      REAL(dp), DIMENSION(3)                             :: displacement
      TYPE(gto_basis_set_type), POINTER                  :: basis_left, basis_right

      NULLIFY (basis_left, basis_right)
      IF (kind_left < 1 .OR. kind_left > SIZE(basis_set_RI) .OR. &
          kind_right < 1 .OR. kind_right > SIZE(basis_set_RI) .OR. &
          kind_left > SIZE(ri_transform_cache) .OR. &
          kind_right > SIZE(ri_transform_cache)) THEN
         CPABORT("Invalid basis kind in LRI-HFXk metric cache")
      END IF
      use_transpose = .FALSE.
      is_zero = .FALSE.
      islot = 0
      displacement(:) = r_right(:) - r_left(:)

      ! Match the shell-distance test in eri_2center before touching either Libint or the cache.
      ! The kind radius bounds every primitive radius, hence this rejects only blocks for which all
      ! primitive pairs are screened by the existing HFX cutoff.
      screen_distance = ri_transform_cache(kind_left)%radius + &
                        ri_transform_cache(kind_right)%radius + &
                        coulomb_cutoff_screen_factor*ri_data%hfx_pot%cutoff_radius
      IF (SUM(displacement**2) > screen_distance*screen_distance) THEN
         metric_cache%screened = metric_cache%screened + 1_int_8
         is_zero = .TRUE.
         RETURN
      END IF

      IF (.NOT. ALLOCATED(metric_cache%entries)) THEN
         ALLOCATE (metric_cache%entries(lri_metric_cache_capacity))
      END IF
      CPASSERT(MOD(SIZE(metric_cache%entries), lri_metric_cache_ways) == 0)
      nsets = SIZE(metric_cache%entries)/lri_metric_cache_ways
      metric_cache%clock = metric_cache%clock + 1_int_8

      ! Four-way set association bounds each lookup to at most eight comparisons (including the
      ! reverse-symmetry key), instead of scanning the whole cache for every auxiliary-center pair.
      direct_set = hfx_lri_metric_hash(kind_left, kind_right, displacement, nsets)
      first_slot = (direct_set - 1)*lri_metric_cache_ways + 1
      last_slot = first_slot + lri_metric_cache_ways - 1
      DO ientry = first_slot, last_slot
         IF (.NOT. metric_cache%entries(ientry)%valid) CYCLE
         IF (metric_cache%entries(ientry)%kind_left /= kind_left .OR. &
             metric_cache%entries(ientry)%kind_right /= kind_right) CYCLE
         IF (MAXVAL(ABS(metric_cache%entries(ientry)%displacement - displacement)) > &
             lri_metric_cache_position_tol) CYCLE
         metric_cache%entries(ientry)%stamp = metric_cache%clock
         metric_cache%hits = metric_cache%hits + 1_int_8
         islot = ientry
         is_zero = metric_cache%entries(ientry)%zero
         RETURN
      END DO

      ! V_AB(R) = transpose(V_BA(-R)); accepting the reverse key increases reuse without
      ! changing or approximating the metric.
      reverse_set = hfx_lri_metric_hash(kind_right, kind_left, -displacement, nsets)
      first_slot = (reverse_set - 1)*lri_metric_cache_ways + 1
      last_slot = first_slot + lri_metric_cache_ways - 1
      DO ientry = first_slot, last_slot
         IF (.NOT. metric_cache%entries(ientry)%valid) CYCLE
         IF (metric_cache%entries(ientry)%kind_left /= kind_right .OR. &
             metric_cache%entries(ientry)%kind_right /= kind_left) CYCLE
         IF (MAXVAL(ABS(metric_cache%entries(ientry)%displacement + displacement)) > &
             lri_metric_cache_position_tol) CYCLE
         metric_cache%entries(ientry)%stamp = metric_cache%clock
         metric_cache%hits = metric_cache%hits + 1_int_8
         islot = ientry
         use_transpose = .TRUE.
         is_zero = metric_cache%entries(ientry)%zero
         RETURN
      END DO

      metric_cache%misses = metric_cache%misses + 1_int_8
      islot = 0
      first_slot = (direct_set - 1)*lri_metric_cache_ways + 1
      last_slot = first_slot + lri_metric_cache_ways - 1
      oldest_slot = first_slot
      oldest_stamp = HUGE(oldest_stamp)
      DO ientry = first_slot, last_slot
         IF (.NOT. metric_cache%entries(ientry)%valid) THEN
            islot = ientry
            EXIT
         END IF
         IF (metric_cache%entries(ientry)%stamp < oldest_stamp) THEN
            oldest_stamp = metric_cache%entries(ientry)%stamp
            oldest_slot = ientry
         END IF
      END DO
      IF (islot == 0) islot = oldest_slot

      basis_left => basis_set_RI(kind_left)%gto_basis_set
      basis_right => basis_set_RI(kind_right)%gto_basis_set
      IF (.NOT. ASSOCIATED(basis_left) .OR. .NOT. ASSOCIATED(basis_right) .OR. &
          .NOT. ALLOCATED(ri_transform_cache(kind_left)%sphi) .OR. &
          .NOT. ALLOCATED(ri_transform_cache(kind_right)%sphi)) THEN
         CPABORT("Missing basis data in LRI-HFXk metric cache")
      END IF

      IF (ALLOCATED(metric_cache%entries(islot)%values)) THEN
         DEALLOCATE (metric_cache%entries(islot)%values)
      END IF
      metric_cache%entries(islot)%zero = .FALSE.
      ALLOCATE (v_raw(SIZE(ri_transform_cache(kind_left)%sphi, 1), &
                      SIZE(ri_transform_cache(kind_right)%sphi, 1)), SOURCE=0.0_dp)
      CALL hfx_lri_fill_2c_block(basis_left, basis_right, r_left, r_right, integral_context, &
                                 ri_data, v_raw, 0, 0)
      ALLOCATE (tmp(SIZE(v_raw, 1), SIZE(ri_transform_cache(kind_right)%sphi, 2)))
      tmp(:, :) = MATMUL(v_raw, ri_transform_cache(kind_right)%sphi)
      ALLOCATE (metric_cache%entries(islot)%values( &
                SIZE(ri_transform_cache(kind_left)%sphi, 2), &
                SIZE(ri_transform_cache(kind_right)%sphi, 2)))
      metric_cache%entries(islot)%values(:, :) = &
         MATMUL(TRANSPOSE(ri_transform_cache(kind_left)%sphi), tmp)
      DEALLOCATE (tmp, v_raw)
      block_norm2 = SUM(metric_cache%entries(islot)%values**2)
      metric_cache%entries(islot)%norm2 = block_norm2
      ! Cache an exactly zero backend result without retaining its dense storage. No numerical
      ! threshold is applied here; EPS_FILTER_2C is evaluated on the assembled pair block below.
      IF (block_norm2 == 0.0_dp) THEN
         DEALLOCATE (metric_cache%entries(islot)%values)
         metric_cache%entries(islot)%zero = .TRUE.
         is_zero = .TRUE.
      END IF
      metric_cache%entries(islot)%valid = .TRUE.
      metric_cache%entries(islot)%kind_left = kind_left
      metric_cache%entries(islot)%kind_right = kind_right
      metric_cache%entries(islot)%displacement(:) = displacement(:)
      metric_cache%entries(islot)%stamp = metric_cache%clock

   END SUBROUTINE hfx_lri_metric_cache_get

! **************************************************************************************************
!> \brief Builds the metric between every auxiliary center of one oriented pair and one fixed center.
!> \param basis_set_RI ...
!> \param pair_left ...
!> \param left_reversed ...
!> \param right_kind ...
!> \param right_position ...
!> \param right_size ...
!> \param ri_data ...
!> \param v_sph ...
!> \param integral_context ...
!> \param ri_transform_cache ...
!> \param metric_cache ...
!> \param metric_table ...
!> \param v_norm2 ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_build_pair_center_v_sph_context( &
      basis_set_RI, pair_left, left_reversed, right_kind, right_position, right_size, &
      ri_data, v_sph, integral_context, ri_transform_cache, metric_cache, metric_table, v_norm2)

      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair_left
      LOGICAL, INTENT(IN)                                :: left_reversed
      INTEGER, INTENT(IN)                                :: right_kind
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: right_position
      INTEGER, INTENT(IN)                                :: right_size
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: v_sph
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: integral_context
      TYPE(lri_basis_transform_type), DIMENSION(:), &
         INTENT(IN)                                      :: ri_transform_cache
      TYPE(hfx_lri_metric_cache_type), INTENT(INOUT)     :: metric_cache
      TYPE(hfx_lri_metric_table_type), INTENT(IN)        :: metric_table
      REAL(dp), INTENT(OUT)                              :: v_norm2

      INTEGER                                            :: cache_slot, ic_left, kind_left, &
                                                            left_offset, left_size, table_entry
      LOGICAL                                            :: cache_transpose, metric_zero, table_found
      REAL(dp), DIMENSION(3)                             :: left_origin, left_position

      IF (right_kind < 1 .OR. right_kind > SIZE(basis_set_RI) .OR. &
          right_kind > SIZE(ri_transform_cache)) THEN
         CPABORT("Invalid fixed-center kind in LRI-HFXk metric build")
      END IF
      IF (.NOT. ALLOCATED(ri_transform_cache(right_kind)%sphi)) THEN
         CPABORT("Missing fixed-center transformation in LRI-HFXk metric build")
      END IF
      IF (right_size /= SIZE(ri_transform_cache(right_kind)%sphi, 2)) THEN
         CPABORT("Inconsistent fixed-center size in LRI-HFXk metric build")
      END IF
      IF (pair_left%naux_centers <= 0 .OR. pair_left%naux_centers > 2) THEN
         CPABORT("Invalid left pair in LRI-HFXk fixed-center metric build")
      END IF

      IF (ALLOCATED(v_sph)) THEN
         IF (ANY(SHAPE(v_sph) /= [pair_left%naux_pair, right_size])) DEALLOCATE (v_sph)
      END IF
      IF (.NOT. ALLOCATED(v_sph)) ALLOCATE (v_sph(pair_left%naux_pair, right_size))
      v_sph(:, :) = 0.0_dp
      v_norm2 = 0.0_dp
      left_origin(:) = 0.0_dp
      IF (left_reversed) left_origin(:) = pair_left%rab(:)

      DO ic_left = 1, pair_left%naux_centers
         kind_left = pair_left%aux_kind(ic_left)
         left_offset = pair_left%aux_offset(ic_left)
         left_size = pair_left%aux_size(ic_left)
         IF (kind_left < 1 .OR. kind_left > SIZE(basis_set_RI) .OR. &
             kind_left > SIZE(ri_transform_cache) .OR. left_size <= 0 .OR. &
             left_offset < 0 .OR. left_offset + left_size > pair_left%naux_pair) THEN
            CPABORT("Invalid left auxiliary center in LRI-HFXk fixed-center metric build")
         END IF
         left_position(:) = pair_left%aux_r(:, ic_left) - left_origin(:)
         CALL hfx_lri_metric_table_find( &
            metric_table, kind_left, right_kind, right_position - left_position, &
            table_entry, cache_transpose, table_found)
         IF (table_found) THEN
            metric_zero = metric_table%entries(table_entry)%zero
            IF (metric_zero) CYCLE
            IF (cache_transpose) THEN
               IF (ANY(SHAPE(metric_table%entries(table_entry)%values) /= &
                       [right_size, left_size])) THEN
                  CPABORT("Inconsistent transposed shared metric block in LRI-HFXk")
               END IF
               v_sph(left_offset + 1:left_offset + left_size, :) = &
                  TRANSPOSE(metric_table%entries(table_entry)%values)
            ELSE
               IF (ANY(SHAPE(metric_table%entries(table_entry)%values) /= &
                       [left_size, right_size])) THEN
                  CPABORT("Inconsistent shared metric block in LRI-HFXk")
               END IF
               v_sph(left_offset + 1:left_offset + left_size, :) = &
                  metric_table%entries(table_entry)%values
            END IF
            v_norm2 = v_norm2 + metric_table%entries(table_entry)%norm2
         ELSE
            ! A changed density sparsity pattern can introduce a geometry key after the immutable
            ! table was built.  Preserve correctness with the old private cache on that rare path.
            CALL hfx_lri_metric_cache_get( &
               basis_set_RI, ri_transform_cache, ri_data, integral_context, &
               kind_left, right_kind, left_position, right_position, metric_cache, &
               cache_slot, cache_transpose, metric_zero)
            IF (metric_zero) CYCLE
            IF (cache_transpose) THEN
               IF (ANY(SHAPE(metric_cache%entries(cache_slot)%values) /= &
                       [right_size, left_size])) THEN
                  CPABORT("Inconsistent transposed block in LRI-HFXk metric fallback")
               END IF
               v_sph(left_offset + 1:left_offset + left_size, :) = &
                  TRANSPOSE(metric_cache%entries(cache_slot)%values)
            ELSE
               IF (ANY(SHAPE(metric_cache%entries(cache_slot)%values) /= &
                       [left_size, right_size])) THEN
                  CPABORT("Inconsistent block in LRI-HFXk metric fallback")
               END IF
               v_sph(left_offset + 1:left_offset + left_size, :) = &
                  metric_cache%entries(cache_slot)%values
            END IF
            v_norm2 = v_norm2 + metric_cache%entries(cache_slot)%norm2
         END IF
      END DO

   END SUBROUTINE hfx_lri_build_pair_center_v_sph_context

! **************************************************************************************************
!> \brief Conservative pair-level form of the primitive distance screening in eri_2center.
!>        False means that every primitive pair in every auxiliary-center block is screened.
!> \param pair_left ...
!> \param pair_right ...
!> \param left_reversed ...
!> \param right_reversed ...
!> \param right_shift ...
!> \param ri_data ...
!> \param ri_transform_cache ...
!> \return ...
! **************************************************************************************************
   LOGICAL FUNCTION hfx_lri_pair_pair_may_couple(pair_left, pair_right, left_reversed, &
                                                 right_reversed, right_shift, &
                                                 ri_data, ri_transform_cache) RESULT(may_couple)

      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair_left, pair_right
      LOGICAL, INTENT(IN)                                :: left_reversed, right_reversed
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: right_shift
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      TYPE(lri_basis_transform_type), DIMENSION(:), &
         INTENT(IN)                                      :: ri_transform_cache

      INTEGER                                            :: ic_left, ic_right, kind_left, kind_right
      REAL(dp)                                           :: screen_distance
      REAL(dp), DIMENSION(3)                             :: displacement, left_origin, right_origin

      may_couple = .FALSE.
      left_origin(:) = 0.0_dp
      right_origin(:) = 0.0_dp
      IF (left_reversed) left_origin(:) = pair_left%rab(:)
      IF (right_reversed) right_origin(:) = pair_right%rab(:)

      DO ic_left = 1, pair_left%naux_centers
         kind_left = pair_left%aux_kind(ic_left)
         IF (kind_left < 1 .OR. kind_left > SIZE(ri_transform_cache)) THEN
            CPABORT("Invalid left kind in LRI-HFXk pair distance screening")
         END IF
         DO ic_right = 1, pair_right%naux_centers
            kind_right = pair_right%aux_kind(ic_right)
            IF (kind_right < 1 .OR. kind_right > SIZE(ri_transform_cache)) THEN
               CPABORT("Invalid right kind in LRI-HFXk pair distance screening")
            END IF
            displacement(:) = pair_right%aux_r(:, ic_right) - right_origin + right_shift - &
                              pair_left%aux_r(:, ic_left) + left_origin
            screen_distance = ri_transform_cache(kind_left)%radius + &
                              ri_transform_cache(kind_right)%radius + &
                              coulomb_cutoff_screen_factor*ri_data%hfx_pot%cutoff_radius
            IF (SUM(displacement*displacement) <= screen_distance*screen_distance) THEN
               may_couple = .TRUE.
               RETURN
            END IF
         END DO
      END DO

   END FUNCTION hfx_lri_pair_pair_may_couple

! **************************************************************************************************
!> \brief Builds a spherical pair-pair auxiliary coupling block from pair metadata.
!> \param basis_set_RI ...
!> \param pair_left ...
!> \param pair_right ...
!> \param ri_data ...
!> \param v_sph ...
!> \param right_shift ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_build_pair_pair_v_sph(basis_set_RI, pair_left, pair_right, ri_data, v_sph, &
                                            right_shift)

      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair_left, pair_right
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: v_sph
      REAL(dp), DIMENSION(3), INTENT(IN), OPTIONAL       :: right_shift

      INTEGER                                            :: icenter, maxl_all, maxl_center
      TYPE(coulomb_integral_context_type)                :: integral_context
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(lri_basis_transform_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: ri_transform_cache

      NULLIFY (basis)
      maxl_all = 0
      DO icenter = 1, pair_left%naux_centers
         basis => basis_set_RI(pair_left%aux_kind(icenter))%gto_basis_set
         IF (.NOT. ASSOCIATED(basis)) RETURN
         CALL get_gto_basis_set(basis, maxl=maxl_center)
         maxl_all = MAX(maxl_all, maxl_center)
      END DO
      DO icenter = 1, pair_right%naux_centers
         basis => basis_set_RI(pair_right%aux_kind(icenter))%gto_basis_set
         IF (.NOT. ASSOCIATED(basis)) RETURN
         CALL get_gto_basis_set(basis, maxl=maxl_center)
         maxl_all = MAX(maxl_all, maxl_center)
      END DO

      CALL init_md_ftable(nmax=2*maxl_all)
      CALL coulomb_integral_init(integral_context, max_am_2c=maxl_all)
      CALL hfx_lri_build_ri_transform_cache(basis_set_RI, ri_transform_cache)
      CALL hfx_lri_build_pair_pair_v_sph_context( &
         basis_set_RI, pair_left, pair_right, ri_data, v_sph, integral_context, &
         ri_transform_cache, right_shift)
      CALL coulomb_integral_cleanup(integral_context)

   END SUBROUTINE hfx_lri_build_pair_pair_v_sph

! **************************************************************************************************
!> \brief Context-reusing implementation of the pair-pair auxiliary coupling block.
!> \param basis_set_RI ...
!> \param pair_left ...
!> \param pair_right ...
!> \param ri_data ...
!> \param v_sph ...
!> \param integral_context ...
!> \param ri_transform_cache ...
!> \param right_shift ...
!> \param metric_cache ...
!> \param left_reversed ...
!> \param right_reversed ...
!> \param v_norm2 ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_build_pair_pair_v_sph_context(basis_set_RI, pair_left, pair_right, ri_data, &
                                                    v_sph, integral_context, ri_transform_cache, &
                                                    right_shift, metric_cache, &
                                                    left_reversed, right_reversed, v_norm2)

      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair_left, pair_right
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: v_sph
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: integral_context
      TYPE(lri_basis_transform_type), DIMENSION(:), &
         INTENT(IN)                                      :: ri_transform_cache
      REAL(dp), DIMENSION(3), INTENT(IN), OPTIONAL       :: right_shift
      TYPE(hfx_lri_metric_cache_type), INTENT(INOUT), &
         OPTIONAL                                        :: metric_cache
      LOGICAL, INTENT(IN), OPTIONAL                      :: left_reversed, right_reversed
      REAL(dp), INTENT(OUT), OPTIONAL                    :: v_norm2

      INTEGER :: cache_slot, ic_left, ic_right, ikind_left, ikind_right, naux_center_left, &
         naux_center_right, naux_left, naux_right, ncart_left, ncart_right, offset_left, &
         offset_right
      LOGICAL                                            :: cache_transpose, left_t, metric_zero, &
                                                            right_t
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: tmp, v_raw
      REAL(dp), DIMENSION(3)                             :: left_origin, output_shift, right_origin
      TYPE(gto_basis_set_type), POINTER                  :: rbas_left, rbas_right

      NULLIFY (rbas_left, rbas_right)

      IF (pair_left%naux_centers <= 0 .OR. pair_right%naux_centers <= 0) RETURN
      IF (pair_left%naux_centers > 2 .OR. pair_right%naux_centers > 2) RETURN
      left_t = .FALSE.
      right_t = .FALSE.
      IF (PRESENT(left_reversed)) left_t = left_reversed
      IF (PRESENT(right_reversed)) right_t = right_reversed
      output_shift(:) = 0.0_dp
      IF (PRESENT(right_shift)) output_shift(:) = right_shift(:)
      left_origin(:) = 0.0_dp
      right_origin(:) = 0.0_dp
      IF (left_t) THEN
         left_origin(:) = pair_left%rab(:)
      END IF
      IF (right_t) THEN
         right_origin(:) = pair_right%rab(:)
      END IF
      IF (PRESENT(v_norm2)) v_norm2 = 0.0_dp

      naux_left = 0
      naux_right = 0

      DO ic_left = 1, pair_left%naux_centers
         rbas_left => basis_set_RI(pair_left%aux_kind(ic_left))%gto_basis_set
         IF (.NOT. ASSOCIATED(rbas_left)) RETURN
         naux_left = naux_left + pair_left%aux_size(ic_left)
      END DO
      DO ic_right = 1, pair_right%naux_centers
         rbas_right => basis_set_RI(pair_right%aux_kind(ic_right))%gto_basis_set
         IF (.NOT. ASSOCIATED(rbas_right)) RETURN
         naux_right = naux_right + pair_right%aux_size(ic_right)
      END DO

      IF (naux_left /= pair_left%naux_pair .OR. naux_right /= pair_right%naux_pair) THEN
         CPABORT("Inconsistent auxiliary metadata in LRI-HFXk pair-pair V build")
      END IF

      IF (ALLOCATED(v_sph)) THEN
         IF (ANY(SHAPE(v_sph) /= [naux_left, naux_right])) DEALLOCATE (v_sph)
      END IF
      IF (.NOT. ALLOCATED(v_sph)) ALLOCATE (v_sph(naux_left, naux_right))
      v_sph(:, :) = 0.0_dp
      DO ic_left = 1, pair_left%naux_centers
         ikind_left = pair_left%aux_kind(ic_left)
         IF (ikind_left < 1 .OR. ikind_left > SIZE(basis_set_RI) .OR. &
             ikind_left > SIZE(ri_transform_cache)) THEN
            CPABORT("Invalid left auxiliary kind in LRI-HFXk pair-pair V build")
         END IF
         rbas_left => basis_set_RI(ikind_left)%gto_basis_set
         ncart_left = hfx_lri_count_cart_functions(rbas_left)
         naux_center_left = pair_left%aux_size(ic_left)
         offset_left = pair_left%aux_offset(ic_left)
         IF (.NOT. ALLOCATED(ri_transform_cache(ikind_left)%sphi)) RETURN
         IF (ANY(SHAPE(ri_transform_cache(ikind_left)%sphi) /= &
                 [ncart_left, naux_center_left])) THEN
            CPABORT("Inconsistent left auxiliary transformation in LRI-HFXk pair-pair V build")
         END IF
         DO ic_right = 1, pair_right%naux_centers
            ikind_right = pair_right%aux_kind(ic_right)
            IF (ikind_right < 1 .OR. ikind_right > SIZE(basis_set_RI) .OR. &
                ikind_right > SIZE(ri_transform_cache)) THEN
               CPABORT("Invalid right auxiliary kind in LRI-HFXk pair-pair V build")
            END IF
            rbas_right => basis_set_RI(ikind_right)%gto_basis_set
            ncart_right = hfx_lri_count_cart_functions(rbas_right)
            naux_center_right = pair_right%aux_size(ic_right)
            offset_right = pair_right%aux_offset(ic_right)
            IF (.NOT. ALLOCATED(ri_transform_cache(ikind_right)%sphi)) RETURN
            IF (ANY(SHAPE(ri_transform_cache(ikind_right)%sphi) /= &
                    [ncart_right, naux_center_right])) THEN
               CPABORT("Inconsistent right auxiliary transformation in LRI-HFXk pair-pair V build")
            END IF

            IF (PRESENT(metric_cache)) THEN
               CALL hfx_lri_metric_cache_get( &
                  basis_set_RI, ri_transform_cache, ri_data, integral_context, &
                  ikind_left, ikind_right, pair_left%aux_r(:, ic_left) - left_origin, &
                  pair_right%aux_r(:, ic_right) - right_origin + output_shift, &
                  metric_cache, cache_slot, &
                  cache_transpose, metric_zero)
               IF (metric_zero) THEN
                  CYCLE
               ELSE IF (cache_transpose) THEN
                  v_sph(offset_left + 1:offset_left + naux_center_left, &
                        offset_right + 1:offset_right + naux_center_right) = &
                     TRANSPOSE(metric_cache%entries(cache_slot)%values)
               ELSE
                  v_sph(offset_left + 1:offset_left + naux_center_left, &
                        offset_right + 1:offset_right + naux_center_right) = &
                     metric_cache%entries(cache_slot)%values
               END IF
               IF (PRESENT(v_norm2)) THEN
                  v_norm2 = v_norm2 + metric_cache%entries(cache_slot)%norm2
               END IF
            ELSE
               ALLOCATE (v_raw(ncart_left, ncart_right), SOURCE=0.0_dp)
               CALL hfx_lri_fill_2c_block(rbas_left, rbas_right, &
                                          pair_left%aux_r(:, ic_left) - left_origin, &
                                          pair_right%aux_r(:, ic_right) - right_origin + output_shift, &
                                          integral_context, ri_data, v_raw, 0, 0)
               ALLOCATE (tmp(ncart_left, naux_center_right))
               tmp(:, :) = MATMUL(v_raw, ri_transform_cache(ikind_right)%sphi)
               v_sph(offset_left + 1:offset_left + naux_center_left, &
                     offset_right + 1:offset_right + naux_center_right) = &
                  MATMUL(TRANSPOSE(ri_transform_cache(ikind_left)%sphi), tmp)
               IF (PRESENT(v_norm2)) THEN
                  v_norm2 = v_norm2 + SUM(v_sph(offset_left + 1:offset_left + naux_center_left, &
                                                offset_right + 1:offset_right + naux_center_right)**2)
               END IF
               DEALLOCATE (tmp, v_raw)
            END IF
         END DO
      END DO
   END SUBROUTINE hfx_lri_build_pair_pair_v_sph_context

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
      DO iao = 1, pair_block%nao_i
         DO jao = 1, pair_block%nao_j
            irow = hfx_lri_product_index(iao, jao, pair_block%nao_j)
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
