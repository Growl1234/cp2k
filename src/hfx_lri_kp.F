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
   USE cell_types, ONLY: cell_type
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
   USE hfx_lri_contract, ONLY: hfx_lri_contract_workspace_type, &
                               hfx_lri_loop3_aux_transform, &
                               hfx_lri_loop3_central_contract, &
                               hfx_lri_loop3_density_c_add, &
                               hfx_lri_loop3_pack_aux, &
                               hfx_lri_loop3_remote_contract, &
                               hfx_lri_product_index
   USE hfx_lri_types, ONLY: hfx_lri_kp_release, &
                            hfx_lri_centered_block_type, &
                            hfx_lri_kp_type, &
                            hfx_lri_pair_block_type
   USE hfx_types, ONLY: hfx_ri_type
   USE input_constants, ONLY: do_potential_short
   USE integral_library_types, ONLY: active_integral_library, &
                                     library_libint
   USE kinds, ONLY: dp, &
                    int_8
   USE kpoint_types, ONLY: get_kpoint_info, &
                           kpoint_type
   USE mathconstants, ONLY: twopi
   USE message_passing, ONLY: mp_para_env_type
   USE orbital_pointers, ONLY: ncoset
   USE particle_types, ONLY: particle_type
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

   INTEGER, PARAMETER :: lri_metric_cache_capacity = 128
   INTEGER, PARAMETER :: lri_metric_cache_ways = 4
   REAL(dp), PARAMETER :: lri_metric_cache_position_tol = 1.0E-12_dp

   TYPE hfx_lri_metric_entry_type
      LOGICAL                                            :: valid = .FALSE., zero = .FALSE.
      INTEGER                                            :: kind_left = 0, kind_right = 0
      INTEGER(int_8)                                     :: stamp = 0_int_8
      REAL(dp), DIMENSION(3)                             :: displacement = 0.0_dp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: values
   END TYPE hfx_lri_metric_entry_type

   TYPE hfx_lri_metric_cache_type
      INTEGER(int_8)                                     :: clock = 0_int_8
      TYPE(hfx_lri_metric_entry_type), ALLOCATABLE, DIMENSION(:) :: entries
   END TYPE hfx_lri_metric_cache_type

   TYPE hfx_lri_metric_table_type
      INTEGER                                            :: capacity = 0, nentries = 0
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: slot
      TYPE(hfx_lri_metric_entry_type), ALLOCATABLE, DIMENSION(:) :: entries
   END TYPE hfx_lri_metric_table_type

   TYPE lri_pair_lookup_type
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_count, atom_offset, atom_pairs
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

   TYPE lri_basis_transform_type
      REAL(dp)                                           :: radius = 0.0_dp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: sphi
   END TYPE lri_basis_transform_type

   TYPE lri_loop3_density_block_type
      REAL(dp)                                           :: norm2 = 0.0_dp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: values
   END TYPE lri_loop3_density_block_type

CONTAINS

! **************************************************************************************************
!> \brief Moves reusable buffers between persistent and thread-private contraction workspaces.
!> \param source ...
!> \param destination ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_contract_workspace_move(source, destination)

      TYPE(hfx_lri_contract_workspace_type), &
         INTENT(INOUT)                                   :: source, destination

      CALL MOVE_ALLOC(source%loop3_s_flat, destination%loop3_s_flat)
      CALL MOVE_ALLOC(source%loop3_t_flat, destination%loop3_t_flat)
      CALL MOVE_ALLOC(source%loop3_u_flat, destination%loop3_u_flat)
      CALL MOVE_ALLOC(source%loop3_block, destination%loop3_block)
      CALL MOVE_ALLOC(source%loop3_s, destination%loop3_s)
      CALL MOVE_ALLOC(source%loop3_t, destination%loop3_t)
      CALL MOVE_ALLOC(source%loop3_u, destination%loop3_u)
      CALL MOVE_ALLOC(source%loop3_w, destination%loop3_w)

   END SUBROUTINE hfx_lri_contract_workspace_move

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

      INTEGER :: handle, iatom, ibasis, ikind, ipair_block, ithread, jatom, jkind, maxl_ao, &
         maxl_ao_kind, maxl_ri, maxl_ri_kind, mepos, mmax, nao_i, nao_j, naux_i, naux_j, &
         naux_pair, nkind, npair_blocks_local, npair_blocks_max, npair_blocks_sum, &
         npair_candidates, npairs_local, npairs_max, npairs_sum, nthreads
      INTEGER(int_8) :: c_nreal_local, c_nreal_max, c_nreal_stored_local, c_nreal_stored_max, &
         c_nreal_stored_sum, c_nreal_sum, sum_nao_pair_local, sum_nao_pair_max, sum_nao_pair_sum, &
         sum_naux_pair_local, sum_naux_pair_max, sum_naux_pair_sum
      INTEGER, DIMENSION(3)                              :: cell_j
      REAL(dp) :: c_mem_mib_local, c_mem_mib_max, c_mem_mib_stored_local, c_mem_mib_stored_max, &
         c_mem_mib_stored_sum, c_mem_mib_sum
      REAL(dp), DIMENSION(3)                             :: rab
      TYPE(coulomb_integral_context_type)                :: integral_context
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_AO, basis_set_RI
      TYPE(gto_basis_set_type), POINTER                  :: obas_i, obas_j, rbas_i, rbas_j
      TYPE(hfx_lri_metric_cache_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: metric_caches
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

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, nkind=nkind, para_env=para_env, &
                      dft_control=dft_control)

      CALL hfx_lri_kp_release(ri_data%lri_kp)
      ALLOCATE (ri_data%lri_kp)
      nthreads = 1
!$    nthreads = omp_get_max_threads()
      ALLOCATE (metric_caches(nthreads))

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

      ! Build the symmetric AO-pair neighbor list with the LRI screening radii.
      CALL setup_neighbor_list(soo_list, basis_set_AO, qs_env=qs_env, mic=.FALSE., symmetric=.TRUE.)

      ! Count first so that every significant pair is retained without a fixed-size buffer.
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

      npairs_local = 0
      sum_nao_pair_local = 0_int_8
      sum_naux_pair_local = 0_int_8
      c_nreal_local = 0_int_8
      mepos = 0
      DO WHILE (neighbor_list_iterate(nl_iterator, mepos) == 0)

         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, &
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

         ipair_block = ri_data%lri_kp%npair_blocks + 1
         CALL hfx_lri_init_pair_block(ri_data%lri_kp%pair_blocks(ipair_block), iatom, jatom, &
                                      ikind, jkind, cell_j, rab, nao_i, nao_j, naux_i, naux_j, &
                                      naux_pair)
         ri_data%lri_kp%npair_blocks = ipair_block

         npairs_local = npairs_local + 1
         sum_nao_pair_local = sum_nao_pair_local + INT(nao_i, int_8)*INT(nao_j, int_8)
         sum_naux_pair_local = sum_naux_pair_local + INT(naux_pair, int_8)
         c_nreal_local = c_nreal_local + &
                         INT(nao_i, int_8)*INT(nao_j, int_8)*INT(naux_pair, int_8)

      END DO

      CALL neighbor_list_iterator_release(nl_iterator)
      CALL release_neighbor_list_sets(soo_list)

      ! The neighbor-list iterator is serial, but every fitted pair is independent. Reuse the
      ! immutable basis transforms and give every worker its own Libint engines.
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP          SHARED(ao_transform_cache, basis_set_AO, basis_set_RI, maxl_ao, maxl_ri, &
      !$OMP                 metric_caches, ri_data, ri_transform_cache) &
      !$OMP          PRIVATE(integral_context, ipair_block, ithread)
      ithread = 1
!$    ithread = omp_get_thread_num() + 1
      metric_caches(ithread)%clock = 0_int_8
      CALL coulomb_integral_init(integral_context, max_am_2c=maxl_ri, &
                                 max_am_3c=MAX(maxl_ao, maxl_ri))
      !$OMP DO SCHEDULE(DYNAMIC, 1)
      DO ipair_block = 1, ri_data%lri_kp%npair_blocks
         CALL hfx_lri_build_pair_coefficients( &
            basis_set_AO, basis_set_RI, ri_data%lri_kp%pair_blocks(ipair_block), ri_data, &
            integral_context, ao_transform_cache, ri_transform_cache, metric_caches(ithread))
      END DO
      !$OMP END DO
      CALL coulomb_integral_cleanup(integral_context)
      !$OMP END PARALLEL
      DEALLOCATE (metric_caches)

      c_nreal_stored_local = 0_int_8
      npair_blocks_local = 0
      DO ipair_block = 1, ri_data%lri_kp%npair_blocks
         IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks(ipair_block)%c_fit)) CYCLE
         npair_blocks_local = npair_blocks_local + 1
         c_nreal_stored_local = c_nreal_stored_local + &
                                INT(SIZE(ri_data%lri_kp%pair_blocks(ipair_block)%c_fit), int_8)
      END DO

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
      END IF

      CALL timestop(handle)

   END SUBROUTINE hfx_lri_pre_scf_kp

! **************************************************************************************************
!> \brief Counts raw Cartesian primitive-function slots in one basis.
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
!> \brief Fills a raw Cartesian two-center integral block in a dense metric.
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
!> \brief Fills a raw Cartesian three-center integral block in a dense J matrix.
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

      CALL hfx_lri_metric_inverse(metric_sph, j_sph, ri_data%eps_eigval, &
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
   SUBROUTINE hfx_lri_metric_inverse(metric_block, j_block, eig_threshold, &
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

   END SUBROUTINE hfx_lri_metric_inverse

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
!> \brief Transfers an allocated K block together with its target-pair metadata.
!> \param k_pair ...
!> \param pair_block ...
!> \param k_block ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_store_k_block(k_pair, pair_block, k_block)

      USE hfx_lri_types, ONLY: hfx_lri_k_block_type
      TYPE(hfx_lri_k_block_type), INTENT(INOUT)          :: k_pair
      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair_block
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: k_block

      IF (ALLOCATED(k_pair%k_block)) DEALLOCATE (k_pair%k_block)
      IF (.NOT. ALLOCATED(k_block)) CPABORT("Cannot store an unallocated LRI-HFXk K block")

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
!> \brief Reduces a lattice vector modulo a regular k mesh and returns its real boundary phase.
!>        Standard even Monkhorst-Pack meshes are anti-periodic and therefore return -1 whenever
!>        an odd number of mesh periods is crossed.
!> \param cell ...
!> \param grid ...
!> \param twist ...
!> \param reduced ...
!> \param factor ...
! **************************************************************************************************
   PURE SUBROUTINE hfx_lri_lattice_reduce(cell, grid, twist, reduced, factor)

      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell, grid, twist
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: reduced
      REAL(dp), INTENT(OUT), OPTIONAL                    :: factor

      INTEGER                                            :: idir, periods
      REAL(dp)                                           :: phase

      phase = 1.0_dp
      DO idir = 1, 3
         reduced(idir) = MODULO(cell(idir), grid(idir))
         periods = (cell(idir) - reduced(idir))/grid(idir)
         IF (twist(idir) == -1 .AND. MODULO(periods, 2) /= 0) phase = -phase
      END DO
      IF (PRESENT(factor)) factor = phase

   END SUBROUTINE hfx_lri_lattice_reduce

! **************************************************************************************************
!> \brief Linear index of a reduced lattice cell, with the first direction varying fastest.
!> \param cell ...
!> \param grid ...
!> \return ...
! **************************************************************************************************
   PURE INTEGER FUNCTION hfx_lri_lattice_cell_index(cell, grid) RESULT(index)

      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell, grid

      index = 1 + cell(1) + grid(1)*(cell(2) + grid(2)*cell(3))

   END FUNCTION hfx_lri_lattice_cell_index

! **************************************************************************************************
!> \brief Decodes a reduced-cell linear index.
!> \param index ...
!> \param grid ...
!> \param cell ...
! **************************************************************************************************
   PURE SUBROUTINE hfx_lri_lattice_cell_from_index(index, grid, cell)

      INTEGER, INTENT(IN)                                :: index
      INTEGER, DIMENSION(3), INTENT(IN)                  :: grid
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: cell

      INTEGER                                            :: work

      work = index - 1
      cell(1) = MODULO(work, grid(1))
      work = work/grid(1)
      cell(2) = MODULO(work, grid(2))
      cell(3) = work/grid(2)

   END SUBROUTINE hfx_lri_lattice_cell_from_index

! **************************************************************************************************
!> \brief Direct index of an oriented atom-pair block on the finite lattice quotient.
!> \param iatom ...
!> \param jatom ...
!> \param cell_index ...
!> \param natom ...
!> \param ncell ...
!> \return ...
! **************************************************************************************************
   PURE INTEGER FUNCTION hfx_lri_lattice_block_index(iatom, jatom, cell_index, natom, ncell) &
      RESULT(index)

      INTEGER, INTENT(IN)                                :: iatom, jatom, cell_index, natom, ncell

      index = cell_index + ncell*((jatom - 1) + natom*(iatom - 1))

   END FUNCTION hfx_lri_lattice_block_index

! **************************************************************************************************
!> \brief Decodes an oriented atom-pair/lattice block index.
!> \param index ...
!> \param natom ...
!> \param ncell ...
!> \param iatom ...
!> \param jatom ...
!> \param cell_index ...
! **************************************************************************************************
   PURE SUBROUTINE hfx_lri_lattice_block_from_index(index, natom, ncell, iatom, jatom, cell_index)

      INTEGER, INTENT(IN)                                :: index, natom, ncell
      INTEGER, INTENT(OUT)                               :: iatom, jatom, cell_index

      INTEGER                                            :: work

      cell_index = MODULO(index - 1, ncell) + 1
      work = (index - 1)/ncell
      jatom = MODULO(work, natom) + 1
      iatom = work/natom + 1

   END SUBROUTINE hfx_lri_lattice_block_from_index

! **************************************************************************************************
!> \brief Establishes the finite lattice quotient and its real Monkhorst-Pack boundary twist.
!> \param kpoints ...
!> \param grid ...
!> \param twist ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_lattice_grid(kpoints, grid, twist)

      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: grid, twist

      CHARACTER(LEN=20)                                  :: kp_scheme
      INTEGER                                            :: idir
      LOGICAL                                            :: gamma_centered
      REAL(dp)                                           :: angle, mesh_shift
      REAL(dp), DIMENSION(3)                             :: kp_shift

      IF (.NOT. ASSOCIATED(kpoints)) CPABORT("LRI loop3 requires k-point metadata")
      CALL get_kpoint_info(kpoints, kp_scheme=kp_scheme, nkp_grid=grid, &
                           kp_shift=kp_shift, gamma_centered=gamma_centered)
      IF (TRIM(kp_scheme) /= "MONKHORST-PACK" .AND. TRIM(kp_scheme) /= "MACDONALD") THEN
         CPABORT("LRI loop3 currently requires a regular Monkhorst-Pack or MacDonald mesh")
      END IF
      IF (ANY(grid <= 0)) CPABORT("Invalid k-point grid in LRI loop3")

      DO idir = 1, 3
         mesh_shift = REAL(grid(idir), dp)*kp_shift(idir)
         IF (MODULO(grid(idir), 2) == 0 .AND. .NOT. gamma_centered) THEN
            mesh_shift = mesh_shift + 0.5_dp
         END IF
         angle = twopi*mesh_shift
         IF (ABS(SIN(angle)) > 2.0E-10_dp) THEN
            CPABORT("LRI loop3 needs a real (+/-1) k-mesh boundary phase")
         END IF
         twist(idir) = NINT(COS(angle))
         IF (ABS(twist(idir)) /= 1) CPABORT("Invalid real boundary phase in LRI loop3")
      END DO

   END SUBROUTINE hfx_lri_lattice_grid

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
      TYPE(hfx_lri_metric_entry_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: new_entries

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
      metric_table%entries(ientry)%kind_left = kind_left
      metric_table%entries(ientry)%kind_right = kind_right
      metric_table%entries(ientry)%displacement(:) = displacement(:)
      metric_table%slot(islot) = ientry

   END SUBROUTINE hfx_lri_metric_table_add_key

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

      INTEGER                                            :: capacity, existing_iatom, &
                                                            existing_jatom, iatom, ipair, islot, &
                                                            nprobe, nrefs, pair_ref, ref_iatom, &
                                                            ref_jatom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cursor, slot
      INTEGER, DIMENSION(3)                              :: existing_cell, ref_cell
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: include_reverse

      CPASSERT(npair_blocks >= 0)
      CPASSERT(npair_blocks <= SIZE(pair_blocks))
      CPASSERT(natom > 0)

      ! The stored neighbor list is symmetric and therefore contains only one representative of
      ! each pair. Add signed reverse aliases to the lookup without copying any C tensor.
      capacity = MAX(17, 4*npair_blocks + 1)
      ALLOCATE (slot(capacity), SOURCE=0)
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
         islot = hfx_lri_pair_hash(ref_iatom, ref_jatom, ref_cell, capacity)
         DO nprobe = 1, capacity
            IF (slot(islot) == 0) THEN
               slot(islot) = pair_ref
               EXIT
            END IF
            CALL hfx_lri_pair_ref_key(pair_blocks, slot(islot), &
                                      existing_iatom, existing_jatom, existing_cell)
            IF (existing_iatom == ref_iatom .AND. existing_jatom == ref_jatom .AND. &
                ALL(existing_cell == ref_cell)) THEN
               CPABORT("Duplicate atom-pair/cell key in the LRI-HFXk pair list")
            END IF
            islot = MOD(islot, capacity) + 1
         END DO
         IF (nprobe > capacity) THEN
            CPABORT("LRI-HFXk pair lookup table is full")
         END IF
      END DO

      ! Add missing reverse orientations as negative references. Self pairs in the origin cell and
      ! any unexpectedly explicit reverse pair already have a direct entry and are not duplicated.
      DO ipair = 1, npair_blocks
         IF (.NOT. ALLOCATED(pair_blocks(ipair)%c_fit)) CYCLE
         pair_ref = -ipair
         CALL hfx_lri_pair_ref_key(pair_blocks, pair_ref, ref_iatom, ref_jatom, ref_cell)
         islot = hfx_lri_pair_hash(ref_iatom, ref_jatom, ref_cell, capacity)
         DO nprobe = 1, capacity
            IF (slot(islot) == 0) THEN
               slot(islot) = pair_ref
               include_reverse(ipair) = .TRUE.
               pair_lookup%atom_count(ref_iatom) = pair_lookup%atom_count(ref_iatom) + 1
               EXIT
            END IF
            CALL hfx_lri_pair_ref_key(pair_blocks, slot(islot), &
                                      existing_iatom, existing_jatom, existing_cell)
            IF (existing_iatom == ref_iatom .AND. existing_jatom == ref_jatom .AND. &
                ALL(existing_cell == ref_cell)) EXIT
            islot = MOD(islot, capacity) + 1
         END DO
         IF (nprobe > capacity) THEN
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
      DEALLOCATE (slot)

   END SUBROUTINE hfx_lri_pair_lookup_create

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
      DEALLOCATE (density_cache%slot)
      density_cache%capacity = 0
   END SUBROUTINE hfx_lri_density_cache_create

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

      TYPE(hfx_lri_metric_entry_type), INTENT(INOUT)     :: metric_entry
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(lri_basis_transform_type), DIMENSION(:), &
         INTENT(IN)                                      :: ri_transform_cache
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: integral_context

      INTEGER                                            :: kind_left, kind_right
      REAL(dp)                                           :: block_norm2, screen_distance
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: tmp, v_raw
      REAL(dp), DIMENSION(3)                             :: r0
      TYPE(gto_basis_set_type), POINTER                  :: basis_left, basis_right

      NULLIFY (basis_left, basis_right)
      kind_left = metric_entry%kind_left
      kind_right = metric_entry%kind_right
      r0(:) = 0.0_dp
      metric_entry%zero = .FALSE.
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
      block_norm2 = SUM(metric_entry%values**2)
      IF (block_norm2 == 0.0_dp) THEN
         DEALLOCATE (metric_entry%values)
         metric_entry%zero = .TRUE.
      END IF

   END SUBROUTINE hfx_lri_metric_table_compute_entry

! **************************************************************************************************
!> \brief Adds one signed physical fitting coefficient to its centred k-mesh quotient block.
!> \param block ...
!> \param pair ...
!> \param reversed ...
!> \param factor ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_centered_add(block, pair, reversed, factor)

      TYPE(hfx_lri_centered_block_type), INTENT(INOUT)   :: block
      TYPE(hfx_lri_pair_block_type), INTENT(IN)          :: pair
      LOGICAL, INTENT(IN)                                :: reversed
      REAL(dp), INTENT(IN)                               :: factor

      INTEGER                                            :: iaux, naux, ncenter, noffset, nremote

      IF (reversed) THEN
         IF (pair%onsite_pair) CPABORT("On-site LRI coefficient cannot have a reverse slice")
         ncenter = pair%nao_j
         nremote = pair%nao_i
         naux = pair%naux_j
         noffset = pair%naux_i
      ELSE
         ncenter = pair%nao_i
         nremote = pair%nao_j
         naux = pair%naux_i
         noffset = 0
      END IF
      IF (noffset + naux > SIZE(pair%c_fit, 3)) THEN
         CPABORT("Invalid auxiliary slice while periodizing LRI coefficients")
      END IF
      IF (.NOT. ALLOCATED(block%values)) THEN
         ALLOCATE (block%values(ncenter, nremote, naux), SOURCE=0.0_dp)
      ELSE IF (ANY(SHAPE(block%values) /= [ncenter, nremote, naux])) THEN
         CPABORT("Inconsistent aliased coefficient dimensions in LRI loop3")
      END IF
      DO iaux = 1, naux
         IF (reversed) THEN
            block%values(:, :, iaux) = block%values(:, :, iaux) + &
                                       factor*TRANSPOSE(pair%c_fit(:, :, noffset + iaux))
         ELSE
            block%values(:, :, iaux) = block%values(:, :, iaux) + &
                                       factor*pair%c_fit(:, :, noffset + iaux)
         END IF
      END DO

   END SUBROUTINE hfx_lri_centered_add

! **************************************************************************************************
!> \brief Builds atom-centred C blocks on the finite, twisted lattice quotient.
!> \param lri_kp ...
!> \param natom ...
!> \param grid ...
!> \param twist ...
!> \param filter_eps ...
!> \param stored_blocks ...
!> \param stored_reals ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_centered_build(lri_kp, natom, grid, twist, filter_eps, &
                                     stored_blocks, stored_reals)

      TYPE(hfx_lri_kp_type), INTENT(INOUT)               :: lri_kp
      INTEGER, INTENT(IN)                                :: natom
      INTEGER, DIMENSION(3), INTENT(IN)                  :: grid, twist
      REAL(dp), INTENT(IN)                               :: filter_eps
      INTEGER, INTENT(OUT)                               :: stored_blocks
      INTEGER(int_8), INTENT(OUT)                        :: stored_reals

      INTEGER                                            :: block_index, cell_index, iatom, iaux, &
                                                            ipair, jatom, nactive, ncell, &
                                                            pair_cell(3), reduced(3)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cursor
      REAL(dp)                                           :: factor

      ncell = PRODUCT(grid)
      lri_kp%lattice_natom = natom
      lri_kp%lattice_ncell = ncell
      lri_kp%lattice_grid(:) = grid(:)
      lri_kp%lattice_twist(:) = twist(:)
      ALLOCATE (lri_kp%atom_kind(natom), lri_kp%atom_nao(natom), lri_kp%atom_naux(natom), &
                lri_kp%atom_ao_offset(natom + 1), SOURCE=0)
      DO ipair = 1, lri_kp%npair_blocks
         iatom = lri_kp%pair_blocks(ipair)%iatom
         jatom = lri_kp%pair_blocks(ipair)%jatom
         lri_kp%atom_kind(iatom) = lri_kp%pair_blocks(ipair)%ikind
         lri_kp%atom_nao(iatom) = lri_kp%pair_blocks(ipair)%nao_i
         lri_kp%atom_naux(iatom) = lri_kp%pair_blocks(ipair)%naux_i
         lri_kp%atom_kind(jatom) = lri_kp%pair_blocks(ipair)%jkind
         lri_kp%atom_nao(jatom) = lri_kp%pair_blocks(ipair)%nao_j
         lri_kp%atom_naux(jatom) = lri_kp%pair_blocks(ipair)%naux_j
      END DO
      IF (ANY(lri_kp%atom_kind <= 0) .OR. ANY(lri_kp%atom_nao <= 0) .OR. &
          ANY(lri_kp%atom_naux <= 0)) THEN
         CPABORT("Incomplete atom metadata while building LRI loop3 coefficients")
      END IF
      lri_kp%atom_ao_offset(1) = 1
      DO iatom = 1, natom
         lri_kp%atom_ao_offset(iatom + 1) = lri_kp%atom_ao_offset(iatom) + lri_kp%atom_nao(iatom)
      END DO
      lri_kp%lattice_total_nao = lri_kp%atom_ao_offset(natom + 1) - 1

      ALLOCATE (lri_kp%centered_blocks(natom*natom*ncell))
      DO ipair = 1, lri_kp%npair_blocks
         IF (.NOT. ALLOCATED(lri_kp%pair_blocks(ipair)%c_fit)) CYCLE

         iatom = lri_kp%pair_blocks(ipair)%iatom
         jatom = lri_kp%pair_blocks(ipair)%jatom
         pair_cell(:) = lri_kp%pair_blocks(ipair)%cell(:)
         CALL hfx_lri_lattice_reduce(pair_cell, grid, twist, reduced, factor)
         cell_index = hfx_lri_lattice_cell_index(reduced, grid)
         block_index = hfx_lri_lattice_block_index(iatom, jatom, cell_index, natom, ncell)
         CALL hfx_lri_centered_add(lri_kp%centered_blocks(block_index), &
                                   lri_kp%pair_blocks(ipair), .FALSE., factor)

         IF (.NOT. lri_kp%pair_blocks(ipair)%onsite_pair) THEN
            pair_cell(:) = -lri_kp%pair_blocks(ipair)%cell(:)
            CALL hfx_lri_lattice_reduce(pair_cell, grid, twist, reduced, factor)
            cell_index = hfx_lri_lattice_cell_index(reduced, grid)
            block_index = hfx_lri_lattice_block_index(jatom, iatom, cell_index, natom, ncell)
            CALL hfx_lri_centered_add(lri_kp%centered_blocks(block_index), &
                                      lri_kp%pair_blocks(ipair), .TRUE., factor)
         END IF
      END DO

      ! A quotient self block may contain both the genuine on-site coefficient and folded
      ! non-zero translations.  Only the former has no distinct "inner" AO orientation.
      ! Retain a small override for that block instead of dropping the complete self orbit.
      DO ipair = 1, lri_kp%npair_blocks
         IF (.NOT. ALLOCATED(lri_kp%pair_blocks(ipair)%c_fit)) CYCLE
         IF (.NOT. lri_kp%pair_blocks(ipair)%onsite_pair) CYCLE
         cell_index = hfx_lri_lattice_cell_index([0, 0, 0], grid)
         block_index = hfx_lri_lattice_block_index( &
                       lri_kp%pair_blocks(ipair)%iatom, lri_kp%pair_blocks(ipair)%iatom, &
                       cell_index, natom, ncell)
         IF (.NOT. ALLOCATED(lri_kp%centered_blocks(block_index)%values)) THEN
            CPABORT("Missing on-site quotient C block in LRI loop3")
         END IF
         IF (.NOT. ALLOCATED(lri_kp%centered_blocks(block_index)%inner_values)) THEN
            ALLOCATE (lri_kp%centered_blocks(block_index)%inner_values, &
                      SOURCE=lri_kp%centered_blocks(block_index)%values)
         END IF
         lri_kp%centered_blocks(block_index)%inner_values = &
            lri_kp%centered_blocks(block_index)%inner_values - &
            lri_kp%pair_blocks(ipair)%c_fit(:, :, 1:lri_kp%pair_blocks(ipair)%naux_i)
      END DO

      stored_blocks = 0
      stored_reals = 0_int_8
      ALLOCATE (lri_kp%centered_atom_count(natom), lri_kp%centered_remote_count(natom), SOURCE=0)
      DO block_index = 1, SIZE(lri_kp%centered_blocks)
         IF (.NOT. ALLOCATED(lri_kp%centered_blocks(block_index)%values)) CYCLE
         lri_kp%centered_blocks(block_index)%norm2 = &
            SUM(lri_kp%centered_blocks(block_index)%values**2)
         IF (lri_kp%centered_blocks(block_index)%norm2 < filter_eps*filter_eps) THEN
            DEALLOCATE (lri_kp%centered_blocks(block_index)%values)
            IF (ALLOCATED(lri_kp%centered_blocks(block_index)%inner_values)) THEN
               DEALLOCATE (lri_kp%centered_blocks(block_index)%inner_values)
            END IF
            lri_kp%centered_blocks(block_index)%inner_active = .FALSE.
            lri_kp%centered_blocks(block_index)%norm2 = 0.0_dp
            CYCLE
         END IF
         ALLOCATE (lri_kp%centered_blocks(block_index)%remote_flat( &
                   SIZE(lri_kp%centered_blocks(block_index)%values, 1)* &
                   SIZE(lri_kp%centered_blocks(block_index)%values, 3), &
                   SIZE(lri_kp%centered_blocks(block_index)%values, 2)))
         DO iaux = 1, SIZE(lri_kp%centered_blocks(block_index)%values, 3)
            lri_kp%centered_blocks(block_index)%remote_flat( &
               (iaux - 1)*SIZE(lri_kp%centered_blocks(block_index)%values, 1) + 1: &
               iaux*SIZE(lri_kp%centered_blocks(block_index)%values, 1), :) = &
               lri_kp%centered_blocks(block_index)%values(:, :, iaux)
         END DO
         IF (ALLOCATED(lri_kp%centered_blocks(block_index)%inner_values)) THEN
            IF (MAXVAL(ABS(lri_kp%centered_blocks(block_index)%inner_values)) == 0.0_dp) THEN
               lri_kp%centered_blocks(block_index)%inner_active = .FALSE.
               DEALLOCATE (lri_kp%centered_blocks(block_index)%inner_values)
            ELSE
               ALLOCATE (lri_kp%centered_blocks(block_index)%inner_remote_flat( &
                         SIZE(lri_kp%centered_blocks(block_index)%inner_values, 1)* &
                         SIZE(lri_kp%centered_blocks(block_index)%inner_values, 3), &
                         SIZE(lri_kp%centered_blocks(block_index)%inner_values, 2)))
               DO iaux = 1, SIZE(lri_kp%centered_blocks(block_index)%inner_values, 3)
                  lri_kp%centered_blocks(block_index)%inner_remote_flat( &
                     (iaux - 1)*SIZE(lri_kp%centered_blocks(block_index)%inner_values, 1) + 1: &
                     iaux*SIZE(lri_kp%centered_blocks(block_index)%inner_values, 1), :) = &
                     lri_kp%centered_blocks(block_index)%inner_values(:, :, iaux)
               END DO
            END IF
         END IF
         CALL hfx_lri_lattice_block_from_index( &
            block_index, natom, ncell, iatom, jatom, cell_index)
         lri_kp%centered_atom_count(iatom) = lri_kp%centered_atom_count(iatom) + 1
         lri_kp%centered_remote_count(jatom) = lri_kp%centered_remote_count(jatom) + 1
         stored_blocks = stored_blocks + 1
         stored_reals = stored_reals + &
                        INT(SIZE(lri_kp%centered_blocks(block_index)%values), int_8) + &
                        INT(SIZE(lri_kp%centered_blocks(block_index)%remote_flat), int_8)
         IF (ALLOCATED(lri_kp%centered_blocks(block_index)%inner_values)) THEN
            stored_reals = stored_reals + &
                           INT(SIZE(lri_kp%centered_blocks(block_index)%inner_values), int_8) + &
                           INT(SIZE(lri_kp%centered_blocks(block_index)%inner_remote_flat), int_8)
         END IF
      END DO

      ALLOCATE (lri_kp%centered_atom_offset(natom + 1), &
                lri_kp%centered_remote_offset(natom + 1), SOURCE=1)
      DO iatom = 1, natom
         lri_kp%centered_atom_offset(iatom + 1) = lri_kp%centered_atom_offset(iatom) + &
                                                  lri_kp%centered_atom_count(iatom)
         lri_kp%centered_remote_offset(iatom + 1) = lri_kp%centered_remote_offset(iatom) + &
                                                    lri_kp%centered_remote_count(iatom)
      END DO
      nactive = lri_kp%centered_atom_offset(natom + 1) - 1
      ALLOCATE (lri_kp%centered_atom_refs(nactive), lri_kp%centered_remote_refs(nactive), SOURCE=0)
      ALLOCATE (cursor(natom))
      cursor(:) = lri_kp%centered_atom_offset(1:natom)
      DO block_index = 1, SIZE(lri_kp%centered_blocks)
         IF (.NOT. ALLOCATED(lri_kp%centered_blocks(block_index)%values)) CYCLE
         CALL hfx_lri_lattice_block_from_index( &
            block_index, natom, ncell, iatom, jatom, cell_index)
         lri_kp%centered_atom_refs(cursor(iatom)) = block_index
         cursor(iatom) = cursor(iatom) + 1
      END DO
      cursor(:) = lri_kp%centered_remote_offset(1:natom)
      DO block_index = 1, SIZE(lri_kp%centered_blocks)
         IF (.NOT. ALLOCATED(lri_kp%centered_blocks(block_index)%values)) CYCLE
         CALL hfx_lri_lattice_block_from_index( &
            block_index, natom, ncell, iatom, jatom, cell_index)
         lri_kp%centered_remote_refs(cursor(jatom)) = block_index
         cursor(jatom) = cursor(jatom) + 1
      END DO
      DEALLOCATE (cursor)

   END SUBROUTINE hfx_lri_centered_build

! **************************************************************************************************
!> \brief Adds a physical metric block to the untwisted lattice quotient.
!> \param lri_kp ...
!> \param iatom ...
!> \param jatom ...
!> \param cell ...
!> \param values ...
!> \param transposed ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_metric_lattice_add(lri_kp, iatom, jatom, cell, values, transposed)

      TYPE(hfx_lri_kp_type), INTENT(INOUT)               :: lri_kp
      INTEGER, INTENT(IN)                                :: iatom, jatom
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: values
      LOGICAL, INTENT(IN)                                :: transposed

      INTEGER                                            :: block_index, cell_index
      INTEGER, DIMENSION(3)                              :: reduced

      CALL hfx_lri_lattice_reduce(cell, lri_kp%lattice_grid, [1, 1, 1], reduced)
      cell_index = hfx_lri_lattice_cell_index(reduced, lri_kp%lattice_grid)
      block_index = hfx_lri_lattice_block_index( &
                    iatom, jatom, cell_index, lri_kp%lattice_natom, lri_kp%lattice_ncell)
      IF (.NOT. ALLOCATED(lri_kp%lattice_metric(block_index)%values)) THEN
         ALLOCATE (lri_kp%lattice_metric(block_index)%values( &
                   lri_kp%atom_naux(iatom), lri_kp%atom_naux(jatom)), SOURCE=0.0_dp)
      END IF
      IF (transposed) THEN
         IF (ANY(SHAPE(lri_kp%lattice_metric(block_index)%values) /= &
                 [SIZE(values, 2), SIZE(values, 1)])) THEN
            CPABORT("Inconsistent transposed metric dimensions in LRI loop3")
         END IF
         lri_kp%lattice_metric(block_index)%values(:, :) = &
            lri_kp%lattice_metric(block_index)%values + TRANSPOSE(values)
      ELSE
         IF (ANY(SHAPE(lri_kp%lattice_metric(block_index)%values) /= SHAPE(values))) THEN
            CPABORT("Inconsistent metric dimensions in LRI loop3")
         END IF
         lri_kp%lattice_metric(block_index)%values(:, :) = &
            lri_kp%lattice_metric(block_index)%values + values
      END IF

   END SUBROUTINE hfx_lri_metric_lattice_add

! **************************************************************************************************
!> \brief Folds one physical two-center metric into every matching atom/cell quotient block.
!> \param lri_kp ...
!> \param metric_entry ...
!> \param particle_set ...
!> \param cell ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_metric_lattice_add_entry(lri_kp, metric_entry, particle_set, cell)

      TYPE(hfx_lri_kp_type), INTENT(INOUT)               :: lri_kp
      TYPE(hfx_lri_metric_entry_type), INTENT(IN)        :: metric_entry
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell

      REAL(dp), PARAMETER                                :: cell_tolerance = 2.0E-8_dp

      INTEGER                                            :: iatom, jatom, nmapped
      INTEGER, DIMENSION(3)                              :: lattice_cell
      REAL(dp), DIMENSION(3)                             :: atom_delta, scaled

      IF (metric_entry%zero .OR. .NOT. ALLOCATED(metric_entry%values)) RETURN
      nmapped = 0
      DO iatom = 1, lri_kp%lattice_natom
         DO jatom = 1, lri_kp%lattice_natom
            atom_delta(:) = particle_set(jatom)%r(:) - particle_set(iatom)%r(:)
            IF (lri_kp%atom_kind(iatom) == metric_entry%kind_left .AND. &
                lri_kp%atom_kind(jatom) == metric_entry%kind_right) THEN
               scaled(:) = MATMUL(cell%h_inv, metric_entry%displacement - atom_delta)
               lattice_cell(:) = NINT(scaled(:))
               IF (MAXVAL(ABS(scaled - REAL(lattice_cell, dp))) < cell_tolerance) THEN
                  CALL hfx_lri_metric_lattice_add( &
                     lri_kp, iatom, jatom, lattice_cell, metric_entry%values, .FALSE.)
                  nmapped = nmapped + 1
               END IF
            END IF
            IF (MAXVAL(ABS(metric_entry%displacement)) <= lri_metric_cache_position_tol .AND. &
                metric_entry%kind_left == metric_entry%kind_right) CYCLE
            IF (lri_kp%atom_kind(iatom) == metric_entry%kind_right .AND. &
                lri_kp%atom_kind(jatom) == metric_entry%kind_left) THEN
               scaled(:) = MATMUL(cell%h_inv, -metric_entry%displacement - atom_delta)
               lattice_cell(:) = NINT(scaled(:))
               IF (MAXVAL(ABS(scaled - REAL(lattice_cell, dp))) < cell_tolerance) THEN
                  CALL hfx_lri_metric_lattice_add( &
                     lri_kp, iatom, jatom, lattice_cell, metric_entry%values, .TRUE.)
                  nmapped = nmapped + 1
               END IF
            END IF
         END DO
      END DO
      IF (nmapped == 0) CPABORT("Physical LRI metric key does not map to the lattice quotient")

   END SUBROUTINE hfx_lri_metric_lattice_add_entry

! **************************************************************************************************
!> \brief Computes physical metric blocks in bounded batches and immediately folds them.
!> \param ri_data ...
!> \param pair_lookup ...
!> \param density_cache ...
!> \param target_pairs ...
!> \param basis_set_RI ...
!> \param ri_transform_cache ...
!> \param particle_set ...
!> \param cell ...
!> \param maxl_ri ...
!> \param filter_eps ...
!> \param unique_entries ...
!> \param evaluated_entries ...
!> \param screened_entries ...
!> \param batch_mib ...
!> \param stored_blocks ...
!> \param stored_reals ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_metric_lattice_build( &
      ri_data, pair_lookup, density_cache, target_pairs, basis_set_RI, ri_transform_cache, &
      particle_set, cell, maxl_ri, filter_eps, unique_entries, evaluated_entries, &
      screened_entries, batch_mib, stored_blocks, stored_reals)

      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(lri_pair_lookup_type), INTENT(IN)             :: pair_lookup
      TYPE(lri_density_cache_type), INTENT(IN)           :: density_cache
      TYPE(hfx_lri_pair_block_type), DIMENSION(:), &
         INTENT(IN)                                      :: target_pairs
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_RI
      TYPE(lri_basis_transform_type), DIMENSION(:), &
         INTENT(IN)                                      :: ri_transform_cache
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: maxl_ri
      REAL(dp), INTENT(IN)                               :: filter_eps
      INTEGER, INTENT(OUT)                               :: unique_entries, evaluated_entries, &
                                                            screened_entries
      REAL(dp), INTENT(OUT)                              :: batch_mib
      INTEGER, INTENT(OUT)                               :: stored_blocks
      INTEGER(int_8), INTENT(OUT)                        :: stored_reals

      INTEGER                                            :: batch_first, batch_last, batch_size, &
                                                            block_index, ientry, nthreads
      INTEGER(int_8)                                     :: batch_peak_reals, batch_reals
      TYPE(coulomb_integral_context_type)                :: integral_context
      TYPE(hfx_lri_metric_entry_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: compact_entries
      TYPE(hfx_lri_metric_table_type)                    :: metric_table

      IF (.NOT. ASSOCIATED(ri_data%lri_kp) .OR. .NOT. ASSOCIATED(particle_set) .OR. &
          .NOT. ASSOCIATED(cell)) THEN
         CPABORT("Missing state or geometry while streaming the LRI metric")
      END IF
      CALL hfx_lri_metric_table_collect( &
         metric_table, ri_data%lri_kp%pair_blocks, pair_lookup, &
         density_cache, target_pairs)
      unique_entries = metric_table%nentries
      IF (unique_entries <= 0) CPABORT("LRI-HFXk collected no metric-table entries")
      DEALLOCATE (metric_table%slot)
      metric_table%capacity = 0
      IF (SIZE(metric_table%entries) /= unique_entries) THEN
         ALLOCATE (compact_entries(unique_entries))
         compact_entries(:) = metric_table%entries(1:unique_entries)
         CALL MOVE_ALLOC(compact_entries, metric_table%entries)
      END IF
      ALLOCATE (ri_data%lri_kp%lattice_metric( &
                ri_data%lri_kp%lattice_natom**2*ri_data%lri_kp%lattice_ncell))

      nthreads = 1
!$    nthreads = omp_get_max_threads()
      batch_size = MAX(1, 4*nthreads)
      batch_peak_reals = 0_int_8
      batch_reals = 0_int_8
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP          SHARED(basis_set_RI, batch_peak_reals, batch_reals, batch_size, cell, &
      !$OMP                 maxl_ri, metric_table, particle_set, ri_data, ri_transform_cache, &
      !$OMP                 unique_entries) &
      !$OMP          PRIVATE(batch_first, batch_last, ientry, integral_context)
      CALL coulomb_integral_init(integral_context, max_am_2c=maxl_ri)
      DO batch_first = 1, unique_entries, batch_size
         batch_last = MIN(unique_entries, batch_first + batch_size - 1)
         !$OMP DO SCHEDULE(DYNAMIC, 1)
         DO ientry = batch_first, batch_last
            CALL hfx_lri_metric_table_compute_entry( &
               metric_table%entries(ientry), basis_set_RI, &
               ri_transform_cache, ri_data, integral_context)
         END DO
         !$OMP END DO
         !$OMP SINGLE
         batch_reals = 0_int_8
         DO ientry = batch_first, batch_last
            IF (ALLOCATED(metric_table%entries(ientry)%values)) THEN
               batch_reals = batch_reals + INT(SIZE( &
                                               metric_table%entries(ientry)%values), int_8)
               CALL hfx_lri_metric_lattice_add_entry( &
                  ri_data%lri_kp, metric_table%entries(ientry), &
                  particle_set, cell)
               DEALLOCATE (metric_table%entries(ientry)%values)
            END IF
         END DO
         batch_peak_reals = MAX(batch_peak_reals, batch_reals)
         !$OMP END SINGLE
      END DO
      CALL coulomb_integral_cleanup(integral_context)
      !$OMP END PARALLEL

      evaluated_entries = 0
      screened_entries = 0
      DO ientry = 1, unique_entries
         IF (metric_table%entries(ientry)%zero) THEN
            screened_entries = screened_entries + 1
         ELSE
            evaluated_entries = evaluated_entries + 1
         END IF
      END DO
      batch_mib = REAL(batch_peak_reals, dp)*8.0_dp/1024.0_dp**2

      stored_blocks = 0
      stored_reals = 0_int_8
      DO block_index = 1, SIZE(ri_data%lri_kp%lattice_metric)
         IF (.NOT. ALLOCATED(ri_data%lri_kp%lattice_metric(block_index)%values)) CYCLE
         ri_data%lri_kp%lattice_metric(block_index)%norm2 = &
            SUM(ri_data%lri_kp%lattice_metric(block_index)%values**2)
         IF (ri_data%lri_kp%lattice_metric(block_index)%norm2 < filter_eps*filter_eps) THEN
            DEALLOCATE (ri_data%lri_kp%lattice_metric(block_index)%values)
            ri_data%lri_kp%lattice_metric(block_index)%norm2 = 0.0_dp
            CYCLE
         END IF
         stored_blocks = stored_blocks + 1
         stored_reals = stored_reals + &
                        INT(SIZE(ri_data%lri_kp%lattice_metric(block_index)%values), int_8)
      END DO

   END SUBROUTINE hfx_lri_metric_lattice_build

! **************************************************************************************************
!> \brief Stores one exact quasiperiodic density representative per atom pair and k-mesh cell.
!> \param lri_kp ...
!> \param density_cache ...
!> \param density_blocks ...
!> \param alias_residual ...
!> \param hermitian_residual ...
!> \param stored_reals ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_density_create(lri_kp, density_cache, density_blocks, &
                                           alias_residual, hermitian_residual, stored_reals)

      TYPE(hfx_lri_kp_type), INTENT(IN)                  :: lri_kp
      TYPE(lri_density_cache_type), INTENT(IN)           :: density_cache
      TYPE(lri_loop3_density_block_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(OUT)                       :: density_blocks
      REAL(dp), INTENT(OUT)                              :: alias_residual, hermitian_residual
      INTEGER(int_8), INTENT(OUT)                        :: stored_reals

      INTEGER                                            :: block_index, cell_index, iatom, iblock, &
                                                            jatom, reverse_cell_index, &
                                                            reverse_index
      INTEGER, DIMENSION(3)                              :: reduced, reverse_cell
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: seen
      REAL(dp)                                           :: block_scale, difference, factor, &
                                                            hermitian_scale, reverse_factor

      ALLOCATE (density_blocks(lri_kp%lattice_natom**2*lri_kp%lattice_ncell))
      ALLOCATE (seen(SIZE(density_blocks)), SOURCE=.FALSE.)
      alias_residual = 0.0_dp
      stored_reals = 0_int_8
      DO iblock = 1, density_cache%nblocks
         IF (.NOT. ALLOCATED(density_cache%blocks(iblock)%values)) CYCLE
         CALL hfx_lri_lattice_reduce( &
            density_cache%blocks(iblock)%cell, lri_kp%lattice_grid, &
            lri_kp%lattice_twist, reduced, factor)
         cell_index = hfx_lri_lattice_cell_index(reduced, lri_kp%lattice_grid)
         block_index = hfx_lri_lattice_block_index( &
                       density_cache%blocks(iblock)%iatom, density_cache%blocks(iblock)%jatom, &
                       cell_index, lri_kp%lattice_natom, lri_kp%lattice_ncell)
         IF (.NOT. seen(block_index)) THEN
            ALLOCATE (density_blocks(block_index)%values( &
                      SIZE(density_cache%blocks(iblock)%values, 1), &
                      SIZE(density_cache%blocks(iblock)%values, 2)))
            density_blocks(block_index)%values(:, :) = &
               factor*density_cache%blocks(iblock)%values
            seen(block_index) = .TRUE.
         ELSE
            IF (ANY(SHAPE(density_blocks(block_index)%values) /= &
                    SHAPE(density_cache%blocks(iblock)%values))) THEN
               CPABORT("Inconsistent aliased density shape in LRI loop3")
            END IF
            difference = MAXVAL(ABS(density_blocks(block_index)%values - &
                                    factor*density_cache%blocks(iblock)%values))
            block_scale = MAX(1.0E-14_dp, &
                              MAXVAL(ABS(density_blocks(block_index)%values)), &
                              MAXVAL(ABS(density_cache%blocks(iblock)%values)))
            alias_residual = MAX(alias_residual, difference/block_scale)
         END IF
      END DO
      DEALLOCATE (seen)

      hermitian_residual = 0.0_dp
      DO block_index = 1, SIZE(density_blocks)
         IF (.NOT. ALLOCATED(density_blocks(block_index)%values)) CYCLE
         density_blocks(block_index)%norm2 = SUM(density_blocks(block_index)%values**2)
         stored_reals = stored_reals + INT(SIZE(density_blocks(block_index)%values), int_8)
         CALL hfx_lri_lattice_block_from_index( &
            block_index, lri_kp%lattice_natom, lri_kp%lattice_ncell, &
            iatom, jatom, cell_index)
         CALL hfx_lri_lattice_cell_from_index(cell_index, lri_kp%lattice_grid, reduced)
         CALL hfx_lri_lattice_reduce( &
            -reduced, lri_kp%lattice_grid, lri_kp%lattice_twist, &
            reverse_cell, reverse_factor)
         reverse_cell_index = hfx_lri_lattice_cell_index(reverse_cell, lri_kp%lattice_grid)
         reverse_index = hfx_lri_lattice_block_index( &
                         jatom, iatom, reverse_cell_index, &
                         lri_kp%lattice_natom, lri_kp%lattice_ncell)
         IF (.NOT. ALLOCATED(density_blocks(reverse_index)%values)) THEN
            CPABORT("Missing Hermitian density partner in LRI loop3")
         END IF
         difference = MAXVAL(ABS(density_blocks(block_index)%values - &
                                 reverse_factor*TRANSPOSE(density_blocks(reverse_index)%values)))
         hermitian_scale = MAX(1.0E-14_dp, &
                               MAXVAL(ABS(density_blocks(block_index)%values)), &
                               MAXVAL(ABS(density_blocks(reverse_index)%values)))
         hermitian_residual = MAX(hermitian_residual, difference/hermitian_scale)
      END DO
      IF (alias_residual > 1.0E-9_dp) THEN
         CPABORT("Density aliases violate the LRI loop3 k-mesh boundary condition")
      END IF
      IF (hermitian_residual > 1.0E-9_dp) THEN
         CPABORT("Density quotient is not Hermitian in LRI loop3")
      END IF

   END SUBROUTINE hfx_lri_loop3_density_create

! **************************************************************************************************
!> \brief Adds one physical contribution to the exact twisted k-mesh quotient of K.
!> \param lri_kp ...
!> \param k_accum ...
!> \param iatom ...
!> \param jatom ...
!> \param cell ...
!> \param block ...
!> \param transposed ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_k_add(lri_kp, k_accum, iatom, jatom, cell, block, transposed)

      TYPE(hfx_lri_kp_type), INTENT(IN)                  :: lri_kp
      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: k_accum
      INTEGER, INTENT(IN)                                :: iatom, jatom
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: block
      LOGICAL, INTENT(IN), OPTIONAL                      :: transposed

      INTEGER                                            :: cell_index, first_i, first_j, last_i, &
                                                            last_j
      INTEGER, DIMENSION(3)                              :: reduced
      LOGICAL                                            :: use_transpose
      REAL(dp)                                           :: factor

      CALL hfx_lri_lattice_reduce( &
         cell, lri_kp%lattice_grid, lri_kp%lattice_twist, reduced, factor)
      cell_index = hfx_lri_lattice_cell_index(reduced, lri_kp%lattice_grid)
      first_i = lri_kp%atom_ao_offset(iatom)
      last_i = lri_kp%atom_ao_offset(iatom + 1) - 1
      first_j = lri_kp%atom_ao_offset(jatom)
      last_j = lri_kp%atom_ao_offset(jatom + 1) - 1
      use_transpose = .FALSE.
      IF (PRESENT(transposed)) use_transpose = transposed
      IF (use_transpose) THEN
         IF (ANY(SHAPE(block) /= [last_j - first_j + 1, last_i - first_i + 1])) THEN
            CPABORT("Inconsistent transposed output block in LRI loop3")
         END IF
         k_accum(first_i:last_i, first_j:last_j, cell_index) = &
            k_accum(first_i:last_i, first_j:last_j, cell_index) + factor*TRANSPOSE(block)
      ELSE
         IF (ANY(SHAPE(block) /= [last_i - first_i + 1, last_j - first_j + 1])) THEN
            CPABORT("Inconsistent output block in LRI loop3")
         END IF
         k_accum(first_i:last_i, first_j:last_j, cell_index) = &
            k_accum(first_i:last_i, first_j:last_j, cell_index) + factor*block
      END IF

   END SUBROUTINE hfx_lri_loop3_k_add

! **************************************************************************************************
!> \brief Builds S(J,B,u)=sum[L,b] D(J,L,u+b) C(B,L,b) on the finite k-mesh quotient.
!> \param lri_kp ...
!> \param density_blocks ...
!> \param group_index ...
!> \param s_values ...
!> \param active ...
!> \param used_terms ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_s_quotient_build(lri_kp, density_blocks, group_index, &
                                             s_values, active, used_terms)

      TYPE(hfx_lri_kp_type), INTENT(IN)                  :: lri_kp
      TYPE(lri_loop3_density_block_type), DIMENSION(:), &
         INTENT(IN)                                      :: density_blocks
      INTEGER, INTENT(IN)                                :: group_index
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: s_values
      LOGICAL, INTENT(OUT)                               :: active
      INTEGER(int_8), INTENT(OUT)                        :: used_terms

      INTEGER                                            :: c_cell_index, c_center, c_index, &
                                                            c_remote, density_cell_index, &
                                                            density_index, group_cell_index, &
                                                            iatom, iref, jatom
      INTEGER, DIMENSION(3)                              :: c_cell, density_cell, group_cell, reduced
      LOGICAL                                            :: first_s
      REAL(dp)                                           :: density_factor

      used_terms = 0_int_8
      CALL hfx_lri_lattice_block_from_index( &
         group_index, lri_kp%lattice_natom, lri_kp%lattice_ncell, &
         iatom, jatom, group_cell_index)
      CALL hfx_lri_lattice_cell_from_index( &
         group_cell_index, lri_kp%lattice_grid, group_cell)
      first_s = .TRUE.
      DO iref = lri_kp%centered_atom_offset(jatom), &
         lri_kp%centered_atom_offset(jatom + 1) - 1
         c_index = lri_kp%centered_atom_refs(iref)
         CALL hfx_lri_lattice_block_from_index( &
            c_index, lri_kp%lattice_natom, lri_kp%lattice_ncell, &
            c_center, c_remote, c_cell_index)
         CPASSERT(c_center == jatom)
         CALL hfx_lri_lattice_cell_from_index( &
            c_cell_index, lri_kp%lattice_grid, c_cell)
         density_cell(:) = group_cell(:) + c_cell(:)
         CALL hfx_lri_lattice_reduce( &
            density_cell, lri_kp%lattice_grid, lri_kp%lattice_twist, &
            reduced, density_factor)
         density_cell_index = hfx_lri_lattice_cell_index(reduced, lri_kp%lattice_grid)
         density_index = hfx_lri_lattice_block_index( &
                         iatom, c_remote, density_cell_index, &
                         lri_kp%lattice_natom, lri_kp%lattice_ncell)
         IF (.NOT. ALLOCATED(density_blocks(density_index)%values)) CYCLE
         CALL hfx_lri_loop3_density_c_add( &
            lri_kp%centered_blocks(c_index)%values, density_blocks(density_index)%values, &
            s_values, density_factor, first_s, .TRUE., &
            lri_kp%centered_blocks(c_index)%remote_flat)
         first_s = .FALSE.
         used_terms = used_terms + 1_int_8
      END DO
      active = .NOT. first_s

   END SUBROUTINE hfx_lri_loop3_s_quotient_build

! **************************************************************************************************
!> \brief Contracts the inner/outer and outer/outer terms for one quotient group.
!> \param lri_kp ...
!> \param s_values ...
!> \param group_index ...
!> \param k_accum ...
!> \param workspace ...
!> \param inner_outer ...
!> \param outer_outer ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_s_quotient_contract(lri_kp, s_values, group_index, k_accum, &
                                                workspace, inner_outer, outer_outer)

      TYPE(hfx_lri_kp_type), INTENT(IN)                  :: lri_kp
      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: s_values
      INTEGER, INTENT(IN)                                :: group_index
      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: k_accum
      TYPE(hfx_lri_contract_workspace_type), &
         INTENT(INOUT)                                   :: workspace
      INTEGER(int_8), INTENT(OUT)                        :: inner_outer, outer_outer

      INTEGER                                            :: c_cell_index, c_center, c_index, &
                                                            c_remote, group_cell_index, iatom, &
                                                            iref, jatom, metric_cell_index, &
                                                            metric_index
      INTEGER, DIMENSION(3)                              :: c_cell, group_cell, output_cell, reduced

      inner_outer = 0_int_8
      outer_outer = 0_int_8
      CALL hfx_lri_lattice_block_from_index( &
         group_index, lri_kp%lattice_natom, lri_kp%lattice_ncell, &
         iatom, jatom, group_cell_index)
      CALL hfx_lri_lattice_cell_from_index( &
         group_cell_index, lri_kp%lattice_grid, group_cell)
      CALL hfx_lri_loop3_pack_aux( &
         s_values, workspace%loop3_s_flat)

      IF (ALLOCATED(lri_kp%lattice_metric(group_index)%values)) THEN
         CALL hfx_lri_loop3_aux_transform( &
            s_values, lri_kp%lattice_metric(group_index)%values, &
            workspace%loop3_u, 1.0_dp, .TRUE.)
         CALL hfx_lri_loop3_pack_aux( &
            workspace%loop3_u, &
            workspace%loop3_u_flat)
         DO iref = lri_kp%centered_atom_offset(iatom), &
            lri_kp%centered_atom_offset(iatom + 1) - 1
            c_index = lri_kp%centered_atom_refs(iref)
            IF (.NOT. lri_kp%centered_blocks(c_index)%inner_active) CYCLE
            CALL hfx_lri_lattice_block_from_index( &
               c_index, lri_kp%lattice_natom, lri_kp%lattice_ncell, &
               c_center, c_remote, c_cell_index)
            CALL hfx_lri_lattice_cell_from_index( &
               c_cell_index, lri_kp%lattice_grid, c_cell)
            output_cell(:) = group_cell(:) - c_cell(:)
            IF (ALLOCATED(lri_kp%centered_blocks(c_index)%inner_values)) THEN
               CALL hfx_lri_loop3_remote_contract( &
                  lri_kp%centered_blocks(c_index)%inner_remote_flat, &
                  workspace%loop3_u_flat, workspace%loop3_block)
            ELSE
               CALL hfx_lri_loop3_remote_contract( &
                  lri_kp%centered_blocks(c_index)%remote_flat, &
                  workspace%loop3_u_flat, workspace%loop3_block)
            END IF
            CALL hfx_lri_loop3_k_add( &
               lri_kp, k_accum, c_remote, jatom, output_cell, &
               workspace%loop3_block)
            CALL hfx_lri_loop3_k_add( &
               lri_kp, k_accum, jatom, c_remote, -output_cell, &
               workspace%loop3_block, transposed=.TRUE.)
            inner_outer = inner_outer + 1_int_8
         END DO
      END IF

      DO iref = lri_kp%centered_remote_offset(iatom), &
         lri_kp%centered_remote_offset(iatom + 1) - 1
         c_index = lri_kp%centered_remote_refs(iref)
         CALL hfx_lri_lattice_block_from_index( &
            c_index, lri_kp%lattice_natom, lri_kp%lattice_ncell, &
            c_center, c_remote, c_cell_index)
         CPASSERT(c_remote == iatom)
         CALL hfx_lri_lattice_cell_from_index( &
            c_cell_index, lri_kp%lattice_grid, c_cell)
         output_cell(:) = group_cell(:) + c_cell(:)
         CALL hfx_lri_lattice_reduce( &
            output_cell, lri_kp%lattice_grid, lri_kp%lattice_twist, reduced)
         metric_cell_index = hfx_lri_lattice_cell_index(reduced, lri_kp%lattice_grid)
         metric_index = hfx_lri_lattice_block_index( &
                        c_center, jatom, metric_cell_index, &
                        lri_kp%lattice_natom, lri_kp%lattice_ncell)
         IF (.NOT. ALLOCATED(lri_kp%lattice_metric(metric_index)%values)) CYCLE
         CALL hfx_lri_loop3_central_contract( &
            lri_kp%centered_blocks(c_index)%values, &
            lri_kp%lattice_metric(metric_index)%values, &
            workspace%loop3_s_flat, workspace%loop3_block, workspace)
         CALL hfx_lri_loop3_k_add( &
            lri_kp, k_accum, c_center, jatom, output_cell, &
            workspace%loop3_block)
         outer_outer = outer_outer + 1_int_8
      END DO

   END SUBROUTINE hfx_lri_loop3_s_quotient_contract

! **************************************************************************************************
!> \brief Contracts the both-inner term for one quotient group.
!> \param lri_kp ...
!> \param density_blocks ...
!> \param group_index ...
!> \param k_accum ...
!> \param workspace ...
!> \param t_terms ...
!> \param inner_inner ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_t_quotient_contract(lri_kp, density_blocks, group_index, k_accum, &
                                                workspace, t_terms, inner_inner)

      TYPE(hfx_lri_kp_type), INTENT(IN)                  :: lri_kp
      TYPE(lri_loop3_density_block_type), DIMENSION(:), &
         INTENT(IN)                                      :: density_blocks
      INTEGER, INTENT(IN)                                :: group_index
      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: k_accum
      TYPE(hfx_lri_contract_workspace_type), &
         INTENT(INOUT)                                   :: workspace
      INTEGER(int_8), INTENT(OUT)                        :: t_terms, inner_inner

      INTEGER :: c_cell_index, c_center, c_index, c_remote, density_cell_index, density_index, &
         group_cell_index, iatom, iref, jatom, metric_index
      INTEGER, DIMENSION(3)                              :: c_cell, density_cell, group_cell, &
                                                            output_cell, reduced
      LOGICAL                                            :: first_t
      REAL(dp)                                           :: density_factor

      t_terms = 0_int_8
      inner_inner = 0_int_8
      CALL hfx_lri_lattice_block_from_index( &
         group_index, lri_kp%lattice_natom, lri_kp%lattice_ncell, &
         iatom, jatom, group_cell_index)
      CALL hfx_lri_lattice_cell_from_index( &
         group_cell_index, lri_kp%lattice_grid, group_cell)
      first_t = .TRUE.
      DO iref = lri_kp%centered_remote_offset(jatom), &
         lri_kp%centered_remote_offset(jatom + 1) - 1
         c_index = lri_kp%centered_remote_refs(iref)
         IF (.NOT. lri_kp%centered_blocks(c_index)%inner_active) CYCLE
         CALL hfx_lri_lattice_block_from_index( &
            c_index, lri_kp%lattice_natom, lri_kp%lattice_ncell, &
            c_center, c_remote, c_cell_index)
         CPASSERT(c_remote == jatom)
         CALL hfx_lri_lattice_cell_from_index( &
            c_cell_index, lri_kp%lattice_grid, c_cell)
         density_cell(:) = group_cell(:) - c_cell(:)
         CALL hfx_lri_lattice_reduce( &
            density_cell, lri_kp%lattice_grid, lri_kp%lattice_twist, &
            reduced, density_factor)
         density_cell_index = hfx_lri_lattice_cell_index(reduced, lri_kp%lattice_grid)
         density_index = hfx_lri_lattice_block_index( &
                         iatom, c_center, density_cell_index, &
                         lri_kp%lattice_natom, lri_kp%lattice_ncell)
         metric_index = hfx_lri_lattice_block_index( &
                        iatom, c_center, density_cell_index, &
                        lri_kp%lattice_natom, lri_kp%lattice_ncell)
         IF (.NOT. ALLOCATED(density_blocks(density_index)%values) .OR. &
             .NOT. ALLOCATED(lri_kp%lattice_metric(metric_index)%values)) CYCLE
         IF (ALLOCATED(lri_kp%centered_blocks(c_index)%inner_values)) THEN
            CALL hfx_lri_loop3_density_c_add( &
               lri_kp%centered_blocks(c_index)%inner_values, &
               density_blocks(density_index)%values, workspace%loop3_w, &
               density_factor, .TRUE., .FALSE.)
         ELSE
            CALL hfx_lri_loop3_density_c_add( &
               lri_kp%centered_blocks(c_index)%values, &
               density_blocks(density_index)%values, workspace%loop3_w, &
               density_factor, .TRUE., .FALSE.)
         END IF
         CALL hfx_lri_loop3_aux_transform( &
            workspace%loop3_w, &
            lri_kp%lattice_metric(metric_index)%values, &
            workspace%loop3_t, 1.0_dp, first_t)
         first_t = .FALSE.
         t_terms = t_terms + 1_int_8
      END DO
      IF (first_t) RETURN
      CALL hfx_lri_loop3_pack_aux( &
         workspace%loop3_t, &
         workspace%loop3_t_flat)

      DO iref = lri_kp%centered_atom_offset(iatom), &
         lri_kp%centered_atom_offset(iatom + 1) - 1
         c_index = lri_kp%centered_atom_refs(iref)
         IF (.NOT. lri_kp%centered_blocks(c_index)%inner_active) CYCLE
         CALL hfx_lri_lattice_block_from_index( &
            c_index, lri_kp%lattice_natom, lri_kp%lattice_ncell, &
            c_center, c_remote, c_cell_index)
         CALL hfx_lri_lattice_cell_from_index( &
            c_cell_index, lri_kp%lattice_grid, c_cell)
         output_cell(:) = group_cell(:) - c_cell(:)
         IF (ALLOCATED(lri_kp%centered_blocks(c_index)%inner_values)) THEN
            CALL hfx_lri_loop3_remote_contract( &
               lri_kp%centered_blocks(c_index)%inner_remote_flat, &
               workspace%loop3_t_flat, workspace%loop3_block)
         ELSE
            CALL hfx_lri_loop3_remote_contract( &
               lri_kp%centered_blocks(c_index)%remote_flat, &
               workspace%loop3_t_flat, workspace%loop3_block)
         END IF
         CALL hfx_lri_loop3_k_add( &
            lri_kp, k_accum, c_remote, jatom, output_cell, &
            workspace%loop3_block)
         inner_inner = inner_inner + 1_int_8
      END DO

   END SUBROUTINE hfx_lri_loop3_t_quotient_contract

! **************************************************************************************************
!> \brief Checks and removes roundoff-level Hermiticity noise from quotient K blocks.
!> \param lri_kp ...
!> \param k_accum ...
!> \param hermitian_residual ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_hermitize(lri_kp, k_accum, hermitian_residual)

      TYPE(hfx_lri_kp_type), INTENT(IN)                  :: lri_kp
      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: k_accum
      REAL(dp), INTENT(OUT)                              :: hermitian_residual

      INTEGER                                            :: block_index, cell_index, first_i, &
                                                            first_j, iatom, jatom, last_i, last_j, &
                                                            reverse_cell_index, reverse_index
      INTEGER, DIMENSION(3)                              :: cell, reverse_cell
      REAL(dp)                                           :: difference, factor, scale
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: symmetric_block

      scale = MAX(1.0E-14_dp, MAXVAL(ABS(k_accum)))
      hermitian_residual = 0.0_dp
      DO block_index = 1, lri_kp%lattice_natom**2*lri_kp%lattice_ncell
         CALL hfx_lri_lattice_block_from_index( &
            block_index, lri_kp%lattice_natom, lri_kp%lattice_ncell, &
            iatom, jatom, cell_index)
         CALL hfx_lri_lattice_cell_from_index(cell_index, lri_kp%lattice_grid, cell)
         CALL hfx_lri_lattice_reduce( &
            -cell, lri_kp%lattice_grid, lri_kp%lattice_twist, reverse_cell, factor)
         reverse_cell_index = hfx_lri_lattice_cell_index(reverse_cell, lri_kp%lattice_grid)
         reverse_index = hfx_lri_lattice_block_index( &
                         jatom, iatom, reverse_cell_index, &
                         lri_kp%lattice_natom, lri_kp%lattice_ncell)
         IF (block_index > reverse_index) CYCLE
         first_i = lri_kp%atom_ao_offset(iatom)
         last_i = lri_kp%atom_ao_offset(iatom + 1) - 1
         first_j = lri_kp%atom_ao_offset(jatom)
         last_j = lri_kp%atom_ao_offset(jatom + 1) - 1
         difference = MAXVAL(ABS( &
                             k_accum(first_i:last_i, first_j:last_j, cell_index) - &
                             factor*TRANSPOSE( &
                             k_accum(first_j:last_j, first_i:last_i, reverse_cell_index))))
         hermitian_residual = MAX(hermitian_residual, difference/scale)
      END DO
      IF (hermitian_residual > 1.0E-8_dp) THEN
         CPABORT("LRI loop3 K contraction failed the Hermiticity check")
      END IF

      DO block_index = 1, lri_kp%lattice_natom**2*lri_kp%lattice_ncell
         CALL hfx_lri_lattice_block_from_index( &
            block_index, lri_kp%lattice_natom, lri_kp%lattice_ncell, &
            iatom, jatom, cell_index)
         CALL hfx_lri_lattice_cell_from_index(cell_index, lri_kp%lattice_grid, cell)
         CALL hfx_lri_lattice_reduce( &
            -cell, lri_kp%lattice_grid, lri_kp%lattice_twist, reverse_cell, factor)
         reverse_cell_index = hfx_lri_lattice_cell_index(reverse_cell, lri_kp%lattice_grid)
         reverse_index = hfx_lri_lattice_block_index( &
                         jatom, iatom, reverse_cell_index, &
                         lri_kp%lattice_natom, lri_kp%lattice_ncell)
         IF (block_index > reverse_index) CYCLE
         first_i = lri_kp%atom_ao_offset(iatom)
         last_i = lri_kp%atom_ao_offset(iatom + 1) - 1
         first_j = lri_kp%atom_ao_offset(jatom)
         last_j = lri_kp%atom_ao_offset(jatom + 1) - 1
         IF (ALLOCATED(symmetric_block)) DEALLOCATE (symmetric_block)
         ALLOCATE (symmetric_block(last_i - first_i + 1, last_j - first_j + 1))
         symmetric_block(:, :) = 0.5_dp*( &
                                 k_accum(first_i:last_i, first_j:last_j, cell_index) + &
                                 factor*TRANSPOSE( &
                                 k_accum(first_j:last_j, first_i:last_i, &
                                         reverse_cell_index)))
         k_accum(first_i:last_i, first_j:last_j, cell_index) = symmetric_block
         k_accum(first_j:last_j, first_i:last_i, reverse_cell_index) = &
            factor*TRANSPOSE(symmetric_block)
      END DO
      IF (ALLOCATED(symmetric_block)) DEALLOCATE (symmetric_block)

   END SUBROUTINE hfx_lri_loop3_hermitize

! **************************************************************************************************
!> \brief Stores one physical representative of each nonzero quotient K orbit.
!> \param lri_kp ...
!> \param target_pairs ...
!> \param k_accum ...
!> \param filter_eps ...
!> \param stored_blocks ...
!> \param target_orbits ...
!> \param filtered_blocks ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_store_k(lri_kp, target_pairs, k_accum, filter_eps, &
                                    stored_blocks, target_orbits, filtered_blocks)

      TYPE(hfx_lri_kp_type), INTENT(INOUT)               :: lri_kp
      TYPE(hfx_lri_pair_block_type), DIMENSION(:), &
         INTENT(IN)                                      :: target_pairs
      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: k_accum
      REAL(dp), INTENT(IN)                               :: filter_eps
      INTEGER, INTENT(OUT)                               :: stored_blocks, target_orbits, &
                                                            filtered_blocks

      INTEGER :: block_iatom, block_index, block_jatom, cell_index, first_i, first_j, ipair, &
         last_i, last_j, reverse_cell_index, reverse_index
      INTEGER, DIMENSION(3)                              :: reduced, reverse_cell
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: selected
      REAL(dp)                                           :: factor, k_norm
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: k_block

      IF (ALLOCATED(lri_kp%k_blocks)) THEN
         DO ipair = 1, SIZE(lri_kp%k_blocks)
            IF (ALLOCATED(lri_kp%k_blocks(ipair)%k_block)) THEN
               DEALLOCATE (lri_kp%k_blocks(ipair)%k_block)
            END IF
         END DO
         DEALLOCATE (lri_kp%k_blocks)
      END IF
      ALLOCATE (lri_kp%k_blocks(lri_kp%lattice_natom**2*lri_kp%lattice_ncell))
      ALLOCATE (selected(SIZE(lri_kp%k_blocks)), SOURCE=.FALSE.)
      stored_blocks = 0
      target_orbits = 0
      filtered_blocks = 0

      DO ipair = 1, SIZE(target_pairs)
         CALL hfx_lri_lattice_reduce( &
            target_pairs(ipair)%cell, lri_kp%lattice_grid, lri_kp%lattice_twist, &
            reduced, factor)
         cell_index = hfx_lri_lattice_cell_index(reduced, lri_kp%lattice_grid)
         block_index = hfx_lri_lattice_block_index( &
                       target_pairs(ipair)%iatom, target_pairs(ipair)%jatom, cell_index, &
                       lri_kp%lattice_natom, lri_kp%lattice_ncell)
         CALL hfx_lri_lattice_reduce( &
            -reduced, lri_kp%lattice_grid, lri_kp%lattice_twist, &
            reverse_cell)
         reverse_cell_index = hfx_lri_lattice_cell_index(reverse_cell, lri_kp%lattice_grid)
         reverse_index = hfx_lri_lattice_block_index( &
                         target_pairs(ipair)%jatom, target_pairs(ipair)%iatom, &
                         reverse_cell_index, lri_kp%lattice_natom, lri_kp%lattice_ncell)
         IF (selected(block_index) .OR. selected(reverse_index)) CYCLE
         selected(block_index) = .TRUE.
         selected(reverse_index) = .TRUE.
         target_orbits = target_orbits + 1

         first_i = lri_kp%atom_ao_offset(target_pairs(ipair)%iatom)
         last_i = lri_kp%atom_ao_offset(target_pairs(ipair)%iatom + 1) - 1
         first_j = lri_kp%atom_ao_offset(target_pairs(ipair)%jatom)
         last_j = lri_kp%atom_ao_offset(target_pairs(ipair)%jatom + 1) - 1
         ALLOCATE (k_block(last_i - first_i + 1, last_j - first_j + 1))
         k_block(:, :) = factor*k_accum(first_i:last_i, first_j:last_j, cell_index)
         k_norm = SQRT(SUM(k_block*k_block))
         IF (k_norm < filter_eps) THEN
            DEALLOCATE (k_block)
            filtered_blocks = filtered_blocks + 1
            CYCLE
         END IF
         stored_blocks = stored_blocks + 1
         CALL hfx_lri_store_k_block( &
            lri_kp%k_blocks(stored_blocks), target_pairs(ipair), k_block)
      END DO

      DO block_index = 1, SIZE(selected)
         IF (selected(block_index)) CYCLE
         CALL hfx_lri_lattice_block_from_index( &
            block_index, lri_kp%lattice_natom, lri_kp%lattice_ncell, &
            block_iatom, block_jatom, cell_index)
         first_i = lri_kp%atom_ao_offset(block_iatom)
         last_i = lri_kp%atom_ao_offset(block_iatom + 1) - 1
         first_j = lri_kp%atom_ao_offset(block_jatom)
         last_j = lri_kp%atom_ao_offset(block_jatom + 1) - 1
         IF (SQRT(SUM(k_accum(first_i:last_i, first_j:last_j, cell_index)**2)) >= filter_eps) THEN
            CPABORT("A nonzero LRI loop3 quotient block has no KS target representative")
         END IF
      END DO
      DEALLOCATE (selected)
      lri_kp%nk_blocks = stored_blocks

   END SUBROUTINE hfx_lri_loop3_store_k

! **************************************************************************************************
!> \brief Exact finite-k loop3 contraction using quotient C, V, D and completed-K folding.
!> \param ri_data ...
!> \param density_cache ...
!> \param target_pairs ...
!> \param density_alias_residual ...
!> \param density_hermitian_residual ...
!> \param k_hermitian_residual ...
!> \param density_mib ...
!> \param working_mib ...
!> \param s_terms ...
!> \param t_terms ...
!> \param inner_outer ...
!> \param outer_outer ...
!> \param inner_inner ...
!> \param stored_blocks ...
!> \param target_orbits ...
!> \param filtered_blocks ...
!> \param k_norm ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_quotient_contract( &
      ri_data, density_cache, target_pairs, density_alias_residual, density_hermitian_residual, &
      k_hermitian_residual, density_mib, working_mib, s_terms, t_terms, inner_outer, &
      outer_outer, inner_inner, stored_blocks, target_orbits, filtered_blocks, k_norm)

      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(lri_density_cache_type), INTENT(IN)           :: density_cache
      TYPE(hfx_lri_pair_block_type), DIMENSION(:), &
         INTENT(IN)                                      :: target_pairs
      REAL(dp), INTENT(OUT)                              :: density_alias_residual, &
                                                            density_hermitian_residual, &
                                                            k_hermitian_residual, density_mib, &
                                                            working_mib
      INTEGER(int_8), INTENT(OUT)                        :: s_terms, t_terms, inner_outer, &
                                                            outer_outer, inner_inner
      INTEGER, INTENT(OUT)                               :: stored_blocks, target_orbits, &
                                                            filtered_blocks
      REAL(dp), INTENT(OUT)                              :: k_norm

      INTEGER                                            :: density_cell_index, density_iatom, &
                                                            density_jatom, igroup, ithread, &
                                                            nblocks, nthreads
      INTEGER(int_8) :: density_reals, inner_inner_local, inner_outer_local, outer_outer_local, &
         s_terms_local, t_terms_local, working_reals
      LOGICAL                                            :: s_active
      REAL(dp)                                           :: exchange_trace
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: k_accum
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :, :)       :: k_thread
      TYPE(hfx_lri_contract_workspace_type)              :: workspace
      TYPE(lri_loop3_density_block_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: density_blocks

      CALL hfx_lri_loop3_density_create( &
         ri_data%lri_kp, density_cache, density_blocks, density_alias_residual, &
         density_hermitian_residual, density_reals)
      density_mib = REAL(density_reals, dp)*8.0_dp/1024.0_dp**2
      nblocks = ri_data%lri_kp%lattice_natom**2*ri_data%lri_kp%lattice_ncell

      nthreads = SIZE(ri_data%lri_kp%thread_workspaces)
      ALLOCATE (k_thread(ri_data%lri_kp%lattice_total_nao, &
                         ri_data%lri_kp%lattice_total_nao, &
                         ri_data%lri_kp%lattice_ncell, nthreads), SOURCE=0.0_dp)
      s_terms = 0_int_8
      t_terms = 0_int_8
      inner_outer = 0_int_8
      outer_outer = 0_int_8
      inner_inner = 0_int_8

      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP          SHARED(density_blocks, k_thread, nblocks, ri_data) &
      !$OMP          PRIVATE(workspace, igroup, inner_inner_local, inner_outer_local, &
      !$OMP                  ithread, outer_outer_local, s_active, s_terms_local, t_terms_local) &
      !$OMP          REDUCTION(+:inner_inner, inner_outer, outer_outer, s_terms, t_terms)
      ithread = 1
!$    ithread = omp_get_thread_num() + 1
      CALL hfx_lri_contract_workspace_move( &
         ri_data%lri_kp%thread_workspaces(ithread), workspace)

      !$OMP DO SCHEDULE(DYNAMIC, 1)
      DO igroup = 1, nblocks
         CALL hfx_lri_loop3_s_quotient_build( &
            ri_data%lri_kp, density_blocks, igroup, &
            workspace%loop3_s, s_active, s_terms_local)
         s_terms = s_terms + s_terms_local
         inner_outer_local = 0_int_8
         outer_outer_local = 0_int_8
         IF (s_active) THEN
            CALL hfx_lri_loop3_s_quotient_contract( &
               ri_data%lri_kp, workspace%loop3_s, igroup, &
               k_thread(:, :, :, ithread), workspace, &
               inner_outer_local, outer_outer_local)
         END IF
         inner_outer = inner_outer + inner_outer_local
         outer_outer = outer_outer + outer_outer_local
      END DO
      !$OMP END DO

      !$OMP DO SCHEDULE(DYNAMIC, 1)
      DO igroup = 1, nblocks
         CALL hfx_lri_loop3_t_quotient_contract( &
            ri_data%lri_kp, density_blocks, igroup, k_thread(:, :, :, ithread), &
            workspace, t_terms_local, inner_inner_local)
         t_terms = t_terms + t_terms_local
         inner_inner = inner_inner + inner_inner_local
      END DO
      !$OMP END DO

      CALL hfx_lri_contract_workspace_move( &
         workspace, ri_data%lri_kp%thread_workspaces(ithread))
      !$OMP END PARALLEL

      working_reals = INT(SIZE(k_thread), int_8) + &
                      INT(ri_data%lri_kp%lattice_total_nao, int_8)**2* &
                      INT(ri_data%lri_kp%lattice_ncell, int_8)
      DO ithread = 1, nthreads
         IF (ALLOCATED(ri_data%lri_kp%thread_workspaces(ithread)%loop3_s_flat)) THEN
            working_reals = working_reals + INT(SIZE( &
                                                ri_data%lri_kp%thread_workspaces(ithread)%loop3_s_flat), int_8)
         END IF
         IF (ALLOCATED(ri_data%lri_kp%thread_workspaces(ithread)%loop3_t_flat)) THEN
            working_reals = working_reals + INT(SIZE( &
                                                ri_data%lri_kp%thread_workspaces(ithread)%loop3_t_flat), int_8)
         END IF
         IF (ALLOCATED(ri_data%lri_kp%thread_workspaces(ithread)%loop3_u_flat)) THEN
            working_reals = working_reals + INT(SIZE( &
                                                ri_data%lri_kp%thread_workspaces(ithread)%loop3_u_flat), int_8)
         END IF
         IF (ALLOCATED(ri_data%lri_kp%thread_workspaces(ithread)%loop3_block)) THEN
            working_reals = working_reals + INT(SIZE( &
                                                ri_data%lri_kp%thread_workspaces(ithread)%loop3_block), int_8)
         END IF
         IF (ALLOCATED(ri_data%lri_kp%thread_workspaces(ithread)%loop3_s)) THEN
            working_reals = working_reals + INT(SIZE( &
                                                ri_data%lri_kp%thread_workspaces(ithread)%loop3_s), int_8)
         END IF
         IF (ALLOCATED(ri_data%lri_kp%thread_workspaces(ithread)%loop3_t)) THEN
            working_reals = working_reals + INT(SIZE( &
                                                ri_data%lri_kp%thread_workspaces(ithread)%loop3_t), int_8)
         END IF
         IF (ALLOCATED(ri_data%lri_kp%thread_workspaces(ithread)%loop3_u)) THEN
            working_reals = working_reals + INT(SIZE( &
                                                ri_data%lri_kp%thread_workspaces(ithread)%loop3_u), int_8)
         END IF
         IF (ALLOCATED(ri_data%lri_kp%thread_workspaces(ithread)%loop3_w)) THEN
            working_reals = working_reals + INT(SIZE( &
                                                ri_data%lri_kp%thread_workspaces(ithread)%loop3_w), int_8)
         END IF
      END DO
      working_mib = REAL(working_reals, dp)*8.0_dp/1024.0_dp**2
      ALLOCATE (k_accum(ri_data%lri_kp%lattice_total_nao, &
                        ri_data%lri_kp%lattice_total_nao, &
                        ri_data%lri_kp%lattice_ncell), SOURCE=0.0_dp)
      DO ithread = 1, nthreads
         k_accum(:, :, :) = k_accum + k_thread(:, :, :, ithread)
      END DO
      DEALLOCATE (k_thread)
      CALL hfx_lri_loop3_hermitize(ri_data%lri_kp, k_accum, k_hermitian_residual)
      k_norm = SQRT(SUM(k_accum*k_accum))

      exchange_trace = 0.0_dp
      DO igroup = 1, SIZE(density_blocks)
         IF (.NOT. ALLOCATED(density_blocks(igroup)%values)) CYCLE
         CALL hfx_lri_lattice_block_from_index( &
            igroup, ri_data%lri_kp%lattice_natom, ri_data%lri_kp%lattice_ncell, &
            density_iatom, density_jatom, density_cell_index)
         exchange_trace = exchange_trace + SUM( &
                          density_blocks(igroup)%values* &
                          k_accum( &
                          ri_data%lri_kp%atom_ao_offset(density_iatom): &
                          ri_data%lri_kp%atom_ao_offset(density_iatom + 1) - 1, &
                          ri_data%lri_kp%atom_ao_offset(density_jatom): &
                          ri_data%lri_kp%atom_ao_offset(density_jatom + 1) - 1, &
                          density_cell_index))
      END DO
      ri_data%lri_kp%exchange_trace_raw = exchange_trace
      CALL hfx_lri_loop3_store_k( &
         ri_data%lri_kp, target_pairs, k_accum, ri_data%filter_eps, &
         stored_blocks, target_orbits, filtered_blocks)

   END SUBROUTINE hfx_lri_loop3_quotient_contract

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

      INTEGER :: c_blocks_loop3, filtered_k_blocks, handle, ikind, ipair, max_k_blocks, maxl_ri, &
         maxl_ri_kind, metric_evaluated, metric_screened, metric_unique, natom, nk_blocks_built, &
         nkind, nthreads, v_blocks_loop3
      INTEGER(int_8) :: c_reals_loop3, inner_inner_loop3, inner_outer_loop3, outer_outer_loop3, &
         s_terms_loop3, t_terms_loop3, v_reals_loop3
      INTEGER, DIMENSION(3)                              :: grid, twist
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: need_geometry
      REAL(dp) :: density_alias_loop3, density_hermitian_loop3, density_mib_loop3, &
         exchange_trace_raw, k_blocks_norm, k_hermitian_loop3, metric_batch_mib, result_mib_loop3
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_AO, basis_set_RI
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(hfx_lri_pair_block_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: target_pairs
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(lri_basis_transform_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: ri_transform_cache
      TYPE(lri_density_cache_type)                       :: density_cache
      TYPE(lri_pair_lookup_type)                         :: pair_lookup
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (basis_set_AO, basis_set_RI, cell, cell_to_index, dft_control, kpoints, particle_set, &
               qs_kind_set, ri_basis, sab_nl)

      IF (.NOT. ASSOCIATED(ri_data%lri_kp)) CPABORT("LRI-HFXk state is not associated")
      IF (.NOT. ALLOCATED(ri_data%lri_kp%pair_blocks)) THEN
         CPABORT("LRI-HFXk fitting-pair blocks are not allocated")
      END IF
      IF (ri_data%lri_kp%npair_blocks <= 0) CPABORT("LRI-HFXk has no fitting-pair blocks")
      IF (.NOT. ASSOCIATED(rho_ao)) CPABORT("LRI-HFXk density matrix is not associated")

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, cell=cell, dft_control=dft_control, kpoints=kpoints, &
                      particle_set=particle_set, qs_kind_set=qs_kind_set, nkind=nkind, natom=natom)
      IF (ASSOCIATED(kpoints)) THEN
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index, sab_nl=sab_nl)
      END IF

      need_geometry = .NOT. ri_data%lri_kp%loop3_ready
      ALLOCATE (basis_set_AO(nkind))
      CALL basis_set_list_setup(basis_set_AO, ri_data%orb_basis_type, qs_kind_set)
      CALL hfx_lri_density_cache_create(rho_ao, cell_to_index, sab_nl, ri_data%filter_eps, density_cache)
      CALL hfx_lri_target_pairs_create(sab_nl, cell_to_index, SIZE(rho_ao, 2), &
                                       basis_set_AO, target_pairs)
      IF (SIZE(target_pairs) == 0) THEN
         CPABORT("LRI-HFXk found no writable target blocks in the k-point neighbor list")
      END IF

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

      IF (need_geometry) THEN
         ALLOCATE (basis_set_RI(nkind))
         CALL basis_set_list_setup(basis_set_RI, ri_data%ri_basis_type, qs_kind_set)
         ! The basis radii are mutable shared metadata. Use the LRI threshold while evaluating
         ! the streamed metric and restore the ordinary QS value before returning.
         DO ikind = 1, nkind
            ri_basis => basis_set_RI(ikind)%gto_basis_set
            IF (ASSOCIATED(ri_basis)) THEN
               CALL init_interaction_radii_orb_basis(ri_basis, ri_data%eps_pgf_orb)
            END IF
         END DO
         CALL hfx_lri_pair_lookup_create(pair_lookup, ri_data%lri_kp%pair_blocks, &
                                         ri_data%lri_kp%npair_blocks, natom)
         CALL hfx_lri_build_ri_transform_cache(basis_set_RI, ri_transform_cache)
         maxl_ri = 0
         DO ikind = 1, nkind
            ri_basis => basis_set_RI(ikind)%gto_basis_set
            IF (.NOT. ASSOCIATED(ri_basis)) CYCLE
            CALL get_gto_basis_set(ri_basis, maxl=maxl_ri_kind)
            maxl_ri = MAX(maxl_ri, maxl_ri_kind)
         END DO
         CALL init_md_ftable(nmax=2*maxl_ri)

         CALL hfx_lri_lattice_grid(kpoints, grid, twist)
         CALL hfx_lri_centered_build( &
            ri_data%lri_kp, natom, grid, twist, ri_data%filter_eps, &
            c_blocks_loop3, c_reals_loop3)

         IF (ri_data%unit_nr > 0) THEN
            WRITE (ri_data%unit_nr, "(/,T2,A)") "LRI-HFXk streamed metric setup"
            WRITE (ri_data%unit_nr, "(T4,A,T28,I10)") "target blocks", SIZE(target_pairs)
            WRITE (ri_data%unit_nr, "(T4,A,T28,I10)") "density blocks", density_cache%nblocks
            FLUSH (ri_data%unit_nr)
         END IF
         CALL hfx_lri_metric_lattice_build( &
            ri_data, pair_lookup, density_cache, target_pairs, basis_set_RI, &
            ri_transform_cache, particle_set, cell, maxl_ri, ri_data%filter_eps_2c, &
            metric_unique, metric_evaluated, metric_screened, metric_batch_mib, &
            v_blocks_loop3, v_reals_loop3)
         ri_data%lri_kp%loop3_ready = .TRUE.

         ! The quotient contains all information needed by subsequent SCF iterations.
         ! Release physical coefficient blocks and the small construction caches now.
         DO ipair = 1, ri_data%lri_kp%npair_blocks
            IF (ALLOCATED(ri_data%lri_kp%pair_blocks(ipair)%c_fit)) THEN
               DEALLOCATE (ri_data%lri_kp%pair_blocks(ipair)%c_fit)
            END IF
         END DO
         IF (ri_data%unit_nr > 0) THEN
            WRITE (ri_data%unit_nr, "(T4,A,T28,I10)") &
               "unique metric keys", metric_unique
            WRITE (ri_data%unit_nr, "(T4,A,T28,I10)") "evaluated metric keys", metric_evaluated
            WRITE (ri_data%unit_nr, "(T4,A,T28,I10)") "distance-screened", metric_screened
            WRITE (ri_data%unit_nr, "(T4,A,T28,F10.2)") "peak metric batch/MiB", metric_batch_mib
            WRITE (ri_data%unit_nr, "(/,T2,A)") "LRI-HFXk finite-k loop3 plan"
            WRITE (ri_data%unit_nr, "(T4,A,T28,3I6)") &
               "quotient grid", ri_data%lri_kp%lattice_grid
            WRITE (ri_data%unit_nr, "(T4,A,T28,3I6)") &
               "boundary twist", ri_data%lri_kp%lattice_twist
            WRITE (ri_data%unit_nr, "(T4,A,T28,I14)") &
               "quotient groups", ri_data%lri_kp%lattice_natom**2*ri_data%lri_kp%lattice_ncell
            WRITE (ri_data%unit_nr, "(T4,A,T28,I14)") "centered C blocks", c_blocks_loop3
            WRITE (ri_data%unit_nr, "(T4,A,T28,F14.2)") &
               "C plus packing/MiB", REAL(c_reals_loop3, dp)*8.0_dp/1024.0_dp**2
            WRITE (ri_data%unit_nr, "(T4,A,T28,I14)") "quotient V blocks", v_blocks_loop3
            WRITE (ri_data%unit_nr, "(T4,A,T28,F14.2)") &
               "quotient V/MiB", REAL(v_reals_loop3, dp)*8.0_dp/1024.0_dp**2
            FLUSH (ri_data%unit_nr)
         END IF
      END IF

      CALL hfx_lri_loop3_quotient_contract( &
         ri_data, density_cache, target_pairs, density_alias_loop3, density_hermitian_loop3, &
         k_hermitian_loop3, density_mib_loop3, result_mib_loop3, s_terms_loop3, t_terms_loop3, &
         inner_outer_loop3, outer_outer_loop3, inner_inner_loop3, nk_blocks_built, &
         max_k_blocks, filtered_k_blocks, k_blocks_norm)
      exchange_trace_raw = ri_data%lri_kp%exchange_trace_raw
      IF (ri_data%unit_nr > 0) THEN
         WRITE (ri_data%unit_nr, "(/,T2,A)") "LRI-HFXk loop3 contraction"
         WRITE (ri_data%unit_nr, "(T4,A,T28,I14)") "S density-C terms", s_terms_loop3
         WRITE (ri_data%unit_nr, "(T4,A,T28,I14)") "T density-C-V terms", t_terms_loop3
         WRITE (ri_data%unit_nr, "(T4,A,T28,I14)") "inner-outer", inner_outer_loop3
         WRITE (ri_data%unit_nr, "(T4,A,T28,I14)") "outer-outer", outer_outer_loop3
         WRITE (ri_data%unit_nr, "(T4,A,T28,I14)") "inner-inner", inner_inner_loop3
         WRITE (ri_data%unit_nr, "(T4,A,T28,I10,T41,I10)") &
            "stored/target orbits", nk_blocks_built, max_k_blocks
         WRITE (ri_data%unit_nr, "(T4,A,T28,I10)") "filtered K blocks", filtered_k_blocks
         WRITE (ri_data%unit_nr, "(T4,A,T28,F10.3,T41,F10.3)") &
            "density/working MiB", density_mib_loop3, result_mib_loop3
         WRITE (ri_data%unit_nr, "(T4,A,T28,ES10.3,T41,ES10.3)") &
            "density alias/Herm", density_alias_loop3, density_hermitian_loop3
         WRITE (ri_data%unit_nr, "(T4,A,T28,ES10.3)") "K Herm residual", k_hermitian_loop3
         WRITE (ri_data%unit_nr, "(T4,A,T28,ES10.3,T41,ES14.6)") &
            "K norm/Tr[D K]", k_blocks_norm, exchange_trace_raw
         FLUSH (ri_data%unit_nr)
      END IF

      IF (ALLOCATED(target_pairs)) DEALLOCATE (target_pairs)
      IF (need_geometry .AND. ASSOCIATED(dft_control)) THEN
         DO ikind = 1, nkind
            ri_basis => basis_set_RI(ikind)%gto_basis_set
            IF (ASSOCIATED(ri_basis)) THEN
               CALL init_interaction_radii_orb_basis( &
                  ri_basis, dft_control%qs_control%eps_pgf_orb)
            END IF
         END DO
      END IF
      IF (ASSOCIATED(basis_set_AO)) DEALLOCATE (basis_set_AO)
      IF (ASSOCIATED(basis_set_RI)) DEALLOCATE (basis_set_RI)
      CALL timestop(handle)

   END SUBROUTINE hfx_lri_build_rank_local_k_blocks

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
         islot = ientry
         use_transpose = .TRUE.
         is_zero = metric_cache%entries(ientry)%zero
         RETURN
      END DO

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

END MODULE hfx_lri_kp
