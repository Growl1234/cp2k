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
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_block_p,&
                                              dbcsr_get_stored_coordinates,&
                                              dbcsr_p_type
   USE hfx_lri,                         ONLY: hfx_lri_build_exchange,&
                                              hfx_lri_collect_pairs,&
                                              hfx_lri_density_cache_add,&
                                              hfx_lri_prepare_pairs,&
                                              lri_density_cache_type
   USE hfx_lri_parallel,                ONLY: hfx_lri_block_type,&
                                              hfx_lri_exchange_blocks
   USE hfx_lri_types,                   ONLY: hfx_lri_pair_block_type,&
                                              hfx_lri_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hfx_lri_update_ks_kp

CONTAINS

! **************************************************************************************************
!> \brief Updates the k-point KS matrices with the pair-local LRI HFX contribution.
!> \param qs_env ...
!> \param lri ...
!> \param ks_matrix ...
!> \param ehfx ...
!> \param rho_ao ...
!> \param geometry_did_change ...
!> \param hf_fraction ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_update_ks_kp(qs_env, lri, ks_matrix, ehfx, rho_ao, &
                                   geometry_did_change, hf_fraction)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_lri_type), INTENT(INOUT)                  :: lri
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      REAL(KIND=dp), INTENT(OUT)                         :: ehfx
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      LOGICAL, INTENT(IN)                                :: geometry_did_change
      REAL(KIND=dp), INTENT(IN)                          :: hf_fraction

      INTEGER                                            :: nkind
      INTEGER, DIMENSION(3)                              :: grid, twist
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: stable_images
      REAL(dp)                                           :: exchange_trace
      TYPE(hfx_lri_pair_block_type), ALLOCATABLE         :: target_pairs(:)
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(lri_density_cache_type)                       :: density_cache
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl

      ehfx = 0.0_dp
      CALL get_qs_env(qs_env, para_env=para_env, kpoints=kpoints, nkind=nkind)
      CALL hfx_lri_lattice_grid(kpoints, grid, twist)
      CALL get_kpoint_info(kpoints, symmetry=stable_images, cell_to_index=cell_to_index, sab_nl=sab_nl)
      IF (.NOT. ASSOCIATED(sab_nl)) CPABORT("LRI-HFXk could not obtain its matrix neighbor list")
      IF (geometry_did_change .OR. .NOT. ALLOCATED(lri%lattice)) THEN
         CALL hfx_lri_prepare_pairs(qs_env, lri, grid, twist, stable_images)
      END IF

      CALL hfx_lri_target_pairs_create(sab_nl, cell_to_index, SIZE(rho_ao, 2), nkind, target_pairs)
      CALL hfx_lri_collect_pairs(target_pairs, para_env)
      CALL hfx_lri_density_cache_create(rho_ao, cell_to_index, target_pairs, para_env, density_cache)
      CALL hfx_lri_build_exchange(qs_env, lri, grid, twist, stable_images, target_pairs, density_cache, exchange_trace)
      ! For an RKS spin-summed density, F_x = -hf_fraction*K_raw/2 and
      ! E_x = Tr[D F_x]/2 = -hf_fraction*Tr[D K_raw]/4.
      ehfx = -0.25_dp*hf_fraction*exchange_trace
      CALL hfx_lri_write_rank_local_k_blocks(cell_to_index, lri, ks_matrix, hf_fraction)

   END SUBROUTINE hfx_lri_update_ks_kp

! **************************************************************************************************
!> \brief Resolve a physical cell, returning zero for an absent or excluded image.
!> \param cell ...
!> \param cell_to_index null for a single-image representation
!> \param nimages number of active real-space matrices
!> \return ...
! **************************************************************************************************
   PURE INTEGER FUNCTION hfx_lri_image_index(cell, cell_to_index, nimages) RESULT(image)
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell
      INTEGER, DIMENSION(:, :, :), INTENT(IN), POINTER   :: cell_to_index
      INTEGER, INTENT(IN)                                :: nimages

      image = 1
      IF (ASSOCIATED(cell_to_index)) THEN
         image = 0
         IF (ANY(cell < LBOUND(cell_to_index)) .OR. ANY(cell > UBOUND(cell_to_index))) RETURN
         image = cell_to_index(cell(1), cell(2), cell(3))
      END IF
      IF (image < 1 .OR. image > nimages) image = 0
   END FUNCTION hfx_lri_image_index

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
!> \brief Reconstructs and caches only the logical density blocks present in CP2K's symmetric
!>        real-space representation. This is the sparse equivalent of get_pmat_images().
!> \param rho_ao ...
!> \param cell_to_index ...
!> \param targets globally collected physical target blocks
!> \param para_env ...
!> \param density_cache ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_density_cache_create(rho_ao, cell_to_index, targets, para_env, density_cache)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(hfx_lri_pair_block_type), DIMENSION(:), &
         INTENT(IN)                                      :: targets
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      TYPE(lri_density_cache_type), INTENT(OUT)          :: density_cache

      INTEGER                                            :: iatom, icol, ikind, image, image_rev, &
                                                            ipair, irow, jatom, jkind, nentry, &
                                                            owner, timing_handle
      INTEGER, DIMENSION(3)                              :: cell
      LOGICAL                                            :: found, source_transposed
      REAL(dp)                                           :: factor
      REAL(dp), DIMENSION(3)                             :: rab
      REAL(dp), DIMENSION(:, :), POINTER                 :: rho_block

      CALL timeset('hfx_lri_density_cache_create', timing_handle)
      NULLIFY (rho_block)
      CPASSERT(ASSOCIATED(rho_ao))

      ! The symmetric list contains at most one stored representative per physical pair. Each
      ! representative contributes its direct orientation and, where present, its Hermitian image.
      nentry = SIZE(targets)

      ALLOCATE (density_cache%blocks(2*nentry))
      density_cache%capacity = MAX(17, 4*nentry + 1)
      ALLOCATE (density_cache%slot(density_cache%capacity), SOURCE=0)
      density_cache%nblocks = 0

      DO ipair = 1, SIZE(targets)
         iatom = targets(ipair)%iatom
         jatom = targets(ipair)%jatom
         ikind = targets(ipair)%ikind
         jkind = targets(ipair)%jkind
         cell = targets(ipair)%cell
         rab = targets(ipair)%rab
         image = hfx_lri_image_index(cell, cell_to_index, SIZE(rho_ao, 2))
         IF (image == 0) CYCLE
         image_rev = hfx_lri_image_index(-cell, cell_to_index, SIZE(rho_ao, 2))
         factor = MERGE(0.5_dp, 1.0_dp, iatom == jatom .AND. image_rev > 0)

         irow = iatom
         icol = jatom
         source_transposed = .FALSE.
         IF (iatom > jatom) THEN
            irow = jatom
            icol = iatom
            source_transposed = .TRUE.
         END IF
         CALL dbcsr_get_stored_coordinates(rho_ao(1, image)%matrix, irow, icol, owner)
         IF (owner /= para_env%mepos) CYCLE
         NULLIFY (rho_block)
         found = .FALSE.
         CALL dbcsr_get_block_p(matrix=rho_ao(1, image)%matrix, row=irow, col=icol, &
                                block=rho_block, found=found)
         IF (.NOT. found) CYCLE

         CALL hfx_lri_density_cache_add(density_cache, iatom, jatom, ikind, jkind, cell, rab, &
                                        rho_block, source_transposed, factor)
         IF (image_rev > 0) THEN
            CALL hfx_lri_density_cache_add(density_cache, jatom, iatom, jkind, ikind, -cell, -rab, &
                                           rho_block,.NOT. source_transposed, factor)
         END IF
      END DO

      CALL timestop(timing_handle)
   END SUBROUTINE hfx_lri_density_cache_create

! **************************************************************************************************
!> \brief Builds the exact symmetric set of AO blocks used by the k-point KS matrices. Fitting-pair
!>        locality must not be used to truncate the output K matrix: exchange can couple two local
!>        fitting pairs whose outer AO centers do not overlap.
!> \param sab_nl ...
!> \param cell_to_index ...
!> \param nimages ...
!> \param nkind ...
!> \param target_pairs ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_target_pairs_create(sab_nl, cell_to_index, nimages, nkind, target_pairs)

      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      INTEGER, INTENT(IN)                                :: nimages, nkind
      TYPE(hfx_lri_pair_block_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(OUT)                       :: target_pairs

      INTEGER                                            :: iatom, ikind, image, ipair, jatom, &
                                                            jkind, mepos, ntarget
      INTEGER, DIMENSION(3)                              :: cell
      REAL(dp), DIMENSION(3)                             :: rab
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      NULLIFY (nl_iterator)
      CPASSERT(ASSOCIATED(sab_nl))
      CPASSERT(nimages > 0)
      CPASSERT(nkind > 0)

      ntarget = 0
      CALL neighbor_list_iterator_create(nl_iterator, sab_nl)
      mepos = 0
      DO WHILE (neighbor_list_iterate(nl_iterator, mepos) == 0)
         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, cell=cell)
         IF (ikind < 1 .OR. ikind > nkind .OR. jkind < 1 .OR. jkind > nkind) CYCLE
         image = hfx_lri_image_index(cell, cell_to_index, nimages)
         IF (image == 0) CYCLE
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
         IF (ikind < 1 .OR. ikind > nkind .OR. jkind < 1 .OR. jkind > nkind) CYCLE
         image = hfx_lri_image_index(cell, cell_to_index, nimages)
         IF (image == 0) CYCLE

         ipair = ipair + 1
         target_pairs(ipair)%iatom = iatom
         target_pairs(ipair)%jatom = jatom
         target_pairs(ipair)%ikind = ikind
         target_pairs(ipair)%jkind = jkind
         target_pairs(ipair)%cell(:) = cell(:)
         target_pairs(ipair)%rab(:) = rab(:)
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)
      CPASSERT(ipair == ntarget)

   END SUBROUTINE hfx_lri_target_pairs_create

! **************************************************************************************************
!> \brief Adds the completed pair-local exchange blocks to the k-point KS matrices.
!> \param cell_to_index ...
!> \param lri ...
!> \param ks_matrix ...
!> \param hf_fraction ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_write_rank_local_k_blocks(cell_to_index, lri, ks_matrix, hf_fraction)
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(hfx_lri_type), INTENT(IN)                     :: lri
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      REAL(dp), INTENT(IN)                               :: hf_fraction

      INTEGER                                            :: icol, img, ipair, irow, timing_handle
      INTEGER, ALLOCATABLE                               :: images(:), owners(:)
      LOGICAL                                            :: found, use_transpose
      REAL(dp)                                           :: scale
      REAL(dp), DIMENSION(:, :), POINTER                 :: ks_block
      TYPE(hfx_lri_block_type), ALLOCATABLE              :: received(:)

      CALL timeset('hfx_lri_write_rank_local_k_blocks', timing_handle)
      NULLIFY (ks_block)
      CPASSERT(ALLOCATED(lri%lattice))
      CPASSERT(ASSOCIATED(ks_matrix))
      CPASSERT(ASSOCIATED(cell_to_index))
      CPASSERT(ALLOCATED(lri%lattice%k_blocks))
      CPASSERT(SIZE(lri%lattice%k_blocks) == SIZE(lri%lattice%target_pairs))
      ALLOCATE (received(SIZE(lri%lattice%k_blocks)))
      ALLOCATE (owners(SIZE(received)), images(SIZE(received)), SOURCE=0)
      DO ipair = 1, SIZE(received)
         img = hfx_lri_image_index(lri%lattice%target_pairs(ipair)%cell, cell_to_index, SIZE(ks_matrix, 2))
         CPASSERT(img >= 1 .AND. img <= SIZE(ks_matrix, 2))
         images(ipair) = img
         irow = MIN(lri%lattice%target_pairs(ipair)%iatom, lri%lattice%target_pairs(ipair)%jatom)
         icol = MAX(lri%lattice%target_pairs(ipair)%iatom, lri%lattice%target_pairs(ipair)%jatom)
         CALL dbcsr_get_stored_coordinates(ks_matrix(1, img)%matrix, irow, icol, owners(ipair))
      END DO
      CALL hfx_lri_exchange_blocks(lri%lattice%k_blocks, owners, received, lri%lattice%parallel%para_env)
      ! One physical KS representative per quotient orbit. No MPI occurs inside this loop.
      ! The input density is spin-summed RKS, so its same-spin density is rho/2.
      scale = -0.5_dp*hf_fraction
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(lri, ks_matrix, received, images, scale) &
      !$OMP PRIVATE(ipair, img, irow, icol, use_transpose, ks_block, found) SCHEDULE(DYNAMIC, 1)
      DO ipair = 1, SIZE(received)
         IF (.NOT. ALLOCATED(received(ipair)%values)) CYCLE
         img = images(ipair)
         irow = MIN(lri%lattice%target_pairs(ipair)%iatom, lri%lattice%target_pairs(ipair)%jatom)
         icol = MAX(lri%lattice%target_pairs(ipair)%iatom, lri%lattice%target_pairs(ipair)%jatom)
         use_transpose = lri%lattice%target_pairs(ipair)%iatom > lri%lattice%target_pairs(ipair)%jatom
         NULLIFY (ks_block)
         CALL dbcsr_get_block_p(matrix=ks_matrix(1, img)%matrix, row=irow, col=icol, block=ks_block, found=found)
         IF (.NOT. found) CPABORT("LRI-HFXk target block is absent on its DBCSR owner")
         IF (use_transpose) THEN
            IF (ANY(SHAPE(ks_block) /= [SIZE(received(ipair)%values, 2), SIZE(received(ipair)%values, 1)])) THEN
               CPABORT("Inconsistent transposed LRI-HFXk KS block dimensions")
            END IF
            ks_block = ks_block + scale*TRANSPOSE(received(ipair)%values)
         ELSE
            CPASSERT(ALL(SHAPE(ks_block) == SHAPE(received(ipair)%values)))
            ks_block = ks_block + scale*received(ipair)%values
         END IF
      END DO
      !$OMP END PARALLEL DO
      CALL timestop(timing_handle)
   END SUBROUTINE hfx_lri_write_rank_local_k_blocks

END MODULE hfx_lri_kp
