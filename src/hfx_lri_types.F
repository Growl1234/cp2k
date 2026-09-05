!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE hfx_lri_types
   USE hfx_lri_contract,                ONLY: hfx_lri_contract_workspace_type
   USE hfx_lri_parallel,                ONLY: hfx_lri_metric_lattice_block_type => hfx_lri_block_type,&
                                              hfx_lri_parallel_type
   USE kinds,                           ONLY: dp

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hfx_lri_centered_block_type, hfx_lri_k_block_type, hfx_lri_kp_type, &
             hfx_lri_metric_lattice_block_type, hfx_lri_pair_block_type, hfx_lri_kp_release

   TYPE hfx_lri_pair_block_type
      LOGICAL                                            :: onsite_pair = .FALSE., fit_active = .FALSE.
      INTEGER                                            :: iatom = 0, ikind = 0, jatom = 0, &
                                                            jkind = 0, nao_i = 0, nao_j = 0, &
                                                            naux_i = 0, naux_j = 0, naux_pair = 0
      INTEGER                                            :: naux_centers = 0
      INTEGER, DIMENSION(3)                              :: cell = 0
      INTEGER, DIMENSION(2)                              :: aux_kind = 0
      REAL(dp), DIMENSION(3)                             :: rab = 0.0_dp
      REAL(dp), DIMENSION(3, 2)                          :: aux_r = 0.0_dp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: c_fit
   END TYPE hfx_lri_pair_block_type

   TYPE hfx_lri_k_block_type
      INTEGER                                            :: iatom = 0, jatom = 0
      INTEGER, DIMENSION(3)                              :: cell = 0
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: k_block
   END TYPE hfx_lri_k_block_type

   ! One coefficient slice with its auxiliary functions fixed on the first atom.  Lattice cells
   ! are reduced to the finite k-mesh quotient before these blocks are stored.
   TYPE hfx_lri_centered_block_type
      LOGICAL                                            :: inner_active = .TRUE.
      LOGICAL                                            :: present = .FALSE., inner_override = .FALSE.
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: values
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: remote_flat
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: inner_values
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: inner_remote_flat
   END TYPE hfx_lri_centered_block_type

   TYPE hfx_lri_kp_type
      TYPE(hfx_lri_parallel_type) :: parallel = hfx_lri_parallel_type()
      INTEGER                                            :: lattice_natom = 0, lattice_ncell = 0, &
                                                            lattice_total_nao = 0
      INTEGER, DIMENSION(3)                              :: lattice_grid = 1, lattice_twist = 1
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_ao_offset, atom_kind, atom_naux, &
                                                            centered_atom_offset, centered_atom_refs, &
                                                            centered_remote_offset, &
                                                            centered_remote_refs
      TYPE(hfx_lri_centered_block_type), ALLOCATABLE, DIMENSION(:) :: centered_blocks
      TYPE(hfx_lri_metric_lattice_block_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: lattice_metric

      TYPE(hfx_lri_pair_block_type), ALLOCATABLE, DIMENSION(:) :: pair_blocks
      TYPE(hfx_lri_pair_block_type), ALLOCATABLE, DIMENSION(:) :: target_pairs
      INTEGER, ALLOCATABLE, DIMENSION(:) :: k_target
      TYPE(hfx_lri_contract_workspace_type), ALLOCATABLE, DIMENSION(:) :: thread_workspaces

      INTEGER :: nk_blocks = 0
      TYPE(hfx_lri_k_block_type), ALLOCATABLE, DIMENSION(:) :: k_blocks
   END TYPE hfx_lri_kp_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param lri_kp ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_kp_release(lri_kp)
      TYPE(hfx_lri_kp_type), POINTER                     :: lri_kp

      IF (ASSOCIATED(lri_kp)) THEN
         IF (ASSOCIATED(lri_kp%parallel%subgroup)) THEN
            CALL lri_kp%parallel%subgroup%free()
            DEALLOCATE (lri_kp%parallel%subgroup)
         END IF
         DEALLOCATE (lri_kp)
      END IF

   END SUBROUTINE hfx_lri_kp_release

END MODULE hfx_lri_types
