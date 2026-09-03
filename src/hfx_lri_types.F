!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE hfx_lri_types
   USE hfx_lri_contract,                ONLY: hfx_lri_contract_workspace_type
   USE kinds,                           ONLY: dp

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hfx_lri_centered_block_type, hfx_lri_k_block_type, hfx_lri_kp_type, &
             hfx_lri_metric_lattice_block_type, hfx_lri_pair_block_type, hfx_lri_kp_release

   TYPE hfx_lri_pair_block_type
      LOGICAL                                            :: onsite_pair = .FALSE.
      INTEGER                                            :: iatom = 0, ikind = 0, jatom = 0, &
                                                            jkind = 0, nao_i = 0, nao_j = 0, &
                                                            naux_i = 0, naux_j = 0, naux_pair = 0
      INTEGER                                            :: naux_centers = 0
      INTEGER, DIMENSION(3)                              :: cell = 0
      INTEGER, DIMENSION(2)                              :: aux_atom = 0, aux_kind = 0, &
                                                            aux_offset = 0, aux_size = 0
      INTEGER, DIMENSION(3, 2)                           :: aux_cell = 0
      REAL(dp), DIMENSION(3)                             :: rab = 0.0_dp
      REAL(dp), DIMENSION(3, 2)                          :: aux_r = 0.0_dp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: c_fit
   END TYPE hfx_lri_pair_block_type

   TYPE hfx_lri_k_block_type
      INTEGER                                            :: iatom = 0, ikind = 0, jatom = 0, &
                                                            jkind = 0, nao_i = 0, nao_j = 0
      INTEGER, DIMENSION(3)                              :: cell = 0
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: k_block
   END TYPE hfx_lri_k_block_type

   ! One coefficient slice with its auxiliary functions fixed on the first atom.  Lattice cells
   ! are reduced to the finite k-mesh quotient before these blocks are stored.
   TYPE hfx_lri_centered_block_type
      LOGICAL                                            :: inner_active = .TRUE.
      REAL(dp)                                           :: norm2 = 0.0_dp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: values
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: remote_flat
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: inner_values
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: inner_remote_flat
   END TYPE hfx_lri_centered_block_type

   ! Periodized two-center metric block on the finite k-mesh lattice quotient.
   TYPE hfx_lri_metric_lattice_block_type
      REAL(dp)                                           :: norm2 = 0.0_dp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: values
   END TYPE hfx_lri_metric_lattice_block_type

   TYPE hfx_lri_kp_type
      LOGICAL :: initialized = .FALSE.
      LOGICAL :: loop3_ready = .FALSE.

      INTEGER :: npair_blocks = 0
      REAL(dp) :: exchange_trace_raw = 0.0_dp

      INTEGER                                            :: lattice_natom = 0, lattice_ncell = 0, &
                                                            lattice_total_nao = 0
      INTEGER, DIMENSION(3)                              :: lattice_grid = 1, lattice_twist = 1
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_ao_offset, atom_kind, atom_nao, &
                                                            atom_naux, centered_atom_count, &
                                                            centered_atom_offset, centered_atom_refs, &
                                                            centered_remote_count, &
                                                            centered_remote_offset, &
                                                            centered_remote_refs
      TYPE(hfx_lri_centered_block_type), ALLOCATABLE, DIMENSION(:) :: centered_blocks
      TYPE(hfx_lri_metric_lattice_block_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: lattice_metric

      TYPE(hfx_lri_pair_block_type), ALLOCATABLE, DIMENSION(:) :: pair_blocks
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

      IF (ASSOCIATED(lri_kp)) DEALLOCATE (lri_kp)

   END SUBROUTINE hfx_lri_kp_release

END MODULE hfx_lri_types
