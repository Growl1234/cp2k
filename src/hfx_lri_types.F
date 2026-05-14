!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE hfx_lri_types
   USE kinds,                           ONLY: dp,&
                                              int_8

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hfx_lri_k_block_type, hfx_lri_kp_type, hfx_lri_pair_block_type, hfx_lri_kp_release

   TYPE hfx_lri_pair_block_type
      INTEGER                                            :: iatom = 0, ikind = 0, jatom = 0, &
                                                            jkind = 0, nao_i = 0, nao_j = 0, &
                                                            naux_pair = 0
      INTEGER, DIMENSION(3)                              :: cell = 0
      REAL(dp), DIMENSION(3)                             :: rab = 0.0_dp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: c_fit
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: v_fit
   END TYPE hfx_lri_pair_block_type

   TYPE hfx_lri_k_block_type
      INTEGER                                            :: pair_index = 0
      INTEGER                                            :: iatom = 0, ikind = 0, jatom = 0, &
                                                            jkind = 0, nao_i = 0, nao_j = 0
      INTEGER, DIMENSION(3)                              :: cell = 0
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: k_block
   END TYPE hfx_lri_k_block_type

   TYPE hfx_lri_kp_type
      LOGICAL :: initialized = .FALSE.

      INTEGER :: npairs = 0
      INTEGER :: npair_blocks = 0

      INTEGER(int_8) :: sum_nao_pair = 0_int_8
      INTEGER(int_8) :: sum_naux_pair = 0_int_8
      INTEGER(int_8) :: c_nreal_est = 0_int_8
      INTEGER(int_8) :: c_nreal_stored = 0_int_8

      REAL(dp) :: c_mem_mib_est = 0.0_dp
      REAL(dp) :: c_mem_mib_stored = 0.0_dp

      TYPE(hfx_lri_pair_block_type), ALLOCATABLE, DIMENSION(:) :: pair_blocks

      TYPE(hfx_lri_k_block_type) :: debug_k_block
   END TYPE hfx_lri_kp_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param lri_kp ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_kp_release(lri_kp)
      TYPE(hfx_lri_kp_type), POINTER                     :: lri_kp

      INTEGER                                            :: ipair

      IF (ASSOCIATED(lri_kp)) THEN
         IF (ALLOCATED(lri_kp%debug_k_block%k_block)) THEN
            DEALLOCATE (lri_kp%debug_k_block%k_block)
         END IF
         IF (ALLOCATED(lri_kp%pair_blocks)) THEN
            DO ipair = 1, SIZE(lri_kp%pair_blocks)
               IF (ALLOCATED(lri_kp%pair_blocks(ipair)%c_fit)) THEN
                  DEALLOCATE (lri_kp%pair_blocks(ipair)%c_fit)
               END IF
               IF (ALLOCATED(lri_kp%pair_blocks(ipair)%v_fit)) THEN
                  DEALLOCATE (lri_kp%pair_blocks(ipair)%v_fit)
               END IF
            END DO
            DEALLOCATE (lri_kp%pair_blocks)
         END IF
         DEALLOCATE (lri_kp)
      END IF
   END SUBROUTINE hfx_lri_kp_release

END MODULE hfx_lri_types
