!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Tensor contractions used by the pair-local RI HFX backend.
! **************************************************************************************************
MODULE hfx_lri_contract
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hfx_lri_contract_exchange, hfx_lri_density_cell

CONTAINS

! **************************************************************************************************
!> \brief Returns the density-cell displacement for an exchange quartet.
!> \param cell_output R_IK
!> \param cell_left R_IJ
!> \param cell_right R_KL
!> \return R_JL = R_IK + R_KL - R_IJ
! **************************************************************************************************
   PURE FUNCTION hfx_lri_density_cell(cell_output, cell_left, cell_right) RESULT(cell_density)

      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell_output, cell_left, cell_right
      INTEGER, DIMENSION(3)                              :: cell_density

      cell_density(:) = cell_output(:) + cell_right(:) - cell_left(:)

   END FUNCTION hfx_lri_density_cell

! **************************************************************************************************
!> \brief Contracts C_IJ, V, C_KL and D_JL into K_IK.
!> \param c_left C_IJ(mu,sigma,X)
!> \param c_right C_KL(nu,lambda,Y)
!> \param v_aux V(X,Y)
!> \param density D_JL(sigma,lambda)
!> \param k_block K_IK(mu,nu)
! **************************************************************************************************
   SUBROUTINE hfx_lri_contract_exchange(c_left, c_right, v_aux, density, k_block)

      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: c_left, c_right
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: v_aux, density
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: k_block

      INTEGER                                            :: iaux, jaux
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: c_right_v, density_right

      IF (SIZE(c_left, 2) /= SIZE(density, 1) .OR. &
          SIZE(c_right, 2) /= SIZE(density, 2)) THEN
         CPABORT("Inconsistent AO dimensions in pair-local RI exchange contraction")
      END IF
      IF (SIZE(c_left, 3) /= SIZE(v_aux, 1) .OR. &
          SIZE(c_right, 3) /= SIZE(v_aux, 2)) THEN
         CPABORT("Inconsistent auxiliary dimensions in pair-local RI exchange contraction")
      END IF

      ! K(mu,nu) = sum[sigma,lambda,X,Y]
      !              C_IJ(mu,sigma,X) V_XY C_KL(nu,lambda,Y) D_JL(sigma,lambda).
      ! Combining V with the right C block first leaves one small GEMM per left auxiliary.
      ALLOCATE (c_right_v(SIZE(c_right, 1), SIZE(c_right, 2)), &
                density_right(SIZE(c_left, 2), SIZE(c_right, 1)), &
                k_block(SIZE(c_left, 1), SIZE(c_right, 1)), SOURCE=0.0_dp)
      DO iaux = 1, SIZE(c_left, 3)
         c_right_v(:, :) = 0.0_dp
         DO jaux = 1, SIZE(c_right, 3)
            c_right_v(:, :) = c_right_v(:, :) + v_aux(iaux, jaux)*c_right(:, :, jaux)
         END DO
         density_right(:, :) = MATMUL(density, TRANSPOSE(c_right_v))
         k_block(:, :) = k_block(:, :) + MATMUL(c_left(:, :, iaux), density_right)
      END DO

      DEALLOCATE (density_right, c_right_v)

   END SUBROUTINE hfx_lri_contract_exchange

END MODULE hfx_lri_contract
