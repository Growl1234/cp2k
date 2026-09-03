!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM hfx_lri_contract_unittest
   USE hfx_lri_contract,                ONLY: hfx_lri_contract_exchange,&
                                              hfx_lri_density_cell
   USE kinds,                           ONLY: dp

   IMPLICIT NONE

   INTEGER, PARAMETER                                  :: nlambda = 2, nmu = 2, nnu = 4, &
                                                          nsigma = 3, nx = 2, ny = 3
   INTEGER                                             :: ilambda, imu, inu, isigma, ix, iy
   INTEGER, DIMENSION(3)                               :: cell_density
   REAL(dp)                                            :: error_max, scale
   REAL(dp), DIMENSION(nmu, nsigma, nx)                :: c_left
   REAL(dp), DIMENSION(nnu, nlambda, ny)               :: c_right
   REAL(dp), DIMENSION(nsigma, nlambda)                :: density
   REAL(dp), DIMENSION(nx, ny)                         :: v_aux
   REAL(dp), DIMENSION(nmu, nnu)                       :: k_reference
   REAL(dp), ALLOCATABLE, DIMENSION(:, :)              :: k_block

   DO ix = 1, nx
      DO isigma = 1, nsigma
         DO imu = 1, nmu
            c_left(imu, isigma, ix) = &
               0.07_dp*REAL(2*imu - 3*isigma + 5*ix, KIND=dp)
         END DO
      END DO
   END DO
   DO iy = 1, ny
      DO ilambda = 1, nlambda
         DO inu = 1, nnu
            c_right(inu, ilambda, iy) = &
               0.05_dp*REAL(-inu + 4*ilambda + 2*iy, KIND=dp)
         END DO
      END DO
   END DO
   DO ilambda = 1, nlambda
      DO isigma = 1, nsigma
         density(isigma, ilambda) = &
            0.11_dp*REAL(3*isigma - 2*ilambda, KIND=dp)
      END DO
   END DO
   DO iy = 1, ny
      DO ix = 1, nx
         v_aux(ix, iy) = 0.09_dp*REAL(2*ix + iy - 2, KIND=dp)
      END DO
   END DO

   CALL hfx_lri_contract_exchange(c_left, c_right, v_aux, density, k_block)

   cell_density(:) = hfx_lri_density_cell([2, -1, 3], [-2, 4, 1], [1, -3, 5])
   IF (ANY(cell_density /= [5, -8, 7])) THEN
      ERROR STOP "pair-local RI exchange density-cell mapping is incorrect"
   END IF

   k_reference(:, :) = 0.0_dp
   DO inu = 1, nnu
      DO imu = 1, nmu
         DO iy = 1, ny
            DO ix = 1, nx
               DO ilambda = 1, nlambda
                  DO isigma = 1, nsigma
                     k_reference(imu, inu) = k_reference(imu, inu) + &
                                             c_left(imu, isigma, ix)*v_aux(ix, iy)* &
                                             c_right(inu, ilambda, iy)*density(isigma, ilambda)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO

   IF (ANY(SHAPE(k_block) /= [nmu, nnu])) THEN
      ERROR STOP "pair-local RI exchange contraction returned the wrong shape"
   END IF
   error_max = MAXVAL(ABS(k_block - k_reference))
   scale = MAX(1.0_dp, MAXVAL(ABS(k_reference)))
   IF (error_max > 1.0E-13_dp*scale) THEN
      ERROR STOP "pair-local RI exchange contraction disagrees with the explicit four-index sum"
   END IF

END PROGRAM hfx_lri_contract_unittest
