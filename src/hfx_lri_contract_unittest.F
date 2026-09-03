!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM hfx_lri_contract_unittest
   USE hfx_lri_contract,                ONLY: hfx_lri_contract_exchange,&
                                              hfx_lri_contract_outer_accumulate,&
                                              hfx_lri_contract_outer_finalize,&
                                              hfx_lri_contract_outer_prepare,&
                                              hfx_lri_contract_workspace_type,&
                                              hfx_lri_density_cell,&
                                              hfx_lri_product_index,&
                                              hfx_lri_right_cell
   USE kinds,                           ONLY: dp

   IMPLICIT NONE

   INTEGER, PARAMETER                                  :: nlambda = 2, nmu = 2, nnu = 4, &
                                                          nsigma = 3, nx = 2, ny = 3
   INTEGER :: iflag_left, iflag_right, ilambda, imu, inu, isigma, ix, iy, &
              nlambda_case, nmu_case, nnu_case, nsigma_case
   INTEGER, DIMENSION(3)                               :: cell_density, cell_right
   LOGICAL                                             :: transpose_left, transpose_right
   REAL(dp)                                            :: error_max, scale
   REAL(dp), DIMENSION(nmu, nsigma, nx)                :: c_left
   REAL(dp), DIMENSION(nnu, nlambda, ny)               :: c_right
   REAL(dp), DIMENSION(nsigma, nlambda)                :: density
   REAL(dp), DIMENSION(nx, ny)                         :: v_aux
   REAL(dp), DIMENSION(nmu, nnu)                       :: k_reference
   REAL(dp), ALLOCATABLE, DIMENSION(:, :)              :: density_case, k_block, k_case, &
                                                          k_outer, k_reference_case
   TYPE(hfx_lri_contract_workspace_type)               :: workspace

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

   IF (hfx_lri_product_index(2, 3, 4) /= 7) THEN
      ERROR STOP "pair-local RI AO-product indexing is incorrect"
   END IF

   cell_density(:) = hfx_lri_density_cell([2, -1, 3], [-2, 4, 1], [1, -3, 5])
   IF (ANY(cell_density /= [5, -8, 7])) THEN
      ERROR STOP "pair-local RI exchange density-cell mapping is incorrect"
   END IF
   cell_right(:) = hfx_lri_right_cell(cell_density, [2, -1, 3], [-2, 4, 1])
   IF (ANY(cell_right /= [1, -3, 5])) THEN
      ERROR STOP "pair-local RI exchange inverse cell mapping is incorrect"
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

   ! A symmetric fitting-pair list stores C_ij only once. Exercise all zero-copy oriented views
   ! used for C_ji by the production contraction.
   DO iflag_left = 0, 1
      transpose_left = iflag_left == 1
      nmu_case = MERGE(nsigma, nmu, transpose_left)
      nsigma_case = MERGE(nmu, nsigma, transpose_left)
      DO iflag_right = 0, 1
         transpose_right = iflag_right == 1
         nnu_case = MERGE(nlambda, nnu, transpose_right)
         nlambda_case = MERGE(nnu, nlambda, transpose_right)
         ALLOCATE (density_case(nsigma_case, nlambda_case))
         DO ilambda = 1, nlambda_case
            DO isigma = 1, nsigma_case
               density_case(isigma, ilambda) = &
                  0.03_dp*REAL(5*isigma - 3*ilambda + iflag_left + 2*iflag_right, KIND=dp)
            END DO
         END DO
         CALL hfx_lri_contract_exchange(c_left, c_right, v_aux, density_case, k_case, &
                                        transpose_left, transpose_right)

         ALLOCATE (k_reference_case(nmu_case, nnu_case), SOURCE=0.0_dp)
         DO inu = 1, nnu_case
            DO imu = 1, nmu_case
               DO iy = 1, ny
                  DO ix = 1, nx
                     DO ilambda = 1, nlambda_case
                        DO isigma = 1, nsigma_case
                           IF (transpose_left) THEN
                              IF (transpose_right) THEN
                                 k_reference_case(imu, inu) = k_reference_case(imu, inu) + &
                                                              c_left(isigma, imu, ix)*v_aux(ix, iy)* &
                                                              c_right(ilambda, inu, iy)*density_case(isigma, ilambda)
                              ELSE
                                 k_reference_case(imu, inu) = k_reference_case(imu, inu) + &
                                                              c_left(isigma, imu, ix)*v_aux(ix, iy)* &
                                                              c_right(inu, ilambda, iy)*density_case(isigma, ilambda)
                              END IF
                           ELSE
                              IF (transpose_right) THEN
                                 k_reference_case(imu, inu) = k_reference_case(imu, inu) + &
                                                              c_left(imu, isigma, ix)*v_aux(ix, iy)* &
                                                              c_right(ilambda, inu, iy)*density_case(isigma, ilambda)
                              ELSE
                                 k_reference_case(imu, inu) = k_reference_case(imu, inu) + &
                                                              c_left(imu, isigma, ix)*v_aux(ix, iy)* &
                                                              c_right(inu, ilambda, iy)*density_case(isigma, ilambda)
                              END IF
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
         error_max = MAXVAL(ABS(k_case - k_reference_case))
         scale = MAX(1.0_dp, MAXVAL(ABS(k_reference_case)))
         IF (error_max > 1.0E-13_dp*scale) THEN
            ERROR STOP "oriented pair-local RI contraction disagrees with the explicit sum"
         END IF

         ! Treat every right auxiliary as belonging to the reusable outer-center slice. This must
         ! reproduce the ordinary contraction, including accumulation of multiple density blocks.
         ALLOCATE (k_outer(nmu_case, nnu_case), SOURCE=0.0_dp)
         CALL hfx_lri_contract_outer_prepare(c_left, v_aux, workspace, transpose_left)
         CALL hfx_lri_contract_outer_accumulate(c_right, density_case, workspace, &
                                                first=.TRUE., transpose_right=transpose_right)
         CALL hfx_lri_contract_outer_finalize(k_outer, workspace)
         error_max = MAXVAL(ABS(k_outer - k_reference_case))
         IF (error_max > 1.0E-13_dp*scale) THEN
            ERROR STOP "outer-center pair-local RI contraction disagrees with the explicit sum"
         END IF
         k_outer(:, :) = 0.0_dp
         CALL hfx_lri_contract_outer_accumulate(c_right, density_case, workspace, &
                                                first=.TRUE., transpose_right=transpose_right)
         CALL hfx_lri_contract_outer_accumulate(c_right, 0.25_dp*density_case, workspace, &
                                                first=.FALSE., transpose_right=transpose_right)
         CALL hfx_lri_contract_outer_finalize(k_outer, workspace)
         error_max = MAXVAL(ABS(k_outer - 1.25_dp*k_reference_case))
         IF (error_max > 1.0E-13_dp*scale) THEN
            ERROR STOP "outer-center pair-local RI accumulation is incorrect"
         END IF
         DEALLOCATE (density_case, k_case, k_outer, k_reference_case)
      END DO
   END DO

END PROGRAM hfx_lri_contract_unittest
