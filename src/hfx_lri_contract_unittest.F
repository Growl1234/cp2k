!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM hfx_lri_contract_unittest
   USE hfx_lri_contract,                ONLY: hfx_lri_contract_workspace_type,&
                                              hfx_lri_loop3_aux_transform,&
                                              hfx_lri_loop3_central_contract,&
                                              hfx_lri_loop3_density_c_add,&
                                              hfx_lri_loop3_pack_aux,&
                                              hfx_lri_loop3_remote_contract,&
                                              hfx_lri_product_index
   USE kinds,                           ONLY: dp

   IMPLICIT NONE

   INTEGER, PARAMETER                                  :: loop3_center = 2, loop3_remote = 3, &
                                                          loop3_first = 4, loop3_output = 4, &
                                                          loop3_p = 2, loop3_q = 3
   INTEGER :: iaux, icenter, ifirst, ioutput, iremote, ix, iy
   REAL(dp)                                            :: error_max, scale
   REAL(dp), DIMENSION(loop3_center, loop3_remote, loop3_p) :: loop3_c
   REAL(dp), DIMENSION(loop3_first, loop3_center)       :: loop3_d_center
   REAL(dp), DIMENSION(loop3_first, loop3_remote)       :: loop3_d_remote
   REAL(dp), DIMENSION(loop3_p, loop3_q)                :: loop3_metric
   REAL(dp), DIMENSION(loop3_center*loop3_p, loop3_remote) :: loop3_packed_c
   REAL(dp), DIMENSION(loop3_center, loop3_output, loop3_p) :: loop3_u
   REAL(dp), DIMENSION(loop3_remote, loop3_output, loop3_q) :: loop3_s
   REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)            :: loop3_result, loop3_reference
   REAL(dp), ALLOCATABLE, DIMENSION(:, :)              :: loop3_block, loop3_block_reference, &
                                                          loop3_aux_flat
   TYPE(hfx_lri_contract_workspace_type)               :: workspace

   IF (hfx_lri_product_index(2, 3, 4) /= 7) THEN
      ERROR STOP "pair-local RI AO-product indexing is incorrect"
   END IF

   ! Exercise the packed BLAS kernels used by the finite-k loop3 implementation.  These tests use
   ! non-square dimensions deliberately, so a hidden transpose or leading-dimension error cannot
   ! pass accidentally.
   DO iaux = 1, loop3_p
      DO iremote = 1, loop3_remote
         DO icenter = 1, loop3_center
            loop3_c(icenter, iremote, iaux) = &
               0.013_dp*REAL(2*icenter - 3*iremote + 5*iaux, KIND=dp)
         END DO
      END DO
   END DO
   DO icenter = 1, loop3_center
      DO ifirst = 1, loop3_first
         loop3_d_center(ifirst, icenter) = &
            0.017_dp*REAL(3*ifirst - 2*icenter, KIND=dp)
      END DO
   END DO
   DO iremote = 1, loop3_remote
      DO ifirst = 1, loop3_first
         loop3_d_remote(ifirst, iremote) = &
            0.019_dp*REAL(2*ifirst + iremote - 4, KIND=dp)
      END DO
   END DO
   DO iy = 1, loop3_q
      DO ix = 1, loop3_p
         loop3_metric(ix, iy) = 0.023_dp*REAL(ix - 2*iy + 1, KIND=dp)
      END DO
   END DO
   DO iaux = 1, loop3_p
      loop3_packed_c((iaux - 1)*loop3_center + 1:iaux*loop3_center, :) = &
         loop3_c(:, :, iaux)
      DO ioutput = 1, loop3_output
         DO icenter = 1, loop3_center
            loop3_u(icenter, ioutput, iaux) = &
               0.029_dp*REAL(icenter + 2*ioutput - iaux, KIND=dp)
         END DO
      END DO
   END DO
   DO iy = 1, loop3_q
      DO ioutput = 1, loop3_output
         DO iremote = 1, loop3_remote
            loop3_s(iremote, ioutput, iy) = &
               0.031_dp*REAL(iremote - ioutput + 3*iy, KIND=dp)
         END DO
      END DO
   END DO

   CALL hfx_lri_loop3_density_c_add( &
      loop3_c, loop3_d_center, loop3_result, 0.7_dp, .TRUE., .FALSE.)
   CALL hfx_lri_loop3_density_c_add( &
      loop3_c, loop3_d_center, loop3_result, -0.2_dp, .FALSE., .FALSE.)
   ALLOCATE (loop3_reference(loop3_first, loop3_remote, loop3_p), SOURCE=0.0_dp)
   DO iaux = 1, loop3_p
      DO iremote = 1, loop3_remote
         DO ifirst = 1, loop3_first
            DO icenter = 1, loop3_center
               loop3_reference(ifirst, iremote, iaux) = &
                  loop3_reference(ifirst, iremote, iaux) + &
                  0.5_dp*loop3_d_center(ifirst, icenter)*loop3_c(icenter, iremote, iaux)
            END DO
         END DO
      END DO
   END DO
   error_max = MAXVAL(ABS(loop3_result - loop3_reference))
   scale = MAX(1.0_dp, MAXVAL(ABS(loop3_reference)))
   IF (error_max > 1.0E-13_dp*scale) THEN
      ERROR STOP "loop3 central-index density-C contraction is incorrect"
   END IF
   DEALLOCATE (loop3_result, loop3_reference)

   CALL hfx_lri_loop3_density_c_add( &
      loop3_c, loop3_d_remote, loop3_result, 1.0_dp, .TRUE., .TRUE., loop3_packed_c)
   ALLOCATE (loop3_reference(loop3_first, loop3_center, loop3_p), SOURCE=0.0_dp)
   DO iaux = 1, loop3_p
      DO icenter = 1, loop3_center
         DO ifirst = 1, loop3_first
            DO iremote = 1, loop3_remote
               loop3_reference(ifirst, icenter, iaux) = &
                  loop3_reference(ifirst, icenter, iaux) + &
                  loop3_d_remote(ifirst, iremote)*loop3_c(icenter, iremote, iaux)
            END DO
         END DO
      END DO
   END DO
   error_max = MAXVAL(ABS(loop3_result - loop3_reference))
   scale = MAX(1.0_dp, MAXVAL(ABS(loop3_reference)))
   IF (error_max > 1.0E-13_dp*scale) THEN
      ERROR STOP "loop3 packed remote-index density-C contraction is incorrect"
   END IF
   DEALLOCATE (loop3_result, loop3_reference)

   CALL hfx_lri_loop3_aux_transform( &
      loop3_s, loop3_metric, loop3_result, 1.0_dp, .TRUE.)
   ALLOCATE (loop3_reference(loop3_remote, loop3_output, loop3_p), SOURCE=0.0_dp)
   DO ix = 1, loop3_p
      DO ioutput = 1, loop3_output
         DO iremote = 1, loop3_remote
            DO iy = 1, loop3_q
               loop3_reference(iremote, ioutput, ix) = &
                  loop3_reference(iremote, ioutput, ix) + &
                  loop3_s(iremote, ioutput, iy)*loop3_metric(ix, iy)
            END DO
         END DO
      END DO
   END DO
   error_max = MAXVAL(ABS(loop3_result - loop3_reference))
   scale = MAX(1.0_dp, MAXVAL(ABS(loop3_reference)))
   IF (error_max > 1.0E-13_dp*scale) THEN
      ERROR STOP "loop3 auxiliary metric contraction is incorrect"
   END IF
   DEALLOCATE (loop3_result, loop3_reference)

   CALL hfx_lri_loop3_pack_aux(loop3_u, loop3_aux_flat)
   CALL hfx_lri_loop3_remote_contract( &
      loop3_packed_c, loop3_aux_flat, loop3_block)
   ALLOCATE (loop3_block_reference(loop3_remote, loop3_output), SOURCE=0.0_dp)
   DO ioutput = 1, loop3_output
      DO iremote = 1, loop3_remote
         DO iaux = 1, loop3_p
            DO icenter = 1, loop3_center
               loop3_block_reference(iremote, ioutput) = &
                  loop3_block_reference(iremote, ioutput) + &
                  loop3_c(icenter, iremote, iaux)*loop3_u(icenter, ioutput, iaux)
            END DO
         END DO
      END DO
   END DO
   error_max = MAXVAL(ABS(loop3_block - loop3_block_reference))
   scale = MAX(1.0_dp, MAXVAL(ABS(loop3_block_reference)))
   IF (error_max > 1.0E-13_dp*scale) THEN
      ERROR STOP "loop3 packed remote contraction is incorrect"
   END IF
   DEALLOCATE (loop3_block, loop3_block_reference, loop3_aux_flat)

   CALL hfx_lri_loop3_pack_aux(loop3_s, loop3_aux_flat)
   CALL hfx_lri_loop3_central_contract( &
      loop3_c, loop3_metric, loop3_aux_flat, loop3_block, workspace)
   ALLOCATE (loop3_block_reference(loop3_center, loop3_output), SOURCE=0.0_dp)
   DO ioutput = 1, loop3_output
      DO icenter = 1, loop3_center
         DO iy = 1, loop3_q
            DO iaux = 1, loop3_p
               DO iremote = 1, loop3_remote
                  loop3_block_reference(icenter, ioutput) = &
                     loop3_block_reference(icenter, ioutput) + &
                     loop3_c(icenter, iremote, iaux)*loop3_metric(iaux, iy)* &
                     loop3_s(iremote, ioutput, iy)
               END DO
            END DO
         END DO
      END DO
   END DO
   error_max = MAXVAL(ABS(loop3_block - loop3_block_reference))
   scale = MAX(1.0_dp, MAXVAL(ABS(loop3_block_reference)))
   IF (error_max > 1.0E-13_dp*scale) THEN
      ERROR STOP "loop3 central contraction is incorrect"
   END IF
   DEALLOCATE (loop3_block, loop3_block_reference, loop3_aux_flat)

END PROGRAM hfx_lri_contract_unittest
