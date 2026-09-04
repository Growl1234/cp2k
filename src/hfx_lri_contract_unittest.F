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
   USE hfx_lri_lattice,                 ONLY: &
        hfx_lri_lattice_block_from_index, hfx_lri_lattice_block_index, &
        hfx_lri_lattice_cell_from_index, hfx_lri_lattice_cell_index, &
        hfx_lri_lattice_density_tolerance, hfx_lri_lattice_displacement_cell, &
        hfx_lri_lattice_pair_representative, hfx_lri_lattice_project_pair, &
        hfx_lri_lattice_project_sample, hfx_lri_lattice_project_self, hfx_lri_lattice_reduce, &
        hfx_lri_reciprocal_diagnostics
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

   CALL test_four_topology_contraction()
   CALL test_displacement_cell()
   CALL test_pair_representatives()
   CALL test_lattice_quotient()
   CALL test_density_projection()
   CALL test_storage_trace_alias()
   CALL test_reciprocal_diagnostics()

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

CONTAINS

! **************************************************************************************************
!> \brief Checks displacement-to-cell mapping across a wrapped atomic coordinate.
! **************************************************************************************************
   SUBROUTINE test_displacement_cell()

      INTEGER, DIMENSION(3)                              :: cell
      REAL(dp)                                           :: residual
      REAL(dp), DIMENSION(3, 3)                          :: h_inv

      h_inv(:, :) = 0.0_dp
      h_inv(1, 1) = 1.0_dp
      h_inv(2, 2) = 1.0_dp
      h_inv(3, 3) = 1.0_dp
      CALL hfx_lri_lattice_displacement_cell( &
         [0.75_dp, -0.25_dp, -0.25_dp], &
         [0.25_dp, 0.25_dp, 0.25_dp], [0.0_dp, 0.0_dp, 0.0_dp], &
         h_inv, cell, residual)
      IF (ANY(cell /= [1, 0, 0]) .OR. residual > 1.0E-14_dp) THEN
         ERROR STOP "wrapped LRI metric displacement maps to the wrong lattice cell"
      END IF

   END SUBROUTINE test_displacement_cell

! **************************************************************************************************
!> \brief Checks that diagonal R/-R pairs are fitted exactly once while other pairs are retained.
! **************************************************************************************************
   SUBROUTINE test_pair_representatives()

      IF (.NOT. hfx_lri_lattice_pair_representative(2, 2, [0, 0, 0])) THEN
         ERROR STOP "on-site LRI pair was rejected"
      END IF
      IF (.NOT. hfx_lri_lattice_pair_representative(2, 2, [1, -3, 0]) .OR. &
          hfx_lri_lattice_pair_representative(2, 2, [-1, 3, 0])) THEN
         ERROR STOP "first diagonal LRI translation orbit has inconsistent representatives"
      END IF
      IF (.NOT. hfx_lri_lattice_pair_representative(2, 2, [0, 2, -4]) .OR. &
          hfx_lri_lattice_pair_representative(2, 2, [0, -2, 4])) THEN
         ERROR STOP "second diagonal LRI translation orbit has inconsistent representatives"
      END IF
      IF (.NOT. hfx_lri_lattice_pair_representative(1, 2, [-4, 0, 1]) .OR. &
          .NOT. hfx_lri_lattice_pair_representative(2, 1, [4, 0, -1])) THEN
         ERROR STOP "off-diagonal LRI pair was unexpectedly rejected"
      END IF

   END SUBROUTINE test_pair_representatives

! **************************************************************************************************
!> \brief Checks quotient indexing and real Monkhorst-Pack boundary phases.
! **************************************************************************************************
   SUBROUTINE test_lattice_quotient()

      INTEGER, DIMENSION(3), PARAMETER                   :: grid = [2, 3, 4], twist = [-1, 1, -1]
      INTEGER, PARAMETER                                 :: natom = 3

      INTEGER                                            :: block_index, cell_index, &
                                                            decoded_cell_index, iatom, &
                                                            iatom_decoded, ix, iy, iz, jatom, &
                                                            jatom_decoded
      INTEGER, DIMENSION(3)                              :: cell, cell_a, cell_b, combined, reduced, &
                                                            reduced_a, reduced_b, &
                                                            reduced_combined, reduced_sum
      REAL(dp)                                           :: factor, factor_a, factor_b, &
                                                            factor_combined, factor_sum

      DO cell_index = 1, PRODUCT(grid)
         CALL hfx_lri_lattice_cell_from_index(cell_index, grid, cell)
         IF (hfx_lri_lattice_cell_index(cell, grid) /= cell_index) THEN
            ERROR STOP "pair-local RI quotient cell indexing is not invertible"
         END IF
      END DO

      DO iatom = 1, natom
         DO jatom = 1, natom
            DO cell_index = 1, PRODUCT(grid)
               block_index = hfx_lri_lattice_block_index( &
                             iatom, jatom, cell_index, natom, PRODUCT(grid))
               CALL hfx_lri_lattice_block_from_index( &
                  block_index, natom, PRODUCT(grid), iatom_decoded, jatom_decoded, &
                  decoded_cell_index)
               IF (iatom_decoded /= iatom .OR. jatom_decoded /= jatom .OR. &
                   decoded_cell_index /= cell_index) THEN
                  ERROR STOP "pair-local RI quotient block indexing is not invertible"
               END IF
            END DO
         END DO
      END DO

      CALL hfx_lri_lattice_reduce([2, 0, 0], grid, twist, reduced, factor)
      IF (ANY(reduced /= [0, 0, 0]) .OR. factor /= -1.0_dp) THEN
         ERROR STOP "pair-local RI positive boundary phase is incorrect"
      END IF
      CALL hfx_lri_lattice_reduce([-1, 0, 0], grid, twist, reduced, factor)
      IF (ANY(reduced /= [1, 0, 0]) .OR. factor /= -1.0_dp) THEN
         ERROR STOP "pair-local RI negative boundary phase is incorrect"
      END IF
      CALL hfx_lri_lattice_reduce([-1, 0, -4], grid, twist, reduced, factor)
      IF (ANY(reduced /= [1, 0, 0]) .OR. factor /= 1.0_dp) THEN
         ERROR STOP "pair-local RI multi-direction boundary phase is incorrect"
      END IF

      DO ix = -5, 5
         DO iy = -4, 4
            DO iz = -3, 3
               cell_a = [ix, iy, iz]
               cell_b = [2 - cell_a(2), cell_a(3) - 1, -1 - cell_a(1)]
               combined = cell_a + cell_b
               CALL hfx_lri_lattice_reduce(cell_a, grid, twist, reduced_a, factor_a)
               CALL hfx_lri_lattice_reduce(cell_b, grid, twist, reduced_b, factor_b)
               CALL hfx_lri_lattice_reduce( &
                  reduced_a + reduced_b, grid, twist, reduced_sum, factor_sum)
               CALL hfx_lri_lattice_reduce( &
                  combined, grid, twist, reduced_combined, factor_combined)
               IF (ANY(reduced_sum /= reduced_combined) .OR. &
                   factor_a*factor_b*factor_sum /= factor_combined) THEN
                  ERROR STOP "pair-local RI quotient reduction violates phase composition"
               END IF
            END DO
         END DO
      END DO

   END SUBROUTINE test_lattice_quotient

! **************************************************************************************************
!> \brief Checks cross-grid quotient averaging, Hermitian projection, and adaptive tolerances.
! **************************************************************************************************
   SUBROUTINE test_density_projection()

      INTEGER                                            :: sample_count
      REAL(dp)                                           :: factor, tolerance_filtered, &
                                                            tolerance_scaled
      REAL(dp), DIMENSION(2, 2)                          :: average, expected, left, right, &
                                                            sample_a, sample_b, self_block

      sample_a(:, :) = RESHAPE([1.0_dp, 0.2_dp, -0.3_dp, 0.7_dp], [2, 2])
      sample_b(:, :) = RESHAPE([-0.4_dp, 0.5_dp, 0.1_dp, -0.8_dp], [2, 2])
      average(:, :) = 0.0_dp
      sample_count = 0
      CALL hfx_lri_lattice_project_sample(average, sample_count, sample_a, 1.0_dp)
      CALL hfx_lri_lattice_project_sample(average, sample_count, sample_b, -1.0_dp)
      expected(:, :) = 0.5_dp*(sample_a - sample_b)
      IF (sample_count /= 2 .OR. MAXVAL(ABS(average - expected)) > 1.0E-15_dp) THEN
         ERROR STOP "cross-grid density samples were not averaged in the twisted quotient"
      END IF

      left(:, :) = RESHAPE([0.8_dp, -0.1_dp, 0.4_dp, 1.2_dp], [2, 2])
      right(:, :) = RESHAPE([0.3_dp, 0.7_dp, -0.2_dp, 0.5_dp], [2, 2])
      factor = -1.0_dp
      expected(:, :) = (2.0_dp*left + factor*TRANSPOSE(right))/3.0_dp
      CALL hfx_lri_lattice_project_pair(left, 2, right, 1, factor)
      IF (MAXVAL(ABS(left - expected)) > 1.0E-15_dp .OR. &
          MAXVAL(ABS(left - factor*TRANSPOSE(right))) > 1.0E-15_dp) THEN
         ERROR STOP "weighted density quotient projection is not Hermitian"
      END IF
      expected(:, :) = left
      CALL hfx_lri_lattice_project_pair(left, 2, right, 1, factor)
      IF (MAXVAL(ABS(left - expected)) > 1.0E-15_dp) THEN
         ERROR STOP "density quotient projection is not idempotent"
      END IF

      self_block(:, :) = RESHAPE([0.6_dp, 0.2_dp, -0.4_dp, 0.9_dp], [2, 2])
      CALL hfx_lri_lattice_project_self(self_block, -1.0_dp)
      IF (MAXVAL(ABS(self_block + TRANSPOSE(self_block))) > 1.0E-15_dp) THEN
         ERROR STOP "self-reverse twisted density block was not projected correctly"
      END IF

      tolerance_scaled = hfx_lri_lattice_density_tolerance(10.0_dp, 0.0_dp)
      tolerance_filtered = hfx_lri_lattice_density_tolerance(0.0_dp, 1.0E-6_dp)
      IF (tolerance_scaled < 10.0_dp*SQRT(EPSILON(1.0_dp)) .OR. &
          tolerance_filtered < 8.0E-6_dp) THEN
         ERROR STOP "density quotient tolerance does not track scale and EPS_FILTER"
      END IF

   END SUBROUTINE test_density_projection

! **************************************************************************************************
!> \brief Checks trace reconstruction when the KS-storage image is a twisted quotient alias.
! **************************************************************************************************
   SUBROUTINE test_storage_trace_alias()

      INTEGER, DIMENSION(3), PARAMETER                   :: grid = [2, 1, 1], twist = [-1, 1, 1]

      INTEGER, DIMENSION(3)                              :: reduced
      REAL(dp)                                           :: factor, quotient_trace, storage_trace
      REAL(dp), DIMENSION(2, 2)                          :: density_quotient, kernel_quotient, &
                                                            kernel_stored

      density_quotient(:, :) = RESHAPE( &
                               [1.1_dp, -0.2_dp, 0.4_dp, 0.7_dp], [2, 2])
      kernel_quotient(:, :) = RESHAPE( &
                              [0.8_dp, 0.3_dp, -0.1_dp, 1.2_dp], [2, 2])
      CALL hfx_lri_lattice_reduce([2, 0, 0], grid, twist, reduced, factor)
      IF (ANY(reduced /= [0, 0, 0]) .OR. factor /= -1.0_dp) THEN
         ERROR STOP "unexpected storage alias in pair-local RI trace test"
      END IF

      ! Both the physical density and the stored physical K acquire the boundary factor.  The
      ! production trace starts from quotient D and stored K, and must therefore restore the
      ! density factor explicitly.  Looking for this exact physical density image would instead
      ! miss the contribution when CP2K selected another symmetric storage representative.
      kernel_stored(:, :) = factor*kernel_quotient
      quotient_trace = SUM(density_quotient*kernel_quotient)
      storage_trace = factor*SUM(density_quotient*kernel_stored)
      IF (ABS(storage_trace - quotient_trace) > &
          1.0E-14_dp*MAX(1.0_dp, ABS(quotient_trace))) THEN
         ERROR STOP "twisted KS-storage alias does not preserve the quotient trace"
      END IF

   END SUBROUTINE test_storage_trace_alias

! **************************************************************************************************
!> \brief Checks reciprocal reconstruction, Hermiticity, spectra, and Parseval trace consistency.
! **************************************************************************************************
   SUBROUTINE test_reciprocal_diagnostics()

      INTEGER, DIMENSION(3), PARAMETER                   :: grid = [2, 1, 1], twist = [-1, 1, 1]

      REAL(dp) :: density_eig_max, density_eig_min, density_hermitian, kernel_eig_max, &
         kernel_eig_min, kernel_hermitian, real_trace, reciprocal_trace
      REAL(dp), DIMENSION(2, 2, 2)                       :: density, kernel

      density(:, :, 1) = RESHAPE([2.0_dp, 0.3_dp, 0.3_dp, 1.5_dp], [2, 2])
      density(:, :, 2) = RESHAPE([0.0_dp, -0.2_dp, 0.2_dp, 0.0_dp], [2, 2])
      kernel(:, :, 1) = RESHAPE([1.2_dp, 0.1_dp, 0.1_dp, 0.9_dp], [2, 2])
      kernel(:, :, 2) = RESHAPE([0.0_dp, -0.15_dp, 0.15_dp, 0.0_dp], [2, 2])
      real_trace = SUM(density*kernel)

      CALL hfx_lri_reciprocal_diagnostics( &
         density, kernel, grid, twist, density_hermitian, kernel_hermitian, &
         density_eig_min, density_eig_max, kernel_eig_min, kernel_eig_max, &
         reciprocal_trace)

      IF (MAX(density_hermitian, kernel_hermitian) > 1.0E-13_dp) THEN
         ERROR STOP "twisted reciprocal matrices are not Hermitian"
      END IF
      IF (density_eig_min <= 0.0_dp .OR. density_eig_max <= density_eig_min) THEN
         ERROR STOP "unexpected density spectrum after reciprocal reconstruction"
      END IF
      IF (kernel_eig_min <= 0.0_dp .OR. kernel_eig_max <= kernel_eig_min) THEN
         ERROR STOP "unexpected kernel spectrum after reciprocal reconstruction"
      END IF
      IF (ABS(reciprocal_trace - real_trace) > &
          1.0E-13_dp*MAX(1.0_dp, ABS(real_trace))) THEN
         ERROR STOP "real- and reciprocal-space traces violate Parseval consistency"
      END IF

   END SUBROUTINE test_reciprocal_diagnostics

! **************************************************************************************************
!> \brief Compares all four pair-local auxiliary-center topologies against explicit sums.
! **************************************************************************************************
   SUBROUTINE test_four_topology_contraction()

      INTEGER, PARAMETER                                 :: nlambda = 2, nmu = 2, nnu = 4, npi = 3, &
                                                            npo = 2, nqi = 2, nqo = 4, nsigma = 3, &
                                                            ntopology = 4

      INTEGER                                            :: ilambda, imu, inu, ip, iq, isigma, &
                                                            itopology
      REAL(dp)                                           :: error_max, reference_scale
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: k_combined, k_delta, k_optimized, &
                                                            k_reference, k_zero
      REAL(dp), DIMENSION(nlambda, nnu, nqi)             :: c_right_inner
      REAL(dp), DIMENSION(nmu, nsigma, npo)              :: c_left_outer
      REAL(dp), DIMENSION(nnu, nlambda, nqo)             :: c_right_outer
      REAL(dp), DIMENSION(npi, nqi)                      :: metric_ii
      REAL(dp), DIMENSION(npi, nqo)                      :: metric_io
      REAL(dp), DIMENSION(npo, nqi)                      :: metric_oi
      REAL(dp), DIMENSION(npo, nqo)                      :: metric_oo
      REAL(dp), DIMENSION(nsigma, nlambda)               :: density, density_delta
      REAL(dp), DIMENSION(nsigma, nmu, npi)              :: c_left_inner

      DO ip = 1, npo
         DO isigma = 1, nsigma
            DO imu = 1, nmu
               c_left_outer(imu, isigma, ip) = &
                  0.013_dp*REAL(2*imu - 3*isigma + 5*ip, KIND=dp)
            END DO
         END DO
      END DO
      DO ip = 1, npi
         DO imu = 1, nmu
            DO isigma = 1, nsigma
               c_left_inner(isigma, imu, ip) = &
                  0.017_dp*REAL(-isigma + 4*imu + 3*ip, KIND=dp)
            END DO
         END DO
      END DO
      DO iq = 1, nqo
         DO ilambda = 1, nlambda
            DO inu = 1, nnu
               c_right_outer(inu, ilambda, iq) = &
                  0.019_dp*REAL(-2*inu + ilambda + 4*iq, KIND=dp)
            END DO
         END DO
      END DO
      DO iq = 1, nqi
         DO inu = 1, nnu
            DO ilambda = 1, nlambda
               c_right_inner(ilambda, inu, iq) = &
                  0.023_dp*REAL(3*ilambda - inu + 2*iq, KIND=dp)
            END DO
         END DO
      END DO
      DO ilambda = 1, nlambda
         DO isigma = 1, nsigma
            density(isigma, ilambda) = &
               0.029_dp*REAL(5*isigma - 2*ilambda - 3, KIND=dp)
            density_delta(isigma, ilambda) = &
               0.011_dp*REAL(isigma + 3*ilambda - 4, KIND=dp)
         END DO
      END DO
      DO iq = 1, nqo
         DO ip = 1, npo
            metric_oo(ip, iq) = 0.031_dp*REAL(2*ip - iq + 1, KIND=dp)
         END DO
         DO ip = 1, npi
            metric_io(ip, iq) = 0.037_dp*REAL(-ip + 3*iq - 2, KIND=dp)
         END DO
      END DO
      DO iq = 1, nqi
         DO ip = 1, npo
            metric_oi(ip, iq) = 0.041_dp*REAL(ip + 2*iq - 1, KIND=dp)
         END DO
         DO ip = 1, npi
            metric_ii(ip, iq) = 0.043_dp*REAL(2*ip + iq - 3, KIND=dp)
         END DO
      END DO

      CALL four_topology_reference( &
         c_left_outer, c_left_inner, c_right_outer, c_right_inner, &
         metric_oo, metric_oi, metric_io, metric_ii, density, k_reference)
      CALL four_topology_optimized( &
         c_left_outer, c_left_inner, c_right_outer, c_right_inner, &
         metric_oo, metric_oi, metric_io, metric_ii, density, k_optimized)
      DO itopology = 1, ntopology
         error_max = MAXVAL(ABS(k_optimized(:, :, itopology) - &
                                k_reference(:, :, itopology)))
         reference_scale = MAX(1.0_dp, MAXVAL(ABS(k_reference(:, :, itopology))))
         IF (error_max > 1.0E-13_dp*reference_scale) THEN
            ERROR STOP "pair-local RI four-topology contraction disagrees with the explicit sum"
         END IF
      END DO

      CALL four_topology_optimized( &
         c_left_outer, c_left_inner, c_right_outer, c_right_inner, &
         metric_oo, metric_oi, metric_io, metric_ii, density_delta, k_delta)
      CALL four_topology_optimized( &
         c_left_outer, c_left_inner, c_right_outer, c_right_inner, &
         metric_oo, metric_oi, metric_io, metric_ii, density + density_delta, k_combined)
      reference_scale = MAX(1.0_dp, MAXVAL(ABS(k_combined)))
      IF (MAXVAL(ABS(k_combined - k_optimized - k_delta)) > &
          1.0E-13_dp*reference_scale) THEN
         ERROR STOP "pair-local RI four-topology contraction is not linear in the density"
      END IF

      CALL four_topology_optimized( &
         c_left_outer, c_left_inner, c_right_outer, c_right_inner, &
         metric_oo, metric_oi, metric_io, metric_ii, 0.0_dp*density, k_zero)
      IF (MAXVAL(ABS(k_zero)) > 1.0E-15_dp) THEN
         ERROR STOP "pair-local RI four-topology contraction is nonzero for zero density"
      END IF

      DEALLOCATE (k_combined, k_delta, k_optimized, k_reference, k_zero)

   END SUBROUTINE test_four_topology_contraction

! **************************************************************************************************
!> \brief Transparent four-topology contraction used only as a correctness oracle.
!> \param c_left_outer C(mu,sigma,P) with P on the outer/target center
!> \param c_left_inner C(sigma,mu,P) with P on the inner/contracted center
!> \param c_right_outer C(nu,lambda,Q) with Q on the outer/output center
!> \param c_right_inner C(lambda,nu,Q) with Q on the inner/contracted center
!> \param metric_oo outer-outer auxiliary metric
!> \param metric_oi outer-inner auxiliary metric
!> \param metric_io inner-outer auxiliary metric
!> \param metric_ii inner-inner auxiliary metric
!> \param density D(sigma,lambda)
!> \param k_terms K(mu,nu,topology), ordered as OO, OI, IO, II
! **************************************************************************************************
   SUBROUTINE four_topology_reference( &
      c_left_outer, c_left_inner, c_right_outer, c_right_inner, &
      metric_oo, metric_oi, metric_io, metric_ii, density, k_terms)

      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: c_left_outer, c_left_inner, &
                                                            c_right_outer, c_right_inner
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: metric_oo, metric_oi, metric_io, &
                                                            metric_ii, density
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: k_terms

      ALLOCATE (k_terms(SIZE(c_left_outer, 1), SIZE(c_right_outer, 1), 4))
      CALL topology_reference( &
         c_left_outer, c_right_outer, metric_oo, density, &
         transpose_left=.FALSE., transpose_right=.FALSE., k_block=k_terms(:, :, 1))
      CALL topology_reference( &
         c_left_outer, c_right_inner, metric_oi, density, &
         transpose_left=.FALSE., transpose_right=.TRUE., k_block=k_terms(:, :, 2))
      CALL topology_reference( &
         c_left_inner, c_right_outer, metric_io, density, &
         transpose_left=.TRUE., transpose_right=.FALSE., k_block=k_terms(:, :, 3))
      CALL topology_reference( &
         c_left_inner, c_right_inner, metric_ii, density, &
         transpose_left=.TRUE., transpose_right=.TRUE., k_block=k_terms(:, :, 4))

   END SUBROUTINE four_topology_reference

! **************************************************************************************************
!> \brief Applies the finite-k loop3 kernels to each physical auxiliary-center topology separately.
!> \param c_left_outer C(mu,sigma,P) with P on the outer/target center
!> \param c_left_inner C(sigma,mu,P) with P on the inner/contracted center
!> \param c_right_outer C(nu,lambda,Q) with Q on the outer/output center
!> \param c_right_inner C(lambda,nu,Q) with Q on the inner/contracted center
!> \param metric_oo outer-outer auxiliary metric
!> \param metric_oi outer-inner auxiliary metric
!> \param metric_io inner-outer auxiliary metric
!> \param metric_ii inner-inner auxiliary metric
!> \param density D(sigma,lambda)
!> \param k_terms K(mu,nu,topology), ordered as OO, OI, IO, II
! **************************************************************************************************
   SUBROUTINE four_topology_optimized( &
      c_left_outer, c_left_inner, c_right_outer, c_right_inner, &
      metric_oo, metric_oi, metric_io, metric_ii, density, k_terms)

      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: c_left_outer, c_left_inner, &
                                                            c_right_outer, c_right_inner
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: metric_oo, metric_oi, metric_io, &
                                                            metric_ii, density
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: k_terms

      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: block, left_inner_flat, s_inner_flat, &
                                                            s_outer_flat, transformed_flat
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: s_inner, s_outer, transformed
      TYPE(hfx_lri_contract_workspace_type)              :: local_workspace

      ALLOCATE (k_terms(SIZE(c_left_outer, 1), SIZE(c_right_outer, 1), 4))

      CALL hfx_lri_loop3_density_c_add( &
         c_right_outer, density, s_outer, 1.0_dp, .TRUE., .TRUE.)
      CALL hfx_lri_loop3_pack_aux(s_outer, s_outer_flat)
      CALL hfx_lri_loop3_central_contract( &
         c_left_outer, metric_oo, s_outer_flat, block, local_workspace)
      k_terms(:, :, 1) = block

      CALL hfx_lri_loop3_density_c_add( &
         c_right_inner, density, s_inner, 1.0_dp, .TRUE., .FALSE.)
      CALL hfx_lri_loop3_pack_aux(s_inner, s_inner_flat)
      CALL hfx_lri_loop3_central_contract( &
         c_left_outer, metric_oi, s_inner_flat, block, local_workspace)
      k_terms(:, :, 2) = block

      CALL hfx_lri_loop3_pack_aux(c_left_inner, left_inner_flat)
      CALL hfx_lri_loop3_aux_transform( &
         s_outer, metric_io, transformed, 1.0_dp, .TRUE.)
      CALL hfx_lri_loop3_pack_aux(transformed, transformed_flat)
      CALL hfx_lri_loop3_remote_contract(left_inner_flat, transformed_flat, block)
      k_terms(:, :, 3) = block

      CALL hfx_lri_loop3_aux_transform( &
         s_inner, metric_ii, transformed, 1.0_dp, .TRUE.)
      CALL hfx_lri_loop3_pack_aux(transformed, transformed_flat)
      CALL hfx_lri_loop3_remote_contract(left_inner_flat, transformed_flat, block)
      k_terms(:, :, 4) = block

      DEALLOCATE (block, left_inner_flat, s_inner, s_inner_flat, s_outer, s_outer_flat, &
                  transformed, transformed_flat)

   END SUBROUTINE four_topology_optimized

! **************************************************************************************************
!> \brief Evaluates one topology by direct summation, without packed views or BLAS.
!> \param c_left left fitting coefficients in stored orientation
!> \param c_right right fitting coefficients in stored orientation
!> \param metric auxiliary-center interaction
!> \param density contracted AO density block
!> \param transpose_left whether the first two left-C indices are stored as (sigma,mu)
!> \param transpose_right whether the first two right-C indices are stored as (lambda,nu)
!> \param k_block resulting exchange block K(mu,nu)
! **************************************************************************************************
   SUBROUTINE topology_reference(c_left, c_right, metric, density, &
                                 transpose_left, transpose_right, k_block)

      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: c_left, c_right
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: metric, density
      LOGICAL, INTENT(IN)                                :: transpose_left, transpose_right
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: k_block

      INTEGER                                            :: ilambda, imu, inu, ip, iq, isigma, &
                                                            nlambda, nmu, nnu, nsigma
      REAL(dp)                                           :: c_left_value, c_right_value

      IF (transpose_left) THEN
         nmu = SIZE(c_left, 2)
         nsigma = SIZE(c_left, 1)
      ELSE
         nmu = SIZE(c_left, 1)
         nsigma = SIZE(c_left, 2)
      END IF
      IF (transpose_right) THEN
         nnu = SIZE(c_right, 2)
         nlambda = SIZE(c_right, 1)
      ELSE
         nnu = SIZE(c_right, 1)
         nlambda = SIZE(c_right, 2)
      END IF
      IF (ANY(SHAPE(metric) /= [SIZE(c_left, 3), SIZE(c_right, 3)]) .OR. &
          ANY(SHAPE(density) /= [nsigma, nlambda]) .OR. &
          ANY(SHAPE(k_block) /= [nmu, nnu])) THEN
         ERROR STOP "inconsistent dimensions in pair-local RI topology reference"
      END IF

      k_block(:, :) = 0.0_dp
      DO inu = 1, nnu
         DO imu = 1, nmu
            DO iq = 1, SIZE(c_right, 3)
               DO ip = 1, SIZE(c_left, 3)
                  DO ilambda = 1, nlambda
                     DO isigma = 1, nsigma
                        IF (transpose_left) THEN
                           c_left_value = c_left(isigma, imu, ip)
                        ELSE
                           c_left_value = c_left(imu, isigma, ip)
                        END IF
                        IF (transpose_right) THEN
                           c_right_value = c_right(ilambda, inu, iq)
                        ELSE
                           c_right_value = c_right(inu, ilambda, iq)
                        END IF
                        k_block(imu, inu) = k_block(imu, inu) + &
                                            c_left_value*metric(ip, iq)*c_right_value* &
                                            density(isigma, ilambda)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE topology_reference

END PROGRAM hfx_lri_contract_unittest
