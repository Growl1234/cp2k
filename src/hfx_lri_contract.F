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

   PUBLIC :: hfx_lri_contract_exchange, hfx_lri_contract_exchange_add, &
             hfx_lri_contract_outer_accumulate, hfx_lri_contract_outer_finalize, &
             hfx_lri_contract_outer_prepare, hfx_lri_contract_workspace_type, &
             hfx_lri_density_cell, hfx_lri_product_index, hfx_lri_right_cell

   TYPE hfx_lri_contract_workspace_type
      REAL(dp), ALLOCATABLE, DIMENSION(:, :) :: density_right, left_coeff, right_coeff, &
                                                weighted_left, weighted_right_flat, left_flat, &
                                                left_metric_flat, outer_density_sum
   END TYPE hfx_lri_contract_workspace_type

CONTAINS

! **************************************************************************************************
!> \brief Linear index of the AO product (ia,ib), with the second AO varying fastest.
! **************************************************************************************************
   PURE INTEGER FUNCTION hfx_lri_product_index(ia, ib, nb) RESULT(index)

      INTEGER, INTENT(IN)                                :: ia, ib, nb

      index = (ia - 1)*nb + ib

   END FUNCTION hfx_lri_product_index

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
!> \brief Returns the right fitting-pair cell that produces a requested density cell.
!> \param cell_density cell of D_JL
!> \param cell_output cell of K_IK
!> \param cell_left cell of the left fitting pair I-J
!> \return cell of the right fitting pair K-L
! **************************************************************************************************
   PURE FUNCTION hfx_lri_right_cell(cell_density, cell_output, cell_left) RESULT(cell_right)

      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell_density, cell_output, cell_left
      INTEGER, DIMENSION(3)                              :: cell_right

      cell_right(:) = cell_density(:) - cell_output(:) + cell_left(:)

   END FUNCTION hfx_lri_right_cell

! **************************************************************************************************
!> \brief Contracts C_IJ, V, C_KL and D_JL into K_IK.
!> \param c_left C_IJ(mu,sigma,X)
!> \param c_right C_KL(nu,lambda,Y)
!> \param v_aux V(X,Y)
!> \param density D_JL(sigma,lambda)
!> \param k_block K_IK(mu,nu)
! **************************************************************************************************
   SUBROUTINE hfx_lri_contract_exchange(c_left, c_right, v_aux, density, k_block, &
                                        transpose_left, transpose_right)

      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: c_left, c_right
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: v_aux, density
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: k_block
      LOGICAL, INTENT(IN), OPTIONAL                      :: transpose_left, transpose_right

      INTEGER                                            :: nmu, nnu
      LOGICAL                                            :: left_t, right_t
      TYPE(hfx_lri_contract_workspace_type)              :: workspace

      left_t = .FALSE.
      right_t = .FALSE.
      IF (PRESENT(transpose_left)) left_t = transpose_left
      IF (PRESENT(transpose_right)) right_t = transpose_right
      nmu = MERGE(SIZE(c_left, 2), SIZE(c_left, 1), left_t)
      nnu = MERGE(SIZE(c_right, 2), SIZE(c_right, 1), right_t)

      ALLOCATE (k_block(nmu, nnu), SOURCE=0.0_dp)
      CALL hfx_lri_contract_exchange_add(c_left, c_right, v_aux, density, k_block, workspace, &
                                         left_t, right_t)

   END SUBROUTINE hfx_lri_contract_exchange

! **************************************************************************************************
!> \brief Adds one C V C^T D contraction to an existing K block while reusing all temporaries.
! **************************************************************************************************
   SUBROUTINE hfx_lri_contract_exchange_add(c_left, c_right, v_aux, density, k_block, workspace, &
                                            transpose_left, transpose_right)

      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: c_left, c_right
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: v_aux, density
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: k_block
      TYPE(hfx_lri_contract_workspace_type), &
         INTENT(INOUT)                                   :: workspace
      LOGICAL, INTENT(IN), OPTIONAL                      :: transpose_left, transpose_right

      INTEGER                                            :: iaux, iy, nlambda, nmu, nnu, nsigma, &
                                                            nsigma_nnu, right_col
      LOGICAL                                            :: left_t, right_t

      left_t = .FALSE.
      right_t = .FALSE.
      IF (PRESENT(transpose_left)) left_t = transpose_left
      IF (PRESENT(transpose_right)) right_t = transpose_right
      nmu = MERGE(SIZE(c_left, 2), SIZE(c_left, 1), left_t)
      nsigma = MERGE(SIZE(c_left, 1), SIZE(c_left, 2), left_t)
      nnu = MERGE(SIZE(c_right, 2), SIZE(c_right, 1), right_t)
      nlambda = MERGE(SIZE(c_right, 1), SIZE(c_right, 2), right_t)

      IF (nsigma /= SIZE(density, 1) .OR. nlambda /= SIZE(density, 2)) THEN
         CPABORT("Inconsistent AO dimensions in pair-local RI exchange contraction")
      END IF
      IF (SIZE(c_left, 3) /= SIZE(v_aux, 1) .OR. &
          SIZE(c_right, 3) /= SIZE(v_aux, 2)) THEN
         CPABORT("Inconsistent auxiliary dimensions in pair-local RI exchange contraction")
      END IF
      IF (ANY(SHAPE(k_block) /= [nmu, nnu])) THEN
         CPABORT("Inconsistent output dimensions in pair-local RI exchange contraction")
      END IF

      ! K(mu,nu) = sum[sigma,lambda,X,Y]
      !              C_IJ(mu,sigma,X) V_XY C_KL(nu,lambda,Y) D_JL(sigma,lambda).
      ! Pack both C tensors so the complete contraction uses three GEMMs instead of one small
      ! MATMUL per auxiliary function. Buffer reshaping/permutation is linear in tensor size and
      ! all buffers survive across pair contractions in the caller.
      nsigma_nnu = nsigma*nnu
      IF (ALLOCATED(workspace%right_coeff)) THEN
         IF (ANY(SHAPE(workspace%right_coeff) /= [nlambda, nnu*SIZE(c_right, 3)])) THEN
            DEALLOCATE (workspace%right_coeff)
         END IF
      END IF
      IF (.NOT. ALLOCATED(workspace%right_coeff)) THEN
         ALLOCATE (workspace%right_coeff(nlambda, nnu*SIZE(c_right, 3)))
      END IF
      DO iy = 1, SIZE(c_right, 3)
         right_col = (iy - 1)*nnu
         IF (right_t) THEN
            workspace%right_coeff(:, right_col + 1:right_col + nnu) = c_right(:, :, iy)
         ELSE
            workspace%right_coeff(:, right_col + 1:right_col + nnu) = &
               TRANSPOSE(c_right(:, :, iy))
         END IF
      END DO

      IF (ALLOCATED(workspace%density_right)) THEN
         IF (ANY(SHAPE(workspace%density_right) /= [nsigma, nnu*SIZE(c_right, 3)])) THEN
            DEALLOCATE (workspace%density_right)
         END IF
      END IF
      IF (.NOT. ALLOCATED(workspace%density_right)) THEN
         ALLOCATE (workspace%density_right(nsigma, nnu*SIZE(c_right, 3)))
      END IF
      IF (ALLOCATED(workspace%weighted_right_flat)) THEN
         IF (ANY(SHAPE(workspace%weighted_right_flat) /= [nsigma_nnu, SIZE(c_left, 3)])) THEN
            DEALLOCATE (workspace%weighted_right_flat)
         END IF
      END IF
      IF (.NOT. ALLOCATED(workspace%weighted_right_flat)) THEN
         ALLOCATE (workspace%weighted_right_flat(nsigma_nnu, SIZE(c_left, 3)))
      END IF
      IF (ALLOCATED(workspace%weighted_left)) THEN
         IF (ANY(SHAPE(workspace%weighted_left) /= [nsigma*SIZE(c_left, 3), nnu])) THEN
            DEALLOCATE (workspace%weighted_left)
         END IF
      END IF
      IF (.NOT. ALLOCATED(workspace%weighted_left)) THEN
         ALLOCATE (workspace%weighted_left(nsigma*SIZE(c_left, 3), nnu))
      END IF
      IF (ALLOCATED(workspace%left_coeff)) THEN
         IF (ANY(SHAPE(workspace%left_coeff) /= [nmu, nsigma*SIZE(c_left, 3)])) THEN
            DEALLOCATE (workspace%left_coeff)
         END IF
      END IF
      IF (.NOT. ALLOCATED(workspace%left_coeff)) THEN
         ALLOCATE (workspace%left_coeff(nmu, nsigma*SIZE(c_left, 3)))
      END IF
      DO iaux = 1, SIZE(c_left, 3)
         IF (left_t) THEN
            workspace%left_coeff(:, (iaux - 1)*nsigma + 1:iaux*nsigma) = &
               TRANSPOSE(c_left(:, :, iaux))
         ELSE
            workspace%left_coeff(:, (iaux - 1)*nsigma + 1:iaux*nsigma) = c_left(:, :, iaux)
         END IF
      END DO

      CALL dgemm("N", "N", nsigma, nnu*SIZE(c_right, 3), nlambda, 1.0_dp, &
                 density, nsigma, workspace%right_coeff, nlambda, 0.0_dp, &
                 workspace%density_right, nsigma)
      ! density_right is stored as (sigma,nu,Y). Its contiguous storage is therefore already the
      ! (sigma*nu,Y) matrix needed below; passing the combined leading dimension avoids a reshape
      ! and a complete tensor copy.
      CALL dgemm("N", "T", nsigma_nnu, SIZE(c_left, 3), SIZE(c_right, 3), 1.0_dp, &
                 workspace%density_right, nsigma_nnu, v_aux, SIZE(v_aux, 1), 0.0_dp, &
                 workspace%weighted_right_flat, nsigma_nnu)
      DO iaux = 1, SIZE(c_left, 3)
         DO iy = 1, nnu
            workspace%weighted_left((iaux - 1)*nsigma + 1:iaux*nsigma, iy) = &
               workspace%weighted_right_flat((iy - 1)*nsigma + 1:iy*nsigma, iaux)
         END DO
      END DO
      CALL dgemm("N", "N", nmu, nnu, nsigma*SIZE(c_left, 3), 1.0_dp, &
                 workspace%left_coeff, nmu, workspace%weighted_left, nsigma*SIZE(c_left, 3), &
                 1.0_dp, k_block, nmu)

   END SUBROUTINE hfx_lri_contract_exchange_add

! **************************************************************************************************
!> \brief Precontracts all auxiliary functions of the left pair with one fixed outer-center metric.
! **************************************************************************************************
   SUBROUTINE hfx_lri_contract_outer_prepare(c_left, v_outer, workspace, transpose_left)

      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: c_left
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: v_outer
      TYPE(hfx_lri_contract_workspace_type), &
         INTENT(INOUT)                                   :: workspace
      LOGICAL, INTENT(IN), OPTIONAL                      :: transpose_left

      INTEGER                                            :: iaux, naux, nmu, nq, nsigma
      LOGICAL                                            :: left_t

      left_t = .FALSE.
      IF (PRESENT(transpose_left)) left_t = transpose_left
      nmu = MERGE(SIZE(c_left, 2), SIZE(c_left, 1), left_t)
      nsigma = MERGE(SIZE(c_left, 1), SIZE(c_left, 2), left_t)
      naux = SIZE(c_left, 3)
      nq = SIZE(v_outer, 2)
      IF (SIZE(v_outer, 1) /= naux) THEN
         CPABORT("Inconsistent auxiliary dimensions in pair-local RI outer preparation")
      END IF

      IF (ALLOCATED(workspace%left_flat)) THEN
         IF (ANY(SHAPE(workspace%left_flat) /= [nmu*nsigma, naux])) THEN
            DEALLOCATE (workspace%left_flat)
         END IF
      END IF
      IF (.NOT. ALLOCATED(workspace%left_flat)) THEN
         ALLOCATE (workspace%left_flat(nmu*nsigma, naux))
      END IF
      IF (ALLOCATED(workspace%left_metric_flat)) THEN
         IF (ANY(SHAPE(workspace%left_metric_flat) /= [nmu*nsigma, nq])) THEN
            DEALLOCATE (workspace%left_metric_flat)
         END IF
      END IF
      IF (.NOT. ALLOCATED(workspace%left_metric_flat)) THEN
         ALLOCATE (workspace%left_metric_flat(nmu*nsigma, nq))
      END IF

      DO iaux = 1, naux
         IF (left_t) THEN
            workspace%left_flat(:, iaux) = &
               RESHAPE(TRANSPOSE(c_left(:, :, iaux)), [nmu*nsigma])
         ELSE
            workspace%left_flat(:, iaux) = RESHAPE(c_left(:, :, iaux), [nmu*nsigma])
         END IF
      END DO
      CALL dgemm("N", "N", nmu*nsigma, nq, naux, 1.0_dp, workspace%left_flat, &
                 nmu*nsigma, v_outer, naux, 0.0_dp, workspace%left_metric_flat, nmu*nsigma)

   END SUBROUTINE hfx_lri_contract_outer_prepare

! **************************************************************************************************
!> \brief Adds D times the outer-center slice of a right C block to a reusable target/left sum.
! **************************************************************************************************
   SUBROUTINE hfx_lri_contract_outer_accumulate(c_right, density, workspace, first, transpose_right)

      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: c_right
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: density
      TYPE(hfx_lri_contract_workspace_type), &
         INTENT(INOUT)                                   :: workspace
      LOGICAL, INTENT(IN)                                :: first
      LOGICAL, INTENT(IN), OPTIONAL                      :: transpose_right

      INTEGER                                            :: iq, nlambda, nnu, nq, right_col
      LOGICAL                                            :: right_t
      REAL(dp)                                           :: beta

      right_t = .FALSE.
      IF (PRESENT(transpose_right)) right_t = transpose_right
      nnu = MERGE(SIZE(c_right, 2), SIZE(c_right, 1), right_t)
      nlambda = MERGE(SIZE(c_right, 1), SIZE(c_right, 2), right_t)
      nq = SIZE(c_right, 3)
      IF (SIZE(density, 2) /= nlambda) THEN
         CPABORT("Inconsistent density dimensions in pair-local RI outer accumulation")
      END IF

      IF (ALLOCATED(workspace%right_coeff)) THEN
         IF (ANY(SHAPE(workspace%right_coeff) /= [nlambda, nnu*nq])) THEN
            DEALLOCATE (workspace%right_coeff)
         END IF
      END IF
      IF (.NOT. ALLOCATED(workspace%right_coeff)) THEN
         ALLOCATE (workspace%right_coeff(nlambda, nnu*nq))
      END IF
      DO iq = 1, nq
         right_col = (iq - 1)*nnu
         IF (right_t) THEN
            workspace%right_coeff(:, right_col + 1:right_col + nnu) = c_right(:, :, iq)
         ELSE
            workspace%right_coeff(:, right_col + 1:right_col + nnu) = &
               TRANSPOSE(c_right(:, :, iq))
         END IF
      END DO

      IF (ALLOCATED(workspace%outer_density_sum)) THEN
         IF (ANY(SHAPE(workspace%outer_density_sum) /= [SIZE(density, 1), nnu*nq])) THEN
            DEALLOCATE (workspace%outer_density_sum)
         END IF
      END IF
      IF (.NOT. ALLOCATED(workspace%outer_density_sum)) THEN
         ALLOCATE (workspace%outer_density_sum(SIZE(density, 1), nnu*nq))
      END IF
      beta = MERGE(0.0_dp, 1.0_dp, first)
      CALL dgemm("N", "N", SIZE(density, 1), nnu*nq, nlambda, 1.0_dp, density, &
                 SIZE(density, 1), workspace%right_coeff, nlambda, beta, &
                 workspace%outer_density_sum, SIZE(density, 1))

   END SUBROUTINE hfx_lri_contract_outer_accumulate

! **************************************************************************************************
!> \brief Completes one fixed target/left outer-center contribution with a single GEMM.
! **************************************************************************************************
   SUBROUTINE hfx_lri_contract_outer_finalize(k_block, workspace)

      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: k_block
      TYPE(hfx_lri_contract_workspace_type), &
         INTENT(INOUT)                                   :: workspace

      INTEGER                                            :: imu_sigma, inu, iq, nmu, nnu, nq, nsigma

      IF (.NOT. ALLOCATED(workspace%left_metric_flat) .OR. &
          .NOT. ALLOCATED(workspace%outer_density_sum)) THEN
         CPABORT("Missing pair-local RI outer-contraction workspace")
      END IF
      nmu = SIZE(k_block, 1)
      nnu = SIZE(k_block, 2)
      IF (MOD(SIZE(workspace%left_metric_flat, 1), nmu) /= 0) THEN
         CPABORT("Inconsistent left dimensions in pair-local RI outer finalization")
      END IF
      nsigma = SIZE(workspace%left_metric_flat, 1)/nmu
      nq = SIZE(workspace%left_metric_flat, 2)
      IF (ANY(SHAPE(workspace%outer_density_sum) /= [nsigma, nnu*nq])) THEN
         CPABORT("Inconsistent right dimensions in pair-local RI outer finalization")
      END IF

      IF (ALLOCATED(workspace%weighted_left)) THEN
         IF (ANY(SHAPE(workspace%weighted_left) /= [nsigma*nq, nnu])) THEN
            DEALLOCATE (workspace%weighted_left)
         END IF
      END IF
      IF (.NOT. ALLOCATED(workspace%weighted_left)) THEN
         ALLOCATE (workspace%weighted_left(nsigma*nq, nnu))
      END IF
      DO iq = 1, nq
         DO inu = 1, nnu
            imu_sigma = (iq - 1)*nsigma
            workspace%weighted_left(imu_sigma + 1:imu_sigma + nsigma, inu) = &
               workspace%outer_density_sum(:, (iq - 1)*nnu + inu)
         END DO
      END DO
      ! left_metric_flat is stored as (mu*sigma,Q). The same contiguous data are the
      ! (mu,sigma*Q) matrix consumed by this GEMM.
      CALL dgemm("N", "N", nmu, nnu, nsigma*nq, 1.0_dp, workspace%left_metric_flat, &
                 nmu, workspace%weighted_left, nsigma*nq, 1.0_dp, k_block, nmu)

   END SUBROUTINE hfx_lri_contract_outer_finalize

END MODULE hfx_lri_contract
