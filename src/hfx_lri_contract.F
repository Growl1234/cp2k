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

   PUBLIC :: hfx_lri_contract_workspace_type, hfx_lri_loop3_aux_transform, &
             hfx_lri_loop3_central_contract, hfx_lri_loop3_density_c_add, &
             hfx_lri_loop3_density_remote_add, hfx_lri_loop3_pack_aux, &
             hfx_lri_loop3_unpack_aux, hfx_lri_loop3_remote_contract, &
             hfx_lri_product_index

   TYPE hfx_lri_contract_workspace_type
      REAL(dp), ALLOCATABLE, DIMENSION(:, :) :: loop3_s_flat, loop3_t_flat, loop3_u_flat, &
                                                loop3_block
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :) :: loop3_s, loop3_t, loop3_u, loop3_w
   END TYPE hfx_lri_contract_workspace_type

CONTAINS

! **************************************************************************************************
!> \brief Linear index of the AO product (ia,ib), with the second AO varying fastest.
!> \param ia ...
!> \param ib ...
!> \param nb ...
!> \return ...
! **************************************************************************************************
   PURE INTEGER FUNCTION hfx_lri_product_index(ia, ib, nb) RESULT(index)

      INTEGER, INTENT(IN)                                :: ia, ib, nb

      index = (ia - 1)*nb + ib

   END FUNCTION hfx_lri_product_index

! **************************************************************************************************
!> \brief Adds D times one atom-centred C tensor to a loop3 intermediate.
!> \param centered C(center,remote,aux)
!> \param density density(first,contracted)
!> \param RESULT ...
!> \param factor boundary-phase factor
!> \param first overwrite rather than add to result
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_density_c_add(centered, density, RESULT, factor, first)

      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: centered
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: density
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: result
      REAL(dp), INTENT(IN)                               :: factor
      LOGICAL, INTENT(IN)                                :: first

      INTEGER                                            :: ncontract, nfirst, nout
      REAL(dp)                                           :: beta

      nfirst = SIZE(density, 1)
      ncontract = SIZE(density, 2)
      nout = SIZE(centered, 2)
      IF (SIZE(centered, 1) /= ncontract) THEN
         CPABORT("Inconsistent central C dimension in LRI loop3 density contraction")
      END IF

      beta = MERGE(0.0_dp, 1.0_dp, first)
      IF (ALLOCATED(RESULT)) THEN
         IF (ANY(SHAPE(RESULT) /= [nfirst, nout, SIZE(centered, 3)])) DEALLOCATE (RESULT)
      END IF
      IF (.NOT. ALLOCATED(RESULT)) THEN
         ALLOCATE (RESULT(nfirst, nout, SIZE(centered, 3)))
         beta = 0.0_dp
      END IF
      ! C(center,remote,aux) is a contiguous center-by-(remote*aux) matrix.
      CALL dgemm("N", "N", nfirst, nout*SIZE(centered, 3), ncontract, factor, &
                 density(1, 1), nfirst, centered(1, 1, 1), ncontract, beta, &
                 RESULT(1, 1, 1), nfirst)

   END SUBROUTINE hfx_lri_loop3_density_c_add

! **************************************************************************************************
!> \brief Adds D times C directly from the persistent (center*aux,remote) layout.
!> \param centered_remote_flat ...
!> \param ncenter ...
!> \param density ...
!> \param RESULT ...
!> \param factor ...
!> \param first ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_density_remote_add(centered_remote_flat, ncenter, density, RESULT, factor, first)
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: centered_remote_flat, density
      INTEGER, INTENT(IN)                                :: ncenter
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :), INTENT(INOUT) :: result
      REAL(dp), INTENT(IN)                               :: factor
      LOGICAL, INTENT(IN)                                :: first

      INTEGER                                            :: naux, nfirst, nremote
      REAL(dp)                                           :: beta

      CPASSERT(ncenter > 0)
      CPASSERT(MODULO(SIZE(centered_remote_flat, 1), ncenter) == 0)
      naux = SIZE(centered_remote_flat, 1)/ncenter
      nremote = SIZE(centered_remote_flat, 2)
      nfirst = SIZE(density, 1)
      CPASSERT(SIZE(density, 2) == nremote)
      beta = MERGE(0.0_dp, 1.0_dp, first)
      IF (ALLOCATED(RESULT)) THEN
         IF (ANY(SHAPE(RESULT) /= [nfirst, ncenter, naux])) DEALLOCATE (RESULT)
      END IF
      IF (.NOT. ALLOCATED(RESULT)) THEN
         ALLOCATE (RESULT(nfirst, ncenter, naux))
         beta = 0.0_dp
      END IF
      CALL dgemm("N", "T", nfirst, ncenter*naux, nremote, factor, &
                 density(1, 1), nfirst, centered_remote_flat(1, 1), &
                 ncenter*naux, beta, RESULT(1, 1, 1), nfirst)
   END SUBROUTINE hfx_lri_loop3_density_remote_add

! **************************************************************************************************
!> \brief Transforms the auxiliary index of a loop3 tensor and optionally accumulates the result.
!> \param source ...
!> \param metric ...
!> \param TARGET ...
!> \param factor ...
!> \param first ...
!> \param transposed metric is stored in its reverse orientation
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_aux_transform(source, metric, TARGET, factor, first, transposed)

      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: source
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: metric
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: target
      REAL(dp), INTENT(IN)                               :: factor
      LOGICAL, INTENT(IN)                                :: first
      LOGICAL, INTENT(IN), OPTIONAL                      :: transposed

      INTEGER                                            :: nflat, nleft, nright
      CHARACTER                                         :: operation
      REAL(dp)                                           :: beta

      nflat = SIZE(source, 1)*SIZE(source, 2)
      nright = SIZE(metric, 2)
      nleft = SIZE(metric, 1)
      operation = "T"
      IF (PRESENT(transposed)) THEN
         IF (transposed) THEN
            nright = SIZE(metric, 1)
            nleft = SIZE(metric, 2)
            operation = "N"
         END IF
      END IF
      IF (SIZE(source, 3) /= nright) THEN
         CPABORT("Inconsistent metric dimension in LRI loop3 auxiliary contraction")
      END IF
      beta = MERGE(0.0_dp, 1.0_dp, first)
      IF (ALLOCATED(TARGET)) THEN
         IF (ANY(SHAPE(TARGET) /= [SIZE(source, 1), SIZE(source, 2), nleft])) THEN
            DEALLOCATE (TARGET)
         END IF
      END IF
      IF (.NOT. ALLOCATED(TARGET)) THEN
         ALLOCATE (TARGET(SIZE(source, 1), SIZE(source, 2), nleft))
         beta = 0.0_dp
      END IF
      CALL dgemm("N", operation, nflat, nleft, nright, factor, source(1, 1, 1), nflat, &
                 metric(1, 1), SIZE(metric, 1), beta, TARGET(1, 1, 1), nflat)

   END SUBROUTINE hfx_lri_loop3_aux_transform

! **************************************************************************************************
!> \brief Packs tensor(first,second,aux) as matrix(first*aux,second).
!> \param source ...
!> \param packed ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_pack_aux(source, packed)

      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: source
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: packed

      INTEGER                                            :: iaux, naux, nfirst, nsecond

      nfirst = SIZE(source, 1)
      nsecond = SIZE(source, 2)
      naux = SIZE(source, 3)
      IF (ALLOCATED(packed)) THEN
         IF (ANY(SHAPE(packed) /= [nfirst*naux, nsecond])) DEALLOCATE (packed)
      END IF
      IF (.NOT. ALLOCATED(packed)) ALLOCATE (packed(nfirst*naux, nsecond))
      DO iaux = 1, naux
         packed((iaux - 1)*nfirst + 1:iaux*nfirst, :) = source(:, :, iaux)
      END DO

   END SUBROUTINE hfx_lri_loop3_pack_aux

! **************************************************************************************************
!> \brief Reorders one C block in reusable scratch, without retaining a second copy of all C.
!> \param packed matrix(first*aux,second)
!> \param nfirst ...
!> \param unpacked tensor(first,second,aux)
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_unpack_aux(packed, nfirst, unpacked)
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: packed
      INTEGER, INTENT(IN)                                :: nfirst
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :), INTENT(INOUT) :: unpacked

      INTEGER                                            :: iaux, naux, nsecond

      CPASSERT(nfirst > 0)
      CPASSERT(MODULO(SIZE(packed, 1), nfirst) == 0)
      naux = SIZE(packed, 1)/nfirst
      nsecond = SIZE(packed, 2)
      IF (ALLOCATED(unpacked)) THEN
         IF (ANY(SHAPE(unpacked) /= [nfirst, nsecond, naux])) DEALLOCATE (unpacked)
      END IF
      IF (.NOT. ALLOCATED(unpacked)) ALLOCATE (unpacked(nfirst, nsecond, naux))
      DO iaux = 1, naux
         unpacked(:, :, iaux) = packed((iaux - 1)*nfirst + 1:iaux*nfirst, :)
      END DO
   END SUBROUTINE hfx_lri_loop3_unpack_aux

! **************************************************************************************************
!> \brief Contracts C(center,remote,aux) with U(center,output,aux), leaving remote/output.
!> \param centered_remote_flat packed C(center*aux,remote)
!> \param intermediate_aux_flat ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_remote_contract(centered_remote_flat, intermediate_aux_flat, block)

      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: centered_remote_flat, &
                                                            intermediate_aux_flat
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: block

      INTEGER                                            :: ncontract, nout, nremote

      ncontract = SIZE(centered_remote_flat, 1)
      nremote = SIZE(centered_remote_flat, 2)
      nout = SIZE(intermediate_aux_flat, 2)
      IF (SIZE(intermediate_aux_flat, 1) /= ncontract) THEN
         CPABORT("Inconsistent packed intermediate in LRI loop3 remote contraction")
      END IF
      IF (ALLOCATED(block)) THEN
         IF (ANY(SHAPE(block) /= [nremote, nout])) DEALLOCATE (block)
      END IF
      IF (.NOT. ALLOCATED(block)) ALLOCATE (block(nremote, nout))
      CALL dgemm("T", "N", nremote, nout, ncontract, 1.0_dp, &
                 centered_remote_flat(1, 1), ncontract, &
                 intermediate_aux_flat(1, 1), ncontract, &
                 0.0_dp, block(1, 1), nremote)

   END SUBROUTINE hfx_lri_loop3_remote_contract

! **************************************************************************************************
!> \brief Contracts C(center,remote,P) V(P,Q) S(remote,output,Q), leaving center/output.
!> \param centered ...
!> \param metric ...
!> \param intermediate_aux_flat ...
!> \param block ...
!> \param workspace ...
!> \param transposed metric is stored in its reverse orientation
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_central_contract(centered, metric, intermediate_aux_flat, block, &
                                             workspace, transposed)

      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: centered
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: metric, intermediate_aux_flat
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: block
      TYPE(hfx_lri_contract_workspace_type), &
         INTENT(INOUT)                                   :: workspace
      LOGICAL, INTENT(IN), OPTIONAL                      :: transposed

      INTEGER                                            :: naux, ncenter, nintermediate_aux, nout, &
                                                            nremote
      CHARACTER                                         :: operation

      ncenter = SIZE(centered, 1)
      nremote = SIZE(centered, 2)
      naux = SIZE(centered, 3)
      nout = SIZE(intermediate_aux_flat, 2)
      nintermediate_aux = SIZE(metric, 2)
      operation = "N"
      IF (PRESENT(transposed)) THEN
         IF (transposed) THEN
            nintermediate_aux = SIZE(metric, 1)
            operation = "T"
         END IF
      END IF
      IF (SIZE(metric, MERGE(1, 2, operation == "N")) /= naux) THEN
         CPABORT("Inconsistent auxiliary dimensions in LRI loop3 central contraction")
      END IF
      IF (ANY(SHAPE(intermediate_aux_flat) /= &
              [nremote*nintermediate_aux, nout])) THEN
         CPABORT("Inconsistent packed intermediate in LRI loop3 central contraction")
      END IF
      IF (ALLOCATED(workspace%loop3_w)) THEN
         IF (ANY(SHAPE(workspace%loop3_w) /= &
                 [ncenter, nremote, nintermediate_aux])) THEN
            DEALLOCATE (workspace%loop3_w)
         END IF
      END IF
      IF (.NOT. ALLOCATED(workspace%loop3_w)) THEN
         ALLOCATE (workspace%loop3_w(ncenter, nremote, nintermediate_aux))
      END IF
      CALL dgemm("N", operation, ncenter*nremote, nintermediate_aux, naux, &
                 1.0_dp, centered(1, 1, 1), ncenter*nremote, metric(1, 1), &
                 SIZE(metric, 1), 0.0_dp, workspace%loop3_w(1, 1, 1), &
                 ncenter*nremote)
      IF (ALLOCATED(block)) THEN
         IF (ANY(SHAPE(block) /= [ncenter, nout])) DEALLOCATE (block)
      END IF
      IF (.NOT. ALLOCATED(block)) ALLOCATE (block(ncenter, nout))
      ! The transformed C tensor is contiguous as (center,remote*intermediate_aux).
      CALL dgemm("N", "N", ncenter, nout, nremote*nintermediate_aux, &
                 1.0_dp, workspace%loop3_w(1, 1, 1), ncenter, &
                 intermediate_aux_flat(1, 1), nremote*nintermediate_aux, &
                 0.0_dp, block(1, 1), ncenter)

   END SUBROUTINE hfx_lri_loop3_central_contract

END MODULE hfx_lri_contract
