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
             hfx_lri_loop3_pack_aux, hfx_lri_loop3_remote_contract, &
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
!> \param transpose_centered contract the remote rather than the central C index
!> \param centered_remote_flat ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_density_c_add(centered, density, RESULT, factor, first, &
                                          transpose_centered, centered_remote_flat)

      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: centered
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: density
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: result
      REAL(dp), INTENT(IN)                               :: factor
      LOGICAL, INTENT(IN)                                :: first, transpose_centered
      REAL(dp), DIMENSION(:, :), INTENT(IN), OPTIONAL    :: centered_remote_flat

      INTEGER                                            :: iaux, ncontract, nfirst, nout
      REAL(dp)                                           :: beta

      nfirst = SIZE(density, 1)
      ncontract = SIZE(density, 2)
      IF (transpose_centered) THEN
         nout = SIZE(centered, 1)
         IF (SIZE(centered, 2) /= ncontract) THEN
            CPABORT("Inconsistent remote C dimension in LRI loop3 density contraction")
         END IF
      ELSE
         nout = SIZE(centered, 2)
         IF (SIZE(centered, 1) /= ncontract) THEN
            CPABORT("Inconsistent central C dimension in LRI loop3 density contraction")
         END IF
      END IF

      IF (ALLOCATED(RESULT)) THEN
         IF (ANY(SHAPE(RESULT) /= [nfirst, nout, SIZE(centered, 3)])) DEALLOCATE (RESULT)
      END IF
      IF (.NOT. ALLOCATED(RESULT)) THEN
         ALLOCATE (RESULT(nfirst, nout, SIZE(centered, 3)), SOURCE=0.0_dp)
      ELSE IF (first) THEN
         RESULT(:, :, :) = 0.0_dp
      END IF
      beta = MERGE(0.0_dp, 1.0_dp, first)
      IF (.NOT. transpose_centered) THEN
         ! C(center,remote,aux) is already a contiguous center-by-(remote*aux) matrix.
         CALL dgemm("N", "N", nfirst, nout*SIZE(centered, 3), ncontract, factor, &
                    density(1, 1), nfirst, centered(1, 1, 1), ncontract, beta, &
                    RESULT(1, 1, 1), nfirst)
      ELSE IF (PRESENT(centered_remote_flat)) THEN
         IF (ANY(SHAPE(centered_remote_flat) /= &
                 [nout*SIZE(centered, 3), ncontract])) THEN
            CPABORT("Inconsistent packed C dimensions in LRI loop3 density contraction")
         END IF
         CALL dgemm("N", "T", nfirst, nout*SIZE(centered, 3), ncontract, factor, &
                    density(1, 1), nfirst, centered_remote_flat(1, 1), &
                    nout*SIZE(centered, 3), beta, RESULT(1, 1, 1), nfirst)
      ELSE
         ! Compatibility path for callers which do not retain the packed C view.
         DO iaux = 1, SIZE(centered, 3)
            CALL dgemm("N", "T", nfirst, nout, ncontract, factor, density, nfirst, &
                       centered(:, :, iaux), nout, beta, RESULT(:, :, iaux), nfirst)
         END DO
      END IF

   END SUBROUTINE hfx_lri_loop3_density_c_add

! **************************************************************************************************
!> \brief Transforms the auxiliary index of a loop3 tensor and optionally accumulates the result.
!> \param source ...
!> \param metric ...
!> \param TARGET ...
!> \param factor ...
!> \param first ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_aux_transform(source, metric, TARGET, factor, first)

      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: source
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: metric
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: target
      REAL(dp), INTENT(IN)                               :: factor
      LOGICAL, INTENT(IN)                                :: first

      INTEGER                                            :: nflat, nleft, nright
      REAL(dp)                                           :: beta

      nflat = SIZE(source, 1)*SIZE(source, 2)
      nright = SIZE(metric, 2)
      nleft = SIZE(metric, 1)
      IF (SIZE(source, 3) /= nright) THEN
         CPABORT("Inconsistent metric dimension in LRI loop3 auxiliary contraction")
      END IF
      IF (ALLOCATED(TARGET)) THEN
         IF (ANY(SHAPE(TARGET) /= [SIZE(source, 1), SIZE(source, 2), nleft])) THEN
            DEALLOCATE (TARGET)
         END IF
      END IF
      IF (.NOT. ALLOCATED(TARGET)) THEN
         ALLOCATE (TARGET(SIZE(source, 1), SIZE(source, 2), nleft), SOURCE=0.0_dp)
      ELSE IF (first) THEN
         TARGET(:, :, :) = 0.0_dp
      END IF
      beta = MERGE(0.0_dp, 1.0_dp, first)
      CALL dgemm("N", "T", nflat, nleft, nright, factor, source(1, 1, 1), nflat, &
                 metric(1, 1), nleft, beta, TARGET(1, 1, 1), nflat)

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
!> \brief Contracts C(center,remote,aux) with U(center,output,aux), leaving remote/output.
!> \param centered_remote_flat ...
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
! **************************************************************************************************
   SUBROUTINE hfx_lri_loop3_central_contract(centered, metric, intermediate_aux_flat, block, &
                                             workspace)

      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: centered
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: metric, intermediate_aux_flat
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: block
      TYPE(hfx_lri_contract_workspace_type), &
         INTENT(INOUT)                                   :: workspace

      INTEGER                                            :: naux, ncenter, nintermediate_aux, nout, &
                                                            nremote

      ncenter = SIZE(centered, 1)
      nremote = SIZE(centered, 2)
      naux = SIZE(centered, 3)
      nout = SIZE(intermediate_aux_flat, 2)
      nintermediate_aux = SIZE(metric, 2)
      IF (SIZE(metric, 1) /= naux) THEN
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
      CALL dgemm("N", "N", ncenter*nremote, nintermediate_aux, naux, &
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
