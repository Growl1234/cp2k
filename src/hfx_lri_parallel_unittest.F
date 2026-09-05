!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Exercise sparse LRI transport, empty ranks, subgroup sums and concurrent block updates.
! **************************************************************************************************
PROGRAM hfx_lri_parallel_unittest
   USE hfx_lri_parallel,                ONLY: &
        hfx_lri_assign_tasks, hfx_lri_block_type, hfx_lri_buffer_add, hfx_lri_buffer_clear, &
        hfx_lri_buffer_init, hfx_lri_buffer_release, hfx_lri_buffer_type, hfx_lri_exchange_blocks, &
        hfx_lri_fetch_blocks, hfx_lri_gather_metadata
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: default_output_unit
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_world_finalize,&
                                              mp_world_init
#include "./base/base_uses.f90"

   IMPLICIT NONE
   TYPE(mp_para_env_type) :: world, subgroup
   TYPE(hfx_lri_buffer_type) :: buffer
   TYPE(hfx_lri_block_type), ALLOCATABLE :: owned(:), fetched(:), partial(:), grouped(:)
   INTEGER, PARAMETER :: nblocks = 9, nupdates = 128
   INTEGER :: key, iter, rank, p, ngroups, group_size, source_rank, nr, nc
   INTEGER :: owners(nblocks), subgroup_owners(nblocks), panels(8)
   INTEGER, ALLOCATABLE :: task_owner(:), local_metadata(:, :), metadata(:, :)
   REAL(dp) :: cost(8), value(2, 3), expected, error
   REAL(dp), ALLOCATABLE :: large(:, :)
   LOGICAL :: requested(nblocks)

   CALL mp_world_init(world)
   rank = world%mepos
   owners = [(MODULO(key - 1, world%num_pe), key=1, nblocks)]
   CALL hfx_lri_buffer_init(buffer, nblocks)
   ALLOCATE (owned(nblocks), fetched(nblocks))
   ! The last rank sends nothing in a multi-rank run; ranks with no owned blocks still participate.
   IF (rank < MAX(1, world%num_pe - 1)) THEN
      value = REAL(rank + 1, dp)
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(buffer, value) PRIVATE(iter)
      DO iter = 1, nupdates
         CALL hfx_lri_buffer_add(buffer, 1, value, 1.0_dp, .FALSE.)
         CALL hfx_lri_buffer_add(buffer, 2, value, -0.5_dp, .TRUE.)
      END DO
      !$OMP END PARALLEL DO
      ! Force several communication batches, including an indivisible block larger than the cap.
      DO key = 3, nblocks
         nr = 193 + key
         nc = 173 + 2*key
         IF (key == nblocks) nr = 801
         ALLOCATE (large(nr, nc), SOURCE=REAL((rank + 1)*key, dp))
         CALL hfx_lri_buffer_add(buffer, key, large, 1.0_dp, .FALSE.)
         DEALLOCATE (large)
      END DO
   END IF
   CALL hfx_lri_exchange_blocks(buffer%blocks, owners, owned, world)
   p = MAX(1, world%num_pe - 1)
   expected = REAL(p*(p + 1)/2, dp)
   error = 0.0_dp
   DO key = 1, nblocks
      CPASSERT(ALLOCATED(owned(key)%values) .EQV. (owners(key) == rank))
      IF (owners(key) /= rank) CYCLE
      SELECT CASE (key)
      CASE (1)
         CPASSERT(ALL(SHAPE(owned(key)%values) == [2, 3]))
         error = MAX(error, MAXVAL(ABS(owned(key)%values - nupdates*expected)))
      CASE (2)
         CPASSERT(ALL(SHAPE(owned(key)%values) == [3, 2]))
         error = MAX(error, MAXVAL(ABS(owned(key)%values + 0.5_dp*nupdates*expected)))
      CASE DEFAULT
         error = MAX(error, MAXVAL(ABS(owned(key)%values - REAL(key, dp)*expected)))
      END SELECT
   END DO
   CALL world%max(error)
   CPASSERT(error < 1.0E-12_dp)

   requested = [(MODULO(key + rank, 2) == 0, key=1, nblocks)]
   CALL hfx_lri_fetch_blocks(owned, owners, requested, fetched, world)
   DO key = 1, nblocks
      CPASSERT(ALLOCATED(fetched(key)%values) .EQV. requested(key))
      IF (.NOT. requested(key)) CYCLE
      SELECT CASE (key)
      CASE (1)
         error = MAX(error, MAXVAL(ABS(fetched(key)%values - nupdates*expected)))
      CASE (2)
         error = MAX(error, MAXVAL(ABS(fetched(key)%values + 0.5_dp*nupdates*expected)))
      CASE DEFAULT
         error = MAX(error, MAXVAL(ABS(fetched(key)%values - REAL(key, dp)*expected)))
      END SELECT
   END DO
   CALL world%max(error)
   CPASSERT(error < 1.0E-12_dp)
   DEALLOCATE (fetched)

   ! Every divisor exercises group-local partial reductions and weighted panel assignment.
   panels = [1, 1, 2, 2, 3, 3, 4, 4]
   cost = [1.0_dp, 20.0_dp, 0.0_dp, 3.0_dp, 8.0_dp, 2.0_dp, 0.0_dp, 0.0_dp]
   DO ngroups = 1, world%num_pe
      IF (MODULO(world%num_pe, ngroups) /= 0) CYCLE
      group_size = world%num_pe/ngroups
      CALL subgroup%from_split(world, rank/group_size)
      ALLOCATE (partial(nblocks), grouped(nblocks))
      subgroup_owners = MODULO(owners, group_size)
      CALL hfx_lri_exchange_blocks(buffer%blocks, subgroup_owners, partial, subgroup)
      CALL hfx_lri_exchange_blocks(partial, owners, grouped, world)
      DO key = 1, nblocks
         CPASSERT(ALLOCATED(grouped(key)%values) .EQV. ALLOCATED(owned(key)%values))
         IF (ALLOCATED(owned(key)%values)) THEN
            error = MAXVAL(ABS(grouped(key)%values - owned(key)%values))
            CPASSERT(error < 1.0E-12_dp)
         END IF
      END DO
      DEALLOCATE (partial, grouped)
      CALL subgroup%free()
      CALL hfx_lri_assign_tasks(cost, panels, ngroups, world%num_pe, task_owner)
      CPASSERT(ALL(PACK(task_owner, cost == 0.0_dp) == -1))
      CPASSERT(ALL(PACK(task_owner, cost > 0.0_dp) >= 0))
      CPASSERT(ALL(PACK(task_owner, cost > 0.0_dp) < world%num_pe))
      CPASSERT(task_owner(1)/group_size == task_owner(2)/group_size)
      CPASSERT(task_owner(5)/group_size == task_owner(6)/group_size)
   END DO
   CALL hfx_lri_assign_tasks(cost(:0), panels(:0), 1, world%num_pe, task_owner)
   CPASSERT(SIZE(task_owner) == 0)
   CALL hfx_lri_assign_tasks(0.0_dp*cost, panels, 1, world%num_pe, task_owner)
   CPASSERT(ALL(task_owner == -1))
   ! Variable and zero-length metadata, in source-rank order.
   ALLOCATE (local_metadata(2, rank))
   local_metadata = rank
   CALL hfx_lri_gather_metadata(local_metadata, metadata, world)
   CPASSERT(SIZE(metadata, 2) == world%num_pe*(world%num_pe - 1)/2)
   key = 0
   DO source_rank = 0, world%num_pe - 1
      DO iter = 1, source_rank
         key = key + 1
         CPASSERT(ALL(metadata(:, key) == source_rank))
      END DO
   END DO
   CALL hfx_lri_buffer_clear(buffer)
   DO key = 1, nblocks
      CPASSERT(.NOT. ALLOCATED(buffer%blocks(key)%values))
   END DO
   ! An entirely empty collective must terminate without entering a payload exchange.
   DEALLOCATE (owned)
   ALLOCATE (owned(nblocks))
   CALL hfx_lri_exchange_blocks(buffer%blocks, owners, owned, world)
   CALL hfx_lri_buffer_release(buffer)
   IF (rank == 0) WRITE (default_output_unit, '(A)') 'LRI MPI/OpenMP transport tests passed.'
   CALL mp_world_finalize()
END PROGRAM hfx_lri_parallel_unittest
