!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Block communication and bounded OpenMP accumulation for pair-local RI.
! **************************************************************************************************
MODULE hfx_lri_parallel
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_para_env_type
   USE omp_lib,                         ONLY: omp_destroy_lock,&
                                              omp_init_lock,&
                                              omp_lock_kind,&
                                              omp_set_lock,&
                                              omp_unset_lock
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hfx_lri_block_type, hfx_lri_buffer_type, hfx_lri_parallel_type, &
             hfx_lri_buffer_init, hfx_lri_buffer_clear, hfx_lri_buffer_release, &
             hfx_lri_buffer_add, hfx_lri_exchange_blocks, hfx_lri_fetch_blocks, &
             hfx_lri_gather_metadata, hfx_lri_assign_tasks

   ! Limit each sender's numerical payload per collective. A single larger block is indivisible.
   INTEGER, PARAMETER :: transfer_elements = 131072

   TYPE hfx_lri_block_type
      REAL(dp), ALLOCATABLE, DIMENSION(:, :) :: values
   END TYPE hfx_lri_block_type

   TYPE hfx_lri_buffer_type
      TYPE(hfx_lri_block_type), ALLOCATABLE, DIMENSION(:) :: blocks
      INTEGER(omp_lock_kind), ALLOCATABLE, DIMENSION(:) :: locks
   END TYPE hfx_lri_buffer_type

   TYPE hfx_lri_parallel_type
      ! The QS environment owns para_env; this object owns only the split communicator.
      TYPE(mp_para_env_type), POINTER :: para_env => NULL(), subgroup => NULL()
      INTEGER :: ngroups = 1, group_size = 1
      INTEGER, ALLOCATABLE, DIMENSION(:) :: owner, tasks_s, tasks_t
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: density_keys
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: need_d, metric_present
   END TYPE hfx_lri_parallel_type

   INTERFACE hfx_lri_gather_metadata
      MODULE PROCEDURE gather_integer_metadata, gather_real_metadata
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Allocate one sparse rank buffer, with a lock for each possible output block.
!> \param buffer ...
!> \param nblocks ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_buffer_init(buffer, nblocks)
      TYPE(hfx_lri_buffer_type), INTENT(OUT)             :: buffer
      INTEGER, INTENT(IN)                                :: nblocks

      INTEGER                                            :: i

      ALLOCATE (buffer%blocks(nblocks), buffer%locks(nblocks))
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(buffer, nblocks) PRIVATE(i)
      DO i = 1, nblocks
         CALL omp_init_lock(buffer%locks(i))
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE hfx_lri_buffer_init

! **************************************************************************************************
!> \brief Clear numerical payloads between task batches; retain the index space and locks.
!> \param buffer ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_buffer_clear(buffer)
      TYPE(hfx_lri_buffer_type), INTENT(INOUT)           :: buffer

      INTEGER                                            :: i

      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(buffer) PRIVATE(i)
      DO i = 1, SIZE(buffer%blocks)
         IF (ALLOCATED(buffer%blocks(i)%values)) DEALLOCATE (buffer%blocks(i)%values)
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE hfx_lri_buffer_clear

! **************************************************************************************************
!> \brief Release the locks after all workers have completed.
!> \param buffer ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_buffer_release(buffer)
      TYPE(hfx_lri_buffer_type), INTENT(INOUT)           :: buffer

      INTEGER                                            :: i

      IF (.NOT. ALLOCATED(buffer%locks)) RETURN
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(buffer) PRIVATE(i)
      DO i = 1, SIZE(buffer%locks)
         CALL omp_destroy_lock(buffer%locks(i))
      END DO
      !$OMP END PARALLEL DO
      DEALLOCATE (buffer%locks, buffer%blocks)
   END SUBROUTINE hfx_lri_buffer_release

! **************************************************************************************************
!> \brief Add a completed BLAS block. Locks protect only allocation and addition, never contraction.
!> \param buffer ...
!> \param key ...
!> \param values ...
!> \param factor ...
!> \param transposed ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_buffer_add(buffer, key, values, factor, transposed)
      TYPE(hfx_lri_buffer_type), INTENT(INOUT)           :: buffer
      INTEGER, INTENT(IN)                                :: key
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: values
      REAL(dp), INTENT(IN)                               :: factor
      LOGICAL, INTENT(IN)                                :: transposed

      INTEGER                                            :: ncol, nrow

      CPASSERT(key >= 1 .AND. key <= SIZE(buffer%blocks))
      nrow = SIZE(values, 1)
      ncol = SIZE(values, 2)
      IF (transposed) THEN
         nrow = SIZE(values, 2)
         ncol = SIZE(values, 1)
      END IF
      CALL omp_set_lock(buffer%locks(key))
      IF (.NOT. ALLOCATED(buffer%blocks(key)%values)) THEN
         ALLOCATE (buffer%blocks(key)%values(nrow, ncol), SOURCE=0.0_dp)
      END IF
      CPASSERT(ALL(SHAPE(buffer%blocks(key)%values) == [nrow, ncol]))
      IF (transposed) THEN
         buffer%blocks(key)%values(:, :) = buffer%blocks(key)%values + factor*TRANSPOSE(values)
      ELSE
         buffer%blocks(key)%values(:, :) = buffer%blocks(key)%values + factor*values
      END IF
      CALL omp_unset_lock(buffer%locks(key))
   END SUBROUTINE hfx_lri_buffer_add

! **************************************************************************************************
!> \brief Gather small integer metadata; numerical tensors use selective block communication.
!> \param local ...
!> \param global ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE gather_integer_metadata(local, global, para_env)
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: local
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: global
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      INTEGER                                            :: i, width
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: counts, offsets, packed, received

      width = SIZE(local, 1)
      CPASSERT(width > 0)
      ALLOCATE (counts(para_env%num_pe), offsets(para_env%num_pe))
      CALL para_env%allgather(SIZE(local), counts)
      CALL block_offsets(counts, offsets)
      ALLOCATE (packed(SIZE(local)), received(SUM(counts)))
      packed(:) = RESHAPE(local, [SIZE(local)])
      CALL para_env%allgatherv(packed, received, counts, offsets)
      ALLOCATE (global(width, SIZE(received)/width))
      global(:, :) = RESHAPE(received, SHAPE(global))
      DO i = 1, SIZE(counts)
         CPASSERT(MODULO(counts(i), width) == 0)
      END DO
   END SUBROUTINE gather_integer_metadata

! **************************************************************************************************
!> \brief Gather positions accompanying the integer pair metadata.
!> \param local ...
!> \param global ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE gather_real_metadata(local, global, para_env)
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: local
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: global
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      INTEGER                                            :: width
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: counts, offsets
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: packed, received

      width = SIZE(local, 1)
      CPASSERT(width > 0)
      ALLOCATE (counts(para_env%num_pe), offsets(para_env%num_pe))
      CALL para_env%allgather(SIZE(local), counts)
      CALL block_offsets(counts, offsets)
      ALLOCATE (packed(SIZE(local)), received(SUM(counts)))
      packed(:) = RESHAPE(local, [SIZE(local)])
      CALL para_env%allgatherv(packed, received, counts, offsets)
      ALLOCATE (global(width, SIZE(received)/width))
      global(:, :) = RESHAPE(received, SHAPE(global))
   END SUBROUTINE gather_real_metadata

! **************************************************************************************************
!> \brief Assign related panels to subgroups, then balance their tasks within each subgroup.
!> \param cost nonnegative task work estimates; zero denotes an absent task
!> \param panel one-based panel identifier for each task
!> \param ngroups ...
!> \param nranks ...
!> \param owner zero-based execution rank for each task, or -1 for an absent task
! **************************************************************************************************
   SUBROUTINE hfx_lri_assign_tasks(cost, panel, ngroups, nranks, owner)
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: cost
      INTEGER, DIMENSION(:), INTENT(IN)                  :: panel
      INTEGER, INTENT(IN)                                :: ngroups, nranks
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: owner

      INTEGER                                            :: first, group, group_size, i, last, p
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: panel_group
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: group_cost, panel_cost, rank_cost

      CPASSERT(ngroups > 0 .AND. nranks >= ngroups)
      CPASSERT(MODULO(nranks, ngroups) == 0)
      CPASSERT(SIZE(cost) == SIZE(panel))
      CPASSERT(ALL(panel > 0) .AND. ALL(cost >= 0.0_dp))
      group_size = nranks/ngroups
      ALLOCATE (owner(SIZE(cost)), SOURCE=-1)
      IF (SIZE(cost) == 0) RETURN
      ALLOCATE (panel_group(MAXVAL(panel)), SOURCE=0)
      ALLOCATE (panel_cost(SIZE(panel_group)), group_cost(ngroups), rank_cost(nranks), SOURCE=0.0_dp)
      DO i = 1, SIZE(cost)
         panel_cost(panel(i)) = panel_cost(panel(i)) + cost(i)
      END DO
      DO p = 1, SIZE(panel_group)
         group = MINLOC(group_cost, 1)
         panel_group(p) = group
         group_cost(group) = group_cost(group) + panel_cost(p)
      END DO
      DO i = 1, SIZE(cost)
         IF (cost(i) == 0.0_dp) CYCLE
         group = panel_group(panel(i))
         first = (group - 1)*group_size + 1
         last = first + group_size - 1
         p = first + MINLOC(rank_cost(first:last), 1) - 1
         owner(i) = p - 1
         rank_cost(p) = rank_cost(p) + cost(i)
      END DO
   END SUBROUTINE hfx_lri_assign_tasks

! **************************************************************************************************
!> \brief Route and sum sparse partial blocks at their unique owners. No numerical filtering here.
!> \param source ...
!> \param owner ...
!> \param destination already allocated index space; existing values are accumulated
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_exchange_blocks(source, owner, destination, para_env)
      TYPE(hfx_lri_block_type), DIMENSION(:), INTENT(IN) :: source
      INTEGER, DIMENSION(:), INTENT(IN)                  :: owner
      TYPE(hfx_lri_block_type), DIMENSION(:), &
         INTENT(INOUT)                                   :: destination
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      INTEGER                                            :: i, n
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: keys, ranks

      CPASSERT(SIZE(source) == SIZE(owner) .AND. SIZE(source) == SIZE(destination))
      n = 0
      DO i = 1, SIZE(source)
         IF (ALLOCATED(source(i)%values)) n = n + 1
      END DO
      ALLOCATE (keys(n), ranks(n))
      n = 0
      DO i = 1, SIZE(source)
         IF (.NOT. ALLOCATED(source(i)%values)) CYCLE
         n = n + 1
         keys(n) = i
         ranks(n) = owner(i)
      END DO
      CALL send_blocks(source, keys, ranks, destination, para_env, .TRUE.)
   END SUBROUTINE hfx_lri_exchange_blocks

! **************************************************************************************************
!> \brief Import only requested blocks from their canonical owners. All ranks participate.
!> \param source ...
!> \param owner ...
!> \param requested ...
!> \param destination requested slots must initially be empty
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_fetch_blocks(source, owner, requested, destination, para_env)
      TYPE(hfx_lri_block_type), DIMENSION(:), INTENT(IN) :: source
      INTEGER, DIMENSION(:), INTENT(IN)                  :: owner
      LOGICAL, DIMENSION(:), INTENT(IN)                  :: requested
      TYPE(hfx_lri_block_type), DIMENSION(:), &
         INTENT(INOUT)                                   :: destination
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      INTEGER                                            :: i, k, p
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cursor, keys, ranks, rcount, received, &
                                                            roff, scount, soff

      CPASSERT(SIZE(source) == SIZE(owner) .AND. SIZE(source) == SIZE(requested))
      ALLOCATE (scount(para_env%num_pe), rcount(para_env%num_pe), &
                soff(para_env%num_pe), roff(para_env%num_pe), cursor(para_env%num_pe))
      scount = 0
      DO i = 1, SIZE(source)
         IF (.NOT. requested(i)) CYCLE
         CPASSERT(owner(i) >= 0 .AND. owner(i) < para_env%num_pe)
         scount(owner(i) + 1) = scount(owner(i) + 1) + 1
      END DO
      CALL para_env%alltoall(scount, rcount, 1)
      CALL block_offsets(scount, soff)
      CALL block_offsets(rcount, roff)
      ALLOCATE (keys(SUM(scount)), received(SUM(rcount)), ranks(SUM(rcount)))
      cursor(:) = soff
      DO i = 1, SIZE(source)
         IF (.NOT. requested(i)) CYCLE
         p = owner(i) + 1
         cursor(p) = cursor(p) + 1
         keys(cursor(p)) = i
      END DO
      CALL para_env%alltoall(keys, scount, soff, received, rcount, roff)
      DO p = 1, para_env%num_pe
         DO k = roff(p) + 1, roff(p) + rcount(p)
            ranks(k) = p - 1
         END DO
      END DO
      CALL send_blocks(source, received, ranks, destination, para_env, .FALSE.)
   END SUBROUTINE hfx_lri_fetch_blocks

! **************************************************************************************************
!> \brief Exchange packed headers and values in bounded batches, including empty senders.
!> \param source ...
!> \param keys ...
!> \param ranks ...
!> \param destination ...
!> \param para_env ...
!> \param add ...
! **************************************************************************************************
   SUBROUTINE send_blocks(source, keys, ranks, destination, para_env, add)
      TYPE(hfx_lri_block_type), DIMENSION(:), INTENT(IN) :: source
      INTEGER, DIMENSION(:), INTENT(IN)                  :: keys, ranks
      TYPE(hfx_lri_block_type), DIMENSION(:), &
         INTENT(INOUT)                                   :: destination
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      LOGICAL, INTENT(IN)                                :: add

      INTEGER                                            :: first, h, i, key, last, more, n, nc, nr, &
                                                            p, v
      INTEGER(int_8)                                     :: total
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: hc, headers, incoming, rh, rho, rv, rvo, &
                                                            sh, sho, sv, svo, vc
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: received, values

      CPASSERT(SIZE(keys) == SIZE(ranks))
      n = para_env%num_pe
      ALLOCATE (sh(n), rh(n), sv(n), rv(n), sho(n), rho(n), svo(n), rvo(n), hc(n), vc(n))
      first = 1
      DO
         more = MERGE(1, 0, first <= SIZE(keys))
         CALL para_env%max(more)
         IF (more == 0) EXIT
         sh = 0
         sv = 0
         total = 0_int_8
         last = first - 1
         DO i = first, SIZE(keys)
            key = keys(i)
            CPASSERT(key >= 1 .AND. key <= SIZE(source))
            CPASSERT(ALLOCATED(source(key)%values))
            n = SIZE(source(key)%values)
            CPASSERT(n > 0)
            IF (i > first .AND. total + INT(n, int_8) > transfer_elements) EXIT
            CPASSERT(total + INT(n, int_8) <= HUGE(0))
            p = ranks(i) + 1
            CPASSERT(p >= 1 .AND. p <= para_env%num_pe)
            sh(p) = sh(p) + 3
            sv(p) = sv(p) + n
            total = total + INT(n, int_8)
            last = i
         END DO
         CALL para_env%alltoall(sh, rh, 1)
         CALL para_env%alltoall(sv, rv, 1)
         CALL block_offsets(sh, sho)
         CALL block_offsets(rh, rho)
         CALL block_offsets(sv, svo)
         CALL block_offsets(rv, rvo)
         ALLOCATE (headers(SUM(sh)), incoming(SUM(rh)), values(SUM(sv)), received(SUM(rv)))
         hc(:) = sho
         vc(:) = svo
         DO i = first, last
            key = keys(i)
            p = ranks(i) + 1
            nr = SIZE(source(key)%values, 1)
            nc = SIZE(source(key)%values, 2)
            headers(hc(p) + 1:hc(p) + 3) = [key, nr, nc]
            values(vc(p) + 1:vc(p) + nr*nc) = RESHAPE(source(key)%values, [nr*nc])
            hc(p) = hc(p) + 3
            vc(p) = vc(p) + nr*nc
         END DO
         CALL para_env%alltoall(headers, sh, sho, incoming, rh, rho)
         CALL para_env%alltoall(values, sv, svo, received, rv, rvo)
         ! Different source ranks may contribute to the same key. Merge in source-rank order.
         DO p = 1, para_env%num_pe
            v = rvo(p)
            DO h = rho(p) + 1, rho(p) + rh(p), 3
               key = incoming(h)
               nr = incoming(h + 1)
               nc = incoming(h + 2)
               CPASSERT(key >= 1 .AND. key <= SIZE(destination))
               CPASSERT(nr > 0 .AND. nc > 0)
               CPASSERT(INT(nr, int_8)*INT(nc, int_8) <= INT(rv(p), int_8))
               n = nr*nc
               CPASSERT(v + n <= rvo(p) + rv(p))
               IF (ALLOCATED(destination(key)%values)) THEN
                  CPASSERT(add)
                  CPASSERT(ALL(SHAPE(destination(key)%values) == [nr, nc]))
                  destination(key)%values(:, :) = destination(key)%values + RESHAPE(received(v + 1:v + n), [nr, nc])
               ELSE
                  ALLOCATE (destination(key)%values(nr, nc))
                  destination(key)%values(:, :) = RESHAPE(received(v + 1:v + n), [nr, nc])
               END IF
               v = v + n
            END DO
            CPASSERT(v == rvo(p) + rv(p))
         END DO
         DEALLOCATE (headers, incoming, values, received)
         first = last + 1
      END DO
   END SUBROUTINE send_blocks

! **************************************************************************************************
!> \brief Form checked zero-based MPI displacements.
!> \param counts ...
!> \param offsets ...
! **************************************************************************************************
   SUBROUTINE block_offsets(counts, offsets)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: counts
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: offsets

      INTEGER                                            :: i
      INTEGER(int_8)                                     :: total

      CPASSERT(SIZE(counts) == SIZE(offsets))
      total = 0_int_8
      DO i = 1, SIZE(counts)
         CPASSERT(counts(i) >= 0)
         offsets(i) = INT(total)
         total = total + INT(counts(i), int_8)
         IF (total > INT(HUGE(0), int_8)) CPABORT("LRI-HFXk MPI buffer exceeds the integer count range")
      END DO
   END SUBROUTINE block_offsets

END MODULE hfx_lri_parallel
