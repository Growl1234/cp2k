!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Finite-lattice quotient helpers used by pair-local RI HFX.
! **************************************************************************************************
MODULE hfx_lri_lattice
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hfx_lri_lattice_block_from_index, hfx_lri_lattice_block_index, &
             hfx_lri_lattice_cell_from_index, hfx_lri_lattice_cell_index, &
             hfx_lri_lattice_displacement_cell, hfx_lri_lattice_pair_representative, &
             hfx_lri_lattice_reduce, &
             hfx_lri_reciprocal_diagnostics

CONTAINS

! **************************************************************************************************
!> \brief Resolves a physical displacement into the lattice cell of two canonical atom images.
!> \param displacement physical vector from the left center to the right center
!> \param position_left canonical position of the left atom
!> \param position_right canonical position of the right atom
!> \param h_inv inverse cell matrix
!> \param cell integer lattice translation of the right atom
!> \param residual distance from the nearest integer in scaled coordinates
! **************************************************************************************************
   PURE SUBROUTINE hfx_lri_lattice_displacement_cell( &
      displacement, position_left, position_right, h_inv, cell, residual)

      REAL(dp), DIMENSION(3), INTENT(IN)                 :: displacement, position_left, &
                                                            position_right
      REAL(dp), DIMENSION(3, 3), INTENT(IN)              :: h_inv
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: cell
      REAL(dp), INTENT(OUT)                              :: residual

      REAL(dp), DIMENSION(3)                             :: scaled

      scaled(:) = MATMUL(h_inv, displacement - (position_right - position_left))
      cell(:) = NINT(scaled(:))
      residual = MAXVAL(ABS(scaled - REAL(cell, dp)))

   END SUBROUTINE hfx_lri_lattice_displacement_cell

! **************************************************************************************************
!> \brief Selects one representative of a Hermitian AO-pair translation orbit.
!> \param iatom first atom
!> \param jatom second atom
!> \param cell lattice translation of the second atom
!> \return true for the representative retained by the LRI pair fitter
!>
!> CP2K's symmetric neighbor-list convention chooses one atom ordering but, for a diagonal atom
!> block, contains both R and -R.  LRI explicitly reconstructs the reverse coefficient, hence
!> retaining both diagonal translations would add the same physical C coefficient twice.
! **************************************************************************************************
   PURE LOGICAL FUNCTION hfx_lri_lattice_pair_representative(iatom, jatom, cell) RESULT(keep)

      INTEGER, INTENT(IN)                                :: iatom, jatom
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell

      INTEGER                                            :: idir

      keep = .TRUE.
      IF (iatom /= jatom) RETURN

      ! Keep the origin and the lexicographically positive member of every {R,-R} pair.
      DO idir = 1, 3
         IF (cell(idir) == 0) CYCLE
         keep = cell(idir) > 0
         RETURN
      END DO

   END FUNCTION hfx_lri_lattice_pair_representative

! **************************************************************************************************
!> \brief Reduces a cell modulo a regular k mesh and returns its real boundary phase.
!> \param cell physical lattice cell
!> \param grid regular k-mesh dimensions
!> \param twist boundary phase for crossing one mesh period in each direction
!> \param reduced representative in the finite lattice quotient
!> \param factor accumulated boundary phase
! **************************************************************************************************
   PURE SUBROUTINE hfx_lri_lattice_reduce(cell, grid, twist, reduced, factor)

      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell, grid, twist
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: reduced
      REAL(dp), INTENT(OUT), OPTIONAL                    :: factor

      INTEGER                                            :: idir, periods
      REAL(dp)                                           :: phase

      phase = 1.0_dp
      DO idir = 1, 3
         reduced(idir) = MODULO(cell(idir), grid(idir))
         periods = (cell(idir) - reduced(idir))/grid(idir)
         IF (twist(idir) == -1 .AND. MODULO(periods, 2) /= 0) phase = -phase
      END DO
      IF (PRESENT(factor)) factor = phase

   END SUBROUTINE hfx_lri_lattice_reduce

! **************************************************************************************************
!> \brief Returns the linear index of a reduced cell, with the first direction varying fastest.
!> \param cell reduced lattice cell
!> \param grid regular k-mesh dimensions
!> \return one-based cell index
! **************************************************************************************************
   PURE INTEGER FUNCTION hfx_lri_lattice_cell_index(cell, grid) RESULT(index)

      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell, grid

      index = 1 + cell(1) + grid(1)*(cell(2) + grid(2)*cell(3))

   END FUNCTION hfx_lri_lattice_cell_index

! **************************************************************************************************
!> \brief Decodes a reduced-cell linear index.
!> \param index one-based cell index
!> \param grid regular k-mesh dimensions
!> \param cell reduced lattice cell
! **************************************************************************************************
   PURE SUBROUTINE hfx_lri_lattice_cell_from_index(index, grid, cell)

      INTEGER, INTENT(IN)                                :: index
      INTEGER, DIMENSION(3), INTENT(IN)                  :: grid
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: cell

      INTEGER                                            :: work

      work = index - 1
      cell(1) = MODULO(work, grid(1))
      work = work/grid(1)
      cell(2) = MODULO(work, grid(2))
      cell(3) = work/grid(2)

   END SUBROUTINE hfx_lri_lattice_cell_from_index

! **************************************************************************************************
!> \brief Returns the direct index of an oriented atom-pair block on the lattice quotient.
!> \param iatom first atom
!> \param jatom second atom
!> \param cell_index reduced-cell index
!> \param natom number of atoms
!> \param ncell number of reduced cells
!> \return one-based block index
! **************************************************************************************************
   PURE INTEGER FUNCTION hfx_lri_lattice_block_index(iatom, jatom, cell_index, natom, ncell) &
      RESULT(index)

      INTEGER, INTENT(IN)                                :: iatom, jatom, cell_index, natom, ncell

      index = cell_index + ncell*((jatom - 1) + natom*(iatom - 1))

   END FUNCTION hfx_lri_lattice_block_index

! **************************************************************************************************
!> \brief Decodes an oriented atom-pair/lattice block index.
!> \param index one-based block index
!> \param natom number of atoms
!> \param ncell number of reduced cells
!> \param iatom first atom
!> \param jatom second atom
!> \param cell_index reduced-cell index
! **************************************************************************************************
   PURE SUBROUTINE hfx_lri_lattice_block_from_index(index, natom, ncell, iatom, jatom, cell_index)

      INTEGER, INTENT(IN)                                :: index, natom, ncell
      INTEGER, INTENT(OUT)                               :: iatom, jatom, cell_index

      INTEGER                                            :: work

      cell_index = MODULO(index - 1, ncell) + 1
      work = (index - 1)/ncell
      jatom = MODULO(work, natom) + 1
      iatom = work/natom + 1

   END SUBROUTINE hfx_lri_lattice_block_from_index

! **************************************************************************************************
!> \brief Reconstructs two real twisted-quotient matrices on the complete reciprocal mesh and
!>        reports Hermiticity, spectral bounds, and their reciprocal-space trace.
!> \param density real-space density on the twisted quotient
!> \param kernel real-space exchange kernel on the same quotient
!> \param grid regular reciprocal-mesh dimensions
!> \param twist boundary phase for crossing one mesh period in each direction
!> \param density_hermitian largest relative Hermiticity residual of D(k)
!> \param kernel_hermitian largest relative Hermiticity residual of K(k)
!> \param density_eig_min smallest eigenvalue of D(k) over the complete mesh
!> \param density_eig_max largest eigenvalue of D(k) over the complete mesh
!> \param kernel_eig_min smallest eigenvalue of K(k) over the complete mesh
!> \param kernel_eig_max largest eigenvalue of K(k) over the complete mesh
!> \param reciprocal_trace mesh average of Tr[D(k) K(k)]
! **************************************************************************************************
   SUBROUTINE hfx_lri_reciprocal_diagnostics( &
      density, kernel, grid, twist, density_hermitian, kernel_hermitian, &
      density_eig_min, density_eig_max, kernel_eig_min, kernel_eig_max, &
      reciprocal_trace)

      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: density, kernel
      INTEGER, DIMENSION(3), INTENT(IN)                  :: grid, twist
      REAL(dp), INTENT(OUT) :: density_hermitian, kernel_hermitian, density_eig_min, &
         density_eig_max, kernel_eig_min, kernel_eig_max, reciprocal_trace

      COMPLEX(dp)                                        :: phase
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:)             :: work
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:, :)          :: density_k, kernel_k, matrix_work
      INTEGER                                            :: cell_index, idir, info, k_index, &
                                                            matrix_size, ncell
      INTEGER, DIMENSION(3)                              :: cell, k_cell
      REAL(dp)                                           :: angle, density_scale, kernel_scale
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: eigenvalues, rwork
      REAL(dp), DIMENSION(3)                             :: k_vector

      IF (ANY(grid <= 0) .OR. ANY(ABS(twist) /= 1)) THEN
         CPABORT("Invalid lattice metadata in LRI reciprocal-space validation")
      END IF
      IF (SIZE(density, 1) /= SIZE(density, 2) .OR. &
          ANY(SHAPE(kernel) /= SHAPE(density))) THEN
         CPABORT("Inconsistent matrices in LRI reciprocal-space validation")
      END IF
      ncell = PRODUCT(grid)
      IF (SIZE(density, 3) /= ncell) THEN
         CPABORT("Inconsistent cell count in LRI reciprocal-space validation")
      END IF

      matrix_size = SIZE(density, 1)
      ALLOCATE (density_k(matrix_size, matrix_size), kernel_k(matrix_size, matrix_size), &
                matrix_work(matrix_size, matrix_size), eigenvalues(matrix_size), &
                work(MAX(1, 2*matrix_size - 1)), rwork(MAX(1, 3*matrix_size - 2)))
      density_hermitian = 0.0_dp
      kernel_hermitian = 0.0_dp
      density_eig_min = HUGE(1.0_dp)
      density_eig_max = -HUGE(1.0_dp)
      kernel_eig_min = HUGE(1.0_dp)
      kernel_eig_max = -HUGE(1.0_dp)
      reciprocal_trace = 0.0_dp

      DO k_index = 1, ncell
         CALL hfx_lri_lattice_cell_from_index(k_index, grid, k_cell)
         DO idir = 1, 3
            k_vector(idir) = (REAL(k_cell(idir), dp) + &
                              MERGE(0.5_dp, 0.0_dp, twist(idir) == -1))/ &
                             REAL(grid(idir), dp)
         END DO
         density_k(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         kernel_k(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO cell_index = 1, ncell
            CALL hfx_lri_lattice_cell_from_index(cell_index, grid, cell)
            angle = twopi*DOT_PRODUCT(k_vector, REAL(cell, dp))
            phase = CMPLX(COS(angle), SIN(angle), KIND=dp)
            density_k(:, :) = density_k + phase*density(:, :, cell_index)
            kernel_k(:, :) = kernel_k + phase*kernel(:, :, cell_index)
         END DO

         density_scale = MAX(TINY(1.0_dp), MAXVAL(ABS(density_k)))
         kernel_scale = MAX(TINY(1.0_dp), MAXVAL(ABS(kernel_k)))
         density_hermitian = MAX( &
                             density_hermitian, &
                             MAXVAL(ABS(density_k - CONJG(TRANSPOSE(density_k))))/ &
                             density_scale)
         kernel_hermitian = MAX( &
                            kernel_hermitian, &
                            MAXVAL(ABS(kernel_k - CONJG(TRANSPOSE(kernel_k))))/ &
                            kernel_scale)
         reciprocal_trace = reciprocal_trace + &
                            REAL(SUM(CONJG(density_k)*kernel_k), KIND=dp)

         matrix_work(:, :) = 0.5_dp*(density_k + CONJG(TRANSPOSE(density_k)))
         CALL zheev('N', 'U', matrix_size, matrix_work, matrix_size, eigenvalues, &
                    work, SIZE(work), rwork, info)
         IF (info /= 0) CPABORT("D(k) diagonalization failed in LRI validation")
         density_eig_min = MIN(density_eig_min, eigenvalues(1))
         density_eig_max = MAX(density_eig_max, eigenvalues(matrix_size))

         matrix_work(:, :) = 0.5_dp*(kernel_k + CONJG(TRANSPOSE(kernel_k)))
         CALL zheev('N', 'U', matrix_size, matrix_work, matrix_size, eigenvalues, &
                    work, SIZE(work), rwork, info)
         IF (info /= 0) CPABORT("K(k) diagonalization failed in LRI validation")
         kernel_eig_min = MIN(kernel_eig_min, eigenvalues(1))
         kernel_eig_max = MAX(kernel_eig_max, eigenvalues(matrix_size))
      END DO
      reciprocal_trace = reciprocal_trace/REAL(ncell, dp)

   END SUBROUTINE hfx_lri_reciprocal_diagnostics

END MODULE hfx_lri_lattice
