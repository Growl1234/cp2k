!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Product-driven Gaussian auxiliary bases for pair-local RI-HFX.
!>        Construct normalized, contracted, same-center ORB products, separately for each allowed L.
!>        Fit their radial parts in a Gaussian space and select products by pivoted Cholesky of their
!>        overlap Gram matrix. The radial L2 error controls basis construction, not the HFX metric.
!>        All m components use the same radial basis. No exchange integrals are evaluated here.
! **************************************************************************************************
MODULE hfx_lri_basis
   USE IEEE_ARITHMETIC,                 ONLY: ieee_is_finite
   USE aux_basis_set,                   ONLY: create_aux_basis
   USE basis_set_types,                 ONLY: gto_basis_set_type,&
                                              sort_gto_basis_set
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE orbital_pointers,                ONLY: init_orbital_pointers
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_lri_basis'
   PUBLIC :: create_hfx_lri_basis

   ! A radial function r**power SUM_i coeff(i) EXP(-zet(i)*r*r).
   TYPE radial_product_type
      INTEGER :: power = 0
      REAL(dp) :: norm2 = 0.0_dp
      REAL(dp), ALLOCATABLE :: zet(:), coeff(:)
   END TYPE radial_product_type

   TYPE radial_basis_type
      INTEGER :: nproduct = 0
      REAL(dp) :: residual = 0.0_dp, absolute_residual = 0.0_dp
      REAL(dp), ALLOCATABLE :: zet(:), coeff(:, :)
   END TYPE radial_basis_type

CONTAINS

! **************************************************************************************************
!> \brief Generate an RI_HFX basis before init_qs_kind_set normalizes ORB and RI_HFX.
!> \param basis Newly allocated basis; explicit RI_HFX bases must bypass this routine.
!> \param orb Uninitialized ORB basis, including its original normalization convention.
!> \param basis_size AUTO_BASIS size: SMALL, MEDIUM, LARGE, HUGE (0 through 3).
!> \param basis_sort Sorting convention used for the other basis sets.
!> \param output_unit Logger unit, negative on non-writing ranks.
! **************************************************************************************************
   SUBROUTINE create_hfx_lri_basis(basis, orb, basis_size, basis_sort, output_unit)
      TYPE(gto_basis_set_type), POINTER                  :: basis, orb
      INTEGER, INTENT(IN)                                :: basis_size, basis_sort, output_unit

      REAL(dp), PARAMETER :: tolerances(0:3) = [1.0E-4_dp, 1.0E-6_dp, 1.0E-8_dp, 1.0E-10_dp]

      CHARACTER(default_string_length)                   :: name
      INTEGER                                            :: handle, i, iset, ish, l, lmax, maxpgf, &
                                                            maxshell, n, nset, nshell
      INTEGER, ALLOCATABLE                               :: angular(:), nl(:, :), npgf(:)
      REAL(dp)                                           :: exponent, tolerance
      REAL(dp), ALLOCATABLE                              :: coeff(:, :, :), zet(:, :)
      TYPE(radial_basis_type), ALLOCATABLE               :: channels(:)
      TYPE(radial_product_type), ALLOCATABLE             :: shells(:)

      CALL timeset('create_hfx_lri_basis', handle)
      CPASSERT(.NOT. ASSOCIATED(basis))
      CPASSERT(ASSOCIATED(orb))
      CPASSERT(basis_size >= 0 .AND. basis_size <= 3)
      IF (orb%nset == 0) CPABORT("Automatic LRI-HFX basis requires a nonempty ORB basis")
      tolerance = tolerances(basis_size)
      nshell = SUM(orb%nshell)
      lmax = 2*MAXVAL(orb%lmax)
      ALLOCATE (shells(nshell), channels(0:lmax))
      i = 0
      DO iset = 1, orb%nset
         DO ish = 1, orb%nshell(iset)
            i = i + 1
            l = orb%l(ish, iset)
            n = orb%npgf(iset)
            shells(i)%power = l
            shells(i)%zet = PACK(orb%zet(1:n, iset), ABS(orb%gcc(1:n, ish, iset)) > 0.0_dp)
            shells(i)%coeff = PACK(orb%gcc(1:n, ish, iset), ABS(orb%gcc(1:n, ish, iset)) > 0.0_dp)
            ! Primitive-dependent factors must match init_orb_basis_set. Constants depending only
            ! on l and the final shell normalization cancel when each radial product is normalized.
            SELECT CASE (orb%norm_type)
            CASE (-1, 2)
               exponent = 0.25_dp*(2*l + 3)
            CASE (1)
               exponent = l + 1.5_dp
            CASE (0, 3)
               exponent = 0.0_dp
            CASE DEFAULT
               CPABORT("Unsupported ORB normalization in automatic LRI-HFX basis")
            END SELECT
            CPASSERT(ALL(shells(i)%zet > 0.0_dp))
            shells(i)%coeff = shells(i)%coeff*shells(i)%zet**exponent
            CALL normalize_radial(shells(i))
         END DO
      END DO

      ! Each rank follows the ordinary replicated basis initialization. Only independent L channels
      ! run in parallel; orbital pointers, the basis container and logging are initialized serially.
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(shells, channels, lmax, tolerance) PRIVATE(l) SCHEDULE(DYNAMIC)
      DO l = 0, lmax
         CALL build_radial_basis(shells, l, tolerance, channels(l))
      END DO
!$OMP END PARALLEL DO

      nset = COUNT(channels%nproduct > 0)
      CPASSERT(nset > 0)
      maxpgf = 0
      maxshell = 0
      DO l = 0, lmax
         IF (channels(l)%nproduct == 0) CYCLE
         maxpgf = MAX(maxpgf, SIZE(channels(l)%zet))
         maxshell = MAX(maxshell, SIZE(channels(l)%coeff, 2))
      END DO
      ALLOCATE (angular(nset), npgf(nset), nl(0:lmax, nset))
      ALLOCATE (zet(maxpgf, nset), coeff(maxpgf, maxshell, nset))
      nl = 0
      zet = 0.0_dp
      coeff = 0.0_dp
      iset = 0
      IF (output_unit > 0) THEN
         WRITE (output_unit, '(A,A)') ' LRI BASIS| ORB products: ', TRIM(orb%name)
         WRITE (output_unit, '(A,ES12.4)') ' LRI BASIS| Maximum relative squared radial error: ', tolerance
         WRITE (output_unit, '(A)') ' LRI BASIS|    L  Products  Primitives  Radial  Max rel. error  Max abs. error'
      END IF
      DO l = 0, lmax
         IF (channels(l)%nproduct == 0) CYCLE
         iset = iset + 1
         angular(iset) = l
         npgf(iset) = SIZE(channels(l)%zet)
         n = SIZE(channels(l)%coeff, 2)
         nl(l, iset) = n
         zet(1:npgf(iset), iset) = channels(l)%zet
         coeff(1:npgf(iset), 1:n, iset) = channels(l)%coeff
         IF (output_unit > 0) THEN
            WRITE (output_unit, '(A,I5,I10,I12,I8,2ES16.6)') ' LRI BASIS|', l, &
               channels(l)%nproduct, npgf(iset), n, channels(l)%residual, channels(l)%absolute_residual
         END IF
      END DO
      CALL init_orbital_pointers(lmax)
      WRITE (name, '(A,I1)') 'LRI-HFX-PRODUCT-', basis_size
      ! These coefficients multiply normalized primitives. Do not apply the heuristic contraction
      ! or diffuse-tail merging used by the general AUTO_BASIS constructor.
      CALL create_aux_basis(basis, name, nset, angular, angular, nl, npgf, zet, gcc=coeff)
      CALL sort_gto_basis_set(basis, basis_sort)
      IF (output_unit > 0) WRITE (output_unit, '(A,I8)') ' LRI BASIS| Spherical auxiliary functions per atom: ', basis%nsgf
      CALL timestop(handle)
   END SUBROUTINE create_hfx_lri_basis

! **************************************************************************************************
!> \brief Integral of r**(power+2) exp(-zeta*r*r), including the radial volume element.
!> \param zeta Positive exponent sum.
!> \param power Sum of radial powers.
!> \return Analytic radial overlap.
! **************************************************************************************************
   ELEMENTAL FUNCTION radial_integral(zeta, power) RESULT(VALUE)
      REAL(dp), INTENT(IN)                               :: zeta
      INTEGER, INTENT(IN)                                :: power
      REAL(dp)                                           :: VALUE

      REAL(dp)                                           :: p

      p = 0.5_dp*(power + 3)
      VALUE = 0.5_dp*GAMMA(p)*EXP(-p*LOG(zeta))
   END FUNCTION radial_integral

! **************************************************************************************************
!> \brief Normalize a radial function without changing relative primitive weights.
!> \param radial Function; norm2 records its squared norm before normalization.
! **************************************************************************************************
   SUBROUTINE normalize_radial(radial)
      TYPE(radial_product_type), INTENT(INOUT)           :: radial

      INTEGER                                            :: i

      radial%norm2 = 0.0_dp
      DO i = 1, SIZE(radial%zet)
         radial%norm2 = radial%norm2 + radial%coeff(i)*SUM(radial%coeff* &
                                                           radial_integral(radial%zet(i) + radial%zet, 2*radial%power))
      END DO
      ! A zero or numerically singular ORB/product must not silently disappear from training.
      CPASSERT(ieee_is_finite(radial%norm2))
      CPASSERT(radial%norm2 > TINY(1.0_dp))
      radial%coeff = radial%coeff/SQRT(radial%norm2)
   END SUBROUTINE normalize_radial

! **************************************************************************************************
!> \brief Form allowed contracted products and fit/compress one complete angular channel.
!> \param shells Normalized ORB radial shells, read only.
!> \param l Angular momentum of the auxiliary functions.
!> \param tolerance Maximum relative squared radial L2 error for each product.
!> \param basis Result for this L, including all radial functions and errors.
! **************************************************************************************************
   SUBROUTINE build_radial_basis(shells, l, tolerance, basis)
      TYPE(radial_product_type), INTENT(IN)              :: shells(:)
      INTEGER, INTENT(IN)                                :: l
      REAL(dp), INTENT(IN)                               :: tolerance
      TYPE(radial_basis_type), INTENT(OUT)               :: basis

      CHARACTER(LEN=256)                                 :: message
      INTEGER                                            :: attempt, i, ip, j, jp, k, m, n, &
                                                            nproduct, nradial, p
      INTEGER, ALLOCATABLE                               :: selected(:)
      REAL(dp)                                           :: zmax, zmin
      REAL(dp), ALLOCATABLE                              :: coefficients(:, :), errors(:), &
                                                            gram(:, :), norms(:), zet(:)
      TYPE(radial_product_type), ALLOCATABLE             :: products(:)

      nproduct = 0
      DO i = 1, SIZE(shells)
         DO j = i, SIZE(shells)
            IF (allowed_product(shells(i)%power, shells(j)%power, l)) nproduct = nproduct + 1
         END DO
      END DO
      basis%nproduct = nproduct
      IF (nproduct == 0) RETURN
      ALLOCATE (products(nproduct), errors(nproduct), norms(nproduct), selected(nproduct))
      k = 0
      DO i = 1, SIZE(shells)
         DO j = i, SIZE(shells)
            IF (.NOT. allowed_product(shells(i)%power, shells(j)%power, l)) CYCLE
            k = k + 1
            n = SIZE(shells(i)%zet)
            m = SIZE(shells(j)%zet)
            products(k)%power = shells(i)%power + shells(j)%power
            ALLOCATE (products(k)%zet(n*m), products(k)%coeff(n*m))
            p = 0
            DO ip = 1, n
               DO jp = 1, m
                  p = p + 1
                  products(k)%zet(p) = shells(i)%zet(ip) + shells(j)%zet(jp)
                  products(k)%coeff(p) = shells(i)%coeff(ip)*shells(j)%coeff(jp)
               END DO
            END DO
            CALL normalize_radial(products(k))
            norms(k) = products(k)%norm2
         END DO
      END DO
      zmin = HUGE(1.0_dp)
      zmax = 0.0_dp
      DO i = 1, nproduct
         ! Lower-L channels include r**(2q). Bracket the product exponent on both sides, so its
         ! Gaussian representation need not approximate exponent derivatives at a grid boundary.
         zmin = MIN(zmin, MINVAL(products(i)%zet)/(1 + (products(i)%power - l)/2))
         zmax = MAX(zmax, MAXVAL(products(i)%zet)*(1 + (products(i)%power - l)/2))
      END DO
      DO attempt = 0, 4
         CALL gaussian_candidates(products, l, zmin, zmax, attempt, zet)
         CALL fit_radial_products(products, l, zet, coefficients, gram)
         DO i = 1, nproduct
            errors(i) = MAX(0.0_dp, 1.0_dp - gram(i, i))
         END DO
         IF (MAXVAL(errors) <= 0.5_dp*tolerance) EXIT
      END DO
      IF (MAXVAL(errors) > 0.5_dp*tolerance) THEN
         WRITE (message, '(A,I3,A,ES12.4,A)') 'Automatic LRI-HFX basis: L=', l, &
            ' Gaussian product residual=', MAXVAL(errors), '; supply an explicit RI_HFX basis'
         CPABORT(TRIM(message))
      END IF
      CALL select_radial_products(gram, norms, tolerance, selected, nradial, errors)
      basis%zet = zet
      ! Keep fitted products, not normalized Cholesky residuals. They span the selected space while
      ! avoiding large cancelling coefficients in nearly dependent radial directions.
      basis%coeff = coefficients(:, selected(1:nradial))
      basis%residual = MAXVAL(errors)
      basis%absolute_residual = MAXVAL(errors*norms)
   END SUBROUTINE build_radial_basis

! **************************************************************************************************
!> \brief Gaunt triangle and parity rules for two same-center spherical harmonics.
!> \param la First ORB angular momentum.
!> \param lb Second ORB angular momentum.
!> \param l Product angular momentum.
!> \return Whether this channel can occur.
! **************************************************************************************************
   PURE FUNCTION allowed_product(la, lb, l) RESULT(allowed)
      INTEGER, INTENT(IN)                                :: la, lb, l
      LOGICAL                                            :: allowed

      allowed = ABS(la - lb) <= l .AND. l <= la + lb .AND. MOD(la + lb - l, 2) == 0
   END FUNCTION allowed_product

! **************************************************************************************************
!> \brief Candidate exponents, refined only if their product representation is insufficient.
!> \param products Normalized radial products.
!> \param l Auxiliary angular momentum.
!> \param zmin Diffuse end of the Gaussian grid.
!> \param zmax Tight end of the Gaussian grid.
!> \param refinement Grid refinement iteration.
!> \param zet Sorted unique candidate exponents.
! **************************************************************************************************
   SUBROUTINE gaussian_candidates(products, l, zmin, zmax, refinement, zet)
      TYPE(radial_product_type), INTENT(IN)              :: products(:)
      INTEGER, INTENT(IN)                                :: l
      REAL(dp), INTENT(IN)                               :: zmin, zmax
      INTEGER, INTENT(IN)                                :: refinement
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: zet(:)

      INTEGER                                            :: i, j, n, ngrid, nunique
      INTEGER, ALLOCATABLE                               :: order(:)
      LOGICAL                                            :: exact
      REAL(dp), ALLOCATABLE                              :: candidates(:)

      exact = ALL(products%power == l)
      ngrid = 0
      IF (.NOT. exact) ngrid = MAX(2, 1 + CEILING(LOG(zmax/zmin)/LOG(1.6_dp)*1.5_dp**refinement))
      n = ngrid
      DO i = 1, SIZE(products)
         IF (exact .OR. SIZE(products(i)%zet) == 1) n = n + SIZE(products(i)%zet)
      END DO
      ALLOCATE (candidates(n), order(n))
      DO i = 1, ngrid
         candidates(i) = zmin*EXP(LOG(zmax/zmin)*REAL(i - 1, dp)/REAL(ngrid - 1, dp))
      END DO
      n = ngrid
      DO i = 1, SIZE(products)
         ! Include single-primitive exponent sums, including diffuse s*p and p*p products, and
         ! centers for the exponent derivatives needed to represent lower-L polynomial factors.
         ! If all products are already ordinary GTOs, no Gaussian fitting grid is needed at all.
         IF (exact .OR. SIZE(products(i)%zet) == 1) THEN
            j = SIZE(products(i)%zet)
            candidates(n + 1:n + j) = products(i)%zet
            n = n + j
         END IF
      END DO
      CALL sort(candidates, n, order)
      nunique = 1
      DO i = 2, n
         IF (candidates(i) - candidates(nunique) <= 16*EPSILON(1.0_dp)*candidates(i)) CYCLE
         nunique = nunique + 1
         candidates(nunique) = candidates(i)
      END DO
      zet = candidates(1:nunique)
   END SUBROUTINE gaussian_candidates

! **************************************************************************************************
!> \brief Fit radial products in normalized Gaussians using analytic ordinary overlaps.
!>        The overlap cutoff only removes numerical null modes; it is independent of EPS_EIGVAL.
!> \param products Normalized radial targets, including lower-L polynomial factors.
!> \param l Auxiliary angular momentum.
!> \param zet Candidate exponents.
!> \param coefficients Fitted product coefficients on normalized primitives.
!> \param gram Overlap Gram matrix of the fitted products.
! **************************************************************************************************
   SUBROUTINE fit_radial_products(products, l, zet, coefficients, gram)
      TYPE(radial_product_type), INTENT(IN)              :: products(:)
      INTEGER, INTENT(IN)                                :: l
      REAL(dp), INTENT(IN)                               :: zet(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: coefficients(:, :), gram(:, :)

      REAL(dp), PARAMETER                                :: overlap_cutoff = 1.0E-14_dp

      INTEGER                                            :: first, i, info, j, k, n, nproduct
      REAL(dp), ALLOCATABLE                              :: eigenvalues(:), norm(:), overlap(:, :), &
                                                            projected(:, :), rhs(:, :), work(:)

      n = SIZE(zet)
      nproduct = SIZE(products)
      ALLOCATE (overlap(n, n), norm(n), coefficients(n, nproduct), gram(nproduct, nproduct))
      norm = 1.0_dp/SQRT(radial_integral(2*zet, 2*l))
      DO i = 1, n
         overlap(:, i) = norm*norm(i)*radial_integral(zet + zet(i), 2*l)
      END DO
      IF (ALL(products%power == l)) THEN
         ! Exact ordinary-GTO products need neither grid fitting nor an ill-conditioned inverse.
         coefficients = 0.0_dp
         DO i = 1, nproduct
            DO j = 1, SIZE(products(i)%zet)
               k = MINLOC(ABS(zet - products(i)%zet(j)), DIM=1)
               coefficients(k, i) = coefficients(k, i) + products(i)%coeff(j)/norm(k)
            END DO
         END DO
         gram = MATMUL(TRANSPOSE(coefficients), MATMUL(overlap, coefficients))
      ELSE
         ALLOCATE (eigenvalues(n), work(MAX(1, 3*n - 1)), rhs(n, nproduct))
         DO i = 1, nproduct
            DO j = 1, n
               rhs(j, i) = norm(j)*SUM(products(i)%coeff* &
                                       radial_integral(zet(j) + products(i)%zet, l + products(i)%power))
            END DO
         END DO
         CALL dsyev('V', 'U', n, overlap, n, eigenvalues, work, SIZE(work), info)
         CPASSERT(info == 0)
         CPASSERT(eigenvalues(n) > 0.0_dp)
         first = COUNT(eigenvalues <= overlap_cutoff*eigenvalues(n)) + 1
         DO i = first, n
            overlap(:, i) = overlap(:, i)/SQRT(eigenvalues(i))
         END DO
         projected = MATMUL(TRANSPOSE(overlap(:, first:n)), rhs)
         coefficients = MATMUL(overlap(:, first:n), projected)
         gram = MATMUL(TRANSPOSE(projected), projected)
      END IF
      CPASSERT(ALL(ieee_is_finite(coefficients)))
      CPASSERT(ALL(ieee_is_finite(gram)))
   END SUBROUTINE fit_radial_products

! **************************************************************************************************
!> \brief Select products until every normalized target meets the squared radial error tolerance.
!> \param gram Gram matrix of the represented, normalized targets (diagonal at most one).
!> \param norms Squared norms of products of normalized ORB shells, used only to break error ties.
!> \param tolerance Maximum relative squared error, including the Gaussian fitting residual.
!> \param selected Indices of retained products in pivot order.
!> \param nselected Number of retained radial functions.
!> \param errors Final per-product squared errors.
! **************************************************************************************************
   SUBROUTINE select_radial_products(gram, norms, tolerance, selected, nselected, errors)
      REAL(dp), INTENT(IN)                               :: gram(:, :), norms(:), tolerance
      INTEGER, INTENT(OUT)                               :: selected(:), nselected
      REAL(dp), INTENT(OUT)                              :: errors(:)

      INTEGER                                            :: i, n, pivot
      REAL(dp)                                           :: maximum, norm2
      REAL(dp), ALLOCATABLE                              :: cholesky(:, :), column(:)

      n = SIZE(gram, 1)
      ALLOCATE (cholesky(n, n), column(n))
      cholesky = 0.0_dp
      errors = 1.0_dp
      nselected = 0
      DO WHILE (MAXVAL(errors) > tolerance .AND. nselected < n)
         maximum = MAXVAL(errors)
         ! Equal-error pivots (in particular the first) favor the more localized product, rather
         ! than an arbitrary ORB set order. All subsequent decisions still use relative errors.
         pivot = MAXLOC(norms, DIM=1, MASK=errors >= maximum - 32*EPSILON(1.0_dp))
         column = gram(:, pivot)
         DO i = 1, nselected
            column = column - cholesky(:, i)*cholesky(pivot, i)
         END DO
         norm2 = column(pivot)
         CPASSERT(norm2 > 0.0_dp)
         nselected = nselected + 1
         selected(nselected) = pivot
         cholesky(:, nselected) = column/SQRT(norm2)
         errors = MAX(0.0_dp, 1.0_dp - SUM(cholesky(:, 1:nselected)**2, DIM=2))
      END DO
      CPASSERT(MAXVAL(errors) <= tolerance)
   END SUBROUTINE select_radial_products
END MODULE hfx_lri_basis
