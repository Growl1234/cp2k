!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Pair-local LRI backend for k-point HFX.
!> \author SY Wang (05.2026)
! **************************************************************************************************
MODULE hfx_lri_kp
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE gamma,                           ONLY: init_md_ftable
   USE hfx_lri_types,                   ONLY: hfx_lri_kp_release
   USE hfx_types,                       ONLY: hfx_ri_type
   USE input_constants,                 ONLY: do_potential_short
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE libint_2c_3c,                    ONLY: eri_2center,&
                                              eri_3center
   USE libint_wrapper,                  ONLY: cp_libint_cleanup_2eri,&
                                              cp_libint_cleanup_3eri,&
                                              cp_libint_init_2eri,&
                                              cp_libint_init_3eri,&
                                              cp_libint_set_contrdepth,&
                                              cp_libint_t
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: ncoset
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_neighbor_lists,               ONLY: setup_neighbor_list
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hfx_lri_update_ks_kp

CONTAINS

! **************************************************************************************************
!> \brief Updates the k-point KS matrices with the pair-local LRI HFX contribution.
!> \param qs_env ...
!> \param ri_data ...
!> \param ks_matrix ...
!> \param ehfx ...
!> \param rho_ao ...
!> \param geometry_did_change ...
!> \param nspins ...
!> \param hf_fraction ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_update_ks_kp(qs_env, ri_data, ks_matrix, ehfx, rho_ao, &
                                   geometry_did_change, nspins, hf_fraction)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      REAL(KIND=dp), INTENT(OUT)                         :: ehfx
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      LOGICAL, INTENT(IN)                                :: geometry_did_change
      INTEGER, INTENT(IN)                                :: nspins
      REAL(KIND=dp), INTENT(IN)                          :: hf_fraction

      ehfx = 0.0_dp

      IF (ri_data%hfx_pot%potential_type /= do_potential_short) THEN
         CPABORT("LRI-HFXk is currently restricted to short-range/HSE exchange")
      END IF

      ! Keep dummy arguments referenced in this scaffold.
      IF (ASSOCIATED(qs_env) .AND. ASSOCIATED(ks_matrix) .AND. ASSOCIATED(rho_ao) .AND. &
          geometry_did_change .AND. nspins > 0 .AND. hf_fraction >= 0.0_dp) THEN
         CONTINUE
      END IF

      IF (geometry_did_change .OR. .NOT. ASSOCIATED(ri_data%lri_kp)) THEN
         CALL hfx_lri_pre_scf_kp(qs_env, ri_data)
      END IF

      CPABORT("LRI-HFXk pair-local contraction is not implemented yet")

   END SUBROUTINE hfx_lri_update_ks_kp

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_pre_scf_kp(qs_env, ri_data)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_lri_pre_scf_kp'

      INTEGER :: dbg_iatom, dbg_ikind, dbg_jatom, dbg_ji_ncart, dbg_jj_ncart, dbg_jkind, &
         dbg_m_ncart, dbg_nao_i, dbg_nao_j, dbg_naux_i, dbg_naux_j, dbg_naux_pair, handle, iatom, &
         ibasis, ikind, ilist, jatom, jkind, jneighbor, mepos, nao_i, nao_j, naux_i, naux_j, &
         naux_pair, nkind, nlist, nneighbor, npairs_local, npairs_max, npairs_sum
      INTEGER(int_8) :: c_nreal, c_nreal_local, c_nreal_max, c_nreal_sum, sum_nao_pair_local, &
         sum_nao_pair_max, sum_nao_pair_sum, sum_naux_pair_local, sum_naux_pair_max, &
         sum_naux_pair_sum
      INTEGER, DIMENSION(3)                              :: cell_j, dbg_cell
      LOGICAL                                            :: dbg_integrals_done, dbg_pair_found
      REAL(dp)                                           :: c_mem_mib_local, c_mem_mib_max, &
                                                            c_mem_mib_sum, dbg_ji_norm, &
                                                            dbg_jj_norm, dbg_m_norm, dbg_r2
      REAL(dp), DIMENSION(3)                             :: dbg_rab, rab
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_AO, basis_set_RI
      TYPE(gto_basis_set_type), POINTER                  :: obas_i, obas_j, rbas_i, rbas_j
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: soo_list
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (qs_kind_set, soo_list, nl_iterator, para_env, dft_control)
      NULLIFY (obas_i, obas_j, rbas_i, rbas_j)
      NULLIFY (basis_set_AO, basis_set_RI)

      dbg_pair_found = .FALSE.
      dbg_iatom = 0
      dbg_jatom = 0
      dbg_ikind = 0
      dbg_jkind = 0
      dbg_nao_i = 0
      dbg_nao_j = 0
      dbg_naux_i = 0
      dbg_naux_j = 0
      dbg_naux_pair = 0
      dbg_cell(:) = 0
      dbg_r2 = HUGE(1.0_dp)
      dbg_rab(:) = 0.0_dp
      dbg_integrals_done = .FALSE.
      dbg_m_ncart = 0
      dbg_ji_ncart = 0
      dbg_jj_ncart = 0
      dbg_m_norm = 0.0_dp
      dbg_ji_norm = 0.0_dp
      dbg_jj_norm = 0.0_dp

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, nkind=nkind, para_env=para_env, &
                      dft_control=dft_control)

      CALL hfx_lri_kp_release(ri_data%lri_kp)
      ALLOCATE (ri_data%lri_kp)

      ALLOCATE (basis_set_AO(nkind), basis_set_RI(nkind))
      CALL basis_set_list_setup(basis_set_AO, ri_data%orb_basis_type, qs_kind_set)
      CALL basis_set_list_setup(basis_set_RI, ri_data%ri_basis_type, qs_kind_set)

      ! Use RI-HFX screening radii while building the LRI-HFXk pair list.
      ! They are restored to the standard QS radii before leaving this routine.
      DO ibasis = 1, SIZE(basis_set_AO)
         obas_i => basis_set_AO(ibasis)%gto_basis_set
         rbas_i => basis_set_RI(ibasis)%gto_basis_set
         IF (ASSOCIATED(obas_i)) THEN
            CALL init_interaction_radii_orb_basis(obas_i, ri_data%eps_pgf_orb)
         END IF
         IF (ASSOCIATED(rbas_i)) THEN
            CALL init_interaction_radii_orb_basis(rbas_i, ri_data%eps_pgf_orb)
         END IF
      END DO

      ! First prototype: AO-AO neighbor list only.
      ! Later this cutoff should be tightened for HSE short-range LRI.
      CALL setup_neighbor_list(soo_list, basis_set_AO, qs_env=qs_env, mic=.FALSE., symmetric=.TRUE.)

      CALL neighbor_list_iterator_create(nl_iterator, soo_list)

      mepos = 0
      DO WHILE (neighbor_list_iterate(nl_iterator, mepos) == 0)

         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, &
                                nlist=nlist, ilist=ilist, nnode=nneighbor, inode=jneighbor, &
                                iatom=iatom, jatom=jatom, r=rab, cell=cell_j)

         obas_i => basis_set_AO(ikind)%gto_basis_set
         obas_j => basis_set_AO(jkind)%gto_basis_set
         rbas_i => basis_set_RI(ikind)%gto_basis_set
         rbas_j => basis_set_RI(jkind)%gto_basis_set

         IF (.NOT. ASSOCIATED(obas_i)) CYCLE
         IF (.NOT. ASSOCIATED(obas_j)) CYCLE
         IF (.NOT. ASSOCIATED(rbas_i)) CYCLE
         IF (.NOT. ASSOCIATED(rbas_j)) CYCLE

         nao_i = obas_i%nsgf
         nao_j = obas_j%nsgf
         naux_i = rbas_i%nsgf
         naux_j = rbas_j%nsgf

         ! Pair-local / two-center auxiliary space:
         ! Aux(ij) = RI(i) union RI(j).
         IF (iatom == jatom .AND. SUM(rab*rab) < 1.0E-20_dp) THEN
            naux_pair = naux_i
         ELSE
            naux_pair = naux_i + naux_j
         END IF

         IF ((.NOT. dbg_pair_found) .OR. SUM(rab*rab) < dbg_r2) THEN
            dbg_pair_found = .TRUE.
            dbg_r2 = SUM(rab*rab)
            dbg_iatom = iatom
            dbg_jatom = jatom
            dbg_ikind = ikind
            dbg_jkind = jkind
            dbg_nao_i = nao_i
            dbg_nao_j = nao_j
            dbg_naux_i = naux_i
            dbg_naux_j = naux_j
            dbg_naux_pair = naux_pair
            dbg_cell(:) = cell_j(:)
            dbg_rab(:) = rab(:)
         END IF

         c_nreal = INT(nao_i, int_8)*INT(nao_j, int_8)*INT(naux_pair, int_8)

         ri_data%lri_kp%npairs = ri_data%lri_kp%npairs + 1
         ri_data%lri_kp%sum_nao_pair = ri_data%lri_kp%sum_nao_pair + INT(nao_i, int_8)*INT(nao_j, int_8)
         ri_data%lri_kp%sum_naux_pair = ri_data%lri_kp%sum_naux_pair + INT(naux_pair, int_8)
         ri_data%lri_kp%c_nreal_est = ri_data%lri_kp%c_nreal_est + c_nreal

      END DO

      CALL neighbor_list_iterator_release(nl_iterator)
      CALL release_neighbor_list_sets(soo_list)

      IF (dbg_pair_found) THEN
         CALL hfx_lri_debug_pair_integrals(basis_set_AO, basis_set_RI, dbg_ikind, dbg_jkind, &
                                           dbg_rab, ri_data, dbg_m_ncart, dbg_ji_ncart, &
                                           dbg_jj_ncart, dbg_m_norm, dbg_ji_norm, dbg_jj_norm)
         dbg_integrals_done = .TRUE.
      END IF

      ! Restore the standard QS screening radii after the LRI-HFXk setup.
      IF (ASSOCIATED(dft_control)) THEN
         DO ibasis = 1, SIZE(basis_set_AO)
            obas_i => basis_set_AO(ibasis)%gto_basis_set
            rbas_i => basis_set_RI(ibasis)%gto_basis_set
            IF (ASSOCIATED(obas_i)) THEN
               CALL init_interaction_radii_orb_basis(obas_i, dft_control%qs_control%eps_pgf_orb)
            END IF
            IF (ASSOCIATED(rbas_i)) THEN
               CALL init_interaction_radii_orb_basis(rbas_i, dft_control%qs_control%eps_pgf_orb)
            END IF
         END DO
      END IF

      IF (ASSOCIATED(basis_set_AO)) DEALLOCATE (basis_set_AO)
      IF (ASSOCIATED(basis_set_RI)) DEALLOCATE (basis_set_RI)

      npairs_local = ri_data%lri_kp%npairs
      sum_nao_pair_local = ri_data%lri_kp%sum_nao_pair
      sum_naux_pair_local = ri_data%lri_kp%sum_naux_pair
      c_nreal_local = ri_data%lri_kp%c_nreal_est
      c_mem_mib_local = REAL(c_nreal_local, KIND=dp)*8.0_dp/1024.0_dp**2

      npairs_sum = npairs_local
      npairs_max = npairs_local
      sum_nao_pair_sum = sum_nao_pair_local
      sum_nao_pair_max = sum_nao_pair_local
      sum_naux_pair_sum = sum_naux_pair_local
      sum_naux_pair_max = sum_naux_pair_local
      c_nreal_sum = c_nreal_local
      c_nreal_max = c_nreal_local

      IF (ASSOCIATED(para_env)) THEN
         CALL para_env%sum(npairs_sum)
         CALL para_env%max(npairs_max)
         CALL para_env%sum(sum_nao_pair_sum)
         CALL para_env%max(sum_nao_pair_max)
         CALL para_env%sum(sum_naux_pair_sum)
         CALL para_env%max(sum_naux_pair_max)
         CALL para_env%sum(c_nreal_sum)
         CALL para_env%max(c_nreal_max)
      END IF

      c_mem_mib_sum = REAL(c_nreal_sum, KIND=dp)*8.0_dp/1024.0_dp**2
      c_mem_mib_max = REAL(c_nreal_max, KIND=dp)*8.0_dp/1024.0_dp**2

      ri_data%lri_kp%c_mem_mib_est = c_mem_mib_local
      ri_data%lri_kp%initialized = .TRUE.

      IF (ri_data%unit_nr > 0) THEN
         WRITE (ri_data%unit_nr, "(/,T2,A)") "LRI-HFXk pair-local setup"
         WRITE (ri_data%unit_nr, "(T4,A,T28,A)") "orbital basis", TRIM(ri_data%orb_basis_type)
         WRITE (ri_data%unit_nr, "(T4,A,T28,A)") "RI basis", TRIM(ri_data%ri_basis_type)
         WRITE (ri_data%unit_nr, "(/,T4,A,T30,A,T43,A,T56,A)") &
            "quantity", "local", "global", "max/rank"
         WRITE (ri_data%unit_nr, "(T4,A,T28,I10,T41,I10,T54,I10)") &
            "AO pairs", npairs_local, npairs_sum, npairs_max
         WRITE (ri_data%unit_nr, "(T4,A,T28,I10,T41,I10,T54,I10)") &
            "AO-pair size", sum_nao_pair_local, sum_nao_pair_sum, sum_nao_pair_max
         WRITE (ri_data%unit_nr, "(T4,A,T28,I10,T41,I10,T54,I10)") &
            "pair aux size", sum_naux_pair_local, sum_naux_pair_sum, sum_naux_pair_max
         WRITE (ri_data%unit_nr, "(T4,A,T28,F10.2,T41,F10.2,T54,F10.2)") &
            "C storage/MiB", c_mem_mib_local, c_mem_mib_sum, c_mem_mib_max
         IF (npairs_sum > 0) THEN
            WRITE (ri_data%unit_nr, "(T4,A,T28,F10.2,T41,F10.2,T54,F10.2)") &
               "avg AO-pair", &
               REAL(sum_nao_pair_local, KIND=dp)/MAX(1.0_dp, REAL(npairs_local, KIND=dp)), &
               REAL(sum_nao_pair_sum, KIND=dp)/REAL(npairs_sum, KIND=dp), &
               REAL(sum_nao_pair_max, KIND=dp)/MAX(1.0_dp, REAL(npairs_max, KIND=dp))
            WRITE (ri_data%unit_nr, "(T4,A,T28,F10.2,T41,F10.2,T54,F10.2)") &
               "avg pair aux", &
               REAL(sum_naux_pair_local, KIND=dp)/MAX(1.0_dp, REAL(npairs_local, KIND=dp)), &
               REAL(sum_naux_pair_sum, KIND=dp)/REAL(npairs_sum, KIND=dp), &
               REAL(sum_naux_pair_max, KIND=dp)/MAX(1.0_dp, REAL(npairs_max, KIND=dp))
         END IF
         IF (dbg_pair_found) THEN
            WRITE (ri_data%unit_nr, "(/,T4,A)") "debug local AO pair"
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
               "atoms", dbg_iatom, dbg_jatom
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
               "kinds", dbg_ikind, dbg_jkind
            WRITE (ri_data%unit_nr, "(T6,A,T28,3I6)") &
               "cell", dbg_cell
            WRITE (ri_data%unit_nr, "(T6,A,T28,F10.4)") &
               "distance", SQRT(dbg_r2)
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10)") &
               "AO functions", dbg_nao_i, dbg_nao_j
            WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,I10,T54,I10)") &
               "RI functions", dbg_naux_i, dbg_naux_j, dbg_naux_pair
            IF (dbg_integrals_done) THEN
               WRITE (ri_data%unit_nr, "(/,T4,A)") "debug pair short-range integrals"
               WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,ES10.3)") &
                  "metric ncart/norm", dbg_m_ncart, dbg_m_norm
               WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,ES10.3)") &
                  "3c aux-i ncart/norm", dbg_ji_ncart, dbg_ji_norm
               WRITE (ri_data%unit_nr, "(T6,A,T28,I10,T41,ES10.3)") &
                  "3c aux-j ncart/norm", dbg_jj_ncart, dbg_jj_norm
            END IF
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE hfx_lri_pre_scf_kp

! **************************************************************************************************
!> \brief ...
!> \param basis_set_AO ...
!> \param basis_set_RI ...
!> \param ikind ...
!> \param jkind ...
!> \param rab ...
!> \param ri_data ...
!> \param metric_ncart ...
!> \param threec_i_ncart ...
!> \param threec_j_ncart ...
!> \param metric_norm ...
!> \param threec_i_norm ...
!> \param threec_j_norm ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_debug_pair_integrals(basis_set_AO, basis_set_RI, ikind, jkind, rab, ri_data, &
                                           metric_ncart, threec_i_ncart, threec_j_ncart, &
                                           metric_norm, threec_i_norm, threec_j_norm)

      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set_AO, basis_set_RI
      INTEGER, INTENT(IN)                                :: ikind, jkind
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rab
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      INTEGER, INTENT(OUT)                               :: metric_ncart, threec_i_ncart, &
                                                            threec_j_ncart
      REAL(dp), INTENT(OUT)                              :: metric_norm, threec_i_norm, threec_j_norm

      INTEGER                                            :: maxl_ai, maxl_aj, maxl_ri, maxl_rj, mmax
      LOGICAL                                            :: onsite_aux
      REAL(dp)                                           :: metric_norm2, threec_i_norm2, &
                                                            threec_j_norm2
      REAL(dp), DIMENSION(3)                             :: r0
      TYPE(cp_libint_t)                                  :: lib2, lib3
      TYPE(gto_basis_set_type), POINTER                  :: obas_i, obas_j, rbas_i, rbas_j

      NULLIFY (obas_i, obas_j, rbas_i, rbas_j)

      obas_i => basis_set_AO(ikind)%gto_basis_set
      obas_j => basis_set_AO(jkind)%gto_basis_set
      rbas_i => basis_set_RI(ikind)%gto_basis_set
      rbas_j => basis_set_RI(jkind)%gto_basis_set

      metric_ncart = 0
      threec_i_ncart = 0
      threec_j_ncart = 0
      metric_norm = 0.0_dp
      threec_i_norm = 0.0_dp
      threec_j_norm = 0.0_dp

      IF (.NOT. ASSOCIATED(obas_i)) RETURN
      IF (.NOT. ASSOCIATED(obas_j)) RETURN
      IF (.NOT. ASSOCIATED(rbas_i)) RETURN
      IF (.NOT. ASSOCIATED(rbas_j)) RETURN

      CALL get_gto_basis_set(obas_i, maxl=maxl_ai)
      CALL get_gto_basis_set(obas_j, maxl=maxl_aj)
      CALL get_gto_basis_set(rbas_i, maxl=maxl_ri)
      CALL get_gto_basis_set(rbas_j, maxl=maxl_rj)

      mmax = MAX(maxl_ri + maxl_rj, maxl_ai + maxl_aj + MAX(maxl_ri, maxl_rj))
      CALL init_md_ftable(nmax=mmax)

      CALL cp_libint_init_2eri(lib2, MAX(maxl_ri, maxl_rj))
      CALL cp_libint_set_contrdepth(lib2, 1)
      CALL cp_libint_init_3eri(lib3, MAX(maxl_ai, maxl_aj, maxl_ri, maxl_rj))
      CALL cp_libint_set_contrdepth(lib3, 1)

      r0(:) = 0.0_dp
      metric_norm2 = 0.0_dp
      threec_i_norm2 = 0.0_dp
      threec_j_norm2 = 0.0_dp
      onsite_aux = (NORM2(rab) <= 1.0E-12_dp .AND. ikind == jkind)

      CALL hfx_lri_accumulate_2c_norm(rbas_i, rbas_i, r0, r0, lib2, ri_data, metric_ncart, &
                                      metric_norm2)
      IF (.NOT. onsite_aux) THEN
         CALL hfx_lri_accumulate_2c_norm(rbas_i, rbas_j, r0, rab, lib2, ri_data, metric_ncart, &
                                         metric_norm2, prefactor=2.0_dp)
         CALL hfx_lri_accumulate_2c_norm(rbas_j, rbas_j, rab, rab, lib2, ri_data, metric_ncart, &
                                         metric_norm2)
      END IF

      CALL hfx_lri_accumulate_3c_norm(obas_i, obas_j, rbas_i, r0, rab, r0, lib3, ri_data, &
                                      threec_i_ncart, threec_i_norm2)
      IF (.NOT. onsite_aux) THEN
         CALL hfx_lri_accumulate_3c_norm(obas_i, obas_j, rbas_j, r0, rab, rab, lib3, ri_data, &
                                         threec_j_ncart, threec_j_norm2)
      END IF

      CALL cp_libint_cleanup_2eri(lib2)
      CALL cp_libint_cleanup_3eri(lib3)

      metric_norm = SQRT(metric_norm2)
      threec_i_norm = SQRT(threec_i_norm2)
      threec_j_norm = SQRT(threec_j_norm2)

   END SUBROUTINE hfx_lri_debug_pair_integrals

! **************************************************************************************************
!> \brief ...
!> \param basis_a ...
!> \param basis_b ...
!> \param ra ...
!> \param rb ...
!> \param lib ...
!> \param ri_data ...
!> \param ncart_total ...
!> \param norm2_total ...
!> \param prefactor ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_accumulate_2c_norm(basis_a, basis_b, ra, rb, lib, ri_data, ncart_total, &
                                         norm2_total, prefactor)

      TYPE(gto_basis_set_type), POINTER                  :: basis_a, basis_b
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rb
      TYPE(cp_libint_t), INTENT(INOUT)                   :: lib
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      INTEGER, INTENT(INOUT)                             :: ncart_total
      REAL(dp), INTENT(INOUT)                            :: norm2_total
      REAL(dp), INTENT(IN), OPTIONAL                     :: prefactor

      INTEGER                                            :: iset, jset, ncoa, ncob, nseta, nsetb
      INTEGER, DIMENSION(:), POINTER                     :: lmax_a, lmax_b, lmin_a, lmin_b, npgf_a, &
                                                            npgf_b
      REAL(dp)                                           :: dab, pfac
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: sij
      REAL(dp), DIMENSION(:, :), POINTER                 :: rpgf_a, rpgf_b, zeta_a, zeta_b

      CALL get_gto_basis_set(basis_a, lmax=lmax_a, lmin=lmin_a, npgf=npgf_a, nset=nseta, &
                             pgf_radius=rpgf_a, zet=zeta_a)
      CALL get_gto_basis_set(basis_b, lmax=lmax_b, lmin=lmin_b, npgf=npgf_b, nset=nsetb, &
                             pgf_radius=rpgf_b, zet=zeta_b)

      pfac = 1.0_dp
      IF (PRESENT(prefactor)) pfac = prefactor
      dab = NORM2(rb - ra)

      DO iset = 1, nseta
         ncoa = npgf_a(iset)*ncoset(lmax_a(iset))
         DO jset = 1, nsetb
            ncob = npgf_b(jset)*ncoset(lmax_b(jset))
            IF (ncoa*ncob == 0) CYCLE
            ALLOCATE (sij(ncoa, ncob), SOURCE=0.0_dp)
            CALL eri_2center(sij, lmin_a(iset), lmax_a(iset), npgf_a(iset), zeta_a(:, iset), &
                             rpgf_a(:, iset), ra, lmin_b(jset), lmax_b(jset), npgf_b(jset), &
                             zeta_b(:, jset), rpgf_b(:, jset), rb, dab, lib, ri_data%hfx_pot)
            ncart_total = ncart_total + SIZE(sij)
            norm2_total = norm2_total + pfac*SUM(sij*sij)
            DEALLOCATE (sij)
         END DO
      END DO

   END SUBROUTINE hfx_lri_accumulate_2c_norm

! **************************************************************************************************
!> \brief ...
!> \param basis_a ...
!> \param basis_b ...
!> \param basis_c ...
!> \param ra ...
!> \param rb ...
!> \param rc ...
!> \param lib ...
!> \param ri_data ...
!> \param ncart_total ...
!> \param norm2_total ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_accumulate_3c_norm(basis_a, basis_b, basis_c, ra, rb, rc, lib, ri_data, &
                                         ncart_total, norm2_total)

      TYPE(gto_basis_set_type), POINTER                  :: basis_a, basis_b, basis_c
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rb, rc
      TYPE(cp_libint_t), INTENT(INOUT)                   :: lib
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data
      INTEGER, INTENT(INOUT)                             :: ncart_total
      REAL(dp), INTENT(INOUT)                            :: norm2_total

      INTEGER                                            :: iset, jset, kset, ncoa, ncob, ncoc, &
                                                            nseta, nsetb, nsetc
      INTEGER, DIMENSION(:), POINTER                     :: lmax_a, lmax_b, lmax_c, lmin_a, lmin_b, &
                                                            lmin_c, npgf_a, npgf_b, npgf_c
      REAL(dp)                                           :: dab, dac, dbc
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: sijk
      REAL(dp), DIMENSION(:, :), POINTER                 :: rpgf_a, rpgf_b, rpgf_c, zeta_a, zeta_b, &
                                                            zeta_c

      CALL get_gto_basis_set(basis_a, lmax=lmax_a, lmin=lmin_a, npgf=npgf_a, nset=nseta, &
                             pgf_radius=rpgf_a, zet=zeta_a)
      CALL get_gto_basis_set(basis_b, lmax=lmax_b, lmin=lmin_b, npgf=npgf_b, nset=nsetb, &
                             pgf_radius=rpgf_b, zet=zeta_b)
      CALL get_gto_basis_set(basis_c, lmax=lmax_c, lmin=lmin_c, npgf=npgf_c, nset=nsetc, &
                             pgf_radius=rpgf_c, zet=zeta_c)

      dab = NORM2(rb - ra)
      dac = NORM2(rc - ra)
      dbc = NORM2(rc - rb)

      DO iset = 1, nseta
         ncoa = npgf_a(iset)*ncoset(lmax_a(iset))
         DO jset = 1, nsetb
            ncob = npgf_b(jset)*ncoset(lmax_b(jset))
            DO kset = 1, nsetc
               ncoc = npgf_c(kset)*ncoset(lmax_c(kset))
               IF (ncoa*ncob*ncoc == 0) CYCLE
               ALLOCATE (sijk(ncoa, ncob, ncoc), SOURCE=0.0_dp)
               CALL eri_3center(sijk, lmin_a(iset), lmax_a(iset), npgf_a(iset), zeta_a(:, iset), &
                                rpgf_a(:, iset), ra, lmin_b(jset), lmax_b(jset), npgf_b(jset), &
                                zeta_b(:, jset), rpgf_b(:, jset), rb, lmin_c(kset), lmax_c(kset), &
                                npgf_c(kset), zeta_c(:, kset), rpgf_c(:, kset), rc, dab, dac, dbc, &
                                lib, ri_data%hfx_pot)
               ncart_total = ncart_total + SIZE(sijk)
               norm2_total = norm2_total + SUM(sijk*sijk)
               DEALLOCATE (sijk)
            END DO
         END DO
      END DO

   END SUBROUTINE hfx_lri_accumulate_3c_norm

END MODULE hfx_lri_kp
