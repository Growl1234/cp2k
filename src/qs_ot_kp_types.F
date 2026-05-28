!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Type definitions for orbital transformation with k-point calculations
! **************************************************************************************************
MODULE qs_ot_kp_types
   USE cp_cfm_types,                    ONLY: cp_cfm_release,&
                                              cp_cfm_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ot_kp_types'

   PUBLIC :: qs_ot_kp_block_type, qs_ot_kp_env_type, qs_ot_kp_env_release

! **************************************************************************************************
!> \brief Per-local-k-point/per-spin work matrices for k-point OT
! **************************************************************************************************
   TYPE qs_ot_kp_block_type
      LOGICAL :: allocated = .FALSE.
      LOGICAL :: has_data = .FALSE.
      TYPE(cp_cfm_type) :: c0
      TYPE(cp_cfm_type) :: gx
      TYPE(cp_cfm_type) :: dx
      TYPE(cp_cfm_type) :: gx_old
      TYPE(cp_cfm_type) :: dx_old
   END TYPE qs_ot_kp_block_type

! **************************************************************************************************
!> \brief Persistent work state for fixed-occupation k-point OT
! **************************************************************************************************
   TYPE qs_ot_kp_env_type
      INTEGER :: kplocal = -1
      INTEGER :: nspin = -1
      LOGICAL :: has_history = .FALSE.
      LOGICAL :: needs_baseline = .FALSE.
      REAL(KIND=dp) :: gnorm = -1.0_dp
      REAL(KIND=dp) :: gnorm_old = -1.0_dp
      REAL(KIND=dp) :: step_size = 1.0E-4_dp
      REAL(KIND=dp) :: last_step_size = 0.0_dp
      TYPE(qs_ot_kp_block_type), DIMENSION(:, :), POINTER :: block => NULL()
   END TYPE qs_ot_kp_env_type

CONTAINS

! **************************************************************************************************
!> \brief Release a k-point OT work environment
!> \param ot_kp_env  k-point OT work environment
! **************************************************************************************************
   SUBROUTINE qs_ot_kp_env_release(ot_kp_env)

      TYPE(qs_ot_kp_env_type), POINTER                  :: ot_kp_env

      INTEGER                                           :: ikp, ispin

      IF (.NOT. ASSOCIATED(ot_kp_env)) RETURN

      IF (ASSOCIATED(ot_kp_env%block)) THEN
         DO ispin = 1, SIZE(ot_kp_env%block, 2)
            DO ikp = 1, SIZE(ot_kp_env%block, 1)
               IF (ot_kp_env%block(ikp, ispin)%allocated) THEN
                  CALL cp_cfm_release(ot_kp_env%block(ikp, ispin)%c0)
                  CALL cp_cfm_release(ot_kp_env%block(ikp, ispin)%gx)
                  CALL cp_cfm_release(ot_kp_env%block(ikp, ispin)%dx)
                  CALL cp_cfm_release(ot_kp_env%block(ikp, ispin)%gx_old)
                  CALL cp_cfm_release(ot_kp_env%block(ikp, ispin)%dx_old)
               END IF
            END DO
         END DO
         DEALLOCATE (ot_kp_env%block)
      END IF

      DEALLOCATE (ot_kp_env)
      NULLIFY (ot_kp_env)

   END SUBROUTINE qs_ot_kp_env_release

END MODULE qs_ot_kp_types
