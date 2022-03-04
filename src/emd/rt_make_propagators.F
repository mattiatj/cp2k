!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for calculating a complex matrix exponential.
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************

MODULE rt_make_propagators

   USE cp_control_types,                ONLY: rtp_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat
   USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_scale,&
                                              dbcsr_type
   USE input_constants,                 ONLY: do_etrs,&
                                              do_pade,&
                                              do_taylor
   USE kinds,                           ONLY: dp
   USE ls_matrix_exp,                   ONLY: bch_expansion_complex_propagator,&
                                              bch_expansion_imaginary_propagator,&
                                              cp_complex_dbcsr_gemm_3,&
                                              taylor_full_complex_dbcsr,&
                                              taylor_only_imaginary_dbcsr
   USE matrix_exp,                      ONLY: arnoldi,&
                                              exp_pade_full_complex,&
                                              exp_pade_only_imaginary,&
                                              taylor_full_complex,&
                                              taylor_only_imaginary
   USE rt_propagation_types,            ONLY: get_rtp,&
                                              rt_prop_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_make_propagators'

   PUBLIC :: propagate_exp, &
             propagate_arnoldi, &
             compute_exponential, &
             compute_exponential_sparse, &
             propagate_exp_density, &
             propagate_bch

CONTAINS
! **************************************************************************************************
!> \brief performs propagations if explicit matrix exponentials are used
!>        ETRS:  exp(i*H(t+dt)*dt/2)*exp(i*H(t)*dt/2)*MOS
!>        EM:    exp[-idt/2H(t+dt/2)*MOS
!> \param rtp ...
!> \param rtp_control ...
!> \param ortho ...
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************

   SUBROUTINE propagate_exp(rtp, rtp_control, ortho)

      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), POINTER                    :: rtp_control
      LOGICAL, INTENT(IN)                                :: ortho

      CHARACTER(len=*), PARAMETER                        :: routineN = 'propagate_exp'
      REAL(KIND=dp), PARAMETER                           :: one = 1.0_dp, zero = 0.0_dp

      INTEGER                                            :: handle, i, im, nao, ndep, nmo, re
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: exp_H_fm, exp_H_fm_red, mos_new, mos_next, &
         mos_old, propagator_fm, propagator_fm_red, propagator_matrix_fm, propagator_matrix_fm_red
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: exp_H_new, exp_H_old, propagator_matrix

      CALL timeset(routineN, handle)

      NULLIFY (mos_new, mos_next, mos_old, exp_H_new, exp_H_old, propagator_matrix)

      IF (.NOT. ortho) THEN
         CALL get_rtp(rtp=rtp, propagator_matrix=propagator_matrix, mos_old=mos_old, mos_new=mos_new, &
                      mos_next=mos_next, exp_H_new=exp_H_new, exp_H_old=exp_H_old)
      ELSE
         CALL get_rtp(rtp=rtp, propagator_matrix_canon=propagator_matrix, mos_ortho_old=mos_old, &
                      mos_ortho_new=mos_new, mos_ortho_next=mos_next, exp_H_canon_new=exp_H_new, exp_H_canon_old=exp_H_old)
      END IF

      NULLIFY (propagator_fm)
      ALLOCATE (propagator_fm(SIZE(exp_H_new)))
      NULLIFY (propagator_matrix_fm)
      ALLOCATE (propagator_matrix_fm(SIZE(propagator_matrix)))
      DO i = 1, SIZE(exp_H_new)
         CALL cp_fm_create(propagator_fm(i)%matrix, &
                           matrix_struct=rtp%ao_ao_fmstruct, &
                           name="prop_fm")
         CALL copy_dbcsr_to_fm(exp_H_new(i)%matrix, propagator_fm(i)%matrix)
         CALL cp_fm_create(propagator_matrix_fm(i)%matrix, &
                           matrix_struct=rtp%ao_ao_fmstruct, &
                           name="prop_mat_fm")
         CALL copy_dbcsr_to_fm(propagator_matrix(i)%matrix, propagator_matrix_fm(i)%matrix)
      END DO

      ! need to reduce dimension of propagator_matrix_fm if linear dependencies in S were removed
      NULLIFY (propagator_fm_red, propagator_matrix_fm_red)
      IF (ortho) THEN
         ALLOCATE (propagator_fm_red(SIZE(propagator_fm)))
         ALLOCATE (propagator_matrix_fm_red(SIZE(propagator_matrix_fm)))
         ndep = rtp_control%nlindep
         CALL cp_fm_get_info(propagator_fm(1)%matrix, nrow_global=nao)
         DO i = 1, SIZE(propagator_fm)
            CALL cp_fm_create(propagator_fm_red(i)%matrix, &
                              matrix_struct=rtp%canon_ortho%oo_oo_fmstruct, &
                              name="prop fm red")
            CALL cp_fm_to_fm_submat(propagator_fm(i)%matrix, propagator_fm_red(i)%matrix, &
                                    nao - ndep, nao - ndep, 1, 1, 1, 1)
            CALL cp_fm_create(propagator_matrix_fm_red(i)%matrix, &
                              matrix_struct=rtp%canon_ortho%oo_oo_fmstruct, &
                              name="prop fm red")
            CALL cp_fm_to_fm_submat(propagator_matrix_fm(i)%matrix, propagator_matrix_fm_red(i)%matrix, &
                                    nao - ndep, nao - ndep, 1, 1, 1, 1)
         END DO
      END IF

      ! Only compute exponential if a new propagator matrix is available
      IF (.NOT. ortho) THEN
         CALL compute_exponential(propagator_fm, propagator_matrix_fm, rtp_control, rtp)
      ELSE
         CALL compute_exponential(propagator_fm_red, propagator_matrix_fm_red, rtp_control, rtp)
         ! copy back to full dimensinal ao matrices in order to be compatible with dbcsr format
         ndep = rtp_control%nlindep
         CALL cp_fm_get_info(propagator_fm(1)%matrix, nrow_global=nao)
         DO i = 1, SIZE(propagator_fm)
            CALL cp_fm_set_all(propagator_fm(i)%matrix, 0.0_dp)
            CALL cp_fm_set_all(propagator_matrix_fm(i)%matrix, 0.0_dp)
            CALL cp_fm_to_fm_submat(propagator_fm_red(i)%matrix, propagator_fm(i)%matrix, &
                                    nao - ndep, nao - ndep, 1, 1, 1, 1)
            CALL cp_fm_to_fm_submat(propagator_matrix_fm_red(i)%matrix, propagator_matrix_fm(i)%matrix, &
                                    nao - ndep, nao - ndep, 1, 1, 1, 1)
         END DO
      END IF

      DO i = 1, SIZE(exp_H_new)
         CALL copy_fm_to_dbcsr(propagator_fm(i)%matrix, exp_H_new(i)%matrix)
         CALL copy_fm_to_dbcsr(propagator_matrix_fm(i)%matrix, propagator_matrix(i)%matrix)
      END DO

      CALL cp_fm_vect_dealloc(propagator_fm)
      CALL cp_fm_vect_dealloc(propagator_matrix_fm)

      IF (.NOT. ortho) THEN
         DO i = 1, SIZE(mos_new)/2
            re = 2*i - 1
            im = 2*i

            CALL cp_fm_get_info(mos_new(re)%matrix, ncol_global=nmo)
            !Save some work by computing the first half of the propagation only once in case of ETRS
            !For EM this matrix has to be the initial matrix, thus a copy is enough
            IF (rtp%iter == 1) THEN
               IF (rtp_control%propagator == do_etrs) THEN
                  CALL cp_dbcsr_sm_fm_multiply(exp_H_old(re)%matrix, mos_old(re)%matrix, &
                                               mos_next(re)%matrix, nmo, alpha=one, beta=zero)
                  CALL cp_dbcsr_sm_fm_multiply(exp_H_old(im)%matrix, mos_old(im)%matrix, &
                                               mos_next(re)%matrix, nmo, alpha=-one, beta=one)
                  CALL cp_dbcsr_sm_fm_multiply(exp_H_old(re)%matrix, mos_old(im)%matrix, &
                                               mos_next(im)%matrix, nmo, alpha=one, beta=zero)
                  CALL cp_dbcsr_sm_fm_multiply(exp_H_old(im)%matrix, mos_old(re)%matrix, &
                                               mos_next(im)%matrix, nmo, alpha=one, beta=one)
               ELSE
                  CALL cp_fm_to_fm(mos_old(re)%matrix, mos_next(re)%matrix)
                  CALL cp_fm_to_fm(mos_old(im)%matrix, mos_next(im)%matrix)
               END IF
            END IF
            CALL cp_dbcsr_sm_fm_multiply(exp_H_new(re)%matrix, mos_next(re)%matrix, &
                                         mos_new(re)%matrix, nmo, alpha=one, beta=zero)
            CALL cp_dbcsr_sm_fm_multiply(exp_H_new(im)%matrix, mos_next(im)%matrix, &
                                         mos_new(re)%matrix, nmo, alpha=-one, beta=one)
            CALL cp_dbcsr_sm_fm_multiply(exp_H_new(re)%matrix, mos_next(im)%matrix, &
                                         mos_new(im)%matrix, nmo, alpha=one, beta=zero)
            CALL cp_dbcsr_sm_fm_multiply(exp_H_new(im)%matrix, mos_next(re)%matrix, &
                                         mos_new(im)%matrix, nmo, alpha=one, beta=one)
         END DO
      ELSE
         ! get dimensions
         CALL dbcsr_get_info(exp_H_new(1)%matrix, nfullrows_total=nao)
         ndep = rtp_control%nlindep

         IF (rtp%iter == 1) THEN
            IF (rtp_control%propagator == do_etrs) THEN
               NULLIFY (exp_H_fm)
               ALLOCATE (exp_H_fm(SIZE(exp_H_old)))
               DO i = 1, SIZE(exp_H_old)
                  CALL cp_fm_create(exp_H_fm(i)%matrix, matrix_struct=rtp%ao_ao_fmstruct, name="exp_H_fm")
                  CALL copy_dbcsr_to_fm(exp_H_old(i)%matrix, exp_H_fm(i)%matrix)
               END DO
               NULLIFY (exp_H_fm_red)
               ALLOCATE (exp_H_fm_red(SIZE(exp_H_fm)))
               DO i = 1, SIZE(exp_H_fm)
                  CALL cp_fm_create(exp_H_fm_red(i)%matrix, matrix_struct=rtp%canon_ortho%oo_oo_fmstruct, name="exp_H_fm_red")
                  CALL cp_fm_to_fm_submat(exp_H_fm(i)%matrix, exp_H_fm_red(i)%matrix, &
                                          nao - ndep, nao - ndep, 1, 1, 1, 1)
               END DO
               CALL cp_fm_vect_dealloc(exp_H_fm)
               DO i = 1, SIZE(mos_old)/2
                  re = 2*i - 1
                  im = 2*i
                  CALL cp_fm_get_info(mos_old(re)%matrix, ncol_global=nmo)
                  CALL cp_gemm("N", "N", nao - ndep, nmo, nao - ndep, one, exp_H_fm_red(re)%matrix, mos_old(re)%matrix, &
                               zero, mos_next(re)%matrix)
                  CALL cp_gemm("N", "N", nao - ndep, nmo, nao - ndep, -one, exp_H_fm_red(im)%matrix, mos_old(im)%matrix, &
                               one, mos_next(re)%matrix)
                  CALL cp_gemm("N", "N", nao - ndep, nmo, nao - ndep, one, exp_H_fm_red(im)%matrix, mos_old(re)%matrix, &
                               zero, mos_next(im)%matrix)
                  CALL cp_gemm("N", "N", nao - ndep, nmo, nao - ndep, one, exp_H_fm_red(re)%matrix, mos_old(im)%matrix, &
                               one, mos_next(im)%matrix)
               END DO

               CALL cp_fm_vect_dealloc(exp_H_fm_red)
            ELSE
               DO i = 1, SIZE(mos_old)/2
                  re = 2*i - 1
                  im = 2*i
                  CALL cp_fm_to_fm(mos_old(re)%matrix, mos_next(re)%matrix)
                  CALL cp_fm_to_fm(mos_old(im)%matrix, mos_next(im)%matrix)
               END DO
            END IF
         END IF

         DO i = 1, SIZE(mos_new)/2
            re = 2*i - 1
            im = 2*i
            CALL cp_fm_get_info(mos_new(re)%matrix, ncol_global=nmo)
            CALL cp_gemm("N", "N", nao - ndep, nmo, nao - ndep, one, propagator_fm_red(re)%matrix, mos_next(re)%matrix, &
                         zero, mos_new(re)%matrix)
            CALL cp_gemm("N", "N", nao - ndep, nmo, nao - ndep, -one, propagator_fm_red(im)%matrix, mos_next(im)%matrix, &
                         one, mos_new(re)%matrix)
            CALL cp_gemm("N", "N", nao - ndep, nmo, nao - ndep, one, propagator_fm_red(im)%matrix, mos_next(re)%matrix, &
                         zero, mos_new(im)%matrix)
            CALL cp_gemm("N", "N", nao - ndep, nmo, nao - ndep, one, propagator_fm_red(re)%matrix, mos_next(im)%matrix, &
                         one, mos_new(im)%matrix)
         END DO
         CALL cp_fm_vect_dealloc(propagator_fm_red)
         CALL cp_fm_vect_dealloc(propagator_matrix_fm_red)
      END IF

      CALL timestop(handle)

   END SUBROUTINE propagate_exp

! **************************************************************************************************
!> \brief Propagation of the density matrix instead of the atomic orbitals
!>        via a matrix exponential
!> \param rtp ...
!> \param rtp_control ...
!> \author Samuel Andermatt (02.2014)
! **************************************************************************************************

   SUBROUTINE propagate_exp_density(rtp, rtp_control)

      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), POINTER                    :: rtp_control

      CHARACTER(len=*), PARAMETER :: routineN = 'propagate_exp_density'
      REAL(KIND=dp), PARAMETER                           :: one = 1.0_dp, zero = 0.0_dp

      INTEGER                                            :: handle, i, im, re
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: exp_H_new, exp_H_old, propagator_matrix, &
                                                            rho_new, rho_next, rho_old
      TYPE(dbcsr_type), POINTER                          :: tmp_im, tmp_re

      CALL timeset(routineN, handle)

      CALL get_rtp(rtp=rtp, propagator_matrix=propagator_matrix, exp_H_new=exp_H_new, &
                   exp_H_old=exp_H_old, rho_old=rho_old, rho_new=rho_new, rho_next=rho_next)

      CALL compute_exponential_sparse(exp_H_new, propagator_matrix, rtp_control, rtp)

      !I could store these matrices in the type
      NULLIFY (tmp_re)
      ALLOCATE (tmp_re)
      CALL dbcsr_create(tmp_re, template=propagator_matrix(1)%matrix, matrix_type="N")
      NULLIFY (tmp_im)
      ALLOCATE (tmp_im)
      CALL dbcsr_create(tmp_im, template=propagator_matrix(1)%matrix, matrix_type="N")

      DO i = 1, SIZE(exp_H_new)/2
         re = 2*i - 1
         im = 2*i
         !Save some work by computing the first half of the propagation only once in case of ETRS
         !For EM this matrix has to be the initial matrix, thus a copy is enough
         IF (rtp%iter == 1) THEN
            IF (rtp_control%propagator == do_etrs) THEN
               CALL cp_complex_dbcsr_gemm_3("N", "N", one, exp_H_old(re)%matrix, exp_H_old(im)%matrix, &
                                            rho_old(re)%matrix, rho_old(im)%matrix, zero, tmp_re, tmp_im, filter_eps=rtp%filter_eps)
               CALL cp_complex_dbcsr_gemm_3("N", "C", one, tmp_re, tmp_im, exp_H_old(re)%matrix, exp_H_old(im)%matrix, &
                                            zero, rho_next(re)%matrix, rho_next(im)%matrix, filter_eps=rtp%filter_eps)
            ELSE
               CALL dbcsr_copy(rho_next(re)%matrix, rho_old(re)%matrix)
               CALL dbcsr_copy(rho_next(im)%matrix, rho_old(im)%matrix)
            END IF
         END IF
         CALL cp_complex_dbcsr_gemm_3("N", "N", one, exp_H_new(re)%matrix, exp_H_new(im)%matrix, &
                                      rho_next(re)%matrix, rho_next(im)%matrix, zero, tmp_re, tmp_im, filter_eps=rtp%filter_eps)
         CALL cp_complex_dbcsr_gemm_3("N", "C", one, tmp_re, tmp_im, exp_H_new(re)%matrix, exp_H_new(im)%matrix, &
                                      zero, rho_new(re)%matrix, rho_new(im)%matrix, filter_eps=rtp%filter_eps)
      END DO

      CALL dbcsr_deallocate_matrix(tmp_re)
      CALL dbcsr_deallocate_matrix(tmp_im)

      CALL timestop(handle)

   END SUBROUTINE propagate_exp_density

! **************************************************************************************************
!> \brief computes U_prop*MOs using arnoldi subspace algorithm
!> \param rtp ...
!> \param rtp_control ...
!> \param ortho ...
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************

   SUBROUTINE propagate_arnoldi(rtp, rtp_control, ortho)
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), POINTER                    :: rtp_control
      LOGICAL, INTENT(IN)                                :: ortho

      CHARACTER(len=*), PARAMETER                        :: routineN = 'propagate_arnoldi'

      INTEGER                                            :: handle, i, im, ispin, nao, ndep, nspin, &
                                                            re
      REAL(dp)                                           :: eps_arnoldi, t
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mos_new, mos_next, mos_old, &
                                                            propagator_matrix_fm, &
                                                            propagator_matrix_fm_red
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: propagator_matrix

      CALL timeset(routineN, handle)

      NULLIFY (mos_new, mos_old, mos_next, propagator_matrix)
      IF (.NOT. ortho) THEN
         CALL get_rtp(rtp=rtp, dt=t, mos_new=mos_new, mos_old=mos_old, &
                      mos_next=mos_next, propagator_matrix=propagator_matrix)
      ELSE
         CALL get_rtp(rtp=rtp, dt=t, mos_ortho_new=mos_new, mos_ortho_old=mos_old, &
                      mos_ortho_next=mos_next, propagator_matrix_canon=propagator_matrix)
      END IF

      nspin = SIZE(mos_new)/2
      eps_arnoldi = rtp_control%eps_exp
      ! except for the first step the further propagated mos_next
      ! must be copied on mos_old so that we have the half propagated mos
      ! ready on mos_old and only need to perform the second half propagatioon
      IF (rtp_control%propagator == do_etrs .AND. rtp%iter == 1) THEN
         DO i = 1, SIZE(mos_new)
            CALL cp_fm_to_fm(mos_next(i)%matrix, mos_old(i)%matrix)
         END DO
      END IF

      NULLIFY (propagator_matrix_fm)
      ALLOCATE (propagator_matrix_fm(SIZE(propagator_matrix)))
      DO i = 1, SIZE(propagator_matrix)
         ALLOCATE (propagator_matrix_fm(i)%matrix)
         CALL cp_fm_create(propagator_matrix_fm(i)%matrix, &
                           matrix_struct=rtp%ao_ao_fmstruct, &
                           name="prop_fm")
         CALL copy_dbcsr_to_fm(propagator_matrix(i)%matrix, propagator_matrix_fm(i)%matrix)
      END DO

      ! need to reduce dimension of propagator_matrix_fm if linear dependencies in S were removed
      NULLIFY (propagator_matrix_fm_red)
      IF (ortho) THEN
         ndep = rtp_control%nlindep
         CALL cp_fm_get_info(propagator_matrix_fm(1)%matrix, nrow_global=nao)
         ALLOCATE (propagator_matrix_fm_red(SIZE(propagator_matrix_fm)))
         DO i = 1, SIZE(propagator_matrix_fm)
            CALL cp_fm_create(propagator_matrix_fm_red(i)%matrix, &
                              matrix_struct=rtp%canon_ortho%oo_oo_fmstruct, &
                              name="prop_fm_red")
            CALL cp_fm_to_fm_submat(propagator_matrix_fm(i)%matrix, propagator_matrix_fm_red(i)%matrix, &
                                    nao - ndep, nao - ndep, 1, 1, 1, 1)
         END DO
      END IF

      DO ispin = 1, nspin
         re = ispin*2 - 1
         im = ispin*2
         IF (.NOT. ortho) THEN
            IF (rtp_control%fixed_ions .AND. .NOT. rtp%do_hfx) THEN
               CALL arnoldi(mos_old(re:im), mos_new(re:im), &
                            eps_arnoldi, Him=propagator_matrix_fm(im)%matrix, &
                            mos_next=mos_next(re:im), narn_old=rtp%narn_old)
            ELSE
               CALL arnoldi(mos_old(re:im), mos_new(re:im), &
                            eps_arnoldi, Hre=propagator_matrix_fm(re)%matrix, &
                            Him=propagator_matrix_fm(im)%matrix, mos_next=mos_next(re:im), &
                            narn_old=rtp%narn_old)
            END IF
         ELSE
            IF (rtp_control%fixed_ions .AND. .NOT. rtp%do_hfx) THEN
               CALL arnoldi(mos_old(re:im), mos_new(re:im), &
                            eps_arnoldi, Him=propagator_matrix_fm_red(im)%matrix, &
                            mos_next=mos_next(re:im), narn_old=rtp%narn_old)
            ELSE
               CALL arnoldi(mos_old(re:im), mos_new(re:im), &
                            eps_arnoldi, Hre=propagator_matrix_fm_red(re)%matrix, &
                            Him=propagator_matrix_fm_red(im)%matrix, mos_next=mos_next(re:im), &
                            narn_old=rtp%narn_old)
            END IF
         END IF
      END DO

!    DO i=1,SIZE(propagator_matrix)
!         CALL copy_fm_to_dbcsr(propagator_matrix_fm(i)%matrix,propagator_matrix(i)%matrix)
!    END DO
      CALL cp_fm_vect_dealloc(propagator_matrix_fm)
      IF (ortho) CALL cp_fm_vect_dealloc(propagator_matrix_fm_red)

      CALL timestop(handle)

   END SUBROUTINE propagate_arnoldi

! **************************************************************************************************
!> \brief  Propagation using the Baker-Campbell-Hausdorff expansion,
!>         currently only works for rtp
!> \param rtp ...
!> \param rtp_control ...
!> \author Samuel Andermatt (02.2014)
! **************************************************************************************************

   SUBROUTINE propagate_bch(rtp, rtp_control)

      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), POINTER                    :: rtp_control

      CHARACTER(len=*), PARAMETER                        :: routineN = 'propagate_bch'

      INTEGER                                            :: handle, im, ispin, re
      REAL(dp)                                           :: dt
      REAL(KIND=dp)                                      :: prefac
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: exp_H_old, propagator_matrix, rho_new, &
                                                            rho_next, rho_old

      CALL timeset(routineN, handle)

      CALL get_rtp(rtp=rtp, propagator_matrix=propagator_matrix, rho_old=rho_old, rho_new=rho_new, &
                   rho_next=rho_next)

      DO ispin = 1, SIZE(propagator_matrix)/2
         re = 2*ispin - 1
         im = 2*ispin

         IF (rtp%iter == 1) THEN
            ! For EM I have to copy rho_old onto rho_next and for ETRS,
            ! this is the first term of the series of commutators that result in rho_next
            CALL dbcsr_copy(rho_next(re)%matrix, rho_old(re)%matrix)
            CALL dbcsr_copy(rho_next(im)%matrix, rho_old(im)%matrix)
            IF (rtp_control%propagator == do_etrs) THEN
               !since we never calculated the matrix exponential the old matrix exponential stores the unscalled propagator
               CALL get_rtp(rtp=rtp, exp_H_old=exp_H_old, dt=dt)
               prefac = -0.5_dp*dt
               CALL dbcsr_scale(exp_H_old(im)%matrix, prefac)
               IF (.NOT. rtp%do_hfx .AND. rtp_control%fixed_ions) THEN
                  CALL bch_expansion_imaginary_propagator( &
                     exp_H_old(im)%matrix, rho_next(re)%matrix, rho_next(im)%matrix, &
                     rtp%filter_eps, rtp%filter_eps_small, rtp_control%eps_exp)
               ELSE
                  CALL dbcsr_scale(exp_H_old(re)%matrix, prefac)
                  CALL bch_expansion_complex_propagator( &
                     exp_H_old(re)%matrix, exp_H_old(im)%matrix, rho_next(re)%matrix, rho_next(im)%matrix, &
                     rtp%filter_eps, rtp%filter_eps_small, rtp_control%eps_exp)
               END IF
            END IF
         END IF
         CALL dbcsr_copy(rho_new(re)%matrix, rho_next(re)%matrix)
         CALL dbcsr_copy(rho_new(im)%matrix, rho_next(im)%matrix)
         IF (.NOT. rtp%do_hfx .AND. rtp_control%fixed_ions) THEN
            CALL bch_expansion_imaginary_propagator( &
               propagator_matrix(im)%matrix, rho_new(re)%matrix, rho_new(im)%matrix, &
               rtp%filter_eps, rtp%filter_eps_small, rtp_control%eps_exp)
         ELSE
            CALL bch_expansion_complex_propagator( &
               propagator_matrix(re)%matrix, propagator_matrix(im)%matrix, rho_new(re)%matrix, rho_new(im)%matrix, &
               rtp%filter_eps, rtp%filter_eps_small, rtp_control%eps_exp)
         END IF

      END DO

      CALL timestop(handle)

   END SUBROUTINE propagate_bch

! **************************************************************************************************
!> \brief decides which type of exponential has to be computed
!> \param propagator ...
!> \param propagator_matrix ...
!> \param rtp_control ...
!> \param rtp ...
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************

   SUBROUTINE compute_exponential(propagator, propagator_matrix, rtp_control, rtp)
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: propagator, propagator_matrix
      TYPE(rtp_control_type), POINTER                    :: rtp_control
      TYPE(rt_prop_type), POINTER                        :: rtp

      INTEGER                                            :: im, ispin, re

      DO ispin = 1, SIZE(propagator)/2
         re = 2*ispin - 1
         im = 2*ispin
         SELECT CASE (rtp_control%mat_exp)
         CASE (do_taylor)
            IF (rtp_control%fixed_ions .AND. .NOT. rtp%do_hfx) THEN
               CALL taylor_only_imaginary(propagator(re:im), propagator_matrix(im)%matrix, &
                                          rtp%orders(1, ispin), rtp%orders(2, ispin))
            ELSE
               CALL taylor_full_complex(propagator(re:im), propagator_matrix(re)%matrix, propagator_matrix(im)%matrix, &
                                        rtp%orders(1, ispin), rtp%orders(2, ispin))
            END IF
         CASE (do_pade)
            IF (rtp_control%fixed_ions .AND. .NOT. rtp%do_hfx) THEN
               CALL exp_pade_only_imaginary(propagator(re:im), propagator_matrix(im)%matrix, &
                                            rtp%orders(1, ispin), rtp%orders(2, ispin))
            ELSE
               CALL exp_pade_full_complex(propagator(re:im), propagator_matrix(re)%matrix, propagator_matrix(im)%matrix, &
                                          rtp%orders(1, ispin), rtp%orders(2, ispin))
            END IF
         END SELECT
      END DO

   END SUBROUTINE compute_exponential

! **************************************************************************************************
!> \brief Sparse versions of the matrix exponentials
!> \param propagator ...
!> \param propagator_matrix ...
!> \param rtp_control ...
!> \param rtp ...
!> \author Samuel Andermatt (02.14)
! **************************************************************************************************

   SUBROUTINE compute_exponential_sparse(propagator, propagator_matrix, rtp_control, rtp)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: propagator, propagator_matrix
      TYPE(rtp_control_type), POINTER                    :: rtp_control
      TYPE(rt_prop_type), POINTER                        :: rtp

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_exponential_sparse'

      INTEGER                                            :: handle, im, ispin, re

      CALL timeset(routineN, handle)

      DO ispin = 1, SIZE(propagator)/2
         re = 2*ispin - 1
         im = 2*ispin
         IF (rtp_control%fixed_ions .AND. .NOT. rtp%do_hfx) THEN
            CALL taylor_only_imaginary_dbcsr(propagator(re:im), propagator_matrix(im)%matrix, &
                                             rtp%orders(1, ispin), rtp%orders(2, ispin), rtp%filter_eps)
         ELSE
            CALL taylor_full_complex_dbcsr(propagator(re:im), propagator_matrix(re)%matrix, propagator_matrix(im)%matrix, &
                                           rtp%orders(1, ispin), rtp%orders(2, ispin), rtp%filter_eps)
         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE compute_exponential_sparse

END MODULE rt_make_propagators
