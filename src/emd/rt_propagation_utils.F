!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines needed for EMD
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************

MODULE rt_propagation_utils
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              rtp_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale_and_add
   USE cp_fm_diag,                      ONLY: cp_fm_ortho_canon
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_binary_read, dbcsr_checksum, dbcsr_copy, dbcsr_copy_into_existing, &
        dbcsr_create, dbcsr_deallocate_matrix, dbcsr_desymmetrize, dbcsr_distribution_type, &
        dbcsr_filter, dbcsr_get_info, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_scale, dbcsr_set, dbcsr_type
   USE input_constants,                 ONLY: use_restart_wfn,&
                                              use_rt_restart
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE mathconstants,                   ONLY: one,&
                                              zero
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_density_matrices,             ONLY: calculate_density_matrix
   USE qs_dftb_matrices,                ONLY: build_dftb_overlap
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_did_change,&
                                              qs_ks_env_type
   USE qs_mo_io,                        ONLY: read_mo_set_from_restart,&
                                              read_rt_mos_from_restart
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE rt_propagation_types,            ONLY: get_rtp,&
                                              rt_prop_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: get_restart_wfn, &
             calc_S_derivs, &
             calc_update_rho, &
             calc_update_rho_sparse, &
             calculate_P_imaginary, &
             rt_canon_ortho, &
             transform_to_co_coeffs, &
             orthonormalize_complex_mos, &
             canon_ortho_matrix_ao_to_co, &
             canon_ortho_dbcsr_to_fm_reduced, &
             canon_ortho_fm_reduced_to_dbcsr

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_utils'

   INTERFACE canon_ortho_matrix_ao_to_co
      MODULE PROCEDURE canon_ortho_matrix_ao_to_co_p_type, canon_ortho_matrix_ao_to_co_type
   END INTERFACE

   INTERFACE transform_to_co_coeffs
      MODULE PROCEDURE transform_to_co_coeffs_p_type, transform_to_co_coeffs_type
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Calculates dS/dR respectily the velocity weighted derivatves
!>        only needed for ehrenfest MD.
!>
!> \param qs_env the qs environment
!> \par History
!>      02.2009 created [Manuel Guidon]
!>      02.2014 switched to dbcsr matrices [Samuel Andermatt]
!> \author Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE calc_S_derivs(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'calc_S_derivs'

      INTEGER                                            :: col_atom, handle, i, j, m, maxder, n, &
                                                            nder, row_atom
      INTEGER, DIMENSION(6, 2)                           :: c_map_mat
      LOGICAL                                            :: return_s_derivatives
      REAL(dp), DIMENSION(:), POINTER                    :: block_values
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: C_mat, S_der, s_derivs
      TYPE(dbcsr_type), POINTER                          :: B_mat, tmp_mat, tmp_mat2
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(rt_prop_type), POINTER                        :: rtp

! REAL(KIND=dp), PARAMETER                           :: one = 1.0_dp, zero = 0.0_dp

      CALL timeset(routineN, handle)

      return_s_derivatives = .TRUE.

      NULLIFY (particle_set)
      NULLIFY (rtp)
      NULLIFY (s_derivs)
      NULLIFY (dft_control)
      NULLIFY (ks_env)

      CALL get_qs_env(qs_env=qs_env, &
                      rtp=rtp, &
                      particle_set=particle_set, &
                      sab_orb=sab_orb, &
                      dft_control=dft_control, &
                      ks_env=ks_env)

      CALL get_rtp(rtp=rtp, B_mat=B_mat, C_mat=C_mat, S_der=S_der)

      nder = 2
      maxder = ncoset(nder)

      NULLIFY (tmp_mat)
      ALLOCATE (tmp_mat)
      CALL dbcsr_create(tmp_mat, template=S_der(1)%matrix, matrix_type="N")

      IF (rtp%iter < 2) THEN
         ! calculate the overlap derivative matrices
         IF (dft_control%qs_control%dftb) THEN
            CALL build_dftb_overlap(qs_env, nder, s_derivs)
         ELSE
            CALL build_overlap_matrix(ks_env, nderivative=nder, matrix_s=s_derivs, &
                                      basis_type_a="ORB", basis_type_b="ORB", sab_nl=sab_orb)
         END IF

         NULLIFY (tmp_mat2)
         ALLOCATE (tmp_mat2)
         CALL dbcsr_create(tmp_mat2, template=S_der(1)%matrix, matrix_type="S")
         DO m = 1, 9
            CALL dbcsr_copy(tmp_mat2, s_derivs(m + 1)%matrix)
            CALL dbcsr_desymmetrize(tmp_mat2, S_der(m)%matrix)
            CALL dbcsr_scale(S_der(m)%matrix, -one)
            CALL dbcsr_filter(S_der(m)%matrix, rtp%filter_eps)
            !The diagonal should be zero
            CALL dbcsr_iterator_start(iter, S_der(m)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
               IF (row_atom == col_atom) block_values = 0.0_dp
            END DO
            CALL dbcsr_iterator_stop(iter)
         END DO
         CALL dbcsr_deallocate_matrix_set(s_derivs)
         CALL dbcsr_deallocate_matrix(tmp_mat2)
      END IF

      !calculate scalar product v(Rb)*<alpha|d/dRb beta> (B_mat), and store the first derivatives

      CALL dbcsr_set(B_mat, zero)
      DO m = 1, 3
         CALL dbcsr_copy(tmp_mat, S_der(m)%matrix)
         CALL dbcsr_iterator_start(iter, tmp_mat)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
            IF (row_atom == col_atom) block_values = 0.0_dp
            block_values = block_values*particle_set(col_atom)%v(m)
         END DO
         CALL dbcsr_iterator_stop(iter)
         CALL dbcsr_add(B_mat, tmp_mat, one, one)
      END DO
      CALL dbcsr_filter(B_mat, rtp%filter_eps)
      !calculate C matrix: v(Rb)*<d/dRa alpha| d/dRb beta>

      c_map_mat = 0
      n = 0
      DO j = 1, 3
         DO m = j, 3
            n = n + 1
            c_map_mat(n, 1) = j
            IF (m == j) CYCLE
            c_map_mat(n, 2) = m
         END DO
      END DO

      DO i = 1, 3
         CALL dbcsr_set(C_mat(i)%matrix, zero)
      END DO
      DO m = 1, 6
         CALL dbcsr_copy(tmp_mat, S_der(m + 3)%matrix)
         DO j = 1, 2
            IF (c_map_mat(m, j) == 0) CYCLE
            CALL dbcsr_add(C_mat(c_map_mat(m, j))%matrix, tmp_mat, one, one)
         END DO
      END DO

      DO m = 1, 3
         CALL dbcsr_iterator_start(iter, C_mat(m)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
            block_values = block_values*particle_set(row_atom)%v(m)
         END DO
         CALL dbcsr_iterator_stop(iter)
         CALL dbcsr_filter(C_mat(m)%matrix, rtp%filter_eps)
      END DO

      CALL dbcsr_deallocate_matrix(tmp_mat)
      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief reads the restart file. At the moment only SCF (means only real)
!> \param qs_env ...
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************

   SUBROUTINE get_restart_wfn(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=default_path_length)                 :: file_name, project_name
      INTEGER                                            :: i, id_nr, im, ispin, ncol, nspin, re, &
                                                            unit_nr
      REAL(KIND=dp)                                      :: alpha, cs_pos
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mos_new, mos_old
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_rmpv, rho_new, rho_old
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_array
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(section_vals_type), POINTER                   :: dft_section, input

      NULLIFY (atomic_kind_set, qs_kind_set, mo_array, particle_set, rho_struct, para_env)

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, &
                      mos=mo_array, &
                      input=input, &
                      rtp=rtp, &
                      dft_control=dft_control, &
                      rho=rho_struct, &
                      para_env=para_env)
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      END IF

      id_nr = 0
      nspin = SIZE(mo_array)
      CALL qs_rho_get(rho_struct, rho_ao=p_rmpv)
      dft_section => section_vals_get_subs_vals(input, "DFT")
      SELECT CASE (dft_control%rtp_control%initial_wfn)
      CASE (use_restart_wfn)
         CALL read_mo_set_from_restart(mo_array, atomic_kind_set, qs_kind_set, particle_set, para_env, &
                                       id_nr=id_nr, multiplicity=dft_control%multiplicity, dft_section=dft_section)
         DO ispin = 1, nspin
            CALL calculate_density_matrix(mo_array(ispin)%mo_set, p_rmpv(ispin)%matrix)
         END DO
         IF (rtp%linear_scaling) THEN
            CALL get_rtp(rtp=rtp, rho_old=rho_old, rho_new=rho_new)
            DO ispin = 1, nspin
               re = 2*ispin - 1
               im = 2*ispin
               CALL cp_fm_get_info(mo_array(ispin)%mo_set%mo_coeff, ncol_global=ncol)
               alpha = 1.0_dp
               IF (SIZE(mo_array) == 1) alpha = 2*alpha
               CALL cp_dbcsr_plus_fm_fm_t( &
                  sparse_matrix=rho_old(re)%matrix, &
                  matrix_v=mo_array(ispin)%mo_set%mo_coeff, matrix_g=mo_array(ispin)%mo_set%mo_coeff, ncol=ncol, &
                  keep_sparsity=.FALSE., alpha=alpha)
            END DO
            DO i = 1, nspin
               CALL dbcsr_copy(rho_new(i)%matrix, rho_old(i)%matrix)
            END DO
            CALL calc_update_rho_sparse(qs_env)
         ELSE
            CALL get_rtp(rtp=rtp, mos_old=mos_old)
            DO i = 1, SIZE(qs_env%mos)
               CALL cp_fm_to_fm(mo_array(i)%mo_set%mo_coeff, mos_old(2*i - 1)%matrix)
               CALL cp_fm_set_all(mos_old(2*i)%matrix, zero, zero)
            END DO
         END IF
      CASE (use_rt_restart)
         IF (rtp%linear_scaling) THEN
            CALL get_rtp(rtp=rtp, rho_old=rho_old, rho_new=rho_new)
            project_name = logger%iter_info%project_name
            DO ispin = 1, nspin
               re = 2*ispin - 1
               im = 2*ispin
               WRITE (file_name, '(A,I0,A)') TRIM(project_name)//"_LS_DM_SPIN_RE", ispin, "_RESTART.dm"
               CALL dbcsr_get_info(rho_old(re)%matrix, distribution=dist)
               CALL dbcsr_binary_read(file_name, distribution=dist, matrix_new=rho_old(re)%matrix)
               cs_pos = dbcsr_checksum(rho_old(re)%matrix, pos=.TRUE.)
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, '(T2,A,E20.8)') "Read restart DM "//TRIM(file_name)//" with checksum: ", cs_pos
               END IF
               WRITE (file_name, '(A,I0,A)') TRIM(project_name)//"_LS_DM_SPIN_IM", ispin, "_RESTART.dm"
               CALL dbcsr_get_info(rho_old(im)%matrix, distribution=dist)
               CALL dbcsr_binary_read(file_name, distribution=dist, matrix_new=rho_old(im)%matrix)
               cs_pos = dbcsr_checksum(rho_old(im)%matrix, pos=.TRUE.)
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, '(T2,A,E20.8)') "Read restart DM "//TRIM(file_name)//" with checksum: ", cs_pos
               END IF
            END DO
            DO i = 1, SIZE(rho_new)
               CALL dbcsr_copy(rho_new(i)%matrix, rho_old(i)%matrix)
            END DO
            CALL calc_update_rho_sparse(qs_env)
         ELSE
            CALL get_rtp(rtp=rtp, mos_old=mos_old, mos_new=mos_new)
            CALL read_rt_mos_from_restart(mo_array, mos_old, atomic_kind_set, qs_kind_set, particle_set, para_env, &
                                          id_nr, dft_control%multiplicity, dft_section)
            DO ispin = 1, nspin
               CALL calculate_density_matrix(mo_array(ispin)%mo_set, &
                                             p_rmpv(ispin)%matrix)
            END DO
         END IF
      END SELECT

   END SUBROUTINE get_restart_wfn

! **************************************************************************************************
!> \brief calculates the density from the complex MOs and passes the density to
!>        qs_env.
!> \param qs_env ...
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************

   SUBROUTINE calc_update_rho(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'calc_update_rho'

      INTEGER                                            :: handle, i, im, ncol, re
      REAL(KIND=dp)                                      :: alpha
      TYPE(cp_fm_p_type)                                 :: mos_occ
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mos_ao, mos_new
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao, rho_ao_im, rho_co, rho_re
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), POINTER                    :: rtp_control

      CALL timeset(routineN, handle)

      NULLIFY (rho, ks_env, mos_new, rtp, dft_control, rtp_control)
      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      rho=rho, &
                      rtp=rtp, &
                      mos=mos, &
                      dft_control=dft_control)

      rtp_control => dft_control%rtp_control
      CALL get_rtp(rtp=rtp, mos_new=mos_new)
      NULLIFY (rho_ao, rho_co, rho_re)
      CALL qs_rho_get(rho_struct=rho, rho_ao=rho_ao)
      IF (rtp_control%canonical_ortho) THEN
         CALL get_rtp(rtp, P_canon=rho_co)
         rho_re => rho_co
      ELSE
         rho_re => rho_ao
      END IF
      DO i = 1, SIZE(mos_new)/2
         re = 2*i - 1; im = 2*i
         alpha = 3*one - REAL(SIZE(mos_new)/2, dp)
         CALL dbcsr_set(rho_re(i)%matrix, zero)
         CALL cp_fm_get_info(mos_new(re)%matrix, ncol_global=ncol)
         ALLOCATE (mos_occ%matrix)
         CALL cp_fm_create(mos_occ%matrix, &
                           matrix_struct=mos_new(re)%matrix%matrix_struct, &
                           name="mos_occ")
         CALL cp_fm_to_fm(mos_new(re)%matrix, mos_occ%matrix)
         CALL cp_fm_column_scale(mos_occ%matrix, mos(i)%mo_set%occupation_numbers/alpha)
         CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=rho_re(i)%matrix, &
                                    matrix_v=mos_occ%matrix, &
                                    ncol=ncol, &
                                    alpha=alpha)
         ! It is actually complex conjugate but i*i=-1 therefore it must be added
         CALL cp_fm_to_fm(mos_new(im)%matrix, mos_occ%matrix)
         CALL cp_fm_column_scale(mos_occ%matrix, mos(i)%mo_set%occupation_numbers/alpha)
         CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=rho_re(i)%matrix, &
                                    matrix_v=mos_occ%matrix, &
                                    ncol=ncol, &
                                    alpha=alpha)
         CALL cp_fm_release(mos_occ%matrix)
         ! calculate density matrix in AO basis from CO basis if required
         IF (rtp_control%canonical_ortho) THEN
            CALL dbcsr_set(rho_ao(i)%matrix, zero)
            CALL canon_ortho_matrix_co_to_ao_type(rtp, rtp_control, rho_re(i)%matrix, rho_ao(i)%matrix, density=.TRUE.)
         END IF
      END DO

      CALL qs_rho_update_rho(rho, qs_env)

      IF (rtp%do_hfx .OR. rtp%track_imag_density) THEN
         CALL qs_rho_get(rho_struct=rho, rho_ao_im=rho_ao_im)
         CALL calculate_P_imaginary(qs_env, rtp, rho_ao_im, keep_sparsity=.TRUE.)
         CALL qs_rho_set(rho, rho_ao_im=rho_ao_im)
      END IF

      ! transform mo_coeffs to AO basis for admm fit
      IF (rtp_control%canonical_ortho) THEN
         IF (dft_control%do_admm) THEN
            NULLIFY (mos_ao)
            CALL get_rtp(rtp, mos_ao=mos_ao)
            CALL transform_to_ao_coeffs(qs_env, mos_new, mos_ao)
         END IF
      END IF
      CALL qs_ks_did_change(ks_env, rho_changed=.TRUE.)

      CALL timestop(handle)

   END SUBROUTINE calc_update_rho

! **************************************************************************************************
!> \brief Copies the density matrix back into the qs_env%rho%rho_ao
!> \param qs_env ...
!> \author Samuel Andermatt (3.14)
! **************************************************************************************************

   SUBROUTINE calc_update_rho_sparse(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calc_update_rho_sparse'

      INTEGER                                            :: handle, ispin
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao, rho_ao_im, rho_new
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), POINTER                    :: rtp_control

! REAL(KIND=dp), PARAMETER                           :: zero = 0.0_dp

      NULLIFY (rho, ks_env, rtp, dft_control)
      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      rho=rho, &
                      rtp=rtp, &
                      dft_control=dft_control)
      rtp_control => dft_control%rtp_control
      CALL get_rtp(rtp=rtp, rho_new=rho_new)
      CALL qs_rho_get(rho_struct=rho, rho_ao=rho_ao)
      IF (rtp%do_hfx .OR. rtp%track_imag_density) CALL qs_rho_get(rho_struct=rho, rho_ao_im=rho_ao_im)
      DO ispin = 1, SIZE(rho_ao)
         CALL dbcsr_set(rho_ao(ispin)%matrix, zero)
         CALL dbcsr_copy_into_existing(rho_ao(ispin)%matrix, rho_new(ispin*2 - 1)%matrix)
         IF (rtp%do_hfx .OR. rtp%track_imag_density) THEN
            CALL dbcsr_copy_into_existing(rho_ao_im(ispin)%matrix, rho_new(ispin*2)%matrix)
         END IF
      END DO

      CALL qs_rho_update_rho(rho, qs_env)
      CALL qs_ks_did_change(ks_env, rho_changed=.TRUE.)

      CALL timestop(handle)

   END SUBROUTINE calc_update_rho_sparse

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rtp ...
!> \param matrix_p_im ...
!> \param keep_sparsity ...
! **************************************************************************************************
   SUBROUTINE calculate_P_imaginary(qs_env, rtp, matrix_p_im, keep_sparsity)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(OUT), &
         POINTER                                         :: matrix_p_im
      LOGICAL, INTENT(IN), OPTIONAL                      :: keep_sparsity

      INTEGER                                            :: i, im, ncol, re
      LOGICAL                                            :: my_keep_sparsity
      REAL(KIND=dp)                                      :: alpha
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mos_new, mos_occ
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p_im_p
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(rtp_control_type), POINTER                    :: rtp_control

      NULLIFY (mos_new, dft_control, mos, rtp_control, matrix_p_im_p)
      my_keep_sparsity = .FALSE.
      IF (PRESENT(keep_sparsity)) my_keep_sparsity = keep_sparsity

      CALL get_rtp(rtp=rtp, mos_new=mos_new)
      CALL get_qs_env(qs_env, mos=mos, dft_control=dft_control)
      rtp_control => dft_control%rtp_control

      IF (.NOT. rtp_control%canonical_ortho) THEN
         matrix_p_im_p => matrix_p_im
      ELSE
         CALL get_rtp(rtp, P_canon_im=matrix_p_im_p)
      END IF

      ALLOCATE (mos_occ(SIZE(mos)*2))
      DO i = 1, SIZE(mos_new)/2
         re = 2*i - 1; im = 2*i
         alpha = 3.0_dp - REAL(SIZE(matrix_p_im), dp)
         ALLOCATE (mos_occ(re)%matrix, mos_occ(im)%matrix)
         CALL cp_fm_create(mos_occ(re)%matrix, &
                           matrix_struct=mos_new(re)%matrix%matrix_struct, &
                           name="mos_occ")
         CALL cp_fm_create(mos_occ(im)%matrix, &
                           matrix_struct=mos_new(im)%matrix%matrix_struct, &
                           name="mos_occ")
         CALL dbcsr_set(matrix_p_im_p(i)%matrix, 0.0_dp)
         CALL cp_fm_get_info(mos_new(re)%matrix, ncol_global=ncol)
         CALL cp_fm_to_fm(mos_new(re)%matrix, mos_occ(re)%matrix)
         CALL cp_fm_column_scale(mos_occ(re)%matrix, mos(i)%mo_set%occupation_numbers/alpha)
         CALL cp_fm_to_fm(mos_new(im)%matrix, mos_occ(im)%matrix)
         CALL cp_fm_column_scale(mos_occ(im)%matrix, mos(i)%mo_set%occupation_numbers/alpha)
         CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=matrix_p_im_p(i)%matrix, &
                                    matrix_v=mos_occ(im)%matrix, &
                                    matrix_g=mos_occ(re)%matrix, &
                                    ncol=ncol, &
                                    keep_sparsity=my_keep_sparsity, &
                                    alpha=alpha)
         ! It is actually complex conjugate not only transposed
         alpha = -alpha
         CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=matrix_p_im_p(i)%matrix, &
                                    matrix_v=mos_occ(re)%matrix, &
                                    matrix_g=mos_occ(im)%matrix, &
                                    ncol=ncol, &
                                    keep_sparsity=my_keep_sparsity, &
                                    alpha=alpha)
         IF (rtp_control%canonical_ortho) THEN
            CALL dbcsr_set(matrix_p_im(i)%matrix, zero)
            CALL canon_ortho_matrix_co_to_ao_type(rtp, rtp_control, matrix_p_im_p(i)%matrix, matrix_p_im(i)%matrix, &
                                                  density=.TRUE.)
         END IF
      END DO
      CALL cp_fm_vect_dealloc(mos_occ)

   END SUBROUTINE calculate_P_imaginary

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calc_trsfmatrices ...
! **************************************************************************************************
   SUBROUTINE rt_canon_ortho(qs_env, calc_trsfmatrices)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calc_trsfmatrices

      INTEGER                                            :: nao, ndep, output_unit
      REAL(kind=dp)                                      :: epseigval
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: mat_tmp, mat_tmp2, s_fm
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), POINTER                    :: rtp_control
      TYPE(section_vals_type), POINTER                   :: input, scf_section

      NULLIFY (matrix_s, para_env, blacs_env, input, scf_section, rtp_control, dft_control, rtp, logger)
      CALL get_qs_env(qs_env, &
                      blacs_env=blacs_env, &
                      dft_control=dft_control, &
                      input=input, &
                      matrix_s=matrix_s, &
                      para_env=para_env, &
                      rtp=rtp)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)

      scf_section => section_vals_get_subs_vals(input, "DFT%SCF")
      CALL section_vals_val_get(scf_section, "EPS_EIGVAL", r_val=epseigval)

      ! create fm struct for overlap matrix
      CALL cp_fm_struct_create(fm_struct_tmp, &
                               nrow_global=nao, &
                               ncol_global=nao, &
                               context=blacs_env, &
                               para_env=para_env)
      rtp_control => dft_control%rtp_control
      CALL cp_fm_create(s_fm, matrix_struct=fm_struct_tmp, name="overlap fm")
      CALL cp_fm_create(mat_tmp, matrix_struct=fm_struct_tmp, name="temporary fm")
      IF (calc_trsfmatrices) CALL cp_fm_create(mat_tmp2, matrix_struct=fm_struct_tmp, name="temporary fm")
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, s_fm)

      IF (.NOT. calc_trsfmatrices) THEN
         CALL cp_fm_ortho_canon(s_fm, mat_tmp, ndep, epseigval)
      ELSE
         CALL cp_fm_ortho_canon(s_fm, mat_tmp, ndep, epseigval, ymatrix=mat_tmp2)
      END IF
      rtp_control%nlindep = ndep
      CALL cp_fm_release(s_fm)

      IF (calc_trsfmatrices) THEN
         CPASSERT(ASSOCIATED(rtp%canon_ortho))
         CALL cp_fm_to_fm(mat_tmp, rtp%canon_ortho%xmatrix, nao - ndep, 1, 1)
         CALL cp_fm_to_fm(mat_tmp2, rtp%canon_ortho%ymatrix, nao - ndep, 1, 1)

         ! copy to dbcsr counter parts
         CALL dbcsr_set(rtp%canon_ortho%xmatrix_dbcsr, zero)
         CALL dbcsr_set(rtp%canon_ortho%ymatrix_dbcsr, zero)
         CALL copy_fm_to_dbcsr(mat_tmp, rtp%canon_ortho%xmatrix_dbcsr)
         CALL copy_fm_to_dbcsr(mat_tmp2, rtp%canon_ortho%ymatrix_dbcsr)
         ! call dbcsr_print(rtp%canon_ortho%ymatrix_dbcsr, matlab_format=.true.)

         logger => cp_get_default_logger()
         scf_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION")
         output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%PROGRAM_RUN_INFO", &
                                            extension=".scfLog")
         IF (output_unit > 0) THEN
            WRITE (output_unit, FMT="(/,(T3,A,T76,I5))") &
               "Number of truncated linear dependent basis functions for RTP:", ndep
         END IF
      END IF

      CALL cp_fm_release(mat_tmp)
      IF (calc_trsfmatrices) CALL cp_fm_release(mat_tmp2)

   END SUBROUTINE rt_canon_ortho

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param coeffs ...
!> \param coeffs_ortho ...
! **************************************************************************************************
   SUBROUTINE transform_to_co_coeffs_p_type(qs_env, coeffs, coeffs_ortho)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: coeffs
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(OUT), &
         POINTER                                         :: coeffs_ortho

      INTEGER                                            :: ispin, j, nao, nmo, noo
      TYPE(cp_fm_type), POINTER                          :: ymatrix
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(rt_prop_type), POINTER                        :: rtp

      NULLIFY (dft_control, matrix_s, rtp, ymatrix)
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s, &
                      rtp=rtp)

      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)
      noo = nao - dft_control%rtp_control%nlindep

      CALL get_rtp(rtp, ymatrix=ymatrix)
      CPASSERT(ASSOCIATED(ymatrix))

      DO ispin = 1, dft_control%nspins ! spins
         CALL cp_fm_get_info(coeffs(2*ispin)%matrix, ncol_global=nmo)
         DO j = 1, 2               ! real and imaginary part
            CALL cp_gemm("T", "N", noo, nmo, nao, 1.0_dp, ymatrix, coeffs(2*(ispin - 1) + j)%matrix, 0.0_dp, &
                         coeffs_ortho(2*(ispin - 1) + j)%matrix)
         END DO
      END DO

      ! normalize as ymatrix does NOT correspond to a unitary transformation
      CALL orthonormalize_complex_mos(qs_env, coeffs_ortho, orthogonal_basis=.TRUE.)

   END SUBROUTINE transform_to_co_coeffs_p_type

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param coeffs ...
!> \param coeffs_ortho ...
! **************************************************************************************************
   SUBROUTINE transform_to_co_coeffs_type(qs_env, coeffs, coeffs_ortho)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: coeffs
      TYPE(cp_fm_type), INTENT(OUT), POINTER             :: coeffs_ortho

      INTEGER                                            :: nao, nmo, noo
      TYPE(cp_fm_type), POINTER                          :: ymatrix
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(rt_prop_type), POINTER                        :: rtp

      NULLIFY (dft_control, matrix_s, rtp, ymatrix)
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s, &
                      rtp=rtp)

      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)
      noo = nao - dft_control%rtp_control%nlindep

      CALL get_rtp(rtp, ymatrix=ymatrix)
      CPASSERT(ASSOCIATED(ymatrix))

      CALL cp_fm_get_info(coeffs, ncol_global=nmo)
      CALL cp_gemm("T", "N", noo, nmo, nao, 1.0_dp, ymatrix, coeffs, 0.0_dp, &
                   coeffs_ortho)

   END SUBROUTINE transform_to_co_coeffs_type

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param coeffs_ortho ...
!> \param coeffs_ao ...
! **************************************************************************************************
   SUBROUTINE transform_to_ao_coeffs(qs_env, coeffs_ortho, coeffs_ao)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: coeffs_ortho
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(OUT), &
         POINTER                                         :: coeffs_ao

      INTEGER                                            :: ispin, j, nao, nmo, noo
      TYPE(cp_fm_type), POINTER                          :: xmatrix
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(rt_prop_type), POINTER                        :: rtp

      NULLIFY (dft_control, matrix_s, rtp)
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s, &
                      rtp=rtp)

      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)
      noo = nao - dft_control%rtp_control%nlindep

      CALL get_rtp(rtp, xmatrix=xmatrix)
      CPASSERT(ASSOCIATED(xmatrix))

      DO ispin = 1, dft_control%nspins ! spins
         CALL cp_fm_get_info(coeffs_ao(2*ispin)%matrix, ncol_global=nmo)
         DO j = 1, 2               ! real and imaginary part
            CALL cp_gemm("N", "N", nao, nmo, noo, 1.0_dp, xmatrix, coeffs_ortho(2*(ispin - 1) + j)%matrix, 0.0_dp, &
                         coeffs_ao(2*(ispin - 1) + j)%matrix)
         END DO
      END DO

      ! normalize as ymatrix does NOT correspond to a unitary transformation
      CALL orthonormalize_complex_mos(qs_env, coeffs_ao, orthogonal_basis=.FALSE.)

   END SUBROUTINE transform_to_ao_coeffs

! **************************************************************************************************
!> \brief orthonormalize complex mos, e. g. after non-unitary transformations using Löwdin's algorithm
!>        refactored  and extended from rt_delta_pulse.F
!> \param qs_env ...
!> \param coeffs ...
!> \param orthogonal_basis set to .TRUE. if the basis functions are orthogonal, i. e. S = identity
! **************************************************************************************************
   SUBROUTINE orthonormalize_complex_mos(qs_env, coeffs, orthogonal_basis)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: coeffs
      LOGICAL, INTENT(IN)                                :: orthogonal_basis

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: eigenvalues_sqrt
      INTEGER                                            :: im, ispin, j, nao, nmo, nspins, re
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type)                                  :: oo_c, oo_v, oo_vt
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: coeffs_tmp
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type)                                   :: oo_1, oo_2, S_fm, tmp
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      NULLIFY (para_env, blacs_env, dft_control, matrix_s, mos)
      CALL get_qs_env(qs_env, &
                      blacs_env=blacs_env, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s, &
                      mos=mos, &
                      para_env=para_env)
      nspins = dft_control%nspins
      CALL cp_fm_get_info(coeffs(1)%matrix, nrow_global=nao)

      ! get overlap matrix
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nao, &
                               context=blacs_env, para_env=para_env)
      CALL cp_fm_create(S_fm, matrix_struct=fm_struct_tmp, name="overlap fm")
      CALL cp_fm_struct_release(fm_struct_tmp)
      IF (.NOT. orthogonal_basis) THEN
         ! copy overlap matrix
         CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, S_fm)
      ELSE
         ! set overlap matrix to identity
         CALL cp_fm_set_all(S_fm, zero, one)
      END IF

      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin)%mo_set, nmo=nmo)
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                                  nrow_global=nmo, ncol_global=nmo)
         CALL cp_fm_create(oo_1, matrix_struct=fm_struct_tmp, name="oo_1")
         CALL cp_fm_create(oo_2, matrix_struct=fm_struct_tmp, name="oo_2")
         CALL cp_fm_struct_release(fm_struct_tmp)

         CALL cp_fm_create(tmp, matrix_struct=coeffs(2*ispin - 1)%matrix%matrix_struct, name="tmp_mat")
         ! get the complex overlap matrix in MO basis
         ! x^T S x + y^T S y + i (-y^TS x+x^T S y)
         CALL parallel_gemm("N", "N", nao, nmo, nao, 1.0_dp, S_fm, coeffs(2*ispin - 1)%matrix, 0.0_dp, tmp)
         CALL parallel_gemm("T", "N", nmo, nmo, nao, 1.0_dp, coeffs(2*ispin - 1)%matrix, tmp, 0.0_dp, oo_1)
         CALL parallel_gemm("T", "N", nmo, nmo, nao, -1.0_dp, coeffs(2*ispin)%matrix, tmp, 0.0_dp, oo_2)

         CALL parallel_gemm("N", "N", nao, nmo, nao, 1.0_dp, S_fm, coeffs(2*ispin)%matrix, 0.0_dp, tmp)
         CALL parallel_gemm("T", "N", nmo, nmo, nao, 1.0_dp, coeffs(2*ispin)%matrix, tmp, 1.0_dp, oo_1)
         CALL parallel_gemm("T", "N", nmo, nmo, nao, 1.0_dp, coeffs(2*ispin - 1)%matrix, tmp, 1.0_dp, oo_2)
         CALL cp_fm_release(tmp)

         ! complex Löwdin
         CALL cp_cfm_create(oo_c, oo_1%matrix_struct)
         CALL cp_cfm_create(oo_v, oo_1%matrix_struct)
         CALL cp_cfm_create(oo_vt, oo_1%matrix_struct)
         oo_c%local_data = CMPLX(oo_1%local_data, oo_2%local_data, KIND=dp)

         ALLOCATE (eigenvalues(nmo))
         ALLOCATE (eigenvalues_sqrt(nmo))
         CALL cp_cfm_heevd(oo_c, oo_v, eigenvalues)
         eigenvalues_sqrt(:) = CMPLX(one/SQRT(eigenvalues(:)), zero, dp)
         CALL cp_cfm_to_cfm(oo_v, oo_vt)
         CALL cp_cfm_column_scale(oo_v, eigenvalues_sqrt)
         DEALLOCATE (eigenvalues)
         DEALLOCATE (eigenvalues_sqrt)
         CALL parallel_gemm('N', 'C', nmo, nmo, nmo, (1.0_dp, 0.0_dp), &
                            oo_v, oo_vt, (0.0_dp, 0.0_dp), oo_c)
         oo_1%local_data = REAL(oo_c%local_data, KIND=dp)
         oo_2%local_data = AIMAG(oo_c%local_data)
         CALL cp_cfm_release(oo_c)
         CALL cp_cfm_release(oo_v)
         CALL cp_cfm_release(oo_vt)

         ! transform coefficients accordingly
         NULLIFY (coeffs_tmp)
         ALLOCATE (coeffs_tmp(2))
         DO j = 1, 2
            NULLIFY (coeffs_tmp(j)%matrix)
            ALLOCATE (coeffs_tmp(j)%matrix)
            CALL cp_fm_create(coeffs_tmp(j)%matrix, matrix_struct=coeffs(2*(ispin - 1) + j)%matrix%matrix_struct)
         END DO

         ! indices for coeffs_tmp
         re = 1
         im = 2
         CALL parallel_gemm("N", "N", nao, nmo, nmo, one, coeffs(2*ispin - 1)%matrix, oo_1, zero, coeffs_tmp(re)%matrix)
         CALL parallel_gemm("N", "N", nao, nmo, nmo, one, coeffs(2*ispin - 1)%matrix, oo_2, zero, coeffs_tmp(im)%matrix)

         CALL parallel_gemm("N", "N", nao, nmo, nmo, -one, coeffs(2*ispin)%matrix, oo_2, zero, coeffs(2*ispin - 1)%matrix)
         CALL cp_fm_scale_and_add(one, coeffs(2*ispin - 1)%matrix, one, coeffs_tmp(re)%matrix)

         CALL parallel_gemm("N", "N", nao, nmo, nmo, one, coeffs(2*ispin)%matrix, oo_1, one, coeffs_tmp(im)%matrix)
         CALL cp_fm_to_fm(coeffs_tmp(im)%matrix, coeffs(2*ispin)%matrix)

         CALL cp_fm_vect_dealloc(coeffs_tmp)
         CALL cp_fm_release(oo_1)
         CALL cp_fm_release(oo_2)
      END DO
      CALL cp_fm_release(S_fm)

   END SUBROUTINE orthonormalize_complex_mos

! **************************************************************************************************
!> \brief ...
!> \param rtp ...
!> \param rtp_control ...
!> \param matrix_ao ...
!> \param matrix_co ...
!> \param density ...
! **************************************************************************************************
   SUBROUTINE canon_ortho_matrix_ao_to_co_p_type(rtp, rtp_control, matrix_ao, matrix_co, density)
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), INTENT(IN), POINTER        :: rtp_control
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: matrix_ao
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(OUT), &
         POINTER                                         :: matrix_co
      LOGICAL, INTENT(IN)                                :: density

      INTEGER                                            :: i, nao, noo
      TYPE(dbcsr_type), POINTER                          :: trsfmat, work, xmatrix_dbcsr, &
                                                            ymatrix_dbcsr

      NULLIFY (xmatrix_dbcsr, ymatrix_dbcsr, trsfmat, work)

      CALL dbcsr_get_info(matrix_ao(1)%matrix, nfullrows_total=nao)
      noo = nao - rtp_control%nlindep
      CALL get_rtp(rtp, xmatrix_dbcsr=xmatrix_dbcsr, ymatrix_dbcsr=ymatrix_dbcsr, work=work)

      CPASSERT(ASSOCIATED(xmatrix_dbcsr) .AND. ASSOCIATED(ymatrix_dbcsr))
      IF (density) THEN
         trsfmat => ymatrix_dbcsr
      ELSE
         trsfmat => xmatrix_dbcsr
      END IF

      DO i = 1, SIZE(matrix_ao)
         CALL dbcsr_set(work, zero)
         CALL dbcsr_set(matrix_co(i)%matrix, zero)
         CALL dbcsr_multiply("N", "N", one, matrix_ao(i)%matrix, trsfmat, zero, work, 1, nao, 1, noo)
         CALL dbcsr_multiply("T", "N", one, trsfmat, work, zero, matrix_co(i)%matrix, 1, noo, 1, noo)
      END DO

   END SUBROUTINE canon_ortho_matrix_ao_to_co_p_type

! **************************************************************************************************
!> \brief ...
!> \param rtp ...
!> \param rtp_control ...
!> \param matrix_ao ...
!> \param matrix_co ...
!> \param density ...
! **************************************************************************************************
   SUBROUTINE canon_ortho_matrix_ao_to_co_type(rtp, rtp_control, matrix_ao, matrix_co, density)
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), INTENT(IN), POINTER        :: rtp_control
      TYPE(dbcsr_type), INTENT(IN), POINTER              :: matrix_ao
      TYPE(dbcsr_type), INTENT(OUT), POINTER             :: matrix_co
      LOGICAL, INTENT(IN)                                :: density

      INTEGER                                            :: nao, noo
      TYPE(dbcsr_type), POINTER                          :: trsfmat, work, xmatrix_dbcsr, &
                                                            ymatrix_dbcsr

      NULLIFY (xmatrix_dbcsr, ymatrix_dbcsr, trsfmat, work)

      CALL dbcsr_get_info(matrix_ao, nfullrows_total=nao)
      noo = nao - rtp_control%nlindep

      CALL get_rtp(rtp, xmatrix_dbcsr=xmatrix_dbcsr, ymatrix_dbcsr=ymatrix_dbcsr, work=work)
      CPASSERT(ASSOCIATED(xmatrix_dbcsr) .AND. ASSOCIATED(ymatrix_dbcsr))
      IF (density) THEN
         trsfmat => ymatrix_dbcsr
      ELSE
         trsfmat => xmatrix_dbcsr
      END IF

      CALL dbcsr_set(work, zero)
      CALL dbcsr_set(matrix_co, zero)
      CALL dbcsr_multiply("N", "N", one, matrix_ao, trsfmat, zero, work, 1, nao, 1, noo)
      CALL dbcsr_multiply("T", "N", one, trsfmat, work, zero, matrix_co, 1, noo, 1, noo)

   END SUBROUTINE canon_ortho_matrix_ao_to_co_type

! **************************************************************************************************
!> \brief ...
!> \param rtp ...
!> \param rtp_control ...
!> \param matrix_co ...
!> \param matrix_ao ...
!> \param density ...
! **************************************************************************************************
   SUBROUTINE canon_ortho_matrix_co_to_ao_p_type(rtp, rtp_control, matrix_co, matrix_ao, density)
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), INTENT(IN), POINTER        :: rtp_control
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: matrix_co
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(OUT), &
         POINTER                                         :: matrix_ao
      LOGICAL, INTENT(IN)                                :: density

      INTEGER                                            :: i, nao, noo
      TYPE(dbcsr_type), POINTER                          :: trsfmat, work, xmatrix_dbcsr, &
                                                            ymatrix_dbcsr

      NULLIFY (xmatrix_dbcsr, ymatrix_dbcsr, trsfmat, work)

      CALL dbcsr_get_info(matrix_ao(1)%matrix, nfullrows_total=nao)
      noo = nao - rtp_control%nlindep
      CALL get_rtp(rtp, xmatrix_dbcsr=xmatrix_dbcsr, ymatrix_dbcsr=ymatrix_dbcsr, work=work)

      CPASSERT(ASSOCIATED(xmatrix_dbcsr) .AND. ASSOCIATED(ymatrix_dbcsr))
      IF (density) THEN
         trsfmat => xmatrix_dbcsr
      ELSE
         trsfmat => ymatrix_dbcsr
      END IF

      DO i = 1, SIZE(matrix_ao)
         CALL dbcsr_set(work, zero)
         CALL dbcsr_multiply("N", "T", one, matrix_co(i)%matrix, trsfmat, zero, work, 1, noo, 1, nao)
         CALL dbcsr_multiply("N", "N", one, trsfmat, work, zero, matrix_ao(i)%matrix, 1, nao, 1, nao)
      END DO

   END SUBROUTINE canon_ortho_matrix_co_to_ao_p_type

! **************************************************************************************************
!> \brief ...
!> \param rtp ...
!> \param rtp_control ...
!> \param matrix_co ...
!> \param matrix_ao ...
!> \param density ...
! **************************************************************************************************
   SUBROUTINE canon_ortho_matrix_co_to_ao_type(rtp, rtp_control, matrix_co, matrix_ao, density)
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), INTENT(IN), POINTER        :: rtp_control
      TYPE(dbcsr_type), INTENT(IN), POINTER              :: matrix_co
      TYPE(dbcsr_type), INTENT(OUT), POINTER             :: matrix_ao
      LOGICAL, INTENT(IN)                                :: density

      INTEGER                                            :: nao, noo
      TYPE(dbcsr_type), POINTER                          :: trsfmat, work, xmatrix_dbcsr, &
                                                            ymatrix_dbcsr

      NULLIFY (xmatrix_dbcsr, ymatrix_dbcsr, trsfmat, work)

      CALL dbcsr_get_info(matrix_ao, nfullrows_total=nao)
      noo = nao - rtp_control%nlindep
      CALL get_rtp(rtp, xmatrix_dbcsr=xmatrix_dbcsr, ymatrix_dbcsr=ymatrix_dbcsr, work=work)

      CPASSERT(ASSOCIATED(xmatrix_dbcsr) .AND. ASSOCIATED(ymatrix_dbcsr))
      IF (density) THEN
         trsfmat => xmatrix_dbcsr
      ELSE
         trsfmat => ymatrix_dbcsr
      END IF

      CALL dbcsr_set(work, zero)
      CALL dbcsr_multiply("N", "T", one, matrix_co, trsfmat, zero, work, 1, noo, 1, nao)
      CALL dbcsr_multiply("N", "N", one, trsfmat, work, zero, matrix_ao, 1, nao, 1, nao)

   END SUBROUTINE canon_ortho_matrix_co_to_ao_type

! **************************************************************************************************
!> \brief transform dbcsr matrix to fm matrix with possibly reduced dimensions
!> \param rtp ...
!> \param rtp_control ...
!> \param dbcsr_matrix ...
!> \param fm_matrix_red ...
! **************************************************************************************************
   SUBROUTINE canon_ortho_dbcsr_to_fm_reduced(rtp, rtp_control, dbcsr_matrix, fm_matrix_red)
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), INTENT(IN), POINTER        :: rtp_control
      TYPE(dbcsr_type), INTENT(IN), POINTER              :: dbcsr_matrix
      TYPE(cp_fm_type), INTENT(OUT), POINTER             :: fm_matrix_red

      INTEGER                                            :: nao, noo
      TYPE(cp_fm_type), POINTER                          :: work_fm

      CALL dbcsr_get_info(rtp%exp_H_new(1)%matrix, nfullrows_total=nao)
      noo = nao - rtp_control%nlindep

      NULLIFY (work_fm)
      IF (noo .NE. nao) THEN
         CALL get_rtp(rtp, work_fm=work_fm)
         CALL cp_fm_set_all(work_fm, 0.0_dp)
         CALL copy_dbcsr_to_fm(dbcsr_matrix, work_fm)
         CALL cp_fm_to_fm_submat(work_fm, fm_matrix_red, noo, noo, 1, 1, 1, 1)
      ELSE
         CALL copy_dbcsr_to_fm(dbcsr_matrix, fm_matrix_red)
      END IF

   END SUBROUTINE canon_ortho_dbcsr_to_fm_reduced

! **************************************************************************************************
!> \brief transform dbcsr matrix to fm matrix with possibly reduced dimensions
!> \param rtp ...
!> \param rtp_control ...
!> \param fm_matrix_red ...
!> \param dbcsr_matrix ...
! **************************************************************************************************
   SUBROUTINE canon_ortho_fm_reduced_to_dbcsr(rtp, rtp_control, fm_matrix_red, dbcsr_matrix)
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), INTENT(IN), POINTER        :: rtp_control
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: fm_matrix_red
      TYPE(dbcsr_type), INTENT(OUT), POINTER             :: dbcsr_matrix

      INTEGER                                            :: nao, noo
      TYPE(cp_fm_type), POINTER                          :: work_fm

      CALL dbcsr_get_info(rtp%exp_H_new(1)%matrix, nfullrows_total=nao)
      noo = nao - rtp_control%nlindep

      NULLIFY (work_fm)
      IF (noo .NE. nao) THEN
         CALL get_rtp(rtp, work_fm=work_fm)
         CALL cp_fm_set_all(work_fm, 0.0_dp)
         CALL cp_fm_to_fm_submat(fm_matrix_red, work_fm, noo, noo, 1, 1, 1, 1)
         CALL copy_fm_to_dbcsr(work_fm, dbcsr_matrix)
      ELSE
         CALL copy_fm_to_dbcsr(fm_matrix_red, dbcsr_matrix)
      END IF

   END SUBROUTINE canon_ortho_fm_reduced_to_dbcsr

END MODULE rt_propagation_utils
