!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Polarizability calculation by dfpt
!>      Initialization of the polar_env,
!>      Perturbation Hamiltonian by application of the Berry phase operator to psi0
!>      Write output
!>      Deallocate everything
!> periodic Raman SL February 2013
!> \note
! **************************************************************************************************
MODULE qs_linres_polar_utils
   USE bibliography,                    ONLY: Luber2014,&
                                              cite_reference
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type, &
        dbcsr_copy, dbcsr_create, dbcsr_dot, dbcsr_print, dbcsr_release, &
        dbcsr_scale, dbcsr_set, dbcsr_type_no_symmetry, dbcsr_type_antisymmetric,&
        dbcsr_desymmetrize, dbcsr_type_symmetric, dbcsr_init_p, dbcsr_multiply, dbcsr_add,&
        dbcsr_get_info, dbcsr_type, dbcsr_deallocate_matrix
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply, &
                                             dbcsr_allocate_matrix_set, &
                                             dbcsr_deallocate_matrix_set   
   USE iterate_matrix,                  ONLY: matrix_sqrt_Newton_Schulz,&
                                              matrix_sqrt_proot                                             
   USE cp_fm_basic_linalg,              ONLY: cp_fm_trace
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type, cp_fm_get_submatrix
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_result_methods,               ONLY: cp_results_erase,&
                                              put_results
   USE cp_result_types,                 ONLY: cp_result_type
   USE qs_dcdr_utils,                   ONLY: dcdr_env_cleanup,&
                                              dcdr_env_init            
   USE ec_env_types,                    ONLY: energy_correction_type
   USE ec_orth_solver,                  ONLY: ec_response_ao, ec_response_ao_vel     
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE machine,                         ONLY: m_flush
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: angstrom
   USE qs_2nd_kernel_ao,                ONLY: build_dm_response
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_linres_methods,               ONLY: linres_read_restart, linres_read_restart_ao, &
                                              linres_solver, &
                                              linres_write_restart, linres_write_restart_ao
   USE qs_linres_types,                 ONLY: get_polar_env,&
                                              linres_control_type,&
                                              polar_env_type,&
                                              set_polar_env, dcdr_env_type
   USE qs_matrix_pools,                 ONLY: qs_matrix_pools_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              get_neighbor_list_set_p,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_p_env_methods,                ONLY: p_env_create,&
                                              p_env_psi0_changed  
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type                                                                                          
   USE qs_p_env_types,                  ONLY: qs_p_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: polar_env_init, polar_polar, polar_print, polar_response, write_polarisability_tensor,&
               polar_response_ao, polar_polar_ao, polar_print_ao

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_polar_utils'

CONTAINS

! **************************************************************************************************
!> \brief Initialize the polar environment
!> \param qs_env ...
!> \par History
!>      06.2018 polar_env integrated into qs_env (MK)
! **************************************************************************************************
SUBROUTINE polar_env_init(qs_env)

   TYPE(qs_environment_type), POINTER                 :: qs_env

   CHARACTER(LEN=*), PARAMETER                        :: routineN = 'polar_env_init'

   INTEGER                                            :: handle, idir, iounit, ispin, m, nao, &
                                                         nmo, nspins
   LOGICAL                                            :: do_raman_vel
   TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
   TYPE(cp_fm_type), POINTER                          :: mo_coeff
   TYPE(cp_logger_type), POINTER                      :: logger
   TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
   TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
   TYPE(dft_control_type), POINTER                    :: dft_control
   TYPE(linres_control_type), POINTER                 :: linres_control
   TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
   TYPE(neighbor_list_set_p_type), DIMENSION(:), &
   POINTER                                            :: sab_orb
   TYPE(polar_env_type), POINTER                      :: polar_env
   TYPE(section_vals_type), POINTER                   :: lr_section, polar_section

   CALL timeset(routineN, handle)

   NULLIFY (dft_control)
   NULLIFY (linres_control)
   NULLIFY (logger)
   NULLIFY (matrix_s)
   NULLIFY (mos)
   NULLIFY (polar_env)
   NULLIFY (lr_section, polar_section, sab_orb)

   logger => cp_get_default_logger()
   lr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES")

   iounit = cp_print_key_unit_nr(logger, lr_section, "PRINT%PROGRAM_RUN_INFO", &
                                 extension=".linresLog")

   IF (iounit > 0) THEN
      WRITE (iounit, "(/,(T2,A))") "POLAR| Starting polarizability calculation", &
         "POLAR| Initialization of the polar environment"
   END IF

   polar_section => section_vals_get_subs_vals(qs_env%input, &
                                               "PROPERTIES%LINRES%POLAR")

   CALL get_qs_env(qs_env=qs_env, &
                   polar_env=polar_env, &
                   dft_control=dft_control, &
                   matrix_s_kp=matrix_s, &
                   matrix_ks=matrix_ks, &
                   linres_control=linres_control, &
                   mos=mos)
   ! Create polar environment if needed
   IF (.NOT. ASSOCIATED(polar_env)) THEN
      ALLOCATE (polar_env)
      CALL set_qs_env(qs_env=qs_env, polar_env=polar_env)
   END IF
   CALL get_polar_env(polar_env=polar_env, do_raman_vel=do_raman_vel)

   nspins = dft_control%nspins

   CALL section_vals_val_get(polar_section, "DO_RAMAN", l_val=polar_env%do_raman)
   CALL section_vals_val_get(polar_section, "PERIODIC_DIPOLE_OPERATOR", l_val=polar_env%do_periodic)
   CALL section_vals_val_get(polar_section, "AO_SOLVER", l_val=polar_env%ao_solver)
   CALL section_vals_val_get(polar_section, "DO_RAMAN_VEL", l_val=polar_env%do_raman_vel)
   CALL section_vals_val_get(polar_section, "DO_RAMAN_WANNIER", l_val=polar_env%do_raman_wannier)
   ! CALL section_vals_val_get(polar_section, "PRECONDITIONER", i_val=polar_env%preconditioner)
   ! CALL section_vals_val_get(polar_section, "PRECONDITIONER", i_val=linres_control%preconditioner_type)

   ! Allocate components of the polar environment if needed
   IF (.NOT. ASSOCIATED(polar_env%polar)) THEN
      ALLOCATE (polar_env%polar(3, 3))
      polar_env%polar(:, :) = 0.0_dp
   END IF
   IF (.NOT. ASSOCIATED(polar_env%dBerry_psi0)) THEN
      ALLOCATE (polar_env%dBerry_psi0(3, nspins))
   ELSE
      ! Remove previous matrices
      DO ispin = 1, nspins
         DO idir = 1, 3
            CALL cp_fm_release(polar_env%dBerry_psi0(idir, ispin))
         END DO
      END DO
   END IF
   IF (.NOT. ASSOCIATED(polar_env%psi1_dBerry)) THEN
      ALLOCATE (polar_env%psi1_dBerry(3, nspins))
   ELSE
      ! Remove previous matrices
      DO ispin = 1, nspins
         DO idir = 1, 3
            CALL cp_fm_release(polar_env%psi1_dBerry(idir, ispin))
         END DO
      END DO
   END IF
   IF (.NOT. ASSOCIATED(polar_env%psi1_dBerry_v)) THEN
      ALLOCATE (polar_env%psi1_dBerry_v(3, nspins))
   ELSE
      ! Remove previous matrices
      DO ispin = 1, nspins
         DO idir = 1, 3
            CALL cp_fm_release(polar_env%psi1_dBerry_v(idir, ispin))
         END DO
      END DO
   END IF
   DO ispin = 1, nspins
      CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff, nao=nao, nmo=nmo)
      CALL cp_fm_get_info(mo_coeff, ncol_global=m, nrow_global=nao)
      NULLIFY (tmp_fm_struct)
      CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nao, &
                               ncol_global=m, &
                               context=mo_coeff%matrix_struct%context)
      DO idir = 1, 3
         CALL cp_fm_create(polar_env%dBerry_psi0(idir, ispin), tmp_fm_struct)
         CALL cp_fm_create(polar_env%psi1_dBerry(idir, ispin), tmp_fm_struct)
         CALL cp_fm_create(polar_env%psi1_dBerry_v(idir, ispin), tmp_fm_struct)
      END DO
      CALL cp_fm_struct_release(tmp_fm_struct)
   END DO

   IF (ASSOCIATED(polar_env%dBerry_psi0_ao)) CALL dbcsr_deallocate_matrix_set(polar_env%dBerry_psi0_ao)
   IF (ASSOCIATED(polar_env%psi1_dBerry_ao)) CALL dbcsr_deallocate_matrix_set(polar_env%psi1_dBerry_ao)

   CALL dbcsr_allocate_matrix_set(polar_env%psi1_dBerry_ao, 3, nspins)
   CALL dbcsr_allocate_matrix_set(polar_env%dBerry_psi0_ao, 3, nspins)
   DO idir = 1, 3
      DO ispin = 1, nspins
         CALL dbcsr_init_p(polar_env%psi1_dBerry_ao(idir, ispin)%matrix)
         CALL dbcsr_create(polar_env%psi1_dBerry_ao(idir, ispin)%matrix, &
                           template=matrix_ks(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(polar_env%psi1_dBerry_ao(idir, ispin)%matrix, matrix_s(1, 1)%matrix)
         CALL dbcsr_set(polar_env%psi1_dBerry_ao(idir, ispin)%matrix, 0.0)

         CALL dbcsr_init_p(polar_env%dBerry_psi0_ao(idir, ispin)%matrix)
         CALL dbcsr_create(polar_env%dBerry_psi0_ao(idir, ispin)%matrix, &
                           template=matrix_ks(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(polar_env%dBerry_psi0_ao(idir, ispin)%matrix, matrix_s(1, 1)%matrix)
         CALL dbcsr_set(polar_env%dBerry_psi0_ao(idir, ispin)%matrix, 0.0)
      END DO
   END DO

   CALL cp_print_key_finished_output(iounit, logger, lr_section, &
                                     "PRINT%PROGRAM_RUN_INFO")

   CALL timestop(handle)

END SUBROUTINE polar_env_init

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \par History
!>      06.2018 polar_env integrated into qs_env (MK)
! **************************************************************************************************
SUBROUTINE polar_polar(qs_env)

   TYPE(qs_environment_type), POINTER                 :: qs_env

   CHARACTER(LEN=*), PARAMETER                        :: routineN = 'polar_polar'
   INTEGER, DIMENSION(2)                              :: nbr_center = -1
   INTEGER                                            :: handle, i, iounit, ispin, nspins, z, &
                                                         nmo, nao, icenter, j, max_center
   LOGICAL                                            :: do_periodic, do_raman, run_stopped, &
                                                         do_raman_vel, do_raman_wannier
   REAL(dp)                                           :: ptmp
   REAL(dp), DIMENSION(:, :), POINTER                 :: polar, polar_tmp
   REAL(dp), ALLOCATABLE,DIMENSION(:,:,:,: )          :: alpha_w
   REAL(dp), DIMENSION(:, :), POINTER                 :: tmp_psi1, tmp_psi0
   TYPE(cell_type), POINTER                           :: cell
   TYPE(cp_fm_type), POINTER                          :: mo_coeff
   TYPE(cp_fm_type), DIMENSION(:, :), POINTER         :: dBerry_psi0, psi1_dBerry, psi1_dBerry_v
   TYPE(cp_logger_type), POINTER                      :: logger
   TYPE(dbcsr_p_type), DIMENSION(:), &
   POINTER                                            :: matrix_s
   ! TYPE(dcdr_env_type)                                :: dcdr_env
   TYPE(dft_control_type), POINTER                    :: dft_control
   TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
   TYPE(polar_env_type), POINTER                      :: polar_env

   CALL timeset(routineN, handle)

   NULLIFY (cell, dft_control, polar, psi1_dBerry, psi1_dBerry_v, logger)
   NULLIFY (mos, dBerry_psi0)
   logger => cp_get_default_logger()
   iounit = cp_logger_get_default_io_unit(logger)

   CALL get_qs_env(qs_env=qs_env, &
                   cell=cell, &
                   dft_control=dft_control, &
                   mos=mos, &
                   matrix_s = matrix_s, &
                   polar_env=polar_env)

   nspins = dft_control%nspins

   CALL get_polar_env(polar_env=polar_env, &
                      do_raman=do_raman, &
                      do_raman_vel=do_raman_vel, &
                      do_raman_wannier=do_raman_wannier, &
                      nbr_center=nbr_center, &
                      run_stopped=run_stopped)

   IF (.NOT. run_stopped .AND. do_raman  .AND. .NOT. do_raman_wannier) THEN

      CALL cite_reference(Luber2014)

      CALL get_polar_env(polar_env=polar_env, &
                         do_periodic=do_periodic, &
                         dBerry_psi0=dBerry_psi0, &
                         polar=polar, &
                         psi1_dBerry=psi1_dBerry)

      ! Initialize
      ALLOCATE (polar_tmp(3, 3))
      polar_tmp(:, :) = 0.0_dp

      DO i = 1, 3 ! directions of electric field
         DO z = 1, 3 !dipole directions
            DO ispin = 1, dft_control%nspins
               !SL compute trace
               ptmp = 0.0_dp
               CALL cp_fm_trace(psi1_dBerry(i, ispin), dBerry_psi0(z, ispin), ptmp)
               polar_tmp(i, z) = polar_tmp(i, z) - 2.0_dp*ptmp
            END DO
         END DO
      END DO !spin

      IF (do_periodic  .AND. .NOT. do_raman_wannier) THEN
         polar(:, :) = MATMUL(MATMUL(cell%hmat, polar_tmp), TRANSPOSE(cell%hmat))/(twopi*twopi)
      ELSE
         polar(:, :) = polar_tmp(:, :)
      END IF 
      !SL evtl maxocc instead?
      IF (dft_control%nspins == 1) THEN
         polar(:, :) = 2.0_dp*polar(:, :)
      END IF

      IF (ASSOCIATED(polar_tmp)) THEN
         DEALLOCATE (polar_tmp)
      END IF

   END IF ! do_raman
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   IF (.NOT. run_stopped .AND. do_raman .AND. do_raman_wannier ) THEN

      CALL cite_reference(Luber2014)

      CALL get_polar_env(polar_env=polar_env, &
                         do_periodic=do_periodic, &
                         dBerry_psi0=dBerry_psi0, &
                         polar=polar, &
                         psi1_dBerry=psi1_dBerry)

      ! Initialize
      ALLOCATE (polar_tmp(3, 3))
      polar_tmp(:, :) = 0.0_dp

      DO ispin=1, nspins
         CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff)
         CALL cp_fm_get_info(mo_coeff, ncol_global=nmo, nrow_global=nao)
      END DO

      max_center = polar_env%nbr_center(1)
      ALLOCATE(alpha_w(3,3,max_center,nspins))
      ALLOCATE(tmp_psi1(nao,1))
      ALLOCATE(tmp_psi0(nao,1))
      alpha_w= 0.0_dp

      DO ispin = 1,nspins
         DO i=1,3
            DO j=1,3
               DO icenter = 1,polar_env%nbr_center(1)

                  !  CALL cp_fm_get_submatrix(psi1_dBerry(i, ispin),tmp_psi1,1,icenter,nao,1)
                  !  CALL cp_fm_get_submatrix(dBerry_psi0(j, ispin),tmp_psi0,1,icenter,nao,1)
                  ! ! alpha_w is -(1/4) final polarizability 
                  ! alpha_w(i,j,icenter,ispin) = DOT_PRODUCT(tmp_psi1(:,1),tmp_psi0(:,1))
                  ! polar(i,j)= polar(i,j) - 4.0_dp*alpha_w(i,j,icenter,ispin)
               END DO
               ! print *, polar(i,j)
            END DO
          END DO
      END DO !spin

      ptmp = 0.0_dp

      DO i = 1, 3 ! directions of electric field
         DO z = 1, 3 !dipole directions
            DO ispin = 1, dft_control%nspins
               !SL compute trace
               CALL cp_fm_trace(psi1_dBerry(i, ispin), dBerry_psi0(z, ispin), ptmp)
               polar_tmp(i, z) = polar_tmp(i, z) - 2.0_dp*ptmp
            END DO
         END DO
      END DO !spin
       
      polar(:, :) = polar_tmp(:, :)

      IF (dft_control%nspins == 1) THEN
         polar(:, :) = 2.0_dp*polar(:, :)
      END IF

      IF (ASSOCIATED(polar_tmp)) THEN
         DEALLOCATE (polar_tmp)
         DEALLOCATE(alpha_w)
         DEALLOCATE(tmp_psi0)
         DEALLOCATE(tmp_psi1)
      END IF

   END IF ! do_raman_wannier
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   IF (.NOT. run_stopped .AND. do_raman_vel) THEN

      CALL cite_reference(Luber2014)

      CALL get_polar_env(polar_env=polar_env, &
                         do_periodic=do_periodic, &
                         dBerry_psi0=dBerry_psi0, &
                         polar=polar, &
                         psi1_dBerry=psi1_dBerry, &
                         psi1_dBerry_v=psi1_dBerry_v)

      ! Initialize
      ALLOCATE (polar_tmp(3, 3))
      polar_tmp(:, :) = 0.0_dp

      DO i = 1, 3 ! directions of electric field
         DO z = 1, 3 !dipole directions
            DO ispin = 1, dft_control%nspins
               !SL compute trace
               ptmp = 0.0_dp
               CALL cp_fm_trace(psi1_dBerry_v(i, ispin), dBerry_psi0(z, ispin), ptmp)
               polar_tmp(i, z) = polar_tmp(i, z) - 2.0_dp*ptmp
            END DO
         END DO
      END DO !spin

      polar(:, :) = polar_tmp(:, :)
      !SL evtl maxocc instead?
      IF (dft_control%nspins == 1) THEN
         polar(:, :) = 2.0_dp*polar(:, :)
      END IF

      IF (ASSOCIATED(polar_tmp)) THEN
         DEALLOCATE (polar_tmp) 
      END IF

   END IF ! do_raman_vel
   CALL timestop(handle)

END SUBROUTINE polar_polar

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \par History
!>      01.2023 by Ravi,
! **************************************************************************************************
   SUBROUTINE polar_polar_ao(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'polar_polar_ao'

      INTEGER                                            :: handle, i, iounit, ispin, nspins, z
      LOGICAL                                            :: do_periodic, do_raman, do_raman_vel, run_stopped,&
                                                            do_raman_wannier
      REAL(dp)                                           :: ptmp_ao
      REAL(dp), DIMENSION(:, :), POINTER                 :: polar, polar_tmp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type), DIMENSION(:, :), POINTER         :: dBerry_psi0, psi1_dBerry
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: dBerry_psi0_ao, psi1_dBerry_ao
      TYPE(dbcsr_type), POINTER                          :: temp_dbcsr
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(polar_env_type), POINTER                      :: polar_env

      CALL timeset(routineN, handle)

      NULLIFY (cell, dft_control, polar, psi1_dBerry, logger)
      NULLIFY (mos, dBerry_psi0)
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      mos=mos, &
                      polar_env=polar_env)

      nspins = dft_control%nspins

      CALL get_polar_env(polar_env=polar_env, &
                         do_raman=do_raman, &
                         do_raman_vel=do_raman_vel, &
                         run_stopped=run_stopped)

      IF (.NOT. run_stopped .AND. qs_env%polar_env%ao_solver .AND. .NOT. qs_env%polar_env%do_raman_vel) THEN

         CALL cite_reference(Luber2014)

         CALL get_polar_env(polar_env=polar_env, &
                            do_periodic=do_periodic, &
                            do_raman_wannier=do_raman_wannier,&
                            dBerry_psi0=dBerry_psi0, &
                            dBerry_psi0_ao=dBerry_psi0_ao, &
                            polar=polar, &
                            psi1_dBerry=psi1_dBerry, &
                            psi1_dBerry_ao=psi1_dBerry_ao)

         NULLIFY (temp_dbcsr)
         ALLOCATE (temp_dbcsr)
         CALL dbcsr_copy(temp_dbcsr, psi1_dBerry_ao(1, 1)%matrix)

         ! Initialize
         ALLOCATE (polar_tmp(3, 3))
         polar_tmp(:, :) = 0.0_dp

         DO i = 1, 3 ! directions of electric field
            DO z = 1, 3 !dipole directions
               DO ispin = 1, dft_control%nspins
                  !SL compute trace
                  ptmp_ao = 0.0_dp
                  CALL dbcsr_dot(psi1_dBerry_ao(i, ispin)%matrix, dBerry_psi0_ao(z, ispin)%matrix, ptmp_ao)
                  polar_tmp(i, z) = polar_tmp(i, z) - 1.0_dp*ptmp_ao
               END DO
            END DO
         END DO !spin

         IF (do_periodic  .AND. .NOT. do_raman_wannier) THEN
            polar(:, :) = MATMUL(MATMUL(cell%hmat, polar_tmp), TRANSPOSE(cell%hmat))/(twopi*twopi)
         ELSE
            polar(:, :) = polar_tmp(:, :)
         END IF 

         !SL evtl maxocc instead?
         IF (dft_control%nspins == 1) THEN
            polar(:, :) = 2.0_dp*polar(:, :)
         END IF

         IF (ASSOCIATED(polar_tmp)) THEN
            DEALLOCATE (polar_tmp)
         END IF
         CALL dbcsr_deallocate_matrix(temp_dbcsr)

      END IF ! do_raman 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF (.NOT. run_stopped .AND. .NOT. qs_env%polar_env%ao_solver .AND. .NOT. qs_env%polar_env%do_raman_vel &
         .AND. do_raman_wannier) THEN

         CALL cite_reference(Luber2014)

         CALL get_polar_env(polar_env=polar_env, &
                            do_periodic=do_periodic, &
                            do_raman_wannier=do_raman_wannier,&
                            dBerry_psi0=dBerry_psi0, &
                            dBerry_psi0_ao=dBerry_psi0_ao, &
                            polar=polar, &
                            psi1_dBerry=psi1_dBerry, &
                            psi1_dBerry_ao=psi1_dBerry_ao)

         ! Initialize
         ALLOCATE (polar_tmp(3, 3))
         polar_tmp(:, :) = 0.0_dp

         DO i = 1, 3 ! directions of electric field
            DO z = 1, 3 !dipole directions
               DO ispin = 1, dft_control%nspins
                  !SL compute trace
                  ptmp_ao = 0.0_dp
                  CALL dbcsr_dot(psi1_dBerry_ao(i, ispin)%matrix, dBerry_psi0_ao(z, ispin)%matrix, ptmp_ao)
                  polar_tmp(i, z) = polar_tmp(i, z) - 1.0_dp*ptmp_ao
               END DO
            END DO
         END DO !spin

         polar(:, :) = polar_tmp(:, :)

         !SL evtl maxocc instead?
         IF (dft_control%nspins == 1) THEN
            polar(:, :) = 2.0_dp*polar(:, :)
         END IF

         IF (ASSOCIATED(polar_tmp)) THEN
            DEALLOCATE (polar_tmp)
         END IF

      END IF ! do_raman 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IF (.NOT. run_stopped .AND. do_raman_vel) THEN

         CALL cite_reference(Luber2014)

         CALL get_polar_env(polar_env=polar_env, &
                            do_periodic=do_periodic, &
                            dBerry_psi0=dBerry_psi0, &
                            dBerry_psi0_ao=dBerry_psi0_ao, &
                            polar=polar, &
                            psi1_dBerry=psi1_dBerry, &
                            psi1_dBerry_ao=psi1_dBerry_ao)

         ! Initialize
         ALLOCATE (polar_tmp(3, 3))
         polar_tmp(:, :) = 0.0_dp

         DO i = 1, 3 ! directions of electric field
            DO z = 1, 3 !dipole directions
               DO ispin = 1, dft_control%nspins
                  ! compute trace
                  ptmp_ao = 0.0_dp
                  CALL dbcsr_dot(polar_env%psi1_dBerry_ao(z, ispin)%matrix, dBerry_psi0_ao(i, ispin)%matrix, ptmp_ao)
                  polar_tmp(i, z) = polar_tmp(i, z) - 1.0_dp*ptmp_ao
               END DO
            END DO
         END DO !spin

         polar(:, :) = polar_tmp(:, :)

         ! evtl maxocc instead?
         IF (dft_control%nspins == 1) THEN
            polar(:, :) = 2.0_dp*polar(:, :)
         END IF

         IF (ASSOCIATED(polar_tmp)) THEN
            DEALLOCATE (polar_tmp)
         END IF

      END IF        
      CALL timestop(handle)

   END SUBROUTINE polar_polar_ao

! **************************************************************************************************
!> \brief Print information related to the polarisability tensor
!> \param qs_env ...
!> \par History
!>      06.2018 polar_env integrated into qs_env (MK)
! **************************************************************************************************
   SUBROUTINE polar_print(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=default_string_length)               :: description
      INTEGER                                            :: iounit, unit_p
      LOGICAL                                            :: do_raman, do_raman_vel, run_stopped
      REAL(dp), DIMENSION(:, :), POINTER                 :: polar
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(polar_env_type), POINTER                      :: polar_env
      TYPE(section_vals_type), POINTER                   :: polar_section

      NULLIFY (logger, dft_control, para_env, results)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      polar_env=polar_env, &
                      results=results, &
                      para_env=para_env)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      polar_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES%POLAR")

      CALL get_polar_env(polar_env=polar_env, &
                         polar=polar, &
                         do_raman=do_raman, &
                         do_raman_vel=do_raman_vel, &
                         run_stopped=run_stopped)

      IF (.NOT. run_stopped .AND. do_raman .AND. .NOT. do_raman_vel) THEN

         description = "[POLAR]"
         CALL cp_results_erase(results, description=description)
         CALL put_results(results, description=description, values=polar(:, :))

         IF (BTEST(cp_print_key_should_output(logger%iter_info, polar_section, &
                                              "PRINT%POLAR_MATRIX"), cp_p_file)) THEN

            unit_p = cp_print_key_unit_nr(logger, polar_section, "PRINT%POLAR_MATRIX", &
                                          extension=".data", middle_name="raman", log_filename=.FALSE.)
            IF (unit_p > 0) THEN
               IF (unit_p /= iounit) THEN
                  WRITE (unit_p, *)
                  WRITE (unit_p, '(T10,A)') 'POLARIZABILITY TENSOR (atomic units):'
                  WRITE (unit_p, '(T10,A,3F15.5)') "xx,yy,zz", polar(1, 1), polar(2, 2), polar(3, 3)
                  WRITE (unit_p, '(T10,A,3F15.5)') "xy,xz,yz", polar(1, 2), polar(1, 3), polar(2, 3)
                  WRITE (unit_p, '(T10,A,3F15.5)') "yx,zx,zy", polar(2, 1), polar(3, 1), polar(3, 2)
                  WRITE (unit_p, '(T10,A)') 'POLARIZABILITY TENSOR (Angstrom^3):'
                  WRITE (unit_p, '(T10,A,3F15.5)') "xx,yy,zz", polar(1, 1)*angstrom**3, &
                     polar(2, 2)*angstrom**3, polar(3, 3)*angstrom**3
                  WRITE (unit_p, '(T10,A,3F15.5)') "xy,xz,yz", polar(1, 2)*angstrom**3, &
                     polar(1, 3)*angstrom**3, polar(2, 3)*angstrom**3
                  WRITE (unit_p, '(T10,A,3F15.5)') "yx,zx,zy", polar(2, 1)*angstrom**3, &
                     polar(3, 1)*angstrom**3, polar(3, 2)*angstrom**3
                  CALL cp_print_key_finished_output(unit_p, logger, polar_section, &
                                                    "PRINT%POLAR_MATRIX")
               END IF
            END IF
         END IF
         IF (iounit > 0) THEN
            WRITE (iounit, '(/,T2,A)') &
               'POLAR| Polarizability tensor [a.u.]'
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| xx,yy,zz', polar(1, 1), polar(2, 2), polar(3, 3)
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| xy,xz,yz', polar(1, 2), polar(1, 3), polar(2, 3)
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| yx,zx,zy', polar(2, 1), polar(3, 1), polar(3, 2)
            WRITE (iounit, '(/,T2,A)') &
               'POLAR| Polarizability tensor [ang^3]'
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| xx,yy,zz', polar(1, 1)*angstrom**3, polar(2, 2)*angstrom**3, polar(3, 3)*angstrom**3
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| xy,xz,yz', polar(1, 2)*angstrom**3, polar(1, 3)*angstrom**3, polar(2, 3)*angstrom**3
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| yx,zx,zy', polar(2, 1)*angstrom**3, polar(3, 1)*angstrom**3, polar(3, 2)*angstrom**3
         END IF
      END IF

      IF (.NOT. run_stopped .AND. do_raman .AND. do_raman_vel) THEN

         description = "[POLAR]"
         CALL cp_results_erase(results, description=description)
         CALL put_results(results, description=description, values=polar(:, :))

         IF (BTEST(cp_print_key_should_output(logger%iter_info, polar_section, &
                                              "PRINT%POLAR_MATRIX"), cp_p_file)) THEN

            unit_p = cp_print_key_unit_nr(logger, polar_section, "PRINT%POLAR_MATRIX", &
                                          extension=".data", middle_name="raman", log_filename=.FALSE.)
            IF (unit_p > 0) THEN
               IF (unit_p /= iounit) THEN
                  WRITE (unit_p, *)
                  WRITE (unit_p, '(T10,A)') 'POLARIZABILITY TENSOR (atomic units):'
                  WRITE (unit_p, '(T10,A,3F15.5)') "xx,yy,zz", polar(1, 1), polar(2, 2), polar(3, 3)
                  WRITE (unit_p, '(T10,A,3F15.5)') "xy,xz,yz", polar(1, 2), polar(1, 3), polar(2, 3)
                  WRITE (unit_p, '(T10,A,3F15.5)') "yx,zx,zy", polar(2, 1), polar(3, 1), polar(3, 2)
                  WRITE (unit_p, '(T10,A)') 'POLARIZABILITY TENSOR (Angstrom^3):'
                  WRITE (unit_p, '(T10,A,3F15.5)') "xx,yy,zz", polar(1, 1)*angstrom**3, &
                     polar(2, 2)*angstrom**3, polar(3, 3)*angstrom**3
                  WRITE (unit_p, '(T10,A,3F15.5)') "xy,xz,yz", polar(1, 2)*angstrom**3, &
                     polar(1, 3)*angstrom**3, polar(2, 3)*angstrom**3
                  WRITE (unit_p, '(T10,A,3F15.5)') "yx,zx,zy", polar(2, 1)*angstrom**3, &
                     polar(3, 1)*angstrom**3, polar(3, 2)*angstrom**3
                  CALL cp_print_key_finished_output(unit_p, logger, polar_section, &
                                                    "PRINT%POLAR_MATRIX")
               END IF
            END IF
         END IF
         IF (iounit > 0) THEN
            WRITE (iounit, '(/,T2,A)') &
               'POLAR| Polarizability tensor [a.u.]'
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| xx,yy,zz', polar(1, 1), polar(2, 2), polar(3, 3)
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| xy,xz,yz', polar(1, 2), polar(1, 3), polar(2, 3)
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| yx,zx,zy', polar(2, 1), polar(3, 1), polar(3, 2)
            WRITE (iounit, '(/,T2,A)') &
               'POLAR| Polarizability tensor [ang^3]'
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| xx,yy,zz', polar(1, 1)*angstrom**3, polar(2, 2)*angstrom**3, polar(3, 3)*angstrom**3
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| xy,xz,yz', polar(1, 2)*angstrom**3, polar(1, 3)*angstrom**3, polar(2, 3)*angstrom**3
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| yx,zx,zy', polar(2, 1)*angstrom**3, polar(3, 1)*angstrom**3, polar(3, 2)*angstrom**3
         END IF
      END IF  !do_raman_vel

   END SUBROUTINE polar_print


! **************************************************************************************************
!> \brief Print information related to the polarisability tensor
!> \param qs_env ...
!> \par History
!>      01.2023, Ravi
! **************************************************************************************************
   SUBROUTINE polar_print_ao(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=default_string_length)               :: description
      INTEGER                                            :: iounit, unit_p
      LOGICAL                                            :: do_raman, do_raman_vel, run_stopped
      REAL(dp), DIMENSION(:, :), POINTER                 :: polar
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(polar_env_type), POINTER                      :: polar_env
      TYPE(section_vals_type), POINTER                   :: polar_section

      NULLIFY (logger, dft_control, para_env, results)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      polar_env=polar_env, &
                      results=results, &
                      para_env=para_env)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      polar_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES%POLAR")

      CALL get_polar_env(polar_env=polar_env, &
                         polar=polar, &
                         do_raman=do_raman, &
                         do_raman_vel=do_raman_vel, &
                         run_stopped=run_stopped)

      IF (.NOT. run_stopped .AND. qs_env%polar_env%ao_solver .AND. .NOT. qs_env%polar_env%do_raman_vel) THEN

         description = "[POLAR]"
         CALL cp_results_erase(results, description=description)
         CALL put_results(results, description=description, values=polar(:, :))

         IF (BTEST(cp_print_key_should_output(logger%iter_info, polar_section, &
                                              "PRINT%POLAR_MATRIX"), cp_p_file)) THEN

            unit_p = cp_print_key_unit_nr(logger, polar_section, "PRINT%POLAR_MATRIX", &
                                          extension=".data", middle_name="raman", log_filename=.FALSE.)
            IF (unit_p > 0) THEN
               IF (unit_p /= iounit) THEN
                  WRITE (unit_p, *)
                  WRITE (unit_p, '(T10,A)') 'POLARIZABILITY TENSOR (atomic units):'
                  WRITE (unit_p, '(T10,A,3F15.5)') "xx,yy,zz", polar(1, 1), polar(2, 2), polar(3, 3)
                  WRITE (unit_p, '(T10,A,3F15.5)') "xy,xz,yz", polar(1, 2), polar(1, 3), polar(2, 3)
                  WRITE (unit_p, '(T10,A,3F15.5)') "yx,zx,zy", polar(2, 1), polar(3, 1), polar(3, 2)
                  WRITE (unit_p, '(T10,A)') 'POLARIZABILITY TENSOR (Angstrom^3):'
                  WRITE (unit_p, '(T10,A,3F15.5)') "xx,yy,zz", polar(1, 1)*angstrom**3, &
                     polar(2, 2)*angstrom**3, polar(3, 3)*angstrom**3
                  WRITE (unit_p, '(T10,A,3F15.5)') "xy,xz,yz", polar(1, 2)*angstrom**3, &
                     polar(1, 3)*angstrom**3, polar(2, 3)*angstrom**3
                  WRITE (unit_p, '(T10,A,3F15.5)') "yx,zx,zy", polar(2, 1)*angstrom**3, &
                     polar(3, 1)*angstrom**3, polar(3, 2)*angstrom**3
                  CALL cp_print_key_finished_output(unit_p, logger, polar_section, &
                                                    "PRINT%POLAR_MATRIX")
               END IF
            END IF
         END IF
         IF (iounit > 0) THEN
            WRITE (iounit, '(/,T2,A)') &
               'POLAR| Polarizability tensor [a.u.]'
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| xx,yy,zz', polar(1, 1), polar(2, 2), polar(3, 3)
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| xy,xz,yz', polar(1, 2), polar(1, 3), polar(2, 3)
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| yx,zx,zy', polar(2, 1), polar(3, 1), polar(3, 2)
            WRITE (iounit, '(/,T2,A)') &
               'POLAR| Polarizability tensor [ang^3]'
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| xx,yy,zz', polar(1, 1)*angstrom**3, polar(2, 2)*angstrom**3, polar(3, 3)*angstrom**3
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| xy,xz,yz', polar(1, 2)*angstrom**3, polar(1, 3)*angstrom**3, polar(2, 3)*angstrom**3
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| yx,zx,zy', polar(2, 1)*angstrom**3, polar(3, 1)*angstrom**3, polar(3, 2)*angstrom**3
         END IF
      END IF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IF (.NOT. run_stopped .AND. do_raman_vel) THEN
         description = "[POLAR]"
         CALL cp_results_erase(results, description=description)
         CALL put_results(results, description=description, values=polar(:, :))

         IF (BTEST(cp_print_key_should_output(logger%iter_info, polar_section, &
                                              "PRINT%POLAR_MATRIX"), cp_p_file)) THEN

            unit_p = cp_print_key_unit_nr(logger, polar_section, "PRINT%POLAR_MATRIX", &
                                          extension=".data", middle_name="raman", log_filename=.FALSE.)
            IF (unit_p > 0) THEN
               IF (unit_p /= iounit) THEN
                  WRITE (unit_p, *)
                  WRITE (unit_p, '(T10,A)') 'POLARIZABILITY TENSOR (atomic units):'
                  WRITE (unit_p, '(T10,A,3F15.5)') "xx,yy,zz", polar(1, 1), polar(2, 2), polar(3, 3)
                  WRITE (unit_p, '(T10,A,3F15.5)') "xy,xz,yz", polar(1, 2), polar(1, 3), polar(2, 3)
                  WRITE (unit_p, '(T10,A,3F15.5)') "yx,zx,zy", polar(2, 1), polar(3, 1), polar(3, 2)
                  WRITE (unit_p, '(T10,A)') 'POLARIZABILITY TENSOR (Angstrom^3):'
                  WRITE (unit_p, '(T10,A,3F15.5)') "xx,yy,zz", polar(1, 1)*angstrom**3, &
                     polar(2, 2)*angstrom**3, polar(3, 3)*angstrom**3
                  WRITE (unit_p, '(T10,A,3F15.5)') "xy,xz,yz", polar(1, 2)*angstrom**3, &
                     polar(1, 3)*angstrom**3, polar(2, 3)*angstrom**3
                  WRITE (unit_p, '(T10,A,3F15.5)') "yx,zx,zy", polar(2, 1)*angstrom**3, &
                     polar(3, 1)*angstrom**3, polar(3, 2)*angstrom**3
                  CALL cp_print_key_finished_output(unit_p, logger, polar_section, &
                                                    "PRINT%POLAR_MATRIX")
               END IF
            END IF
         END IF
         IF (iounit > 0) THEN
            WRITE (iounit, '(/,T2,A)') &
               'POLAR| Polarizability tensor [a.u.]'
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| xx,yy,zz', polar(1, 1), polar(2, 2), polar(3, 3)
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| xy,xz,yz', polar(1, 2), polar(1, 3), polar(2, 3)
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| yx,zx,zy', polar(2, 1), polar(3, 1), polar(3, 2)
            WRITE (iounit, '(/,T2,A)') &
               'POLAR| Polarizability tensor [ang^3]'
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| xx,yy,zz', polar(1, 1)*angstrom**3, polar(2, 2)*angstrom**3, polar(3, 3)*angstrom**3
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| xy,xz,yz', polar(1, 2)*angstrom**3, polar(1, 3)*angstrom**3, polar(2, 3)*angstrom**3
            WRITE (iounit, '(T2,A,T24,3(1X,F18.12))') &
               'POLAR| yx,zx,zy', polar(2, 1)*angstrom**3, polar(3, 1)*angstrom**3, polar(3, 2)*angstrom**3
         END IF
      END IF

   END SUBROUTINE polar_print_ao

! **************************************************************************************************
!> \brief Calculate the polarisability tensor using response theory
!> \param p_env ...
!> \param qs_env ...
!> \par History
!>      06.2018 polar_env integrated into qs_env (MK)
! **************************************************************************************************
   SUBROUTINE polar_response(p_env, qs_env)

      TYPE(qs_p_env_type)                                :: p_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'polar_response'

      INTEGER                                            :: handle, idir, iounit, ispin, nao, nmo, &
                                                            nspins                                                        
      LOGICAL                                            :: do_periodic, do_raman, do_raman_vel, should_stop,&
                                                            do_raman_wannier
      TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: h1_psi0, psi0_order, psi1
      TYPE(cp_fm_type), DIMENSION(:, :), POINTER         :: dBerry_psi0, psi1_dBerry, psi1_dBerry_v
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(polar_env_type), POINTER                      :: polar_env
      TYPE(qs_matrix_pools_type), POINTER                :: mpools
      TYPE(section_vals_type), POINTER                   :: lr_section, polar_section

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, linres_control, lr_section, polar_section)
      NULLIFY (logger, mpools, mo_coeff, para_env)
      NULLIFY (tmp_fm_struct, psi1_dBerry, dBerry_psi0, psi1_dBerry_v)

      logger => cp_get_default_logger()
      lr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES")
      polar_section => section_vals_get_subs_vals(qs_env%input, &
                                                  "PROPERTIES%LINRES%POLAR")

      iounit = cp_print_key_unit_nr(logger, lr_section, "PRINT%PROGRAM_RUN_INFO", &
                                    extension=".linresLog")
      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(T2,A,/)") &
            "POLAR| Self consistent optimization of the response wavefunctions"
      END IF

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      mpools=mpools, &
                      linres_control=linres_control, &
                      mos=mos, &
                      matrix_s=matrix_s, &
                      polar_env=polar_env, &
                      para_env=para_env)

      nspins = dft_control%nspins

      CALL get_polar_env(polar_env=polar_env, do_raman=do_raman, do_raman_vel=do_raman_vel, &
                           do_periodic=do_periodic)

      ! Allocate the vectors
      ALLOCATE (psi0_order(nspins))
      ALLOCATE (psi1(nspins), h1_psi0(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff)
         psi0_order(ispin) = mo_coeff
         CALL cp_fm_get_info(mo_coeff, ncol_global=nmo, nrow_global=nao)
         NULLIFY (tmp_fm_struct)
         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nao, &
                                  ncol_global=nmo, &
                                  context=mo_coeff%matrix_struct%context)
         CALL cp_fm_create(psi1(ispin), tmp_fm_struct)
         CALL cp_fm_create(h1_psi0(ispin), tmp_fm_struct)
         CALL cp_fm_struct_release(tmp_fm_struct)
      END DO

      IF (do_raman .AND. .NOT. do_raman_vel) THEN
         CALL get_polar_env(polar_env=polar_env, &
                            psi1_dBerry=psi1_dBerry, &
                            dBerry_psi0=dBerry_psi0, &
                            do_raman_wannier=do_raman_wannier))
         ! Restart
         IF (linres_control%linres_restart) THEN
            DO idir = 1, 3
               CALL linres_read_restart(qs_env, lr_section, psi1_dBerry(idir, :), idir, "psi1_dBerry")
            END DO
         ELSE 
            DO idir = 1, 3
               DO ispin = 1, nspins
                  CALL cp_fm_set_all(psi1_dBerry(idir, ispin), 0.0_dp)
               END DO
            END DO
         ENDIF

!!!!!!        
!!$OMP PARALLEL DEFAULT(NONE) &
!!$OMP SHARED (nspins, do_periodic, do_raman_wannier, iounit) &
!!$OMP SHARED (psi1_dBerry, dBerry_psi0, linres_control) &
!!$OMP SHARED (psi1, h1_psi0) &
!!$OMP SHARED (p_env, qs_env, psi0_order) &
!!$OMP PRIVATE (idir, ispin, should_stop, lr_section)

!!$OMP DO
         loop_idir: DO idir = 1, 3
            IF (iounit > 0) THEN
               IF (do_periodic .AND. .NOT. do_raman_wannier) THEN
                  WRITE (iounit, "(/,T2,A)") &
                     "POLAR| Response to the perturbation operator Berry phase_"//ACHAR(idir + 119)
               ELSE IF (do_periodic .AND.  do_raman_wannier) THEN
                  WRITE (iounit, "(/,T2,A)") &
                     "POLAR| Response to the perturbation operator MLWF_"//ACHAR(idir + 119)
                  ELSE 
                     WRITE (iounit, "(/,T2,A)") &
                        "POLAR| Response to the perturbation operator R_"//ACHAR(idir + 119)                     
               END IF
            END IF
            ! Do scf cycle to optimize psi1
!!$ single
!!$ master            
            DO ispin = 1, nspins
               CALL cp_fm_to_fm(psi1_dBerry(idir, ispin), psi1(ispin))
               CALL cp_fm_to_fm(dBerry_psi0(idir, ispin), h1_psi0(ispin))
            END DO
            !
            linres_control%lr_triplet = .FALSE. ! we do singlet response
            linres_control%do_kernel = .TRUE. ! we do coupled response
            linres_control%converged = .FALSE.

            CALL linres_solver(p_env, qs_env, psi1, h1_psi0, psi0_order, iounit, should_stop)

            ! Copy the response
            DO ispin = 1, nspins
               CALL cp_fm_to_fm(psi1(ispin), psi1_dBerry(idir, ispin))
            END DO
            !
            ! Write the new result to the restart file
            IF (linres_control%linres_restart) THEN
               CALL linres_write_restart(qs_env, lr_section, psi1_dBerry(idir, :), idir, "psi1_dBerry")
            END IF
         END DO loop_idir
!!$OMP END DO
!!$OMP END PARALLEL 
!!!!
      END IF ! do_raman

      CALL set_polar_env(polar_env, run_stopped=should_stop)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF (do_raman_vel) THEN

         CALL get_polar_env(polar_env=polar_env, &
                            psi1_dBerry=psi1_dBerry, &
                            psi1_dBerry_v=psi1_dBerry_v, &
                            dBerry_psi0=dBerry_psi0)
         ! Restart
         IF (linres_control%linres_restart) THEN
            DO idir = 1, 3
               ! psi1_ptr = psi1_dBerry(idir, :)
               CALL linres_read_restart(qs_env, lr_section, psi1_dBerry(idir, :), idir, "psi1_dBerry")
            END DO
         ELSE
            DO idir = 1, 3
               DO ispin = 1, nspins
                  CALL cp_fm_set_all(psi1_dBerry(idir, ispin), 0.0_dp)
               END DO
            END DO
         END IF
         loop_idir_vel: DO idir = 1, 3
            IF (iounit > 0) THEN
               WRITE (iounit, "(/,T2,A)") &
                  "POLAR| First uncoupled response to the velocity dipole perturbation operator- "//ACHAR(idir + 119)
            END IF

            ! Do scf cycle to optimize psi1
            DO ispin = 1, nspins
               CALL cp_fm_to_fm(psi1_dBerry(idir, ispin), psi1(ispin))
               CALL cp_fm_to_fm(dBerry_psi0(idir, ispin), h1_psi0(ispin))
            END DO
            !
            linres_control%lr_triplet = .FALSE. ! we do singlet response
            linres_control%do_kernel = .FALSE. ! we do uncoupled response
            linres_control%converged = .FALSE.

            CALL linres_solver(p_env, qs_env, psi1, h1_psi0, psi0_order, iounit, should_stop)

            Do ispin=1, nspins
               CALL cp_fm_to_fm(psi1(ispin), psi1_dBerry(idir, ispin))
            END DO
         END DO loop_idir_vel

         loop_idir_vel2: DO idir = 1, 3
            !Initial guess for psi1
              DO ispin = 1,nspins
               CALL cp_fm_set_all(psi1(ispin),0.0_dp)
               CALL cp_fm_set_all(h1_psi0(ispin),0.0_dp)
               CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, psi1_dBerry(idir,ispin),&
                                             psi1(ispin), ncol=nmo, &
                                             alpha=-1.0_dp)
               CALL cp_fm_to_fm(psi1(ispin),h1_psi0(ispin))
               CALL cp_fm_to_fm(psi1(ispin),dBerry_psi0(idir, ispin))
            ENDDO

          ! Restart
            IF (linres_control%linres_restart) THEN
               ! psi1_ptr = psi1_dBerry_v(idir, :)
               CALL linres_read_restart(qs_env, lr_section, psi1_dBerry_v(idir, :), idir, "psi1_dBerry_v")
            ELSE
               DO ispin = 1, nspins
                  CALL cp_fm_set_all(psi1_dBerry_v(idir, ispin), 0.0_dp)
               END DO
            END IF
            IF (iounit > 0) THEN
               WRITE (iounit, "(/,T2,A)") &
                  "POLAR| Second coupled response to the velocity dipole perturbation operator_ "//ACHAR(idir + 119)
            END IF

            ! Do scf cycle to optimize psi1
            DO ispin = 1, nspins
               CALL cp_fm_to_fm(psi1_dBerry_v(idir, ispin), psi1(ispin))
            END DO

            linres_control%lr_triplet = .FALSE. ! we do singlet response
            linres_control%do_kernel = .TRUE. ! we do coupled response
            linres_control%converged = .FALSE.

            ! print *, '', psi1_dBerry_v(idir, 1)%matrix%local_data(:,:)

            CALL linres_solver(p_env, qs_env, psi1, h1_psi0, psi0_order, iounit, should_stop)

            ! Copy the response
            DO ispin = 1, nspins
               CALL cp_fm_to_fm(psi1(ispin), psi1_dBerry_v(idir, ispin))
            END DO
            !
            ! Write the new result to the restart file
            IF (linres_control%linres_restart) THEN
               ! psi1 => psi1_dBerry_v(idir, :)
               CALL linres_write_restart(qs_env, lr_section, psi1_dBerry_v(idir, :), idir, "psi1_dBerry_v")
            END IF
         END DO loop_idir_vel2
      END IF ! do_raman_vel

      CALL set_polar_env(polar_env, run_stopped=should_stop)

      ! Clean up
      DO ispin = 1, nspins
         CALL cp_fm_release(psi1(ispin))
         CALL cp_fm_release(h1_psi0(ispin))
         ! DEALLOCATE (psi1(ispin), h1_psi0(ispin))
      END DO
      DEALLOCATE (psi1, h1_psi0, psi0_order)

      CALL cp_print_key_finished_output(iounit, logger, lr_section, &
                                        "PRINT%PROGRAM_RUN_INFO")

      CALL timestop(handle)

   END SUBROUTINE polar_response

! **************************************************************************************************
!> \brief Calculate the polarisability tensor using AO-based response theory
!> \param qs_env ...
!> \par History
!>      01.2023 by Ravi,   polar_env integrated into qs_env (MK)
! **************************************************************************************************
   SUBROUTINE polar_response_ao(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'polar_response_ao'

      INTEGER                                            :: handle, idir, iounit, ispin, &
                                                            nspins
      LOGICAL                                            :: do_periodic, do_raman_vel, should_stop, &
                                                            do_raman_wannier
      REAL(KIND=dp)                                      :: eps_filter                                                            
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: h1_matrix_p, matrix_hz, matrix_ks, &
                                                            matrix_p1, matrix_s, &
                                                            matrix_wz, rho_ao, h1_matrix_p2, matrix_p2
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: dBerry_psi0_ao, psi1_dBerry_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sab_all
      TYPE(polar_env_type), POINTER                      :: polar_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_matrix_pools_type), POINTER                :: mpools
      TYPE(qs_p_env_type), POINTER                       :: p_env_point
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: lr_section, polar_section

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, linres_control, lr_section, polar_section)
      NULLIFY (logger, mpools, mo_coeff, para_env)
      NULLIFY (matrix_p1, matrix_p2, h1_matrix_p, matrix_s, rho, matrix_ks)

      logger => cp_get_default_logger()
      lr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES")
      polar_section => section_vals_get_subs_vals(qs_env%input, &
                                                  "PROPERTIES%LINRES%POLAR")

      iounit = cp_print_key_unit_nr(logger, lr_section, "PRINT%PROGRAM_RUN_INFO", &
                                    extension=".linresLog")
      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(T2,A,/)") &
            "POLAR| Self consistent optimization of the response wavefunctions"
      END IF

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      mpools=mpools, &
                      linres_control=linres_control, &
                      matrix_s=matrix_s, &
                      polar_env=polar_env, &
                      para_env=para_env, &
                      rho=rho, &
                      sab_orb=sab_orb, &
                      sab_all=sab_all, &
                      ec_env=ec_env, &
                      matrix_ks=matrix_ks)

      nspins = dft_control%nspins

      CALL get_polar_env(polar_env=polar_env, &
                           do_raman_vel=do_raman_vel, do_periodic=do_periodic,do_raman_wannier=do_raman_wannier)

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      ! Allocate the vectors
      ALLOCATE (matrix_p1(nspins), h1_matrix_p(nspins))

      DO ispin = 1, nspins
         ALLOCATE (matrix_p1(ispin)%matrix, h1_matrix_p(ispin)%matrix)
         CALL dbcsr_create(h1_matrix_p(ispin)%matrix, name="h1. density mat", &
                           template=matrix_ks(ispin)%matrix, &
                           matrix_type=dbcsr_type_symmetric)

         CALL dbcsr_create(matrix_p1(ispin)%matrix, name="prtb density mat", &
                           template=matrix_ks(ispin)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)                           
      END DO

      NULLIFY (matrix_hz, matrix_wz)

      CALL dbcsr_allocate_matrix_set(matrix_hz, nspins)
      CALL dbcsr_allocate_matrix_set(matrix_wz, nspins)
      DO ispin = 1, nspins
         ALLOCATE (matrix_wz(ispin)%matrix)
         ALLOCATE (matrix_hz(ispin)%matrix)
         CALL dbcsr_create(matrix_wz(ispin)%matrix, name="Wz MATRIX", &
                           template=matrix_ks(1)%matrix)
         CALL dbcsr_copy(matrix_wz(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(matrix_wz(ispin)%matrix, 0.0_dp)

         CALL dbcsr_create(matrix_hz(ispin)%matrix, template=matrix_ks(1)%matrix)
         CALL dbcsr_copy(matrix_hz(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(matrix_hz(ispin)%matrix, 0.0_dp)
      END DO

      IF (qs_env%polar_env%ao_solver .AND. .NOT. qs_env%polar_env%do_raman_vel) THEN
         CALL get_polar_env(polar_env=polar_env, &
                            psi1_dBerry_ao=psi1_dBerry_ao, &
                            dBerry_psi0_ao=dBerry_psi0_ao)
         ! Restart

         DO idir = 1, 3
            DO ispin = 1, nspins
               CALL dbcsr_set(psi1_dBerry_ao(idir, ispin)%matrix, 0.0_dp)
            END DO
         END DO

         loop_idir: DO idir = 1, 3
            IF (iounit > 0) THEN
               IF (do_periodic .AND. do_raman_wannier) THEN
                  WRITE (iounit, "(/,T2,A)") &
                     "POLAR| Response to the perturbation operator MLWF_"//ACHAR(idir + 119)
               ELSE
                  WRITE (iounit, "(/,T2,A)") &
                     "POLAR| Response to the perturbation operator R_"//ACHAR(idir + 119)
               END IF
            END IF

            ! Do scf cycle to optimize psi1
            DO ispin = 1, nspins
               CALL dbcsr_copy(matrix_p1(ispin)%matrix, psi1_dBerry_ao(idir, ispin)%matrix)
               CALL dbcsr_copy(h1_matrix_p(ispin)%matrix, dBerry_psi0_ao(idir, ispin)%matrix)

            END DO

            linres_control%lr_triplet = .FALSE. ! we do singlet response
            linres_control%do_kernel = .TRUE. ! we do coupled response
            linres_control%converged = .FALSE.

            should_stop = .FALSE.

            ! initialize p_env
            ! Remark: mos environment is needed for this
            IF (ASSOCIATED(ec_env%p_env)) THEN
               CALL p_env_release(ec_env%p_env)
               DEALLOCATE (ec_env%p_env)
               NULLIFY (ec_env%p_env)
            END IF
            ALLOCATE (ec_env%p_env)
            CALL p_env_create(ec_env%p_env, qs_env, orthogonal_orbitals=.TRUE., &
                              linres_control=linres_control)
            CALL set_qs_env(qs_env, linres_control=linres_control)
            CALL p_env_psi0_changed(ec_env%p_env, qs_env)
            ! Total energy overwritten, replace with Etot from energy correction
            CALL get_qs_env(qs_env, energy=energy)
            energy%total = ec_env%etotal
            !
            p_env_point => ec_env%p_env

            CALL ec_response_ao(qs_env=qs_env, &
                                p_env=p_env_point, &
                                matrix_hz=h1_matrix_p, &
                                matrix_pz=matrix_p1, &
                                matrix_wz=matrix_wz, &
                                iounit=iounit, &
                                should_stop=should_stop)

            ! Copy the response
            DO ispin = 1, nspins
               CALL dbcsr_copy(psi1_dBerry_ao(idir, ispin)%matrix, matrix_p1(ispin)%matrix)
               CALL dbcsr_scale(psi1_dBerry_ao(idir, ispin)%matrix, 0.5_dp)
            END DO

         END DO loop_idir
      END IF ! do_raman

      DO ispin = 1, nspins
         CALL dbcsr_release(matrix_p1(ispin)%matrix)
         DEALLOCATE (matrix_p1(ispin)%matrix)
      END DO
      DEALLOCATE (matrix_p1)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      IF (do_raman_vel) THEN

         ALLOCATE (matrix_p1(nspins))
         ALLOCATE (matrix_p2(nspins))
         ALLOCATE (h1_matrix_p2(nspins))
         DO ispin = 1, nspins

            ALLOCATE (matrix_p1(ispin)%matrix)
            ALLOCATE (matrix_p2(ispin)%matrix)
            ALLOCATE (h1_matrix_p2(ispin)%matrix)

            CALL dbcsr_init_p(matrix_p1(ispin)%matrix)
            CALL dbcsr_create(matrix_p1(ispin)%matrix, template=matrix_ks(1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(matrix_p1(ispin)%matrix, sab_orb)
            CALL dbcsr_set(matrix_p1(ispin)%matrix, 0._dp)

            CALL dbcsr_init_p(matrix_p2(ispin)%matrix)
            CALL dbcsr_create(matrix_p2(ispin)%matrix, template=matrix_ks(1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(matrix_p2(ispin)%matrix, sab_orb)
            CALL dbcsr_set(matrix_p2(ispin)%matrix, 0._dp)

            CALL dbcsr_init_p(h1_matrix_p2(ispin)%matrix)
            CALL dbcsr_create(h1_matrix_p2(ispin)%matrix, template=matrix_ks(1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL cp_dbcsr_alloc_block_from_nbl(h1_matrix_p2(ispin)%matrix, sab_all)
            CALL dbcsr_set(h1_matrix_p2(ispin)%matrix, 0._dp)

         END DO

         CALL get_polar_env(polar_env=polar_env, &
                            psi1_dBerry_ao=psi1_dBerry_ao, &
                            dBerry_psi0_ao=dBerry_psi0_ao)

         ! Restart
         DO idir = 1, 3
            DO ispin = 1, nspins
               CALL dbcsr_set(psi1_dBerry_ao(idir, ispin)%matrix, 0.0_dp)
            END DO
         END DO

         loop_idir_vel: DO idir = 1, 3
            IF (iounit > 0) THEN
               WRITE (iounit, "(/,T2,A)") &
                  "POLAR| First uncoupled response to the velocity dipole perturbation operator_ "//ACHAR(idir + 119)
            END IF

            ! Do scf cycle to optimize psi1
            DO ispin = 1, nspins
               CALL dbcsr_desymmetrize(dBerry_psi0_ao(idir, ispin)%matrix, h1_matrix_p(ispin)%matrix)
            END DO

            linres_control%lr_triplet = .FALSE. ! we do singlet response
            linres_control%do_kernel = .FALSE. ! we do uncoupled response
            linres_control%converged = .FALSE.

            should_stop = .FALSE.

            ! initialize p_env
            ! Remark: mos environment is needed for this
            IF (ASSOCIATED(ec_env%p_env)) THEN
               CALL p_env_release(ec_env%p_env)
               DEALLOCATE (ec_env%p_env)
               NULLIFY (ec_env%p_env)
            END IF
            ALLOCATE (ec_env%p_env)
            CALL p_env_create(ec_env%p_env, qs_env, orthogonal_orbitals=.TRUE., &
                              linres_control=linres_control)
            CALL set_qs_env(qs_env, linres_control=linres_control)
            CALL p_env_psi0_changed(ec_env%p_env, qs_env)
            ! Total energy overwritten, replace with Etot from energy correction
            CALL get_qs_env(qs_env, energy=energy)
            energy%total = ec_env%etotal
            !
            p_env_point => ec_env%p_env

            CALL ec_response_ao_vel(qs_env=qs_env, &
                                    p_env=p_env_point, &
                                    matrix_hz=h1_matrix_p, &
                                    matrix_pz=matrix_p1, &
                                    matrix_wz=matrix_wz, &
                                    iounit=iounit, &
                                    should_stop=should_stop)

            ! Copy the response
            DO ispin = 1, nspins
               CALL dbcsr_scale(matrix_p1(ispin)%matrix, -0.5_dp)              
               CALL dbcsr_copy(psi1_dBerry_ao(idir, ispin)%matrix, matrix_p1(ispin)%matrix)          
            END DO

         END DO  loop_idir_vel
      ! END IF ! do_raman

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
         loop_idir_vel2: DO idir = 1, 3

            eps_filter = linres_control%eps_filter

            ! Initial guess for coupled response

            DO ispin = 1,nspins
               CALL dbcsr_set(matrix_p1(ispin)%matrix, 0.0_dp)

               CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_s(1)%matrix, psi1_dBerry_ao(idir, ispin)%matrix, &
                        0.0_dp, h1_matrix_p(ispin)%matrix, filter_eps=eps_filter, retain_sparsity=.TRUE.)

               CALL dbcsr_multiply("N", "N", -1.0_dp, h1_matrix_p(ispin)%matrix, matrix_s(1)%matrix, &
                        0.0_dp, h1_matrix_p2(ispin)%matrix, filter_eps=eps_filter, retain_sparsity=.TRUE.)

               CALL dbcsr_copy(dBerry_psi0_ao(idir, ispin)%matrix, h1_matrix_p2(ispin)%matrix)

            END DO

            IF (iounit > 0) THEN
               WRITE (iounit, "(/,T2,A)") &
                  "POLAR| Second coupled response to the velocity dipole perturbation operator_ "//ACHAR(idir + 119)
            END IF

            linres_control%lr_triplet = .FALSE. ! we do singlet response
            linres_control%do_kernel = .TRUE. ! we do coupled response
            linres_control%converged = .FALSE.

            CALL ec_response_ao_vel(qs_env=qs_env, &
                              p_env=p_env_point, &
                              matrix_hz=h1_matrix_p2, &
                              matrix_pz=matrix_p2, &
                              matrix_wz=matrix_wz, &
                              iounit=iounit, &
                              should_stop=should_stop)

            ! Copy the response
            DO ispin = 1, nspins
               CALL dbcsr_scale(matrix_p2(1)%matrix, 0.5_dp)
               CALL dbcsr_desymmetrize(matrix_p2(1)%matrix, polar_env%psi1_dBerry_ao(idir, ispin)%matrix)
            END DO
         END DO loop_idir_vel2

         DO ispin = 1, nspins
            CALL dbcsr_release(matrix_p1(ispin)%matrix)
            CALL dbcsr_release(matrix_p2(ispin)%matrix)
            CALL dbcsr_release(h1_matrix_p2(ispin)%matrix)
            DEALLOCATE (matrix_p1(ispin)%matrix, matrix_p2(ispin)%matrix, h1_matrix_p2(ispin)%matrix)
         END DO
         DEALLOCATE (matrix_p1, matrix_p2, h1_matrix_p2)         

      END IF ! do_raman_vel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      CALL set_polar_env(polar_env, run_stopped=should_stop)

      ! Clean up
      DO ispin = 1, nspins
         CALL dbcsr_release(h1_matrix_p(ispin)%matrix)
         DEALLOCATE (h1_matrix_p(ispin)%matrix)
      END DO
      DEALLOCATE (h1_matrix_p)

      CALL dbcsr_deallocate_matrix_set(matrix_wz)
      CALL dbcsr_deallocate_matrix_set(matrix_hz)

      CALL cp_print_key_finished_output(iounit, logger, lr_section, &
                                        "PRINT%PROGRAM_RUN_INFO")
      CALL timestop(handle)
   END SUBROUTINE polar_response_ao

! **************************************************************************************************
!> \brief Prints the polarisability tensor to a file during MD runs
!> \param force_env ...
!> \param motion_section ...
!> \param itimes ...
!> \param time ...
!> \param pos ...
!> \param act ...
!> \par History
!>      06.2018 Creation (MK)
!> \author Matthias Krack (MK)
! **************************************************************************************************
   SUBROUTINE write_polarisability_tensor(force_env, motion_section, itimes, time, pos, act)

      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(section_vals_type), POINTER                   :: motion_section
      INTEGER, INTENT(IN)                                :: itimes
      REAL(KIND=dp), INTENT(IN)                          :: time
      CHARACTER(LEN=default_string_length), INTENT(IN), &
         OPTIONAL                                        :: pos, act

      CHARACTER(LEN=default_string_length)               :: my_act, my_pos
      INTEGER                                            :: iounit
      LOGICAL                                            :: do_raman, new_file, run_stopped
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: polar
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(polar_env_type), POINTER                      :: polar_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      NULLIFY (qs_env)

      CALL force_env_get(force_env, qs_env=qs_env)
      IF (ASSOCIATED(qs_env)) THEN
         NULLIFY (polar_env)
         CALL get_qs_env(qs_env=qs_env, polar_env=polar_env)
         IF (ASSOCIATED(polar_env)) THEN
            CALL get_polar_env(polar_env=polar_env, &
                               polar=polar, &
                               do_raman=do_raman, &
                               run_stopped=run_stopped)
            IF (.NOT. run_stopped .AND. do_raman) THEN
               NULLIFY (logger)
               logger => cp_get_default_logger()
               my_pos = "APPEND"
               my_act = "WRITE"
               IF (PRESENT(pos)) my_pos = pos
               IF (PRESENT(act)) my_act = act
               iounit = cp_print_key_unit_nr(logger, motion_section, "PRINT%POLAR_MATRIX", &
                                             extension=".polar", file_position=my_pos, &
                                             file_action=my_act, file_form="FORMATTED", &
                                             is_new_file=new_file)
            ELSE
               iounit = 0
            END IF
            IF (iounit > 0) THEN
               IF (new_file) THEN
                  WRITE (UNIT=iounit, FMT='(A,9(11X,A2," [a.u.]"),6X,A)') &
                     "#   Step   Time [fs]", "xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"
               END IF
               WRITE (UNIT=iounit, FMT='(I8,F12.3,9(1X,F19.8))') itimes, time, &
                  polar(1, 1), polar(1, 2), polar(1, 3), &
                  polar(2, 1), polar(2, 2), polar(2, 3), &
                  polar(3, 1), polar(3, 2), polar(3, 3)
               CALL m_flush(iounit)
               CALL cp_print_key_finished_output(iounit, logger, motion_section, "PRINT%POLAR_MATRIX")
            END IF
         END IF ! polar_env
      END IF ! qs_env

   END SUBROUTINE write_polarisability_tensor

END MODULE qs_linres_polar_utils
