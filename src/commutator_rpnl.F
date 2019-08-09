!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Calculation of the non-local pseudopotential contribution to the core Hamiltonian
!>         <a|V(non-local)|b> = <a|p(l,i)>*h(i,j)*<p(l,j)|b>
!> \par History
!>      - refactered from qs_core_hamiltian [Joost VandeVondele, 2008-11-01]
!>      - full rewrite [jhu, 2009-01-23]
! **************************************************************************************************
MODULE commutator_rpnl
   USE ai_moments,                      ONLY: moment
   USE ai_overlap,                      ONLY: overlap
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE dbcsr_api,                       ONLY: dbcsr_get_block_p,&
                                              dbcsr_p_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE external_potential_types,        ONLY: gth_potential_p_type,&
                                              gth_potential_type,&
                                              sgp_potential_p_type,&
                                              sgp_potential_type
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE sap_kind_types,                  ONLY: alist_type,&
                                              clist_type,&
                                              get_alist,&
                                              release_sap_int,&
                                              sap_int_type,&
                                              sap_sort

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'commutator_rpnl'

   TYPE block_type
      REAL(KIND=dp), DIMENSION(:, :), POINTER              :: block
   END TYPE

   PUBLIC :: build_com_rpnl, build_com_rrpnl, build_com_rpnl_r, build_com_rxrpnl

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param matrix_rv ...
!> \param qs_kind_set ...
!> \param sab_orb ...
!> \param sap_ppnl ...
!> \param eps_ppnl ...
! **************************************************************************************************
   SUBROUTINE build_com_rpnl(matrix_rv, qs_kind_set, sab_orb, sap_ppnl, eps_ppnl, particle_set, cell)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_rv
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      REAL(KIND=dp), INTENT(IN)                          :: eps_ppnl
      TYPE(particle_type), POINTER, DIMENSION(:),&
                                    OPTIONAL             :: particle_set
      TYPE(cell_type), POINTER, OPTIONAL                 :: cell

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_com_rpnl', routineP = moduleN//':'//routineN

      INTEGER :: handle, i, iab, iac, iatom, ibc, icol, ikind, ilist, inode, irow, iset, jatom, &
         jkind, jneighbor, kac, katom, kbc, kkind, l, lc_max, lc_min, ldai, ldsab, lppnl, maxco, &
         maxder, maxl, maxlgto, maxlppnl, maxppnl, maxsgf, mepos, na, nb, ncoa, ncoc, nkind, &
         nlist, nneighbor, nnode, np, nppnl, nprjc, nseta, nsgfa, nthread, prjc, sgfa
      INTEGER, DIMENSION(3)                              :: cell_b, cell_c
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, npgfa, nprj_ppnl, &
                                                            nsgf_seta
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: found, gpot, ppnl_present, spot
      REAL(KIND=dp)                                      :: dac, ppnl_radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: ai_work, sab, work
      REAL(KIND=dp), DIMENSION(1)                        :: rprjc, zetc
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rac, ra, rc !, raf, rcf
      REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha_ppnl, set_radius_a
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cprj, rpgfa, sphi_a, vprj_ppnl, x_block, &
                                                            y_block, z_block, zeta
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: achint, acint, bchint, bcint
      TYPE(alist_type), POINTER                          :: alist_ac, alist_bc
      TYPE(clist_type), POINTER                          :: clist
      TYPE(gth_potential_p_type), DIMENSION(:), POINTER  :: gpotential
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int
      TYPE(sgp_potential_p_type), DIMENSION(:), POINTER  :: spotential
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      CALL timeset(routineN, handle)

      ppnl_present = ASSOCIATED(sap_ppnl)

      IF (ppnl_present) THEN

         nkind = SIZE(qs_kind_set)

         CALL get_qs_kind_set(qs_kind_set, &
                              maxco=maxco, &
                              maxlgto=maxlgto, &
                              maxsgf=maxsgf, &
                              maxlppnl=maxlppnl, &
                              maxppnl=maxppnl)

         maxl = MAX(maxlgto, maxlppnl)
         CALL init_orbital_pointers(maxl+1)

         ldsab = MAX(maxco, ncoset(maxlppnl), maxsgf, maxppnl)
         ldai = ncoset(maxl+1)

         !sap_int needs to be shared as multiple threads need to access this
         ALLOCATE (sap_int(nkind*nkind))
         DO i = 1, nkind*nkind
            NULLIFY (sap_int(i)%alist, sap_int(i)%asort, sap_int(i)%aindex)
            sap_int(i)%nalist = 0
         END DO

         !set up direct access to basis and potential
         ALLOCATE (basis_set(nkind), gpotential(nkind), spotential(nkind))
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
            IF (ASSOCIATED(orb_basis_set)) THEN
               basis_set(ikind)%gto_basis_set => orb_basis_set
            ELSE
               NULLIFY (basis_set(ikind)%gto_basis_set)
            END IF
            CALL get_qs_kind(qs_kind_set(ikind), gth_potential=gth_potential, &
                             sgp_potential=sgp_potential)
            IF (ASSOCIATED(gth_potential)) THEN
               gpotential(ikind)%gth_potential => gth_potential
               NULLIFY (spotential(ikind)%sgp_potential)
            ELSE IF (ASSOCIATED(sgp_potential)) THEN
               spotential(ikind)%sgp_potential => sgp_potential
               NULLIFY (gpotential(ikind)%gth_potential)
            ELSE
               NULLIFY (gpotential(ikind)%gth_potential)
               NULLIFY (spotential(ikind)%sgp_potential)
            END IF
         END DO

         maxder = 4
         nthread = 1
!$       nthread = omp_get_max_threads()

         !calculate the overlap integrals <a|p>
         CALL neighbor_list_iterator_create(nl_iterator, sap_ppnl, nthread=nthread)
!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (nl_iterator, basis_set, spotential, gpotential, maxder, ncoset, &
!$OMP          sap_int, nkind, ldsab,  ldai, nco, cell, particle_set) &
!$OMP PRIVATE (mepos, ikind, kkind, iatom, katom, nlist, ilist, nneighbor, jneighbor, &
!$OMP          cell_c, rac, iac, first_sgfa, la_max, la_min, npgfa, nseta, nsgfa, nsgf_seta, &
!$OMP          sphi_a, zeta, cprj, lppnl, nppnl, nprj_ppnl, &
!$OMP          clist, iset, ncoa, sgfa, prjc, work, sab, ai_work, nprjc,  ppnl_radius, &
!$OMP          ncoc, rpgfa, vprj_ppnl, i, l, gpot, spot, ra, rc, &
!$OMP          set_radius_a, rprjc, dac, lc_max, lc_min, zetc, alpha_ppnl)
         mepos = 0
!$       mepos = omp_get_thread_num()

         ALLOCATE (sab(ldsab, ldsab, maxder), work(ldsab, ldsab, maxder))
         sab = 0.0_dp
         ALLOCATE (ai_work(ldai, ldai, 1))
         ai_work = 0.0_dp

         DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
            CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=kkind, iatom=iatom, &
                                   jatom=katom, nlist=nlist, ilist=ilist, nnode=nneighbor, inode=jneighbor, cell=cell_c, r=rac)
            iac = ikind+nkind*(kkind-1)
            IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
            gpot = ASSOCIATED(gpotential(kkind)%gth_potential)
            spot = ASSOCIATED(spotential(kkind)%sgp_potential)
            IF ((.NOT. gpot) .AND. (.NOT. spot)) CYCLE
            ! get definition of basis set
            first_sgfa => basis_set(ikind)%gto_basis_set%first_sgf
            la_max => basis_set(ikind)%gto_basis_set%lmax
            la_min => basis_set(ikind)%gto_basis_set%lmin
            npgfa => basis_set(ikind)%gto_basis_set%npgf
            nseta = basis_set(ikind)%gto_basis_set%nset
            nsgfa = basis_set(ikind)%gto_basis_set%nsgf
            nsgf_seta => basis_set(ikind)%gto_basis_set%nsgf_set
            rpgfa => basis_set(ikind)%gto_basis_set%pgf_radius
            set_radius_a => basis_set(ikind)%gto_basis_set%set_radius
            sphi_a => basis_set(ikind)%gto_basis_set%sphi
            zeta => basis_set(ikind)%gto_basis_set%zet
            ! get definition of PP projectors
            IF (gpot) THEN
               alpha_ppnl => gpotential(kkind)%gth_potential%alpha_ppnl
               cprj => gpotential(kkind)%gth_potential%cprj
               lppnl = gpotential(kkind)%gth_potential%lppnl
               nppnl = gpotential(kkind)%gth_potential%nppnl
               nprj_ppnl => gpotential(kkind)%gth_potential%nprj_ppnl
               ppnl_radius = gpotential(kkind)%gth_potential%ppnl_radius
               vprj_ppnl => gpotential(kkind)%gth_potential%vprj_ppnl
            ELSEIF (spot) THEN
               CPABORT('SGP not implemented')
            ELSE
               CPABORT('PPNL unknown')
            END IF
!$OMP CRITICAL(sap_int_critical)
            IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) THEN
               sap_int(iac)%a_kind = ikind
               sap_int(iac)%p_kind = kkind
               sap_int(iac)%nalist = nlist
               ALLOCATE (sap_int(iac)%alist(nlist))
               DO i = 1, nlist
                  NULLIFY (sap_int(iac)%alist(i)%clist)
                  sap_int(iac)%alist(i)%aatom = 0
                  sap_int(iac)%alist(i)%nclist = 0
               END DO
            END IF
            IF (.NOT. ASSOCIATED(sap_int(iac)%alist(ilist)%clist)) THEN
               sap_int(iac)%alist(ilist)%aatom = iatom
               sap_int(iac)%alist(ilist)%nclist = nneighbor
               ALLOCATE (sap_int(iac)%alist(ilist)%clist(nneighbor))
               DO i = 1, nneighbor
                  sap_int(iac)%alist(ilist)%clist(i)%catom = 0
               END DO
            END IF
!$OMP END CRITICAL(sap_int_critical)

            ! Take periodic boundary conditions into account
            IF (PRESENT(cell)) THEN
              ra = particle_set(iatom)%r(:)
              rc = particle_set(katom)%r(:)
              rac = pbc(ra, rc, cell)
            ENDIF

            dac = SQRT(SUM(rac*rac))
            clist => sap_int(iac)%alist(ilist)%clist(jneighbor)
            clist%catom = katom
            clist%cell = cell_c
            clist%rac = rac
            ALLOCATE (clist%acint(nsgfa, nppnl, maxder), &
                      clist%achint(nsgfa, nppnl, maxder))
            clist%acint = 0._dp
            clist%achint = 0._dp
            clist%nsgf_cnt = 0
            NULLIFY (clist%sgf_list)

            DO iset = 1, nseta
               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)
               work = 0._dp
               prjc = 1
               DO l = 0, lppnl
                  nprjc = nprj_ppnl(l)*nco(l)
                  IF (nprjc == 0) CYCLE
                  rprjc(1) = ppnl_radius
                  IF (set_radius_a(iset)+rprjc(1) < dac) CYCLE
                  lc_max = l+2*(nprj_ppnl(l)-1)
                  lc_min = l
                  zetc(1) = alpha_ppnl(l)
                  ncoc = ncoset(lc_max)
                  ! Calculate the primitive overlap and dipole moment integrals
                  CALL overlap(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                               lc_max, lc_min, 1, rprjc, zetc, rac, dac, sab(:, :, 1), 0, .FALSE., ai_work, ldai)

                  CALL moment(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                              lc_max, 1, zetc, rprjc, 1, -rac, (/0._dp, 0._dp, 0._dp/), sab(:, :, 2:4))

                  ! *** Transformation step projector functions (cartesian->spherical) ***
                  DO i = 1, maxder
                     CALL dgemm("N", "N", ncoa, nprjc, ncoc, 1.0_dp, sab(1, 1, i), ldsab, &
                                cprj(1, prjc), SIZE(cprj, 1), 0.0_dp, work(1, prjc, i), ldsab)
                  END DO
                  prjc = prjc+nprjc
               END DO
               DO i = 1, maxder
                  ! Contraction step (basis functions)
                  CALL dgemm("T", "N", nsgf_seta(iset), nppnl, ncoa, 1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             work(1, 1, i), ldsab, 0.0_dp, clist%acint(sgfa, 1, i), nsgfa)
                  ! Multiply with interaction matrix(h)
                  CALL dgemm("N", "N", nsgf_seta(iset), nppnl, nppnl, 1.0_dp, clist%acint(sgfa, 1, i), nsgfa, &
                             vprj_ppnl(1, 1), SIZE(vprj_ppnl, 1), 0.0_dp, clist%achint(sgfa, 1, i), nsgfa)
               END DO
            END DO
            clist%maxac = MAXVAL(ABS(clist%acint(:, :, 1)))
            clist%maxach = MAXVAL(ABS(clist%achint(:, :, 1)))
         END DO

         DEALLOCATE (sab, ai_work, work)
!$OMP END PARALLEL
         CALL neighbor_list_iterator_release(nl_iterator)

         ! *** Set up a sorting index
         CALL sap_sort(sap_int)
         ! *** All integrals needed have been calculated and stored in sap_int
         ! *** We now calculate the Hamiltonian matrix elements
         CALL neighbor_list_iterator_create(nl_iterator, sab_orb, nthread=nthread)

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (nl_iterator, basis_set, matrix_rv, &
!$OMP          sap_int, nkind, eps_ppnl ) &
!$OMP PRIVATE (mepos, ikind, jkind, iatom, jatom, nlist, ilist, nnode, inode, cell_b, rab, &
!$OMP          iab, irow, icol, x_block, y_block, z_block, &
!$OMP          found, iac, ibc, alist_ac, alist_bc, acint, bcint, &
!$OMP          achint, bchint, na, np, nb, katom, rac, kkind, kac, kbc, i)

         mepos = 0
!$       mepos = omp_get_thread_num()

         DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
            CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, iatom=iatom, &
                                   jatom=jatom, nlist=nlist, ilist=ilist, nnode=nnode, inode=inode, cell=cell_b, r=rab)
            IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
            IF (.NOT. ASSOCIATED(basis_set(jkind)%gto_basis_set)) CYCLE
            iab = ikind+nkind*(jkind-1)

            ! *** Create matrix blocks for a new matrix block column ***
            ! IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
            ! ELSE
            !    irow = jatom
            !    icol = iatom
            ! END IF
            CALL dbcsr_get_block_p(matrix_rv(1)%matrix, irow, icol, x_block, found)
            CALL dbcsr_get_block_p(matrix_rv(2)%matrix, irow, icol, y_block, found)
            CALL dbcsr_get_block_p(matrix_rv(3)%matrix, irow, icol, z_block, found)

            ! loop over all kinds for projector atom
            IF (ASSOCIATED(x_block) .AND. ASSOCIATED(y_block) .AND. ASSOCIATED(z_block)) THEN
               DO kkind = 1, nkind
                  iac = ikind+nkind*(kkind-1)
                  ibc = jkind+nkind*(kkind-1)
                  IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) CYCLE
                  IF (.NOT. ASSOCIATED(sap_int(ibc)%alist)) CYCLE
                  CALL get_alist(sap_int(iac), alist_ac, iatom)
                  CALL get_alist(sap_int(ibc), alist_bc, jatom)
                  IF (.NOT. ASSOCIATED(alist_ac)) CYCLE
                  IF (.NOT. ASSOCIATED(alist_bc)) CYCLE
                  DO kac = 1, alist_ac%nclist
                     DO kbc = 1, alist_bc%nclist
                        IF (alist_ac%clist(kac)%catom /= alist_bc%clist(kbc)%catom) CYCLE
                        IF (ALL(cell_b+alist_bc%clist(kbc)%cell-alist_ac%clist(kac)%cell == 0)) THEN
                           IF (alist_ac%clist(kac)%maxac*alist_bc%clist(kbc)%maxach < eps_ppnl) CYCLE
                           acint => alist_ac%clist(kac)%acint
                           bcint => alist_bc%clist(kbc)%acint
                           achint => alist_ac%clist(kac)%achint
                           bchint => alist_bc%clist(kbc)%achint
                           na = SIZE(acint, 1)
                           np = SIZE(acint, 2)
                           nb = SIZE(bcint, 1)
!$OMP CRITICAL(h_block_critical)
                           ! IF (iatom <= jatom) THEN
                           ! Vnl*r
                           CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                                      bcint(1, 1, 2), nb, 1.0_dp, x_block, SIZE(x_block, 1))
                           CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                                      bcint(1, 1, 3), nb, 1.0_dp, y_block, SIZE(y_block, 1))
                           CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                                      bcint(1, 1, 4), nb, 1.0_dp, z_block, SIZE(z_block, 1))
                           ! -r*Vnl
                           CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 2), na, &
                                      bcint(1, 1, 1), nb, 1.0_dp, x_block, SIZE(x_block, 1))
                           CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 3), na, &
                                      bcint(1, 1, 1), nb, 1.0_dp, y_block, SIZE(y_block, 1))
                           CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 4), na, &
                                      bcint(1, 1, 1), nb, 1.0_dp, z_block, SIZE(z_block, 1))
                           ! ELSE
                           !    ! Vnl*r
                           !    CALL dgemm("N", "T", nb, na, np, 1.0_dp, bchint(1, 1, 2), nb, &
                           !               acint(1, 1, 1), na, 1.0_dp, x_block, SIZE(x_block, 1))
                           !    CALL dgemm("N", "T", nb, na, np, 1.0_dp, bchint(1, 1, 3), nb, &
                           !               acint(1, 1, 1), na, 1.0_dp, y_block, SIZE(y_block, 1))
                           !    CALL dgemm("N", "T", nb, na, np, 1.0_dp, bchint(1, 1, 4), nb, &
                           !               acint(1, 1, 1), na, 1.0_dp, z_block, SIZE(z_block, 1))
                           !    ! -r*Vnl
                           !    CALL dgemm("N", "T", nb, na, np, -1.0_dp, bchint(1, 1, 1), nb, &
                           !               acint(1, 1, 2), na, 1.0_dp, x_block, SIZE(x_block, 1))
                           !    CALL dgemm("N", "T", nb, na, np, -1.0_dp, bchint(1, 1, 1), nb, &
                           !               acint(1, 1, 3), na, 1.0_dp, y_block, SIZE(y_block, 1))
                           !    CALL dgemm("N", "T", nb, na, np, -1.0_dp, bchint(1, 1, 1), nb, &
                           !               acint(1, 1, 4), na, 1.0_dp, z_block, SIZE(z_block, 1))
                           ! END IF
!$OMP END CRITICAL(h_block_critical)
                           EXIT ! We have found a match and there can be only one single match
                        END IF
                     END DO
                  END DO
               END DO
            ENDIF
         END DO
!$OMP END PARALLEL
         CALL neighbor_list_iterator_release(nl_iterator)

         CALL release_sap_int(sap_int)

         DEALLOCATE (basis_set, gpotential, spotential)

      END IF !ppnl_present

      CALL timestop(handle)

   END SUBROUTINE build_com_rpnl

! **************************************************************************************************
!> \brief ...
!> \param matrix_rv ...
!> \param qs_kind_set ...
!> \param sab_orb ...
!> \param sap_ppnl ...
!> \param eps_ppnl ...
! **************************************************************************************************
   SUBROUTINE build_com_rrpnl(matrix_rv, qs_kind_set, sab_orb, sap_ppnl, eps_ppnl, qs_env, ref_point)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_rv
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      REAL(KIND=dp), INTENT(IN)                          :: eps_ppnl
      TYPE(qs_environment_type), POINTER, OPTIONAL       :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: ref_point

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_com_rrpnl', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, i, iab, iac, iatom, ibc, icol, ikind, ilist, ind, inode, irow, iset, &
         jatom, jkind, jneighbor, kac, katom, kbc, kkind, l, lc_max, lc_min, ldai, ldsab, lppnl, &
         maxco, maxder, maxl, maxlgto, maxlppnl, maxppnl, maxsgf, mepos, na, nb, ncoa, ncoc, &
         nkind, nlist, nneighbor, nnode, np, nppnl, nprjc, nseta, nsgfa, nthread, prjc, sgfa
      INTEGER, DIMENSION(3)                              :: cell_b, cell_c
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, npgfa, nprj_ppnl, &
                                                            nsgf_seta
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: found, gpot, ppnl_present, spot
      REAL(KIND=dp)                                      :: dac, ppnl_radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: ai_work, sab, work
      REAL(KIND=dp), DIMENSION(1)                        :: rprjc, zetc
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rac, rf, ra, rc, raf, rcf
      REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha_ppnl, set_radius_a
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cprj, rpgfa, sphi_a, vprj_ppnl, zeta
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: achint, acint, bchint, bcint
      TYPE(alist_type), POINTER                          :: alist_ac, alist_bc
      TYPE(block_type), DIMENSION(:), POINTER            :: blocks
      TYPE(clist_type), POINTER                          :: clist
      TYPE(gth_potential_p_type), DIMENSION(:), POINTER  :: gpotential
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int
      TYPE(sgp_potential_p_type), DIMENSION(:), POINTER  :: spotential
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
!, &
! x_block, y_block, z_block

      CALL timeset(routineN, handle)

      ppnl_present = ASSOCIATED(sap_ppnl)

      IF (ppnl_present) THEN

         nkind = SIZE(qs_kind_set)

         CALL get_qs_kind_set(qs_kind_set, &
                              maxco=maxco, &
                              maxlgto=maxlgto, &
                              maxsgf=maxsgf, &
                              maxlppnl=maxlppnl, &
                              maxppnl=maxppnl)

         maxl = MAX(maxlgto, maxlppnl)
         CALL init_orbital_pointers(maxl+1)

         ldsab = MAX(maxco, ncoset(maxlppnl), maxsgf, maxppnl)
         ldai = ncoset(maxl+1)

         !sap_int needs to be shared as multiple threads need to access this
         ALLOCATE (sap_int(nkind*nkind))
         DO i = 1, nkind*nkind
            NULLIFY (sap_int(i)%alist, sap_int(i)%asort, sap_int(i)%aindex)
            sap_int(i)%nalist = 0
         END DO

         !set up direct access to basis and potential
         ALLOCATE (basis_set(nkind), gpotential(nkind), spotential(nkind))
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
            IF (ASSOCIATED(orb_basis_set)) THEN
               basis_set(ikind)%gto_basis_set => orb_basis_set
            ELSE
               NULLIFY (basis_set(ikind)%gto_basis_set)
            END IF
            CALL get_qs_kind(qs_kind_set(ikind), gth_potential=gth_potential, &
                             sgp_potential=sgp_potential)
            IF (ASSOCIATED(gth_potential)) THEN
               gpotential(ikind)%gth_potential => gth_potential
               NULLIFY (spotential(ikind)%sgp_potential)
            ELSE IF (ASSOCIATED(sgp_potential)) THEN
               spotential(ikind)%sgp_potential => sgp_potential
               NULLIFY (gpotential(ikind)%gth_potential)
            ELSE
               NULLIFY (gpotential(ikind)%gth_potential)
               NULLIFY (spotential(ikind)%sgp_potential)
            END IF
         END DO

         IF (PRESENT(ref_point)) THEN
            rf(:) = ref_point(:)
         ELSE
            rf(:) = 0._dp
         END IF

         IF (PRESENT(ref_point)) THEN
            CPASSERT(PRESENT(qs_env))
            NULLIFY (cell, particle_set)
            CALL get_qs_env(qs_env,&
                            particle_set = particle_set,&
                            cell = cell)
         ENDIF

         maxder = 10
         nthread = 1
!$       nthread = omp_get_max_threads()

         !calculate the overlap integrals <a|p>
         CALL neighbor_list_iterator_create(nl_iterator, sap_ppnl, nthread=nthread)
!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (nl_iterator, basis_set, spotential, gpotential, maxder, ncoset, &
!$OMP          sap_int, nkind, ldsab,  ldai, nco, rf, cell, particle_set, ref_point) &
!$OMP PRIVATE (mepos, ikind, kkind, iatom, katom, nlist, ilist, nneighbor, jneighbor, &
!$OMP          cell_c, rac, iac, first_sgfa, la_max, la_min, npgfa, nseta, nsgfa, nsgf_seta, &
!$OMP          sphi_a, zeta, cprj, lppnl, nppnl, nprj_ppnl, &
!$OMP          clist, iset, ncoa, sgfa, prjc, work, sab, ai_work, nprjc,  ppnl_radius, &
!$OMP          ncoc, rpgfa, vprj_ppnl, i, l, gpot, spot, raf, rcf, ra, rc, &
!$OMP          set_radius_a, rprjc, dac, lc_max, lc_min, zetc, alpha_ppnl, ind)
         mepos = 0
!$       mepos = omp_get_thread_num()

         ALLOCATE (sab(ldsab, ldsab, maxder), work(ldsab, ldsab, maxder))
         sab = 0.0_dp
         ALLOCATE (ai_work(ldai, ldai, 1))
         ai_work = 0.0_dp

         DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
            CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=kkind, iatom=iatom, &
                                   jatom=katom, nlist=nlist, ilist=ilist, nnode=nneighbor, inode=jneighbor, cell=cell_c, r=rac)
            iac = ikind+nkind*(kkind-1)
            IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
            gpot = ASSOCIATED(gpotential(kkind)%gth_potential)
            spot = ASSOCIATED(spotential(kkind)%sgp_potential)
            IF ((.NOT. gpot) .AND. (.NOT. spot)) CYCLE
            ! get definition of basis set
            first_sgfa => basis_set(ikind)%gto_basis_set%first_sgf
            la_max => basis_set(ikind)%gto_basis_set%lmax
            la_min => basis_set(ikind)%gto_basis_set%lmin
            npgfa => basis_set(ikind)%gto_basis_set%npgf
            nseta = basis_set(ikind)%gto_basis_set%nset
            nsgfa = basis_set(ikind)%gto_basis_set%nsgf
            nsgf_seta => basis_set(ikind)%gto_basis_set%nsgf_set
            rpgfa => basis_set(ikind)%gto_basis_set%pgf_radius
            set_radius_a => basis_set(ikind)%gto_basis_set%set_radius
            sphi_a => basis_set(ikind)%gto_basis_set%sphi
            zeta => basis_set(ikind)%gto_basis_set%zet
            ! get definition of PP projectors
            IF (gpot) THEN
               alpha_ppnl => gpotential(kkind)%gth_potential%alpha_ppnl
               cprj => gpotential(kkind)%gth_potential%cprj
               lppnl = gpotential(kkind)%gth_potential%lppnl
               nppnl = gpotential(kkind)%gth_potential%nppnl
               nprj_ppnl => gpotential(kkind)%gth_potential%nprj_ppnl
               ppnl_radius = gpotential(kkind)%gth_potential%ppnl_radius
               vprj_ppnl => gpotential(kkind)%gth_potential%vprj_ppnl
            ELSEIF (spot) THEN
               CPABORT('SGP not implemented')
            ELSE
               CPABORT('PPNL unknown')
            END IF
!$OMP CRITICAL(sap_int_critical)
            IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) THEN
               sap_int(iac)%a_kind = ikind
               sap_int(iac)%p_kind = kkind
               sap_int(iac)%nalist = nlist
               ALLOCATE (sap_int(iac)%alist(nlist))
               DO i = 1, nlist
                  NULLIFY (sap_int(iac)%alist(i)%clist)
                  sap_int(iac)%alist(i)%aatom = 0
                  sap_int(iac)%alist(i)%nclist = 0
               END DO
            END IF
            IF (.NOT. ASSOCIATED(sap_int(iac)%alist(ilist)%clist)) THEN
               sap_int(iac)%alist(ilist)%aatom = iatom
               sap_int(iac)%alist(ilist)%nclist = nneighbor
               ALLOCATE (sap_int(iac)%alist(ilist)%clist(nneighbor))
               DO i = 1, nneighbor
                  sap_int(iac)%alist(ilist)%clist(i)%catom = 0
               END DO
            END IF
!$OMP END CRITICAL(sap_int_critical)
            dac = SQRT(SUM(rac*rac))
            clist => sap_int(iac)%alist(ilist)%clist(jneighbor)
            clist%catom = katom
            clist%cell = cell_c
            clist%rac = rac
            ALLOCATE (clist%acint(nsgfa, nppnl, maxder), &
                      clist%achint(nsgfa, nppnl, maxder))
            clist%acint = 0._dp
            clist%achint = 0._dp
            clist%nsgf_cnt = 0
            NULLIFY (clist%sgf_list)

            IF (PRESENT(ref_point)) THEN
              ra(:) = pbc(particle_set(iatom)%r(:)-rf, cell)+rf
              rc(:) = pbc(particle_set(katom)%r(:)-rf, cell)+rf
              raf(:) = ra(:)-rf(:)
              rcf(:) = rc(:)-rf(:)
            ENDIF

            DO iset = 1, nseta
               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)
               work = 0._dp
               prjc = 1
               DO l = 0, lppnl
                  nprjc = nprj_ppnl(l)*nco(l)
                  IF (nprjc == 0) CYCLE
                  rprjc(1) = ppnl_radius
                  IF (set_radius_a(iset)+rprjc(1) < dac) CYCLE
                  lc_max = l+2*(nprj_ppnl(l)-1)
                  lc_min = l
                  zetc(1) = alpha_ppnl(l)
                  ncoc = ncoset(lc_max)
                  ! Calculate the primitive overlap and dipole moment integrals
                  CALL overlap(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                               lc_max, lc_min, 1, rprjc, zetc, rac, dac, sab(:, :, 1), 0, .FALSE., ai_work, ldai)
                  IF (.NOT. PRESENT(ref_point)) THEN
                    CALL moment(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lc_max, 1, zetc, rprjc, 2, -rac, (/0._dp, 0._dp, 0._dp/), sab(:, :, 2:maxder))
                  ELSE
                    CALL moment(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lc_max, 1, zetc, rprjc, 2, raf, rcf, sab(:, :, 2:maxder))
                  ENDIF
                  ! *** Transformation step projector functions (cartesian->spherical) ***
                  DO i = 1, maxder
                     CALL dgemm("N", "N", ncoa, nprjc, ncoc, 1.0_dp, sab(1, 1, i), ldsab, &
                                cprj(1, prjc), SIZE(cprj, 1), 0.0_dp, work(1, prjc, i), ldsab)
                  END DO
                  prjc = prjc+nprjc
               END DO
               DO i = 1, maxder
                  ! Contraction step (basis functions)
                  CALL dgemm("T", "N", nsgf_seta(iset), nppnl, ncoa, 1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             work(1, 1, i), ldsab, 0.0_dp, clist%acint(sgfa, 1, i), nsgfa)
                  ! Multiply with interaction matrix(h)
                  CALL dgemm("N", "N", nsgf_seta(iset), nppnl, nppnl, 1.0_dp, clist%acint(sgfa, 1, i), nsgfa, &
                             vprj_ppnl(1, 1), SIZE(vprj_ppnl, 1), 0.0_dp, clist%achint(sgfa, 1, i), nsgfa)
               END DO
            END DO
            clist%maxac = MAXVAL(ABS(clist%acint(:, :, 1)))
            clist%maxach = MAXVAL(ABS(clist%achint(:, :, 1)))
         END DO

         DEALLOCATE (sab, ai_work, work)
!$OMP END PARALLEL
         CALL neighbor_list_iterator_release(nl_iterator)

         ! *** Set up a sorting index
         CALL sap_sort(sap_int)
         ! *** All integrals needed have been calculated and stored in sap_int
         ! *** We now calculate the Hamiltonian matrix elements
         CALL neighbor_list_iterator_create(nl_iterator, sab_orb, nthread=nthread)

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (nl_iterator, basis_set, matrix_rv, &
!$OMP          sap_int, nkind, eps_ppnl ) &
!$OMP PRIVATE (mepos, ikind, jkind, iatom, jatom, nlist, ilist, nnode, inode, cell_b, rab, &
!$OMP          iab, irow, icol, blocks, &
!$OMP          found, iac, ibc, alist_ac, alist_bc, acint, bcint, &
!$OMP          achint, bchint, na, np, nb, katom, rac, kkind, kac, kbc, i, ind)

         mepos = 0
!$       mepos = omp_get_thread_num()
         NULLIFY (blocks)
         ALLOCATE (blocks(6))

         DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
            CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, iatom=iatom, &
                                   jatom=jatom, nlist=nlist, ilist=ilist, nnode=nnode, inode=inode, cell=cell_b, r=rab)
            IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
            IF (.NOT. ASSOCIATED(basis_set(jkind)%gto_basis_set)) CYCLE
            iab = ikind+nkind*(jkind-1)

            ! *** Create matrix blocks for a new matrix block column ***
            ! IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
            ! ELSE
            !    irow = jatom
            !    icol = iatom
            ! END IF

            DO ind = 1, SIZE(matrix_rv)
               CALL dbcsr_get_block_p(matrix_rv(ind)%matrix, irow, icol, blocks(ind)%block, found)
            ENDDO

            ! loop over all kinds for projector atom
            IF (ASSOCIATED(blocks(1)%block) .AND. ASSOCIATED(blocks(2)%block) .AND. ASSOCIATED(blocks(3)%block) &
                .AND. ASSOCIATED(blocks(4)%block) .AND. ASSOCIATED(blocks(5)%block) .AND. ASSOCIATED(blocks(6)%block)) THEN
               DO kkind = 1, nkind
                  iac = ikind+nkind*(kkind-1)
                  ibc = jkind+nkind*(kkind-1)
                  IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) CYCLE
                  IF (.NOT. ASSOCIATED(sap_int(ibc)%alist)) CYCLE
                  CALL get_alist(sap_int(iac), alist_ac, iatom)
                  CALL get_alist(sap_int(ibc), alist_bc, jatom)
                  IF (.NOT. ASSOCIATED(alist_ac)) CYCLE
                  IF (.NOT. ASSOCIATED(alist_bc)) CYCLE
                  DO kac = 1, alist_ac%nclist
                     DO kbc = 1, alist_bc%nclist
                        IF (alist_ac%clist(kac)%catom /= alist_bc%clist(kbc)%catom) CYCLE
                        IF (ALL(cell_b+alist_bc%clist(kbc)%cell-alist_ac%clist(kac)%cell == 0)) THEN
                           IF (alist_ac%clist(kac)%maxac*alist_bc%clist(kbc)%maxach < eps_ppnl) CYCLE
                           acint => alist_ac%clist(kac)%acint
                           bcint => alist_bc%clist(kbc)%acint
                           achint => alist_ac%clist(kac)%achint
                           bchint => alist_bc%clist(kbc)%achint
                           na = SIZE(acint, 1)
                           np = SIZE(acint, 2)
                           nb = SIZE(bcint, 1)
!$OMP CRITICAL(h_block_critical)
                           ! IF (iatom <= jatom) THEN
                           ! Vnl*r
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                           !            bcint(1, 1, 2), nb, 1.0_dp, blocks(1)%block, SIZE(blocks(1)%block, 1))
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                           !            bcint(1, 1, 3), nb, 1.0_dp, blocks(2)%block, SIZE(blocks(2)%block, 1))
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                           !            bcint(1, 1, 4), nb, 1.0_dp, blocks(3)%block, SIZE(blocks(3)%block, 1))
                           ! ! -r*Vnl
                           ! CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 2), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks(1)%block, SIZE(blocks(1)%block, 1))
                           ! CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 3), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks(2)%block, SIZE(blocks(2)%block, 1))
                           ! CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 4), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks(3)%block, SIZE(blocks(3)%block, 1))

                           ! + r_alpha * Vnl * r_beta
                           CALL dgemm("N", "T", na, nb, np, 2._dp, achint(1, 1, 2), na, &
                                      bcint(1, 1, 2), nb, 1.0_dp, blocks(1)%block, SIZE(blocks(1)%block, 1)) ! xVx
                           CALL dgemm("N", "T", na, nb, np, 2._dp, achint(1, 1, 2), na, &
                                      bcint(1, 1, 3), nb, 1.0_dp, blocks(2)%block, SIZE(blocks(2)%block, 1)) ! xVy
                           CALL dgemm("N", "T", na, nb, np, 2._dp, achint(1, 1, 2), na, &
                                      bcint(1, 1, 4), nb, 1.0_dp, blocks(3)%block, SIZE(blocks(3)%block, 1)) ! xVz
                           CALL dgemm("N", "T", na, nb, np, 2._dp, achint(1, 1, 3), na, &
                                      bcint(1, 1, 3), nb, 1.0_dp, blocks(4)%block, SIZE(blocks(4)%block, 1)) ! yVy
                           CALL dgemm("N", "T", na, nb, np, 2._dp, achint(1, 1, 3), na, &
                                      bcint(1, 1, 4), nb, 1.0_dp, blocks(5)%block, SIZE(blocks(5)%block, 1)) ! yVz
                           CALL dgemm("N", "T", na, nb, np, 2._dp, achint(1, 1, 4), na, &
                                      bcint(1, 1, 4), nb, 1.0_dp, blocks(6)%block, SIZE(blocks(6)%block, 1)) ! zVz

                           ! + r_beta * Vnl * r_alpha
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 2), na, &
                           !            bcint(1, 1, 2), nb, 1.0_dp, blocks(1)%block, SIZE(blocks(1)%block, 1))  ! xVx
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 2), na, &
                           !            bcint(1, 1, 3), nb, 1.0_dp, blocks(2)%block, SIZE(blocks(2)%block, 1))  ! xVy
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 2), na, &
                           !            bcint(1, 1, 4), nb, 1.0_dp, blocks(3)%block, SIZE(blocks(3)%block, 1))  ! xVz
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 3), na, &
                           !            bcint(1, 1, 3), nb, 1.0_dp, blocks(4)%block, SIZE(blocks(4)%block, 1))  ! yVy
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 3), na, &
                           !            bcint(1, 1, 4), nb, 1.0_dp, blocks(5)%block, SIZE(blocks(5)%block, 1))  ! yVz
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 4), na, &
                           !            bcint(1, 1, 4), nb, 1.0_dp, blocks(6)%block, SIZE(blocks(6)%block, 1))  ! zVz

                           ! - r_alpha * r_beta * Vnl
                           CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 5), na, &
                                      bcint(1, 1, 1), nb, 1.0_dp, blocks(1)%block, SIZE(blocks(1)%block, 1)) ! xxV
                           CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 6), na, &
                                      bcint(1, 1, 1), nb, 1.0_dp, blocks(2)%block, SIZE(blocks(2)%block, 1)) ! xyV
                           CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 7), na, &
                                      bcint(1, 1, 1), nb, 1.0_dp, blocks(3)%block, SIZE(blocks(3)%block, 1)) ! xzV
                           CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 8), na, &
                                      bcint(1, 1, 1), nb, 1.0_dp, blocks(4)%block, SIZE(blocks(4)%block, 1)) ! yyV
                           CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 9), na, &
                                      bcint(1, 1, 1), nb, 1.0_dp, blocks(5)%block, SIZE(blocks(5)%block, 1)) ! yzV
                           CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 10), na, &
                                      bcint(1, 1, 1), nb, 1.0_dp, blocks(6)%block, SIZE(blocks(6)%block, 1)) ! zzV

                           ! - Vnl * r_alpha * r_beta
                           CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 1), na, &
                                      bcint(1, 1, 5), nb, 1.0_dp, blocks(1)%block, SIZE(blocks(1)%block, 1)) ! Vxx
                           CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 1), na, &
                                      bcint(1, 1, 6), nb, 1.0_dp, blocks(2)%block, SIZE(blocks(2)%block, 1)) ! Vxy
                           CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 1), na, &
                                      bcint(1, 1, 7), nb, 1.0_dp, blocks(3)%block, SIZE(blocks(3)%block, 1)) ! Vxz
                           CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 1), na, &
                                      bcint(1, 1, 8), nb, 1.0_dp, blocks(4)%block, SIZE(blocks(4)%block, 1)) ! Vyy
                           CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 1), na, &
                                      bcint(1, 1, 9), nb, 1.0_dp, blocks(5)%block, SIZE(blocks(5)%block, 1)) ! Vyz
                           CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 1), na, &
                                      bcint(1, 1, 10), nb, 1.0_dp, blocks(6)%block, SIZE(blocks(6)%block, 1)) ! Vzz
                           ! ELSE
                           !    ! Vnl*r
                           !    CALL dgemm("N", "T", nb, na, np, 1.0_dp, bchint(1, 1, 2), nb, &
                           !               acint(1, 1, 1), na, 1.0_dp, x_block, SIZE(x_block, 1))
                           !    CALL dgemm("N", "T", nb, na, np, 1.0_dp, bchint(1, 1, 3), nb, &
                           !               acint(1, 1, 1), na, 1.0_dp, y_block, SIZE(y_block, 1))
                           !    CALL dgemm("N", "T", nb, na, np, 1.0_dp, bchint(1, 1, 4), nb, &
                           !               acint(1, 1, 1), na, 1.0_dp, z_block, SIZE(z_block, 1))
                           !    ! -r*Vnl
                           !    CALL dgemm("N", "T", nb, na, np, -1.0_dp, bchint(1, 1, 1), nb, &
                           !               acint(1, 1, 2), na, 1.0_dp, x_block, SIZE(x_block, 1))
                           !    CALL dgemm("N", "T", nb, na, np, -1.0_dp, bchint(1, 1, 1), nb, &
                           !               acint(1, 1, 3), na, 1.0_dp, y_block, SIZE(y_block, 1))
                           !    CALL dgemm("N", "T", nb, na, np, -1.0_dp, bchint(1, 1, 1), nb, &
                           !               acint(1, 1, 4), na, 1.0_dp, z_block, SIZE(z_block, 1))
                           ! END IF
!$OMP END CRITICAL(h_block_critical)
                           EXIT ! We have found a match and there can be only one single match
                        END IF
                     END DO
                  END DO
               END DO
            ENDIF
         END DO
         DEALLOCATE (blocks)
!$OMP END PARALLEL
         CALL neighbor_list_iterator_release(nl_iterator)

         CALL release_sap_int(sap_int)
         DEALLOCATE (basis_set, gpotential, spotential)

      END IF !ppnl_present

      CALL timestop(handle)

   END SUBROUTINE build_com_rrpnl

! **************************************************************************************************
!> \brief ...
!> \param matrix_rv ...
!> \param qs_kind_set ...
!> \param sab_orb ...
!> \param sap_ppnl ...
!> \param eps_ppnl ...
! **************************************************************************************************
   SUBROUTINE build_com_rxrpnl(qs_env, matrix_rv, qs_kind_set, sab_orb, sap_ppnl, eps_ppnl, ref_point)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_rv
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      REAL(KIND=dp), INTENT(IN)                          :: eps_ppnl
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), OPTIONAL  :: ref_point


      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_com_rxrpnl', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, i, iab, iac, iatom, ibc, icol, ikind, ilist, ind, inode, irow, iset, &
         jatom, jkind, jneighbor, kac, katom, kbc, kkind, l, lc_max, lc_min, ldai, ldsab, lppnl, &
         maxco, maxder, maxl, maxlgto, maxlppnl, maxppnl, maxsgf, mepos, na, nb, ncoa, ncoc, &
         nkind, nlist, nneighbor, nnode, np, nppnl, nprjc, nseta, nsgfa, nthread, prjc, sgfa
      INTEGER, DIMENSION(3)                              :: cell_b, cell_c
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, npgfa, nprj_ppnl, &
                                                            nsgf_seta
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: found, gpot, ppnl_present, spot
      REAL(KIND=dp)                                      :: dac, ppnl_radius!, rac2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: ai_work, sab, work
      REAL(KIND=dp), DIMENSION(1)                        :: rprjc, zetc
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rac, rc, ra, rf, rcf, raf!, rac_new
      REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha_ppnl, set_radius_a
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cprj, rpgfa, sphi_a, vprj_ppnl, zeta
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: achint, acint, bchint, bcint
      TYPE(alist_type), POINTER                          :: alist_ac, alist_bc
      TYPE(block_type), DIMENSION(:), POINTER            :: blocks
      TYPE(clist_type), POINTER                          :: clist
      TYPE(gth_potential_p_type), DIMENSION(:), POINTER  :: gpotential
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int
      TYPE(sgp_potential_p_type), DIMENSION(:), POINTER  :: spotential
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

!, &
! x_block, y_block, z_block

      CALL timeset(routineN, handle)

      ppnl_present = ASSOCIATED(sap_ppnl)

      IF (ppnl_present) THEN

         nkind = SIZE(qs_kind_set)

         CALL get_qs_kind_set(qs_kind_set, &
                              maxco=maxco, &
                              maxlgto=maxlgto, &
                              maxsgf=maxsgf, &
                              maxlppnl=maxlppnl, &
                              maxppnl=maxppnl)

         maxl = MAX(maxlgto, maxlppnl)
         CALL init_orbital_pointers(maxl+1)

         ldsab = MAX(maxco, ncoset(maxlppnl), maxsgf, maxppnl)
         ldai = ncoset(maxl+1)

         !sap_int needs to be shared as multiple threads need to access this
         ALLOCATE (sap_int(nkind*nkind))
         DO i = 1, nkind*nkind
            NULLIFY (sap_int(i)%alist, sap_int(i)%asort, sap_int(i)%aindex)
            sap_int(i)%nalist = 0
         END DO

         !set up direct access to basis and potential
         ALLOCATE (basis_set(nkind), gpotential(nkind), spotential(nkind))
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
            IF (ASSOCIATED(orb_basis_set)) THEN
               basis_set(ikind)%gto_basis_set => orb_basis_set
            ELSE
               NULLIFY (basis_set(ikind)%gto_basis_set)
            END IF
            CALL get_qs_kind(qs_kind_set(ikind), gth_potential=gth_potential, &
                             sgp_potential=sgp_potential)
            IF (ASSOCIATED(gth_potential)) THEN
               gpotential(ikind)%gth_potential => gth_potential
               NULLIFY (spotential(ikind)%sgp_potential)
            ELSE IF (ASSOCIATED(sgp_potential)) THEN
               spotential(ikind)%sgp_potential => sgp_potential
               NULLIFY (gpotential(ikind)%gth_potential)
            ELSE
               NULLIFY (gpotential(ikind)%gth_potential)
               NULLIFY (spotential(ikind)%sgp_potential)
            END IF
         END DO

         IF (PRESENT(ref_point)) THEN
            rf(:) = ref_point(:)
         ELSE
            rf(:) = 0._dp
         END IF

         NULLIFY (cell, particle_set)

         CALL get_qs_env(qs_env,&
                         particle_set = particle_set,&
                         cell = cell)

         maxder = 10
         nthread = 1
!$       nthread = omp_get_max_threads()

         !calculate the overlap integrals <a|p>
         CALL neighbor_list_iterator_create(nl_iterator, sap_ppnl, nthread=nthread)
!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (nl_iterator, basis_set, spotential, gpotential, maxder, ncoset, &
!$OMP          sap_int, nkind, ldsab,  ldai, nco, rf, cell, particle_set, ref_point) &
!$OMP PRIVATE (mepos, ikind, kkind, iatom, katom, nlist, ilist, nneighbor, jneighbor, &
!$OMP          cell_c, rac, iac, first_sgfa, la_max, la_min, npgfa, nseta, nsgfa, nsgf_seta, &
!$OMP          sphi_a, zeta, cprj, lppnl, nppnl, nprj_ppnl, &
!$OMP          clist, iset, ncoa, sgfa, prjc, work, sab, ai_work, nprjc,  ppnl_radius, &
!$OMP          ncoc, rpgfa, vprj_ppnl, i, l, gpot, spot, ra, rc, raf, rcf, &
!$OMP          set_radius_a, rprjc, dac, lc_max, lc_min, zetc, alpha_ppnl, ind)
         mepos = 0
!$       mepos = omp_get_thread_num()

         ALLOCATE (sab(ldsab, ldsab, maxder), work(ldsab, ldsab, maxder))
         sab = 0.0_dp
         ALLOCATE (ai_work(ldai, ldai, 1))
         ai_work = 0.0_dp

         DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
            CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=kkind, iatom=iatom, &
                                   jatom=katom, nlist=nlist, ilist=ilist, nnode=nneighbor, inode=jneighbor, cell=cell_c, r=rac)
            iac = ikind+nkind*(kkind-1)
            IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
            gpot = ASSOCIATED(gpotential(kkind)%gth_potential)
            spot = ASSOCIATED(spotential(kkind)%sgp_potential)
            IF ((.NOT. gpot) .AND. (.NOT. spot)) CYCLE
            ! get definition of basis set
            first_sgfa => basis_set(ikind)%gto_basis_set%first_sgf
            la_max => basis_set(ikind)%gto_basis_set%lmax
            la_min => basis_set(ikind)%gto_basis_set%lmin
            npgfa => basis_set(ikind)%gto_basis_set%npgf
            nseta = basis_set(ikind)%gto_basis_set%nset
            nsgfa = basis_set(ikind)%gto_basis_set%nsgf
            nsgf_seta => basis_set(ikind)%gto_basis_set%nsgf_set
            rpgfa => basis_set(ikind)%gto_basis_set%pgf_radius
            set_radius_a => basis_set(ikind)%gto_basis_set%set_radius
            sphi_a => basis_set(ikind)%gto_basis_set%sphi
            zeta => basis_set(ikind)%gto_basis_set%zet
            ! get definition of PP projectors
            IF (gpot) THEN
               alpha_ppnl => gpotential(kkind)%gth_potential%alpha_ppnl
               cprj => gpotential(kkind)%gth_potential%cprj
               lppnl = gpotential(kkind)%gth_potential%lppnl
               nppnl = gpotential(kkind)%gth_potential%nppnl
               nprj_ppnl => gpotential(kkind)%gth_potential%nprj_ppnl
               ppnl_radius = gpotential(kkind)%gth_potential%ppnl_radius
               vprj_ppnl => gpotential(kkind)%gth_potential%vprj_ppnl
            ELSEIF (spot) THEN
               CPABORT('SGP not implemented')
            ELSE
               CPABORT('PPNL unknown')
            END IF
!$OMP CRITICAL(sap_int_critical)
            IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) THEN
               sap_int(iac)%a_kind = ikind
               sap_int(iac)%p_kind = kkind
               sap_int(iac)%nalist = nlist
               ALLOCATE (sap_int(iac)%alist(nlist))
               DO i = 1, nlist
                  NULLIFY (sap_int(iac)%alist(i)%clist)
                  sap_int(iac)%alist(i)%aatom = 0
                  sap_int(iac)%alist(i)%nclist = 0
               END DO
            END IF
            IF (.NOT. ASSOCIATED(sap_int(iac)%alist(ilist)%clist)) THEN
               sap_int(iac)%alist(ilist)%aatom = iatom
               sap_int(iac)%alist(ilist)%nclist = nneighbor
               ALLOCATE (sap_int(iac)%alist(ilist)%clist(nneighbor))
               DO i = 1, nneighbor
                  sap_int(iac)%alist(ilist)%clist(i)%catom = 0
               END DO
            END IF
!$OMP END CRITICAL(sap_int_critical)
            dac = SQRT(SUM(rac*rac))
            clist => sap_int(iac)%alist(ilist)%clist(jneighbor)
            clist%catom = katom
            clist%cell = cell_c
            clist%rac = rac
            ALLOCATE (clist%acint(nsgfa, nppnl, maxder), &
                      clist%achint(nsgfa, nppnl, maxder))
            clist%acint = 0._dp
            clist%achint = 0._dp
            clist%nsgf_cnt = 0
            NULLIFY (clist%sgf_list)


            ! using PBC here might screw a molecule that fits the box (but e.g. hasn't been shifted by center_molecule)
            ! by folding around the center, such screwing can be avoided for a proper choice of center.
            ra(:) = pbc(particle_set(iatom)%r(:)-rf, cell)+rf
            rc(:) = pbc(particle_set(katom)%r(:)-rf, cell)+rf
            ! we dont use PBC at this point
            ! rac_new(:) = ra(:)-rc(:)
            raf(:) = ra(:)-rf(:)
            rcf(:) = rc(:)-rf(:)
            ! rac2 = rac_new(1)*rac_new(1)+raf(2)*raf(2)+rcf(3)*rcf(3)
            ! dab = SQRT(rab2)

            DO iset = 1, nseta
               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)
               work = 0._dp
               prjc = 1
               DO l = 0, lppnl
                  nprjc = nprj_ppnl(l)*nco(l)
                  IF (nprjc == 0) CYCLE
                  rprjc(1) = ppnl_radius
                  IF (set_radius_a(iset)+rprjc(1) < dac) CYCLE
                  lc_max = l+2*(nprj_ppnl(l)-1)
                  lc_min = l
                  zetc(1) = alpha_ppnl(l)
                  ncoc = ncoset(lc_max)
                  ! Calculate the primitive overlap and dipole moment integrals
                  CALL overlap(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                               lc_max, lc_min, 1, rprjc, zetc, rac, dac, sab(:, :, 1), 0, .FALSE., ai_work, ldai)
                  IF (.NOT. PRESENT(ref_point)) THEN
                    CALL moment(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lc_max, 1, zetc, rprjc, 2, -rac, (/0._dp, 0._dp, 0._dp/), sab(:, :, 2:maxder))
                  ELSE
                    CALL moment(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lc_max, 1, zetc, rprjc, 2, raf, rcf, sab(:, :, 2:maxder))
                  ENDIF
                  ! *** Transformation step projector functions (cartesian->spherical) ***
                  DO i = 1, maxder
                     CALL dgemm("N", "N", ncoa, nprjc, ncoc, 1.0_dp, sab(1, 1, i), ldsab, &
                                cprj(1, prjc), SIZE(cprj, 1), 0.0_dp, work(1, prjc, i), ldsab)
                  END DO
                  prjc = prjc+nprjc
               END DO
               DO i = 1, maxder
                  ! Contraction step (basis functions)
                  CALL dgemm("T", "N", nsgf_seta(iset), nppnl, ncoa, 1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             work(1, 1, i), ldsab, 0.0_dp, clist%acint(sgfa, 1, i), nsgfa)
                  ! Multiply with interaction matrix(h)
                  CALL dgemm("N", "N", nsgf_seta(iset), nppnl, nppnl, 1.0_dp, clist%acint(sgfa, 1, i), nsgfa, &
                             vprj_ppnl(1, 1), SIZE(vprj_ppnl, 1), 0.0_dp, clist%achint(sgfa, 1, i), nsgfa)
               END DO
            END DO
            clist%maxac = MAXVAL(ABS(clist%acint(:, :, 1)))
            clist%maxach = MAXVAL(ABS(clist%achint(:, :, 1)))
         END DO

         DEALLOCATE (sab, ai_work, work)
!$OMP END PARALLEL
         CALL neighbor_list_iterator_release(nl_iterator)

         ! *** Set up a sorting index
         CALL sap_sort(sap_int)
         ! *** All integrals needed have been calculated and stored in sap_int
         ! *** We now calculate the Hamiltonian matrix elements
         CALL neighbor_list_iterator_create(nl_iterator, sab_orb, nthread=nthread)

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (nl_iterator, basis_set, matrix_rv, &
!$OMP          sap_int, nkind, eps_ppnl ) &
!$OMP PRIVATE (mepos, ikind, jkind, iatom, jatom, nlist, ilist, nnode, inode, cell_b, rab, &
!$OMP          iab, irow, icol, blocks, &
!$OMP          found, iac, ibc, alist_ac, alist_bc, acint, bcint, &
!$OMP          achint, bchint, na, np, nb, katom, rac, kkind, kac, kbc, i, ind)

         mepos = 0
!$       mepos = omp_get_thread_num()
         NULLIFY (blocks)
         ALLOCATE (blocks(3))

         DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
            CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, iatom=iatom, &
                                   jatom=jatom, nlist=nlist, ilist=ilist, nnode=nnode, inode=inode, cell=cell_b, r=rab)
            IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
            IF (.NOT. ASSOCIATED(basis_set(jkind)%gto_basis_set)) CYCLE
            iab = ikind+nkind*(jkind-1)

            ! *** Create matrix blocks for a new matrix block column ***
            ! IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
            ! ELSE
            !    irow = jatom
            !    icol = iatom
            ! END IF

            DO ind = 1, SIZE(matrix_rv)
               CALL dbcsr_get_block_p(matrix_rv(ind)%matrix, irow, icol, blocks(ind)%block, found)
            ENDDO

            ! loop over all kinds for projector atom
            IF (ASSOCIATED(blocks(1)%block) .AND. ASSOCIATED(blocks(2)%block) .AND. ASSOCIATED(blocks(3)%block)) THEN
               DO kkind = 1, nkind
                  iac = ikind+nkind*(kkind-1)
                  ibc = jkind+nkind*(kkind-1)
                  IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) CYCLE
                  IF (.NOT. ASSOCIATED(sap_int(ibc)%alist)) CYCLE
                  CALL get_alist(sap_int(iac), alist_ac, iatom)
                  CALL get_alist(sap_int(ibc), alist_bc, jatom)
                  IF (.NOT. ASSOCIATED(alist_ac)) CYCLE
                  IF (.NOT. ASSOCIATED(alist_bc)) CYCLE
                  DO kac = 1, alist_ac%nclist
                     DO kbc = 1, alist_bc%nclist
                        IF (alist_ac%clist(kac)%catom /= alist_bc%clist(kbc)%catom) CYCLE
                        IF (ALL(cell_b+alist_bc%clist(kbc)%cell-alist_ac%clist(kac)%cell == 0)) THEN
                           IF (alist_ac%clist(kac)%maxac*alist_bc%clist(kbc)%maxach < eps_ppnl) CYCLE
                           acint => alist_ac%clist(kac)%acint
                           bcint => alist_bc%clist(kbc)%acint
                           achint => alist_ac%clist(kac)%achint
                           bchint => alist_bc%clist(kbc)%achint
                           na = SIZE(acint, 1)
                           np = SIZE(acint, 2)
                           nb = SIZE(bcint, 1)
!$OMP CRITICAL(h_block_critical)
                           ! x-component
                           CALL dgemm("N", "T", na, nb, np, 1.0_dp, achint(1, 1, 9), na, &
                                      bcint(1, 1, 1), nb, 1.0_dp, blocks(1)%block, SIZE(blocks(1)%block, 1)) ! yzV
                           CALL dgemm("N", "T", na, nb, np, -1.0_dp, achint(1, 1, 3), na, &
                                      bcint(1, 1, 4), nb, 1.0_dp, blocks(1)%block, SIZE(blocks(1)%block, 1)) ! -yVz
                           CALL dgemm("N", "T", na, nb, np, -1.0_dp, achint(1, 1, 9), na, &
                                      bcint(1, 1, 1), nb, 1.0_dp, blocks(1)%block, SIZE(blocks(1)%block, 1)) ! -zyV
                           CALL dgemm("N", "T", na, nb, np, 1.0_dp, achint(1, 1, 4), na, &
                                      bcint(1, 1, 3), nb, 1.0_dp, blocks(1)%block, SIZE(blocks(1)%block, 1)) ! zVy

                           ! y-component
                           CALL dgemm("N", "T", na, nb, np, 1.0_dp, achint(1, 1, 7), na, &
                                      bcint(1, 1, 1), nb, 1.0_dp, blocks(2)%block, SIZE(blocks(2)%block, 1)) ! zxV
                           CALL dgemm("N", "T", na, nb, np, -1.0_dp, achint(1, 1, 4), na, &
                                      bcint(1, 1, 2), nb, 1.0_dp, blocks(2)%block, SIZE(blocks(2)%block, 1)) ! -zVx
                           CALL dgemm("N", "T", na, nb, np, -1.0_dp, achint(1, 1, 7), na, &
                                      bcint(1, 1, 1), nb, 1.0_dp, blocks(2)%block, SIZE(blocks(2)%block, 1)) ! -xzV
                           CALL dgemm("N", "T", na, nb, np, 1.0_dp, achint(1, 1, 2), na, &
                                      bcint(1, 1, 4), nb, 1.0_dp, blocks(2)%block, SIZE(blocks(2)%block, 1)) ! xVz

                           ! z-component
                           CALL dgemm("N", "T", na, nb, np, 1.0_dp, achint(1, 1, 6), na, &
                                      bcint(1, 1, 1), nb, 1.0_dp, blocks(3)%block, SIZE(blocks(3)%block, 1)) ! xyV
                           CALL dgemm("N", "T", na, nb, np, -1.0_dp, achint(1, 1, 2), na, &
                                      bcint(1, 1, 3), nb, 1.0_dp, blocks(3)%block, SIZE(blocks(3)%block, 1)) ! -xVy
                           CALL dgemm("N", "T", na, nb, np, -1.0_dp, achint(1, 1, 6), na, &
                                      bcint(1, 1, 1), nb, 1.0_dp, blocks(3)%block, SIZE(blocks(3)%block, 1)) ! -yxV
                           CALL dgemm("N", "T", na, nb, np, 1.0_dp, achint(1, 1, 3), na, &
                                      bcint(1, 1, 2), nb, 1.0_dp, blocks(3)%block, SIZE(blocks(3)%block, 1)) ! yVx
!$OMP END CRITICAL(h_block_critical)
                           EXIT ! We have found a match and there can be only one single match
                        END IF
                     END DO
                  END DO
               END DO
            ENDIF
         END DO
         DEALLOCATE (blocks)
!$OMP END PARALLEL
         CALL neighbor_list_iterator_release(nl_iterator)

         CALL release_sap_int(sap_int)

         DEALLOCATE (basis_set, gpotential, spotential)

      END IF !ppnl_present

      CALL timestop(handle)

   END SUBROUTINE build_com_rxrpnl

! **************************************************************************************************
!> \brief ...
!> \param matrix_rv ...
!> \param qs_kind_set ...
!> \param sab_nl ...
!> \param sap_ppnl ...
!> \param eps_ppnl ...
!> \param direction_Or ...
! **************************************************************************************************
   SUBROUTINE build_com_rpnl_r(matrix_rv, qs_kind_set, sab_nl, sap_ppnl, eps_ppnl, direction_Or)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_rv
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl, sap_ppnl
      REAL(KIND=dp), INTENT(IN)                          :: eps_ppnl
      LOGICAL, INTENT(IN)                                :: direction_Or

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_com_rpnl_r', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, i, iab, iac, iatom, ibc, icol, ikind, ilist, inode, irow, iset, jatom, &
         jkind, jneighbor, kac, katom, kbc, kkind, l, lc_max, lc_min, ldai, ldsab, lppnl, maxco, &
         maxder, maxl, maxlgto, maxlppnl, maxppnl, maxsgf, mepos, na, nb, ncoa, ncoc, nkind, &
         nlist, nneighbor, nnode, np, nppnl, nprjc, nseta, nsgfa, nthread, prjc, sgfa
      INTEGER, DIMENSION(3)                              :: cell_b, cell_c
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, npgfa, nprj_ppnl, &
                                                            nsgf_seta
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: found, gpot, ppnl_present, spot
      REAL(KIND=dp)                                      :: dac, ppnl_radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: ai_work, sab, work
      REAL(KIND=dp), DIMENSION(1)                        :: rprjc, zetc
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rac
      REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha_ppnl, set_radius_a
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: cprj, rpgfa, sphi_a, vprj_ppnl, xx_block, &
         xy_block, xz_block, yx_block, yy_block, yz_block, zeta, zx_block, zy_block, zz_block
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: achint, acint, bchint, bcint
      TYPE(alist_type), POINTER                          :: alist_ac, alist_bc
      TYPE(clist_type), POINTER                          :: clist
      TYPE(gth_potential_p_type), DIMENSION(:), POINTER  :: gpotential
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int
      TYPE(sgp_potential_p_type), DIMENSION(:), POINTER  :: spotential
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

! x_block, y_block, z_block,

      CALL timeset(routineN, handle)

      ppnl_present = ASSOCIATED(sap_ppnl)

      IF (ppnl_present) THEN

         nkind = SIZE(qs_kind_set)

         CALL get_qs_kind_set(qs_kind_set, &
                              maxco=maxco, &
                              maxlgto=maxlgto, &
                              maxsgf=maxsgf, &
                              maxlppnl=maxlppnl, &
                              maxppnl=maxppnl)

         maxl = MAX(maxlgto, maxlppnl)
         CALL init_orbital_pointers(maxl+1)

         ldsab = MAX(maxco, ncoset(maxlppnl), maxsgf, maxppnl)
         ldai = ncoset(maxl+1)

         !sap_int needs to be shared as multiple threads need to access this
         ALLOCATE (sap_int(nkind*nkind))
         DO i = 1, nkind*nkind
            NULLIFY (sap_int(i)%alist, sap_int(i)%asort, sap_int(i)%aindex)
            sap_int(i)%nalist = 0
         END DO

         !set up direct access to basis and potential
         ALLOCATE (basis_set(nkind), gpotential(nkind), spotential(nkind))
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
            IF (ASSOCIATED(orb_basis_set)) THEN
               basis_set(ikind)%gto_basis_set => orb_basis_set
            ELSE
               NULLIFY (basis_set(ikind)%gto_basis_set)
            END IF
            CALL get_qs_kind(qs_kind_set(ikind), gth_potential=gth_potential, &
                             sgp_potential=sgp_potential)
            IF (ASSOCIATED(gth_potential)) THEN
               gpotential(ikind)%gth_potential => gth_potential
               NULLIFY (spotential(ikind)%sgp_potential)
            ELSE IF (ASSOCIATED(sgp_potential)) THEN
               spotential(ikind)%sgp_potential => sgp_potential
               NULLIFY (gpotential(ikind)%gth_potential)
            ELSE
               NULLIFY (gpotential(ikind)%gth_potential)
               NULLIFY (spotential(ikind)%sgp_potential)
            END IF
         END DO

!ovlp + dipole + quadrupole
         maxder = 10
         nthread = 1
!$       nthread = omp_get_max_threads()

         !calculate the overlap integrals <a|p>
         CALL neighbor_list_iterator_create(nl_iterator, sap_ppnl, nthread=nthread)
!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (nl_iterator, basis_set, spotential, gpotential, maxder, ncoset, &
!$OMP          sap_int, nkind, ldsab,  ldai, nco) &
!$OMP PRIVATE (mepos, ikind, kkind, iatom, katom, nlist, ilist, nneighbor, jneighbor, &
!$OMP          cell_c, rac, iac, first_sgfa, la_max, la_min, npgfa, nseta, nsgfa, nsgf_seta, &
!$OMP          sphi_a, zeta, cprj, lppnl, nppnl, nprj_ppnl, &
!$OMP          clist, iset, ncoa, sgfa, prjc, work, sab, ai_work, nprjc,  ppnl_radius, &
!$OMP          ncoc, rpgfa, vprj_ppnl, i, l, gpot, spot, &
!$OMP          set_radius_a, rprjc, dac, lc_max, lc_min, zetc, alpha_ppnl)
         mepos = 0
!$       mepos = omp_get_thread_num()

         ALLOCATE (sab(ldsab, ldsab, maxder), work(ldsab, ldsab, maxder))
         sab = 0.0_dp
         work = 0.0_dp
         ALLOCATE (ai_work(ldai, ldai, 1))
         ai_work = 0.0_dp

         DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
            CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=kkind, iatom=iatom, &
                                   jatom=katom, nlist=nlist, ilist=ilist, nnode=nneighbor, inode=jneighbor, cell=cell_c, r=rac)
            iac = ikind+nkind*(kkind-1)
            IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
            gpot = ASSOCIATED(gpotential(kkind)%gth_potential)
            spot = ASSOCIATED(spotential(kkind)%sgp_potential)
            IF ((.NOT. gpot) .AND. (.NOT. spot)) CYCLE

            sab = 0.0_dp
            work = 0.0_dp
            ai_work = 0.0_dp

            ! get definition of basis set
            first_sgfa => basis_set(ikind)%gto_basis_set%first_sgf
            la_max => basis_set(ikind)%gto_basis_set%lmax
            la_min => basis_set(ikind)%gto_basis_set%lmin
            npgfa => basis_set(ikind)%gto_basis_set%npgf
            nseta = basis_set(ikind)%gto_basis_set%nset
            nsgfa = basis_set(ikind)%gto_basis_set%nsgf
            nsgf_seta => basis_set(ikind)%gto_basis_set%nsgf_set
            rpgfa => basis_set(ikind)%gto_basis_set%pgf_radius
            set_radius_a => basis_set(ikind)%gto_basis_set%set_radius
            sphi_a => basis_set(ikind)%gto_basis_set%sphi
            zeta => basis_set(ikind)%gto_basis_set%zet
            ! get definition of PP projectors
            IF (gpot) THEN
               alpha_ppnl => gpotential(kkind)%gth_potential%alpha_ppnl
               cprj => gpotential(kkind)%gth_potential%cprj
               lppnl = gpotential(kkind)%gth_potential%lppnl
               nppnl = gpotential(kkind)%gth_potential%nppnl
               nprj_ppnl => gpotential(kkind)%gth_potential%nprj_ppnl
               ppnl_radius = gpotential(kkind)%gth_potential%ppnl_radius
               vprj_ppnl => gpotential(kkind)%gth_potential%vprj_ppnl
            ELSEIF (spot) THEN
               CPABORT('SGP not implemented')
            ELSE
               CPABORT('PPNL unknown')
            END IF
!$OMP CRITICAL(sap_int_critical)
            IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) THEN
               sap_int(iac)%a_kind = ikind
               sap_int(iac)%p_kind = kkind
               sap_int(iac)%nalist = nlist
               ALLOCATE (sap_int(iac)%alist(nlist))
               DO i = 1, nlist
                  NULLIFY (sap_int(iac)%alist(i)%clist)
                  sap_int(iac)%alist(i)%aatom = 0
                  sap_int(iac)%alist(i)%nclist = 0
               END DO
            END IF
            IF (.NOT. ASSOCIATED(sap_int(iac)%alist(ilist)%clist)) THEN
               sap_int(iac)%alist(ilist)%aatom = iatom
               sap_int(iac)%alist(ilist)%nclist = nneighbor
               ALLOCATE (sap_int(iac)%alist(ilist)%clist(nneighbor))
               DO i = 1, nneighbor
                  sap_int(iac)%alist(ilist)%clist(i)%catom = 0
               END DO
            END IF
!$OMP END CRITICAL(sap_int_critical)

            dac = SQRT(SUM(rac*rac))
            clist => sap_int(iac)%alist(ilist)%clist(jneighbor)
            clist%catom = katom
            clist%cell = cell_c
            clist%rac = rac
            ALLOCATE (clist%acint(nsgfa, nppnl, maxder), &
                      clist%achint(nsgfa, nppnl, maxder))
            clist%acint = 0._dp
            clist%achint = 0._dp
            clist%nsgf_cnt = 0
            NULLIFY (clist%sgf_list)
            DO iset = 1, nseta
               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)
               work = 0._dp
               prjc = 1
               DO l = 0, lppnl
                  nprjc = nprj_ppnl(l)*nco(l)
                  IF (nprjc == 0) CYCLE
                  rprjc(1) = ppnl_radius
                  IF (set_radius_a(iset)+rprjc(1) < dac) CYCLE
                  lc_max = l+2*(nprj_ppnl(l)-1)
                  lc_min = l
                  zetc(1) = alpha_ppnl(l)
                  ncoc = ncoset(lc_max)
                  ! Calculate the primitive overlap and dipole moment integrals
                  CALL overlap(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                               lc_max, lc_min, 1, rprjc, zetc, rac, dac, sab(:, :, 1), 0, .FALSE., ai_work, ldai)
                  CALL moment(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                              lc_max, 1, zetc, rprjc, 2, -rac, (/0._dp, 0._dp, 0._dp/), sab(:, :, 2:10))
                  ! *** Transformation step projector functions (cartesian->spherical) ***
                  DO i = 1, maxder
                     CALL dgemm("N", "N", ncoa, nprjc, ncoc, 1.0_dp, sab(1, 1, i), ldsab, &
                                cprj(1, prjc), SIZE(cprj, 1), 0.0_dp, work(1, prjc, i), ldsab)
                  END DO
                  prjc = prjc+nprjc
               END DO
               DO i = 1, maxder
                  ! Contraction step (basis functions)
                  CALL dgemm("T", "N", nsgf_seta(iset), nppnl, ncoa, 1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             work(1, 1, i), ldsab, 0.0_dp, clist%acint(sgfa, 1, i), nsgfa)
                  ! Multiply with interaction matrix(h)
                  CALL dgemm("N", "N", nsgf_seta(iset), nppnl, nppnl, 1.0_dp, clist%acint(sgfa, 1, i), nsgfa, &
                             vprj_ppnl(1, 1), SIZE(vprj_ppnl, 1), 0.0_dp, clist%achint(sgfa, 1, i), nsgfa)
               END DO
            END DO
            clist%maxac = MAXVAL(ABS(clist%acint(:, :, 1)))
            clist%maxach = MAXVAL(ABS(clist%achint(:, :, 1)))

         END DO

         CALL neighbor_list_iterator_release(nl_iterator)

         DEALLOCATE (sab, ai_work, work)
!$OMP END PARALLEL

         ! *** Set up a sorting index
         CALL sap_sort(sap_int)
         ! *** All integrals needed have been calculated and stored in sap_int
         ! *** We now calculate the Hamiltonian matrix elements
         CALL neighbor_list_iterator_create(nl_iterator, sab_nl, nthread=nthread)

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (nl_iterator, basis_set, matrix_rv, &
!$OMP          sap_int, nkind, eps_ppnl, direction_Or) &
!$OMP PRIVATE (mepos, ikind, jkind, iatom, jatom, nlist, ilist, nnode, inode, cell_b, rab, &
!$OMP          iab, irow, icol, xx_block, xy_block, xz_block, yx_block, yy_block, yz_block, &
!$OMP          zx_block, zy_block, zz_block, found, iac, ibc, alist_ac, alist_bc, acint, bcint, &
!$OMP          achint, bchint, na, np, nb, katom, rac, kkind, kac, kbc, i)

         mepos = 0
!$       mepos = omp_get_thread_num()

         DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
            CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, iatom=iatom, &
                                   jatom=jatom, nlist=nlist, ilist=ilist, nnode=nnode, inode=inode, cell=cell_b, r=rab)
            IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
            IF (.NOT. ASSOCIATED(basis_set(jkind)%gto_basis_set)) CYCLE
            iab = ikind+nkind*(jkind-1)

            ! *** Create matrix blocks for a new matrix block column ***
!            IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
!            ELSE
!               irow = jatom
!               icol = iatom
!            END IF
!            IF (iatom <= jatom) THEN
!              print *, 'iatom <= jatom'
!            END IF
            CALL dbcsr_get_block_p(matrix_rv(1)%matrix, irow, icol, xx_block, found)
            CALL dbcsr_get_block_p(matrix_rv(2)%matrix, irow, icol, yx_block, found)
            CALL dbcsr_get_block_p(matrix_rv(3)%matrix, irow, icol, zx_block, found)
            CALL dbcsr_get_block_p(matrix_rv(4)%matrix, irow, icol, xy_block, found)
            CALL dbcsr_get_block_p(matrix_rv(5)%matrix, irow, icol, yy_block, found)
            CALL dbcsr_get_block_p(matrix_rv(6)%matrix, irow, icol, zy_block, found)
            CALL dbcsr_get_block_p(matrix_rv(7)%matrix, irow, icol, xz_block, found)
            CALL dbcsr_get_block_p(matrix_rv(8)%matrix, irow, icol, yz_block, found)
            CALL dbcsr_get_block_p(matrix_rv(9)%matrix, irow, icol, zz_block, found)

            ! loop over all kinds for projector atom
            IF (ASSOCIATED(xx_block) .AND. ASSOCIATED(zz_block)) THEN
               DO kkind = 1, nkind
                  iac = ikind+nkind*(kkind-1)
                  ibc = jkind+nkind*(kkind-1)
                  IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) CYCLE
                  IF (.NOT. ASSOCIATED(sap_int(ibc)%alist)) CYCLE
                  CALL get_alist(sap_int(iac), alist_ac, iatom)
                  CALL get_alist(sap_int(ibc), alist_bc, jatom)
                  IF (.NOT. ASSOCIATED(alist_ac)) CYCLE
                  IF (.NOT. ASSOCIATED(alist_bc)) CYCLE
                  DO kac = 1, alist_ac%nclist
                     DO kbc = 1, alist_bc%nclist
                        IF (alist_ac%clist(kac)%catom /= alist_bc%clist(kbc)%catom) CYCLE
                        IF (ALL(cell_b+alist_bc%clist(kbc)%cell-alist_ac%clist(kac)%cell == 0)) THEN
                           IF (alist_ac%clist(kac)%maxac*alist_bc%clist(kbc)%maxach < eps_ppnl) CYCLE
                           acint => alist_ac%clist(kac)%acint
                           bcint => alist_bc%clist(kbc)%acint
                           achint => alist_ac%clist(kac)%achint
                           bchint => alist_bc%clist(kbc)%achint
                           na = SIZE(acint, 1)
                           np = SIZE(acint, 2)
                           nb = SIZE(bcint, 1)
!$OMP CRITICAL(h_block_critical)
!                           IF (iatom <= jatom) THEN
                           IF (direction_Or) THEN
                              ! Vnl*r_beta*r_alpha
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                                         bcint(1, 1, 5), nb, 1.0_dp, xx_block, SIZE(xx_block, 1))
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                                         bcint(1, 1, 6), nb, 1.0_dp, xy_block, SIZE(xy_block, 1))
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                                         bcint(1, 1, 7), nb, 1.0_dp, xz_block, SIZE(xz_block, 1))

                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                                         bcint(1, 1, 6), nb, 1.0_dp, yx_block, SIZE(yx_block, 1))
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                                         bcint(1, 1, 8), nb, 1.0_dp, yy_block, SIZE(yy_block, 1))
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                                         bcint(1, 1, 9), nb, 1.0_dp, yz_block, SIZE(yz_block, 1))

                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                                         bcint(1, 1, 7), nb, 1.0_dp, zx_block, SIZE(zx_block, 1))
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                                         bcint(1, 1, 9), nb, 1.0_dp, zy_block, SIZE(zy_block, 1))
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                                         bcint(1, 1, 10), nb, 1.0_dp, zz_block, SIZE(zz_block, 1))

                              ! -r_beta*Vnl*r_alpha
                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 2), na, &
                                         bcint(1, 1, 2), nb, 1.0_dp, xx_block, SIZE(xx_block, 1))
                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 3), na, &
                                         bcint(1, 1, 2), nb, 1.0_dp, xy_block, SIZE(xy_block, 1))
                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 4), na, &
                                         bcint(1, 1, 2), nb, 1.0_dp, xz_block, SIZE(xz_block, 1))

                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 2), na, &
                                         bcint(1, 1, 3), nb, 1.0_dp, yx_block, SIZE(yx_block, 1))
                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 3), na, &
                                         bcint(1, 1, 3), nb, 1.0_dp, yy_block, SIZE(yy_block, 1))
                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 4), na, &
                                         bcint(1, 1, 3), nb, 1.0_dp, yz_block, SIZE(yz_block, 1))

                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 2), na, &
                                         bcint(1, 1, 4), nb, 1.0_dp, zx_block, SIZE(zx_block, 1))
                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 3), na, &
                                         bcint(1, 1, 4), nb, 1.0_dp, zy_block, SIZE(zy_block, 1))
                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 4), na, &
                                         bcint(1, 1, 4), nb, 1.0_dp, zz_block, SIZE(zz_block, 1))

                           ELSE
                              ! r_alpha*Vnl*r_beta
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 2), na, &
                                         bcint(1, 1, 2), nb, 1.0_dp, xx_block, SIZE(xx_block, 1))
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 2), na, &
                                         bcint(1, 1, 3), nb, 1.0_dp, xy_block, SIZE(xy_block, 1))
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 2), na, &
                                         bcint(1, 1, 4), nb, 1.0_dp, xz_block, SIZE(xz_block, 1))

                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 3), na, &
                                         bcint(1, 1, 2), nb, 1.0_dp, yx_block, SIZE(yx_block, 1))
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 3), na, &
                                         bcint(1, 1, 3), nb, 1.0_dp, yy_block, SIZE(yy_block, 1))
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 3), na, &
                                         bcint(1, 1, 4), nb, 1.0_dp, yz_block, SIZE(yz_block, 1))

                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 4), na, &
                                         bcint(1, 1, 2), nb, 1.0_dp, zx_block, SIZE(zx_block, 1))
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 4), na, &
                                         bcint(1, 1, 3), nb, 1.0_dp, zy_block, SIZE(zy_block, 1))
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 4), na, &
                                         bcint(1, 1, 4), nb, 1.0_dp, zz_block, SIZE(zz_block, 1))

                              ! -rbeta*r_alpha*Vnl
                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 5), na, &
                                         bcint(1, 1, 1), nb, 1.0_dp, xx_block, SIZE(xx_block, 1))
                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 6), na, &
                                         bcint(1, 1, 1), nb, 1.0_dp, xy_block, SIZE(xy_block, 1))
                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 7), na, &
                                         bcint(1, 1, 1), nb, 1.0_dp, xz_block, SIZE(xz_block, 1))

                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 6), na, &
                                         bcint(1, 1, 1), nb, 1.0_dp, yx_block, SIZE(yx_block, 1))
                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 8), na, &
                                         bcint(1, 1, 1), nb, 1.0_dp, yy_block, SIZE(yy_block, 1))
                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 9), na, &
                                         bcint(1, 1, 1), nb, 1.0_dp, yz_block, SIZE(yz_block, 1))

                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 7), na, &
                                         bcint(1, 1, 1), nb, 1.0_dp, zx_block, SIZE(zx_block, 1))
                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 9), na, &
                                         bcint(1, 1, 1), nb, 1.0_dp, zy_block, SIZE(zy_block, 1))
                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 10), na, &
                                         bcint(1, 1, 1), nb, 1.0_dp, zz_block, SIZE(zz_block, 1))
                           END IF

!$OMP END CRITICAL(h_block_critical)
                           EXIT ! We have found a match and there can be only one single match
                        END IF
                     END DO
                  END DO
               END DO
            ENDIF
         END DO
!$OMP END PARALLEL
         CALL neighbor_list_iterator_release(nl_iterator)

         DEALLOCATE (basis_set, gpotential, spotential)

      END IF !ppnl_present

      CALL timestop(handle)

   END SUBROUTINE build_com_rpnl_r
END MODULE commutator_rpnl
