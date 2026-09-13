!-----------------------------------------------------------------------------------
! MODULE: MultiscaleDemag
!> @brief Computes the demagnetization field for the micromagnetic (FD) region.
!>
!> The magnetization lives directly on the finite-difference mesh nodes.
!> For every real FD node (finiteDiffIndices(i,j,k) > 0) the magnetization vector
!> is read from emomM and placed into Mgrid(1:3, i, j, k)/cell_volume.
!> The demagnetizing field is then computed as:
!>
!> H(i,j,k) = -sum_{i',j',k'} N(i-i', j-j', k-k') . M(i',j',k')
!>
!> using the exact Newell tensor (near-field, R/h < far_threshold) and the point-dipole
!> approximation (far-field, R/h >= far_threshold).
!>
!> The full Newell tensor is PRECOMPUTED once at initialisation and stored in
!> Ntensor(1:6, -(nx-1):(nx-1), -(ny-1):(ny-1), -(nz-1):(nz-1)).
!> Indices: 1=Nxx, 2=Nyy, 3=Nzz, 4=Nxy, 5=Nxz, 6=Nyz.
!> At each timestep calc_multiscale_demag only performs multiply-add lookups,
!> avoiding recalculating and atan/log/sqrt calls inside the simulation loop.
!-------------------------------------------------------------------------------------
module MultiscaleDemag
    use Parameters, only: dblprec
    use Constants,  only: mu0, mry, mub, pi !4pix10^-7 H/m, 2.179872325d-21 J, 9.274009994d-24 J/T, 3.141592653589793_dblprec

    implicit none
    private

    !> Stores the active cell mask and the demag field arrays.
    type :: MultiscaleDemagData
        logical  :: enabled = .false.
        logical  :: output_enabled = .false.
        integer  :: nxc, nyc, nzc   ! FD grid cell counts
        real(dblprec) :: dx, dy, dz    ! Cell sizes
        real(dblprec) :: volume    ! dx*dy*dz cell volume

        integer  :: nxn, nyn, nzn ! FD grid node dimensions

        ! node_index(i,j,k) = atom index from finiteDiffIndices.
        !  > 0 : real continuum FD node
        !  < 0 : ghost / interpolation node
        !  = 0 : empty (no atom)
        integer, allocatable :: node_index(:,:,:)   !(nxn, nyn, nzn)

        ! active_node(i,j,k) = atom index for that FD cell (> 0 = real node).
        ! Stored so we can still look up emomM at each timestep after
        ! finiteDiffIndices has been deallocated in the setup phase.
        ! integer, allocatable :: active_node(:,:,:)   ! (nxc, nyc, nzc)

        !unitCell information
        real(dblprec) :: unitcell_atoms
        real(dblprec) :: unitcell_volume

        ! active_cell(i,j,k) = 1 if all surrounding corner nodes are real, 0 otherwise.
        ! Dimensioned over the cell grid (nxc, nyc, nzc).
        integer, allocatable :: active_cell(:,:,:)  ! (nxc, nyc, nzc)

        ! Cell Centered Magnetization and Demag Field arrays.
        real(dblprec), allocatable :: Mgrid(:,:,:,:)  ! (3, nxc, nyc, nzc)
        real(dblprec), allocatable :: Hgrid(:,:,:,:)  ! (3, nxc, nyc, nzc)

        ! Precomputed Newell tensor kernel.
        ! Ntensor(c, di, dj, dk): component c for displacement (di, dj, dk)
        ! c=1:Nxx, c=2:Nyy, c=3:Nzz, c=4:Nxy, c=5:Nxz, c=6:Nyz
        ! Bounds: di in -(ni-1):(ni-1)
        real(dblprec), allocatable :: Ntensor(:,:,:,:) ! (6, -(nxc-1):(nxc-1), -(nyc-1):(nyc-1), -(nzc-1):(nzc-1))

        ! Final mu0*H mapped back to every atom index for beff accumulation.
        ! Only cells where active_node > 0 are filled; atomistic atoms are zero.
        integer :: Natom
        ! real(dblprec), allocatable :: Bdemag(:,:)    ! (3, Natom)

        ! Per-step energy accumulators (updated by calc_ms_energies).
        real(dblprec) :: E_demag = 0.0_dblprec  !< Demagnetization energy
        real(dblprec) :: E_ani   = 0.0_dblprec  !< Uniaxial anisotropy energy
        real(dblprec) :: E_xc    = 0.0_dblprec  !< Exchange energy
        real(dblprec) :: E_xc_pen= 0.0_dblprec  !< Penalty gradient exchange energy
        real(dblprec) :: E_dm    = 0.0_dblprec  !< DM energy
        real(dblprec) :: E_total = 0.0_dblprec  !< E_demag + E_ani + E_xc + E_dm
    end type MultiscaleDemagData

    type(MultiscaleDemagData), save :: msd

    public :: setup_multiscale_demag,  &
              build_ms_magnetization,  &
              calc_multiscale_demag,   &
              add_multiscale_demag_to_beff, &
              calc_ms_energies, &
              cleanup_multiscale_demag


contains

    ! Must be called before finiteDiffIndices and mesh are deallocated.
    ! @param[in] nx,ny,nz      Cell grid dimensions  (mesh%nrOfBoxes)
    ! @param[in] dx,dy,dz      Cell sizes            (mesh%boxSize)
    ! @param[in] nxn,nyn,nzn   Node grid dimensions  (mesh%nrOfGridPoints)
    ! @param[in] fd_indices     finiteDiffIndices(nxn,nyn,nzn):
    !                             > 0  real FD node  (atom index)
    !                             < 0  ghost/interpolation node
    !                             = 0  empty cell
    ! @param[in] Natom       Total number of atoms/nodes in the system
    ! @param[in] unitCell    Unit cell information
    subroutine setup_multiscale_demag(nxc, nyc, nzc, dx, dy, dz, nxn, nyn, nzn, fd_indices, Natom, unitCell)
        use InputData, only: ms_demag, ms_demag_output, far_threshold, ham_inp
        use AtomGenerator, only: AtomCell
        type(AtomCell), intent(in) :: unitCell
        integer, intent(in) :: nxc, nyc, nzc
        real(dblprec), intent(in) :: dx, dy, dz
        integer, intent(in) :: nxn, nyn, nzn
        integer, intent(in) :: fd_indices(nxn, nyn, nzn)
        integer, intent(in) :: Natom

        integer :: di, dj, dk, i, j, k
        integer :: di2, dj2, dk2, khi
        real(dblprec) :: rx, ry, rz, R, h_max
        real(dblprec), parameter :: pi4inv    = 1.0_dblprec / (4.0_dblprec * pi)
        ! real(dblprec), parameter :: far_thresh = 20.0_dblprec
        logical :: all_real

        ! Demagnetization tensor output file
        integer :: file_unit_N
        character(len=30) :: filn_N

        write(*,*) 'unitcell%nrOfAtoms: ', unitCell%nrOfAtoms, ' unitcell%size: ', unitCell%size, ' unitcell volume: ', unitCell%size(1)*unitCell%size(2)*unitCell%size(3)
        if (unitCell%nrOfAtoms <= 0) then
            write(*,*) 'ERROR: MultiscaleDemag requires a non-empty unit cell. Please check your multiscale.conf file.'
            stop
        end if

        write(*,*) 'DEBUG: setup_multiscale_demag called. ms_demag is: ', ms_demag

        if (ms_demag == 'Y') then
            msd%enabled = .true.
            write(*,*) 'Demag field will be computed.'
            ! return
        end if

        if (ms_demag_output == 'Y') then
            msd%output_enabled = .true.
            write(*,*) 'Continuum energies will be computed.'
            ! return
        end if

        if (.not.(ms_demag == 'Y' .or. ms_demag_output == 'Y')) then
            write(*,*) 'MultiscaleDemag: ms_demag and ms_demag_output are disabled. Continuum energies and demag field will not be computed.'
            return
        end if

        ! Grid and cell sizes
        ! msd%enabled = .true.
        msd%nxc = nxc;  msd%nyc = nyc;  msd%nzc = nzc
        msd%dx = dx;  msd%dy = dy;  msd%dz = dz
        msd%volume = dx * dy * dz
        msd%Natom  = Natom

        ! Node grid
        msd%nxn = nxn;  msd%nyn = nyn;  msd%nzn = nzn

        !unit cell information
        msd%unitcell_atoms = unitCell%nrOfAtoms
        msd%unitcell_volume = unitCell%size(1)*unitCell%size(2)*unitCell%size(3)

        !Warn if standard do_dip is also active
        if (ham_inp%do_dip > 0) then
            write(*,*) 'ms_demag and do_dip are both enabled.'
        end if

        !Store node-index map 
        allocate(msd%node_index(nxn, nyn, nzn))
        msd%node_index = fd_indices

        !Build cell-active map
        ! Cell (i,j,k) is active if all surrounding corner nodes are real (>0).
        ! In 2-D (nzn=1, nzc=1) we check the 4 XY-plane corners (dk2 loop skipped).
        ! In 3-D we check all 8 corners.
        allocate(msd%active_cell(nxc, nyc, nzc))
        msd%active_cell = 0
        khi = merge(0, 1, nzc == 1)   ! 0 or 1 depending on dimensionality
        do k = 1, nzc
            do j = 1, nyc
                do i = 1, nxc
                    all_real = .true.
                    do dk2 = 0, khi
                        do dj2 = 0, 1
                            do di2 = 0, 1
                                if (fd_indices(i+di2, j+dj2, k+dk2) <= 0) then
                                    ! write(*,*) 'DEBUG: ',i+di2, j+dj2, k+dk2, fd_indices(i+di2, j+dj2, k+dk2)
                                    all_real = .false.
                                end if
                            end do
                        end do
                    end do
                    if (all_real) msd%active_cell(i, j, k) = 1
                end do
            end do
        end do

        allocate(msd%Mgrid(3, nxc, nyc, nzc));  msd%Mgrid = 0.0_dblprec
        allocate(msd%Hgrid(3, nxc, nyc, nzc));  msd%Hgrid = 0.0_dblprec
        ! allocate(msd%Bdemag(3, Natom));      msd%Bdemag = 0.0_dblprec

        ! Precompute the full Newell interaction tensor for all displacement
        ! vectors (di, dj, dk) that can occur in the convolution sum.
        ! This is done ONCE here so the simulation loop is just a lookup.
        allocate(msd%Ntensor(6, -(nxc-1):(nxc-1), -(nyc-1):(nyc-1), -(nzc-1):(nzc-1)))
        msd%Ntensor = 0.0_dblprec

        h_max = max(dx, max(dy, dz))

        write(*,'(a)') ' MultiscaleDemag: Precomputing Newell tensor'
        write(*,*) 'Far field threshold is: ', real(far_threshold, dblprec)

        !Writing the output file for the anisotropy energies
        filn_N = "ms_demag_tensor.out"
        open(newunit=file_unit_N, file=trim(filn_N), status='replace') 
        write(file_unit_N,'(7(a8,10x))') 'Nxx', 'Nyy', 'Nzz', 'Nxy', 'Nxz', 'Nyz', 'Ntrace'
        flush(file_unit_N)

        !$omp parallel do collapse(3) default(shared) private(di, dj, dk, rx, ry, rz, R)
        do dk = -(nzc-1), (nzc-1)
            do dj = -(nyc-1), (nyc-1)
                do di = -(nxc-1), (nxc-1)
                    rx = real(di, dblprec) * dx
                    ry = real(dj, dblprec) * dy
                    rz = real(dk, dblprec) * dz
                    R  = sqrt(rx*rx + ry*ry + rz*rz)

                    if (R < 1.0d-12) then
                        ! Self-interaction: Newell self term
                        msd%Ntensor(1, di, dj, dk) = apply_L(.true.,  0.0_dblprec, 0.0_dblprec, 0.0_dblprec, dx, dy, dz)
                        msd%Ntensor(2, di, dj, dk) = apply_L(.true.,  0.0_dblprec, 0.0_dblprec, 0.0_dblprec, dy, dz, dx)
                        msd%Ntensor(3, di, dj, dk) = apply_L(.true.,  0.0_dblprec, 0.0_dblprec, 0.0_dblprec, dz, dx, dy)
                        msd%Ntensor(4, di, dj, dk) = 0.0_dblprec
                        msd%Ntensor(5, di, dj, dk) = 0.0_dblprec
                        msd%Ntensor(6, di, dj, dk) = 0.0_dblprec
                        write(*,*) 'DEBUG: Nxx,Nyy,Nzz, Nxx+Nyy+Nzz: ', msd%Ntensor(1, di, dj, dk), &
                                    msd%Ntensor(2, di, dj, dk), msd%Ntensor(3, di, dj, dk), msd%Ntensor(1, di, dj, dk) &
                                     + msd%Ntensor(2, di, dj, dk) + msd%Ntensor(3, di, dj, dk)

                    else if (R / h_max > real(far_threshold, dblprec)) then
                        ! Far-field: point-dipole approximation
                        msd%Ntensor(1, di, dj, dk) = -msd%volume * pi4inv * (3.0_dblprec*(rx/R)**2 - 1.0_dblprec) / R**3
                        msd%Ntensor(2, di, dj, dk) = -msd%volume * pi4inv * (3.0_dblprec*(ry/R)**2 - 1.0_dblprec) / R**3
                        msd%Ntensor(3, di, dj, dk) = -msd%volume * pi4inv * (3.0_dblprec*(rz/R)**2 - 1.0_dblprec) / R**3
                        msd%Ntensor(4, di, dj, dk) = -msd%volume * pi4inv *  3.0_dblprec*rx*ry / R**5
                        msd%Ntensor(5, di, dj, dk) = -msd%volume * pi4inv *  3.0_dblprec*rx*rz / R**5
                        msd%Ntensor(6, di, dj, dk) = -msd%volume * pi4inv *  3.0_dblprec*ry*rz / R**5

                    else
                        ! Near-field: exact Newell tensor
                        msd%Ntensor(1, di, dj, dk) = apply_L(.true.,  rx, ry, rz, dx, dy, dz)
                        msd%Ntensor(2, di, dj, dk) = apply_L(.true.,  ry, rz, rx, dy, dz, dx)
                        msd%Ntensor(3, di, dj, dk) = apply_L(.true.,  rz, rx, ry, dz, dx, dy)
                        msd%Ntensor(4, di, dj, dk) = apply_L(.false., rx, ry, rz, dx, dy, dz)
                        msd%Ntensor(5, di, dj, dk) = apply_L(.false., rx, rz, ry, dx, dz, dy)
                        msd%Ntensor(6, di, dj, dk) = apply_L(.false., ry, rz, rx, dy, dz, dx)
                    end if

                    write(file_unit_N,'(7(es18.9,2x))') msd%Ntensor(1, di, dj, dk), msd%Ntensor(2, di, dj, dk), msd%Ntensor(3, di, dj, dk), &
                                                    msd%Ntensor(4, di, dj, dk), msd%Ntensor(5, di, dj, dk), msd%Ntensor(6, di, dj, dk), &
                                                    msd%Ntensor(1, di, dj, dk) + msd%Ntensor(2, di, dj, dk) + msd%Ntensor(3, di, dj, dk)
                    flush(file_unit_N)
                end do
            end do
        end do
        !$omp end parallel do

        write(*,'(a)')       ' MultiscaleDemag: FD-based demag solver initialised (cell-centred).'
        write(*,'(a,3i5)')   '   Cell grid  (nrOfBoxes)      : ', nxc,  nyc,  nzc
        write(*,'(a,3i5)')   '   Node grid  (nrOfGridPoints) : ', nxn, nyn, nzn
        write(*,'(a,3f12.6)')'   Cell sizes (Angstrom)       : ', dx,  dy,  dz
        write(*,'(a,i5)')    '   Active cells                : ', sum(msd%active_cell)
    end subroutine setup_multiscale_demag


    ! Build M(r) on the FD grid directly from emomM.
    ! Only real FD nodes contribute.
    ! param[in] emomM  Magnetic moment vector * magnitude, shape (3, Natom, Mensemble)
    ! param[in] ens    Ensemble index to use
    subroutine build_ms_magnetization(emomM, ens)
        real(dblprec), intent(in) :: emomM(3, msd%Natom, *)
        integer,       intent(in) :: ens

        integer       :: i, j, k, di2, dj2, dk2, idx, cnt
        integer       :: khi
        real(dblprec) :: m_avg(3)

        if (.not. msd%enabled) return

        msd%Mgrid = 0.0_dblprec

        khi = merge(0, 1, msd%nzc == 1)
        !$omp parallel do collapse(3) default(shared) private(i, j, k, di2, dj2, dk2, idx, cnt, m_avg)
        do k = 1, msd%nzc
            do j = 1, msd%nyc
                do i = 1, msd%nxc
                    ! idx = msd%active_node(i, j, k)
                    ! if (idx > 0) then ! emomM = m_s * e_hat  (unit: Bohr-magneton equivalent)
                    ! msd%Mgrid(1:3, i, j, k) = emomM(1:3, idx, ens)/msd%volume ! Divide by cell volume to get magnetization density M.
                    ! end if
                    if (msd%active_cell(i, j, k) /= 1) cycle
                    m_avg = 0.0_dblprec
                    cnt   = 0
                    do dk2 = 0, khi
                        do dj2 = 0, 1
                            do di2 = 0, 1
                                idx = msd%node_index(i+di2, j+dj2, k+dk2)
                                if (idx > 0) then
                                    m_avg = m_avg + emomM(1:3, idx, ens)
                                    cnt   = cnt + 1
                                end if
                            end do
                        end do
                    end do
                    if (cnt > 0) then
                        !average moment of surrounding nodes(/cnt), then multiply by atom density to get magnetization density M     (cell volume)
                        msd%Mgrid(1:3, i, j, k) = (m_avg / (real(cnt, dblprec)))*(msd%unitcell_atoms/msd%unitcell_volume) ! * msd%volume
                        ! write(*,*) 'DEBUG: ', cnt
                    end if
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine build_ms_magnetization


    ! Compute the demagnetizing field H on the FD grid
    ! H_i = -sum_j N_ij . M_j
    ! Uses the precomputed Ntensor for O(1) lookup per cell pair
    subroutine calc_multiscale_demag()
        integer       :: i,  j,  k
        integer       :: i2, j2, k2
        real(dblprec) :: Nxx, Nyy, Nzz, Nxy, Nxz, Nyz
        real(dblprec) :: mx, my, mz

        if (.not. msd%enabled) return

        msd%Hgrid = 0.0_dblprec

        ! $omp parallel do collapse(3) default(shared) &
        ! $omp private(i, j, k, i2, j2, k2, mx, my, mz, Nxx, Nyy, Nzz, Nxy, Nxz, Nyz, local_Hx, local_Hy, local_Hz)
        !$omp parallel do collapse(3) default(shared) &
        !$omp private(i, j, k, i2, j2, k2, mx, my, mz, Nxx, Nyy, Nzz, Nxy, Nxz, Nyz)

        do k  = 1, msd%nzc
            do j  = 1, msd%nyc
                do i  = 1, msd%nxc
                    if (msd%active_cell(i, j, k) /= 1) cycle
                    do k2 = 1, msd%nzc
                        do j2 = 1, msd%nyc
                            do i2 = 1, msd%nxc
                                if (msd%active_cell(i2, j2, k2) /= 1) then
                                    ! write(*,*) 'DEBUG: Skipping ghost source cell: ', i2, j2, k2, msd%active_cell(i2, j2, k2)
                                    cycle
                                end if

                                mx = msd%Mgrid(1, i2, j2, k2)
                                my = msd%Mgrid(2, i2, j2, k2)
                                mz = msd%Mgrid(3, i2, j2, k2)

                                if (abs(mx) < 1.0d-12 .and. abs(my) < 1.0d-12 .and. abs(mz) < 1.0d-12) cycle ! Skip empty source cells

                                ! Tensor lookup:
                                Nxx = msd%Ntensor(1, i-i2, j-j2, k-k2)
                                Nyy = msd%Ntensor(2, i-i2, j-j2, k-k2)
                                Nzz = msd%Ntensor(3, i-i2, j-j2, k-k2)
                                Nxy = msd%Ntensor(4, i-i2, j-j2, k-k2)
                                Nxz = msd%Ntensor(5, i-i2, j-j2, k-k2)
                                Nyz = msd%Ntensor(6, i-i2, j-j2, k-k2)

                                msd%Hgrid(1,i,j,k) = msd%Hgrid(1,i,j,k) - (Nxx*mx + Nxy*my + Nxz*mz)
                                msd%Hgrid(2,i,j,k) = msd%Hgrid(2,i,j,k) - (Nxy*mx + Nyy*my + Nyz*mz)
                                msd%Hgrid(3,i,j,k) = msd%Hgrid(3,i,j,k) - (Nxz*mx + Nyz*my + Nzz*mz)

                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine calc_multiscale_demag


    ! Convert Hgrid to mu0*H and add into beff for every real FD node.
    ! Each corner node (in, jn, kn) receives the average H of all adjacent active cells.  
    subroutine add_multiscale_demag_to_beff(beff)
        real(dblprec), intent(inout) :: beff(3, msd%Natom)

        integer       :: in, jn, kn, di2, dj2, dk2, ic, jc, kc, idx, cnt
        integer       :: khi
        real(dblprec) :: H_node(3)

        if (.not. msd%enabled) return

        ! !$omp parallel do collapse(3) default(shared) private(i, j, k, idx)
        ! do k = 1, msd%nz
        !     do j = 1, msd%ny
        !         do i = 1, msd%nx
        !             idx = msd%active_node(i, j, k)
        !             if (idx == 1) then
        !                 beff(1:3, idx) = beff(1:3, idx) + mu0 * (mub / 1.0d-30) * msd%Hgrid(1:3, i, j, k)
        khi = merge(1, 0, msd%nzc == 1)
        !$omp parallel do collapse(3) default(shared) private(in, jn, kn, di2, dj2, dk2, ic, jc, kc, idx, cnt, H_node)
        do kn = 1, msd%nzn
            do jn = 1, msd%nyn
                do in = 1, msd%nxn
                    idx = msd%node_index(in, jn, kn)
                    if (idx <= 0) cycle   ! not a real atom
                    ! write(*,*) 'DEBUG: ',in, jn, kn, msd%node_index(in, jn, kn)
                    H_node = 0.0_dblprec
                    cnt    = 0
                    do dk2 = 0, khi
                        do dj2 = 0, 1
                            do di2 = 0, 1
                                ic = in + di2 - 1
                                jc = jn + dj2 - 1
                                kc = kn + dk2 - 1
                                if (ic < 1 .or. ic > msd%nxc) cycle
                                if (jc < 1 .or. jc > msd%nyc) cycle
                                if (kc < 1 .or. kc > msd%nzc) cycle
                                if (msd%active_cell(ic, jc, kc) /= 1) cycle
                                H_node = H_node + msd%Hgrid(1:3, ic, jc, kc)
                                cnt = cnt + 1
                                ! write(*,*) 'DEBUG: ', cnt
                            end do
                        end do
                    end do
                    if (cnt > 0) then
                        beff(1:3, idx) = beff(1:3, idx) + mu0 * (mub / 1.0d-30) * H_node / real(cnt, dblprec)
                        ! write(*,*) 'DEBUG: ', cnt
                    end if
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine add_multiscale_demag_to_beff

    ! Calculate demagnetization , uniaxial anisotropy , exchange , DMI energy for all
    ! active continuum FD nodes, then write them to ms_energy.<simid>.out.
    subroutine calc_ms_energies(emomM, mstep, nstep, simid, file_unit)  !, file_unit_aniso)
        use HamiltonianData, only : ham
        use InputData,       only : ham_inp

        real(dblprec), intent(in) :: emomM(3, msd%Natom, *)
        integer,       intent(in) :: mstep
        integer,       intent(in) :: nstep
        character(len=8), intent(in) :: simid
        integer, intent(in) :: file_unit
        ! integer, intent(in) :: file_unit_aniso
        ! integer :: cnt_moment
        integer       :: i, j, k, idx
        ! integer       :: file_unit
        real(dblprec) :: costh
        real(dblprec) :: E_demag_loc, E_ani_loc, E_xc_loc, E_dm_loc, E_xc_pen_loc
        real(dblprec) :: beff_xc(3), beff_dm(3)
        real(dblprec) :: diff_x, diff_y, diff_z
        integer       :: jj, idx_j
        character(len=30) :: filn

        if (.not. msd%output_enabled) return

        E_demag_loc = 0.0_dblprec
        E_ani_loc   = 0.0_dblprec
        E_xc_loc    = 0.0_dblprec
        E_xc_pen_loc= 0.0_dblprec
        E_dm_loc    = 0.0_dblprec

        ! write(*,*) 'DEBUG: nstep: ',nstep
        ! Demagnetization energy: E_demag = -(mu0/2) * sum M.H * dV
        !$omp parallel do collapse(3) default(shared) private(i, j, k) reduction(-:E_demag_loc)
        do k = 1, msd%nzc
            do j = 1, msd%nyc
                do i = 1, msd%nxc
                    if (msd%active_cell(i, j, k) /= 1) cycle
                    E_demag_loc = E_demag_loc &
                        - 0.5_dblprec * mu0&
                        * ( msd%Mgrid(1,i,j,k)*msd%Hgrid(1,i,j,k) &
                          + msd%Mgrid(2,i,j,k)*msd%Hgrid(2,i,j,k) &
                          + msd%Mgrid(3,i,j,k)*msd%Hgrid(3,i,j,k)) &
                        * msd%volume
                end do
            end do
        end do
        !$omp end parallel do

        ! cnt_moment  = 0
        ! Uniaxial anisotropy energy: E_ani = -K1*cos^2(th) - K2*cos^4(th)
        ! uses ham%taniso, ham%eaniso, ham%kaniso which are already populated
        ! only processes nodes marked taniso==1 (uniaxial).
        if (ham_inp%do_anisotropy == 1) then
            ! write(*,*) 'DEBUG: Calculating uniaxial anisotropy energy'
            !$omp parallel do collapse(3) default(shared) private(i, j, k, costh, idx) reduction(+:E_ani_loc)
            do k = 1, msd%nzn
                do j = 1, msd%nyn
                    do i = 1, msd%nxn
                        idx = msd%node_index(i, j, k)
                        if (idx <= 0) cycle
                        if (ham%taniso(idx) /= 1) cycle
                        costh = ham%eaniso(1,idx)*emomM(1,idx,1) &
                              + ham%eaniso(2,idx)*emomM(2,idx,1) &
                              + ham%eaniso(3,idx)*emomM(3,idx,1)
                        E_ani_loc = E_ani_loc + (ham%kaniso(1,idx) * costh**2 &   ! 2.179872325d-21 
                            + ham%kaniso(2,idx) * costh**4)

                        ! if (mstep == nstep) then
                        !     cnt_moment = cnt_moment + 1
                        !     write(file_unit_aniso,'(i6,2x,es18.9)') cnt_moment, mub * (ham%kaniso(1,idx) * costh**2 &
                        !                                                                 + ham%kaniso(2,idx) * costh**4)
                        !     flush(file_unit_aniso)
                        ! end if
                    end do
                end do
            end do
            !$omp end parallel do
        end if

        !Exchange energy
        if (ham_inp%do_dm /= 0 .or. ham%max_no_neigh > 0) then
            !$omp parallel do default(shared) private(i, j, k, idx, jj, idx_j, beff_xc, diff_x, diff_y, diff_z) reduction(-:E_xc_loc) reduction(+:E_xc_pen_loc)
            do k = 1, msd%nzn
                do j = 1, msd%nyn
                    do i = 1, msd%nxn
                        idx = msd%node_index(i, j, k)
                        if (idx <= 0) cycle
                        beff_xc = 0.0_dblprec
                        ! write(*,*) 'DEBUG: size(ham%ncoup(jj, ham%aHam(idx), 1))   ', size(ham%ncoup,1), size(ham%ncoup,2), size(ham%ncoup,3)
                        do jj = 1, ham%nlistsize(ham%aHam(idx))
                            idx_j = ham%nlist(jj, idx)
                            beff_xc = beff_xc + ham%ncoup(jj, ham%aHam(idx), 1) * emomM(1:3, idx_j, 1)
                            ! Calculate penalty gradient exchange
                            ! E_pen = 1/4 * sum_{i,j} J_ij * (m_i - m_j)^2
                            diff_x = emomM(1,idx,1) - emomM(1,idx_j,1)
                            diff_y = emomM(2,idx,1) - emomM(2,idx_j,1)
                            diff_z = emomM(3,idx,1) - emomM(3,idx_j,1)
                            E_xc_pen_loc = E_xc_pen_loc + 0.25_dblprec * ham%ncoup(jj, ham%aHam(idx), 1)* (diff_x**2 + diff_y**2 + diff_z**2)
                        end do
                        E_xc_loc = E_xc_loc - 0.5_dblprec * (emomM(1,idx,1)*beff_xc(1) &
                                             + emomM(2,idx,1)*beff_xc(2) &
                                             + emomM(3,idx,1)*beff_xc(3))
                    end do
                end do
            end do
            !$omp end parallel do
        end if

        !DM energy
        if (ham_inp%do_dm == 1) then
            !$omp parallel do default(shared) private(i, j, k, idx, jj, beff_dm) reduction(-:E_dm_loc)
            do k = 1, msd%nzn
                do j = 1, msd%nyn
                    do i = 1, msd%nxn
                        idx = msd%node_index(i, j, k)
                        if (idx <= 0) cycle
                        beff_dm = 0.0_dblprec
                        do jj = 1, ham%dmlistsize(ham%aHam(idx))
                            beff_dm(1) = beff_dm(1) + ham%dm_vect(3,jj,ham%aHam(idx))*emomM(2,ham%dmlist(jj,idx),1) &
                                - ham%dm_vect(2,jj,ham%aHam(idx))*emomM(3,ham%dmlist(jj,idx),1)
                            beff_dm(2) = beff_dm(2) + ham%dm_vect(1,jj,ham%aHam(idx))*emomM(3,ham%dmlist(jj,idx),1) &
                                - ham%dm_vect(3,jj,ham%aHam(idx))*emomM(1,ham%dmlist(jj,idx),1)
                            beff_dm(3) = beff_dm(3) + ham%dm_vect(2,jj,ham%aHam(idx))*emomM(1,ham%dmlist(jj,idx),1) &
                                - ham%dm_vect(1,jj,ham%aHam(idx))*emomM(2,ham%dmlist(jj,idx),1)
                        end do
                        E_dm_loc = E_dm_loc - 0.5_dblprec * (emomM(1,idx,1)*beff_dm(1) &
                                             + emomM(2,idx,1)*beff_dm(2) &
                                             + emomM(3,idx,1)*beff_dm(3))
                    end do
                end do
            end do
            !$omp end parallel do
        end if

        msd%E_demag = E_demag_loc * (mub**2 / 1.0d-30)
        msd%E_ani   = E_ani_loc * mub 
        msd%E_xc    = E_xc_loc * mub 
        msd%E_xc_pen = E_xc_pen_loc * mub
        msd%E_dm    = E_dm_loc * mub
        msd%E_total = msd%E_demag + msd%E_ani + msd%E_xc_pen + msd%E_dm

        ! Write to ms_energy.<simid>.out
        ! Columns: Iter  E_total  E_demag  E_ani  E_xc  E_xc_pen  E_dm
        write(file_unit,'(i6,7(2x,es18.9))') &
            mstep, msd%E_total, msd%E_demag, msd%E_ani, msd%E_xc, msd%E_xc_pen, msd%E_dm, msd%E_xc_pen + msd%E_dm

    end subroutine calc_ms_energies

    subroutine cleanup_multiscale_demag()
        ! if (allocated(msd%active_node)) deallocate(msd%active_node)
        if (allocated(msd%node_index))  deallocate(msd%node_index)
        if (allocated(msd%active_cell)) deallocate(msd%active_cell)
        if (allocated(msd%Mgrid)) deallocate(msd%Mgrid)
        if (allocated(msd%Hgrid)) deallocate(msd%Hgrid)
        ! if (allocated(msd%Bdemag))  deallocate(msd%Bdemag)
        if (allocated(msd%Ntensor)) deallocate(msd%Ntensor)
    end subroutine cleanup_multiscale_demag


    ! Newell tensor precursor functions and helpers
    ! F(x,y,z) - precursor for diagonal Newell components N_xx, N_yy, N_zz.
    pure function newell_F(x, y, z) result(res)
        real(dblprec), intent(in) :: x, y, z
        real(dblprec) :: res, R
        R = sqrt(x*x + y*y + z*z)
        if (R < 1.0d-12) then
            res = 0.0_dblprec
            return
        end if
        res =(1.0_dblprec/6.0_dblprec)*(2.0_dblprec*x*x - y*y - z*z)*R  &
             + 0.5_dblprec*y*(z*z - x*x)*log_safe(y + R)                    &
             + 0.5_dblprec*z*(y*y - x*x)*log_safe(z + R)                    &
             - x*y*z*atan_safe(y*z, x*R)
    end function newell_F

    ! G(x,y,z) - precursor for off-diagonal Newell components N_xy, N_xz, N_yz.
    pure function newell_G(x, y, z) result(res)
        real(dblprec), intent(in) :: x, y, z
        real(dblprec) :: res, R
        R = sqrt(x*x + y*y + z*z)
        if (R < 1.0d-12) then
            res = 0.0_dblprec
            return
        end if
        res = -(1.0_dblprec/3.0_dblprec)*x*y*R                                        &
             + x*y*z*log_safe(z + R)                                                  &
             + (1.0_dblprec/6.0_dblprec)*y*(3.0_dblprec*z*z - y*y)*log_safe(x + R)    &
             + (1.0_dblprec/6.0_dblprec)*x*(3.0_dblprec*z*z - x*x)*log_safe(y + R)    &
             - (1.0_dblprec/6.0_dblprec)*z*z*z*atan_safe(x*y, z*R)                    &
             - 0.5_dblprec*y*y*z*atan_safe(x*z, y*R)                                  &
             - 0.5_dblprec*x*x*z*atan_safe(y*z, x*R)
    end function newell_G

    ! 27-point operator L applied to F or G.
    ! Converts the precursor function into a Newell tensor element via
    ! N = (1/4*pi*hx*hy*hz) * L{Phi}
    ! L{Phi} = sum_{e1,e2,e3 in {-1,0,1}} [8 / (-2)^(|e1|+|e2|+|e3|)] * Phi(r + e.h)
    function apply_L(is_F, x, y, z, hx, hy, hz) result(res)
        logical,       intent(in) :: is_F
        real(dblprec), intent(in) :: x, y, z, hx, hy, hz
        real(dblprec) :: res, w
        integer :: e1, e2, e3
        ! real(dblprec), parameter :: pi = 3.14159265358979323846_dblprec 

        res = 0.0_dblprec
        do e1 = -1, 1
            do e2 = -1, 1
                do e3 = -1, 1
                    w = 8.0_dblprec / ((-2.0_dblprec)**real(abs(e1)+abs(e2)+abs(e3), dblprec))
                    if (is_F) then
                        res = res + w * newell_F(x + e1*hx, y + e2*hy, z + e3*hz)
                    else
                        res = res + w * newell_G(x + e1*hx, y + e2*hy, z + e3*hz)
                    end if
                end do
            end do
        end do
        res = res / (4.0_dblprec * pi * hx * hy * hz)
    end function apply_L

    pure function log_safe(v) result(r)
        real(dblprec), intent(in) :: v
        real(dblprec) :: r
        if (v <= 0.0_dblprec) then
            !returns a massive negative number,which is closer to the limit of log(x) as x approaches 0.
            r = log(tiny(v)) 
        else
            r = log(v)
        end if
    end function log_safe

    ! Mathematically identical to atan(num / den) but avoids division by zero.
    pure function atan_safe(num, den) result(res)
        real(dblprec), intent(in) :: num, den
        real(dblprec) :: res
        real(dblprec), parameter :: pi_2 = 1.5707963267948966_dblprec

        if (abs(den) < 1.0d-12) then
            if (num > 1.0d-12 .and. den > 0.0_dblprec) then
                res = pi_2
            else if (num < -1.0d-12 .and. den > 0.0_dblprec) then
                res = -pi_2
            else if (num > 1.0d-12 .and. den < 0.0_dblprec) then
                res = -pi_2
            else if (num < -1.0d-12 .and. den < 0.0_dblprec) then
                res = pi_2
            else
                res = 0.0_dblprec
            end if
        else
            res = atan(num / den)
        end if
    end function atan_safe

end module MultiscaleDemag