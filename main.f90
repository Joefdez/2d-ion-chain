program twoDChain

  use mpi
  use constants
  use randomModule
  use initcond_generator
  use initialization
  use support_functions

  implicit none
  !-------------------------------------------------------------------------------------------------------------------------------------------------------
  !                                           VARIABLES BLOCK
  !-------------------------------------------------------------------------------------------------------------------------------------------------------

  !real(kind=8), dimension(:), allocatable ::  xx, yy !, zz
  real(kind=8), dimension(:,:), allocatable :: xx_av, yy_av !, zz_av
  real(kind=8), dimension(:,:), allocatable :: xx_avt, yy_avt !, zz_avt
  real(kind=8), dimension(:,:), allocatable :: kinEn_av, kinEn_avt
  real(kind=8), dimension(:,:), allocatable :: xPx_av, xPx_avt, yPy_av, yPy_avt
  real(kind=8), dimension(:,:), allocatable :: Px, Py !, Pz
  real(kind=8), dimension(:,:), allocatable :: Px_av, Py_av !, Pz_av
  real(kind=8), dimension(:),   allocatable :: YYi, AA, AAi
  real(kind=8), dimension(:,:), allocatable :: YY
  real(kind=8), dimension(:,:), allocatable :: Cf1, Cf2, invD1, invD2
  real(kind=8), dimension(:), allocatable   :: xx0, yy0, zz0, eqX, eqY!, eq
  real(kind=8), dimension(:), allocatable   :: kinEN_f, kinEN_ft, kinEn_f_av
  real(kind=8), dimension(:), allocatable   :: xx_f, xx_ft, yy_f, yy_ft, xx_f_av, yy_f_av
  ! Ion chain + Laser characteristics

  real(kind=8) :: omega_0    ! trageted transition frequency
  real(kind=8) :: Gam1       ! natural linewidth of targeted transition
  real(kind=8) :: long_freq, alpha
  real(kind=8) :: temp1 , temp2       ! bath temperatures
  real(kind=8) :: eta1, eta2          ! to be expanded on if more friction coefficients are needed
  real(kind=8) :: D1, D2              ! Diffusion coefficients
  real(kind=8) :: del1, del2          ! detuning parameters
  real(kind=8) :: I1, I2              ! normalized intensities
  real(kind=8) :: omega1, omega2, k1, k2
  real(kind=8) :: char_length             ! characteristic length
  real(kind=8) :: mass1, charge1
  ! Adimensional variables

  real(kind=8) :: aeta1, aeta2        ! adimensional friction coefficients
  real(kind=8) :: aD1, aD2            ! adimensional Diffusion coefficients


  real(kind=8) :: tt, dt, dst, st
  integer :: n_particles, traj, local_traj, nsteps, rem, save_freq, n_ssteps, p_traj
  integer ::n_elems
  real(kind=8), dimension(:), allocatable   :: YYold, stoch_terms, dStoc
  integer                                   :: ii,jj, kk
  real(kind=8)                              :: ic_radius

  ! mpi variables
  integer :: rank, procs, status(MPI_STATUS_SIZE), alloc_err, source, ierr

  ! Neccesary mpi initialization calls
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, procs, ierr)

  !-------------------------------------------------------------------------------------------------------------------------------------------------------
  !                                                  INITIALIZATION BLOCK
  !-------------------------------------------------------------------------------------------------------------------------------------------------------
  ! read in the simulation information
  if(rank .eq. 0) then
    print*, "Reading solution parameters."
    call initialize_system(n_particles, mass1, charge1, tt, dt, traj, save_freq, long_freq, alpha, ic_radius)
  end if

  call mpi_bcast(n_particles, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(traj, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(nsteps, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(save_freq, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)

  call mpi_bcast(mass1, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(charge1, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(tt, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(dt, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(alpha, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(ic_radius, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(long_freq, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)

  call mpi_barrier(mpi_comm_world, ierr)


  mass1  = mass1*uu
  charge1 = charge1 * ee
  long_freq = 2.0d0*pi*long_freq
  nsteps=int(tt/dt)
  char_length = ((charge1*charge1/(4.0d0*pi*ep0))/(mass1*long_freq*long_freq))**(1.0/3.0)

  local_traj = traj/procs
  rem = mod(traj, procs)
  if (rank .lt. rem) local_traj = local_traj + 1
  n_ssteps = int(nsteps/save_freq)
  st = int(0.8*nsteps)
  dst = sqrt(dt)
  n_elems = n_particles*n_ssteps

  if(rank .eq. 0) then
    print*, "Reading laser parameters."
    call initialize_laser_chain(del1, del2, Gam1, omega_0, I1, I2)
  end if


  call mpi_bcast(del1, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(del2, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Gam1, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(omega_0, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(I1, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(I2, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)

  ! Calculate laser frequencies and wavelengths
  Gam1  = 2.0d0 * pi * Gam1
  del1  = del1 * Gam1
  del2  = del2 * Gam1
  omega_0 = 2.0d0 * pi * omega_0
  omega1 = del1 + omega_0
  omega2 = del2 + omega_0

  k1 =  omega1 / cc   ! Not really wavelengths, dimensionally they are wavenumbers
  k2 =  omega2 / cc
  ! Calculate diffusion and friction coefficients
  eta1 = -4.0d0*hbar*k1*k1*I1*(2.0d0*del1/Gam1)/( (1 + 4.0d0*del1*del1/(Gam1*Gam1)) * (1 + 4.0d0*del1*del1/(Gam1*Gam1)) )
  eta2 = -4.0d0*hbar*k2*k2*I2*(2.0d0*del2/Gam1)/( (1 + 4.0d0*del2*del2/(Gam1*Gam1)) * (1 + 4.0d0*del2*del2/(Gam1*Gam1)) )
  D1   = hbar*hbar*k1*k1*I1*(Gam1)/(1.0d0 + 4.0d0*del1*del1/(Gam1*Gam1))
  D2   = hbar*hbar*k2*k2*I2*(Gam1)/(1.0d0 + 4.0d0*del2*del2/(Gam1*Gam1))

  !tt = long_freq * mass1 / eta1
  !dt = tt / 20000.d0
  !nsteps=int(tt/dt)
  !n_ssteps = int(nsteps/save_freq)


  call dimensionless_doppler_values(eta1, D1, mass1, long_freq, char_length, aeta1, aD1)
  call dimensionless_doppler_values(eta2, D2, mass1, long_freq, char_length, aeta2, aD2)

  allocate(xx_av(1:n_particles, n_ssteps))
  xx_av = 0.0d0
  allocate(yy_av(1:n_particles, n_ssteps))
  yy_av = 0.0d0
  allocate(xx_avt(1:n_particles, n_ssteps))
  xx_avt = 0.0d0
  allocate(yy_avt(1:n_particles, n_ssteps))
  yy_avt = 0.0d0
  allocate(kinEN_av(1:n_particles, n_ssteps))
  kinEn_av = 0.0d0
  allocate(kinEN_avt(1:n_particles, n_ssteps))
  kinEn_avt = 0.0d0
  allocate(xx_f(1:n_particles))
  xx_f = 0.0d0
  allocate(xx_f_av(1:n_particles))
  xx_f_av = 0.0d0
  allocate(yy_f(1:n_particles))
  yy_f = 0.0d0
  allocate(yy_f_av(1:n_particles))
  yy_f_av = 0.0d0
  allocate(xx_ft(1:n_particles))
  xx_ft = 0.0d0
  allocate(yy_ft(1:n_particles))
  yy_ft = 0.0d0
  allocate(kinEn_f(1:n_particles))
  kinEn_f = 0.0d0
  allocate(kinEn_f_av(1:n_particles))
  kinEn_f_av = 0.0d0
  allocate(kinEn_ft(1:n_particles))
  kinEN_ft = 0.0d0
  allocate(xPx_av(1:n_particles, n_ssteps))
  xPx_av = 0.0d0
  allocate(yPy_av(1:n_particles, n_ssteps))
  yPy_av = 0.0d0
  allocate(xPx_avt(1:n_particles, n_ssteps))
  xPx_av = 0.0d0
  allocate(yPy_avt(1:n_particles, n_ssteps))
  yPy_av = 0.0d0
  allocate(YY(1:4*n_particles, nsteps))
  YY = 0.0d0
  allocate(Cf1(2*n_particles,2*n_particles))
  Cf1 = 0.0d0
  allocate(Cf2(2*n_particles,2*n_particles))
  Cf2 = 0.0d0
  allocate(invD1(2*n_particles,n_particles))
  invD1 = 0.0d0
  allocate(invD2(2*n_particles,n_particles))
  invD2 = 0.0d0

  allocate(eqX(1:n_particles))
  eqX = 0.0d0
  allocate(eqY(1:n_particles))
  eqY = 0.0d0
  allocate(xx0(1:n_particles))
  xx0 = 0.0d0
  allocate(yy0(1:n_particles))
  yy0 = 0.0d0

  allocate(AA(1:4*n_particles))
  AA = 0.0d0
  allocate(AAi(1:4*n_particles))
  AAi = 0.0d0
  allocate(YYi(1:4*n_particles))
  YYi = 0.0d0
  allocate(dStoc(1:4*n_particles))
  dStoc = 0.0d0
  allocate(stoch_terms(1:4*n_particles))
  stoch_terms = 0.0d0
  stoch_terms((2*n_particles+1):(2*n_particles+3)) = sqrt(2.0d0 * aD1)
  stoch_terms((3*n_particles-2):(3*n_particles))   = sqrt(2.0d0 * aD2)
  stoch_terms((3*n_particles+1):(3*n_particles+3)) = sqrt(2.0d0 * aD1)
  stoch_terms((4*n_particles-2):(4*n_particles))   = sqrt(2.0d0 * aD2)

  open(unit=15, file="eq_pos.dat", action='read')
  do ii = 1, n_particles, 1
    read(15,*) eqX(ii), eqY(ii)
  end do


!-------------------------------------------------------------------------------------------------------------------------------------------------------
!                                                  Resolution block
!-------------------------------------------------------------------------------------------------------------------------------------------------------



  do ii=1, local_traj, 1
    print*,"Process", rank, "on stochastic trajectory ", ii
    YY = 0.0d0
    call icpgen(n_particles, ic_radius, eqX, eqY, xx0(:), yy0(:))
    YY(1:n_particles,1) =  xx0(:)
    YY(n_particles+1:2*n_particles,1) = yy0(:)
    print*, "Proc. ", rank, " on trajectory", ii
    do jj=1, nsteps-1,1
      call stoch_vector(dst, n_particles,  stoch_terms, dStoc)
      call Cforce(YY(1:n_particles,jj), YY(n_particles+1:2*n_particles,jj), n_particles, invD1, Cf1)
      call ddt(YY(:,jj), AA, n_particles, aeta1, aeta2, alpha, Cf1)
      YYi = YY(:,jj) + AA*dt + dStoc
      call Cforce(YYi(1:n_particles), YYi(n_particles+1:2*n_particles), n_particles, invD2, Cf2)
      call ddt(YYi, AAi, n_particles, aeta1, aeta2, alpha, Cf2)
      YY(:,jj+1) = YY(:,jj) + 0.5d0*(AA+AAi)*dt + dStoc
      if(jj .gt. st) then
        xx_f = xx_f + YY(1:n_particles,jj)
        yy_f = yy_f + YY(n_particles+1:2*n_particles,jj)
        kinEN_f = kinEn_f +  0.5d0*YY((2*n_particles+1):3*n_particles,jj)*YY((2*n_particles+1):3*n_particles,jj) +&
                  0.5d0*YY((3*n_particles+1):4*n_particles,jj)*YY((3*n_particles+1):(4*n_particles),jj)
      end if
    end do
    kinEn_av = kinEn_av + 0.5d0*YY((2*n_particles+1):3*n_particles,1::save_freq)*YY((2*n_particles+1):3*n_particles,1::save_freq) +&
              0.5d0*YY((3*n_particles+1):4*n_particles,1::save_freq)*YY((3*n_particles+1):4*n_particles,1::save_freq)
    xx_av   = xx_av + YY(1:n_particles,1::save_freq)
    yy_av   = yy_av + YY(n_particles+1:2*n_particles,1::save_freq)
    xPx_av  = xPx_av + YY(1:n_particles,1::save_freq)*YY((2*n_particles+1):3*n_particles,1::save_freq)
    yPy_av  = yPy_av + YY(n_particles+1:2*n_particles,1::save_freq)*YY((3*n_particles+1):4*n_particles,1::save_freq)
    kinEN_f_av = kinEN_f_av + kinEn_f/(nsteps-st)
    xx_f_av    = xx_f_av + xx_f/(nsteps-st)
    yy_f_av    = yy_f_av + yy_f/(nsteps-st)
    if( ( mod(ii,5) .eq. 0) .and. (ii .lt. int(traj/procs) ) ) then
      call mpi_reduce(kinEn_av, kinEn_avt , n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(xx_av, xx_avt , n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(yy_av, yy_avt , n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(xPx_av, xPx_avt , n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(yPy_av, yPy_avt , n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(kinEn_f, kinEn_ft, n_particles, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(xx_f_av, xx_ft, n_particles, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(yy_f_av, yy_ft, n_particles, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(ii, p_traj, 1, mpi_integer, mpi_sum, 0, mpi_comm_world, ierr)
      if(rank .eq. 0) then
        open(unit=11, file="results/posX.dat")
        open(unit=12, file="results/posy.dat")
        open(unit=13, file="results/kinEN.dat")
        open(unit=14, file="results/xPx.dat")
        open(unit=15, file="results/yPy.dat")
        open(unit=16, file="results/avOvertime.dat")
        do jj=1, n_particles, 1
          write(11,*) xx_avt(jj,:)/p_traj
          write(12,*) yy_avt(jj,:)/p_traj
          write(13,*) kinEn_avt(jj,:)/p_traj
          write(14,*) xPx_avt(jj,:)/p_traj
          write(15,*) yPy_avt(jj,:)/p_traj
        end do
        write(16,*) xx_ft/p_traj
        write(16,*) yy_ft/p_traj
        write(16,*) kinEn_ft/p_traj
        close(unit=11)
        close(unit=12)
        close(unit=13)
        close(unit=14)
        close(unit=15)
        close(unit=16)
      end if
    end if
  end do

  call mpi_barrier(mpi_comm_world, ierr)
  print*, "Out of loop"
  call mpi_reduce(kinEn_av, kinEn_avt , n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(xx_av, xx_avt , n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(yy_av, yy_avt , n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(xPx_av, xPx_avt , n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(yPy_av, yPy_avt , n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(kinEn_f_av, kinEn_ft, n_particles, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(xx_f_av, xx_ft, n_particles, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(yy_f_av, yy_ft, n_particles, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)

  print*,"Everything written"
  if(rank .eq. 0) then
    print*, "Writing final results."
    open(unit=11, file="results/posX.dat")
    open(unit=12, file="results/posy.dat")
    open(unit=13, file="results/kinEN.dat")
    open(unit=14, file="results/xPx.dat")
    open(unit=15, file="results/yPy.dat")
    open(unit=16, file="results/avOvertime.dat")
    do jj=1, n_particles, 1
      write(11,*) xx_avt(jj,:)/traj
      write(12,*) yy_avt(jj,:)/traj
      write(13,*) kinEn_avt(jj,:)/traj
      write(14,*) xPx_avt(jj,:)/traj
      write(15,*) yPy_avt(jj,:)/traj
    end do
    write(16,*) xx_ft/traj
    write(16,*) yy_ft/traj
    write(16,*) kinEn_ft/traj
    close(unit=11)
    close(unit=12)
    close(unit=13)
    close(unit=14)
    close(unit=15)
    close(unit=16)
  end if

  call mpi_finalize(ierr)

end program twoDChain
