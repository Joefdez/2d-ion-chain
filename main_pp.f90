program twoDChain

  use mpi
  use constants
  use initcond_generator
  use initialization
  use randomModule
  use support_functions_twod

  implicit none

  ! variables
  real(kind=8), dimension(:), allocatable               :: xx0, yy0, px0, py0
  real(kind=8), dimension(:), allocatable               :: xxold, yyold, ppxold, ppyold
  real(kind=8), dimension(:), allocatable               :: xxnew, yynew, ppxnew, ppynew
  real(kind=8), dimension(:,:), allocatable             :: xxs, yys, ppxs, ppys
  real(kind=8), dimension(:,:), allocatable             :: xx2s, yy2s, ppx2s, ppy2s
  real(kind=8), dimension(:,:), allocatable             :: xpxs, ypys
  real(kind=8), dimension(:,:), allocatable             :: xx_avt, yy_avt, ppx_avt, ppy_avt
  real(kind=8), dimension(:,:), allocatable             :: xx2_avt, yy2_avt, ppx2_avt, ppy2_avt
  real(kind=8), dimension(:,:), allocatable             :: xpx_avt, ypy_avt
  real(kind=8), dimension(:,:), allocatable             :: xx_av, yy_av, ppx_av, ppy_av
  real(kind=8), dimension(:,:), allocatable             :: xx2_av, yy2_av, ppx2_av, ppy2_av
  real(kind=8), dimension(:,:), allocatable             :: xPx_av, yPy_av
  real(kind=8), dimension(:), allocatable               :: xxi, yyi, ppxi, ppyi
  real(kind=8), dimension(:,:), allocatable             :: fx1, fy1, fx2, fy2
  real(kind=8), dimension(:), allocatable               :: fx, fy
  real(kind=8), dimension(:), allocatable               :: Axx, Axxi, Ayy, Ayyi
  real(kind=8), dimension(:), allocatable               :: Apx, Apxi, Apy, Apyi
  real(kind=8), dimension(:), allocatable               :: dOmx, dOmy, dOmxc, dOmyc
  real(kind=8), dimension(:), allocatable               :: stermsBx, stermsCx, stermsBy, stermsCy

  real(kind=8)                                          :: tt, dt, dst, mass, charge, dist
  real(kind=8)                                          :: alpha, char_length, long_freq
  real(kind=8)                                          :: del1, del2, delC
  real(kind=8)                                          :: omega0, omega1, omega2, omegaC
  real(kind=8)                                          :: I1, I2, IC
  real(kind=8)                                          :: k1, k2, kC
  real(kind=8)                                          :: ic_radius
  real(kind=8)                                          :: Gam
  real(kind=8)                                          :: eta1, aeta1, eta2, aeta2
  real(kind=8)                                          :: D1, aD1, D2, aD2
  real(kind=8)                                          :: etaC, DC
  real(kind=8)                                          :: aetaC, aDC
  integer                                               :: nsteps, savefreq, nssteps, nparticles, nbath, n_elems
  integer                                               :: traj, local_traj, save_freq, rem
  integer                                               :: ii,jj, kk, ll
  real(kind=8)                                          :: seconds, seconds1
  ! mpi variables
  integer :: rank, procs, status(MPI_STATUS_SIZE), alloc_err, source, ierr
  call MPI_INIT(ierr)                                                                               ! Neccesary mpi initialization calls
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, procs, ierr)

  if(rank .eq. 0 ) then
    print*, "Reading solution and system parameters"
    call initialize_system(nparticles, mass, charge, tt, dt, traj, save_freq, long_freq, alpha, ic_radius)
  end if

  call mpi_bcast(nparticles, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(traj, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(nsteps, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(save_freq, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(mass, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(charge, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(tt, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(dt, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(alpha, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(ic_radius, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(long_freq, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)

  call mpi_barrier(mpi_comm_world, ierr)

  mass   = mass*uu
  charge = charge*ee
  long_freq = 2.0d0*pi*long_freq
  nsteps=int(tt/dt)
  char_length = ((charge*charge/(4.0d0*pi*ep0))/(mass*long_freq*long_freq))**(1.0/3.0)
  nssteps = int(nsteps/save_freq)  ! Not saving every single timestep saves memory. Must ask about this

  if(rank .eq. 0) then
    print*, "Reading laser parameters"
    call initialize_laser_chain(del1, del2, delC, Gam, omega0, I1, I2, IC)
  end if

  call mpi_bcast(del1, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(del2, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(delC, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Gam, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(omega0, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(I1, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(I2, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(IC, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)

  Gam  = 2.0d0 * pi * Gam
  del1  = del1 * Gam
  del2  = del2 * Gam
  delC  = delC * Gam
  omega0 = 2.0d0 * pi * omega0
  omega1 = del1 + omega0
  omega2 = del2 + omega0
  omegaC = delC + omega0

  k1 =  omega1 / cc   ! Not really wavelengths, dimensionally they are wavenumbers
  k2 =  omega2 / cc
  kC = omegaC  / cc
 ! Calculate diffusion and friction coefficients
  eta1 = -4.0d0*hbar*k1*k1*I1*(2.0d0*del1/Gam)/( (1 + 4.0d0*del1*del1/(Gam*Gam)) * (1 + 4.0d0*del1*del1/(Gam*Gam)) )
  eta2 = -4.0d0*hbar*k2*k2*I2*(2.0d0*del2/Gam)/( (1 + 4.0d0*del2*del2/(Gam*Gam)) * (1 + 4.0d0*del2*del2/(Gam*Gam)) )
  D1   = hbar*hbar*k1*k1*I1*(Gam)/(1.0d0 + 4.0d0*del1*del1/(Gam*Gam))
  D2   = hbar*hbar*k2*k2*I2*(Gam)/(1.0d0 + 4.0d0*del2*del2/(Gam*Gam))
  etaC = -4.0d0*hbar*kC*kC*IC*(2.0d0*delC/Gam)/( (1 + 4.0d0*delC*delC/(Gam*Gam)) * (1 + 4.0d0*delC*delC/(Gam*Gam)) )
  DC   = hbar*hbar*kC*kC*IC*(Gam)/(1.0d0 + 4.0d0*delC*delC/(Gam*Gam))

  ! Calculate dimensionless Doppler cooling parameters
  call dimensionless_doppler_values(eta1, D1, mass, long_freq, char_length, aeta1, aD1)
  call dimensionless_doppler_values(eta2, D2, mass, long_freq, char_length, aeta2, aD2)
  call dimensionless_doppler_values(etaC, DC, mass, long_freq, char_length, aetaC, aDc)

  local_traj = traj/procs
  rem = mod(traj, procs)
  if (rank .lt. rem) local_traj = local_traj + 1
  n_elems = nssteps*nparticles
  nbath = 3

  allocate(xx0(1:nparticles))
  xx0 = 0.0d0
  allocate(yy0(1:nparticles))
  yy0 = 0.0d0
  allocate(px0(1:nparticles))
  px0 = 0.0d0
  allocate(py0(1:nparticles))
  py0 = 0.0d0
  allocate(xxold(1:nparticles))
  xxold = 0.0d0
  allocate(yyold(1:nparticles))
  yyold = 0.0d0
  allocate(ppxold(1:nparticles))
  ppxold = 0.0d0
  allocate(ppyold(1:nparticles))
  ppyold = 0.0d0
  allocate(xxnew(1:nparticles))
  xxnew = 0.0d0
  allocate(yynew(1:nparticles))
  yynew = 0.0d0
  allocate(ppxnew(1:nparticles))
  ppxnew = 0.0d0
  allocate(ppynew(1:nparticles))
  ppynew = 0.0d0
  allocate(xxs(1:nparticles, 1:nssteps))
  xxs = 0.0d0
  allocate(yys(1:nparticles, 1:nssteps))
  yys = 0.0d0
  allocate(ppxs(1:nparticles, 1:nssteps))
  ppxs = 0.0d0
  allocate(ppys(1:nparticles, 1:nssteps))
  ppys = 0.0d0
  allocate(xx2s(1:nparticles, 1:nssteps))
  xx2s = 0.0d0
  allocate(yy2s(1:nparticles, 1:nssteps))
  yy2s = 0.0d0
  allocate(ppx2s(1:nparticles, 1:nssteps))
  ppx2s = 0.0d0
  allocate(ppy2s(1:nparticles, 1:nssteps))
  ppy2s = 0.0d0
  allocate(xpxs(1:nparticles, 1:nssteps))
  xpxs = 0.0d0
  allocate(ypys(1:nparticles, 1:nssteps))
  ypys = 0.0d0
  allocate(xxi(1:nparticles))
  xxi = 0.0d0
  allocate(yyi(1:nparticles))
  yyi = 0.0d0
  allocate(ppxi(1:nparticles))
  ppxi = 0.0d0
  allocate(ppyi(1:nparticles))
  ppyi = 0.0d0
  allocate(fx1(1:nparticles,1:nparticles))
  fx1 = 0.0d0
  allocate(fx2(1:nparticles,1:nparticles))
  fx2 = 0.0d0
  allocate(fy1(1:nparticles,1:nparticles))
  fy1 = 0.0d0
  allocate(fy2(1:nparticles,1:nparticles))
  fy2 = 0.0d0
  allocate(fx(1:nparticles))
  fx = 0.0d0
  allocate(fy(1:nparticles))
  fy = 0.0d0
  allocate(Axx(1:nparticles))
  Axx = 0.0d0
  allocate(Ayy(1:nparticles))
  Ayy = 0.0d0
  allocate(Apx(1:nparticles))
  Apx = 0.0d0
  allocate(Apy(1:nparticles))
  Apy = 0.0d0
  allocate(Axxi(1:nparticles))
  Axxi = 0.0d0
  allocate(Ayyi(1:nparticles))
  Ayyi = 0.0d0
  allocate(Apxi(1:nparticles))
  Apxi = 0.0d0
  allocate(Apyi(1:nparticles))
  Apyi = 0.0d0
  allocate(dOmx(1:nparticles))
  dOmx = 0.0d0
  allocate(dOmy(1:nparticles))
  dOmy = 0.0d0
  allocate(dOmxc(1:nparticles))
  dOmxc = 0.0d0
  allocate(dOmyc(1:nparticles))
  dOmyc = 0.0d0
  allocate(stermsBx(1:nparticles))
  stermsBx = 0.0d0
  allocate(stermsBy(1:nparticles))
  stermsBy = 0.0d0
  allocate(stermsCx(1:nparticles))
  stermsCx = 0.0d0
  allocate(stermsCy(1:nparticles))
  stermsCy = 0.0d0

  allocate(xx_av(1:nparticles, 1:nssteps))
  xx_av = 0.0d0
  allocate(yy_av(1:nparticles, 1:nssteps))
  yy_av = 0.0d0
  allocate(ppx_av(1:nparticles, 1:nssteps))
  ppx_av = 0.0d0
  allocate(ppy_av(1:nparticles, 1:nssteps))
  ppy_av = 0.0d0
  allocate(xx2_av(1:nparticles, 1:nssteps))
  xx2_av = 0.0d0
  allocate(yy2_av(1:nparticles, 1:nssteps))
  yy2_av = 0.0d0
  allocate(ppx2_av(1:nparticles, 1:nssteps))
  ppx2_av = 0.0d0
  allocate(ppy2_av(1:nparticles, 1:nssteps))
  ppy2_av = 0.0d0
  allocate(xpx_av(1:nparticles, 1:nssteps))
  xpx_av = 0.0d0
  allocate(ypy_av(1:nparticles, 1:nssteps))
  ypy_av = 0.0d0

  allocate(xx_avt(1:nparticles, 1:nssteps))
  xx_avt = 0.0d0
  allocate(yy_avt(1:nparticles, 1:nssteps))
  yy_avt = 0.0d0
  allocate(ppx_avt(1:nparticles, 1:nssteps))
  ppx_avt = 0.0d0
  allocate(ppy_avt(1:nparticles, 1:nssteps))
  ppy_avt = 0.0d0
  allocate(xx2_avt(1:nparticles, 1:nssteps))
  xx2_avt = 0.0d0
  allocate(yy2_avt(1:nparticles, 1:nssteps))
  yy2_avt = 0.0d0
  allocate(ppx2_avt(1:nparticles, 1:nssteps))
  ppx2_avt = 0.0d0
  allocate(ppy2_avt(1:nparticles, 1:nssteps))
  ppy2_avt = 0.0d0
  allocate(xpx_avt(1:nparticles, 1:nssteps))
  xpx_avt = 0.0d0
  allocate(ypy_avt(1:nparticles, 1:nssteps))
  ypy_avt = 0.0d0

  stermsCx = sqrt(2.0d0*aDc)
  stermsCy = sqrt(2.0d0*aDc)
  stermsBx(1:nbath) = sqrt(2.0d0*aD1)
  stermsBy(1:nbath) = sqrt(2.0d0*aD1)
  stermsBx((nparticles-nbath+1):nparticles) = sqrt(2.0d0*aD2)
  stermsBy((nparticles-nbath+1):nparticles) = sqrt(2.0d0*aD2)
  if(rank .eq. 0) then
     seconds = mpi_wtime()
     seconds1 = mpi_wtime()
  end if

  if (rank .eq. 0) then
    open(unit=15, file="eq_pos.dat", action='read')
    do ii = 1, nparticles, 1
      read(15,*) xx0(ii), yy0(ii)
    end do
  end if

  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_bcast(xx0, nparticles, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(yy0, nparticles, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)

  do kk=1, local_traj, 1
    print*, "Proc.", rank, "on trajectory", kk
    call icpgen(nparticles, 0.02d0, xx0, yy0, xxold, yyold)
    call ranseed()
    xxs(1,:) = xxold
    yys(1,:) = yyold
    ll = 1
    do ii=1, nsteps-1, 1
      call coulombM(nparticles, xxold, yyold, fx1, fy1)
      fx = 0.0d0
      fy = 0.0d0
      do jj=1, nparticles, 1
        fx(jj) = sum(fx1(jj,:), 1)
        fy(jj) = sum(fy1(jj,:), 1)
      end do
      call vecA(xxold, yyold, ppxold, ppyold, fx, fy, alpha, aeta1, aeta2, aetaC, nbath, nparticles, Axx, Ayy, Apx, Apy)
      call vecB_edges(dst, nparticles, dOmx, dOmy)
      !call vecB_cool(dst, nparticles, dOmxc, dOmyc)
      xxi = xxold + Axx*dt
      yyi = yyold + Ayy*dt
      ppxi = ppxold + Apx*dt  + stermsBx*dOmx! + stermsCx*dOmxc
      ppyi = ppyold + Apy*dt  + stermsBy*dOmy! + stermsCy*dOmyc
      fx = 0.0d0
      fy = 0.0d0
      call coulombM(nparticles, xxi, yyi, fx2, fy2)
      do jj=1, nparticles, 1
        fx(jj) = sum(fx2(jj,:), 1)
        fy(jj) = sum(fy2(jj,:), 1)
      end do
      call vecA(xxi, yyi, ppxi, ppyi, fx, fy, alpha, aeta1, aeta2, aetaC, nbath, nparticles, Axxi, Ayyi, Apxi, Apyi)
      call vecB_edges(dst, nparticles, dOmx, dOmy)
      !call vecB_cool(dst, nparticles, dOmxc, dOmyc)
      xxnew   = xxold + 0.5d0*(Axx + Axxi)*dt
      yynew   = yyold + 0.5d0*(Ayy + Ayyi)*dt
      ppxnew  = ppxold + 0.5d0*(Apx + Apxi)*dt + stermsBx*dOmx !+ stermsCx*dOmxc
      ppynew  = ppyold + 0.5d0*(Apy + Apyi)*dt + stermsBy*dOmy !+ stermsCy*dOmyc
      if( mod(ii,save_freq) .eq. 0) then
        ll = ll + 1
        xxs(:,ll)   = xxnew
        yys(:,ll)   = yynew
        ppxs(:,ll)  = ppxnew
        ppys(:,ll)  = ppynew
        xx2s(:,ll)  = xxnew*xxnew
        yy2s(:,ll)  = yynew*yynew
        ppx2s(:,ll) = ppxnew*ppxnew
        ppy2s(:,ll) = ppynew*ppynew
        xpxs(:,ll)  = xxnew*ppxnew
        ypys(:,ll)  = yynew*ppynew
      end if
      xxold   = xxnew
      yyold   = yynew
      ppxold  = ppxnew
      ppyold  = ppynew
    end do
    xx_av  = (xx_av*(kk-1) + xxs)/kk
    yy_av  = (yy_av*(kk-1) + yys)/kk
    ppx_av = (ppx_av*(kk-1) + ppxs)/kk
    ppy_av = (ppy_av*(kk-1) + ppys)/kk
    xx2_av  = (xx2_av*(kk-1) + xx2s)/kk
    yy2_av  = (yy2_av*(kk-1) + yy2s)/kk
    ppx2_av = (ppx2_av*(kk-1) + ppx2s)/kk
    ppy2_av = (ppy2_av*(kk-1) + ppy2s)/kk
    xpx_av  = (xpx_av*(kk-1) + xpxs)/kk
    ypy_av  = (ypy_av*(kk-1) + ypys)/kk
    if( ( mod(kk,5) .eq. 0) .and. (kk .lt. int(traj/procs) ) ) then
     print*, "Writing PARTIAL results to files after ", kk, "trajectories."
     call mpi_reduce(xx_av*kk , xx_avt , n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(yy_av*kk , yy_avt , n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(xx2_av*kk, xx2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(yy2_av*kk, yy2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(ppx_av*kk, ppx_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(ppy_av*kk, ppy_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(ppx2_av*kk, ppx2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(ppy2_av*kk, ppy2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(xpx_av*kk, xpx_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(ypy_av*kk, ypy_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(kk, traj, 1, mpi_integer, mpi_sum, 0, mpi_comm_world, ierr)

      if(rank .eq. 0) then
       xx_avt   = xx_avt*char_length/traj
       yy_avt   = yy_avt*char_length/traj
       xx2_avt  = xx2_avt*char_length*char_length/traj
       yy2_avt  = yy2_avt*char_length*char_length/traj
       ppx_avt  = ppx_avt*char_length*mass*long_freq/traj
       ppy_avt  = ppy_avt*char_length*mass*long_freq/traj
       ppx2_avt = ppx2_avt*char_length*char_length*mass*long_freq*long_freq/(2.0d0*kb)/traj ! Convert to temperature in mK
       ppy2_avt = ppy2_avt*char_length*char_length*mass*long_freq*long_freq/(2.0d0*kb)/traj ! Convert to temperature in mK
       xpx_avt  = xpx_avt*char_length*char_length*mass*long_freq/traj
       ypy_avt  = ypy_avt*char_length*char_length*mass*long_freq/traj
       open(unit=11, file="posX.dat")
       open(unit=12, file="posY.dat")
       open(unit=13, file="temperatures.dat")
       do jj=1, nparticles
        write(11,*) xx_avt(jj,:)
        write(12,*) yy_avt(jj,:)
        write(13,*) ppx2_avt(jj,:) + ppy2_avt(jj,:)
       end do
       close(unit=11)
       close(unit=12)
       close(unit=13)
       xx_avt   = 0.0d0
       yy_avt   = 0.0d0
       xx2_avt  = 0.0d0
       yy2_avt  = 0.0d0
       ppx_avt  = 0.0d0
       ppy_avt  = 0.0d0
       ppx2_avt = 0.0d0
       ppy2_avt = 0.0d0
       xpx_avt  = 0.0d0
       ypy_avt  = 0.000
      end if
    end if
  end do
  print*,"Proc ", rank, " finished integrating"
  call mpi_barrier(mpi_comm_world, ierr)
  print*, "sending"
  call mpi_reduce(xx_av*kk , xx_avt , n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  print*, "sending"
  call mpi_reduce(yy_av*kk , yy_avt , n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  print*, "sending"
  call mpi_reduce(xx2_av*kk, xx2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(yy2_av*kk, yy2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(ppx_av*kk, ppx_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(ppy_av*kk, ppy_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(ppx2_av*kk, ppx2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(ppy2_av*kk, ppy2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(xpx_av*kk, xpx_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(ypy_av*kk, ypy_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(local_traj, traj, 1, mpi_integer, mpi_sum, 0, mpi_comm_world, ierr)

  if(rank .eq. 0) then
    seconds = mpi_wtime() - seconds
    print*, "total integration + partial writing time:", seconds, seconds/3600.0d0
    !call date_and_time(values=time)
    !print*, "Integration completed at " ,time(5:)
    xx_avt   = xx_avt*char_length/traj
    yy_avt   = yy_avt*char_length/traj
    xx2_avt  = xx2_avt*char_length*char_length/traj
    yy2_avt  = yy2_avt*char_length*char_length/traj
    ppx_avt  = ppx_avt*char_length*mass*long_freq/traj
    ppy_avt  = ppy_avt*char_length*mass*long_freq/traj
    ppx2_avt = ppx2_avt*char_length*char_length*mass*long_freq*long_freq/(2.0d0*kb)/traj ! Convert to temperature in mK
    ppy2_avt = ppy2_avt*char_length*char_length*mass*long_freq*long_freq/(2.0d0*kb)/traj ! Convert to temperature in mK
    xpx_avt  = xpx_avt*char_length*char_length*mass*long_freq/traj
    ypy_avt  = ypy_avt*char_length*char_length*mass*long_freq/traj
    open(unit=11, file="posX.dat")
    open(unit=12, file="posY.dat")
    open(unit=13, file="temperatures.dat")
    do jj=1, nparticles
     write(11,*) xx_avt(jj,:)
     write(12,*) yy_avt(jj,:)
     write(13,*) ppx2_avt(jj,:) + ppy2_avt(jj,:)
    end do
    close(unit=11)
    close(unit=12)
    close(unit=13)
    print*, "Writing FINAL results to files after ", traj , "trajectories."
    seconds1 = mpi_wtime() - seconds1
    print*, "integration + message + write time:", seconds1

  end if
    ! back to physical units

  call mpi_finalize(ierr)


end program twoDChain
