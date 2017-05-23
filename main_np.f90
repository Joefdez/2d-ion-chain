program twoDChain


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
  integer                                               :: nsteps, savefreq, nssteps, nparticles, nbath, n_elems, fin
  integer                                               :: traj, local_traj, save_freq, rem
  integer                                               :: ii,jj, kk, ll, mm
  real(kind=8)                                          :: seconds, seconds1
  ! mpi variables



  print*, "Reading solution and system parameters"
  call initialize_system(nparticles, mass, charge, tt, dt, traj, save_freq, long_freq, alpha, ic_radius)

  mass   = mass*uu
  charge = charge*ee
  long_freq = 2.0d0*pi*long_freq
  dst = sqrt(dt)
  nsteps=int(tt/dt)
  char_length = ((charge*charge/(4.0d0*pi*ep0))/(mass*long_freq*long_freq))**(1.0/3.0)
  nssteps = int(nsteps/save_freq)  ! Not saving every single timestep saves memory. Must ask about this
  fin = 0.8*nsteps



  print*, "Reading laser parameters"
  call initialize_laser_chain(del1, del2, delC, Gam, omega0, I1, I2, IC)

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

  local_traj = traj
  n_elems = nssteps*nparticles
  nbath = 3

  include 'allocation.f90'

  stermsCx = sqrt(2.0d0*aDc)
  stermsCy = sqrt(2.0d0*aDc)
  stermsBx(1:nbath) = sqrt(2.0d0*aD1)
  stermsBy(1:nbath) = sqrt(2.0d0*aD1)
  stermsBx((nparticles-nbath+1):nparticles) = sqrt(2.0d0*aD2)
  stermsBy((nparticles-nbath+1):nparticles) = sqrt(2.0d0*aD2)



  open(unit=15, file="eq_pos.dat", action='read')
  do ii = 1, nparticles, 1
    read(15,*) xx0(ii), yy0(ii)
  end do


  do kk=1, local_traj, 1
    call icpgen(nparticles, 0.02d0, xx0, yy0, xxold, yyold)
    call ranseed()
    xxs = 0.0d0
    yys = 0.0d0
    xxs(:,1) = xxold
    yys(:,1) = yyold
    ppxold = 0.0d0
    ppyold = 0.0d0
    ll = 1
    mm = 1
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
      ppxi = ppxold + Apx*dt  + stermsBx*dOmx !+ stermsCx*dOmxc
      ppyi = ppyold + Apy*dt  + stermsBy*dOmy !+ stermsCy*dOmyc
      fx = 0.0d0
      fy = 0.0d0
      call coulombM(nparticles, xxi, yyi, fx2, fy2)
      do jj=1, nparticles, 1
        fx(jj) = sum(fx2(jj,:), 1)
        fy(jj) = sum(fy2(jj,:), 1)
      end do
      call vecA(xxi, yyi, ppxi, ppyi, fx, fy, alpha, aeta1, aeta2, aetaC, nbath, nparticles, Axxi, Ayyi, Apxi, Apyi)
      !call vecB_edges(dst, nparticles, dOmx, dOmy)
      !call vecB_cool(dst, nparticles, dOmxc, dOmyc)
      xxnew   = xxold + 0.5d0*(Axx + Axxi)*dt
      yynew   = yyold + 0.5d0*(Ayy + Ayyi)*dt
      ppxnew  = ppxold + 0.5d0*(Apx + Apxi)*dt + stermsBx*dOmx !+ stermsCx*dOmxc
      ppynew  = ppyold + 0.5d0*(Apy + Apyi)*dt + stermsBy*dOmy !+ stermsCy*dOmyc
      if( mod(ii,save_freq) .eq. 0) then
        ll = ll + 1
        xx2s(:,ll)  = xxnew*xxnew
        yy2s(:,ll)  = yynew*yynew
        ppx2s(:,ll) = ppxnew*ppxnew
        ppy2s(:,ll) = ppynew*ppynew
        xpxs(:,ll)  = xxnew*ppxnew
        ypys(:,ll)  = yynew*ppynew
      end if
      if( ii .ge. fin) then
        xxs(:,mm)   = xxnew
        yys(:,mm)   = yynew
        ppxs(:,mm)  = ppxnew
        ppys(:,mm)  = ppynew
        mm = mm + 1
      end if
      xxold   = xxnew
      yyold   = yynew
      ppxold  = ppxnew
      ppyold  = ppynew
    end do
    xx_av  = (xx_av + xxs)
    yy_av  = (yy_av + yys)
    ppx_av = (ppx_av + ppxs)
    ppy_av = (ppy_av + ppys)
    xx2_av  = (xx2_av + xx2s)
    yy2_av  = (yy2_av + yy2s)
    ppx2_av = (ppx2_av + ppx2s)
    ppy2_av = (ppy2_av + ppy2s)
    xpx_av  = (xpx_av + xpxs)
    ypy_av  = (ypy_av + ypys)
    if(mod(local_traj,20) .eq. 0) then
      open(unit=11, file="posX.dat")
      open(unit=12, file="posY.dat")
      open(unit=13, file="temperatures.dat")
      do jj=1, nparticles
        write(11,*) xx_avt(jj,:)/kk
         write(12,*) yy_avt(jj,:)/kk
         write(13,*) (ppx2_avt(jj,:) + ppy2_avt(jj,:))/kk
      end do
      close(unit=11)
      close(unit=12)
      close(unit=13)
    end if
  end do

  xx_avt   = xx_avt/traj!*char_length/traj
  yy_avt   = yy_avt/traj!*char_length/traj
  xx2_avt  = xx2_avt*char_length*char_length/traj
  yy2_avt  = yy2_avt*char_length*char_length/traj
  ppx_avt  = ppx_avt*char_length*mass*long_freq/traj
  ppy_avt  = ppy_avt*char_length*mass*long_freq/traj
  ppx2_avt = ppx2_avt*char_length*char_length*mass*long_freq*long_freq/(2.0d0*kb)/traj ! Convert to temperature in mK
  ppy2_avt = ppy2_avt*char_length*char_length*mass*long_freq*long_freq/(2.0d0*kb)/traj ! Convert to temperature in mK
  xpx_avt  = xpx_avt*char_length*char_length*mass*long_freq/traj
  ypy_avt  = ypy_avt*char_length*char_length*mass*long_freq/traj


end program twoDChain
