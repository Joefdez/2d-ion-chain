program twoDChain

  use mpi
  use constants
  use initcond_generator
  use initialization
  use randomModule
  use support_functions_twod

  implicit none

  include "declarations.f90"

  ! mpi variables
  integer :: rank, procs, status(MPI_STATUS_SIZE), alloc_err, source, ierr
  call MPI_INIT(ierr)                                                                               ! Neccesary mpi initialization calls
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, procs, ierr)

  if(rank .eq. 0 ) then
    print*, "Reading solution and system parameters"
    call initialize_system(nparticles, mass, charge, tt, dt, traj, save_freq, long_freq, alpha, ic_radius, initT)
  end if


  ! Share system data between differente processes
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
  call mpi_bcast(initT, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)

  call mpi_barrier(mpi_comm_world, ierr)              ! Wait for everyone else to get here

  mass   = mass*uu
  charge = charge*ee
  long_freq = 2.0d0*pi*long_freq
  dst = sqrt(dt)
  nsteps=int(tt/dt)
  char_length = ((charge*charge/(4.0d0*pi*ep0))/(mass*long_freq*long_freq))**(1.0/3.0)
  nssteps = int(nsteps/save_freq)  ! Not saving every single timestep saves memory. Must ask about this
  fin = 0.8d0*nsteps

  if(rank .eq. 0) then
    print*, "Reading laser parameters"
    call initialize_laser_chain(del1, del2, delC, Gam, omega0, I1, I2, IC)
  end if

  ! Share the laser parameters between the difference processes
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
  initT = initT * (kb/(char_length*char_length*mass*long_freq*long_freq))
  initSpeed = sqrt(2*initSpeed)

  call doppler_values(k1, Gam, del1, I1, eta1, D1)
  call doppler_values(k2, Gam, del2, I2, eta2, D1)
  !call doppler_vales(kC, Gam, delC, IC etaC, DC)


  ! Calculate dimensionless Doppler cooling parameters
  call dimensionless_doppler_values(eta1, D1, mass, long_freq, char_length, aeta1, aD1)
  call dimensionless_doppler_values(eta2, D2, mass, long_freq, char_length, aeta2, aD2)
  !call dimensionless_doppler_values(etaC, DC, mass, long_freq, char_length, aetaC, aDc)

  local_traj = traj/procs
  rem = mod(traj, procs)
  if (rank .lt. rem) local_traj = local_traj + 1
  n_elems = nssteps*nparticles
  nbath = 3

  include 'allocation.f90'

  stermsBx(1:nbath) = sqrt(2.0d0*aD1)
  stermsBy(1:nbath) = sqrt(2.0d0*aD1)
  stermsBx((nparticles-nbath+1):nparticles) = sqrt(2.0d0*aD2)
  stermsBy((nparticles-nbath+1):nparticles) = sqrt(2.0d0*aD2)
  !stermsCx = sqrt(2.0d0*aDc)
  !stermsCy = sqrt(2.0d0*aDc)

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
  call sleep(rank) ! Delay for getting different seeds

  do kk=1, local_traj, 1
    print*, "Proc.", rank, "on trajectory", kk
    xxs = 0.0d0
    yys = 0.0d0
    ppxold = 0.0d0
    ppyold = 0.0d0
    ll = 0
    mm = 1
    JJix = 0.0d0
    JJiy = 0.0d0
    call icpgen(nparticles, 0.02d0, xx0, yy0, xxold, yyold)
    call icmomgen(nparticles, initSpeed, ppxold, ppyold)
    call ranseed()
!   JJix_av  = 0.0d0
!   JJiy_av = 0.0d0
    do ii=1, nsteps, 1

      include 'solvers/platen_step.f90'     !! This file contains the lines which calculate the time step.


      if( mod(ii,save_freq) .eq. 0) then      ! Save chosen time steps to get time evolution of variables
        ll = ll + 1
        !xx2s(:,ll)  = xxnew*xxnew
        !yy2s(:,ll)  = yynew*yynew
        ppx2s(:,ll) = ppxnew*ppxnew
        ppy2s(:,ll) = ppynew*ppynew
        xpxs(:,ll)  = xxnew*ppxnew
        ypys(:,ll)  = yynew*ppynew
      end if

      if( ii .gt. fin) then                     ! Save last time steps for steady-state calculation
          xxs(:,mm)   = xxnew
          yys(:,mm)   = yynew
          call local_energy(nparticles, alpha, xxold, yyold, invD1, ppxold, ppyold, energy)
          call heat_current(nparticles, fx1, fy1, ppxold, ppyold, hc)
          call current_Flux(hc, energy, xxold, yyold, ppxold, ppyold, nparticles, JJix, JJiy)
          JJix_av = JJix_av + JJix/(nsteps-fin-1)
          JJiy_av = JJiy_av + JJiy/(nsteps-fin-1)
          JJix_av_v(kk) = JJix_av_v(kk) + JJix/(nsteps-fin-1)
          JJiy_av_v(kk) = JJiy_av_v(kk) + JJiy/(nsteps-fin-1)
          mm = mm + 1
      end if
      xxold   = xxnew
      yyold   = yynew
      ppxold  = ppxnew
      ppyold  = ppynew
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!! Calculate correlation matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !xx2s(:,nssteps)  = xxnew*xxnew
    !yy2s(:,nssteps)  = yynew*yynew
    ppx2s(:,nssteps) = ppxnew*ppxnew
    ppy2s(:,nssteps) = ppynew*ppynew
    xpxs(:,nssteps)  = xxnew*ppxnew
    ypys(:,nssteps)  = yynew*ppynew
    call local_energy(nparticles, alpha, xxold, yyold, invD1, ppxold, ppyold, energy)
    call heat_current(nparticles, fx1, fy1, ppxold, ppyold, hc)
    call current_Flux(hc, energy, xxold, yyold, ppxold, ppyold, nparticles, JJix, JJiy)
    JJix_s(nssteps) = JJix
    JJiy_s(nssteps) = JJiy

    ! Write results to file
    if(rank .eq. 0 .and. kk .eq. 1) then
      open(unit=11, file="posX.dat")
      open(unit=12, file="posY.dat")
      print*, "Printing steady-state spatial distribution of "
      do jj=1, nparticles, 1
        write(11,*) xxs(jj,:)
        write(12,*) yys(jj,:)
      end do
      close(unit=11)
      close(unit=12)
    end if
    !xx2_av  = (xx2_av + xx2s)
    !yy2_av  = (yy2_av + yy2s)
    ppx2_av = (ppx2_av + ppx2s)
    ppy2_av = (ppy2_av + ppy2s)
    xpx_av  = (xpx_av + xpxs)
    ypy_av  = (ypy_av + ypys)
    xx_av   = xx_av + xxs
    JJix_sav = JJix_sav + JJix_s
    JJiy_sav = JJiy_sav + JJiy_s
    if( ( mod(kk,5) .eq. 0) .and. (kk .lt. local_traj) ) then
     print*, "Writing PARTIAL results to files after ", kk, "trajectories."

     include 'mpi_comms.f90'
     call mpi_reduce(kk, traj, 1, mpi_integer, mpi_sum, 0, mpi_comm_world, ierr)

     print*, "Finished writing up to", kk
     if(rank .eq. 0) then
      !xx2_avt  = xx2_avt*char_length*char_length/traj
      !yy2_avt  = yy2_avt*char_length*char_length/traj
      ppx2_avt = ppx2_avt*char_length*char_length*mass*long_freq*long_freq/(2.0d0*kb)/traj ! Convert to temperature in K
      ppy2_avt = ppy2_avt*char_length*char_length*mass*long_freq*long_freq/(2.0d0*kb)/traj ! Convert to temperature in K
      xpx_avt  = xpx_avt*char_length*char_length*mass*long_freq/traj
      ypy_avt  = ypy_avt*char_length*char_length*mass*long_freq/traj

      open(unit=11, file="heatflux.dat")
      open(unit=12, file="temperatures.dat")
      write(11,*) JJix_avt/traj
      write(11,*) JJiy_avt/traj
      do jj=1, nparticles
       write(12,*) ppx2_avt(jj,:) + ppy2_avt(jj,:)
      end do
      close(unit=11)
      close(unit=12)
      open(unit=13, file="posXav.dat")
      do jj=1, nparticles, 1
        write(13,*) xx_av/traj
      end do
      close(unit=13)
      open(unit=14, file="JJxt.dat")
      write(14,*) JJix_savt/traj
      close(unit=14)
      open(unit=14, file="JJyt.dat")
      write(14,*) JJiy_savt/traj
      close(unit=14)
     end if
     !xx2_avt  = 0.0d0
     !yy2_avt  = 0.0d0
     ppx2_avt = 0.0d0
     ppy2_avt = 0.0d0
     xpx_avt  = 0.0d0
     ypy_avt  = 0.000
    end if
  end do
  print*,"Proc ", rank, " finished integrating"


  call mpi_barrier(mpi_comm_world, ierr)
  include 'mpi_comms.f90'
  call mpi_reduce(local_traj, traj, 1, mpi_integer, mpi_sum, 0, mpi_comm_world, ierr)


  call mpi_barrier(mpi_comm_world, ierr)

  if(rank .eq. 0) then
    seconds = mpi_wtime() - seconds
    print*, "total integration + partial writing time:", seconds, seconds/3600.0d0
    print*,traj

    !!!!!!!!!!!!!!!!!!!!!!!!!!!! Calculate correlation matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !xx2_avt  = xx2_avt*char_length*char_length/traj
    !yy2_avt  = yy2_avt*char_length*char_length/traj
    ppx2_avt = ppx2_avt*char_length*char_length*mass*long_freq*long_freq/(2.0d0*kb)/traj ! Convert to temperature in K
    ppy2_avt = ppy2_avt*char_length*char_length*mass*long_freq*long_freq/(2.0d0*kb)/traj ! Convert to temperature in K
    xpx_avt  = xpx_avt*char_length*char_length*mass*long_freq/traj
    ypy_avt  = ypy_avt*char_length*char_length*mass*long_freq/traj

    ! Write results to file
    open(unit=11, file="heatflux.dat")
    open(unit=12, file="temperatures.dat")
    write(11,*) JJix_avt/traj
    write(11,*) JJiy_avt/traj
    do jj=1, nparticles
     write(12,*) ppx2_avt(jj,:) + ppy2_avt(jj,:)
    end do
    close(unit=11)
    close(unit=12)
    open(unit=14, file="JJxt.dat")
    write(14,*) JJix_savt/traj
    close(unit=14)
    open(unit=14, file="JJyt.dat")
    write(14,*) JJiy_savt/traj
    close(unit=14)
    print*, "Writing FINAL results to files after ", traj , "trajectories."
    seconds1 = mpi_wtime() - seconds1
    print*, "integration + message + write time:", seconds1
  end if

  call mpi_finalize(ierr)


end program twoDChain
