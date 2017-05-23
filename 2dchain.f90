program twodchain

  use constants
  use initcond_generator
  implicit none

  real(kind=8), dimension(:), allocatable         :: xx0, yy0
  real(kind=8), dimension(:,:), allocatable       :: xx, yy, pxx, pyy
  real(kind=8), dimension(:), allocatable         :: vYY0, vYYn, vYYi
  real(kind=8), dimension(:), allocatable         :: force, AA, AAi !,Cf1, Cf2

  real(kind=8)                                    :: tt, dt
  integer                                         :: nparticles, nsteps, nsteps_s
  real(kind=8)                                    :: alpha, long_freq
  real(kind=8)                                    :: mass, charge

  integer                                          :: ii, jj, kk


  nparticles=30
  tt = 100.0d0
  dt = 1.0d-3
  nsteps=int(tt/dt)
  nsteps_s = nsteps
  alpha=9.0d0


  allocate(xx0(1:nparticles))
  xx0 = 0.0d0
  allocate(yy0(1:nparticles))
  yy0 = 0.0d0
  allocate(xx(1:nparticles,1:nsteps_s))
  xx  = 0.0d0
  allocate(yy(1:nparticles,1:nsteps_s))
  yy  = 0.0d0
  allocate(pxx(1:nparticles,1:nsteps_s))
  pxx = 0.0d0
  allocate(pyy(1:nparticles,1:nsteps_s))
  pyy = 0.0d0
  allocate(vYY0(1:4*nparticles))
  vYY0 = 0.0d0
  allocate(vYYn(1:4*nparticles))
  vYYn = 0.0d0
  allocate(vYYi(1:4*nparticles))
  vYYi = 0.0d0
  allocate(force(1:2*nparticles))
  force = 0.0d0
  allocate(AA(1:4*nparticles))
  AA = 0.0d0
  allocate(AAi(1:4*nparticles))
  AAi = 0.0d0


  open(unit=15, file="eq_pos.dat", action='read')
  do ii = 1, nparticles, 1
    read(15,*) xx0(ii), yy0(ii)
  end do


  call icpgen(nparticles, 0.02d0, xx0, yy0, vYY0(1:nparticles), vYY0((nparticles+1):))
  do ii=1, nsteps_s, 1
    print*,ii
    xx(:,ii) = vYY0(1:nparticles)
    yy(:,ii) = vYY0((nparticles+1):2*nparticles)
    pxx(:,ii) = vYY0((2*nparticles+1):3*nparticles)
    pyy(:,ii) = vYY0((3*nparticles+1):4*nparticles)
    call Cforce(nparticles, vYY0(1:nparticles), vYY0((nparticles+1):2*nparticles), force)
    call vecA(nparticles, alpha, AA, vYY0(1:2*nparticles), vYY0((2*nparticles+1):), force)
    vYYi = vYY0 + AA*dt
    call Cforce(nparticles, vYYi(1:nparticles), vYYi((nparticles+1):), force)
    call vecA(nparticles, alpha, AAi, vYYi(1:2*nparticles), vYYi((2*nparticles+1):), force)
    vYYn = vYY0 + 0.5*(AA+AAi)*dt
    !xx(:,ii) = vYYn(1:nparticles)
    !yy(:,ii) = vYYn((nparticles+1):2*nparticles)
    !pxx(:,ii) = vYYn((2*nparticles+1):3*nparticles)
    !pyy(:,ii) = vYYn((3*nparticles+1):4*nparticles)
    vYY0 = vYYn
  end do

  xx(:,nsteps_s) = vYY0(1:nparticles)
  yy(:,nsteps_s) = vYY0((nparticles+1):2*nparticles)
  pxx(:,nsteps_s) = vYY0((2*nparticles+1):3*nparticles)
  pyy(:,nsteps_s) = vYY0((3*nparticles+1):4*nparticles)

  print*, "writing"
  open(unit=11, file="results/posx.dat")
  open(unit=12, file="results/posy.dat")
  do ii=1, nparticles, 1
    write(11,*) xx(jj,:)
    write(12,*) yy(jj,:)
  end do
  close(unit=11)
  close(unit=12)
  print*, "written"

  contains

    subroutine CForce(npart, xx, yy, ff)
      implicit none
      integer, intent(in)                       :: npart
      real(kind=8), dimension(:), intent(in)    :: xx, yy ! positions
      real(kind=8), dimension(:), intent(inout) :: ff     ! force array
      integer                                   :: ii, jj
      real(kind=8)                              :: dist

      ff = 0.0d0
      do ii=1, npart, 1
          do jj=1, npart, 1
            if(jj .eq. ii) cycle
              dist = sqrt((xx(ii)-xx(jj))*(xx(ii)-xx(jj)) + (yy(ii)-yy(jj))*(yy(ii)-yy(jj)))
              dist = dist*dist*dist
              ff(ii) = ff(ii) + (xx(ii)-xx(jj))/dist
              ff(ii+npart) =  ff(ii+npart) + (yy(ii)-yy(jj))/dist
          end do
      end do
    end subroutine CForce

    subroutine vecA(npart, alpha, AA, pos, mom, force)
      implicit none
      integer, intent(in)                         :: npart
      real(kind=8), intent(in)                    :: alpha
      real(kind=8), dimension(:), intent(in)      :: pos, mom, force
      real(kind=8), dimension(:), intent(inout)                 :: AA

      AA(1:2*npart) = mom
      AA((2*npart+1):(3*npart)) = -1.0d0*pos(1:npart) + force(1:npart)
      AA((3*npart+1):) = -1.0*alpha*alpha*pos((npart+1):) + force((npart+1):)

    end subroutine vecA

end program twodchain
