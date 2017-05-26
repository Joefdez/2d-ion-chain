module support_functions_twod
  use randomModule
  implicit none
  contains

  subroutine coulomb(nparticles, xx, yy, fx, fy)

    implicit none
    real(kind=8), dimension(:), intent(in)              :: xx, yy
    real(kind=8), dimension(:), intent(inout)           :: fx, fy
    integer, intent(in)                                 :: nparticles
    integer                                             :: ii, jj
    real(kind=8)                                        :: dist
    fx = 0.0d0
    fy = 0.0d0
    do ii=1, nparticles, 1
      do jj=ii+1, nparticles, 1
        dist = 0.0d0
        dist = sqrt( (xx(ii)-xx(jj))**2 + (yy(ii)-yy(jj))**2 )
        dist = dist*dist*dist
        fx(ii) = fx(ii) +  (xx(ii)-xx(jj)) / dist
        fx(jj) = fx(jj) -  (xx(ii)-xx(jj)) / dist
        fy(ii) = fy(ii) +  (yy(ii)-yy(jj)) / dist
        fy(ii) = fy(jj) -  (yy(ii)-yy(jj)) / dist
      end do
    end do

  end subroutine coulomb

  subroutine coulombM(nparticles, xx, yy, fx, fy, invD)
    implicit none
    real(kind=8), dimension(:), intent(in)              :: xx, yy
    real(kind=8), dimension(:,:), intent(inout)         :: fx, fy, invD
    integer, intent(in)                                 :: nparticles
    integer                                             :: ii, jj
    real(kind=8)                                        :: dist
    invD = 0.0d0
    fx = 0.0d0
    fy = 0.0d0
    do ii=1, nparticles, 1
      do jj=ii+1, nparticles, 1
        dist = 0.0d0
        dist = sqrt( (xx(ii)-xx(jj))**2 + (yy(ii)-yy(jj))**2 )
        invD(ii,jj) = 1.0d0/dist
        invD(jj,ii) = 1.0d0/dist
        dist = dist*dist*dist
        fx(ii,jj) = (xx(ii)-xx(jj)) / dist
        fx(jj,ii) = -1.0d0*(xx(ii)-xx(jj)) / dist
        fy(ii,jj) = (yy(ii)-yy(jj)) / dist
        fy(jj,ii) = -1.0d0*(yy(ii)-yy(jj)) / dist
      end do
    end do

  end subroutine coulombM

  subroutine vecA(xx, yy, ppx, ppy, cfx, cfy, alpha, eta1, eta2, etaC, nbath, nparticles, Axx, Ayy, Apx, Apy)

    implicit none
    real(kind=8), dimension(:), intent(in)              :: xx, yy, ppx, ppy
    real(kind=8), dimension(:), intent(in)              :: cfx, cfy
    real(kind=8), intent(in)                            :: alpha, eta1, eta2, etaC
    integer, intent(in)                                 :: nparticles, nbath
    real(kind=8), dimension(:), intent(inout)           :: Axx, Ayy, Apx, Apy

    ! Entire chain: harmnonic and coulomb forces, and cooling laser
    Axx(1:nparticles) = ppx(1:nparticles)
    Ayy(1:nparticles) = ppy(1:nparticles)
    Apx(1:nparticles) = -1.0d0*xx(1:nparticles) + cfx(1:nparticles) !-etaC*ppx
    Apy(1:nparticles) = -1.0d0*alpha*alpha*yy(1:nparticles) + cfy(1:nparticles) !- etaC*ppy
    ! Thermal baths at edges
    Apx(1:nbath) = Apx(1:nbath) - eta1*ppx(1:nbath)
    Apx((nparticles-nbath+1):nparticles) = Apx((nparticles-nbath+1):nparticles) - eta2*ppx((nparticles-nbath+1):nparticles)
    Apy(1:nbath) = Apy(1:nbath) - eta1*ppy(1:nbath)
    Apy((nparticles-nbath+1):nparticles) = Apy((nparticles-nbath+1):nparticles) - eta2*ppy((nparticles-nbath+1):nparticles)

  end subroutine

  subroutine vecB_edges(dst, nparticles, dOmx, dOmy)

    implicit none
    real(kind=8), intent(in)                    :: dst
    integer, intent(in)                         :: nparticles
    real(kind=8), dimension(:), intent(inout)   :: dOmx, dOmy
    real(kind=8)                                :: g1, g2
    integer                                     :: ii

    call ranseed()
    do ii=1,3,1
      call bm(g1, g2)
      dOmx(ii) = dst*g1
      dOmx(nparticles+1-ii) = dst*g2
      call bm(g1, g2)
      dOmy(ii) = dst*g1
      dOmy(nparticles+1-ii) = dst*g2
    end do


  end subroutine vecB_edges

  subroutine vecB_cool(dst, nparticles, dOmx, domy)

    implicit none
    real(kind=8), intent(in)                    :: dst
    integer, intent(in)                         :: nparticles
    real(kind=8), dimension(:), intent(inout)   :: dOmx, dOmy
    real(kind=8)                                :: g1, g2
    integer                                     :: ii


    call ranseed()
    do ii=1, nparticles, 1
        call bm(g1, g2)
        dOmx(ii) = dst*g1
        dOmy(ii) = dst*g2
    end do

  end subroutine vecB_cool

  subroutine local_energy(nparticles, alpha, xx, yy, invD, ppx, ppy, energy)
    implicit none
    integer, intent(in)                                   :: nparticles
    real(kind=8), intent(in)                              :: alpha
    real(kind=8), dimension(:), intent(in)                :: xx, yy, ppx, ppy
    real(kind=8), dimension(:,:), intent(in)              :: invD
    real(kind=8), dimension(:), intent(inout)             :: energy
    real(kind=8)                                          :: kin_term
    integer                                               :: ii
    energy = 0.0d0

    ! Add kinetic energy contribution and harmonic potential energy contributions
    energy(1:nparticles) = 0.5*ppx(1:nparticles)*ppx(1:nparticles) + 0.5d0*ppx(1:nparticles)*ppx(1:nparticles) &
              + 0.5d0*xx(1:nparticles)*xx(1:nparticles) + 0.5d0*alpha*alpha*yy(1:nparticles)*yy(1:nparticles)
    do ii=1, nparticles, 1
      energy(ii) = energy(ii) + 0.5d0*sum(invD(ii,:), 1)          ! Add the coulomb potential contribution
    end do

  end subroutine local_energy

  subroutine heat_current(nparticles, cfx, cfy, ppx, ppy, hc)
    implicit none
    integer, intent(in)                                   :: nparticles
    real(kind=8), dimension(:,:), intent(in)              :: cfx, cfy
    real(kind=8), dimension(:), intent(in)                :: ppx, ppy
    real(kind=8), dimension(:,:), intent(inout)              :: hc
    integer                                               :: ii, jj
    hc=0.0d0
    do ii=1, nparticles, 1
      do jj=ii, nparticles, 1
        hc(ii,jj) = 0.5d0*(ppx(ii)+ppx(jj))*cfx(ii,jj) + 0.5d0*(ppy(ii)+ppy(jj))*cfy(ii,jj)
        hc(jj,ii) = -1.0d0*hc(ii,jj)
      end do
    end do

  end subroutine heat_current

  subroutine current_Flux(hc, energy, xx, yy, ppx, ppy, nparticles, JJintx, JJinty)
    implicit none
    real(kind=8), dimension(:,:), intent(in)              :: hc !cfx, cfy
    real(kind=8), dimension(:), intent(in)                :: xx, yy, ppx, ppy, energy
    integer, intent(in)                                   :: nparticles
    real(kind=8), intent(inout)                           :: JJintx, JJinty
    real(kind=8)                                          :: JJ0
    integer                                               :: nn, ll

    JJintx = 0.0d0
    JJinty = 0.0d0



    do nn=1, nparticles, 1
      JJintx = JJintx + energy(nn)*ppx(nn)
      JJinty = JJintx + energy(nn)*ppy(nn)
      do ll=1, nparticles, 1
        if(nn .eq. ll) cycle
        !JJ0 = 0.5d0*(ppx(nn+1)+ppx(ll))*cfx(nn+1,ll) + 0.5d0*(ppy(nn+1)+ppy(ll))*cfy(nn+1,ll)
        JJintx = JJintx + 0.5d0*(xx(nn)-xx(ll))*hc(nn,ll)
        JJinty = JJinty + 0.5d0*(yy(nn)-yy(ll))*hc(nn,ll)
      end do
    end do

  end subroutine current_Flux

end module support_functions_twod
