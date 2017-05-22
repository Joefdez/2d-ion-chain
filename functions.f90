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

  subroutine coulombM(nparticles, xx, yy, fx, fy)
    implicit none
    real(kind=8), dimension(:), intent(in)              :: xx, yy
    real(kind=8), dimension(:,:), intent(inout)           :: fx, fy
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
    Axx = ppx
    Ayy = ppy
    Apx = -1.0d0*xx + cfx !-etaC*ppx
    Apy = -1.0*alpha*alpha*yy + cfy !- etaC*ppy
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

end module support_functions_twod
