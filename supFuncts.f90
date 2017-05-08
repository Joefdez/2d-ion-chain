module support_functions
  use constants
  use randomModule
  implicit none

  contains

  subroutine Cforce(posx, posy, n_particles, invD, Cfx, Cfy)
    implicit none

    real(kind=8), dimension(:), intent(in)      :: posx, posy
    integer, intent(in)                         :: n_particles
    real(kind=8), dimension(:,:), intent(inout) :: invD, Cfx, Cfy
    integer                                     :: ii, jj
    real(kind=8)                                :: dist

    Cfx = 0.0d0
    Cfy = 0.0d0
    invD = 0.0d0
    do ii=1, n_particles, 1
      dist = 0.0d0
      do jj=ii+1, n_particles, 1
        dist = (posx(ii)-posx(jj))*(posx(ii)-posx(jj)) + (posy(ii)-posy(jj))*(posy(ii)-posy(jj))
        dist = sqrt(dist)
        invD(ii,jj) = 1.0d0/dist
        invD(jj,ii) = 1.0d0/dist
        dist = (dist*dist*dist)
        Cfx(ii,jj) = (posx(ii)-posx(jj))/dist
        Cfx(jj,ii) = -1.0d0 * (posx(ii)-posx(jj))/dist
        Cfy(ii, jj) = (posy(ii)-posy(jj))/dist
        Cfy(jj, ii) = -1.0d0 * (posy(ii)-posy(jj))/dist
      end do
    end do

  end subroutine Cforce

  subroutine ddt(YY, AA, n_particles, eta1, eta2, alpha, Cforcex, Cforcey)
    implicit none
    real(kind=8), dimension(:), intent(in)    :: YY
    real(kind=8), dimension(:,:), intent(in)  :: Cforcex, Cforcey
    real(kind=8), intent(in)                  :: eta1, eta2, alpha
    integer, intent(in)                       :: n_particles
    real(kind=8), dimension(:), intent(inout) :: AA
    integer                                   :: ii

    AA = 0.0d0
    AA(1:2*n_particles)  = YY(3*n_particles+1:)
    do ii=1, n_particles, 1
      AA(ii+2*n_particles) = -0.5d0*YY(ii) + sum(Cforcex(ii,1:n_particles),1)
    end do
    do ii=1, n_particles, 1
      AA(ii+3*n_particles) = -0.5d0*alpha*alpha*YY(n_particles+ii) + sum(Cforcey(ii,1:n_particles),1)
    end do

    do ii=1, 3, 1
      AA(2*n_particles+ii) = AA(2*n_particles+ii) - eta1*YY(2*n_particles+ii)
      AA(3*n_particles-(ii-1)) = AA(3*n_particles-(ii-1)) - eta2*YY(3*n_particles-(ii-1))
      AA(3*n_particles+ii) = AA(3*n_particles+ii) - eta1*YY(3*n_particles+ii)
      AA(4*n_particles-(ii-1)) = AA(4*n_particles-(ii-1)) - eta2*YY(4*n_particles-(ii-1))
    end do
  end subroutine ddt

  subroutine stoch_vector(dst, n_particles, stoch_terms, dStoc)
    implicit none
    real(kind=8), intent(in)                  :: dst
    integer, intent(in)                       :: n_particles
    real(kind=8), dimension(:), intent(in)     :: stoch_terms
    real(kind=8), dimension(:), intent(inout) :: dStoc
    real(kind=8), dimension(:), allocatable   :: dOm
    integer                                   :: nn
    real(kind=8)                              :: rg1, rg2

    allocate(dOm(1:4*n_particles))
    dOm = 0.0d0
    do nn=1, 3, 1
      call bm(rg1, rg2)
      dOm(2*n_particles+nn) = rg1*dst
      dOm(3*n_particles-(nn+1)) = rg2*dst
      call bm(rg1, rg2)
      dOm(3*n_particles+1) = rg1*dst
      dOm(4*n_particles-(nn-1)) = rg2*dst
    end do

    dStoc = stoch_terms*dOm

  end subroutine stoch_vector


end module support_functions
