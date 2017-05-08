module initcond_generator

  use constants
  use randomModule

  implicit none

 contains

   subroutine icpgen(n_particles, radius, xx0, yy0, xx, yy)
     ! Introduce uniform noise in equilibrum positions. Requieres file containing equilibrium positions.
     ! Data in file must be columnwise
     implicit none
     real(kind=8), dimension(:), intent(in)     :: xx0, yy0                ! vectors to contain initial conditions
     real(kind=8), dimension(:), intent(inout)  :: xx, yy
     integer, intent(in)                        :: n_particles
     real(kind=8), intent(in)                   :: radius
     integer                                    :: ii
     logical                                    :: calc=.true.
     real(kind=8), dimension(:), allocatable    :: rg
     !open(unit=15, file="eq_pos.dat", action='read')
     call ranseed()
     allocate(rg(1:2))
     do ii = 1, n_particles, 1
     !   read(15,*) qq0((ii*dim+1):(ii*dim+dim))
       do
        call RANDOM_NUMBER(rg)                                        ! generate three uniformly distributed random numbers
        rg=rg*rg*radius*radius                                        ! scale to chosen radius
        if( sqrt(sum(rg)) .le. radius ) exit                          ! if the distance is larger than the radius, repeat
       end do
       !qq0((ii*dim+1):(ii*dim+dim)) = qq0((ii*dim+1):(ii*dim+dim)) + rg           ! add in the noise
       xx(ii) = xx0(ii) + rg(1)
       yy(ii) = yy0(ii) + rg(2)
     end do

   end subroutine icpgen

end module initcond_generator
