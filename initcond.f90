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
     real(kind=8)                               :: ang
     !open(unit=15, file="eq_pos.dat", action='read')
     call ranseed()
     allocate(rg(1:2))
     do ii = 1, n_particles, 1
       call RANDOM_NUMBER(rg)
       ang = 2*pi*rg(1)
       xx(ii) = xx0(ii) + rg(2)*radius*dcos(ang)
       yy(ii) = yy0(ii) + rg(2)*radius*dsin(ang)
     end do

   end subroutine icpgen

   subroutine icmomgen(n_particles, speed, px, py)
     ! Introduce uniform noise in equilibrum positions. Requieres file containing equilibrium positions.
     ! Data in file must be columnwise
     implicit none
     real(kind=8), dimension(:), intent(inout)  :: px, py
     integer, intent(in)                        :: n_particles
     real(kind=8), intent(in)                   :: speed
     integer                                    :: ii
     logical                                    :: calc=.true.
     real(kind=8), dimension(:), allocatable    :: rg
     real(kind=8)                               :: ang
     !open(unit=15, file="eq_pos.dat", action='read')

     call ranseed()
     allocate(rg(1:2))
     do ii=1, n_particles, 1
       call RANDOM_NUMBER(rg)
       ang = 2*pi*rg(1)
       px(ii) = rg(2)*speed*dcos(ang)
       py(ii) = rg(2)*speed*dsin(ang)
     end do

   end subroutine icmomgen

end module initcond_generator
