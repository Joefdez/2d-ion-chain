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
real(kind=8), dimension(:,:), allocatable             :: xx_avo, yy_avo, ppx_avo, ppy_avo
real(kind=8), dimension(:,:), allocatable             :: xx2_avo, yy2_avo, ppx2_avo, ppy2_avo
real(kind=8), dimension(:,:), allocatable             :: xPx_av, yPy_av
real(kind=8), dimension(:,:), allocatable             :: xPx_avo, yPy_avo
real(kind=8), dimension(:), allocatable               :: xxi, yyi, ppxi, ppyi
real(kind=8), dimension(:,:), allocatable             :: fx1, fy1, fx2, fy2
real(kind=8), dimension(:), allocatable               :: fx, fy
real(kind=8), dimension(:), allocatable               :: Axx, Axxi, Ayy, Ayyi
real(kind=8), dimension(:), allocatable               :: Apx, Apxi, Apy, Apyi
real(kind=8), dimension(:), allocatable               :: dOmx, dOmy, dOmxc, dOmyc
real(kind=8), dimension(:), allocatable               :: stermsBx, stermsCx, stermsBy, stermsCy

real(kind=8), dimension(:), allocatable                :: JJix_s, JJix_sav, JJix_savt
real(kind=8), dimension(:), allocatable                :: JJiy_s, JJiy_sav, JJiy_savt

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
integer                                               :: ii,jj, kk, ll, mm, nn
real(kind=8)                                          :: seconds, seconds1
real(kind=8), dimension(:), allocatable               :: energy
real(kind=8)                                          :: JJix, JJiy
!real(kind=8), dimension(:), allocatable               :: JJix_s, JJiy_s
!real(kind=8), dimension(:), allocatable               :: hcx, hcy, hcx_av, hcy_av, hcx_avt, hcy_avt
real(kind=8), dimension(:,:), allocatable             :: hcx, hcy
real(kind=8), dimension(:,:), allocatable             :: invD1, invD2
real(kind=8)                                          :: JJix_av, JJiy_av, JJix_avt,JJiy_avt
real(kind=8), dimension(:), allocatable               :: JJix_av_v, JJiy_av_v
real(kind=8)                                          :: errJJix, errJJiy, errJJix_t, errJJiy_t

real(kind=8)                                          :: initT, initSpeed
