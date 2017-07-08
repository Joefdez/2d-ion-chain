!!!! Time step using Platen weak 2.0 scheme !!!!

  call coulombM(nparticles, xxold, yyold, fx1, fy1, invD1)
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
  call coulombM(nparticles, xxi, yyi, fx2, fy2, invD2)
  do jj=1, nparticles, 1
    fx(jj) = sum(fx2(jj,:), 1)
    fy(jj) = sum(fy2(jj,:), 1)
  end do
  call vecA(xxi, yyi, ppxi, ppyi, fx, fy, alpha, aeta1, aeta2, aetaC, nbath, nparticles, Axxi, Ayyi, Apxi, Apyi)
  xxnew   = xxold + 0.5d0*(Axx + Axxi)*dt
  yynew   = yyold + 0.5d0*(Ayy + Ayyi)*dt
  ppxnew  = ppxold + 0.5d0*(Apx + Apxi)*dt + stermsBx*dOmx !+ stermsCx*dOmxc
  ppynew  = ppyold + 0.5d0*(Apy + Apyi)*dt + stermsBy*dOmy !+ stermsCy*dOmyc
