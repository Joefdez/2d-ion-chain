allocate(xx0(1:nparticles))
xx0 = 0.0d0
allocate(yy0(1:nparticles))
yy0 = 0.0d0
allocate(px0(1:nparticles))
px0 = 0.0d0
allocate(py0(1:nparticles))
py0 = 0.0d0
allocate(xxold(1:nparticles))
xxold = 0.0d0
allocate(yyold(1:nparticles))
yyold = 0.0d0
allocate(ppxold(1:nparticles))
ppxold = 0.0d0
allocate(ppyold(1:nparticles))
ppyold = 0.0d0
allocate(xxnew(1:nparticles))
xxnew = 0.0d0
allocate(yynew(1:nparticles))
yynew = 0.0d0
allocate(ppxnew(1:nparticles))
ppxnew = 0.0d0
allocate(ppynew(1:nparticles))
ppynew = 0.0d0
allocate(xxs(1:nparticles, 1:(nsteps-fin)))
xxs = 0.0d0
allocate(yys(1:nparticles, 1:(nsteps-fin)))
yys = 0.0d0
allocate(ppxs(1:nparticles, 1:(nsteps-fin)))
ppxs = 0.0d0
allocate(ppys(1:nparticles, 1:(nsteps-fin)))
ppys = 0.0d0
allocate(xx2s(1:nparticles, 1:nssteps))
xx2s = 0.0d0
allocate(yy2s(1:nparticles, 1:nssteps))
yy2s = 0.0d0
allocate(ppx2s(1:nparticles, 1:nssteps))
ppx2s = 0.0d0
allocate(ppy2s(1:nparticles, 1:nssteps))
ppy2s = 0.0d0
allocate(xpxs(1:nparticles, 1:nssteps))
xpxs = 0.0d0
allocate(ypys(1:nparticles, 1:nssteps))
ypys = 0.0d0
allocate(xxi(1:nparticles))
xxi = 0.0d0
allocate(yyi(1:nparticles))
yyi = 0.0d0
allocate(ppxi(1:nparticles))
ppxi = 0.0d0
allocate(ppyi(1:nparticles))
ppyi = 0.0d0
allocate(fx1(1:nparticles,1:nparticles))
fx1 = 0.0d0
allocate(fx2(1:nparticles,1:nparticles))
fx2 = 0.0d0
allocate(fy1(1:nparticles,1:nparticles))
fy1 = 0.0d0
allocate(fy2(1:nparticles,1:nparticles))
fy2 = 0.0d0
allocate(fx(1:nparticles))
fx = 0.0d0
allocate(fy(1:nparticles))
fy = 0.0d0
allocate(Axx(1:nparticles))
Axx = 0.0d0
allocate(Ayy(1:nparticles))
Ayy = 0.0d0
allocate(Apx(1:nparticles))
Apx = 0.0d0
allocate(Apy(1:nparticles))
Apy = 0.0d0
allocate(Axxi(1:nparticles))
Axxi = 0.0d0
allocate(Ayyi(1:nparticles))
Ayyi = 0.0d0
allocate(Apxi(1:nparticles))
Apxi = 0.0d0
allocate(Apyi(1:nparticles))
Apyi = 0.0d0
allocate(dOmx(1:nparticles))
dOmx = 0.0d0
allocate(dOmy(1:nparticles))
dOmy = 0.0d0
allocate(dOmxc(1:nparticles))
dOmxc = 0.0d0
allocate(dOmyc(1:nparticles))
dOmyc = 0.0d0
allocate(stermsBx(1:nparticles))
stermsBx = 0.0d0
allocate(stermsBy(1:nparticles))
stermsBy = 0.0d0
allocate(stermsCx(1:nparticles))
stermsCx = 0.0d0
allocate(stermsCy(1:nparticles))
stermsCy = 0.0d0

allocate(xx_av(1:nparticles, 1:(nsteps-fin)))
xx_av = 0.0d0
allocate(yy_av(1:nparticles, 1:(nsteps-fin)))
yy_av = 0.0d0
allocate(ppx_av(1:nparticles, 1:(nsteps-fin)))
ppx_av = 0.0d0
allocate(ppy_av(1:nparticles, 1:(nsteps-fin)))
ppy_av = 0.0d0
allocate(xx2_av(1:nparticles, 1:nssteps))
xx2_av = 0.0d0
allocate(yy2_av(1:nparticles, 1:nssteps))
yy2_av = 0.0d0
allocate(ppx2_av(1:nparticles, 1:nssteps))
ppx2_av = 0.0d0
allocate(ppy2_av(1:nparticles, 1:nssteps))
ppy2_av = 0.0d0
allocate(xpx_av(1:nparticles, 1:nssteps))
xpx_av = 0.0d0
allocate(ypy_av(1:nparticles, 1:nssteps))
ypy_av = 0.0d0

allocate(xx_avt(1:nparticles, 1:(nsteps-fin)))
xx_avt = 0.0d0
allocate(yy_avt(1:nparticles, 1:(nsteps-fin)))
yy_avt = 0.0d0
allocate(ppx_avt(1:nparticles, 1:(nsteps-fin)))
ppx_avt = 0.0d0
allocate(ppy_avt(1:nparticles, 1:(nsteps-fin)))
ppy_avt = 0.0d0
allocate(xx2_avt(1:nparticles, 1:nssteps))
xx2_avt = 0.0d0
allocate(yy2_avt(1:nparticles, 1:nssteps))
yy2_avt = 0.0d0
allocate(ppx2_avt(1:nparticles, 1:nssteps))
ppx2_avt = 0.0d0
allocate(ppy2_avt(1:nparticles, 1:nssteps))
ppy2_avt = 0.0d0
allocate(xpx_avt(1:nparticles, 1:nssteps))
xpx_avt = 0.0d0
allocate(ypy_avt(1:nparticles, 1:nssteps))
ypy_avt = 0.0d0
