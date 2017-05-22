# !usr/bin/python

from numpy import *
from matplotlib.pylab import *
ion()
xx = loadtxt("results/posx.dat")
yy = loadtxt("results/posy.dat")

dims = shape(xx)
frac = int(dims[1]*0.1)
xxf = zeros([frac*30])
yyf = zeros([frac*30])

for ii in range(0,30):
    xxf[ii*frac:(ii+1)*frac] = xx[ii,-1*frac:]
    yyf[ii*frac:(ii+1)*frac] = yy[ii,-1*frac:]

del xx
del yy

gauss = 0.001
nx, ny = 200,  100
xl, xr = -6.0, 6.0
yb, yt = -0.2, 0.2
dx = (xr-xl)/(nx-1)
dy = (yt-yb)/(ny-1)


xedges = arange(xl ,xr+dx, dx)
yedges = arange(yb, yt+dy, dy)

H, xedges, yedges = histogram2d(xxf, yyf, bins=(xedges, yedges))
H=H.T

vals = zeros([nx, ny])


chain = figure("chain")
axes  = chain.add_subplot(111)


axes.imshow(H, interpolation='nearest', origin='low',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])

axes.set_aspect(8)
show()
"""
grid=np.zeros([nx, ny])

for ii in range(nx):
    print ii
    xr = xl + dx*ii
    for jj in range(ny):
        f  = 0
        yr = yb + dy*jj
        cx = (xr-xxf)**2/(2.0*gauss)
        cy = (yr-yyf)**2/(2.0*gauss)
        fs  = sum(exp(-cx-cy))
        grid[ii,jj] = fs

axes.imshow(grid)
"""
