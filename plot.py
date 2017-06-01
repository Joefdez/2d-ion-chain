# !usr/bin/python

from numpy import *
from matplotlib.pylab import *
ion()
#xx = loadtxt("results/posx.dat")
#yy = loadtxt("results/posy.dat")
print("Loading ion coordinate files.")
xx = loadtxt("posX.dat")
yy = loadtxt("posY.dat")

dims = shape(xx)
frac = int(dims[1])
xxs = zeros(30)
xxf = zeros([frac*30])
yyf = zeros([frac*30])


print("Flattenning position arrays to compose bivariate histogram.")
for ii in range(0,30):
    xxf[ii*frac:(ii+1)*frac] = xx[ii,-1*frac:]
    yyf[ii*frac:(ii+1)*frac] = yy[ii,-1*frac:]
    xxs[ii] = mean(xx[ii,:])

del xx
del yy

print("Composing the bivariate histogram")
gauss = 0.001
nx, ny = 400,  200
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

axes.imshow(H, origin='low',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])

axes.set_aspect(8)
show()


print("Loading ion temperature profiles.")
temps=loadtxt("temperatures.dat")
print("Calculating temperature profile")
tf = zeros(30)
ns = shape(temps)
ns = 0.2*ns[1]
for ii in range(30):       # Temperature in K
    tf[ii] = mean(temps[ii,-int(ns):])

print("Plotting temperature profile.")
tb, th = 3.0, 9.0
txl, txr = -6.0, 6.0
temp = figure("Temperature")
axt = temp.add_subplot(111)
axt.set_xlim([txl, txr])
axt.set_ylim([tb, th])
axt.plot(xxs, 1000*tf) # Plot temperature in mK
