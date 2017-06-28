from pylab import *
import random

matplotlib.rcParams.update({'font.size': 1, 'xtick.major.size': 0, 'ytick.major.size': 0, 'xtick.major.width': 0, 'ytick.major.width': 0, 'ytick.minor.size': 0, 'xtick.minor.size': 0, 'axes.linewidth': 0, 'text.usetex': True, 'font.family': 'serif', 'font.serif': 'Times New Roman'})


# Set manually -- this is the length of the box is image is meant for in Mpc/h
# (e.g. set 250 for Bolshoi, 500 for Millennium, 1000 for MDPL...)
compboxlen = 500

# Sets the physical size of the output image.
pixmin = 250 # Default at 250. Consider lowering if compboxlen>1000
pix = pixmin * compboxlen/500.

# Info on your display
xpix = 1920 # horizontal screen pixels
ypix = 1200 # vertical screen pixels
ss = 27 # diagonal screen size in inches




# Read particle data from a pre-run simulation
with open('simulationData.binary', 'rb') as f:
    mass = np.fromfile(f, np.float32, 1)[0]
    Npart = np.fromfile(f, np.int32, 1)[0]
    x = np.fromfile(f, np.float32, Npart)
    y = np.fromfile(f, np.float32, Npart)
    z = np.fromfile(f, np.float32, Npart)
boxlen = max(max(np.max(x),np.max(y)),np.max(z))*1e-6

# Dilute particles for low-res image (redundant particles)
dil_fac = (250./pixmin)**2
Nkeep = int(len(x)/dil_fac)
indices = random.sample(range(len(x)), Nkeep)
x = x[indices]
y = y[indices]
z = z[indices]
mass *= dil_fac

# Randomly choose the centre of the periodic box
x -= np.random.rand()*boxlen*1e6
x[x<0] += boxlen*1e6
y -= np.random.rand()*boxlen*1e6
y[y<0] += boxlen*1e6
z -= np.random.rand()*boxlen*1e6
z[z<0] += boxlen*1e6

# Replicate the read-in simulation to match the appropriate size
print 'replicating box...'
if compboxlen>boxlen:
    Ncopy = int(compboxlen/boxlen) + 1
    if compboxlen>500:
        dil_fac = (compboxlen/500.)**2
        Nkeep = int(len(x)/dil_fac)
    else:
        dil_fac = 1.0
        Nkeep = len(x)
    xnew = np.empty(Nkeep*Ncopy**2)
    ynew = np.empty(Nkeep*Ncopy**2)
    print 'Need >', round(Nkeep*Ncopy**2*8*2*1e-9, 2), 'GB of RAM to produce this plot'
    count = 0
    for cx in xrange(Ncopy):
        for cy in xrange(Ncopy):
            # Give a random translation (with dilution if box size too large)
            indices = random.sample(range(len(x)), Nkeep)
            x_t = x[indices] - np.random.rand()*boxlen*1e6
            x_t[x_t<0] += boxlen*1e6
            y_t = y[indices] - np.random.rand()*boxlen*1e6
            y_t[y_t<0] += boxlen*1e6
            z_t = z[indices] - np.random.rand()*boxlen*1e6
            z_t[z_t<0] += boxlen*1e6
            
            # Centre box
            x_t -= boxlen*5e5
            y_t -= boxlen*5e5
            z_t -= boxlen*5e5
        
            # Random rotations
            theta = np.pi*np.random.random_integers(0,3)/2.
            phi = np.pi*np.random.random_integers(0,3)/2.
            psi = np.pi*np.random.random_integers(0,3)/2.
            xtemp, ytemp = x_t*np.cos(theta) + y_t*np.sin(theta), y_t*np.cos(theta) - x_t*np.sin(theta)
            y_r, ztemp = ytemp*np.cos(phi) + z_t*np.sin(phi), z_t*np.cos(phi) - ytemp*np.sin(phi)
            x_r = xtemp*np.cos(psi) - ztemp*np.sin(psi)

            # Random mirroring or uncentre
            x_r = boxlen*5e5 - x_r if np.random.random_integers(0,1) else boxlen*5e5 + x_r
            y_r = boxlen*5e5 - y_r if np.random.random_integers(0,1) else boxlen*5e5 + y_r

            xnew[count*Nkeep:(count+1)*Nkeep] = x_r + boxlen*cx*1e6
            ynew[count*Nkeep:(count+1)*Nkeep] = y_r + boxlen*cy*1e6
            count += 1

    x = xnew
    del xnew
    y = ynew
    del ynew

# extract relevant portion of the simulation
f = (x*1e-6<compboxlen) * (y*1e-6<compboxlen)
x = x[f]
y = y[f]

print 'histrogramming...'
im, xedge, yedge = np.histogram2d(x, y, bins=pix)
im *= (mass*dil_fac)
width = (xedge[1]-xedge[0])
im2plot = np.log10(im.transpose() / width**2 + 1e-10)

print 'plotting...'
fig = plt.figure()
plt.clf()
fig.subplots_adjust(left=0, bottom=0)
plt.subplot(111)

plt.imshow(im2plot, interpolation='nearest', vmin=-0.2, vmax=3.5, cmap=plt.cm.hot, extent=[0,compboxlen,0,compboxlen], zorder=4)
plt.xlabel('')
plt.xticks([])
plt.yticks([])

c_red = '#840000'
plt.plot([0,compboxlen], [0,0], c_red, lw=2, zorder=5)
plt.plot([0,0], [0,compboxlen], c_red, lw=2, zorder=5)
plt.plot([0,compboxlen], [compboxlen,compboxlen], c_red, lw=2, zorder=5)
plt.plot([compboxlen,compboxlen], [0,compboxlen], c_red, lw=2, zorder=5)

frac = 0.2
plt.plot([0,frac*compboxlen], [compboxlen, (1+0.5*frac)*compboxlen], c_red, lw=1, zorder=1)
plt.plot([compboxlen,(1+frac)*compboxlen], [compboxlen, (1+0.5*frac)*compboxlen], c_red, lw=1, zorder=1)
plt.plot([compboxlen,(1+frac)*compboxlen], [0, 0.5*frac*compboxlen], c_red, lw=1, zorder=1)
plt.plot([frac*compboxlen,(1+frac)*compboxlen], [(1+0.5*frac)*compboxlen,(1+0.5*frac)*compboxlen], c_red, lw=2, zorder=1)
plt.plot([(1+frac)*compboxlen, (1+frac)*compboxlen], [0.5*frac*compboxlen, (1+0.5*frac)*compboxlen], c_red, lw=2, zorder=1)

plt.plot([0,frac*compboxlen], [0, 0.5*frac*compboxlen], c_red, lw=1, zorder=1)
plt.plot([frac*compboxlen,frac*compboxlen], [0.5*frac*compboxlen, (1+0.5*frac)*compboxlen], c_red, lw=2, zorder=1)
plt.plot([frac*compboxlen,(1+frac)*compboxlen], [0.5*frac*compboxlen, 0.5*frac*compboxlen], c_red, lw=2, zorder=1)


# Add volume comparison of Millennium
edge = 500
thick = 7

c_Mill = '#2b8cbe'

plt.plot([0,0], [0,edge], '-', color=c_Mill, lw=1*3, zorder=6)
plt.plot([0,edge], [0,0], '-', color=c_Mill, lw=1*3, zorder=6)
plt.plot([0,edge], [edge,edge], '-', color=c_Mill, lw=1*3, zorder=6)
plt.plot([edge,edge], [edge,0.0], '-', color=c_Mill, lw=1*3, zorder=6)

plt.plot([0,frac*edge], [edge, (1+0.5*frac)*edge], c_Mill, lw=1, zorder=2)
plt.plot([edge,(1+frac)*edge], [edge, (1+0.5*frac)*edge], c_Mill, lw=1, zorder=2)
plt.plot([edge,(1+frac)*edge], [0, 0.5*frac*edge], c_Mill, lw=1, zorder=2)
plt.plot([frac*edge,(1+frac)*edge], [(1+0.5*frac)*edge,(1+0.5*frac)*edge], c_Mill, lw=2, zorder=2)
plt.plot([(1+frac)*edge, (1+frac)*edge], [0.5*frac*edge, (1+0.5*frac)*edge], c_Mill, lw=2, zorder=2)

plt.plot([0,frac*edge], [0, 0.5*frac*edge], c_Mill, lw=1, zorder=2)
plt.plot([frac*edge,frac*edge], [0.5*frac*edge, (1+0.5*frac)*edge], c_Mill, lw=2, zorder=2)
plt.plot([frac*edge,(1+frac)*edge], [0.5*frac*edge, 0.5*frac*edge], c_Mill, lw=2, zorder=2)

edge = max(edge,compboxlen)
plt.axis([-thick,edge*(1+frac)+thick,-thick,edge*(1+0.5*frac)+thick])


# Add observational volume comparison
edge = 240
c_Obs = '#a6bddb'
plt.plot([0,0], [0,edge], '--', color=c_Obs, lw=3, zorder=7)
plt.plot([0,edge], [0,0], '--', color=c_Obs, lw=3, zorder=7)
plt.plot([0,edge], [edge,edge], '--', color=c_Obs, lw=3, zorder=7)
plt.plot([edge,edge], [edge,0.0], '--', color=c_Obs, lw=3, zorder=7)

plt.plot([0,frac*edge], [edge, (1+0.5*frac)*edge], '--', color=c_Obs, lw=1, zorder=3)
plt.plot([edge,(1+frac)*edge], [edge, (1+0.5*frac)*edge], '--', color=c_Obs, lw=1, zorder=3)
plt.plot([edge,(1+frac)*edge], [0, 0.5*frac*edge], '--', color=c_Obs, lw=1, zorder=3)
plt.plot([frac*edge,(1+frac)*edge], [(1+0.5*frac)*edge,(1+0.5*frac)*edge], '--', color=c_Obs, lw=2, zorder=3)
plt.plot([(1+frac)*edge, (1+frac)*edge], [0.5*frac*edge, (1+0.5*frac)*edge], '--', color=c_Obs, lw=2, zorder=3)

plt.plot([0,frac*edge], [0, 0.5*frac*edge], '--', color=c_Obs, lw=1, zorder=3)
plt.plot([frac*edge,frac*edge], [0.5*frac*edge, (1+0.5*frac)*edge], '--', color=c_Obs, lw=2, zorder=3)
plt.plot([frac*edge,(1+frac)*edge], [0.5*frac*edge, 0.5*frac*edge], '--', color=c_Obs, lw=2, zorder=3)


# Finalise output figure
xpixplot=max(pixmin,pix)
ypixplot=max(pixmin,pix)
mydpi = np.sqrt(xpix**2 + ypix**2)/ss  # The dpi of your screen
xinplot = xpixplot*(9./7.)/mydpi
yinplot = ypixplot*(9./7.)/mydpi
fig.set_size_inches(xinplot,yinplot)
fig.set_dpi(mydpi)
imagename = 'BoxImage_'+str(compboxlen)+'.png'
plt.savefig(imagename, dpi=mydpi, bbox_inches='tight', transparent=True)
