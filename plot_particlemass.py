from __future__ import print_function
from pylab import *

try:
    xrange
except NameError:
    xrange = range

matplotlib.rcParams.update({'font.size': 1, 'xtick.major.size': 0, 'ytick.major.size': 0, 'xtick.major.width': 0, 'ytick.major.width': 0, 'ytick.minor.size': 0, 'xtick.minor.size': 0, 'axes.linewidth': 0, 'text.usetex': True, 'font.family': 'serif', 'font.serif': 'Times New Roman'})

# Particle mass of simulation of interest [Msun/h]
pm = 2.64e6

# Set physical size of image
pixmin = 150

# Info on your display
xpix = 1920 # horizontal screen pixels
ypix = 1200 # vertical screen pixels
ss = 15 # diagonal screen size in inches




def sphere2dk(R,Lbin,Nbin):
    # Make a square 2d kernel of a collapsed sphere.
    Nbin = int(Nbin) # Ensure an integer number of bins
    k = np.zeros((Nbin,Nbin)) # k is the convolution kernel to be output
    for i in xrange(Nbin):
        for j in xrange(Nbin):
            r = Lbin*np.sqrt((i - (Nbin-1)/2)**2 + (j - (Nbin-1)/2)**2) # Average distance of the pixel from the centre
            if r<R: k[i,j] = np.sqrt(R**2 - r**2)
    k /= np.sum(k) # Make it normalised
    return k

def circle(r,colour='white',centre=[0,0],lw=3,ls='-'):
    # Plot a circle with radius r centred on coordinates "centre"
    # "Half" indicates whether to do a full circle (f), or just the top (t) or bottom (b) halves
    r = abs(r)
    x = np.linspace(-r, r, 500)
    y1 = np.sqrt(r**2 - x**2)
    y1[0] = 0
    y1[-1] = 0
    y2 = -y1
    x, y1, y2 = x + centre[0], y1 + centre[1], y2 + centre[1]
    arg = np.where(y1==np.max(y1))[0][0]
    xp = np.concatenate((x[arg:], x[::-1], x[:arg]))
    yp = np.concatenate((y1[arg:], y2[::-1], y1[:arg]))
    plt.plot(xp, yp, linewidth=lw, color=colour, linestyle=ls)

pm_comp = 8.6e8 # Comparative particle mass, set as Millennium

rcomp = (pm/pm_comp)**(1./3) # relative size of radius
pix = pixmin * max(1., rcomp)


# Set up figure
print('plotting...')
fig = plt.figure()
plt.clf()
fig.subplots_adjust(left=0, bottom=0)
plt.subplot(111)

mydpi = np.sqrt(xpix**2 + ypix**2)/ss  # The dpi of your screen
xinplot = pix/mydpi
yinplot = pix/mydpi
fig.set_size_inches(xinplot,yinplot)
fig.set_dpi(mydpi)
lwf = 100. / mydpi


im = sphere2dk(1., 2./pix, pix)
vmin, vmax = np.min(im), np.max(im)*1.4
im[im==0] = inf
plt.imshow(im, interpolation='none', cmap=plt.cm.hot, extent=[-1,1,-1,1], vmax=vmax)
circle(1., 'k', lw=lwf*3)

c = '#2b8cbe'
circle(1./rcomp, c, lw=lwf*4)
xedge = max(1., 1./rcomp)*1.05
plt.axis([-xedge,xedge,-xedge,xedge])


if pm<=pm_comp:
    tstr = r'${\bf'+str(round(pm_comp/pm,2))+r' \times }$'
elif pm>pm_comp:
    tstr = r'${\bf \div '+str(round(pm/pm_comp,2))+r'}$'

if pm_comp/pm<10:
    plt.text(0,0, tstr, color=c, fontsize=lwf*15, ha='center', va='center')
else:
    plt.text(0,1.1, tstr, color=c, fontsize=lwf*15, ha='center', va='bottom')


plt.xticks([])
plt.yticks([])

imagename = 'ParticleImage_'+str(int(100*rcomp))+'.png'
plt.savefig(imagename, dpi=mydpi, bbox_inches='tight', transparent=True)
