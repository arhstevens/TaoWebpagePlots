from pylab import *
import sys
import os
sys.path.insert(0, '/Users/astevens/Dropbox/Swinburne Shared/6-MonthProject')
from galprops import galplot as gp
from galprops import galread as gr
from galprops import galcalc as gc

matplotlib.rcParams.update({'font.size': 1, 'xtick.major.size': 0, 'ytick.major.size': 0, 'xtick.major.width': 0, 'ytick.major.width': 0, 'ytick.minor.size': 0, 'xtick.minor.size': 0, 'axes.linewidth': 0, 'text.usetex': True, 'font.family': 'serif', 'font.serif': 'Times New Roman'})

pm = 6.17e9 # Set manually for simulation of interest [Msun/h]

pm_comp = 8.6e8 # Comparative value, set as Millennium

rcomp = (pm/pm_comp)**(1./3) # relative size of radius
pix = 150 * max(1., rcomp)

im = gc.sphere2dk(1., 2./pix, pix)
vmin, vmax = np.min(im), np.max(im)*1.4
im[im==0] = inf
plt.imshow(im, interpolation='none', cmap=plt.cm.hot, extent=[-1,1,-1,1], vmax=vmax)
gp.circle(1., 'k', lw=3)

#gp.circle(1./rcomp, 'g', lw=4, centre=[0, 1./rcomp-1])
#xedge = max(1., 1./rcomp)*1.05
#ymin = -1.05
#ymax = (2*max(1.,1./rcomp)-1)*1.05
#plt.axis([-xedge,xedge,ymin,ymax])

c = '#2b8cbe'
gp.circle(1./rcomp, c, lw=4)
xedge = max(1., 1./rcomp)*1.05
plt.axis([-xedge,xedge,-xedge,xedge])


if pm<=pm_comp:
    tstr = r'${\bf'+str(round(pm_comp/pm,2))+r' \times }$'
elif pm>pm_comp:
    tstr = r'${\bf \div '+str(round(pm/pm_comp,2))+r'}$'
#'#a6bddb'

if pm_comp/pm<10:
    plt.text(0,0, tstr, color=c, fontsize=15, ha='center', va='center')
else:
    plt.text(0,1.1, tstr, color=c, fontsize=15, ha='center', va='bottom')


plt.xticks([])
plt.yticks([])

gp.savepng('ParticleImage_'+str(int(100*rcomp)), xpixplot=pix, ypixplot=pix, transparent=True, xpix=1920, ypix=1200)
