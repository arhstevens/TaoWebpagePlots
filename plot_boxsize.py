from pylab import *
import sys
import os
sys.path.insert(0, '/Users/astevens/Dropbox/Swinburne Shared/6-MonthProject')
from galprops import galplot as gp
from galprops import galread as gr
from galprops import galcalc as gc
import pynbody
from scipy import signal as ss
import random

matplotlib.rcParams.update({'font.size': 1, 'xtick.major.size': 0, 'ytick.major.size': 0, 'xtick.major.width': 0, 'ytick.major.width': 0, 'ytick.minor.size': 0, 'xtick.minor.size': 0, 'axes.linewidth': 0, 'text.usetex': True, 'font.family': 'serif', 'font.serif': 'Times New Roman'})



compboxlen = 1000 # Set manually -- length of box this image is supposed to be for

pixmin = 250
pix = pixmin * compboxlen/500.




# Inform pynbody of how many processor the sim used
pynbody.ramses.multiprocess_num = 8
pynbody.config['number_of_threads'] = 8

f = pynbody.load('/Users/astevens/ramses/trunk/ramses/testcases/cosmo/output_00008')
redshift = 1.0

# Get particle information for DM
x = np.array(f.d['pos'][:,0].in_units('pc'), dtype=np.float32) * (1+redshift)
y = np.array(f.d['pos'][:,1].in_units('pc'), dtype=np.float32) * (1+redshift)
z = np.array(f.d['pos'][:,2].in_units('pc'), dtype=np.float32) * (1+redshift)
mass = f.d['mass'].in_units('Msol')[0]

#boxlen = (np.max(x)-np.min(x)) * 1e-6
boxlen = max(max(np.max(x),np.max(y)),np.max(z))*1e-6

# Dilute particles
dil_fac = 5. * (250./pixmin)**2
Nkeep = int(len(x)/dil_fac)
indices = random.sample(range(len(x)), Nkeep)
x = x[indices]
y = y[indices]
z = z[indices]
mass *= dil_fac

#print 'lengths', len(x), len(y), len(z)
#print 'mins', np.min(x), np.min(y), np.min(z)
#print 'maxs', np.max(x), np.max(y), np.max(z)
#
# Randomly choose the centre of the periodic box
x -= np.random.rand()*boxlen*1e6
x[x<0] += boxlen*1e6
y -= np.random.rand()*boxlen*1e6
y[y<0] += boxlen*1e6
z -= np.random.rand()*boxlen*1e6
z[z<0] += boxlen*1e6

#print 'lengths', len(x), len(y), len(z)
#print 'mins', np.min(x), np.min(y), np.min(z)
#print 'maxs', np.max(x), np.max(y), np.max(z)
#
#
#print len(x), len(y), len(z)


#if compboxlen>boxlen:
#    Ncopy = int(compboxlen/boxlen)
#    print 'Ncopy, sqr', Ncopy, Ncopy*Ncopy
#    x0, y0, z0 = 1.0*x, 1.0*y, 1.0*z
#    for cx in xrange(Ncopy+1):
#        for cy in xrange(Ncopy+1):
#            if cx==0 and cy==0: continue
#            print 'cx, cy', cx, cy
#            x = np.append(x, x0+boxlen*cx*1e6)
#            y = np.append(y, y0+boxlen*cy*1e6)
#            z = np.append(z, z0)

if compboxlen>boxlen:
    Ncopy = int(compboxlen/boxlen) + 1
    if compboxlen>500:
        dil_fac = (compboxlen/500.)**2
        print 'second dilution factor', dil_fac
        Nkeep = int(len(x)/dil_fac)
        print 'Nkeep', Nkeep
    else:
        dil_fac = 1.0
        Nkeep = len(x)
    xnew = np.empty(Nkeep*Ncopy**2)
    ynew = np.empty(Nkeep*Ncopy**2)
#    znew = np.empty(Ncopy**2*len(x))
    print 'copying box', Ncopy**2, 'times'
    print 'leads to', Nkeep*Ncopy**2*8*2 * 1e-9, 'GB for positions'
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
#            znew[count*Nkeep:(count+1)*Nkeep] = z
            count += 1

    print 'finished doing copies'
    x = xnew
    del xnew
    y = ynew
    del ynew
#    z = znew
#    del znew
else:
    # extract relevant portion of the simulation
    f = (x*1e-6<compboxlen) * (y*1e-6<compboxlen) * (z*1e-6<compboxlen)
    x = x[f]
    y = y[f]
#    z = z[f]
print 'finished extraction'

im, xedge, yedge = np.histogram2d(x, y, bins=pix)
im *= (mass*dil_fac)
width = (xedge[1]-xedge[0])
print 'finished histrogramming'

rsmooth = width / (2**14)
Lpix = width/pix
if Lpix < rsmooth:
    k = gc.sphere2dk(rsmooth, Lpix, 2*rsmooth/Lpix)
    im = ss.convolve2d(im,k,mode='same')
im2plot = np.log10(im.transpose() / width**2 + 1e-10)
print 'finished building im2plot'

gp.figure()
#Ncopy = int(compboxlen/boxlen)+1
#for cx in xrange(Ncopy):
#    for cy in xrange(Ncopy+1):
#        plt.imshow(im2plot, interpolation='nearest', vmin=-0.2, vmax=3.5, cmap=plt.cm.hot, extent=[boxlen*cx,boxlen*(cx+1),boxlen*cy,boxlen*(cy+1)])
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



edge = 500
thick = 7

c_Mill = '#2b8cbe'

plt.plot([0,0], [0,edge], '-', color=c_Mill, lw=1*3, zorder=6)
plt.plot([0,edge], [0,0], '-', color=c_Mill, lw=1*3, zorder=6)
plt.plot([0,edge], [edge,edge], '-', color=c_Mill, lw=1*3, zorder=6)
plt.plot([edge,edge], [edge,0.0], '-', color=c_Mill, lw=1*3, zorder=6)

#plt.fill_between([-thick,edge+thick], [0,0], [-thick,-thick], color=c_Mill)
#plt.fill_between([-thick,edge+thick], [edge+thick,edge+thick], [edge,edge], color=c_Mill)
#plt.fill_betweenx([-thick,edge+thick], [0,0], [-thick,-thick], color=c_Mill)
#plt.fill_betweenx([-thick,edge+thick], [edge+thick,edge+thick], [edge,edge], color=c_Mill)

plt.plot([0,frac*edge], [edge, (1+0.5*frac)*edge], c_Mill, lw=1, zorder=2)
plt.plot([edge,(1+frac)*edge], [edge, (1+0.5*frac)*edge], c_Mill, lw=1, zorder=2)
plt.plot([edge,(1+frac)*edge], [0, 0.5*frac*edge], c_Mill, lw=1, zorder=2)
plt.plot([frac*edge,(1+frac)*edge], [(1+0.5*frac)*edge,(1+0.5*frac)*edge], c_Mill, lw=2, zorder=2)
plt.plot([(1+frac)*edge, (1+frac)*edge], [0.5*frac*edge, (1+0.5*frac)*edge], c_Mill, lw=2, zorder=2)

plt.plot([0,frac*edge], [0, 0.5*frac*edge], c_Mill, lw=1, zorder=2)
plt.plot([frac*edge,frac*edge], [0.5*frac*edge, (1+0.5*frac)*edge], c_Mill, lw=2, zorder=2)
plt.plot([frac*edge,(1+frac)*edge], [0.5*frac*edge, 0.5*frac*edge], c_Mill, lw=2, zorder=2)


edge = max(edge,compboxlen)
#plt.axis([edge-pedge,edge,edge-pedge,pedge])
plt.axis([-thick,edge*(1+frac)+thick,-thick,edge*(1+0.5*frac)+thick])

edge = 240 # Now setting for observational volume comparison
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


gp.savepng('BoxImage_'+str(compboxlen), xpixplot=max(pixmin,pix), ypixplot=max(pixmin,pix), transparent=True, xpix=1920, ypix=1200)

# ADDITIONAL FEATURES COULD INCLUDE:
# ADDING A RANDOM ORIENTATION FOR THE IMAGE
# CHOOSING A RANDOM CENTRE FOR THE BOX
# COPYING THE BOX, USING ITS ACTUAL SIZE, TO MORE ACCURATELY REPRESENT THE VOLUME OF INTEREST