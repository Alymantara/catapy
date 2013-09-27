import pyfits as py
import numpy as n
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import cataclysmic as cv
from matplotlib import mpl,pyplot
import matplotlib.cm as cm
import sys
import os as os
import glob
import mynormalize

''' 
     --  CataPy  --
Trailed Spectra Program

Versions:
1.0 JVHS 28/02/2013
'''
lam=777
dell=3000

bins=10
xaxis='wave'
maxim=1.5
minim=0.8
cmaps=cm.Greys_r# #cm.winter_r cm.Blues#cm.gist_stern
medi=17
line_lbl='K I'
saved=n.load('input_data.npy')

if lam < min(saved[0][0]) or lam > max(saved[0][0]):
    print 'Error: input wavelength out of bounds.'
    print 'Must be between '+str(min(saved[0][0]))+' and '+str(max(saved[0][0]))+'.'
    sys.exit()
plt.ion() # Activate update


class DraggableColorbar(object):
    def __init__(self, cbar, mappable):
        self.cbar = cbar
        self.mappable = mappable
        self.press = None
        self.cycle = sorted([i for i in dir(plt.cm) if hasattr(getattr(plt.cm,i),'N')])
        self.index = self.cycle.index(cbar.get_cmap().name)

    def connect(self):
        """connect to all the events we need"""
        self.cidpress = self.cbar.patch.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.cbar.patch.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.cbar.patch.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)
        self.keypress = self.cbar.patch.figure.canvas.mpl_connect(
            'key_press_event', self.key_press)

    def on_press(self, event):
        """on button press we will see if the mouse is over us and store some data"""
        if event.inaxes != self.cbar.ax: return
        self.press = event.x, event.y

    def key_press(self, event):
        if event.key=='down':
            self.index += 1
        elif event.key=='up':
            self.index -= 1
        if self.index<0:
            self.index = len(self.cycle)
        elif self.index>=len(self.cycle):
            self.index = 0
        cmap = self.cycle[self.index]
        self.cbar.set_cmap(cmap)
        self.cbar.draw_all()
        self.mappable.set_cmap(cmap)
        self.mappable.get_axes().set_title(cmap)
        self.cbar.patch.figure.canvas.draw()

    def on_motion(self, event):
        'on motion we will move the rect if the mouse is over us'
        if self.press is None: return
        if event.inaxes != self.cbar.ax: return
        xprev, yprev = self.press
        dx = event.x - xprev
        dy = event.y - yprev
        self.press = event.x,event.y
        #print 'x0=%f, xpress=%f, event.xdata=%f, dx=%f, x0+dx=%f'%(x0, xpress, event.xdata, dx, x0+dx)
        scale = self.cbar.norm.vmax - self.cbar.norm.vmin
        perc = 0.03
        if event.button==1:
            self.cbar.norm.vmin -= (perc*scale)*n.sign(dy)
            self.cbar.norm.vmax -= (perc*scale)*n.sign(dy)
        elif event.button==3:
            self.cbar.norm.vmin -= (perc*scale)*n.sign(dy)
            self.cbar.norm.vmax += (perc*scale)*n.sign(dy)
        self.cbar.draw_all()
        self.mappable.set_norm(self.cbar.norm)
        self.cbar.patch.figure.canvas.draw()


    def on_release(self, event):
        """on release we reset the press data"""
        self.press = None
        self.mappable.set_norm(self.cbar.norm)
        self.cbar.patch.figure.canvas.draw()

    def disconnect(self):
        """disconnect all the stored connection ids"""
        self.cbar.patch.figure.canvas.mpl_disconnect(self.cidpress)
        self.cbar.patch.figure.canvas.mpl_disconnect(self.cidrelease)
        self.cbar.patch.figure.canvas.mpl_disconnect(self.cidmotion)

# Determine the sub-array of each spectra centered on the desired wave
# and +/- delta wavelength
ss=0
for i in n.arange(len(saved[0][0])-1):    
    if lam >= saved[0][0][i]   and lam <=saved[0][0][i+1]:
        ss=i
# Creates the image and sorts by phase the data
x,y = n.ogrid[0+.5/bins:2:1./bins,saved[0][0][ss-dell]:saved[0][0][ss+dell]:(saved[0][0][1]-saved[0][0][0])*1]
zvals = n.transpose(n.random.rand(len(n.transpose(y))-1,len(x))*0.)  # Put for VIS -1: -1,len(x))*0.
tt=[i[0] for i in sorted(enumerate(saved[4]), key=lambda x:x[1])] #sorts by phase
fig=plt.figure(1)
plt.clf()
plt.title('Spectra per bin. Total bins: '+str(bins))
plt.ylabel('Number of spectra per bin')
plt.xlabel('Phase Bin')
fig = pyplot.figure(1,figsize=(3,4))
hst=plt.hist(saved[4],bins,range=[0,1])
hst1=n.concatenate((n.transpose(hst[1]),n.transpose(hst[1])),axis=0)
hst2=n.concatenate((n.transpose(hst[0]),n.transpose(hst[0])),axis=0)
flag=0
'''
os.chdir('/Users/juan/astro/SDSS1433/spectroscopy/NIR')
os.system('rm phase*.txt')
os.chdir('/Users/juan/python')
'''
for i in n.arange(len(zvals)):
    tempfl=n.arange(len(saved[1][0][ss-dell:ss+dell]))*0.0
    lolo=0
    

    #fsock = open('/Users/juan/astro/SDSS1433/spectroscopy/NIR/phase_'+str(hst1[i])+'.txt', 'w')
    for j in n.arange(len(saved[4])):
                
        if hst1[i]==1.0:
            flag=1
        if saved[4][j] >=hst1[i] and saved[4][j]<hst1[i+1] and flag==0:
            lolo=lolo+1
            tempfl=cv.medfilt1(saved[1][j][ss-dell:ss+dell],medi)+tempfl
        if saved[4][j] >=hst1[i+1] and saved[4][j]<hst1[i+2] and flag==1:
            lolo=lolo+1
            tempfl=cv.medfilt1(saved[1][j][ss-dell:ss+dell],medi)+tempfl
            #print >>fsock,saved[6][j][44:-5]+'j.fits'
    #fsock.close()
    zvals[i]=tempfl/lolo/n.median(tempfl/lolo)

'''
fsock = open('/Users/juan/astro/SDSS1433/spectroscopy/NIR/phase_combine.cl', 'w')
files=glob.glob('/Users/juan/astro/SDSS1433/spectroscopy/NIR/phase_*.txt')
for i in files:
    print >>fsock,'scombine @'+i[44:]+' '+i[44:-4]+'.fits'
fsock.close()    
'''    
# %%%%%%%% NO binning, just stacked in phase order
#x,y = n.ogrid[0:2:1./len(saved[0]),saved[0][0][ss-dell]:saved[0][0][ss+dell]:(saved[0][0][1]-saved[0][0][0])*1]
#zvals = n.transpose(n.random.rand(len(n.transpose(y)),len(x))*0.)
#tt=[i[0] for i in sorted(enumerate(saved[4]), key=lambda x:x[1])]
#for i,j in zip(n.arange(len(zvals)*2),tt+tt):
#    zvals[i]=saved[1][j][ss-dell:ss+dell]


# Plots the Trailed spectra
fig = pyplot.figure(2,figsize=(8,10.5),facecolor='w')
#ax = fig.add_subplot(111)
plt.clf()
if xaxis is 'wave':
    img = plt.imshow(zvals,extent=(y.min(), y.max(),x.min(), x.max()),interpolation='nearest', cmap=cmaps,aspect='auto',origin='lower',vmin=minim,vmax=maxim)
    plt.xlabel('Wavelength, $nm$')
    limits=[y.min(), y.max(),x.min(), x.max()]
    plt.axvline(x=lam,linestyle='--',color='k')
if xaxis is 'vel':
    velo=cv.redshift1(saved[0][0][ss-dell:ss+dell],lam)
    img = plt.imshow(zvals,interpolation='nearest', cmap=cmaps,aspect='auto',origin='lower',extent=(min(velo), max(velo),x.min(), x.max()),vmin=minim,vmax=maxim)
    plt.xlabel('Velocity, [km s$^{-1}$]')
    limits=[min(velo), max(velo),x.min(), x.max()]
    plt.axvline(x=0.0,linestyle='--',color='k')
#extent=(min(velo), max(velo),x.min(), x.max())  ,vmin=0.2,vmax=1.3
#cbar = fig.colorbar(img, ticks=[minim, minim +(maxim-minim)/5.,minim +(maxim-minim)*2/5.,minim +(maxim-minim)*3/5.,minim +(maxim-minim)*4/5.,maxim ])
cbar = plt.colorbar(format='%05.2f')
cbar.set_label('Arbitrary Flux')
cbar.set_norm(mynormalize.MyNormalize(vmin=zvals.min(),vmax=zvals.max(),stretch='linear'))
cbar = DraggableColorbar(cbar,img)
cbar.connect()
plt.show()
plt.ylabel('Orbital Phase')
plt.suptitle('CataPy Trailed Spectra')
plt.title('Object: '+saved[2][0]+'\n Line: '+line_lbl+' - $\lambda_0=$'+str(lam)+' nm, $\Delta$pix='+str(dell)+', $\Delta\phi=$'+str(cv.trunc(1./bins,3))+'\n')
phis=[]
for ii in n.arange(len(zvals)/2-1):
    phis.append((hst[1][ii+1]-hst[1][ii])/2.)

n.save('phase_binned.npy',[y,zvals[:bins],x[:bins],phis,lam,dell,bins,xaxis,limits])

