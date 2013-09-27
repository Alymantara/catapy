import cataclysmic as cv
from pyfits import getval
import glob
from numpy import save
from scipy.ndimage.filters import gaussian_filter as gsmooth
from PyAstronomy import funcFit as fuf
from PyAstronomy import pyasl

arm='VIS'
files=glob.glob('/Users/juan/astro/SDSS1433/spectroscopy/'+arm+'/*'+arm+'.fits')
#files=glob.glob('/Users/juan/astro/V458_Vul/old_spec/vis*.fits')
t_0=53858.85689 # BJD/HJD

porb=0.054240679
#porb=0.06805555555555555
output='input_data'
wave,flux,name,mjd,ra,decl,phase,delphase,files1,hjd=[],[],[],[],[],[],[],[],[],[]
for i in files[:]:
    wav,flu=cv.read_xshooter(i)
    t1,t2,t3,t4,t5=getval(i,'OBJECT'),getval(i,'MJD-OBS'),getval(i,'RA'),getval(i,'DEC'),getval(i,'EXPTIME')
    wave.append(wav),flux.append(gsmooth(flu,0.4)),name.append(t1),mjd.append(t2),ra.append(t3),decl.append(t4),delphase.append(t5/3600/24/porb)
    hjd.append(pyasl.helio_jd(t2+0.5,218.324772,10.19017))
    phase.append(cv.phase(hjd[-1],t_0,porb))
    files1.append(i)
    print t1,' ' , t2,' ' ,hjd[-1],' ' ,phase[-1]
save(output,[wave,flux,name,hjd,phase,delphase,files1])
