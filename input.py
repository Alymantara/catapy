def input(files='*.fits',t0=1,porb=0.01,output='input_data')
    from cataclysmic import read_xshooter,phase
    from pyfits import getval
    import glob
    from numpy import save
    files=glob.glob("/Users/juan/astro/SDSS1433/spectroscopy/UVB/*.fits")
    t_0=53858.350922480226 # MJD-UTC
    porb=0.054240679
    wave,flux,name,mjd,ra,decl,phase,delphase=[],[],[],[],[],[],[],[]
    for i in files:
        wav,flu=read_xshooter(i)
        t1,t2,t3,t4,t5=getval(i,'OBJECT'),getval(i,'MJD-OBS'),getval(i,'RA'),getval(i,'DEC'),getval(i,'EXPTIME')
        wave.append(wav),flux.append(flu),name.append(t1),mjd.append(t2),ra.append(t3),decl.append(t4),delphase.append(t5/3600/24/porb)
        phase.append(phase(t2,t_0,porb))
        print t1,t2,phase[-1],t5/3600/24/porb
    save(output,[wave,flux,name,mjd,phase,delphase])
