import numpy as np
import netCDF4 as nc
import datetime as dt
import string
import glob
import pickle
from salishsea_tools import evaltools as et
import sys

def getyr():
    try:
        yr=int(sys.argv[1])
    except:
        print(' no year specified')
        raise
    return yr

def runcalcs(t0,te,saveloc):
    SOGtmaskPath='../../save/SOGtmask.pkl'
    (tmaskSOG,ig0,ig1,jg0,jg1)=pickle.load(open(SOGtmaskPath,'rb'))
    SOGmask=tmaskSOG[:,:,jg0:jg1,ig0:ig1]
    with nc.Dataset('/ocean/eolson/MEOPAR/NEMO-forcing/grid/mesh_mask201702_noLPE.nc') as fm:
        A=np.expand_dims(fm.variables['e1t'][:,jg0:jg1,ig0:ig1]*fm.variables['e2t'][:,jg0:jg1,ig0:ig1],0)
    mod_basedir='/data/eolson/results/MEOPAR/SS36runs/linkHC201812/'
    fformat0='%Y%m%d'
    mod_nam_fmt='nowcast'
    mod_flen=1
    fver='HC201812'
    fliste3t=et.index_model_files(t0,te,mod_basedir,mod_nam_fmt,mod_flen,'carp_T',1)
    flistPP=et.index_model_files(t0,te,mod_basedir,mod_nam_fmt,mod_flen,'prod_T',1)
    flistB=et.index_model_files(t0,te,mod_basedir,mod_nam_fmt,mod_flen,'ptrc_T',1)
    flistT=et.index_model_files(t0,te,mod_basedir,mod_nam_fmt,mod_flen,'grid_T',1)
    content=dict()
    tr={'NO3':'nitrate', 'NH4':'ammonium', 'SIL':'silicon', 'DIAT':'diatoms', 
        'FLAG':'flagellates', 'MRUB':'ciliates', 'MICZ':'microzooplankton', 
        'DON':'dissolved_organic_nitrogen', 'PON':'particulate_organic_nitrogen', 
        'BSI':'biogenic_silicon'}
    for el in tr.keys():
        content[el]=np.empty((int((te-t0).days+1),))
    rates=dict()
    for el in ('PP','NPP','NHtoNO','PPDIAT','SIDIS'):
        rates[el]=np.empty((int((te-t0).days+1),))
    zz_remin_NH = 4.0e-7
    zz_remin_D_bSi = 2.78e-6
    SOGmask24=np.tile(SOGmask,(24,1,1,1))
    times=[t0+dt.timedelta(days=ii) for ii in range(0,int((te-t0).total_seconds()/3600/24)+1)]
    ## calculations
    for iif in range(0,len(flistPP)):
        print(iif,dt.datetime.now()) #progress report
        li0=iif*mod_flen
        li1=(iif+1)*mod_flen
        with nc.Dataset(flistPP.loc[iif,['paths']].values[0]) as fPP, \
              nc.Dataset(flistT.loc[iif,['paths']].values[0]) as fT, \
               nc.Dataset(flistB.loc[iif,['paths']].values[0]) as fB, \
                nc.Dataset(fliste3t.loc[iif,['paths']].values[0]) as fe3t:
            for ili0 in range(0,mod_flen): # mod_flen is number of days in file                
                # content
                for el in tr.keys():
                    content[el][ili0+li0]=np.sum(1e-3*SOGmask*fB.variables[tr[el]][0,:,jg0:jg1,ig0:ig1]\
                                            *A*fe3t.variables['e3t'][0,:,jg0:jg1,ig0:ig1])
                    # mol/mmol * mmol/m3*A*e3t = mol
                # rates
                rates['PP'][ili0+li0]=np.sum(1e-3*3600*SOGmask*(\
                                                  fPP.variables['PPDIAT'][ili0*24:(ili0+1)*24,:,jg0:jg1,ig0:ig1]+\
                                                  fPP.variables['PPPHY'][ili0*24:(ili0+1)*24,:,jg0:jg1,ig0:ig1]+\
                                                  fPP.variables['PPMRUB'][ili0*24:(ili0+1)*24,:,jg0:jg1,ig0:ig1])*\
                                     A*fe3t.variables['e3t'][ili0*24:(ili0+1)*24,:,jg0:jg1,ig0:ig1])
                    # mmol/m3/s *m3 * 3600 s/hr * 1e-3 mol/mmol = mol/hr
                rates['NPP'][ili0+li0]=np.sum(1e-3*3600*SOGmask*(\
                                                  fPP.variables['PPDIATNO3'][ili0*24:(ili0+1)*24,:,jg0:jg1,ig0:ig1]+\
                                                  fPP.variables['PPPHYNO3'][ili0*24:(ili0+1)*24,:,jg0:jg1,ig0:ig1]+\
                                                  fPP.variables['PPMRUBNO3'][ili0*24:(ili0+1)*24,:,jg0:jg1,ig0:ig1])*\
                                     A*fe3t.variables['e3t'][ili0*24:(ili0+1)*24,:,jg0:jg1,ig0:ig1])
                rates['PPDIAT'][ili0+li0]=np.sum(1e-3*3600*SOGmask*(fPP.variables['PPDIAT'][ili0*24:(ili0+1)*24,:,jg0:jg1,ig0:ig1])*\
                                     A*fe3t.variables['e3t'][ili0*24:(ili0+1)*24,:,jg0:jg1,ig0:ig1])
                TQ10=np.where(SOGmask24==1,
                            np.exp(0.07*(fT.variables['votemper'][ili0*24:(ili0+1)*24,:,jg0:jg1,ig0:ig1]-20.0)),0)
                rates['NHtoNO'][ili0+li0]=np.sum(1e-3*3600*zz_remin_NH*SOGmask*\
                                                  fB.variables['ammonium'][ili0*24:(ili0+1)*24,:,jg0:jg1,ig0:ig1]*\
                                                  fB.variables['ammonium'][ili0*24:(ili0+1)*24,:,jg0:jg1,ig0:ig1]*\
                                             TQ10*A*fe3t.variables['e3t'][ili0*24:(ili0+1)*24,:,jg0:jg1,ig0:ig1])
                rates['SIDIS'][ili0+li0]=np.sum(1e-3*3600*zz_remin_D_bSi*SOGmask*\
                                              fB.variables['biogenic_silicon'][ili0*24:(ili0+1)*24,:,jg0:jg1,ig0:ig1]*\
                                             TQ10*A*fe3t.variables['e3t'][ili0*24:(ili0+1)*24,:,jg0:jg1,ig0:ig1])
    times=np.array(times)
    contentPath=saveloc+'SOGcontent'+t0.strftime(fformat0)+'-'+te.strftime(fformat0)+'.pkl'
    pickle.dump((times,content),open(contentPath,'wb'))
    ratesPath=saveloc+'SOGrates'+t0.strftime(fformat0)+'-'+te.strftime(fformat0)+'.pkl'
    pickle.dump((times,rates),open(ratesPath,'wb'))
    return

if __name__ == "__main__":
    yr = getyr();
    print('calcs for',yr)
    t0=dt.datetime(yr,1,1)
    te=dt.datetime(yr+1,1,1)# include next day to facilitate derivative estimates
    saveloc='../../save/'
    runcalcs(t0,te,saveloc);
    print('done')
