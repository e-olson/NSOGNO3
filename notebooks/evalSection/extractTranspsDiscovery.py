import numpy as np
import netCDF4 as nc
import datetime as dt
import string
import glob
import pickle
import matplotlib as mpl
import NorthNut as nn
from salishsea_tools import evaltools as et
import sys

def getyr():
    try:
        yr=int(sys.argv[1])
    except:
        print(' no year specified')
        raise
    return yr

def runcalcs(mod_start,mod_end,fformat0,savepath):
    sdir='/data/eolson/results/MEOPAR/SS36runs/linkHC201812/'

    with nc.Dataset('/ocean/eolson/MEOPAR/NEMO-forcing/grid/mesh_mask201702_noLPE.nc') as fm:
        tmask=np.copy(fm.variables['tmask'])
        umask=np.copy(fm.variables['umask'])
        vmask=np.copy(fm.variables['vmask'])
        navlon=np.copy(fm.variables['nav_lon'])
        navlat=np.copy(fm.variables['nav_lat'])
        e3t_0=np.copy(fm.variables['e3t_0'])
        e3u_0=np.copy(fm.variables['e3u_0'])
        e3v_0=np.copy(fm.variables['e3v_0'])
        e1t=np.copy(fm.variables['e1t'])
        e2t=np.copy(fm.variables['e2t'])
        e1v=np.copy(fm.variables['e1v'])
        e2v=np.copy(fm.variables['e2v'])
        e1u=np.copy(fm.variables['e1u'])
        e2u=np.copy(fm.variables['e2u'])
        gdept_1d=fm.variables['gdept_1d'][0,:]
        e3t_1d=fm.variables['e3t_1d'][0,:]
    tmask24=np.tile(tmask,(24,1,1,1))
    flist=dict()
    for var in ('dian_V','VT2_20'):
        flist[var]=et.index_model_files(mod_start,mod_end,sdir,nam_fmt='nowcast',flen=1,ftype=var,tres=1)
    masks=dict()
    masks['N']=vmask[0,:,792,122:132]
    masks['S']=nn.vmask[:,93,9:18]
    no3T=dict()
    no3T['N']=np.zeros((24*len(flist['dian_V']['paths']),40,10))
    no3T['S']=np.zeros((24*len(flist['dian_V']['paths']),40,9))
    for ii in range(0,len(flist['dian_V']['paths'])):
        if ii%10==0:
            print('day=',ii)
        with nc.Dataset(flist['dian_V']['paths'][ii]) as fd, \
              nc.Dataset(flist['VT2_20']['paths'][ii]) as ft:
            no3T['N'][ii*24:(ii+1)*24,:,:]=ft.variables['NO3_VT'][:,:,0,:10]
            no3T['S'][ii*24:(ii+1)*24,:,:]=fd.variables['NO3_VT'][:,:,93,9:18]
    data=dict()
    data['mod_start']=mod_start
    data['mod_end']=mod_end
    data['no3T']=no3T
    data['masks']=masks
    data['gdept_1d']=gdept_1d
    data['e3t_1d']=e3t_1d
    data['jis']={'N':(792,(122,132)),'S':(nn.jg0+93,(nn.ig0+9,nn.ig0+18))}
    pickle.dump(data,open(savepath,'wb'))
    return

if __name__ == "__main__":
    yr = getyr();
    print('calcs for',yr)
    fformat0='%Y%m%d'
    mod_start=dt.datetime(yr,1,1)
    mod_end=dt.datetime(yr,12,31)
    savepath='../../save/transpDiscovery'+mod_start.strftime(fformat0)+'-'+mod_end.strftime(fformat0)+'.pkl'
    runcalcs(mod_start,mod_end,fformat0,savepath);
    print('done')
