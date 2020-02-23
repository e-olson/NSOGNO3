import numpy as np
import netCDF4 as nc
import datetime as dt
import string
import glob
import pickle
import matplotlib as mpl
from salishsea_tools import evaltools as et
from NorthNut import vvl_interp_T_to_V, vvl_interp_T_to_U;
import NorthNut as nn
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
    vlines={'N0':{'i':(122,132),'j':792},
            'N':{'i':(119,132),'j':762-1},
            'S':{'i':(119,132),'j':735-1}}
    ulines=dict()

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
    for var in ('ptrc_T','grid_V','dian_V','VT2_20','carp_T'):
        flist[var]=et.index_model_files(mod_start,mod_end,sdir,nam_fmt='nowcast',flen=1,ftype=var,tres=1)

    masks=dict()
    for ipath in vlines.keys():
        j=vlines[ipath]['j'];i0=vlines[ipath]['i'][0];i1=vlines[ipath]['i'][1]
        masks[ipath]=vmask[0,:,j,i0:i1]
    volT=dict()
    no3T=dict()
    no3TD=dict()
    no3=dict()
    for ii in range(0,len(flist['ptrc_T']['paths'])):
        if ii%10==0:
            print('day=',ii)
        with nc.Dataset(flist['ptrc_T']['paths'][ii]) as fp, \
              nc.Dataset(flist['dian_V']['paths'][ii]) as fvd, \
               nc.Dataset(flist['VT2_20']['paths'][ii]) as fvd2, \
                nc.Dataset(flist['grid_V']['paths'][ii]) as fv, \
                 nc.Dataset(flist['carp_T']['paths'][ii]) as fe:
            for ipath in vlines.keys():
                j=vlines[ipath]['j'];i0=vlines[ipath]['i'][0];i1=vlines[ipath]['i'][1]
                if not ipath in volT.keys():
                    volT[ipath]=np.zeros((24*len(flist['ptrc_T']['paths']),40,i1-i0))
                    no3T[ipath]=np.zeros((24*len(flist['ptrc_T']['paths']),40,i1-i0))
                    no3TD[ipath]=np.zeros((24*len(flist['ptrc_T']['paths']),40,i1-i0))
                    no3[ipath]=np.zeros((24*len(flist['ptrc_T']['paths']),40,i1-i0))
                    print(ii,volT.keys())
                v_x=fv.variables['vomecrty'][:,:,j,i0:i1]
                e3v=vvl_interp_T_to_V(fe.variables['e3t'],e1v,e2v,e1t,e2t,vmask,e3t_0,e3v_0)
                e3v_x=e3v[:,:,j,i0:i1]
                e1v_x=e1v[:,j,i0:i1]
                volT_x=v_x*e3v_x*e1v_x*vmask[:,:,j,i0:i1]
                no3_x=np.mean(np.ma.masked_where(tmask24[:,:,j:(j+2),i0:i1]==0,
                                   fp.variables['nitrate'][:,:,j:(j+2),i0:i1]),2)
                no3T_x=v_x*e3v_x*e1v_x*no3_x*vmask[:,:,j,i0:i1]
                volT[ipath][ii*24:(ii+1)*24,:,:]=volT_x
                no3T[ipath][ii*24:(ii+1)*24,:,:]=no3T_x
                no3[ipath][ii*24:(ii+1)*24,:,:]=no3_x
                if ipath in ('N','S'):
                    no3TD[ipath][ii*24:(ii+1)*24,:,:]=fvd.variables['NO3_VT'][:,:,j-nn.jg0,(i0-nn.ig0):(i1-nn.ig0)]
                elif ipath=='N0':
                    no3TD[ipath][ii*24:(ii+1)*24,:,:]=fvd2.variables['NO3_VT'][:,:,0,:10]
    data=dict()
    data['mod_start']=mod_start
    data['mod_end']=mod_end
    data['volT']=volT
    data['no3T']=no3T
    data['no3']=no3
    data['masks']=masks
    data['gdept_1d']=gdept_1d
    data['e3t_1d']=e3t_1d
    data['ulines']=ulines
    data['vlines']=vlines
    data['no3TD']=no3TD
    pickle.dump(data,open(savepath,'wb'))
    return

if __name__ == "__main__":
    #yr = getyr();
    #print('calcs for',yr)
    fformat0='%Y%m%d'
    #mod_start=dt.datetime(yr,1,1)
    #mod_end=dt.datetime(yr,12,31)
    mod_start=dt.datetime(2016,5,15) # originally 5/15-8/15,  but changed to even number of fortnights (6, end is included)
    mod_end=dt.datetime(2016,8,20)
    savepath='../../save/transpDiscovery2'+mod_start.strftime(fformat0)+'-'+mod_end.strftime(fformat0)+'.pkl'
    runcalcs(mod_start,mod_end,fformat0,savepath);
    print('done')
