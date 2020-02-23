import numpy as np
import netCDF4 as nc
import datetime as dt
import string
import glob
import pickle
import matplotlib as mpl
from salishsea_tools import evaltools as et
from NorthNut import vvl_interp_T_to_V, vvl_interp_T_to_U;
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
    ulines={'SJDFK':{'i':140,'j':(253,308)},
            'SJDFW':{'i':20,'j':(360,438)},
            'SJDFE':{'i':186,'j':(228,315)},
            'Malaspina':{'i':187,'j':(722,746)}}
    vlines={'Haro':{'i':(213,246),'j':305},
            'SJC':{'i':(258,270),'j':281},
            'Rosario':{'i':(280,311),'j':265},
            'Discovery':{'i':(120,131),'j':737},
            'Sutil':{'i':(138,164),'j':749}}

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
    for var in ('ptrc_T','grid_U','grid_V','carp_T'):
        flist[var]=et.index_model_files(mod_start,mod_end,sdir,nam_fmt='nowcast',flen=1,ftype=var,tres=1)

    masks=dict()
    for ipath in vlines.keys():
        j=vlines[ipath]['j'];i0=vlines[ipath]['i'][0];i1=vlines[ipath]['i'][1]
        masks[ipath]=vmask[0,:,j,i0:i1]
    for ipath in ulines.keys():
        i=ulines[ipath]['i'];j0=ulines[ipath]['j'][0];j1=ulines[ipath]['j'][1]
        masks[ipath]=umask[0,:,j0:j1,i]
    volT=dict()
    no3T=dict()
    no3=dict()
    for ii in range(0,len(flist['ptrc_T']['paths'])):
        if ii%10==0:
            print('day=',ii)
        with nc.Dataset(flist['ptrc_T']['paths'][ii]) as fp, \
              nc.Dataset(flist['grid_U']['paths'][ii]) as fu, \
               nc.Dataset(flist['grid_V']['paths'][ii]) as fv, \
                nc.Dataset(flist['carp_T']['paths'][ii]) as fe:
            for ipath in vlines.keys():
                j=vlines[ipath]['j'];i0=vlines[ipath]['i'][0];i1=vlines[ipath]['i'][1]
                if not ipath in volT.keys():
                    volT[ipath]=np.zeros((24*len(flist['ptrc_T']['paths']),40,i1-i0))
                    no3T[ipath]=np.zeros((24*len(flist['ptrc_T']['paths']),40,i1-i0))
                    no3[ipath]=np.zeros((24*len(flist['ptrc_T']['paths']),40,i1-i0))
                    print(ii,volT.keys())
                v_x=fv.variables['vomecrty'][:,:,j,i0:i1]
                e3v=vvl_interp_T_to_V(fe.variables['e3t'],e1v,e2v,e1t,e2t,vmask,e3t_0,e3v_0)
                e3v_x=e3v[:,:,j,i0:i1]
                e1v_x=e1v[:,j,i0:i1]
                volT_x=v_x*e3v_x*e1v_x*vmask[:,:,j,i0:i1]
                no3_x=np.mean(np.ma.masked_where(tmask24[:,:,j:(j+2),i0:i1]==0,fp.variables['nitrate'][:,:,j:(j+2),i0:i1]),2)
                no3T_x=v_x*e3v_x*e1v_x*no3_x*vmask[:,:,j,i0:i1]
                volT[ipath][ii*24:(ii+1)*24,:,:]=volT_x
                no3T[ipath][ii*24:(ii+1)*24,:,:]=no3T_x
                no3[ipath][ii*24:(ii+1)*24,:,:]=no3_x
            for ipath in ulines.keys():
                i=ulines[ipath]['i'];j0=ulines[ipath]['j'][0];j1=ulines[ipath]['j'][1]
                if not ipath in volT.keys():
                    volT[ipath]=np.zeros((24*len(flist['ptrc_T']['paths']),40,j1-j0))
                    no3T[ipath]=np.zeros((24*len(flist['ptrc_T']['paths']),40,j1-j0))
                    no3[ipath]=np.zeros((24*len(flist['ptrc_T']['paths']),40,j1-j0))
                    print(ii,volT.keys())
                u_x=fu.variables['vozocrtx'][:,:,j0:j1,i]
                e3u=vvl_interp_T_to_U(fe.variables['e3t'],e1u,e2u,e1t,e2t,umask,e3t_0,e3u_0)
                e3u_x=e3u[:,:,j0:j1,i]
                e2u_x=e2u[:,j0:j1,i]
                volT_x=u_x*e3u_x*e2u_x*umask[:,:,j0:j1,i]
                no3_x=np.mean(np.ma.masked_where(tmask24[:,:,j0:j1,i:(i+2)]==0,fp.variables['nitrate'][:,:,j0:j1,i:(i+2)]),3)
                no3T_x=u_x*e3u_x*e2u_x*no3_x*umask[:,:,j0:j1,i]
                volT[ipath][ii*24:(ii+1)*24,:,:]=volT_x
                no3T[ipath][ii*24:(ii+1)*24,:,:]=no3T_x
                no3[ipath][ii*24:(ii+1)*24,:,:]=no3_x
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
    pickle.dump(data,open(savepath,'wb'))
    return

if __name__ == "__main__":
    yr = getyr();
    print('calcs for',yr)
    fformat0='%Y%m%d'
    mod_start=dt.datetime(yr,1,1)
    mod_end=dt.datetime(yr,12,31)
    savepath='../../save/transpLines'+mod_start.strftime(fformat0)+'-'+mod_end.strftime(fformat0)+'.pkl'
    runcalcs(mod_start,mod_end,fformat0,savepath);
    print('done')
