import netCDF4 as nc
import numpy as np
import datetime as dt
from matplotlib.ticker import FormatStrFormatter
import cmocean
from salishsea_tools import evaltools as et
import pickle

ig0=112
ig1=112+97
jg0=644
jg1=644+130
fformat0='%Y%m%d'
print('NorthNut defined variables: ig0,ig1,jg0,jg1,fformat0')

with nc.Dataset('/ocean/eolson/MEOPAR/NEMO-forcing/grid/mesh_mask201702_noLPE.nc') as fmesh:
    vmask=fmesh.variables['vmask'][0,:,jg0:jg1,ig0:ig1]
    umask=fmesh.variables['umask'][0,:,jg0:jg1,ig0:ig1]
    tmask=fmesh.variables['tmask'][0,:,jg0:jg1,ig0:ig1]
    fmask=fmesh.variables['fmask'][0,:,jg0:jg1,ig0:ig1]
    vmask0=vmask[0,:,:]
    umask0=umask[0,:,:]
    gdept=fmesh.variables['gdept_0'][0,:,jg0:jg1,ig0:ig1]
    gdept_1d=fmesh.variables['gdept_1d'][0,:]
    e1t=fmesh.variables['e1t'][0,jg0:jg1,ig0:ig1]
    e2t=fmesh.variables['e2t'][0,jg0:jg1,ig0:ig1]
    e1f=fmesh.variables['e1f'][0,jg0:jg1,ig0:ig1]
    e2f=fmesh.variables['e2f'][0,jg0:jg1,ig0:ig1]
    e12t=fmesh.variables['e1t'][0,jg0:jg1,ig0:ig1].astype(float)*\
        fmesh.variables['e2t'][0,jg0:jg1,ig0:ig1].astype(float)*\
        fmesh.variables['tmask'][0,0,jg0:jg1,ig0:ig1]
    e1v=np.copy(fmesh.variables['e1v'][0,jg0:jg1,ig0:ig1]).astype(float)
    e2u=np.copy(fmesh.variables['e2u'][0,jg0:jg1,ig0:ig1]).astype(float)
    e3t_1d=fmesh.variables['e3t_1d'][0,:]
    e3t_0=fmesh.variables['e3t_0'][0,:,jg0:jg1,ig0:ig1].astype(float)
print('NorthNut defined variables: vmask, vmask0, umask, umask0, tmask, fmask, gdept, ',
        'gdept_1d, e1t, e2t, e12t, e1f, e2f, e1v, e2u, e3t_1d')

boxCol=(.7,.7,.7)
colL=(.8,.6,0)
colR=(.8,0,0.5)
arrowwidth=1
headwidth=5
headlength=3
alen=3
toff=1
apw=dict(width=arrowwidth,             #the width of the arrow in points
             headwidth=headwidth,         #the width of the base of the arrow head in points
             headlength=headlength,       #the length of the arrow head in points
             shrink=0,                    #fraction of total length to 'shrink' from both ends
             color='w')#edgecolor=boxCol,facecolor='w')  
apk=dict(width=arrowwidth,             #the width of the arrow in points
             headwidth=headwidth,         #the width of the base of the arrow head in points
             headlength=headlength,       #the length of the arrow head in points
             shrink=0,                    #fraction of total length to 'shrink' from both ends
             color='k')#edgecolor=boxCol,facecolor='w') 
apk2=dict(width=.5,             #the width of the arrow in points
                 headwidth=headwidth,         #the width of the base of the arrow head in points
                 headlength=2,       #the length of the arrow head in points
                 shrink=0,                    #fraction of total length to 'shrink' from both ends
                 color='k')#edgecolor=boxCol,facecolor='w') 
print('NorthNut defined variables: boxCol, colL, colR, arrowwidth, headwidth, headlength, alen, toff, apw, apk')


def defboxes(k):
    # calc transports: boxes in full model coords
    boxes=dict()
    boxes[0]={'i':(119,132),'j':(735,762)}
    boxes[1]={'i':(118,146),'j':(720,735)}
    boxes[2]={'i':(118,146),'j':(705,720)}
    boxes[3]={'i':(121,150),'j':(690,705)}
    boxes[4]={'i':(126,154),'j':(675,690)}
    boxes[5]={'i':(129,159),'j':(660,675)}
    #boxes[6]={'i':(130,150),'j':(645,660)}
    # boxes in subgrid coords:
    boxesDIAN=dict()
    print('volumes: ')
    for el in boxes.keys():
        boxesDIAN[el]={'i':[ii-ig0 for ii in boxes[el]['i']],'j':[ii-jg0 for ii in boxes[el]['j']]}
        vv=tmask[:k,boxesDIAN[el]['j'][0]:boxesDIAN[el]['j'][1],boxesDIAN[el]['i'][0]:boxesDIAN[el]['i'][1]]*\
                e12t[boxesDIAN[el]['j'][0]:boxesDIAN[el]['j'][1],boxesDIAN[el]['i'][0]:boxesDIAN[el]['i'][1]]*\
                e3t_0[:k,boxesDIAN[el]['j'][0]:boxesDIAN[el]['j'][1],boxesDIAN[el]['i'][0]:boxesDIAN[el]['i'][1]]
        v=np.sum(np.sum(np.sum(vv,2),1),0)
        A_north=np.sum(tmask[:k,boxesDIAN[el]['j'][1],boxesDIAN[el]['i'][0]:boxesDIAN[el]['i'][1]]*\
                e1t[boxesDIAN[el]['j'][1],boxesDIAN[el]['i'][0]:boxesDIAN[el]['i'][1]]*\
                e3t_0[:k,boxesDIAN[el]['j'][1],boxesDIAN[el]['i'][0]:boxesDIAN[el]['i'][1]])
        A_south=np.sum(tmask[:k,boxesDIAN[el]['j'][0],boxesDIAN[el]['i'][0]:boxesDIAN[el]['i'][1]]*\
                e1t[boxesDIAN[el]['j'][0],boxesDIAN[el]['i'][0]:boxesDIAN[el]['i'][1]]*\
                e3t_0[:k,boxesDIAN[el]['j'][0],boxesDIAN[el]['i'][0]:boxesDIAN[el]['i'][1]])
        A_east=np.sum(tmask[:k,boxesDIAN[el]['j'][0]:boxesDIAN[el]['j'][1],boxesDIAN[el]['i'][1]]*\
                e2t[boxesDIAN[el]['j'][0]:boxesDIAN[el]['j'][1],boxesDIAN[el]['i'][1]]*\
                e3t_0[:k,boxesDIAN[el]['j'][0]:boxesDIAN[el]['j'][1],boxesDIAN[el]['i'][1]])
        A_floor=np.sum(tmask[k,boxesDIAN[el]['j'][0]:boxesDIAN[el]['j'][1],boxesDIAN[el]['i'][0]:boxesDIAN[el]['i'][1]]*\
                e12t[boxesDIAN[el]['j'][0]:boxesDIAN[el]['j'][1],boxesDIAN[el]['i'][0]:boxesDIAN[el]['i'][1]])
        print(np.shape(tmask))
        print(el,'vol:', v,'m3')
        print(el,'north face area:',A_north/1e6,'km2')
        print(el,'south face area:',A_south/1e6,'km2')
        print(el,'east face area:',A_east/1e6,'km2')
        print(el,'floor area:',A_floor/1e6,'km2')
        print(el,'floor area:',A_floor/1e6,'km2')
    return boxes, boxesDIAN

def boxAreas(k):
    boxes, boxesDIAN = defboxes(k);
    ABoxes=dict()
    for el in boxes.keys():
        ABoxes[el]=np.sum(tmask[k,boxesDIAN[el]['j'][0]:boxesDIAN[el]['j'][1],boxesDIAN[el]['i'][0]:boxesDIAN[el]['i'][1]]*\
                e12t[boxesDIAN[el]['j'][0]:boxesDIAN[el]['j'][1],boxesDIAN[el]['i'][0]:boxesDIAN[el]['i'][1]]);
    return ABoxes

def defboxesDiscovery():
    # calc transports: boxes in full model coords
    boxes=dict()
    boxes[0]={'i':(116,130),'j':(762,771)}
    boxes[1]={'i':(119,133),'j':(753,762)}
    boxes[2]={'i':(120,132),'j':(744,753)}
    boxes[3]={'i':(119,132),'j':(735,744)}
    # boxes in subgrid coords:
    boxesDIAN=dict()
    for el in boxes.keys():
        boxesDIAN[el]={'i':[ii-ig0 for ii in boxes[el]['i']],'j':[ii-jg0 for ii in boxes[el]['j']]}
    return boxes, boxesDIAN

def boxcoordsT(box):
    xm=0.5*(box['i'][0]+box['i'][1]-1)
    ym=0.5*(box['j'][0]+box['j'][1]-1)
    x0=box['i'][0]-.5
    y0=box['j'][0]-.5
    x1=box['i'][1]-1+.5
    y1=box['j'][1]-1+.5
    return xm, ym, x0, y0, x1, y1

def boxcoordsU(box):
    xm=0.5*(box['i'][0]-1+box['i'][1]-1)
    ym=0.5*(box['j'][0]+box['j'][1]-1)
    x0=box['i'][0]-1
    y0=box['j'][0]-.5
    x1=box['i'][1]-1
    y1=box['j'][1]-1+.5
    return xm, ym, x0, y0, x1, y1

def boxcoordsV(box):
    xm=0.5*(box['i'][0]+box['i'][1]-1)
    ym=0.5*(box['j'][0]-1+box['j'][1]-1)
    x0=box['i'][0]-.5
    y0=box['j'][0]-1
    x1=box['i'][1]-1+.5
    y1=box['j'][1]-1
    return xm, ym, x0, y0, x1, y1

def makebox(boxcoords):
    iii=np.array((boxcoords[2],boxcoords[4],boxcoords[4],boxcoords[2],boxcoords[2]))
    jjj=np.array((boxcoords[3],boxcoords[3],boxcoords[5],boxcoords[5],boxcoords[3]))
    return iii,jjj

def compileZ_k(t0,te,k1,mod_basedir,mod_nam_fmt,mod_flen,fver,saveloc,recalc=False):
    savepath=saveloc+'saveZ_k_'+fver+'_kmax'+str(k1)+'_'+t0.strftime(fformat0)+\
                                                           '-'+te.strftime(fformat0)+'.pkl'
    times=[t0+dt.timedelta(hours=ii) for ii in range(0,int((te-t0).total_seconds()/3600))]
    with nc.Dataset('/ocean/eolson/MEOPAR/NEMO-forcing/grid/mesh_mask201702_noLPE.nc') as mesh:
        tmask=np.copy(mesh.variables['tmask'][:,:,:,:])
    flistC=et.index_model_files(t0,te,mod_basedir,mod_nam_fmt,mod_flen,'carp_T',1)
    Z_i=list()
    if recalc==True:
        ## calculations
        ti=t0
        for iif in range(0,len(flistC)):
            with nc.Dataset(flistC.loc[iif,['paths']].values[0]) as fC:
                Z_i.append(tmask[:,k1-1,jg0:jg1,ig0:ig1]*np.sum(fC.variables['e3t'][:,:k1,jg0:jg1,ig0:ig1],1))
        data=dict()
        data['Z']=np.concatenate(Z_i,axis=0)
        data['tmask']=tmask[:,k1-1,jg0:jg1,ig0:ig1]
        pickle.dump(data,open(savepath,'wb'))
    else:
        data=pickle.load(open(savepath,'rb'))
    Zs=data['Z']
    tmk=data['tmask']
    return Zs, tmk

def calcTransps(t0,te,k1,mod_flen,fver,saveloc,boxes=None,boxesS=None,flistV=None,flistU=None,flistW=None,flistC=None,flistT=None,
                                                    recalc=False):
    # k1 corresponds to index of w grid at bottom of cell; t grid is summed over 0 to k1-1, so total number of t grid cells included is also k1
    savepath=saveloc+'saveBoxes_'+fver+'_kmax'+str(k1)+'_'+t0.strftime(fformat0)+\
                                                           '-'+te.strftime(fformat0)+'.pkl'
    times=[t0+dt.timedelta(hours=ii) for ii in range(0,int((te-t0).total_seconds()/3600)+24)]
    if recalc==True:
        ## calculations
        NBound=dict(); SBound=dict(); EBound=dict(); BBound=dict()
        NBoundMix=dict(); SBoundMix=dict(); EBoundMix=dict(); BBoundMix=dict()
        Content=dict(); Vol=dict(); A_N=dict(); A_S=dict(); A_E=dict()

        for var in (NBound,SBound,EBound,BBound,NBoundMix,SBoundMix,EBoundMix,BBoundMix,Content,Vol,A_N,A_S,A_E):
            for bkey in boxes.keys():
                var[bkey]=np.empty((int((te-t0).days*24+24),k1))
                var[bkey].fill(np.nan)
        ti=t0
        for iif in range(0,len(flistV)):
            fNV=nc.Dataset(flistV.loc[iif,['paths']].values[0])
            fNU=nc.Dataset(flistU.loc[iif,['paths']].values[0])
            fNW=nc.Dataset(flistW.loc[iif,['paths']].values[0])
            fC=nc.Dataset(flistC.loc[iif,['paths']].values[0])
            fT=nc.Dataset(flistT.loc[iif,['paths']].values[0])
            # every file:
            for bkey in boxesS.keys(): #fill for each box
                i0=boxesS[bkey]['i'][0]
                i1=boxesS[bkey]['i'][1]
                j0=boxesS[bkey]['j'][0]
                j1=boxesS[bkey]['j'][1]
                li0=iif*mod_flen*24
                li1=(iif+1)*mod_flen*24
                NBound[bkey][li0:li1,:]=np.sum(vmask[:k1,j1-1,i0:i1]*\
                                                fNV.variables['NO3_VT'][:,:k1,j1-1,i0:i1].astype(float),2)#mmol N/s
                SBound[bkey][li0:li1,:]=np.sum(vmask[:k1,j0-1,i0:i1]*\
                                                fNV.variables['NO3_VT'][:,:k1,j0-1,i0:i1].astype(float),2)#mmol N/s
                EBound[bkey][li0:li1,:]=np.sum(umask[:k1,j0:j1,i1-1]*\
                                                fNU.variables['NO3_UT'][:,:k1,j0:j1,i1-1].astype(float),2)#mmol N/s
                BBound[bkey][li0:li1,:]=np.sum(np.sum(tmask[1:(k1+1),j0:j1,i0:i1]\
                                      *fNW.variables['NO3_WT'][:,1:(k1+1),j0:j1,i0:i1].astype(float),3),2)#mmol N/s
                NBoundMix[bkey][li0:li1,:]=np.sum(vmask[:k1,j1-1,i0:i1]*\
                                       fNV.variables['VLDFNO3'][:,:k1,j1-1,i0:i1].astype(float),2)#mmol N/s
                SBoundMix[bkey][li0:li1,:]=np.sum(vmask[:k1,j0-1,i0:i1]*\
                                       fNV.variables['VLDFNO3'][:,:k1,j0-1,i0:i1].astype(float),2)#mmol N/s
                EBoundMix[bkey][li0:li1,:]=np.sum(umask[:k1,j0:j1,i1-1]*\
                                       fNU.variables['ULDFNO3'][:,:k1,j0:j1,i1-1].astype(float),2)#mmol N/s
                BBoundMix[bkey][li0:li1,:]=np.sum(np.sum(tmask[1:(k1+1),j0:j1,i0:i1]*e12t[j0:j1,i0:i1].astype(float)\
                                       *fNW.variables['VMIXNO3'][:,1:(k1+1),j0:j1,i0:i1].astype(float),3),2)#mmol N/s
                Content[bkey][li0:li1,:]=np.sum(np.sum(tmask[:k1,j0:j1,i0:i1]*\
                                   fT.variables['nitrate'][:,:k1,(j0+jg0):(j1+jg0),(i0+ig0):(i1+ig0)].astype(float)\
                                  *fC.variables['e3t'][:,:k1,(j0+jg0):(j1+jg0),(i0+ig0):(i1+ig0)].astype(float)*\
                                            e12t[j0:j1,i0:i1].astype(float),3),2) #mmol N
                Vol[bkey][li0:li1,:]=np.sum(np.sum(tmask[:k1,j0:j1,i0:i1]\
                                  *fC.variables['e3t'][:,:k1,(j0+jg0):(j1+jg0),(i0+ig0):(i1+ig0)].astype(float)*\
                                            e12t[j0:j1,i0:i1],3),2)#m^3
                A_N[bkey][li0:li1,:]=np.sum(vmask[:k1,j1-1,i0:i1]*\
                                   np.mean(fC.variables['e3t'][:,:k1,(j1-1+jg0):(j1+1+jg0),(i0+ig0):(i1+ig0)],2)*\
                                            e1v[j1-1,i0:i1].astype(float),2)#m**2
                A_S[bkey][li0:li1,:]=np.sum(vmask[:k1,j1-1,i0:i1]*\
                                   np.mean(fC.variables['e3t'][:,:k1,(j0-1+jg0):(j0+1+jg0),(i0+ig0):(i1+ig0)],2)*\
                                            e1v[j0-1,i0:i1].astype(float),2)#m**2                            
                A_E[bkey][li0:li1,:]=np.sum(umask[:k1,j0:j1,i1-1]*\
                                   np.mean(fC.variables['e3t'][:,:k1,(j0+jg0):(j1+jg0),(i1-1+ig0):(i1+1+ig0)],3)*\
                                            e2u[j0:j1,i1-1].astype(float),2)#m**2  
            fNV.close()
            fNU.close()
            fNW.close()
            fC.close()
            fT.close()

        data=dict()
        data['NBound']=NBound
        data['SBound']=SBound
        data['EBound']=EBound
        data['BBound']=BBound
        data['NBoundMix']=NBoundMix
        data['SBoundMix']=SBoundMix
        data['EBoundMix']=EBoundMix
        data['BBoundMix']=BBoundMix
        data['Content']=Content
        data['Vol']=Vol
        data['A_N']=A_N
        data['A_S']=A_S
        data['A_E']=A_E
        data['t0']=t0
        data['te']=te
        data['boxes']=boxes
        pickle.dump(data,open(savepath,'wb'))
    else:
        data=pickle.load(open(savepath,'rb'))
        NBound=data['NBound']
        SBound=data['SBound']
        EBound=data['EBound']
        BBound=data['BBound']
        NBoundMix=data['NBoundMix']
        SBoundMix=data['SBoundMix']
        EBoundMix=data['EBoundMix']
        BBoundMix=data['BBoundMix']
        Content=data['Content']
        Vol=data['Vol']
        A_N=data['A_N']
        A_S=data['A_S']
        A_E=data['A_E']
        t0=data['t0']
        te=data['te']
        boxes=data['boxes']
    return NBound, SBound, EBound, BBound, NBoundMix, SBoundMix, EBoundMix, BBoundMix, \
            Content, Vol, A_N, A_S, A_E, times, boxes

def calcTranspsReduced(t0,te,k1,mod_flen,fver,saveloc,boxes=None,boxesS=None,flistV=None,
                       flistU=None,flistW=None,flistT=None,
                                                    recalc=False):
    # does not include flistC with e3t; remove all calcs depending on that
    # k1 corresponds to index of w grid at bottom of cell; t grid is summed over 0 to k1-1, so total number of t grid cells included is also k1
    savepath=saveloc+'saveBoxesReduced_'+fver+'_kmax'+str(k1)+'_'+t0.strftime(fformat0)+\
                                                           '-'+te.strftime(fformat0)+'.pkl'
    times=[t0+dt.timedelta(hours=ii) for ii in range(0,int((te-t0).total_seconds()/3600)+24)]
    if recalc==True:
        ## calculations
        NBound=dict(); SBound=dict(); EBound=dict(); BBound=dict()
        NBoundMix=dict(); SBoundMix=dict(); EBoundMix=dict(); BBoundMix=dict()
        #Content=dict(); Vol=dict(); A_N=dict(); A_S=dict(); A_E=dict()
        # get start and end file dates from first and last ptrc files
        tf0=flistT.loc[0,['t_0']].values[0]
        tfe=flistT.loc[len(flistT)-1,['t_n']].values[0]-dt.timedelta(days=1)
        iii0=int((t0-tf0).total_seconds()/3600)
        iiie=int((te-tf0).days*24+24)
        for var in (NBound,SBound,EBound,BBound,NBoundMix,SBoundMix,EBoundMix,BBoundMix):
            for bkey in boxes.keys():
                var[bkey]=np.empty((int((tfe-tf0).days*24+24),k1))
                var[bkey].fill(np.nan)
        ti=t0
        for iif in range(0,len(flistV)):
            fNV=nc.Dataset(flistV.loc[iif,['paths']].values[0])
            fNU=nc.Dataset(flistU.loc[iif,['paths']].values[0])
            fNW=nc.Dataset(flistW.loc[iif,['paths']].values[0])
            fT=nc.Dataset(flistT.loc[iif,['paths']].values[0])
            # every file:
            for bkey in boxesS.keys(): #fill for each box
                i0=boxesS[bkey]['i'][0]
                i1=boxesS[bkey]['i'][1]
                j0=boxesS[bkey]['j'][0]
                j1=boxesS[bkey]['j'][1]
                li0=iif*mod_flen*24
                li1=(iif+1)*mod_flen*24
                NBound[bkey][li0:li1,:]=np.sum(vmask[:k1,j1-1,i0:i1]*\
                                                fNV.variables['NO3_VT'][:,:k1,j1-1,i0:i1].astype(float),2)#mmol N/s
                SBound[bkey][li0:li1,:]=np.sum(vmask[:k1,j0-1,i0:i1]*\
                                                fNV.variables['NO3_VT'][:,:k1,j0-1,i0:i1].astype(float),2)#mmol N/s
                EBound[bkey][li0:li1,:]=np.sum(umask[:k1,j0:j1,i1-1]*\
                                                fNU.variables['NO3_UT'][:,:k1,j0:j1,i1-1].astype(float),2)#mmol N/s
                BBound[bkey][li0:li1,:]=np.sum(np.sum(tmask[1:(k1+1),j0:j1,i0:i1]\
                                      *fNW.variables['NO3_WT'][:,1:(k1+1),j0:j1,i0:i1].astype(float),3),2)#mmol N/s
                NBoundMix[bkey][li0:li1,:]=np.sum(vmask[:k1,j1-1,i0:i1]*\
                                       fNV.variables['VLDFNO3'][:,:k1,j1-1,i0:i1].astype(float),2)#mmol N/s
                SBoundMix[bkey][li0:li1,:]=np.sum(vmask[:k1,j0-1,i0:i1]*\
                                       fNV.variables['VLDFNO3'][:,:k1,j0-1,i0:i1].astype(float),2)#mmol N/s
                EBoundMix[bkey][li0:li1,:]=np.sum(umask[:k1,j0:j1,i1-1]*\
                                       fNU.variables['ULDFNO3'][:,:k1,j0:j1,i1-1].astype(float),2)#mmol N/s
                BBoundMix[bkey][li0:li1,:]=np.sum(np.sum(tmask[1:(k1+1),j0:j1,i0:i1]*e12t[j0:j1,i0:i1].astype(float)\
                                       *fNW.variables['VMIXNO3'][:,1:(k1+1),j0:j1,i0:i1].astype(float),3),2)#mmol N/s
            fNV.close()
            fNU.close()
            fNW.close()
            fT.close()
        for bkey in boxesS.keys(): #constrain to correct times
            NBound[bkey]=NBound[bkey][iii0:iiie,...]
            SBound[bkey]=SBound[bkey][iii0:iiie,...]
            EBound[bkey]=EBound[bkey][iii0:iiie,...]
            BBound[bkey]=BBound[bkey][iii0:iiie,...]
            NBoundMix[bkey]=NBoundMix[bkey][iii0:iiie,...]
            SBoundMix[bkey]=SBoundMix[bkey][iii0:iiie,...]
            EBoundMix[bkey]=EBoundMix[bkey][iii0:iiie,...]
            BBoundMix[bkey]=BBoundMix[bkey][iii0:iiie,...]
        data=dict()
        data['NBound']=NBound
        data['SBound']=SBound
        data['EBound']=EBound
        data['BBound']=BBound
        data['NBoundMix']=NBoundMix
        data['SBoundMix']=SBoundMix
        data['EBoundMix']=EBoundMix
        data['BBoundMix']=BBoundMix
        data['t0']=t0
        data['te']=te
        data['boxes']=boxes
        pickle.dump(data,open(savepath,'wb'))
    else:
        data=pickle.load(open(savepath,'rb'))
        NBound=data['NBound']
        SBound=data['SBound']
        EBound=data['EBound']
        BBound=data['BBound']
        NBoundMix=data['NBoundMix']
        SBoundMix=data['SBoundMix']
        EBoundMix=data['EBoundMix']
        BBoundMix=data['BBoundMix']
        t0=data['t0']
        te=data['te']
        boxes=data['boxes']
    return NBound, SBound, EBound, BBound, NBoundMix, SBoundMix, EBoundMix, BBoundMix, \
            times, boxes

def calcProdBoxes(t0,te,k1,mod_flen,fver,saveloc,boxes=None,boxesS=None,flistR=None,flistC=None,
                                                    recalc=False):
    # k1 corresponds to index of w grid at bottom of cell; t grid is summed over 0 to k1-1, so total number of t grid cells included is also k1
    savepath=saveloc+'saveProdBoxes_'+fver+'_kmax'+str(k1)+'_'+t0.strftime(fformat0)+\
                                                           '-'+te.strftime(fformat0)+'.pkl'
    times=[t0+dt.timedelta(hours=ii) for ii in range(0,int((te-t0).total_seconds()/3600)+24)]
    if recalc==True:
        ## calculations
        PP=dict(); NPP=dict(); 

        for var in (PP,NPP):
            for bkey in boxes.keys():
                var[bkey]=np.empty((int((te-t0).days*24+24),))
                var[bkey].fill(np.nan)
        ti=t0
        for iif in range(0,len(flistC)):
            fR=nc.Dataset(flistR.loc[iif,['paths']].values[0])
            fC=nc.Dataset(flistC.loc[iif,['paths']].values[0])
            # every file:
            for bkey in boxesS.keys(): #fill for each box
                i0=boxesS[bkey]['i'][0]
                i1=boxesS[bkey]['i'][1]
                j0=boxesS[bkey]['j'][0]
                j1=boxesS[bkey]['j'][1]
                li0=iif*mod_flen*24
                li1=(iif+1)*mod_flen*24
                NPP[bkey][li0:li1]=np.sum(np.sum(np.sum(tmask[:k1,j0:j1,i0:i1]*\
                                   (fR.variables['PPDIATNO3'][:,:k1,(j0+jg0):(j1+jg0),(i0+ig0):(i1+ig0)].astype(float)\
                                    +fR.variables['PPPHYNO3'][:,:k1,(j0+jg0):(j1+jg0),(i0+ig0):(i1+ig0)].astype(float)\
                                    +fR.variables['PPMRUBNO3'][:,:k1,(j0+jg0):(j1+jg0),(i0+ig0):(i1+ig0)].astype(float))\
                                  *fC.variables['e3t'][:,:k1,(j0+jg0):(j1+jg0),(i0+ig0):(i1+ig0)].astype(float)*\
                                            e12t[j0:j1,i0:i1].astype(float),3),2),1)*1e-3 #mol N/s
                PP[bkey][li0:li1]=np.sum(np.sum(np.sum(tmask[:k1,j0:j1,i0:i1]*\
                                   (fR.variables['PPDIAT'][:,:k1,(j0+jg0):(j1+jg0),(i0+ig0):(i1+ig0)].astype(float)\
                                    +fR.variables['PPPHY'][:,:k1,(j0+jg0):(j1+jg0),(i0+ig0):(i1+ig0)].astype(float)\
                                    +fR.variables['PPMRUB'][:,:k1,(j0+jg0):(j1+jg0),(i0+ig0):(i1+ig0)].astype(float))\
                                  *fC.variables['e3t'][:,:k1,(j0+jg0):(j1+jg0),(i0+ig0):(i1+ig0)].astype(float)*\
                                            e12t[j0:j1,i0:i1].astype(float),3),2),1)*1e-3 #mol N/s                
            fR.close()
            fC.close()

        data=dict()
        data['PP']=PP
        data['NPP']=NPP
        data['t0']=t0
        data['te']=te
        data['boxes']=boxes
        pickle.dump(data,open(savepath,'wb'))
    else:
        data=pickle.load(open(savepath,'rb'))
        PP=data['PP']
        NPP=data['NPP']
        t0=data['t0']
        te=data['te']
        boxes=data['boxes']
    return PP, NPP, times, boxes

def calcProd(t0,te,fver,saveloc,fliste3t=None,flistPP=None,recalc=True):
    savepath=saveloc+'IPP_'+fver+'_'+t0.strftime(fformat0)+'-'+te.strftime(fformat0)+'.pkl'
    times=[t0+dt.timedelta(hours=ii) for ii in range(0,int((te-t0).total_seconds()/3600)+24)]
    if recalc==True:
        ## calculations
        tf0=flistPP.loc[0,['t_0']].values[0]
        tfe=flistPP.loc[len(flistPP)-1,['t_n']].values[0]-dt.timedelta(days=1)
        IPPx=np.empty((int((tfe-tf0).days*24+24),jg1-jg0,ig1-ig0))
        INPPx=np.empty((int((tfe-tf0).days*24+24),jg1-jg0,ig1-ig0))
        iii0=int((t0-tf0).total_seconds()/3600)
        iiie=int((te-tf0).days*24+24)
        ti=t0
        for iif in range(0,len(flistPP)):
            print(iif)
            li0=iif*mod_flen*24
            li1=(iif+1)*mod_flen*24
            with nc.Dataset(flistPP.loc[iif,['paths']].values[0]) as fPP, \
                    nc.Dataset(fliste3t.loc[iif,['paths']].values[0]) as fe3t:
                IPPx[li0:li1,...]=np.sum(tmask2*\
                        (fPP.variables['PPDIAT'][:,:,jg0:jg1,ig0:ig1]+fPP.variables['PPPHY'][:,:,jg0:jg1,ig0:ig1]+\
                         fPP.variables['PPMRUB'][:,:,jg0:jg1,ig0:ig1])*\
                                     fe3t.variables['e3t'][:,:,jg0:jg1,ig0:ig1],1)
                INPPx[li0:li1,...]=np.sum(tmask2*\
                        (fPP.variables['PPDIATNO3'][:,:,jg0:jg1,ig0:ig1]+fPP.variables['PPPHYNO3'][:,:,jg0:jg1,ig0:ig1]+\
                         fPP.variables['PPMRUBNO3'][:,:,jg0:jg1,ig0:ig1])*\
                                     fe3t.variables['e3t'][:,:,jg0:jg1,ig0:ig1],1)
        #constrain to correct times
        IPP=IPPx[iii0:iiie,...]
        INPP=INPPx[iii0:iiie,...]
        data=dict()
        data['IPP']=IPP
        data['INPP']=INPP
        pickle.dump(data,open(savepath,'wb'))
    else:
        data=pickle.load(open(savepath,'rb'))
        IPP=data['IPP']
        INPP=data['INPP']
    return IPP, INPP

def calcFluxFields(flistV,flistU,flistW,flistC,k,saveloc,t0,te,fver,recalc=True):
    savepath=saveloc+'saveFields_'+fver+'_k'+str(k)+'_'+t0.strftime(fformat0)+\
                                                           '-'+te.strftime(fformat0)+'.pkl'
    if recalc==True:
        Ts=dict()
        Ms=dict()
        Ts['V']=np.zeros((len(flistV)*24,jg1-jg0,ig1-ig0))
        Ms['V']=np.zeros((len(flistV)*24,jg1-jg0,ig1-ig0))
        for iif in range(0,len(flistV)):
            with nc.Dataset(flistV.loc[iif,['paths']].values[0]) as fi, \
                    nc.Dataset(flistC.loc[iif,['paths']].values[0]) as fc:
                ie13v=np.sum(0.5*(fc.variables['e3t'][:,:k,jg0:jg1,ig0:ig1]+\
                                  fc.variables['e3t'][:,:k,(jg0+1):(jg1+1),ig0:ig1])*e1v,1)
                Ts['V'][(iif*24):(iif*24+24),:,:]=np.sum(fi.variables['NO3_VT'][:,:k,:,:],1)/ie13v
                Ms['V'][(iif*24):(iif*24+24),:,:]=np.sum(fi.variables['VLDFNO3'][:,:k,:,:],1)/ie13v
        Ts['U']=np.zeros((len(flistU)*24,jg1-jg0,ig1-ig0))
        Ms['U']=np.zeros((len(flistU)*24,jg1-jg0,ig1-ig0))
        for iif in range(0,len(flistU)):
            with nc.Dataset(flistU.loc[iif,['paths']].values[0]) as fi, \
                    nc.Dataset(flistC.loc[iif,['paths']].values[0]) as fc:
                ie23u=np.sum(0.5*(fc.variables['e3t'][:,:k,jg0:jg1,ig0:ig1]+\
                                  fc.variables['e3t'][:,:k,jg0:jg1,(ig0+1):(ig1+1)])*e2u,1)
                Ts['U'][(iif*24):(iif*24+24),:,:]=np.sum(fi.variables['NO3_UT'][:,:k,:,:],1)/ie23u
                Ms['U'][(iif*24):(iif*24+24),:,:]=np.sum(fi.variables['ULDFNO3'][:,:k,:,:],1)/ie23u
        Ts['W']=np.zeros((len(flistW)*24,jg1-jg0,ig1-ig0))
        Ms['W']=np.zeros((len(flistW)*24,jg1-jg0,ig1-ig0))
        for iif in range(0,len(flistW)):
            with nc.Dataset(flistW.loc[iif,['paths']].values[0]) as fi:
                Ts['W'][(iif*24):(iif*24+24),:,:]=fi.variables['NO3_WT'][:,k,:,:]/e12t
                Ms['W'][(iif*24):(iif*24+24),:,:]=fi.variables['VMIXNO3'][:,k,:,:] # already /m2
        meanT=dict()
        meanM=dict()
        meanT['V']=np.mean(Ts['V'],0)
        meanM['V']=np.mean(Ms['V'],0)
        meanT['U']=np.mean(Ts['U'],0)
        meanM['U']=np.mean(Ms['U'],0)
        meanT['W']=np.mean(Ts['W'],0)
        meanM['W']=np.mean(Ms['W'],0)
        data=dict()
        data['meanT']=meanT
        data['meanM']=meanM
        pickle.dump(data,open(savepath,'wb'))
    else:
        data=pickle.load(open(savepath,'rb'))
        meanT=data['meanT']
        meanM=data['meanM']
    return meanT, meanM

def transpConversions(boxes,NBound,SBound,EBound,BBound,NBoundMix,SBoundMix,EBoundMix,BBoundMix,k):
    # units: BBound: mmol N/s, at hourly res;
    # mmol N/s * 24*3600 s/day * 1e-3 mmol/mol * 1e-6 Mmol/mol
    #kk=10 # to 10 m
    print('units now mol/s')
    NBoundC=dict(); EBoundC=dict(); SBoundC=dict(); BBoundC=dict()
    NBoundMixC=dict(); EBoundMixC=dict(); SBoundMixC=dict(); BBoundMixC=dict()
    for el in boxes.keys():
        NBoundC[el]=np.mean(np.sum(NBound[el][:,:k],1))*1e-3 # mol N/s
        SBoundC[el]=np.mean(np.sum(SBound[el][:,:k],1))*1e-3
        EBoundC[el]=np.mean(np.sum(EBound[el][:,:k],1))*1e-3
        BBoundC[el]=np.mean(BBound[el][:,k-1])*1e-3
        NBoundMixC[el]=np.mean(np.sum(NBoundMix[el][:,:k],1))*1e-3
        SBoundMixC[el]=np.mean(np.sum(SBoundMix[el][:,:k],1))*1e-3
        EBoundMixC[el]=np.mean(np.sum(EBoundMix[el][:,:k],1))*1e-3
        BBoundMixC[el]=np.mean(BBoundMix[el][:,k-1])*1e-3
        print(el)
        print(NBoundC[el],BBoundC[el],EBoundC[el])
        print(NBoundMixC[el],BBoundMixC[el],EBoundMixC[el])
    return NBoundC, SBoundC, EBoundC, BBoundC, NBoundMixC, SBoundMixC, EBoundMixC, BBoundMixC

def drawboxesT(ax,boxes,boxCol):
    for el in boxes.keys():
        iii,jjj=makebox(boxcoordsT(boxes[el]))
        ax.plot(iii-ig0,jjj-jg0,'-',color=boxCol,linewidth=1)
def drawboxesU(ax,boxes,boxCol):
    for el in boxes.keys():
        iii,jjj=makebox(boxcoordsU(boxes[el]))
        ax.plot(iii-ig0,jjj-jg0,'-',color=boxCol,linewidth=1)
def drawboxesV(ax,boxes,boxCol):
    for el in boxes.keys():
        iii,jjj=makebox(boxcoordsV(boxes[el]))
        ax.plot(iii-ig0,jjj-jg0,'-',color=boxCol,linewidth=1)

def annotYTranspUpper(iax,boxes,NBoundC,SBoundC,NBoundMixC,SBoundMixC,prec=1):
    xm=dict(); ym=dict(); x0=dict(); y0=dict(); x1=dict(); y1=dict()
    for el in boxes.keys():
        xm[el], ym[el], x0[el], y0[el], x1[el], y1[el]=boxcoordsV(boxes[el])
    for el in (0,1,2,3,4,5):
        val=SBoundC[el]+SBoundMixC[el]
        if not val==0:
            dl=val/np.abs(val)
            ar=iax.annotate("", xy=(xm[el]-ig0, y0[el]+dl*alen-jg0), 
                            xytext=(xm[el]-ig0, y0[el]-dl*alen-jg0),
                            arrowprops=apk)
        an=iax.annotate(str(np.round(np.abs(val),prec)), xy=(x0[el]-ig0-toff, y0[el]-jg0),
                    fontsize=8,verticalalignment='center',ha='right',color='k')
    for el in (0,):
        val=NBoundC[el]+NBoundMixC[el]
        if not val==0:
            dl=val/np.abs(val)
            ar=iax.annotate("", xy=(xm[el]-ig0, y1[el]+dl*alen-jg0), 
                            xytext=(xm[el]-ig0, y1[el]-dl*alen-jg0),
                            arrowprops=apk)
        an=iax.annotate(str(np.round(np.abs(val),prec)), xy=(x0[el]-ig0-toff, y1[el]-jg0),
                    fontsize=8,verticalalignment='center',ha='right',color='k')
    for el in (1,):
        val=NBoundC[1]-SBoundC[0]+NBoundMixC[1]-SBoundMixC[0]
        if not val==0:
            dl=val/np.abs(val)
            ar=iax.annotate("", xy=(.5*(x1[0]+x1[1])-ig0, y1[el]+dl*alen-jg0), 
                            xytext=(.5*(x1[0]+x1[1])-ig0, y1[el]-dl*alen-jg0),
                            arrowprops=apk)
        an=iax.annotate(str(np.round(np.abs(val),prec)), xy=(x1[1]-ig0+toff+3, y1[1]-jg0),
                    fontsize=8,verticalalignment='center',ha='left',color='k')
    for el in (2,3,4,5):
        if x1[el]>x1[el-1]:
            val=NBoundC[el]-SBoundC[el-1]+NBoundMixC[el]-SBoundMixC[el-1]
            if not val==0:
                dl=val/np.abs(val)
                xp=(x1[el]+x1[el-1])/2+.3
                ar=iax.annotate("", xy=(xp-ig0, y1[el]+dl*alen-jg0), 
                                xytext=(xp-ig0, y1[el]-dl*alen-jg0),
                                arrowprops=apk2)
            an=iax.annotate(str(np.round(np.abs(val),prec)), xy=(x1[el]-ig0+toff+3, y1[el]-jg0),
                        fontsize=8,verticalalignment='center',ha='left',color='k')
    # highlight lines
    for el in (1,2,3,4,5):
        iax.plot((x0[el-1]-ig0,x1[el-1]-ig0),(y1[el]-jg0,y1[el]-jg0),'-',color=colL,linewidth=1)
        iax.plot((x1[el-1]-ig0,x1[el]-ig0),(y1[el]-jg0,y1[el]-jg0),'-',color=colR,linewidth=1)
    for el in (0,):
        iax.plot((x0[el]-ig0,x1[el]-ig0),(y1[el]-jg0,y1[el]-jg0),'-',color=colL,linewidth=1)
    for el in (5,):
        iax.plot((x0[el]-ig0,x1[el]-ig0),(y0[el]-jg0,y0[el]-jg0),'-',color=colL,linewidth=1)
    return

def annotXTranspUpper(iax,boxes,EBoundC,EBoundMixC,prec=1):
    xm=dict(); ym=dict(); x0=dict(); y0=dict(); x1=dict(); y1=dict()
    for el in boxes.keys():
        xm[el], ym[el], x0[el], y0[el], x1[el], y1[el]=boxcoordsU(boxes[el])
    for el in (1,2,3,4,5):
        val=EBoundC[el]+EBoundMixC[el]
        if not val==0:
            dl=val/np.abs(val)
            iax.annotate("", xy=(x1[el]+dl*alen-ig0, ym[el]-jg0), 
                         xytext=(x1[el]-dl*alen-ig0, ym[el]-jg0),
                         arrowprops=apk) 
        else:
            dl=3
        iax.annotate(str(np.round(np.abs(val),prec)), xy=(x1[el]-ig0+dl+toff+5, ym[el]-jg0),fontsize=8,
                     ha='left',va='center',color='k')
    ## highlight lines
    #for el in (1,2,3,4,5):
    #    col=colL if el%2==1 else colR 
    #    iax.plot((x1[el]-ig0,x1[el]-ig0),(y0[el]-jg0,y1[el]-jg0),'-',color=col,linewidth=1)
    return

def annotWTTranspUpper(iax,boxes,BBoundC,prec=1):
    xm=dict(); ym=dict(); x0=dict(); y0=dict(); x1=dict(); y1=dict()
    for el in boxes.keys():
        xm[el], ym[el], x0[el], y0[el], x1[el], y1[el]=boxcoordsT(boxes[el])
    for el in (0,1,2,3,4,5):
        val=BBoundC[el]
        #dl=val/np.abs(val)
        iax.annotate(str(np.round(val,prec)), xy=(x0[el]-ig0-toff+1.5, ym[el]-jg0),fontsize=8,
                     ha='right',va='center',color='k')
    return

def annotWMTranspUpper(iax,boxes,BBoundMixC,prec=1):
    xm=dict(); ym=dict(); x0=dict(); y0=dict(); x1=dict(); y1=dict()
    for el in boxes.keys():
        xm[el], ym[el], x0[el], y0[el], x1[el], y1[el]=boxcoordsT(boxes[el])
    for el in (0,1,2,3,4,5):
        val=BBoundMixC[el]
        #dl=val/np.abs(val)
        iax.annotate(str(np.round(val,prec)), xy=(x0[el]-ig0-toff, ym[el]-jg0),fontsize=8,
                     ha='right',va='center',color='k')
    return

def annotPPUpper(iax,boxes,PP,prec=1):
    xm=dict(); ym=dict(); x0=dict(); y0=dict(); x1=dict(); y1=dict()
    for el in boxes.keys():
        xm[el], ym[el], x0[el], y0[el], x1[el], y1[el]=boxcoordsT(boxes[el])
    for el in (0,1,2,3,4,5):
        val=PP[el]
        #dl=val/np.abs(val)
        iax.annotate(str(np.round(val,prec)), xy=(x1[el]-ig0+toff, ym[el]-jg0),fontsize=8,
                     ha='left',va='center',color='k')
    return

def vvl_interp_T_to_W(pe3_in,e3w_0,e3t_0,tmask):
    pe3_out=np.zeros(np.shape(pe3_in))
    pe3_out[:,0,:,:] = e3w_0[:,0,:,:] + pe3_in[:,0,:,:] - e3t_0[:,0,:,:]
    pe3_out[:,1:,:,:] = e3w_0[:,1:,:,:]+(1-.5*tmask[:,1:,:,:])*(pe3_in[:,:-1,:,:]-e3t_0[:,:-1,:,:])+.5*tmask[:,1:,:,:]*(pe3_in[:,1:,:,:]-e3t_0[:,1:,:,:])
    #     ! - ML - The use of mask in this formaula enables the special treatment of the last w- point without indirect adressing
    #     DO jk = 2, jpk
    #        pe3_out(:,:,jk) = e3w_0(:,:,jk) + ( 1.0_wp - 0.5_wp * tmask(:,:,jk) ) * ( pe3_in(:,:,jk-1) - e3t_0(:,:,jk-1) )   &
    #           &                            +            0.5_wp * tmask(:,:,jk)   * ( pe3_in(:,:,jk  ) - e3t_0(:,:,jk  ) )
    #     END DO
    return pe3_out

def vvl_interp_T_to_V(pe3_in,e1v,e2v,e1t,e2t,vmask,e3t_0,e3v_0):
    # keep extra dimensions in mesh variables
    e12v = e1v[0,:,:]*e2v[0,:,:]
    e12t = e1t[0,:,:]*e2t[0,:,:]
    pe3_out=np.zeros(np.shape(pe3_in))
    pe3_out[:,:,:-1,:-1]=.5*vmask[:,:,:-1,:-1]/e12v[:-1,:-1]*(e12t[:-1,:-1]*(pe3_in[:,:,:-1,:-1]-e3t_0[:,:,:-1,:-1])+e12t[1:,:-1]*(pe3_in[:,:,1:,:-1]-e3t_0[:,:,1:,:-1]))
    pe3_out=pe3_out+e3v_0
    #DO jk = 1, jpk
    #        DO jj = 1, jpjm1
    #           DO ji = 1, fs_jpim1   ! vector opt.
    #              pe3_out(ji,jj,jk) = 0.5_wp * vmask(ji,jj,jk) * r1_e12v(ji,jj)                                   &
    #                 &                       * (   e12t(ji,jj  ) * ( pe3_in(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) )     &
    #                 &                           + e12t(ji,jj+1) * ( pe3_in(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) )
    #           END DO
    #        END DO
    #     END DO
    #pe3_out(:,:,:) = pe3_out(:,:,:) + e3v_0(:,:,:)
    return pe3_out

def vvl_interp_T_to_U(pe3_in,e1u,e2u,e1t,e2t,umask,e3t_0,e3u_0):
    e12u = e1u[0,:,:]*e2u[0,:,:]
    e12t = e1t[0,:,:]*e2t[0,:,:]
    pe3_out=np.zeros(np.shape(pe3_in))
    pe3_out[:,:,:-1,:-1]=.5*umask[:,:,:-1,:-1]/e12u[:-1,:-1]*(e12t[:-1,:-1]*(pe3_in[:,:,:-1,:-1]-e3t_0[:,:,:-1,:-1])+e12t[1:,:-1]*(pe3_in[:,:,1:,:-1]-e3t_0[:,:,1:,:-1]))
    pe3_out=pe3_out+e3u_0
    #! horizontal surface weighted interpolation
    #DO jk = 1, jpk
    #    DO jj = 1, jpjm1
    #       DO ji = 1, fs_jpim1   ! vector opt.
    #          pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk) * r1_e12u(ji,jj)                                   &
    #             &                       * (   e12t(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) )     &
    #             &                           + e12t(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) )
    #       END DO
    #    END DO
    #END DO
    #pe3_out(:,:,:) = pe3_out(:,:,:) + e3u_0(:,:,:)
    return pe3_out

def vvl_deptw(e3t,e3w,tmask):
    #wmask=tmask
    dept=np.zeros(np.shape(e3w))
    dept[:,0,...]=0.5*e3w[:,0,...]
    depw=np.zeros(np.shape(e3w))
    zcoef=0*tmask #tmask-wmask
    for jk in range(1,np.shape(e3w)[1]):
        depw[:,jk,...]=depw[:,jk-1,...]+e3t[:,jk-1,...]
        dept[:,jk,...]=zcoef[:,jk,...]*(depw[:,jk,...]+.5*e3w[:,jk,...])\
                  +(1-zcoef[:,jk,...])*(dept[:,jk-1,...]+e3w[:,jk,...])
    mask=np.where(tmask==1,tmask,np.nan)
    return dept*mask, depw*mask

     # fsdept_n(:,:,1) = 0.5_wp * fse3w_n(:,:,1)
     # fsdepw_n(:,:,1) = 0.0_wp
     # fsde3w_n(:,:,1) = fsdept_n(:,:,1) - sshn(:,:)
     # DO jk = 2, jpk
     #    DO jj = 1,jpj
     #       DO ji = 1,jpi
     #         !    zcoef = (tmask(ji,jj,jk) - wmask(ji,jj,jk))   ! 0 everywhere tmask = wmask, ie everywhere expect at jk = mikt
     #                                                            ! 1 for jk = mikt
     #          zcoef = (tmask(ji,jj,jk) - wmask(ji,jj,jk))
     #          fsdepw_n(ji,jj,jk) = fsdepw_n(ji,jj,jk-1) + fse3t_n(ji,jj,jk-1)
     #          fsdept_n(ji,jj,jk) =      zcoef  * ( fsdepw_n(ji,jj,jk  ) + 0.5 * fse3w_n(ji,jj,jk))  &
     #              &                + (1-zcoef) * ( fsdept_n(ji,jj,jk-1) +       fse3w_n(ji,jj,jk)) 
     #          fsde3w_n(ji,jj,jk) = fsdept_n(ji,jj,jk) - sshn(ji,jj)
     #       END DO
     #    END DO
     # END DO




