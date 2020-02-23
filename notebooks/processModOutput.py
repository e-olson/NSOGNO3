import datetime as dt
import subprocess
import numpy as np
import os
import glob
import time
import sys

# process files for upload in archive, extracting only fields and subdomains used in paper
# to define domains, start is same as xml domain def, end is start+duration-1
# therefore T-grid N domain is -d x,112,208 -d y,644,773
# example: 
# ncks -d x,112,208 -d y,644,773 -v e3t,PAR /results2/SalishSea/nowcast-green.201812/01jan15/SalishSea_1h_20150101_20150101_carp_T.nc ./SalishSea_1h_20150101_20150101_carp_T.nc

maxproc=4
saveloc='/data/eolson/MEOPAR/SS36runs/Archive/'
t0=dt.datetime(2015,1,1)
te=dt.datetime(2017,12,31) # included
fdurHC=1 # length of each results file in days
spath='/results2/SalishSea/nowcast-green.201812/'
fstencilC='{0}/SalishSea_1h_{1}_{1}_carp_T.nc'
ostencilC='SalishSea_1h_{1}_{1}_carp_T.nc'
fmax=maxproc*10
fstencilB='{0}/SalishSea_1h_{1}_{1}_ptrc_T.nc'
ostencilB='SalishSea_1h_{1}_{1}_ptrc_T.nc'
fstencilU='{0}/SalishSea_1h_{1}_{1}_grid_U.nc'
ostencilU='SalishSea_1h_{1}_{1}_grid_U.nc'
varnamesU=('vozocrtx',)
fstencilV='{0}/SalishSea_1h_{1}_{1}_grid_V.nc'
ostencilV='SalishSea_1h_{1}_{1}_grid_V.nc'
varnamesV=('vomecrty',)
fstencilW='{0}/SalishSea_1h_{1}_{1}_grid_W.nc'
ostencilW='SalishSea_1h_{1}_{1}_grid_W.nc'
varnamesW=('vovecrtz',)
fstencilT='{0}/SalishSea_1h_{1}_{1}_grid_T.nc'
ostencilT='SalishSea_1h_{1}_{1}_grid_T.nc'
varnamesT=('vosaline','votemper')
fstencilR='{0}/SalishSea_1h_{1}_{1}_prod_T.nc'
ostencilR='SalishSea_1h_{1}_{1}_prod_T.nc'
varnamesR=('PPDIATNO3','PPPHYNO3','PPMRUBNO3','PPDIAT','PPPHY','PPMRUB')

def listFiles(fdur,fstencil,ostencil):
    ifnames=list()
    ofnames=list()
    ffmt='%Y%m%d'
    dfmt='%d%b%y'
    iits=t0
    while iits<=te:
        iite=iits+dt.timedelta(days=(fdur-1))
        iitn=iits+dt.timedelta(days=fdur)
        try:
            iifstr=glob.glob(spath+fstencil.format(iits.strftime(dfmt).lower(),iits.strftime(ffmt),iite.strftime(ffmt)),recursive=True)[0]
            ifnames.append(iifstr)
            ofnames.append(saveloc+ostencil.format(iits.strftime(dfmt).lower(),iits.strftime(ffmt),iite.strftime(ffmt)))
        except:
            print('file does not exist:  '+spath+fstencil.format(iits.strftime(dfmt).lower(),iits.strftime(ffmt),iite.strftime(ffmt)))
            raise
        iits=iitn
    if not len(ifnames)==len(ofnames):
        print('lists are not equal lenght!')
    return ifnames, ofnames

def pidswait(pids):
    for ipid in pids.keys():
        while pids[ipid].poll() is None:
            time.sleep(10)
    plist=list()
    for ipid in pids.keys():
        for line in pids[ipid].stdout:
            print(line)
        for line in pids[ipid].stderr:
            print(line)
        if pids[ipid].returncode!=0:
            print(pids[ipid].returncode)
        pids[ipid].stdout.close()
        pids[ipid].stderr.close()
        plist.append(ipid)
    for el in plist:
        pids.pop(el)
    return pids

def process_carp_T_ndom(ifnames,ofnames):
    pids=dict()
    ii=0
    for ifname,ofname in zip(ifnames,ofnames):
        pids[ii]=subprocess.Popen('ncks -d x,112,208 -d y,644,773 -v e3t,PAR '+ifname+' '+ofname, shell=True,
                                           stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
        if ii%maxproc==(maxproc-1):
            pids[ii].wait() # wait for the last one in set
        if ii%fmax==(fmax-1):
            pids=pidswait(pids)
        ii+=1
    if len(pids)>0:
        pids=pidswait(pids)
    # check that no others are still running, wait for them
    return pids

def process_ptrc_T_ndom(ifnames,ofnames):
    pids=dict()
    ii=0
    for ifname,ofname in zip(ifnames,ofnames):
        pids[ii]=subprocess.Popen('ncks -d x,112,208 -d y,644,773 -v nitrate,diatoms,flagellates,ciliates '+ifname+' '+ofname, shell=True,
                                           stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
        if ii%maxproc==(maxproc-1):
            pids[ii].wait() # wait for the last one in set
        if ii%fmax==(fmax-1):
            pids=pidswait(pids)
        ii+=1
    if len(pids)>0:
        pids=pidswait(pids)
    # check that no others are still running, wait for them
    return pids

def process_ndom(ifnames,ofnames,varnames):
    pids=dict()
    ii=0
    for ifname,ofname in zip(ifnames,ofnames):
        pids[ii]=subprocess.Popen('ncks -d x,112,208 -d y,644,773 -v '+','.join(varnames)+' '+ifname+' '+ofname, shell=True,
                                           stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
        if ii%maxproc==(maxproc-1):
            pids[ii].wait() # wait for the last one in set
        if ii%fmax==(fmax-1):
            pids=pidswait(pids)
        ii+=1
    if len(pids)>0:
        pids=pidswait(pids)
    # check that no others are still running, wait for them
    return pids

if __name__ == "__main__":
    ## uncomment to process carp_T files:
    #ifnamesC,ofnamesC=listFiles(fdurHC,fstencilC,ostencilC);
    #print('done setup')
    #pids=process_carp_T_ndom(ifnamesC,ofnamesC)
    #print('done')
    ## uncomment to process ptrc_T files:
    #ifnamesB,ofnamesB=listFiles(fdurHC,fstencilB,ostencilB);
    #print('done setup')
    #pids=process_ptrc_T_ndom(ifnamesB,ofnamesB)
    #print('done')
    ## uncomment to process velocity files:
    #ifnamesU,ofnamesU=listFiles(fdurHC,fstencilU,ostencilU);
    #pids=process_ndom(ifnamesU,ofnamesU,varnamesU);
    #print('done U')
    #ifnamesV,ofnamesV=listFiles(fdurHC,fstencilV,ostencilV);
    #pids=process_ndom(ifnamesV,ofnamesV,varnamesV);
    #print('done V')
    #ifnamesW,ofnamesW=listFiles(fdurHC,fstencilW,ostencilW);
    #pids=process_ndom(ifnamesW,ofnamesW,varnamesW);
    #print('done W')
    ifnamesR,ofnamesR=listFiles(fdurHC,fstencilR,ostencilR);
    pids=process_ndom(ifnamesR,ofnamesR,varnamesR);
    print('done R')
    ifnamesT,ofnamesT=listFiles(fdurHC,fstencilT,ostencilT);
    pids=process_ndom(ifnamesT,ofnamesT,varnamesT);
    print('done T')
