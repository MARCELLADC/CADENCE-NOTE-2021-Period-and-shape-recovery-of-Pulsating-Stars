#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: silvioleccia
"""
import numpy as np

def meanmag_antilog(mag):
    mag=np.asarray(mag)
    flux=10.**(-mag/2.5)
    return (-2.5)*np.log10(sum(flux)/len(flux))

def ReadTeoSim(filename,dmod=0.,ebv=0.,t0_input=0.):
    
    time_model=[]
    u_model=[]
    g_model=[]
    r_model=[]
    i_model=[]
    z_model=[]
    y_model=[]
    phase_model=[]

    f=open(filename,'r')
    c=-1
    for line in f:
        c=c+1
        if c==0:
            header=line
            continue
        ll=line.split(',')
        if c==1:
            period_model=float(ll[8])*86400.
            if abs(t0_input)<1e-6:
                time_0=float(ll[0])
            else:
                time_0=t0_input
        time_model.append(float(ll[0])-time_0)
        phase_model.append((float(ll[0])-time_0)/period_model % 1)
        u_model.append(float(ll[2])+dmod+1.55607*3.1*ebv)
        g_model.append(float(ll[3])+dmod+1.18379*3.1*ebv)
        r_model.append(float(ll[4])+dmod+1.87075*3.1*ebv)
        i_model.append(float(ll[5])+dmod+0.67897*3.1*ebv)
        z_model.append(float(ll[6])+dmod+0.51683*3.1*ebv)
        y_model.append(float(ll[7])+dmod+0.42839*3.1*ebv)
    f.close()
    
    meanu=meanmag_antilog(u_model)
    meang=meanmag_antilog(g_model)
    meanr=meanmag_antilog(r_model)
    meani=meanmag_antilog(i_model)
    meanz=meanmag_antilog(z_model)
    meany=meanmag_antilog(y_model)
    amplu=max(u_model)-min(u_model)
    amplg=max(g_model)-min(g_model)
    amplr=max(r_model)-min(r_model)
    ampli=max(i_model)-min(i_model)
    amplz=max(z_model)-min(z_model)
    amply=max(y_model)-min(y_model)

    phase_model.append(1.)
    ind_0=phase_model.index(0.)
    u_model.append(u_model[ind_0])
    g_model.append(g_model[ind_0])
    r_model.append(r_model[ind_0])
    i_model.append(i_model[ind_0])
    z_model.append(z_model[ind_0])
    y_model.append(y_model[ind_0])
    
#    return time_model,phase_model,u_model,g_model,r_model,i_model,z_model,y_model
    output={'time':time_model, 'phase':phase_model,'period':period_model,
            'u':u_model, 'g': g_model, 'r': r_model, 'i': i_model, 'z': z_model, 'y': y_model,
            'meanu':meanu,'meang':meang,'meanr':meanr,'meani':meani,'meanz':meanz,'meany':meany,
            'amplu':amplu,'amplg':amplg,'amplr':amplr,'ampli':ampli,'amplz':amplz,'amply':amply}
    return output