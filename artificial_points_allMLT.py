#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 11:38:33 2020

@author: edithb
"""

import pandas as pd
import numpy as np

cols=['01','02','03','04','05','06','07','08','09','10','11','12']
import itertools as it
from scipy.interpolate import Rbf
from scipy import interpolate as inter
import numpy.polynomial.polynomial as poly



def calcMLT(x,y):
    x=np.array(x) ; y=np.array(y)
    Lon = np.rad2deg(np.arctan2(-y, -x)) 
    Lon = np.where(Lon<0,Lon+360,Lon) #   if Lon < 0  : Lon = Lon + 360.  
    MLT = np.where(Lon<180 , 12. + Lon * (24./360.) , (Lon-180.) * (24./360.)) 
    #if Lon < 180: MLT = 12. + Lon * (24./360.)
    #else:         MLT = (Lon-180.) * (24./360.)    
    return MLT


def interp(x1,y1,x2,y2,xj,yj,mlt,dfp):

    coefs = poly.polyfit([x1,x2]+xj,[y1,y2]+yj, 1)
#    coefs = poly.polyfit([x1,x2],[y1,y2], 2)
    p = poly.Polynomial(coefs)   
#    p = inter.interp1d(np.array([x1,x2]+xj),np.array([y1,y2]+yj))

    xi = np.sort(np.random.uniform(low=np.min([x1,x2]), high=np.max([x1,x2]), size=10))
    yi = p(xi)
    ri = np.zeros([xi.size])
    ri = np.sqrt((xi)**2 + (yi)**2)
    MLTi = calcMLT(xi,yi)    
#    print(MLTi)
    
    dfi=pd.DataFrame(columns=cols+['MLTp'])
    for i in range(len(xi)):
        dfi.loc[i]=[xi[i],yi[i]]+9*[0.1]+[ri[i],MLTi[i]] # to fill the 12 columns !!!!! remove TODO!!!!       
#    print(dfi)

    dfp=pd.concat([dfp,dfi]).reset_index(drop=True)
    print('interpolated points added')    
#    print(dfi)
        
    return  dfp



def fill_pp_holes_allMLT(equatorial_plasmapause_cartesian):
    print('looking for points if necessary ...')
#    import numpy.polynomial.polynomial as poly
    poor=[]
    dfpo=pd.DataFrame(equatorial_plasmapause_cartesian,columns=cols)
    
    j = np.where(dfpo['01'] < -1.0E21)    
    dfp=dfpo.loc[j[0][1]+1:j[0][2]-1]
    
#    dfp =dfp.loc[dfp['01']!=-2000000000000000000000.000].reset_index(drop=True)
    dfp =dfp.loc[np.abs(dfp['01']) < 12].reset_index(drop=True) ## CARE!!!!
    
    xp,yp = dfp['01'],dfp['02']
    dfp['MLTp'] = calcMLT(xp,yp)   
#    dfp=dfp.sort_values(by=['MLTp'])
    
    for mlt in it.chain(range(12,24),range(0,12)):  # check if there are points in all MLT
        
        df=pd.DataFrame()
        if list(dfp.MLTp.loc[dfp.MLTp.astype(int)==mlt]) == []:  # check if there are no points 
            print('hole found in MLT = ', mlt,', try interpolating with previous points ...')

            # take the last point from the last filled mlt 
            ang=[]; mlt_t=[]
            for i in list(it.chain(range(12,24),range(0,12))):
                if mlt==0 : mlt_minus=23
                else:       mlt_minus=mlt-1
                if i==mlt_minus:   # take the previous mlt (it should have points), CARE could be MLT=11 and rather distorted!!!! CHECK
                    mlt_t.append(i)
                    mlt_1=i
                    break
                
# patch: hole found at MLT=12 (20150107)
#            
            if list(dfp.MLTp.loc[dfp.MLTp.astype(int)==mlt_1]) == []: # CARE CARE CARE              
                df_1 = dfp.loc[dfp.MLTp.astype(int)==mlt-2]
            else:
                df_1 = dfp.loc[dfp.MLTp.astype(int)==mlt_1]
# patch                
            df = pd.concat([df,df_1])
            ang.append(15*mlt)
            mlt_t.append(mlt)
            
            if mlt_1>12 : x1,y1=df_1['01'].max(),df_1.loc[df_1['01'].idxmax()]['02']
            else:         x1,y1=df_1['01'].min(),df_1.loc[df_1['01'].idxmin()]['02']
            
            positions=it.cycle(list(np.arange(24)))
            starting_at_mlt1 = it.islice(positions, mlt+1, None)
            for i in range(24):
                mltnext=next(starting_at_mlt1) # over all empty mlt
                
                if list(dfp.MLTp.loc[dfp.MLTp.astype(int)==mltnext]) != []:  # next non-empty
                    df1 = dfp.loc[dfp.MLTp.astype(int)==mltnext]
                    mlt_t.append(mltnext)
                    break
                else: 
                    ang.append(15*mltnext)
                    mlt_t.append(mltnext)
            
            df = pd.concat([df,df1])
                        
            if mltnext<12 : x2,y2=df1['01'].max(),df1.loc[df1['01'].idxmax()]['02']
            else:           x2,y2=df1['01'].min(),df1.loc[df1['01'].idxmin()]['02']
                        
            r1 = np.sqrt((x1)**2 + (y1)**2)
            r2 = np.sqrt((x2)**2 + (y2)**2)
#            r = np.min([r1,r2])
#            r = (r1+r2)/2
            
            Dr = (r1-r2)/len(ang)
            
            xj=[]; yj=[] ; r=r1 # create one point for each empty mlt with r-Dr from x1 to x2
            for angu in ang:
                r = r-Dr
                xj.append(r*np.cos(np.deg2rad(angu+0.1)))
                yj.append(r*np.sin(np.deg2rad(angu)))
                
            if len(xj)>1:
                for k in range(len(xj)):
                    dfp = interp(x1,y1,xj[k],yj[k],[],[],mlt,dfp)  # interp mlt's one by one
                    x1=xj[k] ; y1=yj[k]

            else: dfp = interp(x1,y1,x2,y2,xj,yj,mlt,dfp)                
                    
        else:
            NP=len(list(dfp.MLTp.loc[dfp.MLTp.astype(int)==mlt]))
            print('MLT = ',mlt,' contains ',NP,' points')
            if NP < 3: poor.append([mlt,NP])               
        
    print(poor)
    
    dfp=dfp.drop(columns=['MLTp'])
    
    dfpn = pd.DataFrame()
    dfpn = pd.concat([dfpn,dfpo.loc[0:j[0][1]]])  
    dfpn = pd.concat([dfpn,dfp])
    dfpn = pd.concat([dfpn,dfpo.loc[j[0][2]:j[0][3]]])
    
    equatorial_plasmapause_cartesian = np.array(dfpn)
                
    return equatorial_plasmapause_cartesian

#################################################################################################################
