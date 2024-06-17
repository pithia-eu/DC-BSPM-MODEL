#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 07:51:44 2020

@author: edithb
"""

from dateutil import relativedelta
import datetime
import pandas as pd
import numpy as np

def Kp(kpdata,year,month,day,user_path_in):    
    
    date = str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)
    d_1 = (datetime.datetime(year,month,day)- relativedelta.relativedelta(days=1))
    date_1 = str(d_1.year)+'-'+str(d_1.month).zfill(2)+'-'+str(d_1.day).zfill(2)
    url='https://iswa.gsfc.nasa.gov/IswaSystemWebApp/hapi/'
    
    if kpdata=='NOAA' : 
        file='data?id=NOAA_KP_P3H&time.min='+date_1+'T00:00:00.0Z&time.max='+date+'T23:59:59.0Z'
        columns=['time','Kp']
        
    else :
        file='data?id=gfz_obs_geo_3hour_indices&time.min='+date_1+'T00:00:00.0Z&time.max='+date+'T23:59:59.0Z'
        columns=['time','Kp','other']
        
    df=pd.read_csv(url+file,header=None,names=columns)
    
    dfh=pd.DataFrame((list(np.arange(0,22,3)))*2)
    dfd=pd.DataFrame([1]*8+[2]*8)
    
    df['H']=dfh
    df['D']=dfd

    print('preparing Kp input file ...')

    red_date = rem0(str(day))+'-'+rem0(str(month))+'-'+str(year)
    
    with open(user_path_in+"kp" + red_date + ".dat", 'w') as f:  # write a text file to be read by fortran codes
        space='       '
        for i,r in df.iterrows():
            f.write(space+str(r.D)+space+str(r.H)+space+str(r.Kp)+space+'1.00000'+space+'2.00000\n')

        f.write('   -9999      0      0.00000  1.00000  2.00000\n')
           
    return 

#################################################################################################################
def rem0 (a):
    if a[0]=='0': a=str(a[1])
    return a

#################################################################################################################