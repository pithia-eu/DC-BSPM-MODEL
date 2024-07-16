#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sept 18 2020

@author: edithb
"""

import pandas as pd
import numpy as np
import datetime
from dateutil import relativedelta
from artificial_points_allMLT import calcMLT

cols=['01','02','03','04','05','06','07','08','09','10','11','12']

def fill_pp_holes_int(equatorial_plasmapause_cartesian,k_slide,year,month,day,user_path_out,time_plot,date_data):
    
    print('looking for points if big holes ...')
    
    global df,dfi,cols
    
    df=pd.DataFrame()
    dfpo=pd.DataFrame(equatorial_plasmapause_cartesian,columns=cols)
#    for column in dfp.columns: dfp = dfp[dfp[column]!=-2000000000000000000000.000].reset_index(drop=True)

    j = np.where(dfpo['01'] < -1.0E21)    
    dfp=dfpo.loc[j[0][1]+1:j[0][2]-1]

    xp,yp = dfp['01'],dfp['02']
    dfp['MLTp'] = calcMLT(xp,yp)
    dfp=dfp.sort_values(by=['MLTp'])
    
#  Look for a big big hole !!!!!
    if list(dfp.MLTp.loc[(dfp.MLTp>19)&(dfp.MLTp<23.9)]) == []:  # check if there are no points
        try:
            print('big hole found, looking for ancient plasmapause ...')
            
            date1=str(year) + str(month).zfill(2) + str(day).zfill(2) + str(int(k_slide/2)).zfill(2) + '00'
            date1=datetime.datetime.strptime(date1,'%Y%m%d%H%M')
            date_1 = (date1 - relativedelta.relativedelta(hours=3)).strftime("%Y%m%d")        
            time_1 = (date1 - relativedelta.relativedelta(hours=3)).strftime("%Y-%m-%d_%Hh%Mm")
            print('read 3 hours before : ',time_1)        ##### Put a try here!!!!!!
            date_data_1=date_1[:4]+'/'+date_1[4:6].zfill(2)+'/'+date_1[6:8].zfill(2)+'/'
    
            df=pd.read_csv(user_path_out+date_data_1+'PP_'+time_1+'.csv')        
            xp,yp = df['01'],df['02']        
            df['MLTp']= calcMLT(xp,yp) 
            df=df.loc[(df.MLTp>17)&(df.MLTp<23.9)]  # take points from 3 hours earlier
            df=df.sort_values(by=['MLTp']).drop(columns=['MLTp'])
            if list(df.index) != []:
                df=create_points(df)
                dfp=dfp.drop(columns=['MLTp'])
                dfp=pd.concat([dfp,df]).reset_index(drop=True)
                
                dfpn = pd.DataFrame()
                dfpn = pd.concat([dfpn,dfpo.loc[0:j[0][1]]])  
                dfpn = pd.concat([dfpn,dfp])
                dfpn = pd.concat([dfpn,dfpo.loc[j[0][2]:j[0][3]]])
    
                dfpn.to_csv(user_path_out+date_data+'PP_'+time_plot+'.csv', sep=',',columns=cols,float_format='%.4f',index=False) # save updated PP
                equatorial_plasmapause_cartesian = np.array(dfpn)
       
            else: print('file from last hour does not contain points') 
        except: print('CARE, CARE: ancient plasmapause not available!!!')
                     
    else:
        print('plasmapause hole is not so big')
    
    return equatorial_plasmapause_cartesian

#################################################################################################################
def create_points(df):  
    
    p = np.poly1d(np.polyfit(df['01'],df['02'], 3))
    xi = np.sort(np.random.uniform(low=df['01'].min(), high=df['01'].max(), size=50))
    yi = p(xi)
    
    ri = np.zeros([xi.size]); ri = np.sqrt((xi)**2 + (yi)**2)

    dfi=pd.DataFrame(columns=cols)
    for i in range(len(xi)):
        dfi.loc[i]=[xi[i],yi[i]]+9*[0.1]+[ri[i]]
    
    return dfi

#################################################################################################################