#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 07:51:44 2020

@author: edithb
"""

from dateutil import relativedelta
import datetime
import urllib.request as urllib
import pandas as pd

kp={}

def Kp_new(year,month,day,user_path_in):
    date=str(year) + str(month).zfill(2) + str(day).zfill(2)
    site_today='ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-nowcast-archive/tab/'
    site = 'ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/tab/'
    kp_date = 'kp'+date[2:6]+'.tab'
    date_1 = str(int((datetime.datetime.strptime(date,'%Y%m%d').date() - relativedelta.relativedelta(days=1)).strftime('%Y%m%d')))
#    date_2 = str(int((datetime.datetime.strptime(date,'%Y%m%d').date() - relativedelta.relativedelta(days=2)).strftime('%Y%m%d')))

    if datetime.datetime.strptime(date,'%Y%m%d').date()==datetime.datetime.now().date(): 
        url = site_today+kp_date
    else: url = site+kp_date
    print('Kp from : ',url)

#    date_2_found=False    
    date_1_found=False    
    date_found=False
    lines = urllib.urlopen(url).read().decode('utf-8').splitlines()
    k=len(lines)
    for i in range(k):
        line=' '.join(lines[i].split(' ')).split()
        # try:
        #     if line[0]==date_2[2:8]: date_2_found=True ; store_Kp (line)
        # except:pass
        try:
            if line[0]==date_1[2:8]: date_1_found=True ; store_Kp (line) #; print('date-1 :',line)
        except:pass
        try:
            if line[0]==date[2:8]: date_found=True ; store_Kp (line) #; print('date :',line)
        except:pass
            
    if date_found: print(date, ' data found')

    if date_1_found: print(date_1, ' data found')
    else: 
        kp_date = 'kp'+date_1[2:6]+'.tab'
        url = site+kp_date
        lines = urllib.urlopen(url).read().decode('utf-8').splitlines()
        k=len(lines)
        for i in range(k):
            line=' '.join(lines[i].split(' ')).split()
            if line[0]==date_1[2:8]: 
                date_1_found=True ; store_Kp (line) ; break

        if date_1_found: print(date_1, ' data found')
        
    # if date_2_found: print(date_2, ' data found')
    # else: 
    #     kp_date = 'kp'+date_2[2:6]+'.tab'
    #     url = site+kp_date
    #     lines = urllib.urlopen(url).read().decode('utf-8').splitlines()
    #     k=len(lines)
    #     for i in range(k):
    #         line=' '.join(lines[i].split(' ')).split()
    #         if line[0]==date_2[2:8]: 
    #             date_2_found=True ; store_Kp (line) ; break

    #     if date_2_found: print(date_2, ' data found')
        
#    if np.max(list(kp[str(np.max(list(map(int, list(kp.keys())))))].values())) > 5: bigKp = True
    print('preparing Kp input file ...')
    
#    print(kp)

    # for k1,v1 in sorted(kp.items()):
    #     Kp_day=[]; Kpy=[]; time_Kp=[]
    #     m=0
    #     for k2,v2 in v1.items():
    #         time_Kp.append('20'+k1+str(m*3).zfill(2)+'00');Kpy.append(v2)
    #         m=m+1
                                
    #     if (len(Kpy))!=8:
    #         for i in range(8-len(Kpy)):
    #             time_Kp.append('20'+k1+str(m*3).zfill(2)+'00');Kpy.append(0.0)
    #             m=m+1
    #     Kp_day.append(time_Kp); Kp_day.append(Kpy)
        
    # print(Kp_day)
    
    red_date = rem0(str(day))+'-'+rem0(str(month))+'-'+str(year)
    
    with open(user_path_in+"kp" + red_date + ".dat", 'w') as f:  # write a text file to be read by fortran codes
        space='       '
        t=0
        for k1,v1 in sorted(kp.items()):
            t=t+1
            for v in v1.items():
                # for j in range(len(list(v.keys()))):
                #     f.write(space+str(t)+space+str((list((v.keys()))[j]-1)*3)+space+str("%10.8f"%(list((v.values()))[j]))+space+'1.00000'+space+'2.00000\n')
                f.write(space+str(t)+space+str((v[0]-1)*3)+space+str("%10.8f"%(v[1]))+space+'1.00000'+space+'2.00000\n')
        f.write('   -9999      0      0.00000  1.00000  2.00000\n')
           
    return 

#################################################################################################################
def store_Kp (line):
    kp[line[0]]={}
#    print(line)
    for j in range(1,9):
        try:
            if line[j][1]=='o': kp[line[0]][j]=float(line[j][0])
            elif line[j][1]=='-': kp[line[0]][j]=float(line[j][0])-0.33333333
            elif line[j][1]=='+': kp[line[0]][j]=float(line[j][0])+0.33333333
        except: print(j,' not yet available')
    return

#################################################################################################################
def Kp_noaa(year,month,day,user_path_in):
    
    date=str(year) + str(month).zfill(2) + str(day).zfill(2)
    df_read=pd.read_csv(user_path_in + 'kp_forecast.csv',sep=';',header=None,index_col=0,names=['kp'])

    df_read.index = pd.to_datetime(df_read.index)
    df=pd.DataFrame(columns=[0,3,6,9,12,15,18,21])
    for i in range (0,len(df_read),8):
        d=str(df_read.iloc[i].name)[:10]
        df.loc[d]=list(df_read.loc[d].kp)
        
    date_1 = str((datetime.datetime.strptime(date,'%Y%m%d').date() - relativedelta.relativedelta(days=1)).strftime('%Y-%m-%d'))

    print('preparing Kp input file ...')

    red_date = rem0(str(day))+'-'+rem0(str(month))+'-'+str(year)
    date0 = str(year).zfill(2)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)
    
    with open(user_path_in+"kp" + red_date + ".dat", 'w') as f:  # write a text file to be read by fortran codes
        space='       '
        t=0
        for i in [date_1,date0]:
            t=t+1
            for v in range(8):
                f.write(space+str(t)+space+str(list(df.columns)[v])+space+str("%10.0f"%(list(df.loc[i])[v]))+space+'1.00000'+space+'2.00000\n')
        f.write('   -9999      0      0.00000  1.00000  2.00000\n')
           
    return 

#################################################################################################################
def rem0 (a):
    if a[0]=='0': a=str(a[1])
    return a

#################################################################################################################