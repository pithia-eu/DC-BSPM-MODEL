#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Dec 1 2019
@author: edithb

20230522 
- change MLT1 and MLT2 calculations in mer_coord_grid function to get any meridian plane

"""
############ IMPORT MODULES ############
import os, datetime, subprocess, csv, os.path, glob, copy
import urllib.request as urllib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from dateutil import relativedelta
from astropy.time.formats import erfa, TimeFromEpoch
import time
from matplotlib.patches import Wedge
import argparse
from datetime import date
from ctypes import cdll
import pandas as pd
from get_pp import create_plasmapause
from electron_density import get_plasma_ionosphere_density_electrons_19jan2016
#from artificial_points import fill_pp_holes_int,create_points
from scipy.interpolate import griddata
from artificial_points_allMLT import fill_pp_holes_allMLT
from holes import fill_pp_holes_int
#import png2gif

# Ininialize parameters
parser = argparse.ArgumentParser(description='BPIM')
parser.add_argument("--year", type=int, default=datetime.datetime.now().year)
parser.add_argument("--month", type=int, default=datetime.datetime.now().month)
parser.add_argument("--day", type=int, default=datetime.datetime.now().day)
parser.add_argument("--hour", type=int, default=-1)
parser.add_argument("--min", type=int, default=0)
#parser.add_argument("--year", type=int, required=True)
#parser.add_argument("--month", type=int, required=True)
#parser.add_argument("--day", type=int, required=True)
parser.add_argument("--var", default="density", type=str)
parser.add_argument("--number_of_xpoints", default=100, type=int) # grid resolution
parser.add_argument("--add_trough", default=1, type=int) 
parser.add_argument("--verbose", default=1, type=int) # EDITH : it is not used by now !!!!!
parser.add_argument("--rmin", default=1, type=int, choices=[1,2]) # # rmin = 2 DO NOT plot ionosphere ; rmin = 1 # plot ionosphere
parser.add_argument("--rmax", default=12, type=int) # pour plasma trough; change to 7 for the plasmasphere !!!
parser.add_argument("--nspot", default=60, type=int)
parser.add_argument("--MLT", default=12, type=float)
parser.add_argument("--iri_max_height", default=700, type=int) # altitude in km ; maximum 1000km; Precondition: iri_max_height < 2 Earth Radii
parser.add_argument("--ions_max_height", default=1950, type=int) # altitude in km ; maximum 2000 km (2000 km creates zeros at arbitrair locations by IRI)
parser.add_argument("--iri_temp_max_height", default=700, type=int) # maximum 2500 km
parser.add_argument("--r_maxi", default=6, type=int) # define the scale of the rectangle view
parser.add_argument("--ini", default=0, type=int)
parser.add_argument("--nt", default=48, type=int) # set the number of frames: nt
parser.add_argument("--step", default=2, type=int) # default every hour; set to 1 if every 30 min wanted
parser.add_argument("--interp", default=1, type=int) # set to 0 if not wanted
parser.add_argument("--troughParam", default='VAP', type=str) # other: 'Sheeley', 'VAP'
parser.add_argument("--executionid", default=datetime.datetime.now().strftime('%Y-%m-%d'), type=str) # other: 'Sheeley', 'VAP'

args = parser.parse_args()
print(args)

year = args.year; month = args.month; day = args.day
hour = args.hour; minu = args.min
ini = args.ini ; step = args.step
number_of_xpoints = args.number_of_xpoints 
number_of_ypoints = number_of_xpoints 
number_of_zpoints = int(number_of_xpoints/2)
Var = args.var
add_trough = args.add_trough
verbose = args.verbose
rmin = args.rmin ; rmax = args.rmax
nspot = args.nspot
MLT = args.MLT
iri_max_height = args.iri_max_height
ions_max_height = args.ions_max_height
iri_temp_max_height = args.iri_temp_max_height
r_maxi = args.r_maxi
nt = args.nt
interp = args.interp
troughParam = args.troughParam
executionid = args.executionid

start_time = time.time()
directory = os.path.dirname(__file__)
pause_path = os.path.join(directory, 'Plasmapause') # create path to directory "Plasmapause"

user = os.environ.get('USER')
user_path_in = os.path.join(directory, '../in/bspm/'+user+'/')
user_path_out = os.path.join(directory, '../out/bspm/'+user+'/'+executionid+'/')
iri_path = os.path.join(directory, 'Libs/iri2016/') 
#unilib_path = os.path.join(directory, 'Libs/') 
#print(directory)
#print(iri_path)

irilib = cdll[iri_path + 'irilib64.so'] # to communicate by ctypes
#unilib = cdll[unilib_path + 'unilib.so'] # to communicate by ctypes

pause_kp_filename = ""
xmax = r_maxi ; ymax = r_maxi ; zmax = r_maxi / 2.
vari={'density':{'units':' [1/cm\u00b3]','vmin':1e0, 'vmax':1e4},'temperature':{'units':' [K]','vmin':1e2, 'vmax':1e5}}
title='Electron '+ Var + vari[Var]['units']
vmin=vari[Var]['vmin'] ; vmax=vari[Var]['vmax']    
lambda_equatorial = 0.0 # in degrees
cmap = copy.copy(mpl.cm.get_cmap("turbo"))

#################################################################################################################
def update_iri_data(ig_rz,apf107): 
    file1=iri_path+ig_rz
#    mon=str(time.localtime(os.path.getmtime(file1)).tm_mon).zfill(2)
#    year=str(time.localtime(os.path.getmtime(file1)).tm_year)
#    int(year+mon)
    tdy = date.today()
    date_tdy = int(str(tdy.year)+str(tdy.month).zfill(2))
    
    with open(file1) as f: lines=f.read().splitlines()
    f.close()
    line=' '.join(lines[2].split(',')).split()
    date_max = int(str(line[3]+line[2]))
    
    if date_tdy >= date_max: 
        url='https://chain-new.chain-project.net/echaim_downloads/'
        urllib.urlretrieve(url+ig_rz, iri_path+ig_rz)
        urllib.urlretrieve(url+apf107, iri_path+apf107)
        print('ig_rz.dat and apf107.dat have been updated')
    else:
        pass
    
    return

#################################################################################################################
def dual_half_circle(center, radius, angle=0, ax=None, colors=('w','k'), **kwargs):
    """
    Add two half circles to the axes *ax* (or the current axes) with the 
    specified facecolors *colors* rotated at *angle* (in degrees).
    """
    if ax is None:
        ax = plt.gca()
    theta1, theta2 = angle, angle + 180
    w1 = Wedge(center, radius, theta1, theta2, fc=colors[0], linestyle='-', ec= 'k', **kwargs)
    w2 = Wedge(center, radius, theta2, theta1, fc=colors[1], linestyle='-', **kwargs)
    for wedge in [w1, w2]:
        ax.add_artist(wedge)
    return [w1, w2]

#################################################################################################################
class CDFEpoch(TimeFromEpoch):
    name = 'cdf_epoch'
    unit = 1.0 / (erfa.DAYSEC * 1000) # Milliseconds
    epoch_val = '0000-01-01 00:00:00'
    epoch_val2 = None
    epoch_scale = 'utc'
    epoch_format = 'iso'
    
#################################################################################################################
def rem0 (a):
    if a[0]=='0': a=str(a[1])
    return a

##################################################################################################################
def get_nspot(year,month,day,nspot): # revised
    try:
        date=str(year) + '  ' + str(month) 
        datez=str(year) + ' ' + str(month).zfill(2) 
        import urllib.request
        url='http://www.sidc.be/silso/DATA/SN_ms_tot_V2.0.txt' # EDITH: check if this is the right one!!!!!!
        data = urllib.request.urlopen(url)
        for line in data:
            try:
                if str(date) in str(line) or str(datez) in str(line): 
                    nspot_test = ' '.join(str(line).split(' ')).split()[3]
                    if float(nspot_test) > 0.0 : 
                        nspot=nspot_test 
                        print('From: ', url, ', nspots for ', date, ' = ', nspot)
                    else : print('Assume nspots for ', date, ' = ', nspot)
            except: pass
    except: pass
    return nspot

#####################################################################################
def eq_coord_grid ():
# calculate plasmasphere equatorial plane
    x_array_solar_magnetic_eq,y_array_solar_magnetic_eq = np.mgrid[-xmax:xmax:(number_of_xpoints*1j),-ymax:ymax:(number_of_ypoints*1j)]
    
    radius_in_R_E_eq = np.zeros([number_of_xpoints, number_of_ypoints])
    longitude = np.zeros([number_of_xpoints, number_of_ypoints])
    magnetic_local_time_eq = np.zeros([number_of_xpoints, number_of_ypoints])
                   
    for i in range(number_of_xpoints):
        for j in range(number_of_ypoints):
            x_sm_coord = x_array_solar_magnetic_eq[i][j]
            y_sm_coord = y_array_solar_magnetic_eq[i][j]
            radius_in_R_E_eq[i][j] = np.sqrt((x_sm_coord)**2 + (y_sm_coord)**2)
            longitude[i][j] = np.rad2deg(np.arctan2(y_sm_coord, x_sm_coord))  
    
            if longitude[i][j] < 0  : longitude[i][j] = longitude[i][j] + 360.          
            if longitude[i][j] < 180: magnetic_local_time_eq[i][j] = 12. + longitude[i][j] * (24./360.)
            else:                     magnetic_local_time_eq[i][j] = (longitude[i][j]-180.) * (24./360.)
        
    return x_array_solar_magnetic_eq , y_array_solar_magnetic_eq , radius_in_R_E_eq , magnetic_local_time_eq   

#####################################################################################
def mer_coord_grid (MLT): # add MLT if one !=0 in input
    global MLT1,MLT2

# calculate plasmapause meridian plane
# sun is to the left: MLT=12 is to the left, MLT=0 is to the right

    if (MLT>=6.0 and MLT<18.0): 
        MLT1 = MLT
        MLT2 = (MLT+12) % 24.
    else:
        if MLT>=18: MLT1 = MLT-12
        else: MLT1=MLT+12
        MLT2 = MLT
#    print(MLT1,MLT2)
    # latitude: lambda is from -90 to 90 (in degrees)
    # lambda = ((180.*findgen(angpoints/2)) / (angpoints/2)) - 90.
    #  reverse_lambda = -lambda
    
    x_array_solar_magnetic_mer,z_array_solar_magnetic = np.mgrid[-xmax:xmax:(number_of_xpoints*1j),-zmax:zmax:(number_of_zpoints*1j)]
    
    radius_in_R_E_mer = np.zeros([number_of_xpoints, number_of_zpoints])
    latitude_in_degrees_mer = np.zeros([number_of_xpoints, number_of_zpoints])
    magnetic_local_time_mer = np.zeros([number_of_xpoints, number_of_zpoints])

    for i in range(number_of_xpoints):
        for j in range(number_of_zpoints):
            x_sm_coord = x_array_solar_magnetic_mer[i][j]
            z_sm_coord = z_array_solar_magnetic[i][j]
            radius_in_R_E_mer[i][j] = np.sqrt((x_sm_coord)**2 + (z_sm_coord)**2)
        
            if x_sm_coord != 0: latitude_in_degrees_mer[i][j] = np.rad2deg(np.arctan2(z_sm_coord, np.abs(x_sm_coord)))
            else:               latitude_in_degrees_mer[i][j] = 0
            if x_sm_coord >= 0: magnetic_local_time_mer[i][j] = MLT1 # 12; dayside; longitude=0
            else:               magnetic_local_time_mer[i][j] = MLT2 # 0; nightside; longitude=180

    return x_array_solar_magnetic_mer , z_array_solar_magnetic , radius_in_R_E_mer , magnetic_local_time_mer , latitude_in_degrees_mer# end meridian coordinate grid

#####################################################################################
def read_trajdat (pause_traject_filename):
    # Read the file with the plasmapause of the equatorial plane
#    print(pause_traject_filename)
    trajdat = open(os.path.join(user_path_in, pause_traject_filename), "r")
    x = []                                                                      # Initalise List
    reader = csv.reader(trajdat,delimiter=' ',skipinitialspace=True)            # Initialise a reader able to read the different columns separated with a space
    for row in reader: x.append(row)                                            
    trajdat.close()                                                             # Close the file
    trajdat = np.array(x,dtype=float)
    nl = len(trajdat)
    tind = np.where((trajdat[0:nl-1,0] < -1.0E21) & (trajdat[1:nl,0] < -1.0E21)) 
#    print(nl,len(list(tind)),tind)
    return trajdat, tind

#####################################################################################
def calc_plasmapause ():
        #**************************** build the traject.dat file *************************************#
    pause_kp_filename = "kp" + str(day) + "-" + str(month) + "-" + str(year) + ".dat"
    pause_input_filename = 'input_' + str_date + '.dat'
    pause_input_filename_exists = os.path.isfile(os.path.join(user_path_in, pause_input_filename))
    
    if not pause_input_filename_exists:
        command_ed = [os.path.join(pause_path,'ed_input'), os.path.join(user_path_in,pause_input_filename),\
                      os.path.join(pause_path,'ed_input.nml')]
        subprocess.call(command_ed)  # run fortran program (first of the list) with the rest of the list of command_ed as input files
    
#    time.sleep(5)
    pause_fkp_filename = 'fkp_' + str_date + '.dat'     # smoothing kp's !!!!!!!!!
    pause_fkp_filename_exists = os.path.isfile(os.path.join(user_path_in,pause_fkp_filename))
    if not pause_fkp_filename_exists:
        command_main = [os.path.join(pause_path,'fkp_main'), os.path.join(user_path_in, pause_kp_filename),\
                        os.path.join(user_path_in, pause_fkp_filename), os.path.join(pause_path,'fkp_main.nml')]
        subprocess.call(command_main) # run fortran program (first of the list) with the rest of the list of command_main as input files
    
#    time.sleep(5)    
    pause_traject_filename= 'traject_' + str_date + '.dat'
    pause_output_filename= 'output_' + str_date + '.dat'
    pause_reprise_dat_filename= 'reprise_' + str_date + '.dat'
    pause_reprise_sav_filename= 'reprise_' + str_date + '.sav'
    
    pause_traject_filename_exists = os.path.isfile(os.path.join(user_path_in, pause_traject_filename))
    if not pause_traject_filename_exists:
        command_pls = [os.path.join(pause_path, 'pls_kz'), os.path.join(user_path_in, pause_input_filename),\
                       os.path.join(user_path_in, pause_kp_filename), os.path.join(user_path_in, pause_fkp_filename),\
                       os.path.join(user_path_in, pause_traject_filename), os.path.join( user_path_in, pause_output_filename), \
                       os.path.join(user_path_in, pause_reprise_dat_filename), os.path.join(user_path_in, pause_reprise_sav_filename),\
                       os.path.join(pause_path, 'pls_kz.nml')]
        subprocess.call(command_pls) # run fortran program (first of the list) with the rest of the list of command_main as input files
# the traject.dat file is made 
#    time.sleep(5)    
    return pause_kp_filename,pause_traject_filename

#####################################################################################
def animate (radius_in_R_E_eq,magnetic_local_time_eq,radius_in_R_E_mer,magnetic_local_time_mer,latitude_in_degrees_mer,var,xkp,ykp,plasmapause,dfk,MLT1):
    global cs_eq, p_equatorial, cs_mer, p_meridian, tx, pl_p
    
    new_eq = np.intersect1d(np.where(radius_in_R_E_eq >= rmin),np.where(radius_in_R_E_eq <= rmax))   
    new_mer = np.intersect1d(np.where(radius_in_R_E_mer >= rmin),np.where(radius_in_R_E_mer <= rmax))

    for k_slide in range(ini,nt,step):
       
        bar,H=fig_init(xkp,MLT1)
        
        # lambdaE is in degrees    
        float_hour = k_slide / 2.0 # Universal Time
        integer_hour = int(k_slide / 2.0)
        if minu!=0 : integer_minutes = minu
        else : integer_minutes = int((float_hour - integer_hour) * 60.)
        time_plot = datetime.datetime(year,month,day,integer_hour,integer_minutes).strftime("%Y-%m-%d_%Hh%Mm")
        timestamp = datetime.datetime.now().strftime("_at_%Y-%m-%d_%Hh%Mm") 
        
        print(integer_hour,'H',integer_minutes,' UT')
        
        kp=float(dfk.loc[str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+' '+str(integer_hour).zfill(2)+':'+str(integer_minutes).zfill(2)+':00'])
    
    #   EQUATORIAL PLANE
        p_equatorial = np.zeros([number_of_ypoints,number_of_xpoints])
        p_equatorial = p_equatorial.reshape(number_of_xpoints*number_of_ypoints)
            
        for k in list(new_eq):
            p_equatorial[k] = get_plasma_ionosphere_density_electrons_19jan2016 (r=radius_in_R_E_eq[k], lambdaE=lambda_equatorial, \
                                                                                 MLT=magnetic_local_time_eq[k], year=year, month=month, day=day, \
                                                                                 k_slide=k_slide, nspot=nspot, \
                                                                                 iri_max_height=iri_max_height, add_trough=add_trough, a_plasmapause=plasmapause[k_slide][1], \
#                                                                                 a_clean=plasmapause[k_slide][0], unilib=unilib, verbose=verbose, kp=kp, troughParam=troughParam )
                                                                                 a_clean=plasmapause[k_slide][0], troughParam=troughParam )
                  
        
        filename_results_eq = 'dens_eq_' + time_plot + timestamp + '.csv'
        df=pd.DataFrame([radius_in_R_E_eq,magnetic_local_time_eq,p_equatorial]).transpose()
        df.columns=['r[Re]','MLT[hours]','density_eq[1/cm\u00b3]']
        df.to_csv(user_path_out+filename_results_eq, sep=',',float_format='%.3f',index=False)

        p_equatorial = p_equatorial.reshape(number_of_xpoints,number_of_ypoints)
    
    #    # latitude=90 is magnetic north pole
    #    # latitude=-90 is magnetic south pole
    #    # the points are given in the order: Z pole -> MLT 12 -> N pole -> MLT 0 -> Z pole
    
    #   MERIDIAN PLANE    
        p_meridian = np.zeros([number_of_xpoints,number_of_zpoints])
        p_meridian = p_meridian.reshape(number_of_xpoints*number_of_zpoints)
         
        for k in list(new_mer):
            p_meridian[k] = get_plasma_ionosphere_density_electrons_19jan2016 (r=radius_in_R_E_mer[k], lambdaE=latitude_in_degrees_mer[k], \
                                                                               MLT=magnetic_local_time_mer[k], year=year, month=month, day=day, \
                                                                               k_slide=k_slide, nspot=nspot, \
                                                                               iri_max_height=iri_max_height, add_trough=add_trough, a_plasmapause=plasmapause[k_slide][1], \
#                                                                               a_clean=plasmapause[k_slide][0], unilib=unilib, verbose=verbose, kp=kp, troughParam=troughParam )
                                                                               a_clean=plasmapause[k_slide][0], troughParam=troughParam )
            
        filename_results_mer='dens_mer_' + time_plot + timestamp + '.csv'
        df=pd.DataFrame([radius_in_R_E_mer,magnetic_local_time_mer,latitude_in_degrees_mer,p_meridian]).transpose()
        df.columns=['r[Re]','MLT[hours]','Lat[Deg]','density_mer[1/cm\u00b3]']
        df.to_csv(user_path_out+filename_results_mer, sep=',',float_format='%.3f',index=False)

        p_meridian = p_meridian.reshape(number_of_xpoints,number_of_zpoints)
#        p_meridian = np.ma.masked_where((-0.05 < x_array_solar_magnetic_mer.any() < 0.05) & \
#                                        (z_array_solar_magnetic.any() < -1) & \
#                                        (z_array_solar_magnetic > 1).any(), p_meridian)  # small correction

        x,y=np.mgrid[-xmax:xmax:(4000*1j),-ymax:ymax:(4000*1j)]        
        grid_z2=griddata((-xxx.flatten(),yyy.flatten()),p_equatorial.flatten(),(x,y),method='linear')
        grid_z2 = np.ma.masked_where(grid_z2 == 0.0,grid_z2)        
        cs_eq = ax1.imshow(grid_z2.T,extent=[-xmax, xmax, -ymax, ymax],norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax),cmap=cmap)
        
        x,y=np.mgrid[-xmax:xmax:(4000*1j),-zmax:zmax:(2000*1j)]        
        grid_z2=griddata((-xxxx.flatten(),zzzz.flatten()),p_meridian.flatten(),(x,y),method='linear')
        grid_z2 = np.ma.masked_where(grid_z2 == 0.0,grid_z2)        
        cs_mer = ax2.imshow(grid_z2.T,extent=[-xmax, xmax, -zmax, zmax],norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax),cmap=cmap)
        
#        cs_eq = ax1.scatter(-yyy,-xxx,c=p_equatorial.T,s=msize,marker='s',edgecolors='none',norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax),cmap='jet')
#        cs_mer = ax2.scatter(-xxxx,zzzz,c=p_meridian,s=msize,marker='s',edgecolors='none',norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax),cmap='jet')

#        cs_eq = ax1.pcolormesh(-yyy,-xxx,p_equatorial.T,norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax),shading='gouraud',cmap='jet')
#        cs_mer = ax2.pcolormesh(-xxxx,zzzz,p_meridian,norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax),shading='gouraud',cmap='jet')
        
#        cs_eq = ax1.contourf(-yyy,-xxx,p_equatorial.T,levels=np.linspace(vmin,vmax,200),locator=ticker.LogLocator(),cmap=mpl_cm.get_cmap('jet'))
#        cs_mer = ax2.contourf(-xxxx,zzzz,p_meridian,levels=np.linspace(vmin,vmax,200),locator=ticker.LogLocator(),cmap=mpl_cm.get_cmap('jet'))

        cb = fig.colorbar(cs_eq,orientation='horizontal',cax=cb_ax)
        cb.ax.tick_params(labelsize=12)
        
        xp,yp = plasmapause[k_slide][2][:,0],plasmapause[k_slide][2][:,1]
        pl_p = ax1.scatter(xp,yp,facecolors='none', edgecolors='k')        # plot pp

        # date_data_1=str(year)+'/'+str(month).zfill(2)+'/'+str(day-1).zfill(2)+'/'
        # PP_v=pd.read_csv(user_path_out+date_data_1+'PP_2015-03-25_12h00m.csv', sep=',')
        # PP_v.drop(PP_v.loc[(PP_v['01']>-5)&(PP_v['01']<-4)&(PP_v['02']>-0.5)&(PP_v['02']<2.5)].index, inplace=True)
        # PP_v.drop(PP_v.loc[(PP_v['01']>1.5)&(PP_v['01']<3.5)&(PP_v['02']>-4)&(PP_v['02']<-2.5)].index, inplace=True)
        # PP_v.drop(PP_v.loc[(PP_v['01']>-4)&(PP_v['01']<-2)&(PP_v['02']>-5)&(PP_v['02']<-2)].index, inplace=True)
        # pl_p_v = ax1.scatter(PP_v['01'],PP_v['02'],facecolors='none', edgecolors='grey') 

        # dfPP1 = pd.DataFrame(plasmapause[k_slide][2])
        # dfPP1.drop(dfPP1.loc[dfPP1[1]<=0].index, inplace=True)

        # dfPP2 = pd.DataFrame(plasmapause[k_slide][2])
        # dfPP2.drop(dfPP2.loc[dfPP2[1]>=0].index, inplace=True)

        # # print(dfPP1)
        # # print(dfPP2)
        
        # x_new=np.linspace(-6,6,100)        
        # import numpy.polynomial.polynomial as poly
        # coefs = poly.polyfit(dfPP1[0], dfPP1[1], 5)
        # ffit = poly.Polynomial(coefs)    
        # ax1.plot(x_new, ffit(x_new),color='k',linewidth=3)

        
# Temporary to check PP        
        rp = np.zeros([xp.size]); Lon = np.zeros([xp.size]); MLTp = np.zeros([xp.size])
        rp = np.sqrt((xp)**2 + (yp)**2)
        Lon = np.rad2deg(np.arctan2(-yp, -xp))  
        for p in range(xp.size):
            if Lon[p] < 0  : Lon[p] = Lon[p] + 360.          
            if Lon[p] < 180: MLTp[p] = 12. + Lon[p] * (24./360.)
            else:            MLTp[p] = (Lon[p]-180.) * (24./360.)
        
        filename_PP='PP_' + time_plot + timestamp + '.csv'
        dfp=pd.DataFrame([xp,yp,rp,MLTp]).transpose()
        dfp.columns=['XPgsm','YPgsm','Rp','MLTp']
        for column in dfp.columns: dfp = dfp[dfp[column]!=-2000000000000000000000.000]
        dfp=dfp.sort_values(by=['MLTp'])        
#        dfp.to_csv(user_path_out+filename_PP, sep=',',float_format='%.3f',index=False)
# Temporary to check PP                      

        ax4.text(0.5, 0.5,time_plot,fontsize=18,fontweight="bold")
        ax5.text(0.3, 0.5, title ,fontsize=18,fontweight="bold")
        ax6.set_title('Axes units in Re',size=14)        
        ax7.set_title('Belgian SWIFF Plasmasphere Model v.2021',size=16, fontstyle='italic')
        
        
# Print Kp bars         
#         ykpcut=[0]*len(xkp)
#         a=list(np.arange(0,25,3))
# #        for i in range(len(xkp[2*H:])):
#         for i in range(len(xkp[H:])):
#             if (integer_hour >= a[i]) & (integer_hour < a[i+1]):
# #                ykpcut[:2*H+i+1] = [ykp[j] for j in range(2*H+i+1)]     
#                 ykpcut[:H+i+1] = [ykp[j] for j in range(H+i+1)]     
#                 break
            
#         print(xkp,ykp)    
#        ax3.bar(xkp,ykpcut,bar,color='k',align='edge')
        ax3.bar(xkp,ykp,bar,color='k',align='edge')
        ax3.set_xlim(0,72,3)
# End Print Kp bars
            
#        timestamp = datetime.datetime.now().strftime("_at_%d-%b-%Y_%Hh%Mm")
        plt.savefig(user_path_out+'dens_'+time_plot+timestamp+'.png',dpi=75,bbox_inches='tight')
        plt.close(fig)
        
#    png2gif.all(user_path_out)
           
#    print(time.time() - start_time, ' seconds')
    return cs_eq, cs_mer, cb

#####################################################################################
def fig_init (xkp,MLT1):
    global ticksz, ticksx, msize, ax1, ax2, ax3, ax4, ax5, ax6, ax7, cb_ax, fig
    fig = plt.figure(figsize=(16,10),constrained_layout=False)
    gs = fig.add_gridspec(15,8)

    ticksx=np.arange(-xmax,xmax+1)
    if isinstance(zmax,int): ticksz=np.arange(-zmax,zmax+1)
    else: ticksz=np.arange(-zmax,zmax,dtype=int)

    circles = [1.0,2.0,4.0]
    msize=200
    
    ax1 = fig.add_subplot(gs[3:,:-4])
    ax1.set_title('Equatorial plane',fontsize=18,fontweight='bold')
    for c in circles:
        ax1.add_artist(plt.Circle((0,0), c, color='k', linewidth=0.5, linestyle=':', fill=False))
    dual_half_circle((0.0), radius=1.0, angle=90, ax=ax1)
    plt.xlim(-xmax,xmax); plt.ylim(-ymax,ymax)
    plt.xticks(ticksx);plt.yticks(ticksx)
    plt.grid(linestyle='--')
    ax1.tick_params(axis='both', which='major', labelsize=12) 
    
    
    ax2 = fig.add_subplot(gs[6:12,-4:])
    if MLT1==12: ax2.set_title('Meridian plane',fontsize=18,fontweight='bold')
    else: ax2.set_title('Meridian plane (MLT='+str(MLT1)+')',fontsize=18,fontweight='bold')
    for c in circles:
        ax2.add_artist(plt.Circle((0,0), c, color='k', linewidth=0.5, linestyle=':', fill=False))
    dual_half_circle((0.0), radius=1.0, angle=90, ax=ax2)
    plt.xlim(-xmax,xmax); plt.ylim(-zmax,zmax)
    plt.xticks(ticksx);plt.yticks(ticksz)
    
    plt.grid(linestyle='--')
    ax2.tick_params(axis='both', which='major', labelsize=12) 

    ax3 = fig.add_subplot(gs[:2,3:7])
    today = datetime.datetime(year,month,day).strftime("%d-%b-%Y")
    yesday = (datetime.datetime(year,month,day) - relativedelta.relativedelta(days=1)).strftime("%d-%b-%Y")
    tomday = (datetime.datetime(year,month,day) + relativedelta.relativedelta(days=1)).strftime("%d-%b-%Y")
#    byesday = (datetime.datetime(year,month,day) - relativedelta.relativedelta(days=2)).strftime("%d-%b-%Y")
    ax3.set_title(yesday+'                 '+today+'                 '+tomday)
#    ax3.set_title(byesday+'                 '+yesday+'                 '+today)
    ax3.set_ylabel('Kp',fontsize=10)
#    
    if xkp[1]-xkp[0]==1: plt.xticks(np.arange(0,72,3),3*list(np.arange(0,24,3)),ha='left'); bar=0.8 ; H=24
    else: plt.xticks(np.arange(0,70,3),3*list(np.arange(0,24,3)),ha='left'); bar=2.8 ; H=8
    plt.ylim(0,8); plt.yticks(np.arange(0,9,2)); ax3.set_yticklabels([0,2,4,6,8])    
    plt.grid(linestyle='--')
    
    ax4 = fig.add_subplot(gs[:2,:2]); ax4.axis('off')
    ax5 = fig.add_subplot(gs[13:,-4:]); ax5.axis('off')
    ax6 = fig.add_subplot(gs[5:6,-7:]); ax6.axis('off') 
    ax7 = fig.add_subplot(gs[4:5,-3:]); ax7.axis('off')

    fig.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9,wspace=0.2, hspace=0)
#    fig.tight_layout()
#    plt.subplots_adjust(bottom=0.1)
#    cb_ax = fig.add_axes([0.17, 0.035, 0.7, 0.04])
    cb_ax = fig.add_axes([0.508, 0.11, 0.39, 0.015])
    
#    plt.suptitle(title, fontsize=16)    

    return bar,H

#############################    PROGRAM BEGINS HERE    ###############################
print('Current dir: ',os.getcwd())
dire=os.getcwd()
if not os.path.exists(user_path_in): os.makedirs(user_path_in)
if not os.path.exists(user_path_out): os.makedirs(user_path_out)

str_date = "-".join([str(year),str(month), str(day)]) 
str_date_new = "-".join([str(year),str(month).zfill(2), str(day).zfill(2)])
#red_date = rem0(str(day))+'-'+rem0(str(month))+'-'+str(year)

date_data=str(year)+'/'+str(month).zfill(2)+'/'+str(day).zfill(2)+'/'
if not os.path.exists(user_path_out+date_data): os.makedirs(user_path_out+date_data)

update_iri_data('ig_rz.dat','apf107.dat') # This is checked at every run !!!!!

from get_kp import Kp_new; Kp_new(year,month,day,user_path_in); #user_path_kp = user_path_in+'kp_POTSDAM/'
#from get_kp import Kp_noaa; Kp_noaa(year,month,day,user_path_in); #user_path_kp = os.path.join(directory, 'Results/kp_NOAA/')

nspot=get_nspot(year,month,day,nspot) # Read nspot  # EDITH: CHECK!!!!!!!!!
pause_kp_filename,pause_traject_filename = calc_plasmapause() 
trajdat, tind = read_trajdat (pause_traject_filename)
cols=['01','02','03','04','05','06','07','08','09','10','11','12']
# count the number of TRAJECT*.DAT files --> this will be the number of pictures
# nt = len(tind[0][:]) - 1 # in an animation nt is 48 (48 halfhours in a day)


##################################################################################################################
###### Calculate PP
plasmapause={}

#datee = str(year)+str(month).zfill(2) + str(day).zfill(2)+'/'
# if os.path.exists(user_path_kp+datee):   # PP already stored for this date
#     for k_slide in range(0,nt,2):
#         float_hour = k_slide / 2.0 # Universal Time
#         integer_hour = int(k_slide / 2.0)
#         integer_minutes = int((float_hour - integer_hour) * 60.)
#         time_plot = datetime.datetime(year,month,day,integer_hour,integer_minutes).strftime("%Y-%m-%d_%Hh%Mm")
#         df_aclean=pd.read_csv(user_path_kp+datee+'ACLEAN_'+time_plot+'.csv')
#         df_app=pd.read_csv(user_path_kp+datee+'APP_'+time_plot+'.csv')
#         df_eppc=pd.read_csv(user_path_kp+datee+'EPPC_'+time_plot+'.csv')
#         plasmapause[k_slide]=[np.array(df_aclean), np.array(df_app), np.array(df_eppc)]
        
#else:    # PP not available, so create PP
#os.makedirs(user_path_kp+str(year)+str(month).zfill(2) + str(day).zfill(2))
df=pd.DataFrame()
fill=False

if hour>=0 : ini = hour * 2 ; nt = ini + 1 ; step = 1

datee = datetime.datetime(year,month,day).strftime("%Y%m%d")
datenow = str(datetime.datetime.now().year)+str(datetime.datetime.now().month).zfill(2)+str(datetime.datetime.now().day).zfill(2)

if datee==datenow and hour==-1 : nt= datetime.datetime.now().hour *2 # to finish the loop of hours now

print(ini,nt,step)

for k_slide in range(ini, nt, step):
    if k_slide == 0: equatorial_plasmapause_cartesian = trajdat[0:tind[0][0]+2,:] 
    else:            equatorial_plasmapause_cartesian = trajdat[tind[0][k_slide-1]+2:tind[0][k_slide]+2,:]     

    # equatorial_plasmapause_polar(0) is theta
    # equatorial_plasmapause_polar(1) is r
    # equatorial_plasmapause_cartesian(0) is x
    # equatorial_plasmapause_cartesian(1) is y
    float_hour = k_slide / 2.0 # Universal Time
    integer_hour = int(k_slide / 2.0)
    if minu!=0 : integer_minutes = minu
    else : integer_minutes = int((float_hour - integer_hour) * 60.)
    time_plot = datetime.datetime(year,month,day,integer_hour,integer_minutes).strftime("%Y-%m-%d_%Hh%Mm")
    
    print(integer_hour,'H',integer_minutes,' UT')
    
    dfp=pd.DataFrame(equatorial_plasmapause_cartesian,columns=cols)    
#        dfp.to_csv(user_path_out+date_data+'PP_ori_'+time_plot+'.csv', sep=',',columns=cols,float_format='%.4f',index=False) 
    
########### pass here to fill holes by dynamic interpolation
    if interp == 1:        
        equatorial_plasmapause_cartesian=fill_pp_holes_int(equatorial_plasmapause_cartesian,k_slide,year,month,day,user_path_out,time_plot,date_data)
        equatorial_plasmapause_cartesian=fill_pp_holes_allMLT(equatorial_plasmapause_cartesian)
###########    
    ## Here the 12 columns of trajdat are passed !!!!! TODO: change !!!!!!
    equatorial_plasmapause_polar, a_clean, a_plasmapause = create_plasmapause (equatorial_plasmapause_cartesian=np.transpose(equatorial_plasmapause_cartesian)) 
    print('plasmapause created')    
    plasmapause[k_slide]=[a_clean, a_plasmapause, equatorial_plasmapause_cartesian]

#### Save PP 
    dfp=pd.DataFrame(equatorial_plasmapause_cartesian,columns=cols)    
    dfp.to_csv(user_path_out+date_data+'PP_'+time_plot+'.csv', sep=',',columns=cols,float_format='%.4f',index=False) 
# #### OTHER for routine  !!!!! Check if Ok
#     dfp=pd.DataFrame(equatorial_plasmapause_cartesian)    
#     dfp.to_csv(user_path_kp+datee+'EPPC_'+time_plot+'.csv', sep=',',float_format='%.4f',index=False) 
#     dfp=pd.DataFrame(a_clean)    
#     dfp.to_csv(user_path_kp+datee+'ACLEAN_'+time_plot+'.csv', sep=',',float_format='%.4f',index=False) 
#     dfp=pd.DataFrame(a_plasmapause)    
#     dfp.to_csv(user_path_kp+datee+'APP_'+time_plot+'.csv', sep=',',float_format='%.4f',index=False) 
#     ####    
###### END calculate PP and store
    
#remove unuseful files in /data/in/bpim
for f in set(glob.glob(user_path_in+"*.dat"))-set(glob.glob(user_path_in+"kp*")): os.remove(f)
for f in set(glob.glob(user_path_in+"*.sav")): os.remove(f)

# Kp
with open(os.path.join(user_path_in,pause_kp_filename)) as f: lines = f.read().splitlines()
f.close()
xkp=[];ykp=[];kpp=[]

date1=str(year)+str(month).zfill(2)+str(day).zfill(2)
date_plus1 = str(int((datetime.datetime.strptime(date1,'%Y%m%d').date() + relativedelta.relativedelta(days=1)).strftime('%Y%m%d')))
date_plus1=date_plus1[:4]+'-'+date_plus1[4:6]+'-'+date_plus1[6:]+' 00:00:00'

for i in range(len(lines)-1):
    line=lines[i].split(' ')
    line=' '.join(line).split()
    xkp.append(float(line[1])) ; ykp.append(float(line[2]))
    if line[0]=='2':  # current day
#    if line[0]=='3':  # current day
        kpp.append([str(str_date_new+' '+str(line[1].zfill(2)+':00:00')),line[2]])
kpp.append([str(date_plus1),line[2]])

#print(kpp)

dfk=pd.DataFrame(kpp)
dfk.index=pd.to_datetime(dfk[0])
dfk=dfk.drop(columns=[0])    
dfk=dfk.rename(columns={"1": "kpp"})
dfk=dfk.resample('30min').last().fillna(method='ffill')[:-1]

if xkp[1]-xkp[0]==3 : xkp=np.array([i for i in range(0, 70) if i % 3 == 0])
else : xkp=np.array([i for i in range(0, 72)])
xkp=xkp[:len(ykp)]
# end Kp

x_array_solar_magnetic_eq , y_array_solar_magnetic_eq , radius_in_R_E_eq , magnetic_local_time_eq = eq_coord_grid()   # calculate equatorial coordinate grid
x_array_solar_magnetic_mer , z_array_solar_magnetic , radius_in_R_E_mer , magnetic_local_time_mer , latitude_in_degrees_mer = mer_coord_grid(MLT)  # calculate meridian coordinate grid
bar,H=fig_init(xkp,MLT1)


radius_in_R_E_eq = radius_in_R_E_eq.reshape(number_of_xpoints*number_of_ypoints)            
magnetic_local_time_eq = magnetic_local_time_eq.reshape(number_of_xpoints*number_of_ypoints)
xxx , yyy = x_array_solar_magnetic_eq , y_array_solar_magnetic_eq
p_equatorial = np.zeros([number_of_ypoints,number_of_xpoints])
cs_eq = ax1.scatter(yyy,xxx,c=p_equatorial,s=msize,edgecolors='none',norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax),cmap=cmap)
cb = fig.colorbar(cs_eq,orientation='horizontal',cax=cb_ax)

latitude_in_degrees_mer = latitude_in_degrees_mer.reshape(number_of_xpoints*number_of_zpoints)
radius_in_R_E_mer = radius_in_R_E_mer.reshape(number_of_xpoints*number_of_zpoints)            
magnetic_local_time_mer = magnetic_local_time_mer.reshape(number_of_xpoints*number_of_zpoints)
xxxx,zzzz = x_array_solar_magnetic_mer,z_array_solar_magnetic   
p_meridian = np.zeros([number_of_xpoints,number_of_zpoints])  
cs_mer = ax2.scatter(xxxx,zzzz,c=p_meridian,s=msize,edgecolors='none',norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax),cmap=cmap)

animate(radius_in_R_E_eq,magnetic_local_time_eq,radius_in_R_E_mer,magnetic_local_time_mer,latitude_in_degrees_mer,Var,xkp,ykp,plasmapause,dfk,MLT1)

timestamp = datetime.datetime.now().strftime("_at_%d-%b-%Y_%Hh%Mm")
#ani.save(user_path_movies+'elect_'+Var[:4]+'_'+str(year)+'-'+str(month)+'-'+str(day)+'_#frames_'+str(nt)+timestamp+'.gif',writer='imagemagick')
                    
print(time.time() - start_time, ' seconds')    
exit()    




