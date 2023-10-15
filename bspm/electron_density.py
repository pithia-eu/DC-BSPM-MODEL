#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 15:32:56 2020

@author: edithb
"""

import numpy as np
from ctypes import cdll, byref, c_int, c_double, c_char, POINTER, Structure
from numpy.ctypeslib import ndpointer
import os, datetime
from datetime import timezone
from astropy.time import Time
from spacepy import coordinates as coord
from spacepy.time import Ticktock

directory = os.path.dirname(__file__)
#unilib_path = os.path.join(directory, 'Libs/') 
iri_path = os.path.join(directory, 'Libs/iri2016/') 
irilib = cdll[iri_path + 'irilib64.so'] # to communicate by ctypes
#unilib = cdll[unilib_path + 'unilib.so'] # to communicate by ctypes

#################################################################################################################
def get_plasma_ionosphere_density_electrons_19jan2016(r, lambdaE, MLT, year, month, day, \
                                                      k_slide, nspot, iri_max_height, add_trough, \
                                                      a_plasmapause, a_clean, troughParam ):
    """
    r is in Earth radii
    lambdaE is latitude from -90 to 90 degrees
    MLT is magnetic local time in hours
    1 Earth Radius = 6371 km
    2000 km = 0.313922 Earth Radii  
    """
#    print('edith:',r)
    if r < 1.0: density = 0.0 #; print('dens=0')
    else:
        iri_max_height_earth_radii = iri_max_height / 6371.

        if r <= (1. + iri_max_height_earth_radii):
            density = density_ionosphere_electrons_irbem_jul2021(r=r, lambdaE=lambdaE, MLT=MLT, k_slide=k_slide, nspot=nspot, \
                                                                    year=year, month=month, day=day)
        else: 
# radius_height_iono = 1.313922 Earth Radii
            radius_height_iono = 1. + iri_max_height_earth_radii
            radius_height_plasma = 2.0
# interpolate when 1. + iri_max_height_earth_radii <= r <= 2.0 along the magnetic field line of a dipole
            if r < 2.0: 
                l_shell = r / np.cos(np.deg2rad(lambdaE))**2

                if l_shell > 2.0:
                    magnetic_latitude_field_line_crossing_plasma = np.rad2deg(np.arccos(np.sqrt( 2.0/l_shell ))) 
                    # the output of arccos is between 0 and pi
                    # give magnetic_latitude_field_line_crossing_plasma the same sign as lambdaE
                    if lambdaE < 0: magnetic_latitude_field_line_crossing_plasma = -magnetic_latitude_field_line_crossing_plasma
                    if (magnetic_latitude_field_line_crossing_plasma < -90 or magnetic_latitude_field_line_crossing_plasma > 90): 
                        print('latitude not between -90 and 90')
                        
                    n_plasmasphere = density_plasmasphere(r=radius_height_plasma, lambdaE=magnetic_latitude_field_line_crossing_plasma, MLT=MLT, \
                                                          a_plasmapause=a_plasmapause, a_clean=a_clean, nspot=nspot, year=year, month=month, day=day, \
                                                          add_trough=add_trough, troughParam=troughParam )
                    
                    magnetic_latitude_field_line_crossing_iono = np.rad2deg(np.arccos(np.sqrt(radius_height_iono/l_shell )))
    
                    if lambdaE < 0: magnetic_latitude_field_line_crossing_iono = -magnetic_latitude_field_line_crossing_iono
                    if (magnetic_latitude_field_line_crossing_iono < -90 or magnetic_latitude_field_line_crossing_iono > 90): 
                        print('latitude not between -90 and 90')
                        
                    n_ionosphere = density_ionosphere_electrons_irbem_jul2021(r=radius_height_iono, lambdaE=magnetic_latitude_field_line_crossing_iono, \
                                                                                 MLT=MLT, k_slide=k_slide, nspot=nspot, \
                                                                                 year=year, month=month, day=day)
                    
                    # logarithmic interpolation
                    if n_plasmasphere != 0: 
                        density = 10**((((np.log10(n_plasmasphere)-np.log10(n_ionosphere))/  \
                                         (radius_height_plasma-radius_height_iono))*(r-radius_height_iono)) + np.log10(n_ionosphere))
                    else : density = 0.

                elif l_shell <= 2.0:
                  
                    density_plasmasphere_equator = density_plasmasphere(r=radius_height_plasma, lambdaE=0, MLT=MLT, a_plasmapause=a_plasmapause, \
                                                                        a_clean=a_clean, nspot=nspot, year=year, month=month, day=day, add_trough=add_trough, troughParam=troughParam) 
                    
                    density_ionosphere_equator = density_ionosphere_electrons_irbem_jul2021(r=radius_height_iono, lambdaE=0, MLT=MLT, k_slide=k_slide, \
                                                                                               nspot=nspot, \
                                                                                               year=year, month=month, day=day)
                    
                    interpolation_density_l_shell = 10**((((np.log10(density_plasmasphere_equator)-np.log10(density_ionosphere_equator))/ \
                                                          (radius_height_plasma-radius_height_iono))*(l_shell-radius_height_iono)) + np.log10(density_ionosphere_equator))
                  
                    magnetic_latitude_field_line_crossing_iono = np.rad2deg(np.arccos(np.sqrt( radius_height_iono/l_shell )))
                    
                    if lambdaE < 0: magnetic_latitude_field_line_crossing_iono = -magnetic_latitude_field_line_crossing_iono
                    if (magnetic_latitude_field_line_crossing_iono < -90 or magnetic_latitude_field_line_crossing_iono > 90): 
                        print('latitude not between -90 and 90')
                  
                    density_iono = density_ionosphere_electrons_irbem_jul2021(r=radius_height_iono, lambdaE=magnetic_latitude_field_line_crossing_iono, MLT=MLT, \
                                                                                 k_slide=k_slide, nspot=nspot, \
                                                                                 year=year, month=month, day=day)
                    
    # logarithmic interpolation
                    density = 10**((((np.log10(interpolation_density_l_shell)-np.log10(density_iono))/ \
                                     (l_shell-radius_height_iono))*(r-radius_height_iono)) + np.log10(density_iono))
#r >= 2.0
            else:
                density = density_plasmasphere(r=r, lambdaE=lambdaE, MLT=MLT, a_plasmapause=a_plasmapause, a_clean=a_clean, nspot=nspot, \
                                               year=year, month=month, day=day, add_trough=add_trough, troughParam=troughParam )

    return density

#################################################################################################################
def density_ionosphere_electrons_irbem_jul2021(r, lambdaE, MLT, k_slide, nspot, year, month, day ):
          
# k_slide is heure in 2 * Universal Time
    float_hour = k_slide / 2.0 # Universal Time
    integer_hour = int(k_slide / 2.0)
    integer_minutes = int((float_hour - integer_hour) * 60.)
            
    r_sm = r
    xlat_degree_sm = lambdaE # in degrees
    xlat_sm = np.deg2rad(xlat_degree_sm)
    MLT_sm = MLT
#      ; local noon 12:00 MLT = geomagnetic longitude 0 degrees
#      ; local dusk 18:00 MLT = geomagnetic longitude 90 degrees
#      ; local midnight 00:00 MLT = geomagnetic longitude 180 degrees
#      ; local dawn 6:00 MLT = geomagnetic longitude 270 degrees
#      ; xlon_degree_sm is from 0 to 360 degrees
    xlon_degree_sm = ( ( MLT_sm * 360. ) / 24. + 180. ) % 360. #; in degrees
    xlon_sm = np.deg2rad(xlon_degree_sm)
      
#      ; r_sm, lambda_sm, MLT_sm are SM coordinates
#      ; transform polar r_sm, lambda_sm, MLT_sm to cartesian x_sm, y_sm, z_sm
            
    x_sm = r_sm * np.cos(xlon_sm) * np.cos(xlat_sm) #; in Earth radii
    y_sm = r_sm * np.sin(xlon_sm) * np.cos(xlat_sm) #; in Earth radii
    z_sm = r_sm * np.sin(xlat_sm) #; in Earth radii

    ##### IRBEM COORD TRANSFORMATION    
    cvals = coord.Coords([[x_sm,y_sm,z_sm]], 'SM', 'car')
    d=str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'T'+str(integer_hour).zfill(2)+':'+str(integer_minutes).zfill(2)+':00'
    cvals.ticks = Ticktock([d], 'ISO') # add ticks
    c = cvals.convert('MAG', 'car')
    x_mag = list(c.x)[0]    
    y_mag = list(c.y)[0]
    z_mag = list(c.z)[0]    
    ##### IRBEM COORD TRANSFORMATION    

    r_mag = np.sqrt(x_mag**2 + y_mag**2 + z_mag**2)

    xlon_mag = np.arctan2(y_mag, x_mag) #; arctan only gives results between -pi and pi so we need to add 2pi when it's negative
    if xlon_mag <= 0: xlon_mag = xlon_mag + 2 * np.pi
    xlon_degree_mag = np.rad2deg(xlon_mag)
    xlat_mag = np.arctan2(z_mag, np.sqrt(x_mag**2 + y_mag**2)) #; atan only gives results between -pi and pi and that is not ok, because it should be between -pi/2 and pi/2
    xlat_degree_mag = np.rad2deg(xlat_mag)
      
#      ; IRI expects a height in km from the surface of the Earth
    hx = (r_mag - 1.0) * 6371.

    altitude = c_double(hx) #; in km
    latitude = c_double(xlat_degree_mag) #    ; IRI expects a latitude between -90 and 90 degrees
    longitude = c_double(xlon_degree_mag) #    ; IRI expects a longitude between 0 and 360 degrees
    iyear=c_int(year)
    imd = month * 100 + day
    mmdd = c_int(imd)
    
#    ; ionosphere corotates with earth
    fhour=c_double(float_hour)     #; Universal Time
    
#    ;C*** THE ALTITUDE LIMITS ARE:  LOWER (DAY/NIGHT)  UPPER        ***
#    ;C***     ELECTRON DENSITY         60/80 KM       1000 KM       ***
#    ;C***     TEMPERATURES              120 KM        2500/3000 KM  ***
#    ;C***     ION DENSITIES             100 KM        1000 KM       ***

    res=c_double(0.0) #; 64 bit
    
    old_dir = os.getcwd()
    os.chdir(iri_path)
    irilib.elecdens7jun2016_(byref(altitude),byref(latitude),byref(longitude),byref(iyear),byref(mmdd),byref(fhour),byref(res))
    os.chdir(old_dir)
    
    density = res.value #      ; electron density per cm**-3
    if density < 0: density=0

    return density

#################################################################################################################
def density_plasmasphere (r, lambdaE, MLT, a_plasmapause, a_clean, nspot, \
                          year, month, day, add_trough, troughParam ):
    """    
      ; input: r is radial distance,
      ; lambda is magnetic latitude in degrees: [-90, 90]
      ; MLT is Magnetic Local Time (Longitude) from 0 to 24
      ; input: a_plasmapause and a_clean (see create_plasmapause)
      ; input: heure, nspot
      if (n_elements(add_trough) le 0) then add_trough=1
      
      if (n_elements(p_plasma) le 0) then p_plasma=1 ; power law coefficient between 0 and 1
      if (n_elements(p_trough) le 0) then p_trough=2 ; power law coefficient between 2 and 3
    """    
    artificial_cut_off=10
    p_plasma=1 
    p_trough=2
    epsilon = 0.001
    n=0.
    
    if np.cos(np.deg2rad(lambdaE)) > epsilon:
        is_inside = inside_plasmasphere(r=r, lambdaE=lambdaE, MLT=MLT, a_plasmapause=a_plasmapause, a_clean=a_clean)
        
#            ;has_plume_points=has_plume_points, a_plume_points=a_plume_points, plume_cartesian=plume_cartesian, $
#            ;    equatorial_plasmapause_cartesian=equatorial_plasmapause_cartesian, $
#            ;parameter_angle_delta=parameter_angle_delta, parameter_radius_epsilon=parameter_radius_epsilon)

        l_shell = r / np.cos(np.deg2rad(lambdaE))**2

        if is_inside == 1:
            density_electron_equatorial = density_electron_plasmasphere_equatorial(r=l_shell, nspot=nspot, year=year, month=month, day=day)    
            n = density_electron_equatorial * (np.cos(np.deg2rad(lambdaE)))**(-2*p_plasma)      
        else:
            if add_trough == 1 and l_shell <= artificial_cut_off:  # outside plasmapause 
                if troughParam == 'CA':      n = trough_CA(MLT,l_shell,lambdaE,p_trough)
                if troughParam == 'Sheeley': n = trough_Sheeley(MLT,l_shell,lambdaE,p_trough)
                if troughParam == 'VAP':     n = trough_VAP(MLT,l_shell,lambdaE,p_trough)

    if n > 1.e10: n=0.

    return n

#################################################################################################################
def trough_VAP(MLT,L,Lat,p_trough):
    
#    print('PlasmaTrough parameters from VAP fitting')
    
    switch_plus3 = {
# new : eq dens + L > Lpp (from SPM simulations) a=2 (p_trough) !!!!        
    ( 0 , 6 ):[3945,  227],
    ( 6 , 12 ):[ 300,  848],
    ( 12 , 18 ):[ 301, 1476],
    ( 18 , 24 ):[1578,  527],
# new : eq dens + L > Lpp (from SPM simulations) a=1 !!!!
    # ( 0 , 6 ):[4655,   33],
    # ( 6 , 12 ):[ 101,  897],
    # ( 12 , 18 ):[5914, 1142],
    # ( 18 , 24 ):[ 102,  612],
# new : eq dens + L > Lpp (from SPM linear equation Kp dependant)    
    # ( 0 , 6 ):[ 986, 5648],
    # ( 6 , 12 ):[ 100, 4193],
    # ( 12 , 18 ):[ 421, 4049],
    # ( 18 , 24 ):[ 185, 1768],
# old : <100cm-3 & L>3    
    # ( 0 , 6 ):[1382,  -76],
    # ( 6 , 12 ):[ 114,  159],
    # ( 12 , 18 ):[ 113,  177],
    # ( 18 , 24 ):[ 104,  103],
                    }

    mlt1,mlt2,[a,b] = switch(switch_plus3,MLT)
    e = -4.0 # old : e = -3
# mlt1p,mlt2p are not used    
    if MLT>=18 : mlt1p,mlt2p,[c,d] = switch(switch_plus3,0) ; d=0.0  
    else:        mlt1p,mlt2p,[c,d] = switch(switch_plus3,MLT+6)
    
    return ( ((a + b * (MLT)) * L**(e) + (1.0 - np.exp(-(L-2.)/10.0))) * (np.cos(np.deg2rad(Lat)))**(-2*p_trough) ) * ((mlt2-MLT)/(mlt2-mlt1)) + \
           ( ((c + d * (MLT)) * L**(e) + (1.0 - np.exp(-(L-2.)/10.0))) * (np.cos(np.deg2rad(Lat)))**(-2*p_trough) ) * ((MLT-mlt1)/(mlt2-mlt1))
 

def switch(table, val):
    for (mlt1, mlt2) in table:
        if mlt1 <= val < mlt2:
            return mlt1, mlt2, table[(mlt1, mlt2)]     

#################################################################################################################
def trough_CA(MLT,l_shell,lambdaE,p_trough):
    
#    print('PlasmaTrough parameters from CA')
    
    if MLT >= 0 and MLT < 6:
#                  ; from midnight to morning (dawn)
        n = ((5800.+300.*MLT)*l_shell**(-4.5)+(1.-np.exp(-(l_shell-2.)/10.0)))*(np.cos(np.deg2rad(lambdaE)))**(-2*p_trough)

    if MLT >= 6 and MLT < 15:
#                  ; from morning (dawn) to afternoon
        n = ((-800.+1400.*MLT)*l_shell**(-4.5)+(1.-np.exp(-(l_shell-2.)/10.0)))*(np.cos(np.deg2rad(lambdaE)))**(-2*p_trough)

    if MLT >= 15 and MLT <= 24:
#                  ; from afternoon to midnight
# Changed on April 29 2020 to ensure continuity at MLT=15 and MLT=24
#                    n = (((-800.+1400.*15.)*l_shell**(-4.5)+(1.-np.exp(-(l_shell-2.)/10.0)))*(np.cos(np.deg2rad(lambdaE)))**(-2*p_trough)) * ((24.-MLT)/(24.-15.)) + \
        n = (((44200.-1600.*MLT)*l_shell**(-4.5)+(1.-np.exp(-(l_shell-2.)/10.0)))*(np.cos(np.deg2rad(lambdaE)))**(-2*p_trough)) * ((24.-MLT)/(24.-15.)) + \
            ((5800.*l_shell**(-4.5)+(1.-np.exp(-(l_shell-2.)/10.0)))*(np.cos(np.deg2rad(lambdaE)))**(-2*p_trough)) * ((MLT-15.)/(24.-15.))

    return n

#################################################################################################################
def trough_Sheeley(MLT,l_shell,lambdaE,p_trough):
    
#    print('PlasmaTrough parameters from Sheeley')
 
    n = 124 * (3/l_shell)**4 + 36 * (3/l_shell)**3.5 * np.cos(((MLT) - (7.7 * (3/l_shell)**2 + 12)) * np.pi/12)
#    npm = 78 * (3/l_shell)**4.72 + 17 * (3/l_shell)**3.75 * np.cos(((MLT)-22) * np.pi/ 12)   # minus/plus to define region

    return n

#################################################################################################################
def inside_plasmasphere (r, lambdaE, MLT, a_plasmapause, a_clean):
    """
    ;has_plume_points=has_plume_points, a_plume_points=a_plume_points, plume_cartesian=plume_cartesian, $
    ;equatorial_plasmapause_cartesian=equatorial_plasmapause_cartesian, $
    ;parameter_angle_delta=parameter_angle_delta, parameter_radius_epsilon=parameter_radius_epsilon)
    
    ; input: r is radial distance
    ; input: lambda is Magnetic Latitude in degrees [-90, 90]
    ; input: MLT is Magnetic Local Time (Longitude)
    ; input: a_plasmapause and a_clean (calculated in create_plasmapause)
    ; (are in Cartesian Coordinates and are the plasmapause points in the equatorial plane)
    ; output: is_inside=1 inside; is_inside=0 outside
    """    
    epsilon = 0.001
    is_inside = 1
      
    theta = (MLT*2*np.pi)/ 24.0
    x = r*np.cos(theta) ; y = r*np.sin(theta)
    
    if np.cos(np.deg2rad(lambdaE)) > epsilon:
        inside_plasmapause = inside( x, y, list(np.transpose(a_plasmapause[0,:]*(np.cos(np.deg2rad(lambdaE)))**2)), \
                                           list(np.transpose(a_plasmapause[1,:]*(np.cos(np.deg2rad(lambdaE)))**2)) )
        
        if not inside_plasmapause :
            inside_plasmaplume = inside( x, y, list(np.transpose(a_clean[:,0]*(np.cos(np.deg2rad(lambdaE)))**2)), \
                                               list(np.transpose(a_clean[:,1]*(np.cos(np.deg2rad(lambdaE)))**2)) )
            if not inside_plasmaplume : is_inside = 0

    else: is_inside = 0 # points on the z axis are outside the plasmapause
    
    return is_inside
#####################################################################################
def inside ( x, y, px, py):
    """
     ;  x - The x coordinate of the point.
     ;  y - The y coordinate of the point.
     ; px - The x coordinates of the polygon.
     ; py - The y coordinates of the polygon.
     ;
     ; The return value of the function is 1 if the point is inside the
     ; polygon and 0 if it is outside the polygon.
     EDITH: the sum does not give exactly the same value as IDL !!!!!!
    """
    sx = len(px)
    sy = len(py)
    if sx != 1: NX=sx 
    else: return -1    # Error if px not a vector
    if sy != 1: NY=sy 
    else: return -1    # Error if py not a vector
    
    if NX == NY : N = NX 
    else: return -1        # Incompatible dimensions
    
    tmp_px = px + [px[0]]                             # Close Polygon in x
    tmp_py = py + [py[0]]                             # Close Polygon in y
     
    i  = np.arange(N,dtype=int)                       # Counter (0:NX-1)
    ip = np.arange(1,N+1,dtype=int)                   # Counter (1:nx)    
    
    X1=[]; Y1=[]; X2=[]; Y2=[]
    
    for ii in list(i):
        X1.append(tmp_px[ii]  - x)
        Y1.append(tmp_py[ii]  - y)
    for iip in list(ip):
        X2.append(tmp_px[iip] - x)
        Y2.append(tmp_py[iip] - y)     
    
    X1=np.array(X1)
    Y1=np.array(Y1)
    X2=np.array(X2)
    Y2=np.array(Y2)
        
    dp = X1*X2 + Y1*Y2                               # Dot-product
    cp = X1*Y2 - Y1*X2                               # Cross-product
    theta = np.arctan2(cp,dp)
    
    if np.abs(theta.sum()) > np.pi : return 1 
    else: return 0
################################################################################################################

def compute(datetimes, to_np: bool = False):  
#    if not isinstance(datetimes[0], list):
    datetimes = [datetimes]
    cdf_time = []
    for d in datetimes:
#        print(d)
        unix_seconds = datetime.datetime(d[0], d[1], d[2], d[3], d[4], d[5]).replace(tzinfo=timezone.utc).timestamp()
        if len(d) == 7:
            remainder_seconds = (d[6]/1000.0)
            astrotime = Time(unix_seconds, remainder_seconds, format='unix', precision=9)
            cdf_time.append(astrotime.cdf_epoch)
        if len(d) == 9:
            remainder_seconds = (d[6]/1000.0) + (d[7]/1000000.0) + (d[8]/1000000000.0)
            astrotime = Time(unix_seconds, remainder_seconds, format='unix', precision=9)
            cdf_time.append(astrotime.cdf_tt2000)
        if len(d) == 10:
            remainder_seconds = (d[6]/1000.0) + (d[7]/1000000.0) + (d[8]/1000000000.0) + (d[9]/1000000000000.0)
            astrotime = Time(unix_seconds, remainder_seconds, format='unix', precision=9)
            cdf_time.append(astrotime.cdf_epoch16)
    if to_np:
        return np.array(cdf_time)
    else:
        return cdf_time
#################################################################################################################
def density_electron_plasmasphere_equatorial (r, nspot, year, month, day): 
#  ; cm**-3
    
#      ; date is 'the number of the day'
      date = datetime.datetime(year,month,day).timetuple().tm_yday 
      
      density_electron =  10**( (-0.3145*float(r)+3.9043) + \
                                (0.15*(np.cos(2.*np.pi*(float(date)+9.)/365.) - 0.5*np.cos(4.*np.pi*(float(date)+9.)/365.) + \
                                 0.00127*float(nspot)-0.0635))* np.exp(-(float(r)-2.0)/(1.5)) )
        
      return density_electron
