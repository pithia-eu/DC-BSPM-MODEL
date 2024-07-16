#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 15:42:43 2020

@author: edithb
"""
import numpy as np
import statistics as stat

####################################################################################################
def create_plasmapause (equatorial_plasmapause_cartesian):
    
#; input: equatorial_plasmapause_cartesian
#; output: equatorial_plasmapause_polar
#; equatorial_plasmapause_polar(0) is theta
#; equatorial_plasmapause_polar(1) is r
#; equatorial_plasmapause_cartesian(0) is x
#; equatorial_plasmapause_cartesian(1) is y
#; optional input: parameters
#;
#; necessary for function inside_plasmasphere:
#; output: a_plasmapause (the plasmapause points without plume points) (cartesian coordinates)
#; output: a_clean (all the plasmapause points without outliers) (cartesian coordinates)
    
    print('in create_plasmapause ...')

    j = np.where(equatorial_plasmapause_cartesian[0,:] > -1.0E21)
    nj = len(j[0])
       
    equatorial_plasmapause_polar = np.empty([2,nj])
  
#;  equatorial_plasmapause_cartesian(0,*) = -equatorial_plasmapause_cartesian(0,*)
#;  equatorial_plasmapause_cartesian(1,*) = -equatorial_plasmapause_cartesian(1,*)

#;;;DATA:POLAR COORDINATES;;;
    for m in range(nj):
        equatorial_plasmapause_polar[0,m] = np.arctan( equatorial_plasmapause_cartesian[1,j[0][m]]/equatorial_plasmapause_cartesian[0,j[0][m]] )
        equatorial_plasmapause_polar[1,m] = np.sqrt( ( equatorial_plasmapause_cartesian[0,j[0][m]])**2 + (equatorial_plasmapause_cartesian[1,j[0][m]])**2 )
#;;;ATAN functions only between -pi/2 and pi/2;;;
        if equatorial_plasmapause_cartesian[1,j[0][m]] >= 0:
            if equatorial_plasmapause_cartesian[0,j[0][m]] <= 0 : equatorial_plasmapause_polar[0,m] = equatorial_plasmapause_polar[0,m]+np.pi

        if equatorial_plasmapause_cartesian[1,j[0][m]] < 0:
            if equatorial_plasmapause_cartesian[0,j[0][m]] <  0 : equatorial_plasmapause_polar[0,m] = equatorial_plasmapause_polar[0,m]-np.pi
  
#; [0,2*pi]
    for m in range(nj):
        if equatorial_plasmapause_polar[0,m] < 0: equatorial_plasmapause_polar[0,m] = equatorial_plasmapause_polar[0,m]+2*np.pi

#; equatorial_plasmapause_polar are calculated

#;reversing again the a array
#;  equatorial_plasmapause_cartesian(0,*)=-equatorial_plasmapause_cartesian(0,*)
#;  equatorial_plasmapause_cartesian(1,*)=-equatorial_plasmapause_cartesian(1,*)
  
    b_plume_points = np.zeros(2) ; a_plume_points = np.zeros(2)
    has_plume_points = 0
    b_no_plume_points = np.zeros(2) ; a_no_plume_points = np.zeros(2)
  
#; the plasmapause without outliers is calculated
    a_clean = remove_plume_outliers(equatorial_plasmapause_cartesian[:,j[0][:]], equatorial_plasmapause_polar)
    
    for s in range(nj):
        is_outlier     = point_is_outlier(equatorial_plasmapause_polar, index=s)
        is_plume_point = point_is_plume_point(equatorial_plasmapause_polar, s=s)
#        print(s, is_outlier, is_plume_point)

        if ( (is_outlier == 1) or (is_plume_point == 1) ) :
            b_plume_points = np.row_stack([b_plume_points, [equatorial_plasmapause_polar[0,s], equatorial_plasmapause_polar[1,s]]])
            a_plume_points = np.row_stack([a_plume_points, [equatorial_plasmapause_cartesian[0,j[0][s]], equatorial_plasmapause_cartesian[1,j[0][s]]]])
            has_plume_points = 1
        else:
            b_no_plume_points = np.row_stack([b_no_plume_points, [equatorial_plasmapause_polar[0,s], equatorial_plasmapause_polar[1,s]]])
            a_no_plume_points = np.row_stack([a_no_plume_points, [equatorial_plasmapause_cartesian[0,j[0][s]], equatorial_plasmapause_cartesian[1,j[0][s]]]])
    
    if has_plume_points == 1:
        b_plume_points = b_plume_points[1:, :]
        a_plume_points = a_plume_points[1:, :]      

    b_no_plume_points = b_no_plume_points[1:, :]
    a_no_plume_points = a_no_plume_points[1:, :]
      
    number_a_no_plume_points = len(a_no_plume_points)
#    print('number_a_no_plume_points =',number_a_no_plume_points)

#    ; the plasmapause points without plume points are calculated
    a_plasmapause = np.empty([2,number_a_no_plume_points])
    a_plasmapause[0] = a_no_plume_points[np.argsort(b_no_plume_points[:,0]),0]
    a_plasmapause[1] = a_no_plume_points[np.argsort(b_no_plume_points[:,0]),1]

    return equatorial_plasmapause_polar, a_clean, a_plasmapause 

#################################################################################################################
def distance (x0, y0, x1, y1):
    d = np.sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1))
    return d

####################################################################################################
def remove_plume_outliers (a, b):
    """    
      ; a point is a plume outlier if it is not a plasmapause point
      ; and (if the distance with the previous point is too large
      ; or if the distance with the next point is too large)
      ;
      ; plume_outlier_exclusion=0, remove no outliers.
      ; plume_outlier_exclusion=1, removes some outliers.
      ; plume_outlier_exclusion>1, remove many outliers.
      ; default: plume_outlier_exclusion=6
    
        if len(plume_outlier_exclusion) <= 0 : plume_outlier_exclusion=6
        if plume_outlier_exclusion == 0 : return a 
        else: exclusion = plume_outlier_exclusion
    """
    exclusion=6 

    delta = 2.0 / exclusion 
    na = len(a[0, :])
    a_clean = np.arange(12) #; traject.dat has 12 columns

#  ;i=0    
    add=1 
    
    if (point_is_outlier(b,0)     == 1 or \
        point_is_plume_point(b,0) == 1) :
        
#    ; point is not a plasmapause point, but a plume point    
        if distance(a[0, na-1], a[1, na-1], a[0, 0], a[1, 0]) > delta : add = 0
        if distance(a[0, 0   ], a[1, 0   ], a[0, 1], a[1, 1]) > delta : add = 0
        
    if add == 1 : a_clean = np.row_stack([a_clean, a[:,0]])
    
    for i in range(1,na-1):
        add=1
    
        if (point_is_outlier(b,i)     == 1 or \
            point_is_plume_point(b,i) == 1) :

#      ; point is not a plasmapause point, but a plume point      
            if distance(a[0, i-1], a[1, i-1], a[0, i  ], a[1, i  ]) > delta: add = 0
            if distance(a[0, i  ], a[1, i  ], a[0, i+1], a[1, i+1]) > delta: add = 0

        if add == 1: a_clean = np.row_stack([a_clean, a[:,i]])

#  ;i=na-1
    add=1         
    if (point_is_outlier(b, index=na-1) == 1 or \
        point_is_plume_point(b, s=na-1) == 1) :

#  ; point is not a plasmapause point, but a plume point    
        if distance(a[0, na-2], a[1, na-2], a[0, na-1], a[1, na-1]) > delta: add = 0
        if distance(a[0, na-1], a[1, na-1], a[0, 0   ], a[1, 0   ]) > delta: add = 0
            
    if add == 1: a_clean = np.row_stack([a_clean, a[:,na-1]])
 
    a_clean = a_clean[1:,:]
    n_a_clean = len(a_clean)
    a_clean=np.transpose(a_clean)
    a_smooth = np.arange(12) #; traject.dat has 12 columns
  
    dis = 2.5
  
#  ; k=0
    add = 1
  
    if distance(a_clean[0, n_a_clean-1], a_clean[1, n_a_clean-1], a_clean[0, 0], a_clean[1, 0]) > delta and \
       distance(a_clean[0, 0          ], a_clean[1, 0          ], a_clean[0, 1], a_clean[1, 1]) > delta :
        
        radius1 = np.sqrt( a_clean[0, n_a_clean-1]*a_clean[0, n_a_clean-1] + a_clean[1, n_a_clean-1]*a_clean[1, n_a_clean-1] )
        radius2 = np.sqrt( a_clean[0, 1          ]*a_clean[0, 1          ] + a_clean[1, 1          ]*a_clean[1, 1          ] )
    
        radius = max([radius1, radius2])
    
        k_radius = np.sqrt( a_clean[0, 0]*a_clean[0, 0] + a_clean[1, 0]*a_clean[1, 0] )
        if k_radius > radius + dis : add = 0

    a_clean=np.transpose(a_clean)
    if add == 1 : a_smooth = np.row_stack([a_smooth, a_clean[0,:]])
    
    a_clean=np.transpose(a_clean)
    
    for k in range(1, n_a_clean-1):
        add = 1
        if distance(a_clean[0, k-1], a_clean[1, k-1], a_clean[0, k  ], a_clean[1, k  ] > delta) and \
           distance(a_clean[0, k  ], a_clean[1, k  ], a_clean[0, k+1], a_clean[1, k+1] > delta) :
          
            radius1 = np.sqrt( a_clean[0, k-1]*a_clean[0, k-1] + a_clean[1, k-1]*a_clean[1, k-1] )
            radius2 = np.sqrt( a_clean[0, k+1]*a_clean[0, k+1] + a_clean[1, k+1]*a_clean[1, k+1] )
              
            radius = max([radius1, radius2])
              
            k_radius = np.sqrt( a_clean[0, k]*a_clean[0, k] + a_clean[1, k]*a_clean[1, k] )
            if k_radius > radius + dis: add = 0
            
        if add == 1 : a_smooth = np.row_stack([a_smooth, a_clean[:,k]])
    

#  ; k=n_a_clean-1
    add = 1
    
    if distance(a_clean[0, n_a_clean-2], a_clean[1, n_a_clean-2], a_clean[0, n_a_clean-1], a_clean[1, n_a_clean-1]) > delta and \
       distance(a_clean[0, n_a_clean-1], a_clean[1, n_a_clean-1], a_clean[0, 0          ], a_clean[1, 0          ]) > delta: 
    
        radius1 = np.sqrt( a_clean[0, n_a_clean-2]*a_clean[0, n_a_clean-2] + a_clean[1, n_a_clean-1]*a_clean[1, n_a_clean-1] )
        radius2 = np.sqrt( a_clean[0, n_a_clean-1]*a_clean[0, n_a_clean-1] + a_clean[1, 0          ]*a_clean[1, 0          ] )
        
        radius = max([radius1, radius2])
        
        k_radius = np.sqrt( a_clean[0, n_a_clean-1]*a_clean[0, n_a_clean-1] + a_clean[1, n_a_clean-1]*a_clean[1, n_a_clean-1] )
        if k_radius > radius + dis: add = 0
    
    a_clean=np.transpose(a_clean)    
    if add == 1 : a_smooth = np.row_stack([a_smooth, a_clean[n_a_clean-1,:]])
    
    return a_smooth[1:, :]
    
#################################################################################################################
def point_is_outlier (b, index):
    """    
      ; outlier_exclusion=0, remove no outliers.
      ; default outlier_exclusion=1, removes some outliers.
      ; outlier_exclusion>1, remove many outliers.
    
    if len(plasmapause_outlier_exclusion) <= 0 : plasmapause_outlier_exclusion=1
    if plasmapause_outlier_exclusion == 0 : return 0 
    else : exclusion = plasmapause_outlier_exclusion
    """
    exclusion = 1
    
    delta_radius = 20.0 / exclusion
    mean_radius = np.mean(b[1,:])
    variance_radius = stat.variance(b[1,:]) 

    maximum_radius = mean_radius + delta_radius * variance_radius
 
    is_outlier = 0
    if b[1,index] >= maximum_radius : is_outlier = 1
    
#    if index == 91:
#        print('edith en point_is_outlier: maximum_radius ', maximum_radius, mean_radius, delta_radius, variance_radius, b[1,index], is_outlier)
    
    return is_outlier

#################################################################################################################
def point_is_plume_point (b, s):
    """
  ; we consider the traject.dat points having a random order.

  ; returns 1 if the plasmapause point with index 's' is a point of a plume.
  ; returns 0 if the plasmapause point is not a point of the plume.

  ; a plume only exists at a distance larger than 6

    if len(parameter_angle_delta) <= 0 : angle_delta=0.2 
    else: angle_delta = parameter_angle_delta
    if len(parameter_radius_epsilon) <= 0 : radius_epsilon=0.5 
    else: radius_epsilon = parameter_radius_epsilon
    """

    angle_delta=0.2 ; radius_epsilon=0.5

    nb = len(b[0,:])

    is_plume_point = 0
      
    first_angle  = b[0,s] - angle_delta
    second_angle = b[0,s] + angle_delta
    
    first_radius = b[1,s] - radius_epsilon
    if first_radius < 0. : first_radius = 0.

    second_radius = b[1,s] + radius_epsilon
  
    for n in range(nb-1):
        if n != s:
            if (b[0,n] >= first_angle and b[0,n] <= second_angle):
                if (b[1,n] < first_radius or b[1,n] > second_radius):
#          ; plasmapause has large plume
#          ; point with index 's' is a plume point is radius(s) > radius(n); b(1,s) > b(1,n)
                    if b[1,s] > b[1,n]: is_plume_point = 1 ; break
      
            if (b[0,n]-2*np.pi >= first_angle and b[0,n]-2*np.pi <= second_angle):
                if (b[1,n] < first_radius or b[1,n] > second_radius):
    #          ; plasmapause has large plume
    #          ; point with index 's' is a plume point is radius(s) > radius(n); b(1,s) > b(1,n)
                    if b[1,s] > b[1,n]: is_plume_point = 1 ; break
          
            if (b[0,n]+2*np.pi >= first_angle and b[0,n]+2*np.pi <= second_angle):
                if (b[1,n] < first_radius or b[1,n] > second_radius):
    #          ; plasmapause has large plume
    #          ; point with index 's' is a plume point is radius(s) > radius(n); b(1,s) > b(1,n)
                    if b[1,s] > b[1,n]: is_plume_point = 1 ; break

    return is_plume_point

##########################################################################################################
