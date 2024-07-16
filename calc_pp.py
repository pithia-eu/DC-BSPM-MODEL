#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 16:47:17 2020

@author: edithb
"""
import os, subprocess

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