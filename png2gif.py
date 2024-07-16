#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 10:07:41 2020

@author: edithb
"""

from PIL import Image
import glob,os

os.getcwd()
os.chdir('../out/spm/')
dire=os.getcwd()+'/'
print(dire)
# Create the frames
frames = []
imgs = glob.glob(dire+'/*.png')
print(imgs)
for i in sorted(imgs):
    new_frame = Image.open(i)
    frames.append(new_frame)

# Save into a GIF file that loops forever
# a=imgs[0].split('/'); b=imgs[-1].split('/')
# for i in range(len(a)):
#     if a[i]=='spm': 
#         d1=a[i+1]+a[i+2]+a[i+3]
#         d2=b[i+1]+b[i+2]+b[i+3]
#     else:pass

#frames[0].save(dire+d1+'-'+d2+'.gif', format='GIF',
frames[0].save(dire+'test.gif', format='GIF',
               append_images=frames[1:],
               save_all=True,
               duration=100, loop=0)
