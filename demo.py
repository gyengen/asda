#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Monday Nov 25 2019

Name: demo.py

Purpose: A demo of ASDA

@author: Jaijia Liu at University of Sheffield
"""
__author__ = 'Jiajia Liu'
__copyright__ = 'Copyright 2019, The Solar Physics and Space Plasma ' + \
                'Research Center (SP2RC)'
__license__ = 'GPLv3'
__date__ = '2019/11/25'
__maintainor__ = 'Jiajia Liu'
__email__ = 'jj.liu@sheffield.ac.uk'


import numpy as np
from datetime import datetime
import asda


# def artificial_vortex():
#     '''
#     Generate an artificial vortex using the Lamb_Oseen class in asda, then
#     perform the vortex detection and visualisation
#     '''
beg_time = datetime.today()
# Generate an artificial vortex
vmax = 2.0 # rotating speed
rmax = 50 # radius
ratio = 0.2 # ratio of expanding speed over rotating speed
lo = asda.Lamb_Oseen(vmax=vmax, rmax=rmax, ratio_vradial=ratio)
# Generate vx and vy
vx, vy = lo.get_vxvy(x_range=[-100, 100], y_range=[-100, 100])
# Visualise it
lo.visual_vortex()

# perform vortex detection
gamma = lo.gamma_values()
lo.visual_gamma(gamma2=True)
center_edge = lo.center_edge()
(ve, vr, vc, ia) = lo.vortex_property()

# time used
end_time = datetime.today()
print('Time used ', end_time-beg_time)

print('Detected vortex:')
print('Center Location', cetner_edge['center'])


# vel_file = 'vel_demo.npz'  # file in which velocity field will be stored
# vxvy = np.load(vel_file)
# vx = vxvy['vx']
# vy = vxvy['vy']

# # Perform swirl detection
# factor = 1
# # Gamma1 and Gamma2 
# (gamma1, gamma2) = gamma_values(vx, vy, factor=factor)
# # Store tau1 and tau2
# vcimageout((gamma1, gamma2), gamma_file)
# # Determine Swirls
# center, edge, points, peak, radius = center_edge(gamma1, gamma2,
#                                                  factor=factor)
# # Properties of Swirls
# ve, vr, vc, ia = vortex_property(center, edge, points, vx, vy,
#                                  data0)
# # Save results
# np.savez('vortex_demo.npz', center=center, edge=edge,
#          points=points, peak=peak, radius=radius, ia=ia,
#          ve=ve, vr=vr, vc=vc)

# # Caculate time consumption
# end_time = datetime.today()
# print('Time used ', end_time-beg_time)

# if __name__ == '__main__':
#     artificial_vortex()
