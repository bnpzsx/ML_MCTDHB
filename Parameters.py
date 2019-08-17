# coding=utf-8
import numpy as np
import math
#import QDTK
from sdvr_new import Sindvr
##################################################################################
################################ input parameters ################################
# tree:                         
#                        (+)     
#                        ||| n      
#                         O -m     
#                       /   \
#                      /     \
#                     O -mx   O -my
#                     |       |
#                     O -npx  O -npy

n  = 2
m  = 2

mx = 2
my = 2

###############################################################################
# DVRs in x- and y- direction
x_min = -5.0
x_max =  5.0
npx   =  10

y_min = -5.0
y_max =  5.0
npy   =  10

#Qx_dvr = QDTK.sdvr(npx,x_min,x_max)
#Qy_dvr = QDTK.sdvr(npy,y_min,y_max)
#print(Qx_dvr)
Qx_dvr = Sindvr(npx,x_min,x_max)
Qy_dvr = Sindvr(npy,y_min,y_max)
#print(Qx_dvr)
