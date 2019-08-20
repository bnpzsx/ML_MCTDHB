# coding=utf-8
import numpy as np
from utils import *
import Wavefunction_2d as wfn
import Parameters as Par
###############################################################################
# DVRs in x- and y- direction
Qx_dvr = Par.Qx_dvr
Qy_dvr = Par.Qy_dvr

###############################################################################
# tree:                         
#                        (+)     
#                        ||| n      
#                         O -m     
#                       /   \
#                      /     \
#                     O -mx   O -my
#                     |       |
#                     O -npx  O -npy

n     = Par.n
m     = Par.m
mx    = Par.mx
my    = Par.my
npx   = Par.npx
npy   = Par.npy

tape = (-10,
         n, 1, m, 
        -1, 1,
         2, 0, mx, my,
        -1, 1,
         1, 0, npx,
         0,-1, 2,   
         1, 0, npy,  
        -2)
###############################################################################
# initial conditions:
# initial number states of the species:
ns = np.zeros(m,int)

ns[0] = n-1
ns[m-1] = 1

spfx1 = np.exp(-0.5*(Qx_dvr.x)**2)
spfx1 = spfx1/np.sqrt(sum(spfx1*spfx1))
spfx2 = np.exp(-0.5*(Qx_dvr.x)**2)*Qx_dvr.x
spfx2 = spfx2/np.sqrt(sum(spfx2*spfx2))
spfy1 = Qy_dvr.x/np.sqrt(sum(Qy_dvr.x**2))
spfy2 = spfy1
#spfy2 = Qy_dvr.x**2/np.sqrt(sum(Qy_dvr.x**4))

#必须给出所有spf，并将两个自由度分别打包在一个列表
spfx = [spfx1,spfx2]
spfy = [spfy1,spfy2]

SPF1  = [spfx1,spfy1]
SPF2  = [spfx2,spfy2]

startSPF   = [SPF1,SPF2]
###############################################################################
# accuracies for regarding state as normalized / states as orthonormal
eps_norm = 10**(-13)
eps_over = 10**(-13)
###############################################################################
wfn = wfn.Wavefunction(tape=tape)
#wfn = QDTK.Wfn(tape=tape)
wfn.init_coef_sing_spec_B(ns,startSPF,spfx,spfy)
wfn.createWfnFile('restart.ini')

