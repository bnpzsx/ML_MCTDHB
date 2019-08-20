# coding=utf-8
import numpy as np
#import pylab as pyl
#import os
#import sys
#import math
#import copy
import QDTK
import QDTK.Operatorb as Operb
import QDTK.Operator  as Oper
#import QDTK.Tools.Conversion as conv
import QDTK.PES.Potfit as pfit
import Parameters as Par

########## system parameters #####################################################
num_spec = 1            # number of species

n     = Par.n
m     = Par.m
mx    = Par.mx
my    = Par.my
npx   = Par.npx
npy   = Par.npy
Qx_dvr = Par.Qx_dvr
Qy_dvr = Par.Qy_dvr

########## operator generaltion###################################################
## define Hamiltonians ##
def returnOperator():
    OPER = Operb.Operatorb()
    OPER.define_dofs_and_grids((1,1),(Qx_dvr,Qy_dvr))    
    # Labels
    OPER.addLabel('dq2_x',Operb.OTerm(Qx_dvr.d2dvr))  
    OPER.addLabel('dq2_y',Operb.OTerm(Qy_dvr.d2dvr))  
    OPER.addLabel('q^2_x',Operb.OTerm(Qx_dvr.x**2))  
    OPER.addLabel('q^2_y',Operb.OTerm(Qy_dvr.x**2))  
    OPER.addLabel('q_x',Operb.OTerm(Qx_dvr.x))
    OPER.addLabel('q_y',Operb.OTerm(Qy_dvr.x))
    OPER.addLabel('a_0',Operb.OCoef(1.0))


    tab="""
a_0          |1 dq2_x
a_0          |2 dq2_y
a_0          |1 q^2_x
a_0          |2 q^2_y
a_0          |1 q_x  |2 q_y
"""
    OPER.readTableb(tab)
    return OPER

###############################################################################
op = returnOperator()     
op.createOperatorFileb('dynamic-operator')
