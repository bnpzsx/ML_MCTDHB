'''
对一个指定的二维波色子系统,定义每个node的展开系数
'''
import numpy as np
from constant import m,n,m1,mx,my,npx,npy,type_complex
from utils import *
# tree:
#                         O 
#                         |
#                        (+) -m    
#                        ||| n      
#                         O -m1     
#                       /   \
#                      /     \
#                     O -mx   O -my
#                     |       |
#                     O -npx  O -npy
class tree:
    '定义每个node的展开系数'
    def __init__(self):
        l=self.layers={}
        zeros=lambda shape:np.zeros(shape,type_complex)
        l[1]=zeros(m)
        l[2]=zeros((m,num_combination(m1+n-1,n)))
        l[3]=zeros((m1,mx,my))
        l[4]=np.array([zeros((mx,npx)),zeros((my,npy))])
        s=0
        for k in l:
            for i in l[k]:
                s+=i.size
        self.psilen=s

if __name__=='__main__':
    t=tree()
    l=t.layers
