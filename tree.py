'''
对一个指定的二维波色子系统,定义每个node的展开系数
'''
import numpy as np
from constant import m,n,m1,mx,my,npx,npy,type_complex
from node import node
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
        l0=self.tree=node(0)
        l1=l0.addsubnode(node(1,m,l0,l0,1,n))
        l2=l1.addsubnode(node(2,m1,l0,l1))
        l3x=l2.addsubnode(node(3,mx,l0,l2))
        l3y=l2.addsubnode(node(3,my,l0,l2))
        l4x=l3x.addsubnode(node(4,npx,l0,l3x))
        l4y=l3y.addsubnode(node(4,npy,l0,l3y))
        self.tree.makewavefunction()
        
        l=self.layers={i:[] for i in range(5)} #按每层索引node
        for i in self.tree:
                l[i.layer].append(i)
        s=0
        for k in l:
            for i in l[k]:
                for j in i._psi:
                    s+=j.size
        self.psilen=s

if __name__=='__main__':
    t=tree()
    r=t.tree
    r.print()
    l=t.layers
