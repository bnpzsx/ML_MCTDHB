'''
对一个指定的二维波色子系统,定义每个node的展开系数
'''
import numpy as np
from constant import m,n,mx,my,npx,npy,type_complex
from node import node
from utils import *
from reader import psi_reader
# tree:
#                        (+)
#                        ||| n
#                         O -m
#                       /   \
#                      /     \
#                     O -mx   O -my
#                     |       |
#                     O -npx  O -npy
class tree:
    '定义每个node的展开系数'
    def __init__(self):
        l1=self.tree=(node(1,1,1,n))
        l2=l1.addsubnode(node(2,m))
        l3x=l2.addsubnode(node(3,mx))
        l3y=l2.addsubnode(node(3,my))
        l4x=l3x.addsubnode(node(4,npx))
        l4y=l3y.addsubnode(node(4,npy))
        self.tree.makewavefunction()
        
        l=self.layers={i:[] for i in range(1,5)} #按每层索引node
        for i in self.tree:
            l[i.layer].append(i)
        s=0
        for i in self.tree:
            if i.num_subnode!=0: # 叶子节点没有波函数
                s+=i._psi.size
        self.psilen=s

    def init_from_psi(self,path='restart.ini'):
        '使用QDTK的restart文件给波函数赋初值'
        f=psi_reader(path)
        psi=f.psi[0]
        assert self.psilen==len(psi),'psi file does not match the tree.'
        n=0 #psi 读取位置
        p_node=[] #存放primitive node
        for i in self.tree:
            if i.num_subnode>0:
                if i.subnode[0].num_subnode==0:
                    p_node.append(i)
                    continue #跳过primitive node
            else:
                continue #跳过叶子节点
            l=i._psi.size
            i._psi=np.reshape(psi[n:n+l],(i.num_SPFs,i.dim_SPF))
            n+=l
        for i in p_node:
            l=i._psi.size
            i._psi=np.reshape(psi[n:n+l],(i.num_SPFs,i.dim_SPF))
            n+=l

if __name__=='__main__':
    t=tree()
    r=t.tree
    r.print()
    t.init_from_psi('restart.ini')
