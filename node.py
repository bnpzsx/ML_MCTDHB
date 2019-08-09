import numpy as np
from constant import type_complex
from utils import num_combination

class node:
    def __init__(self,layer=1,num_SPFs=1,symmetry=0,num_particles=1):
        '初始化单个node'
        assert symmetry in [0,1,-1],'Symmetry must be 0, 1 or -1 '
        self.topnode=None
        self.parent=None
        self.layer=layer
        self.symmetry=symmetry #(0 -> no, 1 -> bos., -1 -> ferm. sym.)
        self.num_particles=num_particles #only for symmetry!=0
        self.num_SPFs=num_SPFs
        self.subnode=[]
        self.subdim=[] #SPFs of subnodes
        self.num_subnode=0 # 0 leaf node
        self.index=0 #which node respect to parent node
        self.dim_SPF=0 #length of a SPF
        self._psi=[] #wavefunction of this node
        self.psi=[] #wavefunction of this node and all subnodes
        
    def addsubnode(self,node):
        '添加子node'
        self.subnode.append(node)
        self.subdim.append(node.num_SPFs)
        node.parent=self
        node.topnode=self if self.topnode is None else self.topnode
        node.index=self.num_subnode
        self.num_subnode+=1
        return node

    def makewavefunction(self):
        '定义每一层波函数(展开系数)'
        assert self.psi==[],'The wavefunction has been defined already'
        zeros=lambda shape:np.zeros(shape,type_complex)
        if self.subnode:
            if self.symmetry==0:
                self.dim_SPF=np.prod(self.subdim)
                self._psi=zeros((self.num_SPFs,self.dim_SPF))
            elif self.symmetry==1:
                m=self.num_particles
                n=self.subdim[0]+m-1
                self.dim_SPF=num_combination(n,m)
                self._psi=zeros((self.num_SPFs,self.dim_SPF))
            else:
                raise OSError('Method for fermionic system is not implemented')
            self.psi.append(self._psi)
            subpsi=[]
            for i in self.subnode:
                i.makewavefunction()
                if i.subnode:
                    subpsi.append(i.psi)
            self.psi.append(subpsi)

    def print(self,prefix='├'):
        '打印树图'
        print('%sl%d.%d SPF:%d type:%d psi:%d'%(prefix,self.layer,self.index,self.num_SPFs,
                                         self.symmetry,np.size(self._psi)))
        for i in self.subnode:
            i.print('│'+prefix)
    def __iter__(self):
        '使tree可以被迭代'
        yield self
        for i in self.subnode:
            yield from i
