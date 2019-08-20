# coding=utf-8
import math
import numpy as np

#class Dvr(object):
class Dvr:#python3
    """
    DVR顶级类，属性和方法由实际DVR继承
    
    属性:
      -  x：位置数组
      -  X：位置数组作为对角矩阵
      -  d1dvr：一阶导数DVR矩阵
      -  d2dvr：二阶导数DVR矩阵
      -  d1vbr：一阶导数VBR矩阵
      -  d2vbr：二阶导数VBR矩阵
      -  umat（）：DVR  - > VBR转换矩阵
      -  umatT（）：VBR  - > DVR转换矩阵
    
    方法（Dvr类和子类）：
      -  _calcumat（self，args）：计算变换矩阵，umat。Q = UXU+
      -  _cald1mat（self，args）：在DVR的基础上计算一阶导数的矩阵。D1_（dvr）= U + .D1_（vbr）.U
      -  _calcd2mat（self，args）：与_cald1mat相同
    """
    def __init__(self,npoints):
        self.npoints = npoints
        self.x = np.zeros(self.npoints,float)
        #self.X = None
        self.X = np.zeros((self.npoints,self.npoints),float)
        self.d1dvr = np.zeros((self.npoints,self.npoints),float)
        self.d1vbr = np.zeros((self.npoints,self.npoints),float)
        self.d2dvr = np.zeros((self.npoints,self.npoints),float)
        self.d2vbr = np.zeros((self.npoints,self.npoints),float)
        self.d4dvr = np.zeros((self.npoints,self.npoints),float)
        self.umat = np.zeros((self.npoints,self.npoints),float)
        self.weights = np.zeros(self.npoints,float)

    def delta_w(self):
        """ returns the inverse of the local DVR weight;
            for implementation of exact delta function """
        
        return 1.0/self.weights
        #array([inf, inf, inf])
        
#documentation P71
class Sindvr(Dvr):
    """
    sin DVR class
    """

    def __init__(self,npoints,qmin,qmax,**opts):
        Dvr.__init__(self,npoints,**opts)
        self.minmax = [qmin,qmax]                     #区间[x0,xn]
        self.l = self.minmax[1] - self.minmax[0]      #区间长度l
        self.Nf = float(self.npoints)                 #区间格点数n
        self.L = self.l * (self.Nf+1.0)/(self.Nf-1.0) #L=((n+1)/(n-1))*(xn-x0)
        self.dx = self.l/(self.Nf-1)                  #格点长度dx
        self.xo = self.minmax[0] - self.dx            #x0=x1-L/(n+1)
        self.weights = self.dx                        #w=L/(n+1)=((n+1)/(n-1))*(xn-x0)/(n+1)=(xn-x0)/(n-1)=dx

        self._calcumat()                              #在类中直接执行方法
        self._calcd1dvr()
        self._calcd2dvr()
        self._calcweights()
        
    def dkron(self,i,j):
        if isinstance(i,int) and isinstance(j,int):
            if i==j:
                return 1
            else:
                return 0
        else:
            return 0
        
    def sortEigValsVecs(self,evals,evecs):
        """
        Return the pair (eigvals,eigvecs) sorted from low to high values
        """
        evecs[:] = evecs.take(evals.argsort(),axis=1) # axis=1: eigenvectors are expected in columns
        evals[:] = evals.take(evals.argsort())
        
    def _calcweights(self):
        self.weights = np.zeros(self.npoints,float)
        self.weights[:] = self.dx                    #设置格点权重

    def _calcumat(self):
        # Generate the matrix Q(z), where z=cos(pi(x-xo)/L)
        Qz = np.zeros((self.npoints,self.npoints),float)
        for i in range(self.npoints):
            for j in range(self.npoints):
                Qz[i,j] = 0.5 * float( self.dkron(i,j+1) + self.dkron(i,j-1) )
        evals,evecs = np.linalg.eig(Qz)
        #print(evals)
        evals[:] = (self.L/math.pi)*np.arccos(evals)
        self.sortEigValsVecs(evals,evecs)
        self.umat = evecs
        # change sign of eigenvectors to ensure positive weights
        for i in range(self.npoints):
            if self.umat[0,i] < 0.0:
                self.umat[:,i] = -1.0*self.umat[:,i]
        self.umatT = np.transpose(self.umat)
        self.x = evals + self.xo
        #print(self.x)
        for i in range(self.npoints):
            self.X[i,i] = self.x[i]
        #self.X[i,i] = np.diag(self.x[i])

    def _calcd1dvr(self):
        for j in range(self.npoints):
            j1=j+1
            for k in range(self.npoints):
                k1=k+1
                if j != k:
                    self.d1vbr[j,k] = float((j1-k1)%2) * (4.0/self.L) * j1*k1/(j1**2 - k1**2)
        self.d1dvr = np.dot(self.umatT,np.dot(self.d1vbr,self.umat))

    def _calcd2dvr(self):
        for j in range(self.npoints):
            j1=j+1
            for k in range(self.npoints):
                k1=k+1
                self.d2vbr[j,k] = - self.dkron(j1,k1) * (j1*math.pi/self.L)**2
        self.d2dvr = np.dot(self.umatT,np.dot(self.d2vbr,self.umat))

#test
if __name__=='__main__':
    x_min = -5.0
    x_max =  5.0
    npx   =  10
    Qx_dvr = Sindvr(npx,x_min,x_max)
    print(Qx_dvr.x)
    #Qx_dvr_1=Qx_dvr._calcumat()
    print(Qx_dvr.x)
