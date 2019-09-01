from pylab import *
from utils import *
import constant as C
from tree import tree
# tree:
#                        (+)
#                        ||| n
#                         O -m
#                       /   \
#                      /     \
#                     O -mx   O -my
#                     |       |
#                     O -npx  O -npy
a=-0.5
b=1
c=1

n=C.n
m=C.m
mp=[C.mx,C.my]
np=[C.npx,C.npy]
x=arange(-5,5,0.1)
y=arange(-5,5,0.1)
dvr=[x,y]
sdof=2

def Trimat(md,a=a):
    '三对角矩阵展开在底层'
    temp=zeros((md,md),C.type_complex)
    for i in range(1,md-1):
        temp[i,i]=-2*a
        temp[i,i+1]=temp[i,i-1]=a
    temp[0,0]=temp[md-1,md-1]=-2*a
    temp[0,1]=temp[md-1,md-2]=a
    return temp

def diagmat(x,c,b=b):
    '''对角矩阵与位置c次方关系：维度、指数、系数
        diag([b*x[i]**c for i in range(md)])'''
    temp=diag([b*i**c for i in x])
    return temp

def h_operator(tv,dof,psi):
    '单体算符在所有node上形式'
    temp={i:{} for i in range(3)}
    for i in range(sdof):
        if i==dof:
            if tv==0:
                temp[0][i]=Trimat(np[dof])
            else:
                temp[0][i]=diagmat(dvr[dof],2)
            temp[1][i]=Rt(psi[i],temp[0][i])
        else:
            temp[0][i]=eye(np[i])
            temp[1][i]=eye(mp[i])
    temp1=temp[1][0]
    for i in range(1,sdof):
        temp1=kron(temp1,temp[1][i])
    temp[2][0]=Rt(psi[2],temp1)
    return temp

def W_operator(dof1,dof2,psi):
    '两体算符在所有node上的形式'
    temp={i:{} for i in range(3)}
    for i in range(sdof):
        if i==dof1:
            temp[0][i]=diagmat(dvr[dof1],1)
            temp[1][i]=Rt(psi[i],temp[0][i])
        else:
            temp[0][i]=eye(np[i])
            temp[1][i]=eye(mp[i])
        if i==dof2:
            temp[0][sdof+i]=diagmat(dvr[dof2],1)
            temp[1][sdof+i]=Rt(psi[i],temp[0][sdof+i])
        else:
            temp[0][sdof+i]=eye(np[i])
            temp[1][sdof+i]=eye(mp[i])
    temp1=temp[1][0]
    temp2=temp[1][sdof]
    for i in range(1,sdof):
        temp1=kron(temp1,temp[1][i])
        temp2=kron(temp2,temp[1][i])
    temp[2][0]=Rt(psi[2],temp1)
    temp[2][1]=Rt(psi[2],temp2)
    return temp

def Rt(A,H):
    '<A|H|A>'
    temp=matmul(conj(A),H)
    temp=matmul(temp,A.T)
    return temp

def P_operator(psi,i=2):
    '投影算符在sub_particle上'
    temp=matmul(psi[i].T,conj(psi[i]))
    return temp

def P_sub_operator(psi):
    '当定义单粒子函数为一个自由度的函数时投影算符在底层上的形式'
    temp={}
    for i in range(sdof):
        temp[i]=P_operator(psi,i)
    return temp

def n_base(psi,a):
    '求出波函数任意一个number state分量并展开在截断层基矢上'
    n=sum(a) #粒子数
    m=len(a) #轨道数
    temp=zeros((n,prod(mp)),dtype=C.type_complex)
    index=ns_index_(a)
    num_2=prod(factorial(a))
    num_2=sqrt(num_2/factorial(n))
    for i in range(m):
        for j in range(a[i]):
            for k in range(prod(mp)):
                temp[sum(a[0:i])+j,k]=psi[2][i,k]
    temp3=temp[0,:]
    for i in range(1,n):
        temp3=kron(temp3,temp[i])
    temp3=psi[3][index]*num_2*temp3
    return temp3

def PSI(psi):
    '总的波函数在截断层的展开'
    temp=0
    for i in ns_distrubution_(n,m):
        temp=temp+n_base(psi,i)
    return temp

def shf_sub(s,dof,PSI):
    '当定义单粒子函数为一个自由度的函数时single hole function在截断层形式'
    if dof==0:
        d1=prod(mp)**s
    else:
        d1=prod(mp)**s*prod(mp[0:dof])
    if dof==sdof-1:
        d2=prod(mp)**(n-s-1)
    else:
        d2=prod(mp)**(n-s-1)*prod(mp[dof+1:sdof])
    temp=zeros((mp[dof],d1*d2),dtype=C.type_complex)
    for i in range(mp[dof]):
        for j in range(d1):
            temp[i,j*d2:j*d2+d2]=PSI[(j*mp[dof]+i)*d2:(j*mp[dof]+i+1)*d2]
    return temp

def pho_sub(shf):
    '当定义单粒子函数为一个自由度的函数时密度矩阵形式'
    pass

def ns_distrubution_(n,m):
    '存储N个粒子m个态体系的二次量子化波矢'
    return ns_distribution(n,m)

def ns_index_(a):
    '算一个二次量子化的波矢对应的psi系数的位置'
    return ns_index(a)

def shf_operator(psi):
    'single hole fuction在占据数表象的形式(1.15)'
    ns=ns_distrubution_(n-1,m)
    temp=zeros(ns.shape,C.type_complex).T
    for i,k in enumerate(ns):
        for j in range(m):
            temp2=k.copy()
            temp2[j]=temp2[j]+1
            temp[j][i]=sqrt((k[j]+1)/n)*psi[3][ns_index_(temp2)]
    return temp

def pho_operator(shf):
    '约化密度矩阵在二次量子化的形式1.22'
    temp=matmul(conj(shf),shf.T)
    return temp

def h1_meanoperator(htable):
    #相互作用前一项在展开在subspf上
    l=sdof**2*num_combination(n,2)
    temp={}
    for i in range(l):
        temp1=htable[i][1][0]
        for j in range(1,sdof):
            temp1=kron(temp1,htable[2*sdof*n+i][1][j])
        temp[i]=temp1
    return temp

def h2_meanoperator(htable,top_psi):
    '(1.23)'
    ns=ns_distrubution_(n-2,m)
    l=sdof**2*num_combination(n,2)
    def P(_ns,i,j):
        'ns[(i,j)]+=1\nP(i,j)=sqrt(ni(nj-delta(i,j)))'
        ni=_ns[i]+1
        nj=_ns[p]+1
        P=sqrt(ni*(nj - (1 if i==j else 0)))
        return P
    def A_index(_ns,i,j):
        ns=_ns.copy()
        ns[[i,j]]+=1
        return ns_index_(ns)
    Hv=zeros([l,m,m],complex)
    for v in range(l):
        for i in range(m):
            for j in range(m):
                s=0
                for _i,_j in enumerate(ns):
                    for p in range(m):
                        for q in range(m):
                            P_ip=P(_j,i,p)
                            P_jq=P(_j,j,q)
                            A_ip=top_psi[A_index(_j,i,p)]
                            A_jq=top_psi[A_index(_j,j,q)]
                            p_v_q=htable[2*sdof*n+v][2][1][p,q]
                            s+=P_ip*P_jq*A_ip*A_jq*p_v_q
                Hv[v,i,j]=s
    return Hv

def drives(t,y):
    'psi[2] 1.21'
    l=sdof**2*num_combination(n,2)
    n_spfs=prod(mp)
    dy=zeros_like(y)
    def rh(_n,_i):
        s=0
        for v in range(l):
            a=matmul(eye(100)-P,H1[v])
            b=matmul(linalg.inv(pho),H2[v])
            for _j in range(n_spfs):
                for _m in range(m):
                    s1=a[_i,_j]
                    s2=b[_n,_m]
                    s+=s1*s2*y[_m,_j]
        return s        
    for _n in range(m):
        for _i in range(n_spfs):
            dy[_n,_i]=rh(_n,_i)             
    return dy

if __name__=='__main__':
    tree=tree()
    tree.init_from_psi('dvr/restart.ini')
    t=tree.layers
    psi=[t[3][0]._psi,t[3][1]._psi,t[2][0]._psi,reshape(t[1][0]._psi,(-1))+1e-10]
    l=2*sdof*n+sdof**2*num_combination(n,2)
    htable={}#{i:{j:{} for j in range(3)} for i in range(l)} #三个维度，第一维表示第几个算符，第二维表示层数，第三维表示同一层的第几个
    #动能项
    for i in range(n):
        for j in range(sdof):
            htable[i*sdof+j]=h_operator(0,j,psi)
    #势能项
    for i in range(n):
        for j in range(sdof):
            htable[sdof*n+i*sdof+j]=h_operator(1,j,psi)
    #相互作用能
    for s1 in range(n-1):
        for s2 in range(s1+1,n):
            for dof1 in range(sdof):
                for dof2 in range(sdof):
                    htable[2*sdof*n+(s1*n+s2-s1)*(dof1*sdof+dof2)]=W_operator(dof1,dof2,psi)
    #投影算符
    P=P_operator(psi)
    pho=pho_operator(shf_operator(psi))
    #平均场算符
    H1=h1_meanoperator(htable)  
    H2=h2_meanoperator(htable,psi[3])
    dy=drives(0,psi[2])
    PSI=PSI(psi)
    shf=shf_sub(0,1,PSI)
    
    
    
    

