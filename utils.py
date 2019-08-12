'''
工具模块:包含数学函数和输入输出函数
'''
import numpy as np
from functools import wraps,reduce
from itertools import combinations_with_replacement

def map_list(iter_=False):
    '根据提供的函数对指定序列做映射\n仅用于单参数函数'
    def decorator(func):
        @wraps(func)
        def main(a):
            if hasattr(a,'__iter__'):
                r=map(func,a)
                if iter_:
                    return r
                else:
                    return list(r)
            else:
                return func(a)
        return main
    return decorator

@map_list()
def factorial(a):
    '返回一个整数或一个整数列表的阶乘'
    return np.math.factorial(a)

def num_permutation(n,m):
    'A(n,m) 从n个不同元素中, 任取m个不同的元素按照一定的顺序排成一列'
    assert m<=n,'组合参数不符合要求, %d>%d'%(m,n)
    return factorial(n)//factorial(n-m)

def num_combination(n,m):
    'C(n,m) n 元集合的 m 元集合数'
    assert m<=n,'组合参数不符合要求, %d>%d'%(m,n)
    if m==0 or m==n:
        return 1
    return num_permutation(n,m)//factorial(m)

def num_combination_repetition(n,m):
    'H(n,m) 从n个不同的元素取出m个元素的允许重复的组合总数'
    assert n>=1 ,'组合参数不符合要求, %d<1'%(n)
    return num_combination(n+m-1,m)

def ns_distribution(n,m):
    'n个粒子m个轨道的 number state in a bosonic node'
    #`n+1`个位置插`m-1`个板,板位置差为ns
    pos=combinations_with_replacement(range(n+1),r=m-1)
    pos=np.array(list(pos))
    c=pos[:,1:]-pos[:,:-1]
    t=n*np.ones((1,len(pos)),int)-pos[:,-1]
    ns=np.hstack((pos[:,0:1],c,t.T))
    return ns[::-1]

def ns_index(ns):
    '得到`ns`在一组number state中的序号'
    #先逆变换为插板位置,再变为序号
    n=sum(ns)
    m=len(ns)
    num_ns=num_combination(n+m-1,n)
    t=0
    index=0
    o=0
    for i in range(m-1):
        t=t+ns[i]
        for j in range(t-o):
            w=num_combination_repetition(n+1-j-o,m-2-i)
            index+=w
        o=t
    return num_ns-index-1
