'''
工具模块:包含数学函数和输入输出函数
'''
import numpy as np
from functools import wraps,reduce

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
