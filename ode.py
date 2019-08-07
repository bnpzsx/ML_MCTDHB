'''
solver类:用于求解常微分方程
'''
import numpy as np
from scipy.integrate import ode
#from scipy.integrate import complex_ode

__all__ = ['solver']
class solver():
    '''使用scipy库的积分工具求解微分方程
    Parameters
    ----------
    f : callable ``f(t, y, *f_args)``
        Right-hand side of the differential equation. t is a scalar,
            ``y.shape == (n,)``.
            ``f_args`` is set by calling ``set_f_params(*args)``.
            `f` should return a scalar, array or list (not a tuple).
    jac : callable ``jac(t, y, *jac_args)``, optional
        Jacobian of the right-hand side, ``jac[i,j] = d f[i] / d y[j]``.
        ``jac_args`` is set by calling ``set_jac_params(*args)``.  
    integrator : the integrator used to numeric integrators,
        Available integrators are listed below.
        "vode", "zvode", "lsoda", "dopri5", "dop853".
    stiff : weather the problem is a stiff problem,
        this determine the method to use, Adams (non-stiff) or BDF (stiff).
        
    other parameters :
      - atol : float or sequence absolute tolerance for solution
      - rtol : float or sequence relative tolerance for solution
      - first_step : float
      - min_step : float
      - max_step : float
      - safety : float
        Safety factor on new step selection (default 0.9)
      - ifactor : float
      - dfactor : float
        Maximum factor to increase/decrease step size by in one step
      - beta : float
        Beta parameter for stabilised step size control.
    '''

    def __init__(self,f,jac=None,integrator='zvode',stiff=False,**para):
        self._f=f
        self._jac=jac
        self.integrator=integrator
        self.stiff=stiff
        self.solver=ode(f,jac)
        method='adams' if not stiff else 'bdf'
        para['method']=method
        self.solver.set_integrator(integrator,**para)

    def set_initial_value(self, y0, t0):
        '设置初值'
        self.solver.set_initial_value(y0, t0)

    def set_params(self,f_args=[],jac_args=[]):
        '设置f, jac函数调用时需要的额外参数'
        self.solver.set_f_params(*f_args)
        self.solver.set_jac_params(*jac_args)

    def run(self,t,dt):
        '积分至t, 间隔dt返回结果\n返回列表 t, y'
        r_t=[]
        a_t=r_t.append
        r_y=[]
        a_y=r_y.append
        s=self.solver
        _integrator=s._integrator
        while _integrator.success == 1 and s.t < t:
            a_y(s.integrate(s.t+dt))
            a_t(s.t)
        return r_t, r_y

    def run_relax(self,t):
        '直接积分到t'
        return self.solver.integrate(t,relax=True)

#测试代码
if __name__=='__main__':
    y0, t0 = [3], 0
    def f(x,y,*arg):
        return -x*x*y*y
    def y(x):
        return 3/(1+x**3)
    dt = 0.1
    T = 1.5
    s=solver(f,integrator='dopri5',atol=1e-10,rtol=1e-10)
    s.set_initial_value(y0,t0)
    t,r=s.run(T,dt)
    y_=list(map(y,t))
    error=np.reshape(r,(-1,))-y_
    print(error)
