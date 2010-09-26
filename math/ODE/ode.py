

"""
A general inderface for solving ODEs


"""


class ode:
    """
    Solving y'(t) = RHS(t,y); y(t0) = y0
    
    parameters:
        RHS the right hand side of the equation
        t0 the initial point (start)
        y0 the initial condition
    
    """
    
    def __init__(self, **kwargs):
        """
        parameters:
        [rhs, solver, step_control, output, terminator]
        [step_size, stop]
        """
        self.rhs = None
        self.solver = None
        self.step_control = None
        self.output = None
        self.terminator = None
        
        self.step_size = 1.0
        self.stop = 1.0
        
        for key in kwargs:
            setattr(self, key, kwargs[key])
    
    def run(self, initial_condition ):
        """
        Solve the ode.
        
        initial_condition is a pair, (start, initial_value)
        
        """
        
        if self.terminator is None:
            self.terminator = fixedTerminator(self )
        
        t = initial_condition[0]
        y=None
        try:
            y = initial_condition[1].copy()
        except:
            y = initial_condition[1]
        
        dt = self.step_size
        
        if self.output is not None:
            self.output.output(t, y, start=True, end=False)
        
        continue_flag=True
        
        while continue_flag:
            #print x,y, dx
            
            if self.step_control is not None:
                dt = self.step_control.get_step_size(t, y)
            else:
                t = round((t-initial_condition[0])/dt)*dt + initial_condition[0]
            
            dt, continue_flag, kill_flag = self.terminator.terminate(t, y, dt)
            
            if kill_flag:
                break
            
            if not continue_flag and self.step_control is not None:
                self.step_control.set_step_size(dt)
            
            t, y = self.solver.step(self.rhs, t, dt, y )
            
            if self.output is not None:
                self.output.output(t, y, start=False, end=not continue_flag)
        
        return t, y
    
    def solve(self, initial_condition):
        return self.run(initial_condition)
    
    def __call__(self, initial_condition):
        return self.run(initial_condition )




class fixedTerminator:
    
    def __init__(self, ode_instance ):
        self.ode_instance = ode_instance
    
    def terminate(self, t, y, dt):
        if t+dt > self.ode_instance.stop:
            dt = self.ode_instance.stop-t
            return dt, False, False
        return dt, True, False

class RK4solver:
    
    def step(self, rhs, t, dt, y):
        
        K1 = rhs(t, y)
        K2 = rhs(t + dt/2.0, y+K1*(dt/2.0))
        K3 = rhs(t + dt/2.0, y+K2*(dt/2.0))
        K4 = rhs(t+dt, y+K3*dt)
        
        y += (K1/6.0 + K2/3.0 + K3/3.0 + K4/6.0)*dt
        t += dt
        
        return t, y


class DP54(ode):
    
    """
    solver:
    def step(self, rhs, t, dt, y):
        return t, y
    
    step_control:
    def get_step_size(self, t, y):
        return dt
    def set_step_size(self, dt):
    
    output:
    sef output(self, t, y, start=True, end=False)
    
    terminator:
    def terminate(t, y, dt):
        return dt, continue_flag, kill_flag
    
    """
    
    def __init__(self, **kwargs):
        self.sample_points = None
        
        if self.output is not None:
            self.output.output(t, y, start=True, end=False) 
        
    """
    def step(self, rhs, t, dt, y):
        

    def get_step_size(self, t, y):
        
        self.dt = self.dt_next
        if self.dt < dt_min:
            self.dt = dt_min
        if self.dt > self.dt_max:
            self.dt = self.dt_max
        
        error = (self.eps+1)*2
        
        while error > self.eps:
            
            self.eval_rhs(t, y)
            
            temp  = self.K0*(35.0/384.0-5179.0/57600.0)
            temp += self.K2*(500.0/1113.0-7571.0/16695.0)
            temp += self.K3*(125.0/192.0-393.0/640.0)
            temp += self.K4*(-2187.0/6784.0+92097.0/339200.0)
            temp += self.K5*(11.0/84.0-187.0/2100.0)
            temp += self.K6*(-1.0/40.0)
            
            try:
                error = norm(temp)
            except:
                try:
                    error = temp.norm()
                except:
                    print("Error: norm not defined on vector type.")
            
            if error > self.eps:
                self.dt /= 2.0
        """
            
        



# Unit Tests:

import unittest
import math




class TestSequenceFunctions(unittest.TestCase):
    
    def test_RK4solver(self):
        foo = ode()
        foo.solver = RK4solver()
        foo.step_size = 0.01
        foo.stop = 10.0
        foo.rhs = lambda t,y: -3.0*y + 6.0*t + 5.0
        solution = lambda t: 2.0*math.exp(-3.0*t) + 2.0*t + 1.0
        
        xy = (0.0, 3.0)
        acc_err = 0
        acc_mag = 0
        
        while xy[0] < 10.0:
            foo.stop = xy[0] + 0.1
            xy = foo(xy)
            acc_err += abs(xy[1] - solution(xy[0]))
            acc_mag += abs(solution(xy[0]))
        
        err = acc_err/acc_mag
        self.assertTrue(err < 1E-9)

        




if __name__ == '__main__':
    unittest.main()



