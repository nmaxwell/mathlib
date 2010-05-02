

"""
A general inderface for solving ODEs


"""



class fixedTerminator:
    
    def __init__(self, ode_instance ):
        self.ode_instance = ode_instance
    
    def terminate(self, x, y, dx):
        if x+dx > self.ode_instance.stop:
            dx = self.ode_instance.stop-x
            return dx, False, False
        return dx, True, False


class ode:
    """
    Solving y'(x) = RHS(x,y); y(x0) = y0
    
    parameters:
        RHS the right hand side of the equation
        x0 the initial point
        y0 the initial condition
    
    """
    
    def __init__(self, ):
        self.rhs = None
        self.solver = None
        self.step_control = None
        self.output = None
        self.terminator = None
        
        self.rhs = None
        self.step_size = 1.0
        self.stop = 1.0
    
    
    def run(self, initial_condition ):
        
        if self.terminator is None:
            self.terminator = fixedTerminator(self )
        
        x = initial_condition[0]
        y = initial_condition[1]
        dx = self.step_size
        
        if self.output is not None:
            self.output.output(x, y, start=True, end=False)
        
        continue_flag=True
        
        while continue_flag:
            #print x,y, dx
            
            if self.step_control is not None:
                dx = self.step_control.get_step_size(x,y)
            else:
                x = round((x-initial_condition[0])/dx)*dx + initial_condition[0]
            
            dx, continue_flag, kill_flag = self.terminator.terminate(x, y, dx)
            
            if kill_flag:
                break
            
            if not continue_flag and self.step_control is not None:
                self.step_control.set_step_size(dx)
            
            x, y = self.solver.step(self.rhs, x, dx, y )
            
            if self.output is not None:
                self.output.output(x, y, start=False, end=not continue_flag)
        
        return x, y
    
    def __call__(self, initial_condition):
        return self.run(initial_condition )



class RK4solver:
    
    def step(self, rhs, x, dx, y):
        
        K1 = rhs(x, y)
        K2 = rhs(x + dx/2.0, y+K1*(dx/2.0))
        K3 = rhs(x + dx/2.0, y+K2*(dx/2.0))
        K4 = rhs(x+dx, y+K3*dx)
        
        y += (K1/6.0 + K2/3.0 + K3/3.0 + K4/6.0)*dx
        x += dx
        
        return x, y








