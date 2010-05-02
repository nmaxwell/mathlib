
from math import *
from ode import *

foo = ode()
foo.solver = RK4solver()
foo.step_size = 0.01
foo.stop = 10.0
foo.rhs = lambda t,y: -3.0*y + 6.0*t + 5.0
solution = lambda t: 2.0*exp(-3.0*t) + 2.0*t + 1.0

xy = (0.0, 3.0)

while xy[0] < 10.0:
    foo.stop = xy[0] + 0.1
    xy = foo(xy)
    
    print xy, solution(xy[0]), foo.stop
    
    




"""
xy = foo(xy)

solution = lambda t: 2.0*exp(-3.0*t) + 2.0*t + 1.0

print xy, solution(xy[0])
"""

"""

Y = [ xy = foo(xy) k in range(100) ]





#x,y = foo.run()

"""






















