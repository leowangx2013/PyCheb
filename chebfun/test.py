import PyCheb
import numpy as np

# Linear example
# f = PyCheb.ODEsolver("diff(u,2)", "math.exp(x)", [-1,1], 0, 0, precision=1e-15)
# Using N instead of precision
#f = PyCheb.ODEsolver("diff(u,2)", "math.exp(x)", [-1,1], 0, 0, N=5)
# f.solve(1)

# Nonlinear example
f = PyCheb.ODEsolver("diff(u,2)", "math.exp(2*u)", [-1,1], 0, 0,precision=1e-15)
f.solve(0)