import PyCheb
import numpy as np

# transform from u to eye()
'''
	polyfuum.Polyfun(lfunc, rfunc, domain, lvalue, rvalue, N)
'''
f = PyCheb.ODEsolver("diff(u,2)", "math.exp(x)", [-1,1], 0, 0, precision=1e-11)
f.solve(1)