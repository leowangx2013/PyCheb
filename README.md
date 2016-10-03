# PyCheb

**This is a Python package for solving ODEs using spectral methods**

## Background

**Differential equations** are used to describe phenomena of states and processes. Solutions of these problems explain patterns of them, thus people are eager to seek the solutions of these equations for describing states and making predictions of the future. An **ordinary differential equation (ODE)** is a differential equation containing one function (as the variable of the equation) of one independent variable (of the function) and its derivatives. Solving ODEs is comparatively easy but useful for scientists and engineers. And this is why we are interested in it and make such a Python package to solve it.


## Spectral Methods
Spectral methods are a class of techniques used in applied mathematics to solve differential equations numerically. The idea is to write the solution of the differential equation as a sum of certain "base function" (for example, as a Fourier series which is a sum of sinusoids) and then to choose the coefficients in the sum in order to satisfy the differential equation at any given accuracy. Spectral methods can be used to solve ordinary differential equations (ODEs), partial differential equations (PDEs) and eigenvalue problems involving differential equations. 
Compared to the traditional ODE solving methods, spectral methods are naturally with the advantage of super-fast convergence rate when the target function is smooth enough. As for more details for spectral methods, please check [the offical website of Chebfun](http://www.chebfun.org/publications/). It lists a bibliography for understanding spectral methods and the MATLAB project **Chebfun** which we will specifically talk about later.
## Related WorksAn object-oriented MATLAB system named [Chebfun](http://www.chebfun.org/) was created by a group of developers leading by Prof. Trefethen in University of Oxford in 2002. This system allows users to naturally input the ordinary differential equations in MATLAB code and get the solutions by using spectral methods. Chebfun has updated five versions from 2002 to 2012. And in 2013, it extended to Chebfun2 which moves the computations to multiple dimensions. Now, Chebfun has grown into a large open-source project on Github (<https://github.com/chebfun/chebfun>) and it has many related projects. 
[Pychebfun](https://github.com/olivierverdier/pychebfun) is a partial implementation of Chebfun in Python. This Python package incorporates functions like Chebyshev polynomial expansions, Lagrange interpolation, Clenshaw-Curtis quadrature and so on. However, it is unable to solve ODEs. 
Another Python partial implementation is [chebpy](https://github.com/chebpy/chebpy). It develops fast and is more powerful. 

## About PyChebOur project PyCheb is based on the concepts of Chebfun and its related projects. The reason why we implement spectral collocation method in Python is that first, Python is a high-level programming language and has a maturing ecosystem of scientific libraries, so it powerful enough for the implementation of this ordinary differential equation solving algorithm. Besides, Python is a general-purpose language, so it is not only used in academic settings but also in industry. Therefore, compared to the MATLAB implementation, it can be more widely used, by both scientists and engineers. Also, we wish that our implementation could be simpler and more easily to use. Furthermore, we expect to accelerate our system in the future by taking the advantages of GPU computing. This is important as we plan to expand our algorithm to solving partial differential equations, which is more complicated in computing. Python modules like [PyCUDA](https://mathema.tician.de/software/pycuda) and [Theano](http://deeplearning.net/software/theano) can help us achieve our goal. Thus, the Python implementation has a high extensibility. 
PyCheb can automatically solve the ordinary differential equations using spectral collocation method and plot the graph of the high-precision approximate solution. 

To use PyCheb, first download the PyCheb.py and import it into your project: 

	>> import PyCheb

For example, we want to solve the linear ODE: 	$$u''=e^{2x}, -1<x<1,u(\pm1)=0$$

Use 20 cheb points to solve this ODE:	>> f = PyCheb.ODEsolver("diff(u,2)", "math.exp(x)", [-1,1], 0, 0, N=20)
	>> f.solve(isLinear=1)
This will plot the following fig using the cheb points and the corresponding values produced by barycentric interpolation:
![fig1](https://raw.githubusercontent.com/leowangx2013/PyCheb/master/img/fig1.png)
A [scipy.interpolate.BarycentricInterpolator](http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.BarycentricInterpolator.html) instance will be returned, which is the barycentric interpolation result as well as the approximate solution of this ODE.
PyCheb also can automatically adapt the number of cheb points according to a given precision:

	>> f = PyCheb.ODEsolver("diff(u,2)", "math.exp(x)", [-1,1], 0, 0, precision=1e-15)
	>> f.solve(isLinear=1)PyCheb will recursively try to approximate the precision by increasing the number of cheb points in a Fibonacci sequence manner. In this example, 34 cheb points are used to approximate the solution.

As for nonlinear ODEs, we set the **isLinear** parameter for **solve()** method to 0. For example, we want to solve the nonlinear ODE:
	$$u''=e^{2u}, -1<x<1,u(\pm1)=0$$

We use the following code:
	>> f = PyCheb.ODEsolver("diff(u,2)", "math.exp(2*u)", [-1,1], 0, 0,precision=1e-15)
	>> f.solve(0)
We will get this fig:
![fig2](https://raw.githubusercontent.com/leowangx2013/PyCheb/master/img/fig2.png)
A [scipy.interpolate.BarycentricInterpolator](http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.BarycentricInterpolator.html) instance will be returned just like solving a linear ODE. The instance is the approximate solution of this nonlinear ODE.
