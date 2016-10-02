import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib
import pprint
import sympy
import scipy.interpolate 
from numpy.linalg import inv
from numpy import linalg as LA
import math
import theano 
import theano.tensor as T
from theano import pp
import itertools
import re

def findStrAfterStr(myString, searchText, afterIndex):
    return myString.find(searchText, afterIndex)

def replaceAtIndex(oldString, newString, startIndex, endIndex):
	return oldString[:startIndex] + newString + oldString[endIndex+1:]

def polyFit(xs, ys):
	degree = len(xs)
	coeff = np.polyfit(xs, ys, degree)

	# Construct a complete polynomial
	poly = ""
	degree = len(coeff) - 1
	for i in range(0, len(coeff)):
		poly = poly + str(coeff[i]) + " * x**" + str(degree)
		if i != len(coeff) - 1:
			poly = poly + " + "
		degree = degree - 1
	print 'polyval = ', np.polyval(coeff, 10)
	return poly

def plot(xs, ys, title):

	plt.plot(xs, ys, '-')
	plt.title(title)
	plt.show()

def barycentricInterpolate(xs, ys):
	k = len(xs)
	# Make ws
	ws = []
	for i in range(0, k):
		w = 1
		for j in range(0, k):
			for l in range(0, k):
				if l != j:
					w = w * (xs[j] - xs[l])

		w = w ** (-1)
		ws.append(w)
	# Compute l(x)
	l = ''
	upper_part = '('
	for i in range(0, k):
		upper_part = upper_part + str(ws[i]*ys[i]) + '/(x - ' + str(xs[i]) + ')'
		if i != k-1:
			upper_part = upper_part + ' + '
	upper_part = upper_part + (')')

	lower_part = '('
	for i in range(0, k):
		lower_part = lower_part + str(ws[i]) + '/(x - ' + str(xs[i]) + ')'
		if i != k-1:
			lower_part = lower_part + ' + '
	lower_part = lower_part + ')'

	# print 'upper_part:', upper_part
	# print 'lower_part:', lower_part

	l = upper_part + " / " + lower_part
	x = 10
	print 'in barycentric interpolation: ', eval(l)
	file = open('/Users/tianshiwang/Desktop/f.txt', 'w')
	file.write(l)

	return l

class ODEsolver(object):
	lfunc = ""
	rfunc = ""
	domain = []
	lvalue = 0
	rvalue = 0
	N = 0

	def __init__(self, lfunc, rfunc, domain, lvalue, rvalue, precision = 0, N = 8):
		"""
        Init an Polyfun object from the given six parameters .
        values: Interpolation values
        vscale: The actual vscale; computed automatically if not given
        """

		self.lfunc = lfunc
		self.rfunc = rfunc
		self.domain = domain
		self.lvalue = lvalue
		self.rvalue = rvalue
		self.N = N
		self.precision = precision

		###########################

		# gy = T.grad(y, x)
		# f = theano.function([x], gy)
		# f = theano.function([x], gy)
		# print type(y)


		# fu = self.interpolatorToPoly(fu, domain, N)

		'''
			For auto adaption N
		'''
		# Replace "diff(...)"s with D**(n...n)
		'''
		old_lfunc = lfunc
		for token in tokens:
			m = re.findall(r'\d+', lfunc[token[0]: token[1]])
			#print "m[0]:", m[0]

			if (len(m) <= 1):
				for i in m[0]
					u_sym = T.dscalar('u')
					u_function = eval(fu)
					gy = T.grad(u_function, u_sym)
					temp = gy
					for j in range(0, m[0]-1)
						gy = T.grad(temp, x)
						temp = gy
		'''

	def getChebPoints(self):
		domain = self.domain
		N = self.N
		nodes = []
		center = (domain[0] + domain[1]) / 2.0
		multiple = abs(domain[1] - domain[0]) / 2.0

		for i in range(0, N+1):
			p = np.cos(i*np.pi/N)
			nodes.append(p * multiple + center)
		return nodes

	def getChebMatrix(self):
		N = self.N
		if N == 0:
			D = 0
			x = 1
			return [D, x]
		x = []
		for i in range(0, N+1):
			x.append(np.cos(np.pi * i/N))
		x = np.matrix(x).T

		c = []
		for i in range(0, N+1):
			if i == 0:
				c.append(2*(-1)**(i))
			elif i == N:
				c.append(2*(-1)**(i))
			else:
				c.append((-1)**(i))
		c = np.matrix(c).T
		X = np.matlib.repmat(x, 1, N+1)
		dX = X - np.matrix(X).T

		D = (c*(1.0/c).T) / (dX+(np.identity(N+1)))
		D = D - np.diag(np.squeeze(np.asarray(D.T.sum(axis=0))))

		return D

	'''
		def findStrAfterStr(myString, searchText, afterIndex):
			return myString.find(searchText, afterIndex)

		def replaceAtIndex(oldString, newString, startIndex, endIndex):
			return oldString[:startIndex] + newString + oldString[endIndex+1:]
	'''

	def getL(self, D):
		lfunc = self.lfunc

		N = len(D)
		tokens = []
		L = str(lfunc)
		start = 0
		end = 0

		# Replace 'diff(...)' with D**n
		while True:
			token = []
			start = findStrAfterStr(L, 'diff(', 0)
			end = findStrAfterStr(L, ')', start)

			if (start != -1):
				num = re.findall(r'\d+', L[start: end])
				if num != []:
					L = replaceAtIndex(L, 'D**' + str(num[0]), start, end)
				else:
					L = replaceAtIndex(L, 'D', start, end)
			else:
				break
		# Replace 'u' with eye(n)
		while True:
			start = findStrAfterStr(L, 'u', 0)
			if start != -1:
				L = replaceAtIndex(L, 'np.eye(' + str(N) + ')', start, start + 1)
			else:
				break

		L = eval(L)
		return L

	# us = self.solveLinearODE(ys, domain[0], lvalue, domain[1], rvalue, D, N)
	def solve(self, isLinear):

		domain = self.domain
		precision = self.precision
		rfunc = self.rfunc
		lfunc = self.lfunc
		rvalue = self.rvalue
		lvalue = self.lvalue

		p = 100
		pre_N = 5

		poly = self.testPrecision(isLinear)
		old_poly = poly
		if precision != 0:

			while(p > precision):

				temp = self.increaseN(pre_N, self.N)
				pre_N = self.N
				self.N = temp
				print 'pre_N = ', pre_N
				print 'N = ', self.N
				poly = self.testPrecision(isLinear)
				p = self.getMaxDifference(old_poly, poly, self.getChebPoints())
				print 'p = ', p
				print 'precision = ', precision
				self.N = temp

				print 'N = ', self.N
				old_poly = poly


		title = lfunc + ' = ' + rfunc + ', (' + str(domain[0]) + ', ' \
				+ str(domain[1]) + '), '  + 'u(' + str(domain[0]) + \
				') = ' + str(lvalue) + ', u(' + str(domain[1]) + \
				') = ' + str(rvalue)

		xs = self.getChebPoints()

		us = poly(xs)
		plot(xs, us, title)

		# pprint.pprint(np.polyval(fu, range(-1, 1, 0.01)))
		return poly

	def testPrecision(self, isLinear):
		domain = self.domain
		precision = self.precision
		rfunc = self.rfunc
		lfunc = self.lfunc
		rvalue = self.rvalue
		lvalue = self.lvalue
		N = self.N

		xs = self.getChebPoints()

		ys = []

		for x in xs:
			ys.append(eval(rfunc))
		print("ys=", ys)
		f = scipy.interpolate.BarycentricInterpolator(xs, ys)

		xnew = np.arange(domain[0], domain[1], 0.005)
		ynew = f(xnew)

		# cheb_series = np.polynomial.chebyshev.chebfit(x, y, degree)
		D = self.getChebMatrix()

		L = self.getL(D)

		# DN = eval(lfunc)
		# print "DN:", DN

		# uD = inv(Ld) * np.asarray(y).reshape(N+1,1)
		# Test solveLinearODE()
		# print self.solveLinearODE(y, -1, 0, 1, 0, D, N)
		# y, lbound, lvalue, rbound, rvalue, D, N

		# def solveNonlinearODE(self, L, y, lbound, lvalue, rbound, rvalue, N):

		if isLinear:
			us = solveLinearODE(L, ys, domain[0], lvalue, domain[1], rvalue, N)
		else:
			us = solveNonlinearODE(L, rfunc, domain[0], lvalue, domain[1], rvalue, N)

		# Test my barycentric interpolation
		# fu = barycentricInterpolate(xs, us)
		fu = scipy.interpolate.BarycentricInterpolator(xs, us)
		return fu

	def getMaxDifference(self, poly1, poly2, xs):
		us1 = poly1(xs)
		us2 = poly2(xs)
		print 'us1:'
		print us1
		print 'us2'
		print us2
		print '//////////'

		max_diff = 0
		for i in range(0, len(us1)):
			diff = us1[i] - us2[i]
			if diff > max_diff:
				max_diff = diff
				print 'max_diff:', max_diff
		return max_diff


	def increaseN(self, n1, n2):
		return n1 + n2

def solveLinearODE(L, y, lbound, lvalue, rbound, rvalue, N):

	# D2 = np.asarray(D * D)
	# D2 = D2[1:N, 1:N]
	L = L[1:N, 1:N]
	y = y[1:N]

	u = np.dot(inv(L), np.asarray(y).reshape(N-1,1)).tolist()
	u = np.insert(u, 0, lvalue)
	u = np.append(u, [rvalue])
	return u

# fix bug
def solveNonlinearODE(L, rfunc, lbound, lvalue, rbound, rvalue, N):
	# D2 = np.asarray(D * D)
	# D2 = D2[1:N, 1:N]

	L = L[1:N, 1:N]
	print("L=", L)
	u = []
	x = 0
	for i in range(0, N-1):
		u.append(0)

	# u = np.asarray(u_values).reshape(N-1, 1)

	change = 1
	time = 0
	while change > 1e-15 and time < 50:
		rfunc_values = []
		for i in u:
			x = i
			rfunc_values.append(eval(rfunc))
			print("rfunc_values=", rfunc_values)
		unew = np.dot(inv(L), np.asarray(rfunc_values).reshape(N-1, 1))
		change = LA.norm(unew - u , np.inf)
		u = unew
		time = time + 1
	u = np.squeeze(np.asarray(u))
	u = np.insert(u, 0, lvalue)
	u = np.append(u, rvalue)
	print("u: ", u)
	return u

def calDifference(poly1, poly2,):
	pass

def interpolatorToPoly(interpolator, domain, N):
	xnew = np.arange(domain[0], domain[1], 0.005)
	ynew = interpolator(xnew)

	return np.polyfit(xnew, ynew, N)


