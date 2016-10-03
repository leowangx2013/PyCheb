import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib
import pprint
import sympy
import scipy.interpolate 
from numpy.linalg import inv
from numpy import linalg as LA
import math
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

	l = upper_part + " / " + lower_part
	print 'in barycentric interpolation: ', eval(l)

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

	def getL(self, D):
		lfunc = self.lfunc

		N = len(D)
		L = str(lfunc)

		while True:
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
				poly = self.testPrecision(isLinear)
				p = self.getMaxDifference(old_poly, poly, self.getChebPoints())
				self.N = temp

				old_poly = poly
		print "N=", self.N

		title = lfunc + ' = ' + rfunc + ', (' + str(domain[0]) + ', ' \
				+ str(domain[1]) + '), '  + 'u(' + str(domain[0]) + \
				') = ' + str(lvalue) + ', u(' + str(domain[1]) + \
				') = ' + str(rvalue)

		xs = self.getChebPoints()

		us = poly(xs)
		plot(xs, us, title)

		return poly

	def testPrecision(self, isLinear):
		domain = self.domain
		rfunc = self.rfunc
		rvalue = self.rvalue
		lvalue = self.lvalue
		N = self.N

		xs = self.getChebPoints()

		ys = []



		D = self.getChebMatrix()

		L = self.getL(D)

		if isLinear:
			for x in xs:
				ys.append(eval(rfunc))
			print("ys=", ys)
			us = solveLinearODE(L, ys, domain[0], lvalue, domain[1], rvalue, N)
		else:
			us = solveNonlinearODE(L, rfunc, domain[0], lvalue, domain[1], rvalue, N)

		fu = scipy.interpolate.BarycentricInterpolator(xs, us)
		return fu

	def getMaxDifference(self, poly1, poly2, xs):
		us1 = poly1(xs)
		us2 = poly2(xs)
		print 'us1:'
		print us1
		print 'us2'
		print us2

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

	L = L[1:N, 1:N]
	y = y[1:N]

	u = np.dot(inv(L), np.asarray(y).reshape(N-1,1)).tolist()
	u = np.insert(u, 0, lvalue)
	u = np.append(u, [rvalue])
	return u

def solveNonlinearODE(L, rfunc, lbound, lvalue, rbound, rvalue, N):

	L = L[1:N, 1:N]
	print("L=", L)
	us = []
	x = 0
	for i in range(0, N-1):
		us.append(0)

	change = 1
	time = 0
	while change > 1e-15 and time < 50:
		rfunc_values = []
		for u in us:
			rfunc_values.append(eval(rfunc))
		unew = np.dot(inv(L), np.asarray(rfunc_values).reshape(N-1, 1))
		change = LA.norm(unew - us , np.inf)
		us = unew
		time = time + 1
	us = np.squeeze(np.asarray(us))
	us = np.insert(us, 0, lvalue)
	us = np.append(us, rvalue)
	print("us: ", us)
	return us

def calDifference(poly1, poly2,):
	pass

def interpolatorToPoly(interpolator, domain, N):
	xnew = np.arange(domain[0], domain[1], 0.005)
	ynew = interpolator(xnew)

	return np.polyfit(xnew, ynew, N)


