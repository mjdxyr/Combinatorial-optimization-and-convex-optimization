import numpy as np
from sympy import beta


def f(x):
    return x[0][0]**2+x[1][0]**2-x[0][0]*x[1][0]-10*x[0][0]-4*x[1][0]+60


def fGradient(x):
    i = 2 * x[0][0]-x[1][0]-10
    j = 2 * x[1][0]-x[0][0]-4
    return np.array([[i, j]]).T


def phiDerivative1(x, lambdda, d):
    return 2*x[0][0]*d[0][0]+2*lambdda*d[0][0]**2+2*x[1][0]*d[1][0]+2*lambdda*d[1][0]**2-(x[0][0]*d[1][0]+x[1][0]*d[0][0]+2*lambdda*d[0][0]*d[1][0])-10*d[0][0]-4*d[1][0]


def phiDerivative2(x, lambdda, d):
    return 2*(d[0][0]**2+d[1][0]**2-d[0][0]*d[1][0])


def newton(phiDerivative1, phiDerivative2, x, d, epsilon1=0.0001, epsilon2=0.00001):
    lambdda = 1
    preLambdda = 0

    while True:
        if (abs(phiDerivative1(x, lambdda, d)) < epsilon1):
            return lambdda
        if (phiDerivative2(x, lambdda, d) <= 0):
            exit("失效")
        preLambdda = lambdda
        lambdda = lambdda - \
            phiDerivative1(x, lambdda, d) / phiDerivative2(x, lambdda, d)
        if (abs(lambdda - preLambdda) < epsilon2):
            return lambdda


def fr(f, initialX, H, epsilon=0.000001):
    x = initialX

    d = -fGradient(x)
    while True:
        if(np.linalg.norm(fGradient(x)) < epsilon):
            return x
        lambdda = newton(phiDerivative1, phiDerivative2, x, d)
        preX = x
        x = x+lambdda*d
        beta = np.dot(fGradient(x).T, fGradient(x)) / \
            np.dot(fGradient(preX).T, fGradient(preX))
        d = -fGradient(x)+beta*d


print(fr(f, np.array([[0, 0]]).T, np.eye(2)))
