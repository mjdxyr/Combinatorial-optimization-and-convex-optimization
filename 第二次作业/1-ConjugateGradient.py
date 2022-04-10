import numpy as np


def solveLambda(f, fGradient, x, d):
    # Wolfe-Powell法求解lambda值
    rho = 1/10
    sigma = 9/10
    alpha = 3/2
    beta = 1/2
    lambdda = 1
    phi1 = f(x)
    phi1Derivative = np.dot(fGradient(x).T, d)
    while True:
        phi2 = f(x+lambdda*d)
        if(phi2 <= phi1+rho*phi1Derivative*lambdda):
            phi2Derivative = np.dot(fGradient(x+lambdda*d).T, d)
            if(phi2Derivative >= sigma*phi1Derivative*lambdda):
                return lambdda
            else:
                lambdda = alpha*lambdda
        else:
            lambdda = beta*lambdda


def f(x):
    return x[0][0]-x[1][0]+2*x[0][0]**2+2*x[0][0]*x[1][0]+x[1][0]**2


def fGradient(x):
    x1 = 1+4*x[0][0]+2*x[1][0]
    x2 = -1+2*x[0][0]+2*x[1][0]
    return np.array([[x1, x2]]).T


def conjugateGradient(f, fGradient, initialPoint, epsilon):
    x = initialPoint
    preX = 0

    d = -fGradient(initialPoint)

    while True:
        if np.linalg.norm(fGradient(x)) < epsilon:
            return x
        else:
            lambdda = solveLambda(f, fGradient, x, d)
            preX = x
            x = x+lambdda*d
            beta = np.dot(fGradient(x).T, fGradient(x)-fGradient(preX)
                          )/np.dot(fGradient(preX).T, fGradient(preX))
            d = -fGradient(x)+beta*d


print(conjugateGradient(f, fGradient, np.array([[0, 0]]).T, 10**-6))
