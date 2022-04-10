import numpy as np


def f(x):
    return 100*(x[1][0]-x[0][0]**2)**2+(1-x[0][0])**2


def fGradient(x):
    i = 100*(-4*x[0][0]*x[1][0]+4*x[0][0]**3)+(-2+2*x[0][0])
    j = 100*(2*x[1][0]-2*x[0][0]**2)
    return np.array([[i, j]]).T


def wolfePowell(f, fGradient, x, d):
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


print(wolfePowell(f, fGradient, np.array([[-1, 1]]).T, np.array([[1, 1]]).T))
