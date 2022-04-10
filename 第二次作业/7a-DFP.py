import numpy as np


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


def newton(phiDerivative1, phiDerivative2, x, d, epsilon1=0.0000000001, epsilon2=0.00000000001):
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


def dfp(f, initialX, H, epsilon=0.00000000001):
    x = initialX

    while True:
        fGradient1 = fGradient(x)
        if (np.linalg.norm(fGradient1) < epsilon):
            return x
        d = -np.dot(H, fGradient1)
        lambdda = newton(phiDerivative1, phiDerivative2, x, d)
        s = lambdda * d

        x = x + s
        y = fGradient(x) - fGradient1
        H = H + np.dot(s, s.T) / np.dot(s.T, y) - np.dot(H,
                                                         np.dot(y, np.dot(y.T, H))) / np.dot(y.T, np.dot(H, y))


print(dfp(f, np.array([[0, 0]]).T, np.eye(2)))
