import numpy as np
import math


def f(x):
    return (x[0][0]+2*x[1][0]-7)**2+(2*x[0][0]+x[1][0]-5)**2


def fGradient(x):
    i = 10*x[0][0]+8*x[1][0]-34
    j = 8*x[0][0]+10*x[1][0]-38
    return np.array([[i, j]]).T


def sgd(fGradient, initialX=np.array([[0, 0]]).T, learningRate=0.1, epsilon=0.00000001):
    x = initialX
    while True:
        gradient = fGradient(x)
        if(np.linalg.norm(gradient) < epsilon):
            return x
        x = x-learningRate*gradient


def momentum(fGradient, initialX=np.array([[0, 0]]).T, learningRate=0.1, alpha=0.1, epsilon=0.0000000001):
    x = initialX
    v = np.array([[0, 0]]).T
    while True:
        gradient = fGradient(x)
        v = -learningRate*gradient+alpha*v
        x = x+v
        if(np.linalg.norm(v) < epsilon):
            return x


def nag(fGradient, initialX=np.array([[0, 0]]).T, learningRate=0.1, alpha=0.1, epsilon=0.000000001):
    x = initialX
    v = np.array([[0, 0]]).T
    while True:
        v = alpha*v-learningRate*fGradient(x+alpha*v)
        x = x+v
        if(np.linalg.norm(v) < epsilon):
            return x


def adagrad(fGradient, initialX=np.array([[0, 0]]).T, learningRate=0.1,  delta=0.00000001, epsilon=0.000000000001):
    x = initialX
    s = 0
    while True:
        gradient = fGradient(x)
        s = s+gradient*gradient
        v = -learningRate/np.sqrt(s+delta)*gradient
        x = x+v
        if(np.linalg.norm(v) < epsilon):
            return x


def rmsprop(fGradient, initialX=np.array([[0, 0]]).T, learningRate=0.1,  delta=0.00000001, gamma=0.9, epsilon=0.000001):
    x = initialX
    s = 0
    while True:
        gradient = fGradient(x)
        s = gamma*s+(1-gamma)*gradient*gradient
        v = -learningRate/np.sqrt(s+delta)*gradient
        x = x+v
        if(np.linalg.norm(v) < epsilon):
            return x


def adadelta(fGradient, initialX=np.array([[0, 0]]).T,  delta=0.00000001, rho=0.9, epsilon=0.0000001):
    x = initialX
    s = 0
    deltaX = 0
    while True:
        gradient = fGradient(x)
        s = rho*s+(1-rho)*gradient*gradient
        gradient1 = np.sqrt((deltaX+delta)/(s+delta))*gradient
        x = x-gradient1
        deltaX = rho*deltaX-(1-rho)*gradient1*gradient1

        if(np.linalg.norm(gradient) < epsilon):
            return x


def adam(fGradient, initialX=np.array([[0, 0]]).T, beta1=0.9, beta2=0.999, learningRate=0.1, delta=0.00000001,epsilon=0.0000000001):
    x = initialX
    v = 0
    s = 0
    t = 1
    while True:
        gradient = fGradient(x)
        v = beta1*v+(1-beta1)*gradient
        s = beta2*s+(1-beta2)*gradient*gradient
        v_corr = v/(1-beta1**t)
        s_corr = s/(1-beta2**t)
        gradient1 = learningRate*v_corr/(np.sqrt(s_corr)+delta)
        x = x-gradient1
        if(np.linalg.norm(gradient) < epsilon):
            return x


print(sgd(fGradient))
print(momentum(fGradient))
print(nag(fGradient))
print(adagrad(fGradient))
print(rmsprop(fGradient))
print(adadelta(fGradient))
print(adam(fGradient))
