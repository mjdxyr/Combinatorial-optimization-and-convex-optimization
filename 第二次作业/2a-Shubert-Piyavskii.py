def phi(x):
    return 2*x*x-x-1


def solveIntersection(phi, a, b, L):
    x = ((phi(a)-phi(b))/L+a+b)/2
    y = (L*a-L*b+phi(a)+phi(b))/2
    return [x, y]


def shubert(phi, interval, epsilon):
    alpha = interval[0]
    beta = interval[1]
    L = 100

    x = (alpha+beta)/2
    xL = [alpha, x, beta]
    xIntersection = [solveIntersection(phi, alpha, x, L)[0],
                     solveIntersection(phi, x, beta, L)[0]]
    yIntersection = [solveIntersection(phi, alpha, x, L)[1],
                     solveIntersection(phi, x, beta, L)[1]]

    while True:
        index = yIntersection.index(min(yIntersection))
        if(phi(xIntersection[index])-yIntersection[index] < epsilon):
            return xIntersection[index]
        xL.insert(index+1, xIntersection[index])
        x1 = solveIntersection(
            phi, xL[index], xL[index+1], L)[0]
        y1 = solveIntersection(
            phi, xL[index], xL[index+1], L)[1]
        x2 = solveIntersection(
            phi, xL[index+1], xL[index+2], L)[0]
        y2 = solveIntersection(
            phi, xL[index+1], xL[index+2], L)[1]
        xIntersection[index] = x1
        xIntersection.insert(index+1, x2)
        yIntersection[index] = y1
        yIntersection.insert(index+1, y2)


print(shubert(phi, [-1, 1], 0.06))
