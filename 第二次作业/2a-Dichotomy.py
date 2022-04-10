def phi(x):
    return 2*x*x-x-1


def phiDerivative(x):
    return 4*x-1


def dichotomy(phiDerivative, interval, epsilon):
    alpha = interval[0]
    beta = interval[1]

    while((beta-alpha) > epsilon):
        lambdda = (alpha+beta)/2
        if(phiDerivative(lambdda) > 0):
            beta = lambdda
        elif(phiDerivative(lambdda) < 0):
            alpha = lambdda
        else:
            return lambdda
    return (beta+alpha)/2


print(dichotomy(phiDerivative, [-1, 1], 0.06))
