def phi(x):
    return 3*x*x-21.6*x-1


def phiDerivative(x):
    return 6*x-21.6


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


print(dichotomy(phiDerivative, [0, 25], 0.08))
