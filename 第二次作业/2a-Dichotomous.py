def phi(x):
    return 2*x*x-x-1


def dichotomous(phi, interval, epsilon):
    alpha = interval[0]
    beta = interval[1]

    while((beta-alpha) > 2 * epsilon*(1+epsilon/(beta-alpha))):
        lambdda = (alpha+beta)/2-epsilon
        mu = (alpha+beta)/2+epsilon
        if(phi(lambdda) < phi(mu)):
            beta = mu
        else:
            alpha = lambdda

    return (beta+alpha)/2


print(dichotomous(phi, [-1, 1], 0.06))
