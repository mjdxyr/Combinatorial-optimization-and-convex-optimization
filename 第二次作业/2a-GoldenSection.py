from math import sqrt


def phi(x):
    return 2*x*x-x-1


def glodenSection(phi, interval, epsilon):
    alpha = interval[0]
    beta = interval[1]
    gamma = (sqrt(5)-1)/2

    while True:
        lambdda = alpha+(1-gamma)*(beta-alpha)
        mu = alpha+gamma*(beta-alpha)

        while True:
            if(beta-alpha < epsilon):
                solve = (beta+alpha)/2
                return solve
            else:
                compare = phi(lambdda)-phi(mu)
                if(compare > 0):
                    alpha = lambdda
                    lambdda = mu
                    mu = alpha+gamma*(beta-alpha)
                elif(compare < 0):
                    beta = mu
                    mu = lambdda
                    lambdda = alpha+(1-gamma)*(beta-alpha)
                else:
                    alpha = lambdda
                    beta = mu
                    break


print(glodenSection(phi, [-1, 1], 0.06))
