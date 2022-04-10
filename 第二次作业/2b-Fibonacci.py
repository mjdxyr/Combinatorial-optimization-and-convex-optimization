def fibonacci(n):
    if(n == 1 or n == 0):
        return 1
    else:
        return fibonacci(n-1)+fibonacci(n-2)


def phi(x):
    return 3*x*x-21.6*x-1


def fibonacciSearch(phi, interval, epsilon):
    alpha = interval[0]
    beta = interval[1]

    n = 2
    Fn = fibonacci(n)

    while(Fn < (beta-alpha)/epsilon):
        n += 1
        Fn = fibonacci(n)

    k = 1
    while True:
        lambdda = alpha+(fibonacci(n-k-1)/fibonacci(n-k+1))*(beta-alpha)
        mu = alpha+(fibonacci(n-k)/fibonacci(n-k+1))*(beta-alpha)

        while True:
            k += 1
            compare = phi(lambdda)-phi(mu)
            if(compare > 0):
                alpha = lambdda
                lambdda = mu
                mu = alpha+(fibonacci(n-k)/fibonacci(n-k+1))*(beta-alpha)
            elif(compare < 0):
                beta = mu
                mu = lambdda
                lambdda = alpha+(fibonacci(n-k-1) /
                                 fibonacci(n-k+1))*(beta-alpha)
            else:
                alpha = lambdda
                beta = mu
                break
            if(k == n-2):
                solve = (beta+alpha)/2
                return solve


print(fibonacciSearch(phi, [0, 25], 0.08))
