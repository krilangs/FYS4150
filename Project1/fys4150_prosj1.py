import numpy as np
import matplotlib.pyplot as plt
import time

def f(x):
    """Source term f(x)"""
    return 100.*np.exp(-10.*x)


def gausselim(n):
    a = np.ones(n-1)*(-1)
    b = np.ones(n)*2
    c = np.ones(n-1)*(-1)
    
    x = np.linspace(0, 1., n+2)
    h = 1./(n+1)
    
    #f = 100.*np.exp(-10.*x)*h**2
    b_tilde = np.zeros(n)
    f_tilde = np.zeros(n)
    v = np.zeros(n+2)
    f_vec = np.zeros(n)
    for i in range(0, n):
        f_vec[i] = f(x[i+1])*h**2
    
    b_tilde[0] = b[0]
    f_tilde[0] = f_vec[0]

    # Forward substitution
    t0 = time.time()
    for i in range(1, n):
        b_tilde[i] = b[i] - (a[i-1]*c[i-1])/b_tilde[i-1]
        f_tilde[i] = f_vec[i] - (f_tilde[i-1]*a[i-1])/b_tilde[i-1]

    # Backward substitution
    v[n] = f_tilde[n-1]/b_tilde[n-1]
    for i in range(n-1, 0, -1):
        v[i] = (f_tilde[i-1] - c[i-1]*v[i+1])/b_tilde[i-1]
    
    t1 = time.time()
    timer = t1 - t0
    print("Time used: %s s" %timer )
    
    return v

def special_algorithm(n):
    x = np.linspace(0, 1., n)
    h = 1./(n+1)
    
    f = 100.*np.exp(-10.*x)*h**2
    f_tilde = np.zeros(n)
    v = np.zeros(n)
    
    f_tilde[0] = f[0]
    
    # Forward substitution
    for i in range(1, n-1):
        f_tilde[i] = f[i] - (i-1)*f_tilde[i-1]/i
    
    # Backward substitution
    v[n-1] = f_tilde[n-1]/b_tilde[n-1]
    for i in range(n-2, 1):
        v[i] = (f_tilde[i] + v[i+1])*(i-1)/i
        
    return v
        
def analytic(x):
    return 1 - (1 - np.exp(-10)) * x - np.exp(-10 * x)
    

if __name__=="__main__":
    N = [10, 100, 1000]
    
    for n in N:
        plt.figure()
        x = np.linspace(0, 1., n+2)
        plt.plot(x, gausselim(n), label="Numerical for n=%s" %n)
        plt.plot(x, analytic(x), "--r", label="Analytical for n=%s" %n)
        plt.legend()
    
