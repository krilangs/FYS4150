from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import time

def f(x):
    """Source term f(x)"""
    return 100.*np.exp(-10.*x)

def analytic(x):
    return 1 - (1 - np.exp(-10)) * x - np.exp(-10 * x)

def gausselim(n):
    a = np.ones(n-1)*(-1)
    b = np.ones(n)*2
    c = np.ones(n-1)*(-1)

    x = np.linspace(0, 1., n+2)
    h = 1./(n+1)

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

    return v, timer

def special_algorithm(n):
    x = np.linspace(0, 1., n+2)
    h = 1./(n+1)

    b_tilde = np.zeros(n)
    f_tilde = np.zeros(n)
    v = np.zeros(n+2)
    f_vec = np.zeros(n)
    for i in range(0, n):
        f_vec[i] = f(x[i+1])*h**2

    b_tilde[0] = b_tilde[-1] = 2.
    f_tilde[0] = f_vec[0]

    # Forward substitution
    t0 = time.time()
    for i in range(1, n-1):
        b_tilde[i] = (i+1.)/i
    for i in range(1, n):
        f_tilde[i] = f_vec[i] + f_tilde[i-1]/b_tilde[i-1]

    # Backward substitution
    v[n] = f_tilde[n-1]/b_tilde[n-1]
    for i in range(n-1, 0, -1):
        v[i] = (f_tilde[i-1] + v[i+1])/b_tilde[i-1]

    t1 = time.time()
    timer = t1 - t0
    #print("Time used special: %s s" %timer )

    return v, timer

def rel_error(v, u):
    error = np.max(np.abs((v[1:-1] - u[1:-1])/u[1:-1]))
    return error

def task_b():
    N = [10, 100, 1000]
    for n in N:
        plt.figure()
        x = np.linspace(0, 1., n+2)
        v, time1 = gausselim(n)
        print("Time used general: %s s" %time1)
        plt.plot(x, v, label="Numerical for n=%s" %n)
        plt.plot(x, analytic(x), "--r", label="Analytical for n=%s" %n)
        plt.legend()

def task_c():
    N = [10, 1e2, 1e3, 1e4, 1e5, 1e6]
    for n in N:
        n = int(n)
        time1 = gausselim(n)[1]
        time2 = special_algorithm(n)[1]
        print("n=%.1e" %n)
        print("Time used general: %.3e s" %time1)
        print("Time used special: %.3e s" %time2)
        time_diff = time1 - time2
        print("Time diff: %.3e" %time_diff)

def task_d():
    N = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]
    for n in N:
        n = int(n)
        x = np.linspace(0, 1., n+2)
        u_analytic = analytic(x)
        v_gen, time1 =  gausselim(n)
        error_gen = rel_error(v_gen, u_analytic)
        print("n=%.1e: Error general= %.8e" %(n, error_gen))

    for n in N:
        n = int(n)
        x = np.linspace(0, 1., n+2)
        u_analytic = analytic(x)
        v_spec, time2 = special_algorithm(n)
        error_spec = rel_error(v_spec, u_analytic)
        print("n=%.1e: Error special= %.8e" %(n, error_spec))

def task_e():
    ...

if __name__ == "__main__":
    #task_b()
    #task_c()
    task_d()
    #task_e()
