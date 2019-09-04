from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as scpl
import time

def f(x):
    """Source term f(x)"""
    return 100.*np.exp(-10*x)

def analytic(x):
    return 1 - (1 - np.exp(-10))*x - np.exp(-10*x)

def gausselim(n):
    a = np.ones(n-1)*(-1)
    b = np.ones(n)*2
    c = np.ones(n-1)*(-1)

    h = 1./(n+1)

    x = np.zeros(n+2)
    b_tilde = np.zeros(n)
    f_tilde = np.zeros(n)
    f_vec = np.zeros(n)
    v = np.zeros(n+2)

    for i in range(0, n+2):
        x[i] = i*h

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

    return x, v, timer

def special_algorithm(n):
    h = 1./(n)

    b_tilde = np.zeros(n+1)
    f_tilde = np.zeros(n+1)
    v = np.zeros(n+1)
    x = np.zeros(n+1)

    b_tilde[0] = b_tilde[-1] = 2.

    for i in range(1, n):
        b_tilde[i] = (i+1.)/i

    for i in range(0, n+1):
        x[i] = i*h
        f_tilde[i] = f(i*h)*h**2

    # Forward substitution
    t0 = time.time()
    for i in range(2, n):
        f_tilde[i] = f_tilde[i] + f_tilde[i-1]/b_tilde[i-1]

    # Backward substitution
    v[n-1] = f_tilde[n-1]/b_tilde[n-1]
    for i in range(n-2, 0, -1):
        v[i] = (f_tilde[i] + v[i+1])/b_tilde[i]

    t1 = time.time()
    timer = t1 - t0
    #print("Time used special: %s s" %timer )

    return x, v, timer

def rel_error(v, u):
    error = np.max(np.abs((v[1:-1] - u[1:-1])/u[1:-1]))
    return error

def LU(n):
    """LU-decomposition"""
    matrix = np.array([[0 for i in range(n)] for j in range(n)])
    h = 1 / (n+1)
    x = np.array([i * h for i in range(n)])
    f = 100*np.exp(-10*x)*h**2

    for i in range(n):
        matrix[i][i] = 2
        if i != 0:
            matrix[i][i-1] = -1
        if i != n-1:
            matrix[i][i+1] = -1

    t0 = time.time()
    A = scpl.lu_solve(scpl.lu_factor(matrix), f)
    t1 = time.time()
    timer = t1 - t0

    return A, timer

def task_b():
    N = [10, 100, 1000]
    for n in N:
        plt.figure()
        x, v, time1 = gausselim(n)
        print("Time used general: %s s" %time1)
        plt.plot(x, v, label="Numerical for n=%s" %n)
        plt.plot(x, analytic(x), "--r", label="Analytical for n=%s" %n)
        plt.xlabel("x"); plt.ylabel("u(x)")
        plt.legend()

def task_c():
    N = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6]
    for n in N:
        n = int(n)
        time_gen = gausselim(n)[2]
        time_spec = special_algorithm(n)[2]
        print("n=%.1e" %n)
        print("Time used general: %.3e s" %time_gen)
        print("Time used special: %.3e s" %time_spec)
        time_diff = time_gen - time_spec
        print("Time diff: %.3e" %time_diff)
    """
    N = [10, 100, 1000]
    for n in N:
        plt.figure()
        x, v, time2 = special_algorithm(n)
        plt.plot(x, v, label="Numerical for n=%s" %n)
        plt.plot(x, analytic(x), "--r", label="Analytical for n=%s" %n)
        plt.xlabel("x"); plt.ylabel("u(x)")
        plt.legend()
    """

def task_d():
    """
    Compute the relative error between the analytical solution and the Thomas
    algorithm, and the special algorithm (Thöpnitz?).
    """
    N = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]
    error_gen1 = []
    h_list = []
    error_spec1 = []
    for n in N:
        n = int(n)
        h = 1./(n+1)
        h_list.append(h)
        x, v_gen, time1 = gausselim(n)
        u_analytic = analytic(x)
        error_gen = rel_error(v_gen, u_analytic)
        error_gen1.append(error_gen)
        print("n=%.1e: Error general= %.8e" %(n, error_gen))

    for n in N:
        n = int(n)
        x, v_spec, time2 = special_algorithm(n)
        u_analytic = analytic(x)
        error_spec = rel_error(v_spec, u_analytic)
        error_spec1.append(error_spec)
        print("n=%.1e: Error special= %.8e" %(n, error_spec))

    error_arr_gen = np.array(error_gen1)
    error_arr_spec = np.array(error_spec1)
    h_arr = np.array(h_list)

    plt.figure()
    plt.loglog(h_arr, error_arr_gen)
    plt.xlabel("$log_{10}(h)$")
    plt.ylabel("$\epsilon$")

    plt.figure()
    plt.loglog(h_arr, error_arr_spec)
    plt.xlabel("$log_{10}(h)$")
    plt.ylabel("$\epsilon$")

def task_e(n):
    N = [1e1, 1e2, 1e3]
    for n in N:
        n = int(n)

if __name__ == "__main__":
    #task_b()
    #task_c()
    #task_d()
    task_e()

