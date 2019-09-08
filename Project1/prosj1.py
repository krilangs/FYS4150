from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as scpl
import time

def f(x):
    """Source term f(x)"""
    return 100.*np.exp(-10*x)

def analytic(x):
    """Analytic solution"""
    return 1 - (1 - np.exp(-10))*x - np.exp(-10*x)

def gausselim(n):
    """Solve a general tridiagonal matrix using the Thomas alogrithm for
    Gaussian elimination of a set of linear equations,
    and take the CPU time of the algorithm."""
    a = np.ones(n-1)*(-1)
    b = np.ones(n)*2
    c = np.ones(n-1)*(-1)

    h = 1./(n)  # Stepsize

    x = np.array([i * h for i in range(n+1)])
    d = np.zeros(n)
    y_tilde = np.zeros(n)
    f_vec = np.zeros(n)
    v = np.zeros(n+1)

    for i in range(0, n):
        f_vec[i] = f(x[i+1])*h**2

    # Initial conditions
    d[0] = b[0]
    y_tilde[0] = f_vec[0]

    t0 = time.time()   # Start timer
    # Forward substitution
    for i in range(1, n):
        factor = a[i-1]/d[i-1]
        d[i] = b[i] - (c[i-1]*factor)
        y_tilde[i] = f_vec[i] - (y_tilde[i-1]*factor)

    # Backward substitution
    v[n-1] = y_tilde[n-1]/d[n-1]
    for i in range(n-2, 0, -1):
        v[i] = (y_tilde[i-1] - c[i-1]*v[i+1])/d[i-1]

    t1 = time.time()   # Stop timer
    timer = t1 - t0    # Time of running the algorithm [s]

    return x, v, timer

def special_algorithm(n):
    """Solve the special tridiagonal Töplitz matrix using Gaussian elimination,
    and take the CPU time of the alogrithm."""
    h = 1./(n)    # Step size

    d = np.zeros(n+1)
    y_tilde = np.zeros(n+1)
    v = np.zeros(n+1)
    x = np.array([i * h for i in range(n+1)])

    d[0] = d[-1] = 2.

    for i in range(1, n):
        d[i] = (i+1.)/i

    for i in range(0, n+1):
        y_tilde[i] = f(x[i])*h**2

    t0 = time.time()   # Start timer
    # Forward substitution
    for i in range(2, n):
        y_tilde[i] = y_tilde[i] + y_tilde[i-1]/d[i-1]

    # Backward substitution
    v[n-1] = y_tilde[n-1]/d[n-1]
    for i in range(n-2, 0, -1):
        v[i] = (y_tilde[i] + v[i+1])/d[i]

    t1 = time.time()    # Stop timer
    timer = t1 - t0     # Time of running the algorithm [s]
    #print("Time used special: %s s" %timer )

    return x, v, timer

def rel_error(v, u):
    """Calculate the relative error in the given data set
    as function of log10(h)."""
    error = np.max(np.abs((v[1:-1] - u[1:-1])/u[1:-1]))
    return error

def LU(n):
    """Calculate the LU-decomposition of the tridiagonal matrix,
    and take the CPU time of the calculation."""
    # Create a matrix filled with zeros
    matrix = np.array([[0 for i in range(n-1)] for j in range(n-1)])
    h = 1 / (n)   # Step size
    x = np.array([i * h for i in range(n-1)])
    func = f(x)*h**2
    # Fill in the diagonal and the elements along the diagonal in the matrix
    for i in range(0, n-1):
        matrix[i][i] = 2
        if i != 0:
            matrix[i][i-1] = -1
        if i != n-2:
            matrix[i][i+1] = -1

    t0 = time.time()    # Start timer
    A = scpl.lu_solve(scpl.lu_factor(matrix), func)  # LU-decomposition solver
    t1 = time.time()    # Stop timer
    timer = t1 - t0     # Time of running the algorithm [s]

    v = np.zeros(n+1)
    v[1:-1] = A

    return x, v, timer

def task_b():
    """
    Compare the solution from the general tridiagonal matrix with the
    analytical solution, and plot for different sizes of the matrices and
    number of grid points.
    """
    N = [10, 100, 1000]
    for n in N:
        plt.figure()
        x, v, time1 = gausselim(n)
        print("Time used general: %s s" %time1)
        plt.title("Gaussian elimination with different grid points\n \
        compared with the analytical solution", size=14)
        plt.plot(x, v, label="Numerical for n=%s" %n, markersize=14)
        plt.plot(x, analytic(x), "--r", label="Analytical for n=%s" %n,\
                 markersize=14)
        plt.xlabel("x", size=14); plt.ylabel("u(x)", size=14)
        plt.legend()
        #plt.savefig("Gaussian_n=%s.png" %n)
    plt.show()

def task_c():
    """
    Compare the CPU time of the special algorithm with the CPU time of the
    Thomas algorithm up to n=1e6 grid points.
    """
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

def task_d():
    """
    Compute the relative error between the analytical solution and the Thomas
    algorithm, and the special algorithm of the Töplitz matrix.
    Then plot the relative errors.
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
    plt.title("Relative error with use of the Thomas algorithm", size=14)
    plt.loglog(h_arr, error_arr_gen)
    plt.xlabel("$log_{10}(h)$", size=14)
    plt.ylabel("$\epsilon=log_{10}(|\\frac{v-u}{u}|)$", size=14)
    #plt.savefig("Error_general.png")
    plt.figure()
    plt.title("Relative error with use of the specialized algorithm", size=14)
    plt.loglog(h_arr, error_arr_spec)
    plt.xlabel("$log_{10}(h)$", size=14)
    plt.ylabel("$\epsilon=log_{10}(|\\frac{v-u}{u}|)$", size=14)
    #plt.savefig("Error_special.png")
    plt.show()

def task_e():
    """
    Calcualte the LU decomposition of the matrix for different grid points,
    and compare the CPU time and relative error with the previously used
    algorithms.
    """
    N = [10, 100, 1000]
    error_gen1 = []
    error_LU1 = []
    time_gen1 = []
    time_LU1 = []
    for n in N:
        x, v, time_gen = gausselim(n)
        u_analytic = analytic(x)
        error_gen = rel_error(v, u_analytic)
        error_gen1.append(error_gen)
        time_gen1.append(time_gen)

        x_, A, time_LU = LU(n)
        u_analytic = analytic(x)
        error_LU = rel_error(A, u_analytic)
        error_LU1.append(error_LU)
        time_LU1.append(time_LU)
        plt.figure()
        plt.title("Gaussian elimination vs LU decomposition\n \
                  with different grid points", size=14)
        plt.plot(x, v, "--r", label="Numerical n=%s" %n, markersize=14)
        plt.plot(x, A, label="LU n=%s" %n, markersize=14)
        plt.legend()
        #plt.savefig("Gaussian_vs_LU_n=%s.png" %n)

    print(error_gen1)
    print(error_LU1)
    print(time_gen1)
    print(time_LU1)
    plt.show()

if __name__ == "__main__":
    #task_b()
    #task_c()
    #task_d()
    task_e()
