from __future__ import division
import numpy as np
from numba import jit
import matplotlib.pyplot as plt
import time


def Matrix(n, omega, max_rho, pot):
    """Construct the matrix for b)"""
    rho_0 = 0
    rho_max = max_rho
    rho = np.linspace(rho_0, rho_max, n+1)
    h = rho_max/(n+1)
    V = np.zeros(n+1)
    V[1:] = potentials(rho[1:], omega, pot)

    d = 2./h**2 + V
    a = -1./h**2
    A = np.zeros(shape=(n, n), dtype=np.float64)

    A[range(n), range(n)] = d[1:]
    A[range(1, n), range(n-1)] = a
    A[range(n-1), range(1, n)] = a

    return A, rho

@jit(nopython=True)
def Jacobi_rotate(matrix, vec, k, l):
    """Finding the new matrix elements by a Jacobi rotation"""
    n = len(matrix)
    tau = (matrix[l,l] - matrix[k,k])/(2.*matrix[k,l])

    if tau >= 0:   # tan(theta)
        t = 1./(tau + np.sqrt(1. + tau**2))
    else:
        t = -1./(-tau + np.sqrt(1. + tau**2))

    c = 1./np.sqrt(1. + t**2)   # cos(theta)
    s = t*c                     # sin(theta)

    a_kk = matrix[k,k]
    a_ll = matrix[l,l]
    matrix[k,k] = a_kk*c**2 - 2.*matrix[k,l]*c*s + a_ll*s**2
    matrix[l,l] = a_kk*s**2 + 2.*matrix[k,l]*c*s + a_ll*c**2
    matrix[k,l] = 0
    matrix[l,k] = 0

    for i in range(n):
        if i != k and i != l:
            a_ik = matrix[i,k]
            a_il = matrix[i,l]
            matrix[i,k] = a_ik*c - a_il*s
            matrix[k,i] = matrix[i,k]
            matrix[i,l] = a_il*c + a_ik*s
            matrix[l,i] = matrix[i,l]

        # Eigenvectors
        r_ik = vec[i,k]
        r_il = vec[i,l]
        vec[i,k] = r_ik*c - r_il*s
        vec[i,l] = r_il*c + r_ik*s

    return matrix, vec

@jit(nopython=True)
def offdiag(matrix):
    """Finding the maximum non-diagonal element """
    n = len(matrix)
    max_elem = 0
    for i in range(n):
        for j in range(i+1, n):
            if abs(matrix[i,j]) > max_elem:
                max_elem = abs(matrix[i,j])
                max_k = i
                max_l = j
            else:
                pass

    return max_k, max_l

def potentials(rho, omega, pot):
    """Define the different potentials used in the different exercises"""
    if pot == "pot1":  # Exercise b) and c)
        V = 0
    if pot == "pot2":  # Exercise d)
        V = rho**2
    if pot == "pot3":  # Exercise e)
        V = omega**2*rho**2 + 1./rho
    return V
        

def solve(matrix, tol, time_take=False):
    """
    Find the eigenvalues and eigenvectors of the given matrix.
    Also calculate the number of similarity transformations needed to get
    all non-diagonal elements to zero within a tolerance.
    Can also check time of for the operations if time_take=True.
    """
    M = np.copy(matrix)
    n = len(M)
    print("Dim %i x %i:" %(n,n))
    iterations = 0     # Number of transformations
    eig_vec = np.identity(n)  # Eigenvectors

    k, l = offdiag(M)

    if time_take is True:   # Start timer
        start_time = time.perf_counter()

    while M[k,l]**2 >= tol:
        M, eig_vec = Jacobi_rotate(M, eig_vec, k, l)  # Do the transformations
        k, l = offdiag(M)
        iterations += 1

    if time_take is True:   # Stop timer, and calculate the time of operations
        stop_time = time.perf_counter()
        final_time = stop_time - start_time
        print("Time used on the algorithm: %s s" %final_time)

    eig_val = np.zeros(n)  # Eigenvalues
    for i in range(n):
        eig_val[i] = M[i,i]

    print("Eigenpairs are found after %i similarity transformations" %iterations)
    
    if time_take is True:
        return eig_val, eig_vec, iterations, final_time
    else:
        return eig_val, eig_vec

def figsetup(title, xlab, ylab, fname, show=False):
    """
    Sets up and saves figure for usage in report
    usage:
    plot(...)
    plot(...)
    figsetup("filename")
    """
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(fname)
    plt.tight_layout()
    plt.title(title)
    plt.legend()
    #plt.savefig("../figs/" + fname + ".png", dpi=250)
    if show is False:
        plt.close()
    else:
        plt.show()

def ex_c(show):
    plot = show
    tol = 1e-8
    list_n = [5, 10, 20, 50, 70, 100, 150, 200, 300, 400, 500]
    n_iterations = []
    time_taken = []
    
    for n in list_n:
        M, rho = Matrix(n, 1, 1, pot="pot1")
        M_val, M_vec, iterations, time = solve(M, tol, time_take=True)
        n_iterations = np.append(n_iterations, iterations)
        time_taken = np.append(time_taken, time)
    list_n = np.array(list_n)

    plt.figure(figsize=[5, 5])
    plt.semilogy(list_n, time_taken, "o--")
    figsetup(title="Time taken for solution of N-step Buckling beam\n Semilogy-plot",
            xlab="N", ylab="Time elapsed [s]", fname="q2c_time", show=plot)

    plt.figure(figsize=[5, 5])
    plt.semilogy(list_n, n_iterations, "o--")
    figsetup(title="No. Iterations for solution of N-step Buckling beam\n Semilogy-plot",
             xlab="N", ylab="No. Iterations", fname="q2c_count", show=plot)

    plt.figure(figsize=[5, 5])
    plt.loglog(list_n, time_taken, "o--")
    figsetup(title="Time taken for solution of N-step Buckling beam\n Loglog-plot",
             xlab="N", ylab="Time elapsed [s]", fname="q2c_timeloglog", show=plot)

    plt.figure(figsize=[5, 5])
    plt.loglog(list_n, n_iterations, "o--")
    figsetup(title="No. Iterations for solution of N-step Buckling beam\n Loglog-plot",
             xlab="N", ylab="No. Iterations", fname="q2c_countloglog", show=plot)

    plt.figure(figsize=[5, 5])
    plt.plot(list_n, time_taken, "o--")
    figsetup(title="Time taken for solution of N-step Buckling beam\n Normal plot",
             xlab="N", ylab="Time elapsed [s]", fname="q2c_timenormal", show=plot)

    plt.figure(figsize=[5, 5])
    plt.plot(list_n, n_iterations, "o--")
    figsetup(title="No. Iterations for solution of N-step Buckling beam\n Normal plot",
             xlab="N", ylab="No. Iterations", fname="q2c_countnormal", show=plot)

def ex_d(show):
    N = 600
    rho_max = 10
    plot = show
    
    M, rho = Matrix(N, 1, rho_max, pot="pot2")
    M_val, M_vec = solve(M, tol=1e-8, time_take=False)
    
    permute = M_val.argsort()
    M_val = M_val[permute]
    M_vec = M_vec[:, permute]
    
    plt.figure(figsize=[5, 5])

    for n in [0, 1, 2]:
        plt.plot(rho[1:], M_vec[:, n]**2, label="$\\lambda=$%.4f" % M_val[n])
        #plt.axis([0, 15, 0.0, 0.025])

    figsetup(title="Dimensionless wavefunction for first 3 eigenstates",
             xlab="$\\rho$", ylab="$u(\\rho)$", fname="question2d%i" % N,
             show=plot)
    
    print(M_val[:4])

def ex_e(show):
    N = 400
    rho_max = 10
    plot = show
    omega_r = [0.01, 0.5, 1., 5.]

    plt.figure(figsize=[5, 5])

    for w in omega_r:
        M, rho = Matrix(N, w, rho_max, pot="pot3")
        M_val, M_vec = solve(M, tol=1e-8, time_take=False)
    
        permute = M_val.argsort()
        M_val = M_val[permute]
        M_vec = M_vec[:, permute]
        
        #print(M_val[:4])
        
        plt.plot(rho[1:], M_vec[:, 0], label="$\\omega=$%.2f" %w)
    
    figsetup(title="Dimensionless wavefunction for first eigenstates",
             xlab="$\\rho$", ylab="$u(\\rho)$", fname="question2d%i" % N,
             show=plot)
    
if __name__ == "__main__":
    #ex_c(show=True)   
    #ex_d(show=True)
    ex_e(show=True)
    
