from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from numba import jit
import time

def Matrix(n, omega, max_rho, pot):
    """Construct the tri-diagonal matrix."""
    rho_0 = 0
    rho_max = max_rho
    rho = np.linspace(rho_0, rho_max, n+1)
    h = rho_max/(n+1)   # Step size
    V = np.zeros(n+1)   # Potential which is exercises dependent
    V[1:] = potentials(rho[1:], omega, pot) 

    d = 2./h**2 + V     # Diagonal elements
    a = -1./h**2        # Non-diagonal around the diagonal
    A = np.zeros(shape=(n, n), dtype=np.float64)  # Create matrix with zeros
    # Fill inn the diagonal and the above and below non-diagonal elements
    A[range(n), range(n)] = d[1:]
    A[range(1, n), range(n-1)] = a
    A[range(n-1), range(1, n)] = a

    return A, rho

@jit(nopython=True)
def Jacobi_rotate(matrix, vec, k, l):
    """Finding the new matrix elements by a Jacobi rotation."""
    n = len(matrix)
    tau = (matrix[l,l] - matrix[k,k])/(2.*matrix[k,l])  # cot(2*theta)

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
    """
    Finding the maximum non-diagonal element of a symmetrical matrix.
    Returns the position indices of the maximum element.
    """
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
    """Define the different potentials used in the different exercises."""
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
    Can also check time of the operations if time_take=True.
    """
    M = np.copy(matrix)
    n = len(M)
    print("Dim %i x %i:" %(n,n))
    iterations = 0     # Number of transformations
    eig_vec = np.identity(n)  # Eigenvectors

    k, l = offdiag(M)   # Indices of the maximum non-diagonal element

    if time_take is True:   # Start timer
        start_time = time.perf_counter()

    while M[k,l]**2 >= tol:
        M, eig_vec = Jacobi_rotate(M, eig_vec, k, l)  # Do the transformations
        k, l = offdiag(M)  # Calculate new indices for the new maximum element
        iterations += 1   # Add iterations

    if time_take is True:   # Stop timer, and calculate the time of operations
        stop_time = time.perf_counter()
        final_time = stop_time - start_time
        print("Time used on the algorithm: %.10s s" %final_time)

    eig_val = np.zeros(n)  # Eigenvalues
    for i in range(n):
        eig_val[i] = M[i,i]

    print("Eigenpairs are found after %i similarity transformations" %iterations)
    
    if time_take is True:
        return eig_val, eig_vec, iterations, final_time
    else:
        return eig_val, eig_vec

def figsetup(title, xlabel, ylabel, fname, show=False):
    """
    Sets up and saves figure for usage in report with the option to not plot.
    """
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.legend()
    plt.savefig(fname + ".png", dpi=250)
    if show is False:
        plt.close()   # Does not plot
    else:
        plt.show()    # Plots the figures

def ex_c(show):
    plot = show
    tol = 1e-8
    list_n = [5, 10, 20, 50, 70, 100, 150, 200, 300, 400, 500, 800]
    n_iterations = []   # List for number of iterations
    time_taken = []     # List for time taken for the algorithm
    
    for n in list_n:
        M, rho = Matrix(n, 1, 1, pot="pot1")
        M_val, M_vec, iterations, time = solve(M, tol, time_take=True)
        n_iterations = np.append(n_iterations, iterations)
        time_taken = np.append(time_taken, time)
    list_n = np.array(list_n)

    plt.figure(figsize=[5, 5])
    plt.semilogy(list_n, time_taken, "o--")
    figsetup(title="Time taken for solution of N-step Buckling beam\n Semilogy-plot",
            xlabel="N", ylabel="Time elapsed [s]", fname="Time_semilogy", show=plot)

    plt.figure(figsize=[5, 5])
    plt.semilogy(list_n, n_iterations, "o--")
    figsetup(title="# Iterations for solution of N-step Buckling beam\n Semilogy-plot",
             xlabel="N", ylabel="# Iterations", fname="Iterations_semilogy", show=plot)

    plt.figure(figsize=[5, 5])
    plt.loglog(list_n, time_taken, "o--")
    figsetup(title="Time taken for solution of N-step Buckling beam\n Loglog-plot",
             xlabel="N", ylabel="Time elapsed [s]", fname="Time_loglog", show=plot)

    plt.figure(figsize=[5, 5])
    plt.loglog(list_n, n_iterations, "o--")
    figsetup(title="# Iterations for solution of N-step Buckling beam\n Loglog-plot",
             xlabel="N", ylabel="# Iterations", fname="Iterations_loglog", show=plot)

    plt.figure(figsize=[5, 5])
    plt.plot(list_n, time_taken, "o--")
    figsetup(title="Time taken for solution of N-step Buckling beam\n Normal plot",
             xlabel="N", ylabel="Time elapsed [s]", fname="Time_plot", show=plot)

    plt.figure(figsize=[5, 5])
    plt.plot(list_n, n_iterations, "o--")
    figsetup(title="# Iterations for solution of N-step Buckling beam\n Normal plot",
             xlabel="N", ylabel="# Iterations", fname="Iterations_plot", show=plot)

def ex_d(show):
    N = 800
    rho_max = 50
    plot = show
    
    M, rho = Matrix(N, 1, rho_max, pot="pot2")
    M_val, M_vec, iterations, time = solve(M, tol=1e-8, time_take=True)
    
    permute = M_val.argsort()
    M_val = M_val[permute]
    M_vec = M_vec[:, permute]
    
    plt.figure(figsize=[5, 5])

    for n in [0, 1, 2]:
        plt.plot(rho[1:], M_vec[:, n]**2, label="$\\lambda=$%.4f" % M_val[n])
        plt.axis([0, 8, 0.0, 0.06])

    figsetup(title="Dimensionless wavefunction for first 3 eigenstates",
             xlabel="$\\rho$", ylabel="$u(\\rho)$", fname="Eigenvalues", show=plot)
    
    print(M_val[:4])

def ex_e(show):
    N = 500
    rho_max = 50
    plot = show
    omega_r = [0.01, 0.5, 1., 5.]

    plt.figure(figsize=[5, 5])

    for w in omega_r:
        M, rho = Matrix(N, w, rho_max, pot="pot3")
        M_val, M_vec = solve(M, tol=1e-8, time_take=False)
    
        permute = M_val.argsort()
        M_val = M_val[permute]
        M_vec = M_vec[:, permute]
        
        np_val, np_vec = np.linalg.eig(M)
        np_val = np.sort(np_val)
        np_vec = np.sort(np_vec)

        plt.plot(rho[1:], M_vec[:, 0], label="$\\omega=$%.2f" %w)
        print(M_val[0])
        print(np_val[0])

    figsetup(title="Dimensionless wavefunction for first eigenstates",
             xlabel="$\\rho$", ylabel="$u(\\rho)$", fname="Frequency",
             show=plot)
    
if __name__ == "__main__":
    #ex_c(show=True)   
    #ex_d(show=True)
    ex_e(show=False)   # Maa ha analytiske foerst
    
