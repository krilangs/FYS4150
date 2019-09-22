from __future__ import division
import numpy as np
import time


def matrix_b(n, omega, max_rho):
    """Construct the matrix for b)"""
    rho_0 = 0
    rho_max = max_rho
    rho = np.linspace(rho_0, rho_max, n+1)
    h = rho[1] - rho[0]
    V = np.zeros(n+1)
    V[1:] = omega**2*rho[1:]**2
    d = 2./h**2 + V
    a = -1./h**2
    A = np.zeros(shape=(n, n), dtype=np.float64)

    A[range(n), range(n)] = d[1:]
    A[range(1, n), range(n-1)] = a
    A[range(n-1), range(1, n)] = a
    print(A)

    return A, rho

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

    a_kk = matrix[k][k]
    a_ll = matrix[l][l]
    matrix[k,k] = a_kk*c**2 - 2.*matrix[k,l]*c*s + a_ll*s**2
    matrix[l,l] = a_kk*c**2 + 2.*matrix[k,l]*c*s + a_ll*s**2
    matrix[k,l] = 0
    matrix[l,k] = 0

    for i in range(n):
        if i != k & i != l:
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

def offdiag(matrix):
    """Finding the maximum non-diagonal element """
    n = len(matrix)
    max_elem = 1e-12
    for i in range(n):
        for j in range(i+1, n):
            if abs(matrix[i,j]) > max_elem:
                max_elem = abs(matrix[i,j])
                max_k = i
                max_l = j
            else:
                pass

    return max_k, max_l


def solve(matrix, tol, time_take=False):
    """
    Find the eigenvalues and eigenvectors of the given matrix.
    Also calculate the number of similarity transformations needed to get
    all non-diagonal elements to zero within a tolerance.
    Can also check time of for the operations if time_take=True.
    """
    M = np.copy(matrix)
    n = len(M)
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

    eig_val = np.zeros(n)  # Eigenvalues
    for i in range(n):
        eig_val[i] = M[i,i]

    print("Eigenpairs are found after %i similarity transformations" %iterations)
    
    if time_take is True:
        return eig_val, eig_vec, iterations, final_time
    else:
        return eig_val, eig_vec

