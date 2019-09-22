from __future__ import division
import numpy as np



def matrix_b(n, omega, N):
    rho_0 =0
    rho_max = N
    rho = np.linspace(rho_0, rho_max, n+1)
    h = rho[1] - rho[0]
    V = np.zeros(n+1)
    V[1:] = omega**2*rho[1:]**2
    d = 2./h**2 + V
    a = -1./h**2
    matrix = np.array([[0 for i in range(n-1)] for j in range(n-1)])
    # Fill in the diagonal and the elements along the diagonal in the matrix
    for i in range(0, n-1):
        matrix[i][i] = d[i]
        if i != 0:
            matrix[i][i-1] = a
        if i != n-2:
            matrix[i][i+1] = a
        
    print(matrix)
    EigValues, EigVectors = np.linalg.eig(matrix)
    print(EigValues)
    A = np.zeros(shape=(n,n), dtype=np.float64)

    rho_0 = 0
    rho_n = N
    rho = np.linspace(rho_0, rho_n, n+1)  # quickfix
    h = rho[1]-rho[0]
    V = np.zeros(n+1)
    V[1:] = omega**2*rho[1:]**2
    d = 2/h**2 + V
    e = -1/h**2

    A[range(n), range(n)] = d[1:]
    A[range(1, n), range(n-1)] = e
    A[range(n-1), range(1, n)] = e
    print(A)
    
    return A, rho

def Jacobi_rotate(matrix, vec, k, l):
    n = len(matrix)
    tau = (matrix[l,l] - matrix[k,k])/(2.*matrix[k,l])
    
    if tau >= 0:   # tan(theta)
        t = 1./(tau + np.sqrt(1. + tau**2))
    else:
        t = 1./(-tau + np.sqrt(1. + tau**2))
        
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
    
matrix_b(5, 1., 5)
