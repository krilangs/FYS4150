from __future__ import division
import numpy as np
from numpy.lib.stride_tricks import as_strided
import main

def test_offdiag():
    """
    Test that the function for finding the maximum off-diagonal element works.
    """
    matrix = np.array([[7, 2, 1, 6],
                       [9, 9, 12, 3],
                       [5, 9, 14, 7],
                       [2, 7, 8, 5]])
    
    i, j = main.offdiag(matrix)
    n = matrix.shape[0]
    max_off_diag = np.max(as_strided(matrix, (n-1,n+1), (matrix.itemsize*(n+1), matrix.itemsize))[:,1:])
    
    if i == j:
        raise ValueError("Got a diagonal element")
    if matrix[i,j] != max_off_diag:
        raise ValueError("Got incorrect element %i, should have got %i"\
                         %(matrix[i,j], max_off_diag))
    
def test_eigenpairs():
    """
    Check that the matrix is constructed correctly such that numpy gives
    us the same eigenvalues as the analytical within a tolerance.
    """
    n = 8
    M, rho = main.matrix_b(n, 1, n)
    
    matrix_val, matrix_vec = np.linalg.eig(M)
    
    h = 1./(n+1)
    d = 2./h**2
    a = -1./h**2
    j = np.linspace(1,n,n)
    
    # Analytical eigenvalues
    lamb_j = d + 2.*a*np.cos(j*np.pi/(n + 1)) 

    for i in range(n):
        if abs(matrix_val[i] - lamb_j[i]) > 1e-8:
            raise ValueError("Numerical eigenvalues does not match the analytical")
test_eigenpairs()
