from __future__ import division
import numpy as np
from numpy.lib.stride_tricks import as_strided
import main

def test_offdiag():
    """
    Test that the function for finding the maximum off-diagonal element works.
    """
    print("Test maximum off-diagonal")
    matrix = np.array([[7, 2, 1, 6],
                       [9, 9, 12, 3],
                       [5, 9, 14, 7],
                       [2, 7, 8, 5]])
    
    i, j = main.offdiag(matrix)
    n = matrix.shape[0]
    max_off_diag = np.max(as_strided(matrix, (n-1,n+1),\
                        (matrix.itemsize*(n+1), matrix.itemsize))[:,1:])
    
    if i == j:
        raise ValueError("Got a diagonal element")
    if matrix[i,j] != max_off_diag:
        raise ValueError("Got incorrect element %i, should have got %i"\
                         %(matrix[i,j], max_off_diag))
    
def test_eigenvalues():
    """
    Check that the matrix is constructed correctly such that numpy gives
    us the same eigenvalues as the analytical within a tolerance.
    """
    print("Test eigenvalues of matrix")
    n = 8
    M, rho = main.Matrix(n, 1, pot="pot1")
    
    num_val, num_vec = np.linalg.eig(M)
    num_val = np.sort(num_val)
    
    h = 1./(n+1)
    d = 2./h**2
    a = -1./h**2
    j = np.linspace(1,n,n)
    
    # Analytical eigenvalues
    lamb_j = d + 2.*a*np.cos(j*np.pi/(n + 1))

    for i in range(n):
        if abs(num_val[i] - lamb_j[i]) > 1e-10:
            raise ValueError("Numerical eigenvalues does not match the analytical")

def test_Jacobi():
    """
    Test that the Jacobi solver gives the same results as for numpy:
    1) Test the eigenpairs
    2) Test the orthogonality of eigenvectors
    """
    print("Test Jacobi:")
    n = 5
    tol = 1e-8
    time_take = False
    M, rho = main.Matrix(n, 1, pot="pot1")
    
    if time_take is True:
        num_val, num_vec, iterations, final_time = main.solve(M, tol, time_take)
    else: 
        num_val, num_vec = main.solve(M, tol, time_take)

    ind = num_val.argsort()
    num_val = num_val[ind]
    num_vec = num_vec[:, ind]
    
    np_val, np_vec = np.linalg.eig(M)
    np_val = np.sort(np_val)
    np_vec = np.sort(np_vec)
    
    # Check the numerical eigenpairs against the analytical eigenpairs
    for i in range(n):
        if abs(num_val[i] - np_val[i]) > tol:
            raise ValueError("Numerical eigenvalues does not match the analytical")
 
    # Check that the eigenvectors are orthogonal
    for i in range(n):
        for j in range(n):
            inner_prod = np.matmul(np.transpose(num_vec[:,i]), num_vec[:,j])
            if i != j and inner_prod > tol:
                raise ValueError("Orthogonality not preserved")
            if i == j and abs(inner_prod - 1) > tol:
                raise ValueError("Orthogonality not preserved")

if __name__ == "__main__":
    """
    Running all the tests.
    If there is an error, the error is printed and the program stops.
    """
    print("Running the tests:")
    test_offdiag()
    test_eigenvalues()
    test_Jacobi()    # Test for several n (dims)
    print("All the tests are passed.")
