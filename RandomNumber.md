## Overview

In this project, we are gonna solve exact-diagonalization problem, which is basically focusing on solving a Hermitian matrix (Hamiltonian), so `lanczos algorithm` would be a good approach to do it.

## Lanczos algorithm

The Lanczos algorithm is an iterative method devised by Cornelius Lanczos that is an adaptation of power methods to find the m "most useful"(tending towards extreme highest/lowest)eigenvalues and eigenvectors of an n x n Hermitian matrix.
<br />
<br />
Practically,we use a random-number generator to select each element of the starting vector and look at the sequence of subspaces, Repeated applications of the matrix "amplify" the components of the vector corresponding to the largest eigenvalues, so by the time  subspace will contain eigenvectors corresponding to large eigenvalues. Simply said, the eigenvectors would be converged.

## Sample Data Format

Here is part of the code, where we define Lanczos:
```
def lanczos(matrix, precision=1e-12, max_iterations=1000):
    '''
    Basic Lanczos algorithm without reorthogonalization

    Args:
        matrix (scipy.sparse.csr_matrix): Matrix for Lanczos iterations
        precision (float): precision for convergence of minimal eigenvalue
        max_iterations (int): maximum amount of Lanczos steps

    Returns:
        (t_eigs, alphas, betas) where:
        t_eigs (numpy.array): final eigenvalues of T-matrix
        alphas (numpy.array): Diagonal elements of the T-matrix
        betas (numpy.array): Off-diagonal elements of the T-Matrix
    '''
```
<br />
Here is the part where we generate the random starting vectors.

```
linear_size = matrix.get_shape()[0]

    alphas = []
    betas = [0.]

    v1 = np.random.rand(linear_size) + 1j*np.random.rand(linear_size)
    v1 /= np.linalg.norm(v1)
    v0 = np.zeros(linear_size, dtype=complex)
    w = np.zeros(linear_size, dtype=complex)

    alpha = 0.
    beta = 0.
```
<br />

And this part perform basic Lanczos iteration step:

```
w = matrix.dot(v1)
        alpha = np.real(np.dot(np.conj(v1),w))
        w = w - alpha*v1 - beta*v0
        v0 = np.copy(v1)
        beta = np.real(np.sqrt(np.dot(np.conj(w),w)))
        v1 = 1/beta*w
        alphas.append(alpha)
        betas.append(beta)

```


