

## Overview

**Lambda Lanczos** calculates the smallest or largest eigenvalue and
the corresponding eigenvector of a symmetric (Hermitian) matrix.


```c++
#include <lambda_lanczos/lambda_lanczos.hpp>
using lambda_lanczos::LambdaLanczos;
/* Some include and using declarations are omitted */

void sample() {
  const int n = 3;
  double matrix[n][n] = { {2.0, 1.0, 1.0},
                          {1.0, 2.0, 1.0},
                          {1.0, 1.0, 2.0} };
  // Its eigenvalues are {4, 1, 1}

  /* Prepare matrix-vector multiplication routine used in Lanczos algorithm */
  auto mv_mul = [&](const vector<double>& in, vector<double>& out) {
    for(int i = 0;i < n;i++) {
      for(int j = 0;j < n;j++) {
        out[i] += matrix[i][j]*in[j];
      }
    }
  };


  /* Execute Lanczos algorithm */
  LambdaLanczos<double> engine(mv_mul, n, true, 1); // Find 1 maximum eigenvalue
  vector<double> eigenvalues;
  vector<vector<double>> eigenvectors;
  engine.run(eigenvalues, eigenvectors);
  //// If C++17 is available, the following notation does the same thing:
  // auto [eigenvalues, eigenvectors] = engine.run()


  /* Print result */
  cout << "Eigenvalue: " << setprecision(16) << eigenvalues[0] << endl;
  cout << "Eigenvector:";
  for(int i = 0; i < n; i++) {
    cout << eigenvectors[0][i] << " ";
  }
  cout << endl;
}
```



### 1. Eigenvalue problem
`LambdaLanczos` class computes maximum (minimum) eigenvalues and
corresponding eigenvectors. Degenerate eigenvalues are taken into account.


### 2. Exponentiation
`Exponentiator` class computes the following type of matrix-vector multiplication:

$$\boldsymbol{v}'=e^{\delta A} \boldsymbol{v},$$

where $A$ is a symmetric (Hermitian) matrix and $\delta$ is a scalar parameter.
This class is based on the same theory of the Lanczos algorithm (Krylov subspace method).

As an application, this class may be used for
time evolution of a quantum many-body system:

$$ \ket{\psi(t+\Delta t)} = e^{-iH\Delta t} \ket{\psi(t)},$$
