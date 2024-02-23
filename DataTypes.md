# Data Types
In this computational project, I will calculate the exact-diagonalization based on hubbard model. `Integer`,`double`  and `floats` will be the major datatypes in my program.

## Overview

First of all, We will have to biuld the hamiltonian by creating its matrix elements; the datatypes of those matrix elements would be `Integer` or `double`.

## Sample Data Format

Here is an example, which is part of the code, of the data format:
```
 //build the hopping matrix
  for(int i=0;i<neighbors.size();i++){
    auto nlist = neighbors[i];
    for(auto n = nlist.begin();n!=nlist.end();n++){
      HH(i,*n)=-t;
    }
  } 
  //HH is hermitian
  for(int i=0;i<L;i++){
    for(int j=0;j<L;j++){
      assert( std::abs(HH(i,j)-HH(j,i)) < 1e-6);
    }
  } 
```
## Data Types for the results
Since the results of diagonaliztion should include eigenvalues and eigenvectors, I also need to take care of the datatypes of them. the datatypes of eigenvectors would be `double`
### Sample Data for results Format
```
// call of eigenvectors function follows:
  std::cout << "\nHead of Eigenvector for the lowest eigenvalue:\n\n";
  std::vector<double>::iterator start = eigen.begin();
  std::vector<double>::iterator end = eigen.begin()+1;
  std::vector<Vector> eigenvectors; // for storing the eigen vectors.
  ietl::Info<double> info; // (m1, m2, ma, eigenvalue, residualm, status).
```
