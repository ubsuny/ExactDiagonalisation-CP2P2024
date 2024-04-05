# Modify my data into tensorflow

## INTRODUCTION
My main code is writen in fortron, so the output data are all `.dat` files.
What I want to do is modify these `.dat` files into tensorflow format, which can be used and analysised in python!
## MODEL AND METHODS


### Model
We are interested in an interacting disordered system that
can be described in equilibrium by the single-band Anderson-
Hubbard model defined by
<br />
![image](https://github.com/ubsuny/ExactDiagonalisation-CP2P2024/assets/50903294/dc190fcd-1480-4c03-b872-03789f924f9d)

The first term represents the kinetic energy, the second term
the interaction U between electrons, V describes the random
disorder potential, and μ is the chemical potential.U is the Coulomb
interaction at a doubly occupied site, and Vi is the local onsite disorder potential randomly distributed
according to a probability distribution P(Vi ).We use a box distribution
P(Vi ) = (1/2W)*<b>θ</b>(W − |Vi|). 
<br />
My analysis of the Anderson-Hubbard model will be
done on the Bethe lattice with a large coordination number, while We will study this system at half-filling.


### Nonequilibrium formalism
For the nonequilibrium many-body formalism, starting at
an initial time t<sub>min</sub>, the system is evolved forward in time to
times of physical interest up to a maximal time t<sub>max</sub> and then
backwards again to the initial time t<sub>min</sub>.
<br />
![image](https://github.com/ubsuny/ExactDiagonalisation-CP2P2024/assets/50903294/1bfb6f1f-be54-4c63-99e3-767738dd874a)
<br />
Notice that the virtical path is along the imaginary time -iβ, where temperature T = 1/β.This is so called the Kadanoff-Baym-Keldysh contour; the interaction quench with the interaction being switched from U1 = 0 to a finite value U2 = U occurs at time tquench. The disorder strength W is held fixed.
<br />
To describe such system, the formalism can be formulated either explicitly
in terms of the different Green’s functions.
![image](https://github.com/ubsuny/ExactDiagonalisation-CP2P2024/assets/50903294/9030b336-abb9-4e6e-a251-dab2733472de)

![image](https://github.com/ubsuny/ExactDiagonalisation-CP2P2024/assets/50903294/82b1377f-d71e-4521-87d0-218fc1f1cede)
<br />
So, simply said, the imaginary part of this retarded green's function is the density of states. and that's our main target of this project.
<br />
![image](https://github.com/ubsuny/ExactDiagonalisation-CP2P2024/assets/50903294/e4183f44-62d9-4aaf-85c6-5ff7f7e32aeb)



## PROGRAMMING AND METHODS


### Algorithm of self-consistency loop

The nonequilibrium DMFT+CPA algorithm follows the
self-consistency loop illustrated in figure.
![image](https://github.com/ubsuny/ExactDiagonalisation-CP2P2024/assets/50903294/3ff89b77-c0cd-41df-b062-cfc7afc86fcf)
```
do while(.not. converged .and. (iter .le. iter_max)) !self-consistent loop
