# Particle Physics and Event Reconstruction

## INTRODUCTION
I use a solution for the nonequilibrium dynamics of an interacting disordered system. The approach
adapts the combination of the equilibrium dynamical mean-field theory and the equilibrium coherent potential
approximation methods to the nonequilibrium many-body formalism, using the Kadanoff-Baym-Keldysh complex
time contour.We observe via the time dependence of the potential, kinetic, and total energies
the effect of disorder on the relaxation of the system as a function of final interaction strength. The real-time
approach has the potential to shed light on the fundamental role of disorder in the nonequilibrium dynamics of
interacting quantum systems.

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
Notice that the virtical path is along the imaginary time -iβ, where temperature T = 1/β.
<br />
This is so called the Kadanoff-Baym-Keldysh contour; the interaction quench with the interaction
being switched from U1 = 0 to a finite value U2 = U occurs at
time tquench. The disorder strength W is held fixed.
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


### Algorithm and numerical implementation

The nonequilibrium DMFT+CPA algorithm follows the
self-consistency loop illustrated in figure.
![image](https://github.com/ubsuny/ExactDiagonalisation-CP2P2024/assets/50903294/3ff89b77-c0cd-41df-b062-cfc7afc86fcf)

The loop is started
by setting the hybridization Δ(t, t') to an initial guess (we
use an infinitesimal imaginary number) for the first calculation
of the noninteracting Green’s function on the impurity
in the equation below, which can be considered as the `seed` of the loop:
![image](https://github.com/ubsuny/ExactDiagonalisation-CP2P2024/assets/50903294/7b5d6f8d-90fa-4f92-9424-5f038c07855b)
<br />
And at the end of the process, for each subsequent
iteration, the new hybridization is calculated from the average
Green’s function by the equation:
<br />
![image](https://github.com/ubsuny/ExactDiagonalisation-CP2P2024/assets/50903294/8cc7a039-34fa-4c33-90ac-ff0c4bd7d749)


### Particle Identification
The code identifies particles based on their properties, such as probabilities of being certain particle types (e.g., kaons, pions), charges, and momenta. Particle identification is essential for isolating specific particle decays and studying their properties.

### Invariant Mass Reconstruction
The code reconstructs the invariant masses of particles produced in specific decay channels using their measured momenta. Invariant mass reconstruction allows for precise determination of particle masses and the study of particle decay processes.

### Background Suppression
The code applies selection criteria to suppress background contributions from events not corresponding to the decay channels of interest. Background suppression techniques enhance the signal-to-background ratio, improving the sensitivity and accuracy of the analysis.

## Results
On loading and running the C++ macro on ROOT software, the following output is generated:

![Screenshot 2024-03-25 115722](https://github.com/ubsuny/MLppCollision_CP2P2024/assets/143828394/fa7d0b15-e9cd-4771-b30b-305ad16b4b75)


![Screenshot 2024-03-25 120256](https://github.com/ubsuny/MLppCollision_CP2P2024/assets/143828394/3cbd237b-506d-4479-9a5a-1e24f1304672)

The results successfully reproduce the mass of charged B-mesons i.e 5.279 GeV and also indicate a difference in the number of B+ and B- mesons, the reason for this asymmetry is still not fully understood but is a well-observed phenomenon.
Measured CP Asymmetry is -0.0359549 with an uncertainty of 0.00751989.
