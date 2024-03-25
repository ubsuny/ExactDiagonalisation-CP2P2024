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
My analysis of the Anderson-Hubbard model will be
done on the Bethe lattice with a large coordination number, while We will study this system at half-filling.

### Nonequilibrium formalism
For the nonequilibrium many-body formalism, starting at
an initial time t<sub>min</sub>, the system is evolved forward in time to
times of physical interest up to a maximal time tmax and then
backwards again to the initial time tmin.

### Detectors
Detectors surrounding the collision points are designed to detect and measure the properties of particles produced in the collisions. These detectors are complex instruments composed of various sub-detectors, each optimized for detecting different types of particles and measuring specific properties such as energy, momentum, and charge.

## Event Reconstruction Process
Event reconstruction is the process of identifying and analyzing particles produced in collider collisions based on the signals recorded by detectors. The event reconstruction process typically involves several steps:

1. **Signal Detection**: Detectors record signals produced by particles interacting with detector materials. These signals include energy deposits, ionization trails, and electromagnetic radiation emitted by particles.

2. **Particle Identification**: Based on the recorded signals, detectors identify the types of particles produced in the collisions. Different types of particles leave distinct signatures in detectors, allowing for particle identification based on their properties such as energy loss, trajectory, and interaction patterns.

3. **Trajectory Reconstruction**: Detectors reconstruct the trajectories of particles based on the signals they produce as they pass through the detector layers. Trajectory reconstruction allows scientists to determine the paths particles took after the collision and infer information about their origins and interactions.

4. **Energy Measurement**: Detectors measure the energy deposited by particles as they interact with detector materials. Energy measurements are crucial for determining the total energy of particles and identifying high-energy particles that may indicate the presence of new physics phenomena.

5. **Interaction Analysis**: Scientists analyze the detected particles' interactions to study fundamental physics phenomena, such as particle decays, particle scattering, and the production of new particles. Analysis of particle interactions provides insights into the underlying physics principles governing particle behavior and interactions.

## Present analysis
The code performs event reconstruction and analysis of data recorded by detectors in collider experiments. It applies selection criteria to identify specific particle decays, reconstructs their invariant masses, and analyzes their properties to study fundamental physics phenomena such as CP Asymmetry.

![Screenshot 2024-03-25 115940](https://github.com/ubsuny/MLppCollision_CP2P2024/assets/143828394/6e31efc0-accc-40da-ae21-985c86a70171)

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
