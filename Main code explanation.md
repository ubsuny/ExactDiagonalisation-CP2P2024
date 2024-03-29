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


### Algorithm of self-consistency loop

The nonequilibrium DMFT+CPA algorithm follows the
self-consistency loop illustrated in figure.
![image](https://github.com/ubsuny/ExactDiagonalisation-CP2P2024/assets/50903294/3ff89b77-c0cd-41df-b062-cfc7afc86fcf)
```
do while(.not. converged .and. (iter .le. iter_max)) !self-consistent loop


    if(myid .eq. 0 .and. iter .gt. 1)then
      call check_convergence(converged, g_ave, g_old, errr2, er_min)
      write(*,*) 'iter =', iter, 'error=', errr2
    end if
    
        g_old = g_ave
    call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_BCAST(delta_ii,N_c**2,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
```
The loop is started
by setting the hybridization Δ(t, t') to an initial guess (we
use an infinitesimal imaginary number) for the first calculation
of the noninteracting Green’s function on the impurity
in the equation below, which can be considered as the `seed` of the loop:
![image](https://github.com/ubsuny/ExactDiagonalisation-CP2P2024/assets/50903294/7b5d6f8d-90fa-4f92-9424-5f038c07855b)
<br />
inside the while loop above, the first step is included and the below is the code 
<br />
```
do i1 = idbegin, idend !impurity solver !1, N_d!
        
       V = disorder(i1)

        temp = mjk - delta_ii - deltaij*V!&
!             &0.5*(deltaij_i+deltaij_j)*V !this is inverse of g-script
        
        do j = 1, N_c !prepare to invert temp to get g-script
            do k = 1, N_c
                temp2(j,k) = temp(j,k)*w_i(j)*w_i(k) !temp2 will hold g-script (need to keep inverse around)
            end do
        end do
        
        ipiv = 0.0
        info = 0
        workarray = 0.0
        call zgetrf(n_c, n_c, temp2, n_c, ipiv, info)
        if(info .ne. 0) write(*,*) 'in get_inverse zgetrf returns ', info
        call zgetri(n_c, temp2, n_c, ipiv, workarray, n_c, info)
        
        if(.true.)then
```
<br />

### 
And at the end of the process, for each subsequent
iteration, the new hybridization is calculated from the average
Green’s function by the equation: 
<br />
![image](https://github.com/ubsuny/ExactDiagonalisation-CP2P2024/assets/50903294/e523f87b-405b-437b-bc5c-944c8189b789)
<br />
and the code of this step is below:
```
    !collect all impurity GFs into the average GF
    call MPI_Reduce(g_imp, g_ave, N_c*N_c, MPI_DOUBLE_COMPLEX, MPI_SUM, &
      0, MPI_COMM_WORLD, ierr)
      
        
        g_imp = 0

    !the new hybridization is t^2 * g_ave
    if(myid .eq. 0) delta_ii = t_hop**2 * g_ave
       
    iter = iter + 1
end do !self-consistent loop end
```
And then, we check if the result is converged by the code:
```
subroutine check_convergence(converge, g_in, g_prev, errave, er_minimum)

   logical,intent(out) :: converge
   complex(kind=wp), intent(in) :: g_in(1:N_c, 1:N_c), g_prev(1:N_c, 1:N_c)
   real(kind=wp) :: errr, errrdenom
   real(kind=wp), intent(in) :: er_minimum
   real(kind=wp), intent(out) :: errave

   ! integer :: k

   errr = 0.0
   errrdenom = 0.0
                           

   do i = 1, N_c
      do j = 1, N_c
   !    if((abs(sigma(k)) .ge. eps2) .and. (abs(sigma_old(k)) .ge. eps2))then                                                                            
      errr = errr + abs((g_in(i,j)) - (g_prev(i,j)))
      errrdenom = errrdenom + abs(g_in(i,j) + g_prev(i,j))
   !    end if                                                                                                                                           
      end do
   end do

   errr = 2.0*errr/errrdenom

   errave = errr
   
!    write(*,*) 'averaged error =', errr, 'er_minimum =', er_minimum

   if((errr - er_minimum) .le. 0.0)then
      converge = .true.
!      write(*,*) 'Converged'
end if
 
end subroutine check_convergence
```

### Results and plots

