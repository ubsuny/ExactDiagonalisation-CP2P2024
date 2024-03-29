program main

use mpi

implicit none

integer, parameter :: dp = kind(1.0d0)
integer, parameter :: sp = kind(1.0)
integer, parameter :: wp = dp
real(kind=wp), parameter :: pi = 3.14159265358979324_wp
complex(kind=wp), parameter :: ii = (0.0_wp, 1.0_wp)

real(kind=wp) :: delta_t, t_max, t_min, delta_tau, beta, mu, er_min, er_min2,&
&W, V, Pv, delta_d, alpha1, t_hop, omega_min, omega_max, t_ave, t_rel, U1, U2,&
&delta_eps, emin, emax, errr2, xA, xB, wA, wB, t2
real(kind=wp), dimension(:), allocatable :: disorder, d_weight, NIDOS, omega, U,&
&eps, Ekin, Epot, Ekinfinal

integer :: i, j, k, i1, i2, N_t, N_tau, N_c, info, iter, N_d, N_omega,&
&numprocs, myid, ierr, idbegin, idend, exitint, iave, irel, i_offset, iter_max, N_e, SCSOiter&
&, itermaxU
integer, dimension(:), allocatable :: ipiv

character(LEN=64) :: filename, strA, strB, strC, strD
character(LEN=4) :: str1, str4
character(LEN=2) :: str2, str3

complex(kind=wp) :: integral_temp, int_temp2, temp1
complex(kind=wp), dimension(:), allocatable :: tc, w_i, workarray, NIDOS_t, DOS
complex(kind=wp), dimension(:,:), allocatable :: mjk, temp, delta_ii, g_imp, g_ret, g_old, &
&deltaij_i, deltaij, g_ave, temp2, Sigma, deltaij_j, gamma_ave, vterm, vtermfin, g_less, temp3, distfn, distfntotal

logical converged, SCSOconverged, binarysolver

!MPI initialization
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

if(myid .eq. 0) write(*,*) 'numprocs =', numprocs
!write(*,*) 'myid =', myid

!call readin
if(myid .eq. 0) call readin 
if(myid .eq. 0) exitint = 0
 
!make sure that number of processes is compatible with N_d

if(myid .eq. 0)then
     if(numprocs .gt. N_d**2)then
 
         write(*,*) 'Error! Numprocs greater than N_d^2'
         exitint = 1
 
     end if
 
     if(mod(N_d**2,numprocs) .ne. 0)then
 
         write(*,*) 'Error! Number of processes does not evenly divide N_d'
         write(*,*) 'Mod(max(N_d**2,numprocs),min(N_d**2,numprocs)) =', mod(max(N_d**2,numprocs),min(N_d**2,numprocs))
 
         exitint = 1
 
     end if

     if(exitint .ne. 0) call MPI_Abort(MPI_COMM_WORLD,exitint,ierr)
end if

! if(exitint .ne. 0) call MPI_Abort(MPI_COMM_WORLD,exitint, ierr)

alpha1 = 0.7
omega_min = -4.0
omega_max = 4.0
N_omega = 200
iter_max = 500

!set up file name strings -- this is currently not optimal, needs to include both U1 and U2

if(myid .eq. 0)then

   write(strA,*) N_t
   write(strB,*) 10*U2
   write(strC,*) 10*W
   write(strD,*) N_tau

   strA = adjustl(strA)
   strB = adjustl(strB)
   strC = adjustl(strC)
   strD = adjustl(strD)

   str1 = trim(strA(1:4))
   str2 = trim(strB(1:2))
   str3 = trim(strC(1:2))
   str4 = trim(strD(1:4))

   if(W .eq. 1.0) str3 = '1'
   if(U2 .eq. 1.0) str2 = '1'
   str2 = trim(str2)
   if(U2 .lt. 1.0) str2 = '0'//str2
   if(W .lt. 1.0) str3 = '0'//str3
   if(abs(U2 - 0.125) .lt. 0.01) str2 = '0125'
   if(abs(W - 0.125) .lt. 0.01) str3 = '0125'

   write(*,*) str1
   write(*,*) str2
   write(*,*) str3
   write(*,*) str4

   write(*,*) 'N_t =', N_t, 'U1 =', U1, 'U2 =', U2, 'W =', W, 'beta =', beta
end if


!broadcast input data
call MPI_BCAST(t_min,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(t_max,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(N_t,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(N_tau,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(beta,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(t_hop,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(W,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(N_d,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(U1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(U2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(er_min2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(itermaxU, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
N_e = N_d

!set up elements derived from input data
delta_t = (t_max-t_min)/N_t
delta_tau = beta/N_tau
delta_d = W/real(N_d)
i_offset = int(-t_min/N_t)
N_c = 2*N_t + N_tau
mu = 0.0
! mu = 0.5*U !doing this at half-filling, so mu = 0. In order to maintain particle number, we simply let the hartree-fock term cancel the chemical potential.

allocate(tc(1:2*N_t + N_tau), mjk(1:2*N_t + N_tau, 1:2*N_t + N_tau), w_i(1:2*N_t + N_tau),&
& workarray(1:N_c), ipiv(1:N_c), temp(1:N_c, 1:N_c), delta_ii(1:N_c, 1:N_c), g_imp(1:N_c, 1:N_c),&
&disorder(1:N_d), g_ret(1:N_t, 1:N_t), d_weight(1:N_d), g_old(1:N_c, 1:N_c), deltaij_i(1:N_c, 1:N_c),&
&deltaij(1:N_c, 1:N_c), NIDOS_t(1:N_t), omega(1:N_omega), temp3(1:N_c, 1:N_c)&!NIDOS(1:N_omega), 
&, temp2(1:N_c, 1:N_c), Sigma(1:N_c, 1:N_c), DOS(1:N_omega), U(1:N_c), deltaij_j(N_c, N_c),&
&gamma_ave(N_c,N_c), vterm(1:N_c, 1:N_c), vtermfin(1:N_c, 1:N_c), NIDOS(1:N_e), eps(1:N_e))

allocate(g_ave(1:N_c, 1:N_c), Ekin(1:N_t), Epot(1:N_t), g_less(1:N_t, 1:N_t), Ekinfinal(1:N_t),&
&distfn(1:N_e, 1:N_t), distfntotal(1:N_e, 1:N_t))


!bethe lattice DOS setup

emax = 2*t_hop
emin = -2*t_hop
delta_eps = (emax - emin)/N_e

do i = 1, N_e

    eps(i) = -0.5 + delta_eps*(i - 0.5)

    NIDOS(i) = sqrt((2*t_hop)**2 - eps(i)**2)/(2*pi*t_hop**2)

end do

!disorder setup
do i = 1, N_d
    disorder(i) = -0.5*W + (i - 0.5)*delta_d
end do

do i = 1, N_d
    d_weight(i) = 1.0/real(N_d)
end do

! write(*,*) d_weight(N_d/2)

do i = 1, N_omega
    omega(i) = omega_min + i*(omega_max - omega_min)/N_omega
end do

!build contour time grid and delta functions
do i = 1, N_t
    tc(i) = t_min + (i-1.0)*delta_t
    w_i(i) = delta_t
end do

do i = N_t + 1, 2*N_t 
    tc(i) = t_max - (i-N_t-1.0)*delta_t!
    w_i(i) = -1.0_wp*delta_t
end do

do i = 2*N_t + 1, N_c
    tc(i) = t_min - ii*delta_tau*(i-2*N_t-1.0)
    w_i(i) = -1.0_wp*ii*delta_tau
end do

deltaij_i = 0
deltaij_j = 0
deltaij = 0

do i = 1, N_c - 1
    deltaij_i(i+1, i) = 1.0/w_i(i+1)
    deltaij_j(i, i+1) = 1.0/w_i(i)
    deltaij(i, i) = 1.0/w_i(i)
end do
deltaij_i(1, N_c) = (1.0/delta_t)*mu
deltaij_j(1, N_c) = (1.0/delta_t)*mu
deltaij(1, N_c) = (1.0/delta_t)*mu!1.0/w_i(N_c)
!write(*,*) 'grids built'


! Here we set up the array for U(t) on the contour. 
! If the real part of the contour time is positive, we take the second value, otherwise we take the first

t2 = 5

U = 0

if(.true.)then
do i = 1, N_c

    if(real(tc(i)) .ge. 0)then
        U(i) = U2
        else
        U(i) = U1
    end if

end do
else

do i = 1, N_c


!if(i .lt. 2*N_t)then
   U(i) = U2*exp(((real(tc(i)) - 5)**2)/(2))
!end if
!   if(real(tc(i)) .ge. 0)then
!      if(real(tc(i)) .le. t2)then
!         U(i) = U2/t2*real(tc(i)) 
!         else
!            U(i) = U2
!      end if
!      else
!         U(i) = U1
!  end if

end do

end if


!build the derivative of the delta function -- leave mu out because we're at half-filling
mjk = 0.0

do i = 1, N_c
    mjk(i,i) = 1.0
end do

do i = 1, N_c - 1
    mjk(i+1, i) = -1.0! - ii*w_i(i+1)*mu
!    mjk(i, i+1) = 1.0
end do

mjk(1,N_c) = 1! + ii*delta_t*mu

do i = 1, N_c
    do j = 1, N_c
    
        mjk(i,j) = ii*mjk(i, j)*(1.0/w_i(i))*(1.0/w_i(j))
    
    end do
end do


!define the starting and end points for each process in the disorder grid
idbegin = int((myid)*((N_d) / numprocs) + 1)
idend = int((myid+1)*((N_d) / numprocs))

write(*,*) 'myid =', myid, 'idbegin and idend =', idbegin, idend

g_ave = 0
vterm = 0.0
Sigma = 0.0
gamma_ave = 0.0
!gamma_i = 0
vtermfin = 0.0

iter = 1

if(myid .eq. 0)delta_ii = 0.00001*ii
converged  = .false.

do while(.not. converged .and. (iter .le. iter_max)) !self-consistent loop


    if(myid .eq. 0 .and. iter .gt. 1)then
      call check_convergence(converged, g_ave, g_old, errr2, er_min)
      write(*,*) 'iter =', iter, 'error=', errr2
    end if
    
        g_old = g_ave
    call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_BCAST(delta_ii,N_c**2,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
    
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
        do j = 1, N_c!build second-order Hubbard self-energy
            do k = 1, N_c
                Sigma(j,k) = U(j)*U(k)*temp2(j,k)*temp2(j,k)*temp2(k,j)!
            end do
        end do
        
        temp2 = temp - Sigma!- deltaij_i*V!temp holds inverse of g-script


        !invert to get the impurity Green's function for this value of the disorder, held in temp2
        do j = 1, N_c
            do k = 1, N_c
                temp2(j,k) = temp2(j,k)*w_i(j)*w_i(k)
            end do
        end do
        
        ipiv = 0.0
        info = 0
        workarray = 0.0
        call zgetrf(n_c, n_c, temp2, n_c, ipiv, info)
        if(info .ne. 0) write(*,*) 'zgetrf returns ', info
        call zgetri(n_c, temp2, n_c, ipiv, workarray, n_c, info)
         
         end if
        
        g_imp = g_imp + temp2*d_weight(i1)
        !g_imp holds the sum of all impurity GFs for this process
 
    end do

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    !collect all impurity GFs into the average GF
    call MPI_Reduce(g_imp, g_ave, N_c*N_c, MPI_DOUBLE_COMPLEX, MPI_SUM, &
      0, MPI_COMM_WORLD, ierr)
      
        
        g_imp = 0

    !the new hybridization is t^2 * g_ave
!     if(myid .eq. 0) delta_ii = (1.0-alpha1)*delta_ii + alpha1*(t_hop*t_hop*g_ave)
    if(myid .eq. 0) delta_ii = t_hop**2 * g_ave
       
    iter = iter + 1
end do !self-consistent loop end
! 
if(myid .eq. 0)then 

    temp2 = g_ave

    do j = 1, N_c !invert GF because G = (i dt - Sigma)^-1                                                                                            
        do k = 1, N_c
            temp2(j,k) = temp2(j,k)*w_i(j)*w_i(k)
            end do
    end do

    ipiv = 0.0
    info = 0
    workarray = 0.0
    call zgetrf(n_c, n_c, temp2, n_c, ipiv, info)
    if(info .ne. 0) write(*,*) 'zgetrf returns ', info
    call zgetri(n_c, temp2, n_c, ipiv, workarray, n_c, info)
    if(info .ne. 0) write(*,*) 'zgetri returns ', info

    sigma = mjk - temp2 - delta_ii !sigma now holds the full self-energy, including disorder effects

end if

call MPI_BCAST(sigma,N_c**2,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
Ekin = 0.0
Ekinfinal = 0.0
distfn = 0.0
distfntotal = 0.0

idbegin = int((myid)*((N_e) / numprocs) + 1)
idend = int((myid+1)*((N_e) / numprocs))
!write(*,*) 'generating Ekin, myid =', myid, 'delta_eps =', delta_eps


    do k = idbegin, idend
        
            temp2 = mjk - Sigma - deltaij*eps(k)!(deltaij_i+deltaij_j)*eps(k)

            do i = 1, N_c
                do j = 1, N_c
                    temp2(i,j) = temp2(i,j)*w_i(i)*w_i(j)
                end do
            end do
            
            ipiv = 0.0
            info = 0
            workarray = 0.0
            call zgetrf(n_c, n_c, temp2, n_c, ipiv, info)
            if(info .ne. 0) write(*,*) 'in get_inverse zgetrf returns ', info
            call zgetri(n_c, temp2, n_c, ipiv, workarray, n_c, info)
            if(info .ne. 0) write(*,*) 'in get_inverse zgetri returns ', info
            
            do i = 1, N_t 
            
                Ekin(i) = Ekin(i) + 2.0*NIDOS(k)*eps(k)*delta_eps*aimag(temp2(i, 2*N_t + 2 - i))!2.0*
                distfn(k,i) = -ii*temp2(i,2*N_t + 2 - i)
            end do
        
    end do

    call MPI_REDUCE(Ekin, Ekinfinal, N_t, MPI_DOUBLE, MPI_SUM, &
      &0, MPI_COMM_WORLD, ierr)

    call MPI_REDUCE(distfn, distfntotal, N_t*N_e, MPI_DOUBLE_COMPLEX, MPI_SUM,&
         &0, MPI_COMM_WORLD, ierr)

    if(myid .eq. 0)then

    do i = 1, N_c !Build G*\Sigma
        do j = 1, N_c
            temp(i,j) = w_i(i)*sigma(i,j)
        end do
    end do

    call zgemm('N', 'N', n_c, n_c, n_c, (1.0_wp, 0.0_wp), g_ave, n_c, temp, n_c, (0.0_wp, 0.0_wp), gamma_ave, n_c)

    if(.false.)then
    temp2 = g_ave

    do i = 1, N_c
       do j = 1, N_c

          temp2(i,j) = w_i(i)*w_i(j)*temp2(i,j)

          end do
          end do

    ipiv = 0.0
    info = 0
    workarray = 0.0
    call zgetrf(n_c, n_c, temp2, n_c, ipiv, info)
    if(info .ne. 0) write(*,*) 'in get_inverse zgetrf returns ', info
    call zgetri(n_c, temp2, n_c, ipiv, workarray, n_c, info)
    if(info .ne. 0) write(*,*) 'in get_inverse zgetri returns ', info
    g_ret = 0.0

    do i = 1, N_t
        do j = 1, N_t
!            g_ret(i,j) = 0.5*(g_ave(i,j) - g_ave(i, 2*N_t -j + 2)&                                                                                   
!                            &+g_ave(2*N_t - i + 2, j) -&                                                                                             
!                            &g_ave(2*N_t - i + 2, 2*N_t - j + 2))                                                                                    
            if(i .ge. j)g_ret(i,j) = temp2(2*N_t - i + 2,j) - temp2(i, 2*N_t - j + 2) !there are multiple ways to do this                             

            end do
    end do


    filename = 'g_ret_inv'//'U'//&
        &trim(adjustl(str2))//'W'//trim(adjustl(str3))//'Nt'//trim(adjustl(str1))//'.dat'
    filename = filename

    open(unit = 37, file = trim(adjustl(filename)), status='unknown')
    do i = 1, N_t
        do j = 1, N_t
            write(37,*) i, j, real(g_ret(i,j)), aimag(g_ret(i,j))
        end do
    end do
    close(37)

    
    end if

    if(.true.)then

    g_ret = 0.0

    do i = 1, N_t
        do j = 1, N_t
!            g_ret(i,j) = 0.5*(g_ave(i,j) - g_ave(i, 2*N_t -j + 2)&                                                                                   
!                            &+g_ave(2*N_t - i + 2, j) -&                                                                                             
!                            &g_ave(2*N_t - i + 2, 2*N_t - j + 2))                                                                                    
            if(i .ge. j)g_ret(i,j) = g_ave(2*N_t - i + 2,j) - g_ave(i, 2*N_t - j + 2) !there are multiple ways to do this                             

            end do
         end do


    !build filename and export local retarded GF
    filename = 'g_ret'//'U'//&
        &trim(adjustl(str2))//'W'//trim(adjustl(str3))//'Nt'//trim(adjustl(str1))//'.dat'
    filename = filename

    open(unit = 37, file = trim(adjustl(filename)), status='unknown')
    do i = 1, N_t
        do j = 1, N_t
            write(37,*) i, j, real(g_ret(i,j)), aimag(g_ret(i,j))
        end do
    end do
    close(37)
    
    write(*,*) 'g_ret exported'

    filename = 'g_c'//'U'//&
        &trim(adjustl(str2))//'W'//trim(adjustl(str3))//'Nt'//trim(adjustl(str1))//'.dat'
    filename = filename

    open(unit = 37, file = trim(adjustl(filename)), status='unknown')
    do i = 1, N_c
        do j = 1, N_c
            write(37,*) i, j, real(g_ave(i,j)), aimag(g_ave(i,j))
        end do
    end do
    close(37)

    do i = 1, N_t
       do j = 1, N_t


!          if(i .gt. j)then
!             temp2(i,j) = g_ave(2*N_t - i + 2, 2*N_t - j + 2)
!             else
!             temp2(i,j) = g_ave(i,j)
!             end if
!          temp2(i,j) = 0.5*(g_ave(i, 2*N_t + 2 - j) - conjg(g_ave(j, 2*N_t + 2 - i)))
          temp2(i,j) =(g_ave(i, 2*N_t + 2 -j))!  - conjg(g_ave(2*N_t + 2 - i, j)))*0.5

          end do
          end do

    filename = 'g_less'//'U'//&
         &trim(adjustl(str2))//'W'//trim(adjustl(str3))//'Nt'//trim(adjustl(str1))//'.dat'
    filename = filename

    open(unit = 37, file = trim(adjustl(filename)), status='unknown')
    do i = 1, N_t
       do j = 1, N_t

          write(37,*), i, j, real(temp2(i,j)), aimag(temp2(i,j))

       end do
    end do
    close(37)

    filename = 'g_greater'//'U'//&
         &trim(adjustl(str2))//'W'//trim(adjustl(str3))//'Nt'//trim(adjustl(str1))//'.dat'
    filename = filename

    open(unit=37, file = trim(adjustl(filename)), status='unknown')

    do i = 1, N_t
       do j = 1, N_t

          write(37,*) i, j, real(g_ave(2*N_t + 2 - i, j)), aimag(g_ave(2*N_t + 2 - i, j))

       end do
    end do


    close(37)


    end if

    do i = 1, N_t
        do j = 1, N_t
    
        temp2(i, j) = gamma_ave(i, 2*N_t - j + 2)! + 0.5*U(i)*g_ave(i, 2*N_t - j + 2)!U(i)*
    
        end do
    end do
    
    do i = 1, N_t

        Epot(i) = aimag(temp2(i,i)) + U(i)*0.25

    end do
    
    
     filename = 'Epot'//'U'//&
      &trim(adjustl(str2))//'W'//trim(adjustl(str3))//'Nt'//trim(adjustl(str1))//'Ntau'//trim(adjustl(str4))//'.dat'
     filename = filename
     
     open(unit = 37, file = trim(adjustl(filename)), status='unknown')
     do i = 1, N_t
             write(37,*) real(tc(i)), real(Epot(i))!, aimag(Epot(i))! + U(i)*0.25! 
     end do
     close(37)    

     filename = 'distfn'//'U'//&
      &trim(adjustl(str2))//'W'//trim(adjustl(str3))//'Nt'//trim(adjustl(str1))//'Ntau'//trim(adjustl(str4))//'.dat'
     filename = filename

     open(unit = 37, file = trim(adjustl(filename)), status='unknown')

     do i = 1 , N_t
        do k = 1, N_e

           write(37,*) real(tc(i)), eps(k), real(distfntotal(k,i)), aimag(distfntotal(k,i))

        end do
     end do

     close(37)
        
    filename = 'Ekin'//'U'//&
     &trim(adjustl(str2))//'W'//trim(adjustl(str3))//'Nt'//trim(adjustl(str1))//'Ntau'//trim(adjustl(str4))//'.dat'
    filename = filename
    
    open(unit = 37, file = filename, status='unknown')
    do i = 1, N_t
            write(37,*) real(tc(i)), Ekinfinal(i)!
    end do
    close(37)
    
    filename = 'Etot'//'U'//&
     &trim(adjustl(str2))//'W'//trim(adjustl(str3))//'Nt'//trim(adjustl(str1))//'Ntau'//trim(adjustl(str4))//'.dat'
    filename = filename
    
    open(unit = 37, file = filename, status='unknown')
    do i = 1, N_t
            write(37,*) real(tc(i)), Ekinfinal(i) + Epot(i)!
    end do
    close(37)
    
    end if

    if(myid .eq. 0)then
    
!     vtermfin = vterm
if(.false.)then
    do i = 3, N_t - 3

       workarray(i) = ((1.0/12.0)*g_less(i-2,i) - (2.0/3.0)*g_less(i-1,i) + (2.0/3.0)*g_less(i+1,i) - (1.0/12.0)*g_less(i+2,i)&
           &) / (delta_t) + U(i)*0.25

   end do

   open(unit = 37, file = 'Etotderiv.dat', status='unknown')

      do i = 3, N_t - 3

         write(37, *) real(tc(i)), real(workarray(i)), aimag(workarray(i))

      end do

   close(37)
end if

if(.false.)then!switch for large file exports
   filename = 'gc'//'U'//trim(adjustl(str2))//'W'//trim(adjustl(str3))//'Nt'//&
        &trim(adjustl(str1))//'.dat'
    filename = filename
   open(unit=37, file = filename, status='unknown')
   do i = 1, N_c
      do j = 1, N_c
         write(37,*) i, j, real(g_ave(i,j)), aimag(g_ave(i,j))
      end do
   end do
   close(37)

   filename = 'sigma_loc'//'U'//trim(adjustl(str2))//'W'//trim(adjustl(str3))//'Nt'&
        &//trim(adjustl(str1))//'.dat'
    filename = filename

    open(unit=37, file = filename, status='unknown')
   do i = 1, N_c
      do j= 1, N_c
         write(37,*) i, j, real(sigma(i,j)), aimag(sigma(i,j))
      end do
   end do
   close(37)
end if
end if
call MPI_FINALIZE(ierr)

contains

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

subroutine readin
    integer :: input
    input = 42
    open(unit=input,file='indat',status='old')
    read(input,*) t_min,t_max,N_t, N_tau, N_d, beta, t_hop, W, U1, U2, er_min, er_min2, alpha1&
&,itermaxU
    close(unit=input)
!     mu = U/2.0

end subroutine readin

end program main
