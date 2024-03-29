program processGLGRGG

implicit none

integer, parameter :: dp = kind(1.0d0)
integer, parameter :: sp = kind(1.0)
integer, parameter :: wp = dp
real(kind=wp), parameter :: pi = 3.14159265358979324_wp
complex(kind=wp), parameter :: ii = (0.0_wp, 1.0_wp)

complex(kind=wp), allocatable, dimension(:,:) :: glessw, gretw, ggreatw,glessw2, gretw2, ggreatw2&
&,glessw3, gretw3, ggreatw3
! complex(kind=wp), allocatable, dimension(:,:) :: glessw, gretw, ggreatw
complex(kind=wp), allocatable, dimension(:,:) :: grett, glesst, ggreatt, GLrel, GRrel, GGrel&
&,grett2, glesst2, ggreatt2, grett3, glesst3, ggreatt3, GLrel2, GLrel3, GRrel2, GRrel3, GGrel2, GGrel3&
&,glesswI1, glesswI2, gretwI1, gretwI2, ggreatwI1, ggreatwI2
real(kind=wp), allocatable, dimension(:) :: time, time2, time3, w&
&, time_ave, time_rel, time_ave2, time_rel2, time_ave3, time_rel3
integer :: Nt, i, j, k, k0, k1, k2, Nw, i1, i2, i3, i_offset, i_slice, i_ave, Ntd, i_rel, i_slice_ave&
&,Nt2, Nt3
integer, allocatable, dimension(:,:,:) :: map, map_2ave
real(kind=wp) :: tmin, tmax, xx, yy, zz, deltat, deltat2, deltat3, L0, L1, L2, dt&
&,rtemp, itemp, wmin, trel, t_slice, tminin, tmaxin
complex(kind=wp), allocatable, dimension(:) :: fw, fw1, fw2, fw3, temparray1, temparray2, temparrayX
complex(kind=wp) :: scatemp
real :: T1, T2

call cpu_time(T1)

tminin = -10.0!-5.0!
tmaxin = 15.0!20.0!

tmin = -abs((tmaxin - tminin)*0.5)!tminin!-abs((tmaxin - tminin)*0.5)!
tmax = abs(tmin)!tmaxin!

Nt = 800
Nt2 = 1000
Nt3 = 1200
Nw = 200
! Ntd = 2*Nt

wmin = -1.0
! i_slice = Nt/2
t_slice = 0.0

allocate(&
         &glesst(1:Nt, 1:Nt), grett(1:Nt, 1:Nt), ggreatt(1:Nt, 1:Nt),&
         &time(1:Nt), time2(1:Nt2), time3(1:Nt3), w(1:Nw), fw(1:Nw))

allocate(&
         &glesst2(1:Nt2, 1:Nt2), ggreatt2(1:Nt2, 1:Nt2), grett2(1:Nt2, 1:Nt2)&
         &, glesst3(1:Nt3, 1:Nt3), ggreatt3(1:Nt3, 1:Nt3), grett3(1:Nt3, 1:Nt3))

allocate(GLrel(1:2*Nt, 1:2*Nt), GRrel(1:2*Nt, 1:2*Nt), GGrel(1:2*Nt, 1:2*Nt),&
         &GLrel2(1:2*Nt2, 1:2*Nt2), GRrel2(1:2*Nt2, 1:2*Nt2), GGrel2(1:2*Nt2, 1:2*Nt2),&
         &GLrel3(1:2*Nt3, 1:2*Nt3), GRrel3(1:2*Nt3, 1:2*Nt3), GGrel3(1:2*Nt3, 1:2*Nt3),&
          &time_ave(1:2*Nt), time_rel(1:2*Nt),&
          &time_ave2(1:2*Nt2), time_rel2(1:2*Nt2),&
          &time_ave3(1:2*Nt3), time_rel3(1:2*Nt3),&
          &map(1:Ntd, 1:Ntd, 1:2), map_2ave(1:Nt, 1:Nt, 1:2))

! allocate(glessw(1:Nt,1:Nw), gretw(1:Nt,1:Nw), ggreatw(1:Nt,1:Nw))
allocate(glessw(1:2*Nt, 1:Nw), gretw(1:2*Nt, 1:Nw), ggreatw(1:2*Nt, 1:Nw),&
         &glessw2(1:2*Nt2, 1:Nw), gretw2(1:2*Nt2, 1:Nw), ggreatw2(1:2*Nt2, 1:Nw),&
         &glessw3(1:2*Nt3, 1:Nw), gretw3(1:2*Nt3, 1:Nw), ggreatw3(1:2*Nt3, 1:Nw),&
         &glesswI1(1:2*Nt3, 1:Nw),glesswI2(1:2*Nt3, 1:Nw), gretwI1(1:2*Nt3, 1:Nw),gretwI2(1:2*Nt3, 1:Nw),&
         &ggreatwI1(1:2*Nt3, 1:Nw),ggreatwI2(1:2*Nt3, 1:Nw))

allocate(temparray1(1:2*Nt), temparray2(1:2*Nt2), temparrayx(1:2*Nt3),&
         &fw1(1:Nw), fw2(1:Nw), fw3(1:Nw))!, fw1(1:Nw), fw2(1:Nw), fw3(1:Nw))

do i = 1, Nw

    w(i) = wmin + (i-0.5)*(-2*wmin/real(Nw))

!     write(*,*) i, w(i)

end do

deltat = real(tmax - tmin)/real(Nt)
deltat2 = real(tmax - tmin)/real(Nt2)
deltat3 = real(tmax - tmin)/real(Nt3)

write(*,*) 'deltat, deltat2, deltat3,', deltat, deltat2, deltat3

L0 = deltat2*deltat3/((deltat - deltat2)*(deltat - deltat3))
L1 = deltat*deltat3/((deltat2 - deltat)*(deltat2 - deltat3))
L2 = deltat*deltat2/((deltat3 - deltat)*(deltat3 - deltat2))

write(*,*) 'extrapolation weights,', L0, L1, L2

do i = 1, Nt
    time(i) = tmin + (i-0.5)*deltat
end do
do i = 1, Nt2
    time2(i) = tmin + (i-0.5)*deltat2
end do
do i = 1, Nt3
    time3(i) = tmin + (i-0.5)*deltat3
end do
do i = 1, 2*Nt
    time_ave(i) = tmin + (i)*0.5*deltat
    time_rel(i) = 2.0*tmin + (i)*deltat
end do
do i = 1, 2*Nt2
    time_ave2(i) = tmin + (i)*0.5*deltat2
    time_rel2(i) = 2.0*tmin + (i)*deltat2
end do
do i = 1, 2*Nt3
    time_ave3(i) = tmin + (i)*0.5*deltat3
    time_rel3(i) = 2.0*tmin + (i)*deltat3
end do

! open(unit=16,file='data/params/t_ave.dat',status='old')
!     read(16,*) t_slice
! close(unit=16)

open(unit=16,file='data/params/tave.dat',status='unknown')
    do i1 = 1, 2*Nt3
        write(16,*) i1, time_ave3(i1)
    end do
close(unit=16)

open(unit=16,file='data/params/omega.dat',status='unknown')
    do i1 = 1, Nw
        write(16,*) i1, w(i1)
    end do
close(unit=16)

! Prepare the slice integer from t_ave selection
do i = 1, Nt3
    if(abs(time3(i) - t_slice) .lt. deltat3)then
        i_slice = i
        exit
    end if
end do

i_slice = Nt3/2
! i_slice_ave = Nt3 - 1

! i_slice = i_slice

i_slice_ave = nint(((time3(i_slice) - tmin)/(0.5*deltat3)))

write(*,*) i_slice, i_slice_ave, time3(i_slice), time_ave3(i_slice_ave)


!

!import glesser, ggreater, gretarded for Nt, Nt2, Nt3

! open(unit=16,file='data/gr1.dat',status='old')
open(unit=16,file='data/g_ret1.dat',status='old')
do i1 = 1, Nt
    do i2 = 1, Nt
    read(16, *) j,k, rtemp, itemp
        if(j .ge. k)grett(j,k) =  ii*itemp! 
    end do
end do
close(unit=16)

! open(unit=16,file='data/gl1.dat',status='old')
open(unit=16,file='data/g_less1.dat',status='old')
do i1 = 1, Nt
    do i2 = 1, Nt
    read(16, *) j,k, rtemp, itemp
     glesst(j,k) = rtemp + ii*itemp! 
    end do
end do
close(unit=16)

! open(unit=16,file='data/gg1.dat',status='old')
open(unit=16,file='data/g_greater1.dat',status='old')
do i1 = 1, Nt
    do i2 = 1, Nt
    read(16, *) j,k, rtemp, itemp
     ggreatt(j,k) = rtemp + ii*itemp! 
    end do
end do
close(unit=16)

! open(unit=16,file='data/gr2.dat',status='old')
open(unit=16,file='data/g_ret2.dat',status='old')
do i1 = 1, Nt2
    do i2 = 1, Nt2
    read(16, *) j,k, rtemp, itemp
        if(j .ge. k)grett2(j,k) =  ii*itemp!  
    end do
end do
close(unit=16)

! open(unit=16,file='data/gl2.dat',status='old')
open(unit=16,file='data/g_less2.dat',status='old')
do i1 = 1, Nt2
    do i2 = 1, Nt2
    read(16, *) j,k, rtemp, itemp
     glesst2(j,k) = rtemp + ii*itemp! 
    end do
end do
close(unit=16)

! open(unit=16,file='data/gg2.dat',status='old')
open(unit=16,file='data/g_greater2.dat',status='old')
do i1 = 1, Nt2
    do i2 = 1, Nt2
    read(16, *) j,k, rtemp, itemp
     ggreatt2(j,k) = rtemp + ii*itemp! 
    end do
end do
close(unit=16)

! open(unit=16,file='data/gr3.dat',status='old')
open(unit=16,file='data/g_ret3.dat',status='old')
do i1 = 1, Nt3
    do i2 = 1, Nt3
    read(16, *) j,k, rtemp, itemp
     if(j .ge. k)grett3(j,k) =  ii*itemp! 
    end do
end do
close(unit=16)

! open(unit=16,file='data/gl3.dat',status='old')
open(unit=16,file='data/g_less3.dat',status='old')
do i1 = 1, Nt3
    do i2 = 1, Nt3
    read(16, *) j,k, rtemp, itemp
     glesst3(j,k) = rtemp + ii*itemp! 
    end do
end do
close(unit=16)

! open(unit=16,file='data/gg3.dat',status='old')
open(unit=16,file='data/g_greater3.dat',status='old')
do i1 = 1, Nt3
    do i2 = 1, Nt3
    read(16, *) j,k, rtemp, itemp
     ggreatt3(j,k) = rtemp + ii*itemp! 
    end do
end do
close(unit=16)

! call Wigner transform here
call Wignergridtransform(Nt, time, deltat, glesst, GLrel)
call Wignergridtransform(Nt, time, deltat, ggreatt, GGrel)
call Wignergridtransform(Nt, time, deltat, grett, GRrel)

call Wignergridtransform(Nt2, time2, deltat2, glesst2, GLrel2)
call Wignergridtransform(Nt2, time2, deltat2, ggreatt2, GGrel2)
call Wignergridtransform(Nt2, time2, deltat2, grett2, GRrel2)

call Wignergridtransform(Nt3, time3, deltat3, glesst3, GLrel3)
call Wignergridtransform(Nt3, time3, deltat3, ggreatt3, GGrel3)
call Wignergridtransform(Nt3, time3, deltat3, grett3, GRrel3)

open(unit=16,file='data/grettatrTEST.dat',status='unknown')
i1 = Nt3
do i_ave = 1, 2*Nt3

    write(16,*) time_ave3(i_ave), real(gRrel3(i_ave, i1)), aimag(gRrel3(i_ave,i1))!

end do
close(unit=16)

write(*,*) 'Green''s functions transformed to Wigner time'
if(.false.)then
! open(unit=16,file='data/glesstrel1.dat',status='unknown')
! do i1 = 1, 2*Nt
!     write(16,*) time_rel(i1), real(GLrel(Nt, i1)), aimag(GLrel(Nt, i1))
! end do
! close(unit=16)
! 
! open(unit=16,file='data/ggreattrel1.dat',status='unknown')
! do i1 = 1, 2*Nt
!     write(16,*) time_rel(i1), real(GGrel(Nt, i1)), aimag(GGrel(Nt, i1))
! end do
! close(unit=16)
! 
! open(unit=16,file='data/grettrel1.dat',status='unknown')
! do i1 = 1, 2*Nt
!     write(16,*) time_rel(i1), real(GRrel(Nt, i1)), aimag(GRrel(Nt, i1))
! end do
! close(unit=16)
! 
! open(unit=16,file='data/glesstrel2.dat',status='unknown')
! do i1 = 1, 2*Nt2
!     write(16,*) time_rel2(i1), real(GLrel2(Nt2, i1)), aimag(GLrel2(Nt2, i1))
! end do
! close(unit=16)
! 
! open(unit=16,file='data/ggreattrel2.dat',status='unknown')
! do i1 = 1, 2*Nt2
!     write(16,*) time_rel2(i1), real(GGrel2(Nt2, i1)), aimag(GGrel2(Nt2, i1))
! end do
! close(unit=16)
! 
! open(unit=16,file='data/grettrel2.dat',status='unknown')
! do i1 = 1, 2*Nt2
!     write(16,*) time_rel2(i1), real(GRrel2(Nt2, i1)), aimag(GRrel2(Nt2, i1))
! end do
! close(unit=16)
! 
! open(unit=16,file='data/glesstrel3.dat',status='unknown')
! do i1 = 1, 2*Nt3
!     write(16,*) time_rel3(i1), real(GLrel3(Nt3, i1)), aimag(GLrel3(Nt3, i1))
! end do
! close(unit=16)
! 
! open(unit=16,file='data/ggreattrel3.dat',status='unknown')
! do i1 = 1, 2*Nt3
!     write(16,*) time_rel3(i1), real(GGrel3(Nt3, i1)), aimag(GGrel3(Nt3, i1))
! end do
! close(unit=16)
! 
! open(unit=16,file='data/grettrel3.dat',status='unknown')
! do i1 = 1, 2*Nt3
!     write(16,*) time_rel3(i1), real(GRrel3(Nt3, i1)), aimag(GRrel3(Nt3, i1))
! end do
! close(unit=16)
! stop
end if
i_offset = int(time(1)/deltat)

glessw = 0
gretw = 0
ggreatw = 0

do i_ave = 1, 2*Nt

    call Ftransform1D(i_ave, 2*Nt, time_rel, deltat, Nw, w, GLrel(i_ave, :), fw)
    glessw(i_ave, :) = fw(:)

    call Ftransform1D(i_ave, 2*Nt, time_rel, deltat, Nw, w, GGrel(i_ave, :), fw)
    ggreatw(i_ave, :) = fw(:)
    
    call Ftransform1D(i_ave, 2*Nt, time_rel, deltat, Nw, w, GRrel(i_ave, :), fw)
    gretw(i_ave, :) = fw(:)

end do

do i_ave = 1, 2*Nt2

    call Ftransform1D(i_ave, 2*Nt2, time_rel2, deltat2, Nw, w, GLrel2(i_ave, :), fw)
    glessw2(i_ave, :) = fw(:)

    call Ftransform1D(i_ave, 2*Nt2, time_rel2, deltat2, Nw, w, GGrel2(i_ave, :), fw)
    ggreatw2(i_ave, :) = fw(:)
    
    call Ftransform1D(i_ave, 2*Nt2, time_rel2, deltat2, Nw, w, GRrel2(i_ave, :), fw)
    gretw2(i_ave, :) = fw(:)

end do

do i_ave = 1, 2*Nt3

    call Ftransform1D(i_ave, 2*Nt3, time_rel3, deltat3, Nw, w, GLrel3(i_ave, :), fw)
    glessw3(i_ave, :) = fw(:)

    call Ftransform1D(i_ave, 2*Nt3, time_rel3, deltat3, Nw, w, GGrel3(i_ave, :), fw)
    ggreatw3(i_ave, :) = fw(:)
    
    call Ftransform1D(i_ave, 2*Nt3, time_rel3, deltat3, Nw, w, GRrel3(i_ave, :), fw)
    gretw3(i_ave, :) = fw(:)

end do
i1 = Nw/2
open(unit=16,file='data/grettwTEST5.dat',status='unknown')

do i_ave = 1, 2*Nt3

    write(16,*) time_ave3(i_ave), real(gretw3(i_ave, i1)), aimag(gretw3(i_ave,i1))!

end do
close(unit=16)

open(unit=16,file='data/glesstwTEST5.dat',status='unknown')

do i_ave = 1, 2*Nt3

    write(16,*) time_ave3(i_ave), real(glessw3(i_ave, i1)), aimag(glessw3(i_ave,i1))!

end do
close(unit=16)

open(unit=16,file='data/ggreattwTEST5.dat',status='unknown')

do i_ave = 1, 2*Nt3

    write(16,*) time_ave3(i_ave), real(ggreatw3(i_ave, i1)), aimag(ggreatw3(i_ave,i1))!

end do
close(unit=16)


write(*,*) 'Fourier transforms in relative time complete'

!for each frequency, interpolate in average time
! 
do j = 1, Nw
   
        temparray1(:) = glessw(:, j)
        call interpolatequad1d(0.5*deltat, time_ave, time_ave3, temparray1, temparrayX, 2*Nt, 2*Nt3)
        glesswI1(:, j) = temparrayX(:)

        temparray2(:) = glessw2(:, j)
        call interpolatequad1d(0.5*deltat2, time_ave2, time_ave3, temparray2, temparrayX, 2*Nt2, 2*Nt3)
        glesswI2(:, j) = temparrayX(:)

        temparray1(:) = gretw(:, j)
        call interpolatequad1d(0.5*deltat, time_ave, time_ave3, temparray1, temparrayX, 2*Nt, 2*Nt3)
        gretwI1(:, j) = temparrayX(:)

        temparray2(:) = gretw2(:, j)
        call interpolatequad1d(0.5*deltat2, time_ave2, time_ave3, temparray2, temparrayX, 2*Nt2, 2*Nt3)
        gretwI2(:, j) = temparrayX(:)

        temparray1(:) = ggreatw(:, j)
        call interpolatequad1d(0.5*deltat, time_ave, time_ave3, temparray1, temparrayX, 2*Nt, 2*Nt3)
        ggreatwI1(:, j) = temparrayX(:)

        temparray2(:) = ggreatw2(:, j)
        call interpolatequad1d(0.5*deltat2, time_ave2, time_ave3, temparray2, temparrayX, 2*Nt2, 2*Nt3)
        ggreatwI2(:, j) = temparrayX(:)

end do

write(*,*) 'Interpolation in average time complete'





if(.false.)then !old Fourier transform code
do i = 1, 2*Nt
    do j = 1, Nw
        do i_rel = 1, 2*Nt

            scatemp = exp(ii*w(j)*time_rel(i_rel))*deltat*GLrel(i, i_rel)
            glessw(i,j) = glessw(i,j) + scatemp
            
            scatemp = exp(ii*w(j)*time_rel(i_rel))*deltat*GGrel(i, i_rel)
            ggreatw(i,j) = ggreatw(i,j) + scatemp
            
            scatemp = exp(ii*w(j)*time_rel(i_rel))*deltat*GRrel(i, i_rel)
            gretw(i,j) = gretw(i,j) + scatemp
            
        end do
    end do
end do
end if

write(*,*) 't_ave value exported:', time_ave3(i_slice_ave)
write(*,*) 'center t_ave value:', time3(Nt3/2), time3(Nt3/2) + 2.5
write(*,*) 'difference:', abs(time3(Nt3/2) - time3(i_slice))
call cpu_time(T2)

write(*,*) 'cpu times T1 and T2 ', T1, T2
write(*,*) 'elapsed time in minutes', (T2 - T1)/60

do i1 = 1, Nw

    fw1(i1) = -0.25*(aimag(glesswI1(i_slice_ave, i1)) - aimag(ggreatwI1(i_slice_ave, Nw-i1+1)))&!
              &/aimag(gretwI1(i_slice_ave, i1))
                              
    fw2(i1) = -0.25*(aimag(glesswI2(i_slice_ave, i1)) - aimag(ggreatwI2(i_slice_ave, Nw-i1+1)))&!
                &/aimag(gretwI2(i_slice_ave, i1))
!                 
    fw3(i1) = -0.25*(aimag(glessw3(i_slice_ave, i1)) - aimag(ggreatw3(i_slice_ave, Nw-i1+1)))&!
                &/aimag(gretw3(i_slice_ave, i1))!fw3
                              
end do

fw = L0*fw1 + L1*fw2 + L2*fw3

open(unit=16,file='data/fw.dat',status='unknown')
    do i1 = Nw/4, Nw-Nw/4
        write(16,*) w(i1), real(fw(i1))
    end do
close(unit=16)
open(unit=16,file='data/fw1.dat',status='unknown')
    do i1 = Nw/4, Nw-Nw/4
        write(16,*) w(i1), real(fw1(i1))
    end do
close(unit=16)
open(unit=16,file='data/fw2.dat',status='unknown')
    do i1 = Nw/4, Nw-Nw/4
        write(16,*) w(i1), real(fw2(i1))
    end do
close(unit=16)
open(unit=16,file='data/fw3.dat',status='unknown')
    do i1 = Nw/4, Nw-Nw/4
        write(16,*) w(i1), real(fw3(i1))
    end do
close(unit=16)

gretw3 = L0*gretwI1 + L1*gretwI2 + L2*gretw3
glessw3 = L0*glesswI1 + L1*glesswI2 + L2*glessw3
ggreatw3 = L0*ggreatwI1 + L1*ggreatwI2 + L2*ggreatw3

do i1 = 1, Nw

    fw(i1) = -0.25*(aimag(glessw3(i_slice_ave, i1)) - aimag(ggreatw3(i_slice_ave, Nw-i1+1)))&!
                &/aimag(gretw3(i_slice_ave, i1))!fw3
                              
end do

open(unit=16,file='data/fwalt.dat',status='unknown')
    do i1 = Nw/4, Nw-Nw/4
        write(16,*) w(i1), real(fw(i1))
    end do
close(unit=16)
! 
! open(unit=16,file='ggreatw.dat',status='unknown')
! do i1 = 1, Nw
!     write(16,*) w(i1), real(ggreatw(Nw - i1 + 1)), -aimag(ggreatw(Nw +1 - i1))
! end do
! close(unit=16)
! ! 
open(unit=16,file='data/ggreatw.dat',status='unknown')
do i_ave = 1, 2*Nt3
do i1 = 1, Nw
    write(16,*) i_ave, i1, real(ggreatw3(i_ave, i1)), aimag(ggreatw3(i_ave,i1))!
end do
end do
close(unit=16)
! 
! i1 = Nw/2
! open(unit=16,file='data/ggreatwTEST1.dat',status='unknown')
! 
! do i_ave = 1, 2*Nt3
! 
!     write(16,*) time_ave3(i_ave), real(ggreatwI1(i_ave, i1)), aimag(ggreatwI1(i_ave,i1))!
! 
! end do
! close(unit=16)
! ! i1 = Nw - Nw/4
! open(unit=16,file='data/ggreatwTEST2.dat',status='unknown')
! 
! do i_ave = 1, 2*Nt3
! 
!     write(16,*) time_ave3(i_ave), real(ggreatwI2(i_ave, i1)), aimag(ggreatwI2(i_ave,i1))!
! 
! end do
! close(unit=16)
! 
! open(unit=16,file='data/ggreatwTEST3.dat',status='unknown')
! 
! do i_ave = 1, 2*Nt3
! 
!     write(16,*) time_ave3(i_ave), real(ggreatw3(i_ave, i1)), aimag(ggreatw3(i_ave,i1))!
! 
! end do
! close(unit=16)

! ggreatw3 = L0*ggreatwI1 + L1*ggreatwI2 + L2*ggreatw3
! 
! open(unit=16,file='data/ggreatwTEST4.dat',status='unknown')
! 
! do i_ave = 1, 2*Nt3
! 
!     write(16,*) time_ave3(i_ave), real(ggreatw3(i_ave, i1)), aimag(ggreatw3(i_ave,i1))!
! 
! end do
! 
! open(unit=16,file='data/glesswTEST1.dat',status='unknown')
! 
! do i_ave = 1, 2*Nt3
! 
!     write(16,*) time_ave3(i_ave), real(glesswI1(i_ave, i1)), aimag(glesswI1(i_ave,i1))!
! 
! end do
! close(unit=16)
! ! i1 = Nw - Nw/4
! open(unit=16,file='data/glesswTEST2.dat',status='unknown')
! 
! do i_ave = 1, 2*Nt3
! 
!     write(16,*) time_ave3(i_ave), real(glesswI2(i_ave, i1)), aimag(glesswI2(i_ave,i1))!
! 
! end do
! close(unit=16)
! 
! open(unit=16,file='data/glesstwTEST3.dat',status='unknown')
! 
! do i_ave = 1, 2*Nt3
! 
!     write(16,*) time_ave3(i_ave), real(glessw3(i_ave, i1)), aimag(glessw3(i_ave,i1))!
! 
! end do
! close(unit=16)
! 
! ! glessw3 = L0*glesswI1 + L1*glesswI2 + L2*glessw3
! 
! open(unit=16,file='data/glesstwTEST4.dat',status='unknown')
! 
! do i_ave = 1, 2*Nt3
! 
!     write(16,*) time_ave3(i_ave), real(glessw3(i_ave, i1)), aimag(glessw3(i_ave,i1))!
! 
! end do
! close(unit=16)
! 
! open(unit=16,file='data/gretwTEST1.dat',status='unknown')
! 
! do i_ave = 1, 2*Nt3
! 
!     write(16,*) time_ave3(i_ave), real(gretwI1(i_ave, i1)), aimag(gretwI1(i_ave,i1))!
! 
! end do
! close(unit=16)
! ! i1 = Nw - Nw/4
! open(unit=16,file='data/gretwTEST2.dat',status='unknown')
! 
! do i_ave = 1, 2*Nt3
! 
!     write(16,*) time_ave3(i_ave), real(gretwI2(i_ave, i1)), aimag(gretwI2(i_ave,i1))!
! 
! end do
! close(unit=16)
! 
! open(unit=16,file='data/grettwTEST3.dat',status='unknown')
! 
! do i_ave = 1, 2*Nt3
! 
!     write(16,*) time_ave3(i_ave), real(gretw3(i_ave, i1)), aimag(gretw3(i_ave,i1))!
! 
! end do
! close(unit=16)
! 
! ! gretw3 = L0*gretwI1 + L1*gretwI2 + L2*gretw3
! 
! open(unit=16,file='data/grettwTEST4.dat',status='unknown')
! 
! do i_ave = 1, 2*Nt3
! 
!     write(16,*) time_ave3(i_ave), real(gretw3(i_ave, i1)), aimag(gretw3(i_ave,i1))!
! 
! end do
! close(unit=16)
! close(unit=16)

! ! ! 
open(unit=16,file='data/glessw.dat',status='unknown')
do i_ave = 1, 2*Nt3
do i1 = 1, Nw
    write(16,*) i_ave, i1, real(glessw3(i_ave,i1)), aimag(glessw3(i_ave,i1))!time_ave3(i), 
end do
end do
close(unit=16)

open(unit=16,file='data/gretw.dat',status='unknown')
do i_ave =1, 2*Nt3
do i1 = 1, Nw
    write(16,*) i_ave, i1, real(gretw3(i_ave,i1)), aimag(gretw3(i_ave,i1))!time_ave3(i), 
end do
end do
close(unit=16)

open(unit=16, file='data/params/params.dat', status='unknown')

    write(16,*) Nt3, Nw
    write(16,*) 'Nt3,','Nw'

close(unit=16)

! open(unit=16,file='data/gretw2.dat',status='unknown')
! ! do i =1, 2*Nt3
! do i1 = 1, Nw
!     write(16,*) w(i1), real(gretwI2(i_slice_ave,i1)), aimag(gretwI2(i_slice_ave,i1))!time_ave3(i), 
! end do
! ! end do
! close(unit=16)
! 
! open(unit=16,file='data/gretw3.dat',status='unknown')
! ! do i =1, 2*Nt3
! do i1 = 1, Nw
!     write(16,*) w(i1), real(gretw3(i_slice_ave,i1)), aimag(gretw3(i_slice_ave,i1))!time_ave3(i), 
! end do
! ! end do
! close(unit=16)

contains

! subroutine Ftransform1D(i_ave, Nt_in, trel_in, deltat_in, Nw, w, input, output)
! complex(kind=wp), intent(in), dimension(:,:) :: input
! complex(kind=wp), intent(out), dimension(:) :: output
! real(kind=wp), intent(in), dimension(:) :: trel_in, w
! real(kind=wp), intent(in) :: deltat_in
! integer, intent(in) :: Nt_in, Nw
! complex(kind=wp) :: scatemp
! integer :: i, j, i_rel, i_ave
! 
! output = 0
! scatemp = 0
! 
! write(*,*) deltat_in
! 
! do i_ave = 1, Nt_in
!     do j = 1, Nw
!         do i_rel = 1, Nt_in
! 
!             scatemp = exp(ii*w(j)*trel_in(i_rel))*deltat_in*input(i_ave, i_rel)
!             output(i_ave,j) = output(i_ave, j) + scatemp
!             
!         end do
!     end do
! end do
! 
! 
! end subroutine Ftransform1D

subroutine Ftransform1D(i_ave, Nt_in, trel_in, deltat_in, Nw, w, input, output)
complex(kind=wp), intent(in), dimension(:) :: input
complex(kind=wp), intent(out), dimension(:) :: output
real(kind=wp), intent(in), dimension(:) :: trel_in, w
real(kind=wp), intent(in) :: deltat_in
integer, intent(in) :: Nt_in, Nw, i_ave
complex(kind=wp) :: scatemp
integer :: i, j, i_rel

output = 0
scatemp = 0

! write(*,*) deltat_in

! do i_ave = 1, Nt_in
    do j = 1, Nw
        do i_rel = 1, Nt_in - 1

            scatemp = exp(ii*w(j)*trel_in(i_rel))*deltat_in*input(i_rel)!i_ave, 
            output(j) = output(j) + scatemp
            
        end do
    end do
! end do


end subroutine Ftransform1D

subroutine Wignergridtransform(Nt_in, t_in, deltat, input, output)

complex(kind=wp), intent(in), dimension(:,:) :: input
real(kind=wp), intent(in), dimension(:) :: t_in
real(kind=wp), intent(in) :: deltat
integer, intent(in):: Nt_in
complex(kind=wp), intent(out), dimension(:,:) :: output
integer :: k0, k1, k2, k3, i, j, k, Ntd, i_ave, i_rel, i1, i2
integer, allocatable, dimension(:,:,:) :: map!, map_2ave
real(kind=wp) :: t_min, t_max,  xx, yy
real(kind=wp), allocatable, dimension(:) :: t_ave, t_rel

t_min = t_in(1)
t_max = t_in(Nt_in)

! write(*,*) 'testing...', t_min, t_max

Ntd = 2*Nt_in

allocate(t_ave(1:Ntd), t_rel(1:Ntd))

do i = 1, Ntd
    t_ave(i) = t_min + (i)*0.5*deltat
    t_rel(i) = 2.0*t_min + (i)*deltat
end do

allocate(map(1:Ntd, 1:Ntd, 1:2))!, map_2ave(1:Nt_in, 1:Nt_in, 1:2))

map = -1
map_2ave = -1

do i_ave=1, Ntd
    do i_rel=1, Ntd
        xx = t_ave(i_ave) + (0.5*t_rel(i_rel)) !t
        yy = t_ave(i_ave) - (0.5*t_rel(i_rel)) !t'
        
        i1 = nint(((xx-t_min)/(deltat)) + 1.0)
        i2 = nint(((yy-t_min)/(deltat)) + 1.0)
        
        
        if(i1 .gt. 0 .and. i1 .le. Nt_in .and. i2 .gt. 0 .and. i2 .le. Nt_in)then
        
            
!         write(*,*) xx-time(i1), yy-time(i2)
        
            map(i_ave, i_rel, 1) = i1
            map(i_ave, i_rel, 2) = i2
!             
!             map_2ave(i1, i2, 1) = i_ave
!             map_2ave(i1, i2, 2) = i_rel
        
        end if
        
    end do
end do

output = 0

do i_ave = 1, Ntd
    do i_rel = 1, Ntd

        i1 = map(i_ave, i_rel, 1)
        i2 = map(i_ave, i_rel, 2)
        
        if(i1 .gt. 0 .and. i2 .gt. 0)then
        
            output(i_ave, i_rel) = input(i1, i2)
            
        end if

    end do
end do

end subroutine Wignergridtransform

subroutine interpolatequad1d(deltatA, tA, tB, input, output, NtA, NtB)

complex(kind=wp), intent(in), dimension(:) :: input
real(kind=wp), intent(in), dimension(:) :: tA, tB
real(kind=wp), intent(in) :: deltatA
integer, intent(in):: NtA, NtB
complex(kind=wp), intent(out), dimension(:) :: output

integer :: k0, k1, k2, k3, i, j, k
output = 0
do i = 2, NtB

!find nearest point
    do j = 1, NtA
        if(abs(tA(j) - tB(i)) .lt. deltatA)then
            k = j !k is now index of nearest point
!             write(*,*) i,j,k
!             if(k .ge. Nt1 - 2) write(*,*) 'error!', k, i
            exit
        end if
    end do
    
    
    
    k0 = k
    k1 = k + 1
    k2 = k + 2
    
    if(k .ge. NtA-2)then
    
        k0 = NtA-2
        k1 = NtA-1
        K2 = NtA
    
    end if
    
    if(k0 .eq. 0) write(*,*) 'k0'
    
    L0 = (tB(i) - tA(k1))*(tB(i) - tA(k2))/((tA(k0) - tA(k1))*(tA(k0) - tA(k2)))
    L1 = (tB(i) - tA(k0))*(tB(i) - tA(k2))/((tA(k1) - tA(k0))*(tA(k1) - tA(k2)))
    L2 = (tB(i) - tA(k0))*(tB(i) - tA(k1))/((tA(k2) - tA(k0))*(tA(k2) - tA(k1)))
    
    output(i) = input(k0)*L0 + input(k1)*L1 + input(k2)*L2
    
end do
end subroutine interpolatequad1d

end program processGLGRGG
