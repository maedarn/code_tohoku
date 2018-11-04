program err
  implicit none
  integer :: ndx , i , k,loop
  DOUBLE PRECISION, dimension(:), allocatable :: Phi, Phiexa , gosa , dx,Phi2,gosaa
  character(1) num
  double precision a,gosa1,Lbox,c,d,e,gosa2,ll
  real(4) b

  !allocate(Phi(1:ndx-2))
  !allocate(Phiexa(1:ndx-2))
  loop=5
  Lbox=1.0d2
  allocate(gosa(1:loop))
  allocate(gosaa(1:loop))
  allocate(dx(1:loop))
  gosa(:) = 0.0d0
  ndx=34

  do k=1,loop
     ndx=(ndx-2)*2 + 2
     write(*,*) ndx
     dx(k)=dble(ndx)
     allocate(Phi(1:ndx-2))
     allocate(Phi2(1:ndx-2))
     allocate(Phiexa(-1:ndx))
     write(num,'(i1.1)') k
     open(18,file='/Users/maeda/Desktop/kaiseki/testcode6/phi'//num//'00000.dat')
     open(19,file='/Users/maeda/Desktop/kaiseki/testcode6/phi'//num//'exact.DAT')
     !if(k==5) then
     !   close(18)
     !   open(18,file='/Users/maeda/Desktop/kaiseki/testcode3/phi24999.dat')
     !end if
     gosa1=0.0d0
     gosa2=0.0d0
     do i=1,ndx-2
        read(18,*) a , Phi(i),c,d,e,Phi2(i)
     end do
     close(18)
     do i=-1,ndx
        read(19,*) b,Phiexa(i)
     end do
     close(19)

     do i=1,ndx-2
        gosa1 = (Phi(i)-Phiexa(i)) ** 2 + gosa1
     end do
     do i=1,ndx-2
        gosa2 = (Phi2(i)-Phiexa(i)) ** 2 + gosa2
     end do
     if(loop==5) then
        open(30,file='/Users/maeda/Desktop/kaiseki/testcode6/katamuki.dat')
        ll=Lbox/(dx(5)-2.0d0)/2.0d0
        do i=2,ndx-3
           !if(i==2) then
           !   ll=ll+Lbox/(dx(5)-2.0d0)
           !   write(30,*) ll, (-Phi2(i)+Phi2(i+1))/2.0d0/(Lbox/(dx(5)-2.0d0)),-Phi2(i),Phi2(i+1),&
           !        (-Phi2(i)+Phiexa(i-1)+Phi2(i+1)-Phiexa(i+1))/2.0d0/(Lbox/(dx(5)-2.0d0)),&
           !        Phi2(i)-Phiexa(i),-Phi2(i)+Phi2(i+1)!-Phiexa(i)
           !elseif(i== ndx-3) then
           !   ll=ll+Lbox/(dx(5)-2.0d0)
           !   write(30,*) ll, (-Phi2(i-1)+Phi2(i))/2.0d0/(Lbox/(dx(5)-2.0d0)),-Phi2(i-1),Phi2(i),&
           !        (-Phi2(i-1)+Phiexa(i-1)+Phi2(i)-Phiexa(i+1))/2.0d0/(Lbox/(dx(5)-2.0d0)),&
           !        Phi2(i)-Phiexa(i),-Phi2(i-1)+Phi2(i)!-Phiexa(i)
              !ll=ll+Lbox/(dx(5)-2.0d0)
           !else
              ll=ll+Lbox/(dx(5)-2.0d0)
              write(30,*) ll, (-Phi2(i-1)+Phi2(i+1))/2.0d0/(Lbox/(dx(5)-2.0d0)),-Phi2(i-1),Phi2(i+1),&
                   (-Phi2(i-1)+Phiexa(i-1)+Phi2(i+1)-Phiexa(i+1))/2.0d0/(Lbox/(dx(5)-2.0d0)),&
                   Phi2(i)-Phiexa(i),-Phi2(i-1)+Phi2(i+1)!-Phiexa(i)
              !ll=ll+Lbox/(dx(5)-2.0d0)
           !end if
        end do
        close(30)
     end if

     gosa1=dsqrt(gosa1)
     gosa2=dsqrt(gosa2)
     gosa1=gosa1/dble(ndx-2)
     gosa2=gosa2/dble(ndx-2)
     gosa(k)=gosa1
     gosaa(k)=gosa2
     deallocate(Phi,Phiexa,Phi2)
  end do


  open(20,file='/Users/maeda/Desktop/kaiseki/testcode6/gosa.dat')
  do i=1,loop
     write(20,*) Lbox/(dx(i)-2.0d0) , gosa(i) , gosaa(i)
  end do
  close(20)
  deallocate(dx,gosa,gosaa)
end program err
