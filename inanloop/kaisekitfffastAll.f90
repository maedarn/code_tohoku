program main
  implicit none
  real(4),allocatable,dimension(:,:,:,:) :: U
  real(4),allocatable,dimension(:) :: shd
  integer i,j,k,nx,ny,nz,val,time,lstime,is,maxshd,timejump
  integer intime
  character(3) timech
  real(4) ,allocatable,dimension(:) :: ttt
  real(4) ,allocatable,dimension(:,:) :: tff,Mass,inttff!,intMass
  real(4)  Lbox,Mass1,tff1,dl,Msun,dt
  real(4), parameter :: G=1.11142e-4

  nx=512
  ny=512
  nz=512
  val=9
  lstime=178
  Lbox=100.e0
  Msun=1.473e-2
  maxshd=15
  timejump=4
  intime=2
  dt=1.e0

  dl=Lbox/real(nx)

  allocate(U(nx,ny,nz,val))
  allocate(tff(lstime,maxshd),Mass(lstime,maxshd),ttt(lstime),inttff(0:lstime,maxshd))
  allocate(shd(maxshd))
  tff(:,:)=0.e0
  Mass(:,:)=0.e0
  inttff(:,:)=0.e0
  !intMass(:,:)=0.e0
  ttt(:)=0.e0

  do i=1,maxshd
     shd(i)=1.e4*real(i)
  end do

  do time=intime,lstime,timejump
     write(*,*) time,'timestep'
     ttt(time)=0.25e0*real(time)
     write(timech,'(i3.3)') time
     open(110,file='Allnewbigtime'//timech//'.DAT',FORM='UNFORMATTED')
     U(:,:,:,:)=0.e0
     do k=1,nz
        do j=1,ny
           do i=1,nz
              read(110) U(i,j,k,1), U(i,j,k,2), U(i,j,k,3), U(i,j,k,4), U(i,j,k,5),&
                   U(i,j,k,6), U(i,j,k,7), U(i,j,k,8), U(i,j,k,9)
           end do
        end do
        write(*,*) k, U(5,5,k,2),'read-phyval'
     end do
     close(110)

     do is=1,maxshd
        tff1=0.e0
        Mass1=0.e0
        do k=1,nz
           do j=1,ny
              do i=1,nz
                 if(U(i,j,k,1)>shd(is))then
                    !tff1=1/dsqrt(G*dble(U(i,j,k,1)))+tff1
                    tff1=tff1+U(i,j,k,1)*sqrt(G*U(i,j,k,1))*dl*dl*dl
                    !cc=cc+1
                    Mass1=U(i,j,k,1)*dl*dl*dl+Mass1
                    !mxrho1=amax1(mxrho1,v(i,j,k,1))
                 end if
              end do
           end do
           write(*,*) k,'caltff'
        end do
        tff1=tff1*Msun
        Mass1=Mass1*Msun
        tff(time,is)=tff1
        Mass(time,is)=Mass1
        inttff(time,is)=tff1*dt+inttff(time-timejump,is) !timejump is important
        !intMass(time,is)=Mass1
     end do
  end do

  open(190,file='tffint.DAT',FORM='FORMATTED')
  do time=intime,lstime,timejump
     write(190,*) ttt(time),(inttff(time,is),is=1,maxshd)
  end do
  close(190)

  open(290,file='tffdt.DAT',FORM='FORMATTED')
  do time=intime,lstime,timejump
     write(290,*) ttt(time),(tff(time,is),is=1,maxshd)
  end do
  close(290)

  open(390,file='GasMasstff.DAT',FORM='FORMATTED')
  do time=intime,lstime,timejump
     write(390,*) ttt(time),(Mass(time,is),is=1,maxshd)
  end do
  close(390)

  deallocate(U)
  deallocate(tff,Mass,inttff,ttt,shd)
end program main
