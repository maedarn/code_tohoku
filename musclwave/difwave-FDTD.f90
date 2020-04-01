module comvar
  implicit none
  integer, parameter :: ndx=514,laststep=12000,ist=1,ien=2 !preiodic:ist=1,ien=2 , kotei:ist=2,ien=3 : ndx=130
  !double precision, parameter :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0
  DOUBLE PRECISION :: cg = 10.0d1 , dx != Lbox/dble(ndx-2) !, bcphi1 , bcphi2
  double precision :: Lbox=1.0d2 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0 , ch = 1.d0 , meanrho! , h=10.0d0
  !double precision :: G=1.11142d-4, G4pi=12.56637d0*G , coeff=0.90d0 ,  kappa=1.0d0/3.0d0
  double precision ::  G4pi=12.56637d0*1.11142d-4 , coeff=0.5d0 !,  kappa=1.0d0/3.0d0
  DOUBLE PRECISION , dimension(1:3) :: bcphi1 , bcphi2 ,bcphigrd1 , bcphigrd2
  character(54) :: dir='/Users/maeda/Desktop/Dropbox/kaiseki-desktpo/testwave/'
  integer :: bckd=3
end module comvar

module grvvar
  implicit none
  integer, parameter :: ndx2=514 !パラメータ属性必要
  DOUBLE PRECISION , dimension(-1:ndx2) :: x,rho
  DOUBLE PRECISION , dimension(-1:ndx2) :: Phidt,Phiexa,Phiv
  !DOUBLE PRECISION , dimension(-1:ndx2) :: Phidt
end module grvvar


program main
  use comvar
  use grvvar
  implicit none
  integer i,j,mode
  !integer , parameter :: ndx=514
  DOUBLE PRECISION, dimension(-1:ndx) :: Phidummy
  double precision nu2,dt,h
  integer :: sv=0,svmod=100

  mode=bckd

  call INIT()
  call time(dt)

  Phiv(:)=0.d0
  Phidummy(:)=0.d0

  !Phiv(0) = Phiexa(0)
  !Phiv(-1) = Phiexa(-1)
  !Phiv(ndx-1) = Phiexa(ndx-1)
  !Phiv(ndx) = Phiexa(ndx)
  Phiv(:) = Phiexa(:)

  nu2=dt*dt*cg*cg/dx/dx

  h=0.d0
  do j=1,laststep
     write(*,*)'step=', j
     h= ch*dt+h
     !call exa(h,1)
     call exa(h,mode)
     if((j==1).or.(j==2))then
        Phidummy(:)=Phiv(:)
        Phiv(:)=Phiexa(:)
        goto 300
     end if

     !call BC(1) !固定
     call BC(mode) !周期3

     Phidummy(:)=Phiv(:)
     do i=ist,ndx-ien
        Phiv(i)= 2.0d0*Phidummy(i) - Phidt(i) + nu2 * (Phidummy(i-1) + Phidummy(i+1)  &
             - w1 * Phidummy(i)) - rho(i) * G4pi * dt * dt * cg * cg
     end do

     300 continue
     Phidt(:)=Phidummy(:)
     if(mod(j,svmod)==1)then
        call saveu(sv)
     end if
  end do
end program main

subroutine INIT()
  use comvar
  use grvvar
  integer :: i
  double precision :: amp,pi=3.1415926535d0,haba

  dinit1 = 2.0d0/G4pi/90.d0

  !----------x--------------
  dx = Lbox/dble(ndx-2)
  x(1) = dx/2.0d0
  x(0) = x(1) - dx
  x(-1) = x(0) - dx
  do i=2,ndx
     x(i) = x(i-1) + dx
  end do
  !----------x--------------


  !-------Phidt-----------
  Phidt(:)=0.0d0
  !-------Phdt-----------

  !---------rho-------------
  rho(:)=0.d0
  goto 199
  do i = -1,ndx
     if( dabs(x(i) - hcen) .le. h) then
        rho(i) = dinit1
        !rho(i) = 0.0d0
     else
        rho(i) = 0.0d0
        !rho(i) = dinit1
        rho(i) = dinit1*1.d-2
     end if
  end do
  199 continue
  !---------rho-------------


  !--------Phiexa-----------
  goto 200
  open(142,file=dir//'phiexact.DAT')
  open(143,file=dir//'INIden.DAT')
  open(144,file=dir//'phigrd.DAT')
  do i= -1,ndx
     if( dabs(x(i) - hcen) .le. h ) then
        Phiexa(i) = G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2
        write(142,*) sngl(x(i)) ,  sngl(G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2)
     else
        Phiexa(i) = G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
        write(142,*) sngl(x(i)) , sngl(G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2)
     end if
     write(143,*) sngl(rho(i))
  end do
  200 continue
  !--------Phiexa-----------

end subroutine INIT

subroutine time(dt)
  use comvar
  use grvvar
  double precision :: dt
  dt = dx/cg * coeff
  write(*,*) 'time cg' , dt
end subroutine time

subroutine saveu(in1)
  use comvar
  use grvvar
  integer :: i,in1
  character(5) name

  write(name,'(i5.5)') in1
  open(21,file=dir//'phi'//name//'.dat')
  do i=-1,ndx
        write(21,*) x(i), Phiv(i) ,Phiexa(i) , rho(i)
  end do
  close(21)
  in1=in1+1
end subroutine saveu


     !---------rho-------------
     do i = -1,ndx
        if( dabs(x(i) - hcen) .le. h) then
           rho(i) = dinit1
           !rho(i) = 0.0d0
        else
           rho(i) = 0.0d0
           !rho(i) = dinit1
           !rho(i) = dinit1*1.d-2
        end if
     end do

     meanrho=0.d0
     do i = 1,ndx-2
        meanrho=meanrho+rho(i)
     end do
     meanrho=meanrho/dble(ndx-2)

     do i = -1,ndx
           rho(i) = rho(i)-meanrho
     end do
     !---------rho-------------
  end if


  !--------Phiexa-----------
  !goto 200
!  open(142,file=dir//'phiexact.DAT')
!  open(143,file=dir//'INIden.DAT')
!  open(144,file=dir//'phigrd.DAT')
  do i= -1,ndx
     if( dabs(x(i) - hcen) .le. h ) then
        Phiexa(i) = G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2
        write(142,*) sngl(x(i)) ,  sngl(G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2)
     else
        Phiexa(i) = G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
        write(142,*) sngl(x(i)) , sngl(G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2)
     end if
     write(143,*) sngl(rho(i))
  end do
  !--------Phiexa-----------
end subroutine exa

subroutine BC(mode)
  use comvar
  use grvvar
  integer :: mode

  if(mode==1)then
     Phiv(0) = Phiexa(0)
     Phiv(-1) = Phiexa(-1)
     Phiv(ndx-1) = Phiexa(ndx-1)
     Phiv(ndx) = Phiexa(ndx)
  end if

  if(mode==3)then
     Phiv(0) = Phiexa(ndx-2)
     Phiv(-1) = Phiexa(ndx-3)
     Phiv(ndx-1) = Phiexa(1)
     Phiv(ndx) = Phiexa(2)
  end if

end subroutine BC
