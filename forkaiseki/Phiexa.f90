program main
  implicit none
  integer i,j,k,nx,ny,nz,core,n
  double precision , allocatable , dimension(:,:,:) :: Phiexa
  real(4) , allocatable , dimension(:,:,:) :: Phi , dphi
  real(4) , allocatable , dimension(:) :: x,y,z
  real(8) , allocatable , dimension(:,:,:,:) :: NPhi
  DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G
  double precision censh ,Hsheet ,rho , minexa
  real(4) minphi,sa
  character(3) NPI

  !*************
  core=64
  nx=256
  ny=256
  nz=256
  censh = 50.0d0
  Hsheet =10.0d0
  rho = 2.0d0/G4pi/90.d0
  !*************
  allocate(Phi(1:nx,1:ny,1:nz))
  allocate(dphi(1:nx,1:ny,1:nz))
  allocate(Phiexa(1:nx,1:ny,1:nz))
  allocate(x(1:nx))
  allocate(y(1:ny))
  allocate(z(1:nz))
  allocate(Phi(1:nx+1,1:ny+1,1:nz+1,1:core))

  open(140,file='cdnt.DAT')
  do i=1,nx
     read(140,*)  x(i)
  end do
  close(140)

  open(141,file='final.DAT')
  do k=1,nz
     do j=1,ny
        do i=1,nx
           read(141,*)  Phi(i,j,k)
        end do
     end do
  end do
  close(141)

  do n=0,core-1
     write(NPI,'(i3.3)') n
     open(170+n,file='final'//NPI//'.DAT')
     do k=1,nz+1
        do j=1,ny+1
           do i=1,nx+1
              read(170+n,*)  NPhi(i,j,k,n+1)
           end do
        end do
     end do
     close(170+n)
  end do

  minphi=1.0d2
  do k=1,nz
     do j=1,ny
        do i=1,nx
          minphi = amin1(minphi,Phi(i,j,k))
        end do
     end do
  end do

  minexa = 1.0d2
  open(142,file='phiexact.DAT')
  do k=1,nz
     do j=1,ny
        do i= 1,nx
           if( dabs(dble(x(i)) - censh ) .le. Hsheet ) then
              Phiexa(i,j,k) = G4pi/2.0d0 * rho * (dble(x(i)) - censh )**2
           else
              Phiexa(i,j,k) = G4pi * rho * Hsheet * dabs(dble(x(i)) - censh ) - G4pi/2.0d0 * rho * Hsheet**2
           end if
           write(142,*) sngl(Phiexa(i,j,k))
           minexa=dmin1(minexa,Phiexa(i,j,k))
        end do
     end do
  end do
  close(142)

  sa=minphi - sngl(minexa)

  WRITE(*,*)  Phi(1,1,1) - sngl( Phiexa(1,1,1) ) , Phi(nx/2,ny/2,nz/2) - sngl( Phiexa(nx/2,ny/2,nz/2) ) , &
       Phi(nx/2+1,ny/2+1,nz/2+1) - sngl( Phiexa(nx/2+1,ny/2+1,nz/2+1) ),sa

  open(143,file='dphi.DAT')
  do k=1,nz
     do j=1,ny
        do i= 1,nx
           dphi(i,j,k) = Phi(i,j,k) - sngl(Phiexa(i,j,k)) - sa!( Phi(nx/2,ny/2,nz/2) - sngl( Phiexa(nx/2,ny/2,nz/2) ) ) !Phi(nx/2,ny/2,nz/2)
           write(143,*) sngl(dphi(i,j,k))
        end do
     end do
  end do
  close(143)


end program main
