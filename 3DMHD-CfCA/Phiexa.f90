program main
  implicit none
  integer i,j,k,nx,ny,nz
  double precision , allocatable , dimension(:,:,:) :: Phiexa
  real(4) , allocatable , dimension(:,:,:) :: Phi , dphi
  real(4) , allocatable , dimension(:) :: x,y,z
  DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G
  double precision censh ,Hsheet ,rho

  !*************
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

  write(*,*) 'p1'

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
        end do
     end do
  end do
  close(142)

  write(*,*) 'p2'

  open(143,file='dphi.DAT')
  do k=1,nz
     do j=1,ny
        do i= 1,nx
           dphi(i,j,k) = Phi(i,j,k) - sngl(Phiexa(i,j,k)) - Phi(nx/2,ny/2,nz/2)
           write(143,*) sngl(dphi(i,j,k))
        end do
     end do
  end do
  close(143)

  write(*,*) 'p3'

end program main
