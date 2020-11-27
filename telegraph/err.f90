module comvar
  implicit none
  integer, parameter :: ndx=128,ndy=128,laststep=2000,istx=1,ienx=2,isty=1,ieny=2,svnum=20, dim=4 !preiodic:ist=1,ien=2 , kotei:ist=2,ien=3 : ndx=130
  integer, parameter :: itel=1,iwv=1
  !double precision, parameter :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0
  integer :: iwx,iwy,iwz,bndx0=4,bndy0=4,bndx1=3,bndy1=3 !odd:x, even:y, 1,2:periodic, 3,4:exact, 5,6:exact+free
  DOUBLE PRECISION :: cg = 1.0d0, Tdiff=1.0d1 , dx,dy,ratiodini=1.d0 != Lbox/dble(ndx-2) !, bcphi1 , bcphi2
  double precision :: Lbox=1.0d2 , h=50.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0 ,rsph=10.d0, rch=1.d0,Cnst=0.d0
  DOUBLE PRECISION , dimension(-1:ndx,-1:ndy) :: Phidt,Phiexa,Phicrr
  DOUBLE PRECISION , dimension(-1-1:ndx+1,-1-1:ndy+1) :: Phiexa2,rho, Qgr
  DOUBLE PRECISION , dimension(-1:ndx,-1:ndy,dim) :: Phigrd
  !double precision :: G=1.11142d-4, G4pi=12.56637d0*G , coeff=0.90d0 ,  kappa=1.0d0/3.0d0
  double precision ::  G4pi=12.56637d0*1.11142d-4 , coeff=1.d0 ,meanrho,meanphiexa!,  kappa=1.0d0/3.0d0
  DOUBLE PRECISION , dimension(1:3) :: bcphi1 , bcphi2 ,bcphigrd1 , bcphigrd2
  !character(63) :: dir='/Users/maeda/Desktop/Dropbox/analysis/telegraph-test-2D-1/simu/'
  character(79) :: dir='/Users/maeda/Desktop/Dropbox/analysis/telegraph-test-2D-1/cg1-T10-rho1-cen-128/'
end module comvar

module grvvar
  implicit none
  integer, parameter :: ndx2=130 , ndy2=130 , dim2=4 !パラメータ属性必要
  DOUBLE PRECISION , dimension(-1-1:ndx2+1) :: x
  DOUBLE PRECISION , dimension(-1-1:ndy2+1) :: y
  DOUBLE PRECISION , dimension(-1:ndx2,dim2) :: flmtel1,flmtel2,flmtel3,flmtel4
  DOUBLE PRECISION , dimension(-1:ndx2,-1:ndy2) :: Phidtn  , fx , fy
  DOUBLE PRECISION , dimension(-1:ndx2,-1:ndy2,dim2) :: wp1 ,wp2, wppre1 ,wppre2 !wp->telegraph
  DOUBLE PRECISION , dimension(-1:ndx2,-1:ndy2,dim2) :: fp1 ,fp2 !fp->wave
  DOUBLE PRECISION , dimension(-1:ndx2,-1:ndy2,2) :: source ,sourcedt,sourcedt2
end module grvvar

program main
  use comvar
  use grvvar
  integer :: i,in1,j,k,time=100
  character(5) name
  double precision meanphi1,meanphi2,errphi1,errphi2,dm,phimax,phimin,dm2
  
open(850,file=dir//'phierr.dat',FORM='FORMATTED',position='append')

do k=0,time
  meanphi1=0.d0
  meanphi2=0.d0
  meanphiexa=0.d0
  phimax=-1.d10
  phimin=1.d10
  write(name,'(i5.5)') k
  open(21,file=dir//'phi2D'//name//'.dat')
  do j=-1,ndy
     do i=-1,ndx
        read(21,*) x(i),y(j),wp1(i,j,1),wp1(i,j,2),wp1(i,j,3),wp1(i,j,4) &
        ,wp2(i,j,1),wp2(i,j,2),wp2(i,j,3),wp2(i,j,4)&
        ,fp1(i,j,1),fp1(i,j,2),fp1(i,j,3),fp1(i,j,4) &
        ,fp2(i,j,1),fp2(i,j,2),fp2(i,j,3),fp2(i,j,4)&
        ,dm2,rho(i,j),Phicrr(i,j),Phiexa(i,j),&
        Phigrd(i,j,1),Phigrd(i,j,2),Phigrd(i,j,3),Phigrd(i,j,4),dm2,dm2

        meanphi1=0.25d0*(wp2(i,j,1)+wp2(i,j,2)+wp2(i,j,3)+wp2(i,j,4))+meanphi1
        meanphi2=wp2(i,j,1)+meanphi2
        meanphiexa=Phiexa(i,j)+meanphiexa
        phimax=dmax1(phimax,Phiexa(i,j))
        phimin=dmin1(phimin,Phiexa(i,j))
     end do
     !write(*,*) 'ok'
     read (21,*) dm
  end do
  close(21)

  do j=1,ndy-2
     do i=1,ndx-2
     errphi1=(Phiexa(i,j)+meanphi1-meanphiexa)**2
     errphi2=(Phiexa(i,j)+meanphi2-meanphiexa)**2
     enddo
  end do
  errphi1=dsqrt(errphi1)
  errphi2=dsqrt(errphi2)
  write(850,*) i*svnum,',',errphi1,',', errphi2,',', meanphiexa,',',&
   meanphi1,',', meanphi2,',', errphi1/meanphi1,',', errphi2/meanphi2

end do

close(850)

end program main

