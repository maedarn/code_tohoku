program main
  implicit none
  DOUBLE PRECISION, parameter :: mH=1.d0, mHe=4.d0, mH2=2.d0, mC=12.d0, mCO=28.d0
  DOUBLE PRECISION, parameter :: G0=1.d0, xc=1.4d-4, xo=3.2d-4, dv=2.d0, Tgr=5.d-3
  DOUBLE PRECISION, dimension(:,:,:)  , allocatable :: ndp,ndH,ndH2,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot
  DOUBLE PRECISION, dimension(:,:,:,:), allocatable :: Ntot,NH2,NnC,NCO,tCII
  DOUBLE PRECISION  :: ndpmin,ndHmin,ndH2min,ndHemin,ndHepmin,ndCmin,ndCpmin,ndCOmin
  integer :: Np1x, Np2x, Np1y, Np2y, Np1z, Np2z, nunit, ix, jy, kz
  double precision ::  ql1x,ql2x,ql1y,ql2y,ql1z,ql2z,dinit1,dinit2,pinit1,pinit2, &
       vinitx1,vinitx2,vinity1,vinity2,vinitz1,vinitz2,           &
       binitx1,binitx2,binity1,binity2,binitz1,binitz2
  double precision, dimension(:), allocatable :: x_i,y_i,z_i,dx_i,dy_i,dz_i
  double precision :: theta,amp,xpi,ypi,zpi,phase1,phase2,phase3,kx,ky,kzz,kw,v_bai,b_bai
  double precision :: Hini,pini,H2ini,Heini,Hepini,Cini,COini,Cpini,dBC,kb,vx,by,b,pi

  pi=3.1415926535897932384626d0
!*********parameter**********
  b_bai=1.0d0
  v_bai=4.0d0
  theta = 0.25d0*pi
!*********parameter**********
  kb=8.63359d0
  vx=9.778d0 !1km/s
  by=2.2307d0 !1μG


  pinit1=8.810807d3*kb*1.d-3; pinit2=pinit1
  Hini=0.9219098d0; pini=0.9503446d-2; H2ini=0.9465513d-8; Heini=0.9155226d-1; Hepini=0.5655353d-3
  Cini=0.1565848d-8; COini=0.2202631d-20; Cpini=0.1433520d-3
  dinit1=mH*Hini+mH*pini+mH2*H2ini+mHe*Heini+mHe*Hepini; dinit2=dinit1


  vinitx1=vx * v_bai
  b=by * b_bai
!write(*,*) pi
  binity1= b * dcos(theta)
  binitx1= b * dsin(theta)
!write(*,*) dcos(theta),dsin(theta) ,theta

  write(*,*) dinit1
  write(*,*) vinitx1
  write(*,*) binitx1
  write(*,*) binity1
end program main
