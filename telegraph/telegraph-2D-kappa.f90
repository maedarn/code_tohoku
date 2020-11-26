module comvar
  implicit none
  integer, parameter :: ndx=66,ndy=66,laststep=2000,istx=1,ienx=2,isty=1,ieny=2,svnum=20, dim=4 !preiodic:ist=1,ien=2 , kotei:ist=2,ien=3 : ndx=130
  integer, parameter :: itel=1,iwv=1
  !double precision, parameter :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0
  integer :: iwx,iwy,iwz,bndx0=4,bndy0=4,bndx1=3,bndy1=3 !odd:x, even:y, 1,2:periodic, 3,4:exact, 5,6:exact+free
  DOUBLE PRECISION :: cg = 1.0d0, kappa=0.50d0 , dx,dy != Lbox/dble(ndx-2) !, bcphi1 , bcphi2
  double precision :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0 ,rsph=10.d0, rch=1.d0,Cnst=0.d0
  DOUBLE PRECISION , dimension(-1:ndx,-1:ndy) :: Phidt,Phiexa,Phicrr
  DOUBLE PRECISION , dimension(-1-1:ndx+1,-1-1:ndy+1) :: Phiexa2,rho, Qgr
  DOUBLE PRECISION , dimension(-1:ndx,-1:ndy,dim) :: Phigrd
  !double precision :: G=1.11142d-4, G4pi=12.56637d0*G , coeff=0.90d0 ,  kappa=1.0d0/3.0d0
  double precision ::  G4pi=12.56637d0*1.11142d-4 , coeff=0.5d0 ,meanrho,meanphiexa!,  kappa=1.0d0/3.0d0
  DOUBLE PRECISION , dimension(1:3) :: bcphi1 , bcphi2 ,bcphigrd1 , bcphigrd2
  character(63) :: dir='/Users/maeda/Desktop/Dropbox/analysis/telegraph-test-2D-1/simu/'
end module comvar

module grvvar
  implicit none
  integer, parameter :: ndx2=66 , ndy2=66 , dim2=4 !パラメータ属性必要
  DOUBLE PRECISION , dimension(-1-1:ndx2+1) :: x
  DOUBLE PRECISION , dimension(-1-1:ndy2+1) :: y
  DOUBLE PRECISION , dimension(-1:ndx2,-1:ndy2) :: Phidtn  , fx , fy
  DOUBLE PRECISION , dimension(-1:ndx2,-1:ndy2,dim2) :: wp1 ,wp2, wppre1 ,wppre2 !wp->telegraph
  DOUBLE PRECISION , dimension(-1:ndx2,-1:ndy2,dim2) :: fp1 ,fp2 !fp->wave
  DOUBLE PRECISION , dimension(-1:ndx2,-1:ndy2,2) :: source ,sourcedt,sourcedt2
end module grvvar

program muscl1D
  !implicit none :: まちがった位置
  use comvar
  use grvvar
  implicit none
  DOUBLE PRECISION :: dt=0.0d0
  integer :: i,sv=0,iws,ws=2,j,k
  integer :: n,m


  call INITIAL()
  !call BC(bndx,bndx)
  !call BC(bndy,bndy)
  !call muslcslv1D(Phi,Phi1step,dt,13)
  write(*,*) 'here'
  !call saveu(sv)
  do i=1,laststep

     if(mod(i,svnum)==0) then
     call saveu(sv)
     end if

     call time(dt)
     !call timesource(Phidtn,rho,dt,3)
     write(*,*) i ,dt,'step'

    !goto 4545
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iwx=1;iwy=0
    !call BC(wp2(:,:,1),3,2,0,0,1)
    !call BC(wp2(:,:,2),2,3,0,0,2)
    !call BC(wp2(:,:,3),3,2,0,0,3)
    !call BC(wp2(:,:,4),2,3,0,0,4)
    if(itel==1) then
    call BC(wp2(:,:,1),3,3,0,0,1)
    call BC(wp2(:,:,2),3,3,0,0,2)
    call BC(wp2(:,:,3),3,3,0,0,3)
    call BC(wp2(:,:,4),3,3,0,0,4)
    endif

    if(iwv==1) then
    call BC(fp2(:,:,1),13,13,0,0,1)
    call BC(fp2(:,:,2),13,13,0,0,2)
    call BC(fp2(:,:,3),13,13,0,0,3)
    call BC(fp2(:,:,4),13,13,0,0,4)
    endif

    if(itel==1) then
    call muslcslv1D(wp2(:,:,1),rho,dt*0.25d0,1,2)
    call muslcslv1D(wp2(:,:,2),rho,dt*0.25d0,2,2)
    call muslcslv1D(wp2(:,:,3),rho,dt*0.25d0,1,2)
    call muslcslv1D(wp2(:,:,4),rho,dt*0.25d0,2,2)
    endif

    if(iwv==1) then
    call muslcslv1D(fp2(:,:,1),rho,dt*0.25d0,1,2)
    call muslcslv1D(fp2(:,:,2),rho,dt*0.25d0,2,2)
    call muslcslv1D(fp2(:,:,3),rho,dt*0.25d0,1,2)
    call muslcslv1D(fp2(:,:,4),rho,dt*0.25d0,2,2)
    endif

    !call BCtest()
    iwx=0;iwy=1
    !call BC(wp2(:,:,1),0,0,3,2,1)
    !call BC(wp2(:,:,2),0,0,2,3,2)
    !call BC(wp2(:,:,3),0,0,2,3,3)
    !call BC(wp2(:,:,4),0,0,3,2,4)

    if(itel==1) then
    call BC(wp2(:,:,1),0,0,3,3,1)
    call BC(wp2(:,:,2),0,0,3,3,2)
    call BC(wp2(:,:,3),0,0,3,3,3)
    call BC(wp2(:,:,4),0,0,3,3,4)
    endif

    if(iwv==1) then
    call BC(fp2(:,:,1),0,0,13,13,1)
    call BC(fp2(:,:,2),0,0,13,13,2)
    call BC(fp2(:,:,3),0,0,13,13,3)
    call BC(fp2(:,:,4),0,0,13,13,4)
    endif

    if(itel==1) then
    call muslcslv1D(wp2(:,:,1),rho,dt*0.25d0,1,2)
    call muslcslv1D(wp2(:,:,2),rho,dt*0.25d0,2,2)
    call muslcslv1D(wp2(:,:,3),rho,dt*0.25d0,2,2)
    call muslcslv1D(wp2(:,:,4),rho,dt*0.25d0,1,2)
    endif

    if(iwv==1) then
    call muslcslv1D(fp2(:,:,1),rho,dt*0.25d0,1,2)
    call muslcslv1D(fp2(:,:,2),rho,dt*0.25d0,2,2)
    call muslcslv1D(fp2(:,:,3),rho,dt*0.25d0,2,2)
    call muslcslv1D(fp2(:,:,4),rho,dt*0.25d0,1,2)
    endif
    !call BCtest()

    if(itel==1) then
    do k=1,ndy-2
      do j=1,ndx-2
wp2(j,k,1) = wp2(j,k,1)*dexp(-kappa * 0.5d0 * dt)+wp1(j,k,1)*(1.d0-dexp(-kappa * 0.5d0 * dt))
wp2(j,k,2) = wp2(j,k,2)*dexp(-kappa * 0.5d0 * dt)+wp1(j,k,2)*(1.d0-dexp(-kappa * 0.5d0 * dt))
wp2(j,k,3) = wp2(j,k,3)*dexp(-kappa * 0.5d0 * dt)+wp1(j,k,3)*(1.d0-dexp(-kappa * 0.5d0 * dt))
wp2(j,k,4) = wp2(j,k,4)*dexp(-kappa * 0.5d0 * dt)+wp1(j,k,4)*(1.d0-dexp(-kappa * 0.5d0 * dt))
      enddo
    enddo
    endif
    
    if(iwv==1) then
    do k=1,ndy-2
        do j=1,ndx-2
    fp2(j,k,1) = fp2(j,k,1)+cg*dt*fp1(j,k,1)*0.5d0
    fp2(j,k,2) = fp2(j,k,2)+cg*dt*fp1(j,k,2)*0.5d0
    fp2(j,k,3) = fp2(j,k,3)+cg*dt*fp1(j,k,3)*0.5d0
    fp2(j,k,4) = fp2(j,k,4)+cg*dt*fp1(j,k,4)*0.5d0
        enddo
    enddo
    endif
!do k=1,ndy-2
!  do j=1,ndx-2
!  wp2(j,k,1) = wp2(j,k,1)*dexp(-kappa * dt*0.5d0) &
!  +(wp1(j,k,1) - cg*cg*4.d0*1.d0/kappa/2.d0*1.d0/kappa/2.d0*G4pi*rho(j,k))*(1.d0-dexp(-kappa * dt*0.5d0))
!  wp2(j,k,2) = wp2(j,k,2)*dexp(-kappa * dt*0.5d0) &
!  +(wp1(j,k,2) - cg*cg*4.d0*1.d0/kappa/2.d0*1.d0/kappa/2.d0*G4pi*rho(j,k))*(1.d0-dexp(-kappa * dt*0.5d0))
!  wp2(j,k,3) = wp2(j,k,3)*dexp(-kappa * dt*0.5d0) &
!  +(wp1(j,k,3) - cg*cg*4.d0*1.d0/kappa/2.d0*1.d0/kappa/2.d0*G4pi*rho(j,k))*(1.d0-dexp(-kappa * dt*0.5d0))
!  wp2(j,k,4) = wp2(j,k,4)*dexp(-kappa * dt*0.5d0) &
!  +(wp1(j,k,4) - cg*cg*4.d0*1.d0/kappa/2.d0*1.d0/kappa/2.d0*G4pi*rho(j,k))*(1.d0-dexp(-kappa * dt*0.5d0))
!  enddo
!enddo
    !call BC(wp2(:,:,1),3,3,3,3,1)
    !call BC(wp2(:,:,2),3,3,3,3,2)
    !call BC(wp2(:,:,3),3,3,3,3,3)
    !call BC(wp2(:,:,4),3,3,3,3,4)
    !call BCtest()
    iwx=0;iwy=1
    !call BC(wp2(:,:,1),0,0,3,2,1)
    !call BC(wp2(:,:,2),0,0,2,3,2)
    !call BC(wp2(:,:,3),0,0,2,3,3)
    !call BC(wp2(:,:,4),0,0,3,2,4)
    if(itel==1) then
    call BC(wp2(:,:,1),0,0,3,3,1)
    call BC(wp2(:,:,2),0,0,3,3,2)
    call BC(wp2(:,:,3),0,0,3,3,3)
    call BC(wp2(:,:,4),0,0,3,3,4)
    endif

    if(iwv==1) then
    call BC(fp2(:,:,1),0,0,13,13,1)
    call BC(fp2(:,:,2),0,0,13,13,2)
    call BC(fp2(:,:,3),0,0,13,13,3)
    call BC(fp2(:,:,4),0,0,13,13,4)
    endif

    if(itel==1) then
    call muslcslv1D(wp2(:,:,1),rho,dt*0.25d0,1,2)
    call muslcslv1D(wp2(:,:,2),rho,dt*0.25d0,2,2)
    call muslcslv1D(wp2(:,:,3),rho,dt*0.25d0,2,2)
    call muslcslv1D(wp2(:,:,4),rho,dt*0.25d0,1,2)
    endif


    if(iwv==1) then
    call muslcslv1D(fp2(:,:,1),rho,dt*0.25d0,1,2)
    call muslcslv1D(fp2(:,:,2),rho,dt*0.25d0,2,2)
    call muslcslv1D(fp2(:,:,3),rho,dt*0.25d0,2,2)
    call muslcslv1D(fp2(:,:,4),rho,dt*0.25d0,1,2)
    endif
    !call BCtest()
    iwx=1;iwy=0
    !call BC(wp2(:,:,1),3,2,0,0,1)
    !call BC(wp2(:,:,2),2,3,0,0,2)
    !call BC(wp2(:,:,3),3,2,0,0,3)
    !call BC(wp2(:,:,4),2,3,0,0,4)
    if(itel==1) then
    call BC(wp2(:,:,1),3,3,0,0,1)
    call BC(wp2(:,:,2),3,3,0,0,2)
    call BC(wp2(:,:,3),3,3,0,0,3)
    call BC(wp2(:,:,4),3,3,0,0,4)
    endif

    if(iwv==1) then
    call BC(fp2(:,:,1),13,13,0,0,1)
    call BC(fp2(:,:,2),13,13,0,0,2)
    call BC(fp2(:,:,3),13,13,0,0,3)
    call BC(fp2(:,:,4),13,13,0,0,4)
    endif

    if(itel==1) then
    call muslcslv1D(wp2(:,:,1),rho,dt*0.25d0,1,2)
    call muslcslv1D(wp2(:,:,2),rho,dt*0.25d0,2,2)
    call muslcslv1D(wp2(:,:,3),rho,dt*0.25d0,1,2)
    call muslcslv1D(wp2(:,:,4),rho,dt*0.25d0,2,2)
    endif

    if(iwv==1) then
    call muslcslv1D(fp2(:,:,1),rho,dt*0.25d0,1,2)
    call muslcslv1D(fp2(:,:,2),rho,dt*0.25d0,2,2)
    call muslcslv1D(fp2(:,:,3),rho,dt*0.25d0,1,2)
    call muslcslv1D(fp2(:,:,4),rho,dt*0.25d0,2,2)
    endif
    !call BCtest()
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    iwx=1;iwy=0
    !call BC(wp1(:,:,1),2,4,0,0,1)
    !call BC(wp1(:,:,2),4,2,0,0,2)
    !call BC(wp1(:,:,3),2,4,0,0,3)
    !call BC(wp1(:,:,4),4,2,0,0,4)
    if(itel==1) then
    call BC(wp1(:,:,1),4,4,0,0,1)
    call BC(wp1(:,:,2),4,4,0,0,2)
    call BC(wp1(:,:,3),4,4,0,0,3)
    call BC(wp1(:,:,4),4,4,0,0,4)
    endif

    if(iwv==1) then
    call BC(fp1(:,:,1),14,14,0,0,1)
    call BC(fp1(:,:,2),14,14,0,0,2)
    call BC(fp1(:,:,3),14,14,0,0,3)
    call BC(fp1(:,:,4),14,14,0,0,4)
    endif

    if(itel==1) then
    call muslcslv1D(wp1(:,:,1),rho,dt*0.5d0,2,2)
    call muslcslv1D(wp1(:,:,2),rho,dt*0.5d0,1,2)
    call muslcslv1D(wp1(:,:,3),rho,dt*0.5d0,2,2)
    call muslcslv1D(wp1(:,:,4),rho,dt*0.5d0,1,2)
    endif

    if(iwv==1) then
    call muslcslv1D(fp1(:,:,1),rho,dt*0.5d0,2,2)
    call muslcslv1D(fp1(:,:,2),rho,dt*0.5d0,1,2)
    call muslcslv1D(fp1(:,:,3),rho,dt*0.5d0,2,2)
    call muslcslv1D(fp1(:,:,4),rho,dt*0.5d0,1,2)
    endif
    !call BCtest()
    iwx=0;iwy=1
    !call BC(wp1(:,:,1),0,0,2,4,1)
    !call BC(wp1(:,:,2),0,0,4,2,2)
    !call BC(wp1(:,:,3),0,0,4,2,3)
    !call BC(wp1(:,:,4),0,0,2,4,4)
    if(itel==1) then
    call BC(wp1(:,:,1),0,0,4,4,1)
    call BC(wp1(:,:,2),0,0,4,4,2)
    call BC(wp1(:,:,3),0,0,4,4,3)
    call BC(wp1(:,:,4),0,0,4,4,4)
    endif

    if(iwv==1) then
    call BC(fp1(:,:,1),0,0,14,14,1)
    call BC(fp1(:,:,2),0,0,14,14,2)
    call BC(fp1(:,:,3),0,0,14,14,3)
    call BC(fp1(:,:,4),0,0,14,14,4)
    endif

    if(itel==1) then
    call muslcslv1D(wp1(:,:,1),rho,dt*0.5d0,2,2)
    call muslcslv1D(wp1(:,:,2),rho,dt*0.5d0,1,2)
    call muslcslv1D(wp1(:,:,3),rho,dt*0.5d0,1,2)
    call muslcslv1D(wp1(:,:,4),rho,dt*0.5d0,2,2)
    endif

    if(iwv==1) then
    call muslcslv1D(fp1(:,:,1),rho,dt*0.5d0,2,2)
    call muslcslv1D(fp1(:,:,2),rho,dt*0.5d0,1,2)
    call muslcslv1D(fp1(:,:,3),rho,dt*0.5d0,1,2)
    call muslcslv1D(fp1(:,:,4),rho,dt*0.5d0,2,2)
    endif
    !call BCtest()

    if(itel==1) then
    do k=1,ndy-2
      do j=1,ndx-2
wp1(j,k,1) = wp1(j,k,1)*dexp(-kappa * dt)+(wp2(j,k,1) &
-4.d0*1.d0/kappa/2.d0*cg*cg*(wp2(j+1,k+1,1)-wp2(j+1,k-1,1)-wp2(j-1,k+1,1)+wp2(j-1,k-1,1))/dx/dy &
-cg*cg*4.d0*1.d0/kappa/2.d0*1.d0/kappa/2.d0*G4pi*rho(j,k))*(1.d0-dexp(-kappa * dt))
wp1(j,k,2) = wp1(j,k,2)*dexp(-kappa * dt)+(wp2(j,k,2) &
-4.d0*1.d0/kappa/2.d0*cg*cg*(wp2(j+1,k+1,2)-wp2(j+1,k-1,2)-wp2(j-1,k+1,2)+wp2(j-1,k-1,2))/dx/dy &
-cg*cg*4.d0*1.d0/kappa/2.d0*1.d0/kappa/2.d0*G4pi*rho(j,k))*(1.d0-dexp(-kappa * dt))
wp1(j,k,3) = wp1(j,k,3)*dexp(-kappa * dt)+(wp2(j,k,3) &
+4.d0*1.d0/kappa/2.d0*cg*cg*(wp2(j+1,k+1,3)-wp2(j+1,k-1,3)-wp2(j-1,k+1,3)+wp2(j-1,k-1,3))/dx/dy &
-cg*cg*4.d0*1.d0/kappa/2.d0*1.d0/kappa/2.d0*G4pi*rho(j,k))*(1.d0-dexp(-kappa * dt))
wp1(j,k,4) = wp1(j,k,4)*dexp(-kappa * dt)+(wp2(j,k,4) &
+4.d0*1.d0/kappa/2.d0*cg*cg*(wp2(j+1,k+1,4)-wp2(j+1,k-1,4)-wp2(j-1,k+1,4)+wp2(j-1,k-1,4))/dx/dy &
-cg*cg*4.d0*1.d0/kappa/2.d0*1.d0/kappa/2.d0*G4pi*rho(j,k))*(1.d0-dexp(-kappa * dt))
      enddo
    enddo
    endif


    if(iwv==1) then
    do k=1,ndy-2
      do j=1,ndx-2
fp1(j,k,1) = fp1(j,k,1)-cg*G4pi*rho(j,k)*dt &
-0.5d0*dt*cg*(fp2(j+1,k+1,1)-fp2(j+1,k-1,1)-fp2(j-1,k+1,1)+fp2(j-1,k-1,1))/dx/dy
fp1(j,k,2) = fp1(j,k,2)-cg*G4pi*rho(j,k)*dt &
-0.5d0*dt*cg*(fp2(j+1,k+1,2)-fp2(j+1,k-1,2)-fp2(j-1,k+1,2)+fp2(j-1,k-1,2))/dx/dy
fp1(j,k,3) = fp1(j,k,3)-cg*G4pi*rho(j,k)*dt &
+0.5d0*dt*cg*(fp2(j+1,k+1,3)-fp2(j+1,k-1,3)-fp2(j-1,k+1,3)+fp2(j-1,k-1,3))/dx/dy
fp1(j,k,4) = fp1(j,k,4)-cg*G4pi*rho(j,k)*dt &
+0.5d0*dt*cg*(fp2(j+1,k+1,4)-fp2(j+1,k-1,4)-fp2(j-1,k+1,4)+fp2(j-1,k-1,4))/dx/dy
      enddo
    enddo
    endif
    !call BC(wp1(:,:,1),4,4,4,4,1)
    !call BC(wp1(:,:,2),4,4,4,4,2)
    !call BC(wp1(:,:,3),4,4,4,4,3)
    !call BC(wp1(:,:,4),4,4,4,4,4)
    !call BCtest()
    !call BC(wp1(:,:,1),0,0,2,4,1)
    !call BC(wp1(:,:,2),0,0,4,2,2)
    !call BC(wp1(:,:,3),0,0,4,2,3)
    !call BC(wp1(:,:,4),0,0,2,4,4)
    iwx=0;iwy=1
    if(itel==1) then
    call BC(wp1(:,:,1),0,0,4,4,1)
    call BC(wp1(:,:,2),0,0,4,4,2)
    call BC(wp1(:,:,3),0,0,4,4,3)
    call BC(wp1(:,:,4),0,0,4,4,4)
    endif

    if(iwv==1) then
    call BC(fp1(:,:,1),0,0,14,14,1)
    call BC(fp1(:,:,2),0,0,14,14,2)
    call BC(fp1(:,:,3),0,0,14,14,3)
    call BC(fp1(:,:,4),0,0,14,14,4)
    endif

    if(itel==1) then
    call muslcslv1D(wp1(:,:,1),rho,dt*0.5d0,2,2)
    call muslcslv1D(wp1(:,:,2),rho,dt*0.5d0,1,2)
    call muslcslv1D(wp1(:,:,3),rho,dt*0.5d0,1,2)
    call muslcslv1D(wp1(:,:,4),rho,dt*0.5d0,2,2)
    endif

    if(iwv==1) then
    call muslcslv1D(fp1(:,:,1),rho,dt*0.5d0,2,2)
    call muslcslv1D(fp1(:,:,2),rho,dt*0.5d0,1,2)
    call muslcslv1D(fp1(:,:,3),rho,dt*0.5d0,1,2)
    call muslcslv1D(fp1(:,:,4),rho,dt*0.5d0,2,2)
    endif
    !call BCtest()
    iwx=1;iwy=0
    !call BC(wp1(:,:,1),2,4,0,0,1)
    !call BC(wp1(:,:,2),4,2,0,0,2)
    !call BC(wp1(:,:,3),2,4,0,0,3)
    !call BC(wp1(:,:,4),4,2,0,0,4)

    if(itel==1) then
    call BC(wp1(:,:,1),4,4,0,0,1)
    call BC(wp1(:,:,2),4,4,0,0,2)
    call BC(wp1(:,:,3),4,4,0,0,3)
    call BC(wp1(:,:,4),4,4,0,0,4)
    endif

    if(iwv==1) then
    call BC(fp1(:,:,1),14,14,0,0,1)
    call BC(fp1(:,:,2),14,14,0,0,2)
    call BC(fp1(:,:,3),14,14,0,0,3)
    call BC(fp1(:,:,4),14,14,0,0,4)
    endif

    if(itel==1) then
    call muslcslv1D(wp1(:,:,1),rho,dt*0.5d0,2,2)
    call muslcslv1D(wp1(:,:,2),rho,dt*0.5d0,1,2)
    call muslcslv1D(wp1(:,:,3),rho,dt*0.5d0,2,2)
    call muslcslv1D(wp1(:,:,4),rho,dt*0.5d0,1,2)
    endif

    if(iwv==1) then
    call muslcslv1D(fp1(:,:,1),rho,dt*0.5d0,2,2)
    call muslcslv1D(fp1(:,:,2),rho,dt*0.5d0,1,2)
    call muslcslv1D(fp1(:,:,3),rho,dt*0.5d0,2,2)
    call muslcslv1D(fp1(:,:,4),rho,dt*0.5d0,1,2)
    endif
    !call BCtest()


       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       iwx=1;iwy=0
       !call BC(wp2(:,:,1),3,2,0,0,1)
       !call BC(wp2(:,:,2),2,3,0,0,2)
       !call BC(wp2(:,:,3),3,2,0,0,3)
       !call BC(wp2(:,:,4),2,3,0,0,4)
       if(itel==1) then
       call BC(wp2(:,:,1),3,3,0,0,1)
       call BC(wp2(:,:,2),3,3,0,0,2)
       call BC(wp2(:,:,3),3,3,0,0,3)
       call BC(wp2(:,:,4),3,3,0,0,4)
       endif

       if(iwv==1) then
       call BC(fp2(:,:,1),13,13,0,0,1)
       call BC(fp2(:,:,2),13,13,0,0,2)
       call BC(fp2(:,:,3),13,13,0,0,3)
       call BC(fp2(:,:,4),13,13,0,0,4)
       endif

       if(itel==1) then
       call muslcslv1D(wp2(:,:,1),rho,dt*0.25d0,1,2)
       call muslcslv1D(wp2(:,:,2),rho,dt*0.25d0,2,2)
       call muslcslv1D(wp2(:,:,3),rho,dt*0.25d0,1,2)
       call muslcslv1D(wp2(:,:,4),rho,dt*0.25d0,2,2)
       endif

       if(iwv==1) then
       call muslcslv1D(fp2(:,:,1),rho,dt*0.25d0,1,2)
       call muslcslv1D(fp2(:,:,2),rho,dt*0.25d0,2,2)
       call muslcslv1D(fp2(:,:,3),rho,dt*0.25d0,1,2)
       call muslcslv1D(fp2(:,:,4),rho,dt*0.25d0,2,2)
       endif

       !call BCtest()
       iwx=0;iwy=1
       !call BC(wp2(:,:,1),0,0,3,2,1)
       !call BC(wp2(:,:,2),0,0,2,3,2)
       !call BC(wp2(:,:,3),0,0,2,3,3)
       !call BC(wp2(:,:,4),0,0,3,2,4)

       if(itel==1) then
       call BC(wp2(:,:,1),0,0,3,3,1)
       call BC(wp2(:,:,2),0,0,3,3,2)
       call BC(wp2(:,:,3),0,0,3,3,3)
       call BC(wp2(:,:,4),0,0,3,3,4)
       endif

       if(iwv==1) then
       call BC(fp2(:,:,1),0,0,13,13,1)
       call BC(fp2(:,:,2),0,0,13,13,2)
       call BC(fp2(:,:,3),0,0,13,13,3)
       call BC(fp2(:,:,4),0,0,13,13,4)
       endif

       if(itel==1) then
       call muslcslv1D(wp2(:,:,1),rho,dt*0.25d0,1,2)
       call muslcslv1D(wp2(:,:,2),rho,dt*0.25d0,2,2)
       call muslcslv1D(wp2(:,:,3),rho,dt*0.25d0,2,2)
       call muslcslv1D(wp2(:,:,4),rho,dt*0.25d0,1,2)
       endif

       if(iwv==1) then
       call muslcslv1D(fp2(:,:,1),rho,dt*0.25d0,1,2)
       call muslcslv1D(fp2(:,:,2),rho,dt*0.25d0,2,2)
       call muslcslv1D(fp2(:,:,3),rho,dt*0.25d0,2,2)
       call muslcslv1D(fp2(:,:,4),rho,dt*0.25d0,1,2)
       endif
       !call BCtest()

       if(itel==1) then
       do k=1,ndy-2
         do j=1,ndx-2
   wp2(j,k,1) = wp2(j,k,1)*dexp(-kappa * 0.5d0 * dt)+wp1(j,k,1)*(1.d0-dexp(-kappa * 0.5d0 * dt))
   wp2(j,k,2) = wp2(j,k,2)*dexp(-kappa * 0.5d0 * dt)+wp1(j,k,2)*(1.d0-dexp(-kappa * 0.5d0 * dt))
   wp2(j,k,3) = wp2(j,k,3)*dexp(-kappa * 0.5d0 * dt)+wp1(j,k,3)*(1.d0-dexp(-kappa * 0.5d0 * dt))
   wp2(j,k,4) = wp2(j,k,4)*dexp(-kappa * 0.5d0 * dt)+wp1(j,k,4)*(1.d0-dexp(-kappa * 0.5d0 * dt))
         enddo
       enddo
       endif
       
       if(iwv==1) then
       do k=1,ndy-2
           do j=1,ndx-2
       fp2(j,k,1) = fp2(j,k,1)+cg*dt*fp1(j,k,1)*0.5d0
       fp2(j,k,2) = fp2(j,k,2)+cg*dt*fp1(j,k,2)*0.5d0
       fp2(j,k,3) = fp2(j,k,3)+cg*dt*fp1(j,k,3)*0.5d0
       fp2(j,k,4) = fp2(j,k,4)+cg*dt*fp1(j,k,4)*0.5d0
           enddo
       enddo
       endif
   !do k=1,ndy-2
   !  do j=1,ndx-2
   !  wp2(j,k,1) = wp2(j,k,1)*dexp(-kappa * dt*0.5d0) &
   !  +(wp1(j,k,1) - cg*cg*4.d0*1.d0/kappa/2.d0*1.d0/kappa/2.d0*G4pi*rho(j,k))*(1.d0-dexp(-kappa * dt*0.5d0))
   !  wp2(j,k,2) = wp2(j,k,2)*dexp(-kappa * dt*0.5d0) &
   !  +(wp1(j,k,2) - cg*cg*4.d0*1.d0/kappa/2.d0*1.d0/kappa/2.d0*G4pi*rho(j,k))*(1.d0-dexp(-kappa * dt*0.5d0))
   !  wp2(j,k,3) = wp2(j,k,3)*dexp(-kappa * dt*0.5d0) &
   !  +(wp1(j,k,3) - cg*cg*4.d0*1.d0/kappa/2.d0*1.d0/kappa/2.d0*G4pi*rho(j,k))*(1.d0-dexp(-kappa * dt*0.5d0))
   !  wp2(j,k,4) = wp2(j,k,4)*dexp(-kappa * dt*0.5d0) &
   !  +(wp1(j,k,4) - cg*cg*4.d0*1.d0/kappa/2.d0*1.d0/kappa/2.d0*G4pi*rho(j,k))*(1.d0-dexp(-kappa * dt*0.5d0))
   !  enddo
   !enddo
       !call BC(wp2(:,:,1),3,3,3,3,1)
       !call BC(wp2(:,:,2),3,3,3,3,2)
       !call BC(wp2(:,:,3),3,3,3,3,3)
       !call BC(wp2(:,:,4),3,3,3,3,4)
       !call BCtest()
       iwx=0;iwy=1
       !call BC(wp2(:,:,1),0,0,3,2,1)
       !call BC(wp2(:,:,2),0,0,2,3,2)
       !call BC(wp2(:,:,3),0,0,2,3,3)
       !call BC(wp2(:,:,4),0,0,3,2,4)
       if(itel==1) then
       call BC(wp2(:,:,1),0,0,3,3,1)
       call BC(wp2(:,:,2),0,0,3,3,2)
       call BC(wp2(:,:,3),0,0,3,3,3)
       call BC(wp2(:,:,4),0,0,3,3,4)
       endif

       if(iwv==1) then
       call BC(fp2(:,:,1),0,0,13,13,1)
       call BC(fp2(:,:,2),0,0,13,13,2)
       call BC(fp2(:,:,3),0,0,13,13,3)
       call BC(fp2(:,:,4),0,0,13,13,4)
       endif

       if(itel==1) then
       call muslcslv1D(wp2(:,:,1),rho,dt*0.25d0,1,2)
       call muslcslv1D(wp2(:,:,2),rho,dt*0.25d0,2,2)
       call muslcslv1D(wp2(:,:,3),rho,dt*0.25d0,2,2)
       call muslcslv1D(wp2(:,:,4),rho,dt*0.25d0,1,2)
       endif


       if(iwv==1) then
       call muslcslv1D(fp2(:,:,1),rho,dt*0.25d0,1,2)
       call muslcslv1D(fp2(:,:,2),rho,dt*0.25d0,2,2)
       call muslcslv1D(fp2(:,:,3),rho,dt*0.25d0,2,2)
       call muslcslv1D(fp2(:,:,4),rho,dt*0.25d0,1,2)
       endif
       !call BCtest()
       iwx=1;iwy=0
       !call BC(wp2(:,:,1),3,2,0,0,1)
       !call BC(wp2(:,:,2),2,3,0,0,2)
       !call BC(wp2(:,:,3),3,2,0,0,3)
       !call BC(wp2(:,:,4),2,3,0,0,4)
       if(itel==1) then
       call BC(wp2(:,:,1),3,3,0,0,1)
       call BC(wp2(:,:,2),3,3,0,0,2)
       call BC(wp2(:,:,3),3,3,0,0,3)
       call BC(wp2(:,:,4),3,3,0,0,4)
       endif

       if(iwv==1) then
       call BC(fp2(:,:,1),13,13,0,0,1)
       call BC(fp2(:,:,2),13,13,0,0,2)
       call BC(fp2(:,:,3),13,13,0,0,3)
       call BC(fp2(:,:,4),13,13,0,0,4)
       endif

       if(itel==1) then
       call muslcslv1D(wp2(:,:,1),rho,dt*0.25d0,1,2)
       call muslcslv1D(wp2(:,:,2),rho,dt*0.25d0,2,2)
       call muslcslv1D(wp2(:,:,3),rho,dt*0.25d0,1,2)
       call muslcslv1D(wp2(:,:,4),rho,dt*0.25d0,2,2)
       endif

       if(iwv==1) then
       call muslcslv1D(fp2(:,:,1),rho,dt*0.25d0,1,2)
       call muslcslv1D(fp2(:,:,2),rho,dt*0.25d0,2,2)
       call muslcslv1D(fp2(:,:,3),rho,dt*0.25d0,1,2)
       call muslcslv1D(fp2(:,:,4),rho,dt*0.25d0,2,2)
       endif
       !call BCtest()
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   wppre1(:,:,:)=wp1(:,:,:)
   wppre2(:,:,:)=wp2(:,:,:)


   !  if(mod(i,svnum)==1) then
   !     call saveu(sv)
   !  end if
  end do
  call saveu(sv)
end program muscl1D

subroutine INITIAL()
  use comvar
  use grvvar
  integer :: i,j
  double precision :: amp,pi=3.1415926535d0,haba,meanphi
  double precision :: rdmy=0.d0,rdmyx=0.d0,rdmyy=0.d0

  dinit1 = 2.0d0/G4pi/90.d0

  !----------x--------------
  dx = Lbox/dble(ndx-2)
  x(1) = dx/2.0d0
  x(0) = x(1) - dx
  x(-1) = x(0) - dx
  x(-2) = x(-1) - dx
  do i=2,ndx+1
     x(i) = x(i-1) + dx
  end do
  !----------x--------------

  !----------y--------------
  dy = Lbox/dble(ndy-2)
  y(1) = dy/2.0d0
  y(0) = y(1) - dy
  y(-1) = y(0) - dy
  y(-2) = y(-1) - dy
  do i=2,ndy+1
     y(i) = y(i-1) + dy
  end do
  !----------y--------------


  !---------Phi-------------
  !Phi(:,:)=0.0d0
  !---------Phi-------------


  !-------Phidtn-----------
  Phidtn(:,:)=0.0d0
  Phiexa(:,:)=0.d0
  Phiexa2(:,:)=0.d0
  !-------Phdtn-----------


  !----------wp--------------
  wp1(:,:,:)=0.d0
  wp2(:,:,:)=0.d0
  wppre1(:,:,:)=0.d0
  wppre2(:,:,:)=0.d0
  !----------wp--------------

  !----------fp--------------
  fp1(:,:,:)=0.d0
  fp2(:,:,:)=0.d0
  !----------fp--------------

   !----------f--------------
  fx(:,:)=0.d0
  fy(:,:)=0.d0
  rho(:,:)=0.d0
  !----------f--------------



  !---------rho-------------
  goto 2060
  do i = -1,ndx
     if( dabs(x(i) - hcen) .le. h) then
        rho(i,:) = dinit1
        !rho(i) = 0.0d0
     else
        rho(i,:) = 0.0d0
        !rho(i) = dinit1
        !rho(i) = dinit1*1.d-2
     end if
  end do

  meanrho=0.d0
  do j = 1,ndy-2
     do i = 1,ndx-2
        meanrho=meanrho+rho(i,j)
     end do
  end do
  meanrho=meanrho/dble(ndx-2)/dble(ndy-2)
  
  !rho(:,:) = rho(:,:) - meanrho

  do j = -1,ndy
     do i = -1,ndx
        rho(i,j)=rho(i,j)!-meanrho
     end do
  end do
  2060 continue
  !goto 2061

  meanrho=0.d0
  do i = -1-1,ndx+1
  do j = -1-1,ndy+1
     if( dsqrt((x(i) - hcen)**2+(y(j) - hcen)**2) .le. rsph) then
        rho(i,j) = dinit1
        !rho(i) = 0.0d0
        meanrho=meanrho+dinit1
     else
        rho(i,j) = 0.0d0
        !rho(i) = dinit1
        !rho(i) = dinit1*1.d-2
     end if
     !meanrho=meanrho!+rho(i,j)
  end do
  end do

  !do i = -1-1,ndx+1
  !do j = -1-1,ndy+1
  !if(i==j)then
  !Qgr(i,j)=(3.d0*(x(i)-hcen)*(y(j)-hcen)-(x(i)-hcen)**2-(y(j)-hcen)**2)*rho(i,j)*dx*dy
  !endif
  !if(i.ne.j) then
  !Qgr(i,j)=(3.d0*(x(i)-hcen)*(y(j)-hcen))*rho(i,j)*dx*dy
  !endif
  !end do
  !end do

!  meanrho=meanrho/dble(ndx+2)/dble(ndy+2)
  meanrho=0.d0
  !rho(:,:)=rho(:,:)-meanrho
!  do j = -1,ndy
!     do i = -1,ndx
!        rho(i,j)=rho(i,j)-meanrho
!     end do
!  end do
!  write(*,*)'meanrho',meanrho
  !2061 continue
  !---------rho-------------



  !--------Phiexa-----------
  !goto 200
  !dinit1=dinit1-meanrho
  meanphi=0.d0
  meanphiexa=0.d0
  open(142,file=dir//'phiexact.DAT')
  open(143,file=dir//'INIden.DAT')
  open(144,file=dir//'phigrd.DAT')
  goto 3060
!  do j= -1,ndy
!  do i= -1,ndx
!     if( dabs(x(i) - hcen) .le. h ) then
!        Phiexa(i,j) = G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2
!        write(142,*) sngl(x(i)) ,  sngl(G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2)
!        meanphi=meanphi+G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2
!     else
!        Phiexa(i,j) = G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
!        write(142,*) sngl(x(i)) , sngl(G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2)
!        meanphi=meanphi+G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
!     end if
!     write(143,*) sngl(rho(i,1))
!  end do
!  end do
!  meanphi=meanphi/dble(ndx+2)/dble(ndy+2)

  !Phiexa(:,:)=Phiexa(:,:)-meanphi
!  do j=-1,ndx
!  do i=0,ndx-1
!     Phigrd(i,j,1)=(-Phiexa(i-1,j)+Phiexa(i+1,j))*0.5d0/dx
     !write(144,*) sngl(x(i)) , Phigrd(i) , Phiexa(i-1),Phiexa(i+1)
!  end do
!  Phigrd(-1,j)=(-Phiexa(0,j)+Phiexa(1,j))/dx
!  Phigrd(ndx,j)=(Phiexa(ndx-1,j)-Phiexa(ndx-2,j))/dx
!  end do
  3060 continue

do j= -1,ndy
do i= -1,ndx
   if( dsqrt((x(i) - hcen)**2+(y(j) - hcen)**2) .le. rsph) then
      Phicrr(i,j) = G4pi/4.0d0 * rho(i,j) * ((x(i) - hcen)**2+(y(j) - hcen)**2)+Cnst !pi*G*rho*r^2
      write(142,*) sngl(x(i)),sngl(y(j)) ,  sngl(Phiexa(i,j))
      meanphi=meanphi+G4pi/2.0d0 * rho(i,j) * (x(i) - hcen )**2
   else
      Phicrr(i,j) = G4pi/2.0d0 * dinit1 * rsph **2 *dlog(dsqrt((x(i) - hcen)**2+(y(j) - hcen)**2)) &
           + ( G4pi/4.0d0 * dinit1 * rsph**2 - G4pi/2.0d0 * dinit1 * rsph**2 * dlog(rsph))+Cnst
      write(142,*) sngl(x(i)),sngl(y(j)) , sngl(Phiexa(i,j))
      meanphi=meanphi+G4pi * rho(i,j) * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
   end if
   write(143,*) sngl(x(i)),sngl(y(j)),sngl(rho(i,j))
end do
write(142,*)
write(143,*)
end do

dinit1=dinit1*rsph**3/(rch*rsph)**3
rsph=rch*rsph

  do j= -1-1,ndy+1
  do i= -1-1,ndx+1
     if( dsqrt((x(i) - hcen)**2+(y(j) - hcen)**2) .le. rsph) then
        Phiexa2(i,j) = G4pi/4.0d0 * rho(i,j) * ((x(i) - hcen)**2+(y(j) - hcen)**2)+Cnst !pi*G*rho*r^2
        !write(142,*) sngl(x(i)),sngl(y(j)) ,  sngl(Phiexa(i,j))
        meanphi=meanphi+G4pi/2.0d0 * rho(i,j) * (x(i) - hcen )**2
     else
        Phiexa2(i,j) = G4pi/2.0d0 * dinit1 * rsph **2 *dlog(dsqrt((x(i) - hcen)**2+(y(j) - hcen)**2)) &
             + ( G4pi/4.0d0 * dinit1 * rsph**2 - G4pi/2.0d0 * dinit1 * rsph**2 * dlog(rsph))+Cnst
        !write(142,*) sngl(x(i)),sngl(y(j)) , sngl(Phiexa(i,j))
        meanphi=meanphi+G4pi * rho(i,j) * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
     end if
     !write(143,*) sngl(x(i)),sngl(y(j)),sngl(rho(i,j))
  end do
  !write(142,*)
  !write(143,*)
  end do
  
 !Point source
 !do j= -1-1,ndy+1
 !do i= -1-1,ndx+1
 !Phiexa2(i,j) = -G4pi / 3.d0 * dinit1 * rsph**3.d0 / (dsqrt((x(i) - hcen)**2+(y(j) - hcen)**2+1.0d-10))!+1.d0
 !end do
 !end do

  !πGmλ/d
!  do j= -1-1,ndy+1
!  do i= -1-1,ndx+1
!Phiexa2(i,j) = 2.d0 * G4pi * 3.141592d0 * rsph * rsph *  dinit1 * dlog(dsqrt((x(i) - hcen)**2+(y(j) - hcen)**2+1.0d-10)) &
!-2.d0 * G4pi * 3.141592d0 * rsph * rsph *  dinit1 * dlog(dsqrt((0.d0 - hcen)**2+(0.d0 - hcen)**2+1.0d-10))!+1.d0
!Phiexa2(i,j)=G4pi/2.0d0 * dinit1 * rsph **2 *dlog(dsqrt((x(i) - hcen)**2+(y(j) - hcen)**2))
!  end do
!  end do

  !do j= -1-1,ndy+1
  !do i= -1-1,ndx+1
  !rdmy=dsqrt((x(i) - hcen)**2+(y(j) - hcen)**2)
  !rdmyx=x(i) - hcen
  !rdmyy=y(j) - hcen
  !Phiexa2(i,j) = -G4pi / 3.d0 * dinit1 * rsph**3.d0 / (rdmy**2+1.0d-10)- &
  !G4pi/8.d0/3.141592d0*(Qgr(i,i)*rdmyx*rdmyx/(rdmy**5)+Qgr(i,j)*rdmyx*rdmyy/(rdmy**5)&
  !+Qgr(j,i)*rdmyy*rdmyx/(rdmy**5)+Qgr(j,j)*rdmyy*rdmyy/(rdmy**5))
  !end do
  !end do
  
  do j=-1,ndy
  do i=-1,ndx
     Phigrd(i,j,1)= (-Phiexa2(i-1,j)+Phiexa2(i+1,j))*0.5d0/dx+(-Phiexa2(i,j-1)+Phiexa2(i,j+1))*0.5d0/dy
     Phigrd(i,j,2)=-(-Phiexa2(i-1,j)+Phiexa2(i+1,j))*0.5d0/dx-(-Phiexa2(i,j-1)+Phiexa2(i,j+1))*0.5d0/dy
     Phigrd(i,j,3)= (-Phiexa2(i-1,j)+Phiexa2(i+1,j))*0.5d0/dx-(-Phiexa2(i,j-1)+Phiexa2(i,j+1))*0.5d0/dy
     Phigrd(i,j,4)=-(-Phiexa2(i-1,j)+Phiexa2(i+1,j))*0.5d0/dx+(-Phiexa2(i,j-1)+Phiexa2(i,j+1))*0.5d0/dy
     Phiexa(i,j)=Phiexa2(i,j)
  end do
  end do
  do j=1,ndy-2
  do i=1,ndx-2
     meanphiexa=Phiexa(i,j)+meanphiexa
  end do
  end do
  meanphiexa=meanphiexa/dble(ndx-2)/dble(ndy-2)



  close(142)
  close(143)
  close(144)
  !200 continue
  !--------Phiexa-----------


  !---------wave--------
  goto 201
  !do i = -1, ndx
  !   amp = 1.d-3
  !   Phi(i) =  amp*dsin(2.d0*pi*x(i)/Lbox)
  !   Phi1step(i) =  amp*dsin(2.d0*pi*x(i)/Lbox)
  !end do


!  do i = -1, ndx
!     amp = 1.d-3
!     haba=10.0d0
     !Phi(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
!     Phi1step(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
!  end do
  201 continue
  !---------wave--------
end subroutine INITIAL

subroutine BCtest()
  use comvar
  use grvvar

wp1(:, 0,1)    =Phigrd(:, 0,1)
wp1(:,-1,1)    =Phigrd(:,-1,1)
wp1(:,ndy2,1)  =Phigrd(:,ndy2,1)
wp1(:,ndy2-1,1)=Phigrd(:,ndy2-1,1)
wp2(:, 0,1)    =Phiexa(:,0)
wp2(:,-1,1)    =Phiexa(:,-1)
wp2(:,ndy2  ,1)=Phiexa(:,ndy2)
wp2(:,ndy2-1,1)=Phiexa(:,ndy2-1)
wp1( 0,:,1)    =Phigrd( 0,:,1)
wp1(-1,:,1)    =Phigrd(-1,:,1)
wp1(ndx2  ,:,1)=Phigrd(ndx2,:,1)
wp1(ndx2-1,:,1)=Phigrd(ndx2-1,:,1)
wp2( 0,:,1)    =Phiexa(0,:)
wp2(-1,:,1)    =Phiexa(-1,:)
wp2(ndx2  ,:,1)=Phiexa(ndx2,:)
wp2(ndx2-1,:,1)=Phiexa(ndx2-1,:)
end subroutine BCtest


subroutine BC(U,modex0,modex1,modey0,modey1,numwp2)
  use comvar
!  use grvvar
  integer :: i,modex0,modey0,modex1,modey1,numwp2
  !double precision , dimension(1:2) :: pl,pr
  DOUBLE PRECISION , dimension(-1:ndx,-1:ndy) :: U
  
!**********periodic*******************
      if(modex0==1) then
             U( 0,:)    =U(ndx-2,:)
             U(-1,:)    =U(ndx-3,:)
       endif
       if(modex1==1) then
             U(ndx  ,:)=U(2,:)
             U(ndx-1,:)=U(1,:)
       end if
       if(modey0==1) then
              U(:, 0)    =U(:,ndy-2)
              U(:,-1)    =U(:,ndy-3)
       endif
       if(modey1==1) then
              U(:,ndy  )=U(:,2)
              U(:,ndy-1)=U(:,1)
       end if
!**********periodic*******************

!*************free********************
if(modex0==2)then
U(-1,:)    =U(0,:)
U( 0,:)    =U(1,:)
endif
if(modex1==2)then
U(ndx  ,:)=U(ndx-1,:)
U(ndx-1,:)=U(ndx-2,:)
endif
if(modey0==2)then
U(:,-1)    =U(:,0)
U(:, 0)    =U(:,1)
endif
if(modey1==2)then
U(:,ndy  )=U(:,ndy-1)
U(:,ndy-1)=U(:,ndy-2)
endif
!*************free********************

!************exact(wp->telegraph)********************
!%%%%%phi%%%%
if(modex0==3)then
    U( 0,:)    =Phiexa(0,:)
    U(-1,:)    =Phiexa(-1,:)
endif
if(modex1==3)then
    U(ndx  ,:)=Phiexa(ndx,:)
    U(ndx-1,:)=Phiexa(ndx-1,:)
endif
if(modey0==3)then
U(:, 0)    =Phiexa(:,0)
U(:,-1)    =Phiexa(:,-1)
endif
if(modey1==3)then
U(:,ndy  )=Phiexa(:,ndy  )
U(:,ndy-1)=Phiexa(:,ndy-1)
endif
!%%%%%phi%%%%
!%%%%%%f(pair of phi)%%%%%
if(modex0==4)then
  U( 0,:)    =cg*2.d0*1.d0/kappa/2.d0*Phigrd( 0,:,numwp2)+Phiexa( 0,:)
  U(-1,:)    =cg*2.d0*1.d0/kappa/2.d0*Phigrd(-1,:,numwp2)+Phiexa(-1,:)
endif
if(modex1==4)then
  U(ndx  ,:)=cg*2.d0*1.d0/kappa/2.d0*Phigrd(ndx  ,:,numwp2)+Phiexa(ndx  ,:)
  U(ndx-1,:)=cg*2.d0*1.d0/kappa/2.d0*Phigrd(ndx-1,:,numwp2)+Phiexa(ndx-1,:)
endif
if(modey0==4)then
   U(:, 0)    =cg*2.d0*1.d0/kappa/2.d0*Phigrd(:, 0,numwp2)+Phiexa(:, 0)
   U(:,-1)    =cg*2.d0*1.d0/kappa/2.d0*Phigrd(:,-1,numwp2)+Phiexa(:,-1)
endif
if(modey1==4)then
   U(:,ndy)  =cg*2.d0*1.d0/kappa/2.d0*Phigrd(:,ndy  ,numwp2)+Phiexa(:,ndy  )
   U(:,ndy-1)=cg*2.d0*1.d0/kappa/2.d0*Phigrd(:,ndy-1,numwp2)+Phiexa(:,ndy-1)
endif
!%%%%%%f(pair of phi)%%%%%
!************exact********************


!************exact(fp->wave)********************
!%%%%%phi%%%%
if(modex0==13)then
    U( 0,:)    =Phiexa(0,:)
    U(-1,:)    =Phiexa(-1,:)
endif
if(modex1==13)then
    U(ndx  ,:)=Phiexa(ndx,:)
    U(ndx-1,:)=Phiexa(ndx-1,:)
endif
if(modey0==13)then
U(:, 0)    =Phiexa(:,0)
U(:,-1)    =Phiexa(:,-1)
endif
if(modey1==13)then
U(:,ndy  )=Phiexa(:,ndy  )
U(:,ndy-1)=Phiexa(:,ndy-1)
endif
!%%%%%phi%%%%
!%%%%%%f(pair of phi)%%%%%
if(modex0==14)then
  U( 0,:)    =Phigrd( 0,:,numwp2)
  U(-1,:)    =Phigrd(-1,:,numwp2)
endif
if(modex1==14)then
  U(ndx  ,:)=Phigrd(ndx  ,:,numwp2)
  U(ndx-1,:)=Phigrd(ndx-1,:,numwp2)
endif
if(modey0==14)then
   U(:, 0)    =Phigrd(:, 0,numwp2)
   U(:,-1)    =Phigrd(:,-1,numwp2)
endif
if(modey1==14)then
   U(:,ndy)  =Phigrd(:,ndy  ,numwp2)
   U(:,ndy-1)=Phigrd(:,ndy-1,numwp2)
endif
!%%%%%%f%%%%%
!************exact********************
end subroutine BC


subroutine time(dt)
  use comvar
  use grvvar
  double precision :: dt
  dt = dx/cg * coeff
  write(*,*) 'time cg' , dt
end subroutine time



subroutine timesource(Phiv,source,dt,mode)
  use comvar
  !use grvver
  integer i,mode,j
  double precision :: dt,sdt,mindt,maxdt,epsl = 1.0d-4
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy) :: Phiv,source

  if(mode==3)then
     maxdt=0.d0
     mindt=1.d10
     do j=1,ndy-2
        do i=1,ndx-2
           if((source(i,j) .ne. 0.0d0) .and. (Phiv(i,j) .ne. 0.0d0))then
              sdt = dabs(0.5d0 * Phiv(i,j) /( G4pi * cg * cg  * source(i,j) ))
              mindt=dmin1(mindt,sdt)
           end if
        end do
     end do
     if( (mindt < dt) .and. (mindt .ne. 0.0d0)) then
        dt=dmin1(mindt,dt)
        !dt = sdt
     end if
  end if

  write(*,*) 'time source' , dt
end subroutine timesource


subroutine muslcslv1D(Phiv,source,dt,mode,hazi)
  use comvar
  double precision :: nu2 , w=6.0d0 , dt2 , dt , deltap,deltam  !kappa -> comver  better?
  integer :: direction , mode , invdt , loopmode , dloop,cnt=0
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy) :: Phigrad,Phipre,fluxphi&
       ,Phiv,Phi2dt,Phiu,sourcepre,sourcepri
  DOUBLE PRECISION, dimension(-1-1:ndx+1,-1-1:ndy+1) :: source
  character(5) name
  integer Ncell,Ncm,Ncl,ix,jy,kz,Lnum,Mnum,hazi,is,ie,idm
  !DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G

 if(iwx.eq.1) then; Ncell = ndx; Ncm = ndy; endif!  BT1 = 2; BT2 = 3; VN = 2; end if
    if(iwy.eq.1) then; Ncell = ndy; Ncm = ndx; endif! BT1 = 3; BT2 = 1; VN = 3; end if

  !----kyoukai-----
   if(hazi==1)then
      is = 2
      ie = Ncell-3
   end if
   if(hazi==2)then
      is = 1
      ie = Ncell-2
   end if
  !----kyoukai-----
  nu2 = cg * dt / dx
  Phipre(:,:) = Phiv(:,:)
  !write(name,'(i5.5)') cnt
  !------------ul.solver.+cg-------------
  if(mode==1) then
     call fluxcal(Phipre,Phipre,Phiu,0.0d0,1.d0/3.0d0,10,is,ie)
     !call fluxcal(Phipre,Phipre,Phiu,0.0d0,0.0d0,10)
     write(*,*) Phiu(1,1)
     !------------calcurate dt/2------------
!     DO Lnum = 1, Ncl-2
        DO Mnum = 1, Ncm-2
           do i = is-1,ie+1
              ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
              jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Mnum! + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
     !do i=ist-1,ndx-ien+1 !一次なので大丈夫
              Phi2dt(ix,jy) = Phipre(ix,jy)- 0.5d0 * nu2 * ( Phiu(ix,jy) - Phiu(ixm,jym))
           end do
        end DO
!     end DO
     !------------calcurate dt/2------------
     call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,1.d0/3.0d0,1,is,ie)
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,1)
     !write(*,*) Phiu(127),'127-2'
     !do i = ist , ndx-ien
!      DO Lnum = 1, Ncl-2
        DO Mnum = 1, Ncm-2
           do i = is,ie
              ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
              jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Mnum! + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              Phiv(ix,jy) = Phipre(ix,jy) - nu2 * (Phiu(ix,jy) - Phiu(ixm,jym))
           end do
        end DO
!     end DO
  end if
  !------------ul.solver.+cg-------------



  !------------ul.solver.-cg-------------
  if(mode==2) then

     call fluxcal(Phipre,Phipre,Phiu,0.0d0,1.d0/3.0d0,11,is,ie)
     !call fluxcal(Phipre,Phipre,Phiu,0.0d0,0.0d0,11)
     !------------calcurate dt/2------------
!     DO Lnum = 1, Ncl-2
        DO Mnum = 1, Ncm-2
           do i = is-1,ie+1
              ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
              jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Mnum !+ iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i=ist-1,ndx-ien+1
              Phi2dt(ix,jy) = Phipre(ix,jy) + 0.5d0 * nu2 * ( Phiu(ixp,jyp) - Phiu(ix,jy))
           end do
        end DO
!     end DO
     !------------calcurate dt/2------------
     call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,1.d0/3.0d0,4,is,ie)
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,4)

     !do i = ist , ndx-ien
!     DO Lnum = 1, Ncl-2
        DO Mnum = 1, Ncm-2
           do i = is,ie
              ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
              jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Mnum !+ iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              Phiv(ix,jy) = Phipre(ix,jy) + nu2 * (Phiu(ixp,jyp) - Phiu(ix,jy))
           end do
        end DO
!     end DO

     !do i=-1,ndx
     !   write(202,*) i, Phiv(i)
     !end do

  end if
  !------------ul.solver.-cg-------------


  !--------------source------------------
  if(mode==3) then
!     do k = 1,ndz-2
        do j = 1,ndy-2
           do i = 1,ndx-2
              !Phiv(i,j) =  cg * G4pi * source(i,j) * dt * phiratio + Phipre(i,j)
               Phiv(i,j) =  cg * source(i,j) * dt + Phipre(i,j)
           end do
        end do
!     end do
  end if

  if(mode==4) then
!     do k = 1,ndz-2
        do j = 1,ndy-2
           do i = 1,ndx-2
              !Phiv(i,j) = -cg * source(i,j) * dt + Phipre(i,j)
              Phiv(i,j) = -cg * source(i,j) * dt + Phipre(i,j)
           end do
        end do
!     end do
  end if
  !--------------source------------------

!  close(201)
!  close(202)
  cnt=cnt+2
end subroutine muslcslv1D

!subroutine vanalbada(fg,gradfg,iwx,iwy,iwz)
subroutine vanalbada(Mnum,Lnum,Phipre,Phigrad,i_sta,i_end,dmein)
  use comvar
  double precision :: delp , delm ,flmt,eps=1.0d-10
  !integer :: i , ip , im , flmt ,eps=1.0d-10
  integer :: Mnum,Lnum,Ncell,i_sta,i_end,k,dmein
  integer ix,jy,kz,ixp,jyp,kzp,ixm,jym,kzm
  integer :: i , ip , im
  !DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phigrad,Phipre
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy) :: Phipre
  DOUBLE PRECISION, dimension(-1:dmein) :: Phigrad


  !if(iwx.eq.1) Ncell = ndx
  !if(iwy.eq.1) Ncell = ndy
  !if(iwz.eq.1) Ncell = ndz

  do i = i_sta-1 , i_end+1
     ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
     jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!     kz  = iwx*Lnum + iwy*Mnum + iwz*i
     ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
     jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!     kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
     ixm = iwx*(i-1)+ iwy*Mnum! + iwz*Mnum
     jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!     kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)

     delp = Phipre(ixp,jyp)-Phipre(ix,jy)
     delm = Phipre(ix,jy)-Phipre(ixm,jym)
     flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
     !Phigrad(ix,jy,kz) = flmt
     Phigrad(i) = flmt
  end do

end subroutine vanalbada


subroutine fluxcal(preuse,pre,uin,ep,eta,mode,is,ie)
  use comvar
  double precision :: ep , eta
  DOUBLE PRECISION , dimension(-1:ndx,-1:ndy) :: ul,ur,pre,preuse,uin
  DOUBLE PRECISION , dimension(-1:ndx) :: slop  !------------- need allocation --------------
  integer :: i,mode,Ncell,Ncl,Ncm,j,k,Lnum,Mnum
  integer ix,jy,kz,ixp,jyp,kzp,ixm,jym,kzm,is,ie
  !DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G
  !uin(:)=0.0d0
!  if(iwx.eq.1) then; Ncell = ndx; Ncm = ndy; Ncl = ndz;  end if
!     if(iwy.eq.1) then; Ncell = ndy; Ncm = ndz; Ncl = ndx;  end if
!        if(iwz.eq.1) then; Ncell = ndz; Ncm = ndx; Ncl = ndy;  end if
  if(iwx.eq.1) then; Ncell = ndx; Ncm = ndy;  end if
     if(iwy.eq.1) then; Ncell = ndy; Ncm = ndx;  end if

           !call vanalbada(pre,slop)
           if(mode==1) then
!              DO Lnum = 1, Ncl-2
              DO Mnum = 1, Ncm-2
              call vanalbada(Mnum,Lnum,pre,slop,is,ie,Ncell)
              do i = is-1,ie+1
              ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
              jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Mnum !+ iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !call vanalbada(pre,slop)
              !do i = is,ie
              ul(ix,jy) = preuse(ix,jy) + 0.25d0 * ep * slop(i) &
                   * ((1.0d0-slop(i)*eta)*(pre(ix,jy)-pre(ixm,jym)) + &
                   (1.0d0+slop(i)*eta)*(pre(ixp,jyp) - pre(ix,jy))) !i+1/2
              uin(ix,jy)=ul(ix,jy)
              end do
              end DO
!              end DO
              !write(*,*) slop(127),'127slop'
              !uin(:)=ul(:)
           end if


           if(mode==4) then
!              DO Lnum = 1, Ncl-2
              DO Mnum = 1, Ncm-2
              call vanalbada(Mnum,Lnum,pre,slop,is,ie,Ncell)
              do i = is-1,ie+1
              ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
              jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Mnum! + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i = ist-1,ndx-ien+1
              ur(ix,jy) = preuse(ix,jy) - 0.25d0 * ep * slop(i) &
                   * ((1.0d0+slop(i)*eta)*(pre(ix,jy)-pre(ixm,jym)) + &
                   (1.0d0-slop(i)*eta)*(pre(ixp,jyp) - pre(ix,jy))) !i-1/2
              uin(ix,jy)=ur(ix,jy)
              end do
              end DO
!              end DO
              !write(*,*) slop(127),'127slop'
              !write(*,*) slop(ndx-ien),ndx-ien,slop(ndx-ien+1)
              !write(*,*) u(2)
              !uin(:)=ur(:)
           end if

           if(mode==10) then
!              DO Lnum = 1, Ncl-2
              DO Mnum = 1, Ncm-2
              !call vanalbada(pre,slop)
              do i = is-2,ie+2
              ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
              jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Mnum! + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i = ist-2,ndx-ien+2
              ul(ix,jy) = preuse(ix,jy)
              uin(ix,jy)=ul(ix,jy)
              end do
              end DO
!              end DO
           end if

           if(mode==11) then
!              DO Lnum = 1, Ncl-2
              DO Mnum = 1, Ncm-2
              !call vanalbada(pre,slop)
              do i = is-2,ie+2
              ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
              jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Mnum !+ iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i = is,ie
              ur(ix,jy) = preuse(ix,jy)
              uin(ix,jy)=ur(ix,jy)
              end do
              end DO
!              end DO
           end if


end subroutine fluxcal

subroutine saveu(in1)
  use comvar
  use grvvar
  integer :: i,in1,j
  character(5) name
  double precision meanphi1,meanphi2,errphi1,errphi2
  
  meanphi1=0.d0
  meanphi2=0.d0
  write(*,*) 'here2'!,in1,x(1),y(1),wp1(1,1),wp2(1,1),rho(1,1)
  write(name,'(i5.5)') in1
  open(21,file=dir//'phi2D'//name//'.dat')
  do j=-1,ndy
     do i=-1,ndx
      !  write(21,*) x(i),',',y(j),',',wp1(i,j,1),',',wp1(i,j,2),',',wp1(i,j,3),',',wp1(i,j,4) &
      !,',',wp2(i,j,1),',',wp2(i,j,2),',',wp2(i,j,3),',',wp2(i,j,4)&
      !,',',0.25d0*(wp2(i,j,1)+wp2(i,j,2)+wp2(i,j,3)+wp2(i,j,4)),',',rho(i,j),',',Phiexa(i,j),',',Phigrd(i,j,1),',',&
      !Phigrd(i,j,2),',',Phigrd(i,j,3),',',Phigrd(i,j,4),',',wp2(i,j,1)-Phiexa(i,j),',',wp1(i,j,1)-Phigrd(i,j,1)
        write(21,*) x(i),y(j),wp1(i,j,1),wp1(i,j,2),wp1(i,j,3),wp1(i,j,4) &
        ,wp2(i,j,1),wp2(i,j,2),wp2(i,j,3),wp2(i,j,4)&
        ,fp1(i,j,1),fp1(i,j,2),fp1(i,j,3),fp1(i,j,4) &
        ,fp2(i,j,1),fp2(i,j,2),fp2(i,j,3),fp2(i,j,4)&
        ,0.25d0*(wp2(i,j,1)+wp2(i,j,2)+wp2(i,j,3)+wp2(i,j,4)),rho(i,j),Phicrr(i,j),Phiexa(i,j),Phigrd(i,j,1),&
        Phigrd(i,j,2),Phigrd(i,j,3),Phigrd(i,j,4),wp2(i,j,1)-Phiexa(i,j),wp1(i,j,1)-Phigrd(i,j,1)
     end do
     write(21,*)
     meanphi1=0.25d0*(wp2(i,j,1)+wp2(i,j,2)+wp2(i,j,3)+wp2(i,j,4))+meanphi1
     meanphi2=wp2(i,j,1)+meanphi2
  end do
  close(21)
  !write(*,*) 'here3'
  open(22,file=dir//'xphi2D'//name//'.dat')
  j=ndy/2
  do i=1,ndx-2
        write(22,*) x(i),wp1(i,j,1),wp2(i,j,1),rho(i,j),Phiexa(i,j),Phigrd(i,j,1)!,wp2(j+1,k+1)-wp2(j+1,k-1)-wp2(j-1,k+1)+wp2(j-1,k-1)
  end do
  close(22)
  
  do j=1,ndy-2
     do i=1,ndx-2
     errphi1=(Phiexa(i,j)+meanphi1-meanphiexa)**2
     errphi2=(Phiexa(i,j)+meanphi2-meanphiexa)**2
     enddo
  end do
  errphi1=dsqrt(errphi1)
  errphi2=dsqrt(errphi2)
  
  open(850,file=dir//'err.dat',FORM='FORMATTED',position='append')
  write(850,*) in1, errphi1, errphi2, meanphiexa, meanphi1, meanphi2, errphi1/meanphi1, errphi2/meanphi2

  write(*,*) 'save step : ',in1
  in1=in1+1
end subroutine saveu

