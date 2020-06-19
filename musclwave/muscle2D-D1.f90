module comvar
  implicit none
  integer, parameter :: ndx=66,ndy=66,laststep=500,istx=1,ienx=2,isty=1,ieny=2,svnum=1 !preiodic:ist=1,ien=2 , kotei:ist=2,ien=3 : ndx=130
  !double precision, parameter :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0
  integer :: iwx,iwy,iwz
  DOUBLE PRECISION :: cg = 1.0d0 , dx,dy != Lbox/dble(ndx-2) !, bcphi1 , bcphi2
double precision :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0 ,rsph=5.d0
  !double precision :: G=1.11142d-4, G4pi=12.56637d0*G , coeff=0.90d0 ,  kappa=1.0d0/3.0d0
  double precision ::  G4pi=12.56637d0*1.11142d-4 , coeff=0.5d0 ,meanrho!,  kappa=1.0d0/3.0d0
  DOUBLE PRECISION , dimension(1:3) :: bcphi1 , bcphi2 ,bcphigrd1 , bcphigrd2
  character(51) :: dir='/Users/maeda/Desktop/Dropbox/analysis/muscle2Dtest/'
end module comvar

module grvvar
  implicit none
  integer, parameter :: ndx2=66 , ndy2=66 !パラメータ属性必要
  DOUBLE PRECISION , dimension(-1:ndx2) :: x
  DOUBLE PRECISION , dimension(-1:ndy2) :: y
  DOUBLE PRECISION , dimension(-1:ndx2,-1:ndy2) :: Phidtn , rho , fx , fy ,wp1 ,wp2
  DOUBLE PRECISION , dimension(-1:ndx2,-1:ndy2) :: Phidt,Phigrd,Phiexa
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
  call BC(1)
  !call muslcslv1D(Phi,Phi1step,dt,13)
  write(*,*) 'here'
  !call saveu(sv)
  do i=1,laststep
     call time(dt)
     !call timesource(Phidtn,rho,dt,3)
     write(*,*) i ,dt,'step'
     
    
    iwx=1;iwy=0
    call muslcslv1D(wp1,rho,dt,2,2)
    call BC(1)
    if(mod(i,svnum)==0) then
       call saveu(sv)
    end if
    iwx=0;iwy=1
    call muslcslv1D(wp1,rho,dt,2,2)
    call BC(2)
    do k=1,ndy-2
      do j=1,ndx-2
        wp1(j,k) = wp1(j,k)-cg*G4pi*rho(j,k)*dt &
                   -0.5d0*dt*cg*(wp2(j+1,k+1)-wp2(j+1,k-1)-wp2(j-1,k+1)+wp2(j-1,k-1))/dx/dy
      enddo
    enddo
    call BC(1)
    call BC(2)

    iwx=1;iwy=0
    call muslcslv1D(wp2,rho,dt,1,2)
    call BC(1)
    iwx=0;iwy=1
    call muslcslv1D(wp2,rho,dt,1,2)
    call BC(2)
    do k=1,ndy-2
      do j=1,ndx-2
        wp2(j,k) = wp2(j,k)+cg*dt*wp1(j,k)
      enddo
    enddo
    call BC(1)
    call BC(2)

     

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

  dinit1 = 2.0d0/G4pi/90.d0*1.d-10

  !----------x--------------
  dx = Lbox/dble(ndx-2)
  x(1) = dx/2.0d0
  x(0) = x(1) - dx
  x(-1) = x(0) - dx
  do i=2,ndx
     x(i) = x(i-1) + dx
  end do
  !----------x--------------

  !----------y--------------
  dy = Lbox/dble(ndy-2)
  y(1) = dy/2.0d0
  y(0) = y(1) - dy
  y(-1) = y(0) - dy
  do i=2,ndy
     y(i) = y(i-1) + dy
  end do
  !----------y--------------


  !---------Phi-------------
  !Phi(:,:)=0.0d0
  !---------Phi-------------


  !-------Phidtn-----------
  Phidtn(:,:)=0.0d0
  !-------Phdtn-----------


  !----------wp--------------
  wp1(:,:)=0.d0
  wp2(:,:)=0.d0
  !----------wp--------------

   !----------f--------------
  fx(:,:)=0.d0
  fy(:,:)=0.d0
  !----------f--------------



  !---------rho-------------
  !goto 2060
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

  do j = -1,ndy
     do i = -1,ndx
        rho(i,j)=rho(i,j)!-meanrho
     end do
  end do
  !2060 continue
  goto 2061
  do i = -1,ndx
  do j = -1,ndy
     if( dsqrt((x(i) - hcen)**2+(y(j) - hcen)**2) .le. rsph) then
        rho(i,j) = dinit1
        !rho(i) = 0.0d0
     else
        rho(i,j) = 0.0d0
        !rho(i) = dinit1
        !rho(i) = dinit1*1.d-2
     end if
  end do
  end do
  2061 continue
  !---------rho-------------



  !--------Phiexa-----------
  !goto 200
  !dinit1=dinit1-meanrho
  meanphi=0.d0
  open(142,file=dir//'phiexact.DAT')
  open(143,file=dir//'INIden.DAT')
  open(144,file=dir//'phigrd.DAT')
  do j= -1,ndy
  do i= -1,ndx
     if( dabs(x(i) - hcen) .le. h ) then
        Phiexa(i,j) = G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2
        write(142,*) sngl(x(i)) ,  sngl(G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2)
        meanphi=meanphi+G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2
     else
        Phiexa(i,j) = G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
        write(142,*) sngl(x(i)) , sngl(G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2)
        meanphi=meanphi+G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
     end if
     write(143,*) sngl(rho(i,1))
  end do
  end do
  meanphi=meanphi/dble(ndx+2)/dble(ndy+2)

  !Phiexa(:,:)=Phiexa(:,:)-meanphi

  do j=-1,ndx
  do i=0,ndx-1
     Phigrd(i,j)=(-Phiexa(i-1,j)+Phiexa(i+1,j))*0.5d0/dx
     !write(144,*) sngl(x(i)) , Phigrd(i) , Phiexa(i-1),Phiexa(i+1)
  end do
  Phigrd(-1,j)=(-Phiexa(0,j)+Phiexa(1,j))/dx
  Phigrd(ndx,j)=(Phiexa(ndx-1,j)-Phiexa(ndx-2,j))/dx
  end do
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



subroutine BC(mode)
  use comvar
  use grvvar
  integer :: i,mode
  double precision , dimension(1:2) :: pl,pr

  if(mode==1) then
    wp1( 0,:)    =wp1( 1,:)
    wp1(-1,:)    =wp1( 0,:)
    wp1(ndx2  ,:)=wp1(ndx2-1,:)
    wp1(ndx2-1,:)=wp1(ndx2-2,:)
    
    wp2( 0,:)    =wp2( 1,:)
    wp2(-1,:)    =wp2( 0,:)
    wp2(ndx2  ,:)=wp2(ndx2-1,:)
    wp2(ndx2-1,:)=wp2(ndx2-2,:)
  end if

  if(mode==2) then
        wp1(:, 0)    =wp1(:,ndy2-2)
        wp1(:,-1)    =wp1(:,ndy2-3)
        wp1(:,ndy2  )=wp1(:,2)
        wp1(:,ndy2-1)=wp1(:,1)
        
        wp2(:, 0)    =wp2(:,ndy2-2)
        wp2(:,-1)    =wp2(:,ndy2-3)
        wp2(:,ndy2  )=wp2(:,2)
        wp2(:,ndy2-1)=wp2(:,1)
     !wp1(:,-1)=0.d0
     !wp1(:,0)=0.d0
     !wp1(:,ndy2)=0.d0
     !wp1(:,ndy2-1)=0.d0
     
     !wp2(:,-1)=0.d0
     !wp2(:, 0)=0.d0
     !wp2(:,ndy2)=0.d0
     !wp2(:,ndy2-1)=0.d0
  end if

  if(mode==3)then
wp1( 0,:)    =Phigrd( 0,:)
wp1(-1,:)    =Phigrd(-1,:)
wp1(ndx2  ,:)=Phigrd(ndx2,:)
wp1(ndx2-1,:)=Phigrd(ndx2-1,:)

wp2( 0,:)    =Phiexa(0,:)
wp2(-1,:)    =Phiexa(-1,:)
wp2(ndx2  ,:)=Phiexa(ndx2,:)
wp2(ndx2-1,:)=Phiexa(ndx2-1,:)
endif
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
       ,Phiv,source,Phi2dt,Phiu,sourcepre,sourcepri
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


subroutine fluxcal(preuse,pre,uin,ep,kappa,mode,is,ie)
  use comvar
  double precision :: ep , kappa
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
                   * ((1.0d0-slop(i)*kappa)*(pre(ix,jy)-pre(ixm,jym)) + &
                   (1.0d0+slop(i)*kappa)*(pre(ixp,jyp) - pre(ix,jy))) !i+1/2
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
                   * ((1.0d0+slop(i)*kappa)*(pre(ix,jy)-pre(ixm,jym)) + &
                   (1.0d0-slop(i)*kappa)*(pre(ixp,jyp) - pre(ix,jy))) !i-1/2
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
  write(*,*) 'here2'!,in1,x(1),y(1),wp1(1,1),wp2(1,1),rho(1,1)
  write(name,'(i5.5)') in1
  open(21,file=dir//'phi'//name//'.dat')
  do j=1,ndy-2
     do i=1,ndx-2
        write(21,*) x(i),y(j),wp1(i,j),wp2(i,j),rho(i,j)
     end do
     write(21,*)
  end do
  close(21)
  write(*,*) 'here3'
  open(22,file=dir//'xphi'//name//'.dat')
  j=ndy/2
  do i=1,ndx-2
        write(22,*) x(i),wp1(i,j),wp2(i,j),rho(i,j),wp2(j+1,k+1)-wp2(j+1,k-1)-wp2(j-1,k+1)+wp2(j-1,k-1)
  end do
  close(22)
  write(*,*) 'save step : ',in1
  in1=in1+1
end subroutine saveu

