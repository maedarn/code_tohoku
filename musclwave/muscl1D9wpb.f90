module comvar
  implicit none
  integer, parameter :: ndx=66,ndy=66,ndz=66,laststep=2000,ist=2,ien=3 , istp=1,ienp=2!preiodic:ist=1,ien=2 , kotei:ist=2,ien=3 :: ndx=130
  !double precision, parameter :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0
  DOUBLE PRECISION :: cg = 1.0d0 , dx ,dy,dz!= Lbox/dble(ndx-2) !, bcphi1 , bcphi2
  integer :: iwx,iwy,iwz,ifEVO,ifEVO2 , svv=15
  double precision :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0
  !double precision :: G=1.11142d-4, G4pi=12.56637d0*G , coeff=0.90d0 ,  kappa=1.0d0/3.0d0
  double precision ::  G4pi=12.56637d0*1.11142d-4 , coeff=0.1d0 !,  kappa=1.0d0/3.0d0
  DOUBLE PRECISION , dimension(1:3,-1:ndy,-1:ndz) :: bcphi1 , bcphi2 ,bcphigrd1 , bcphigrd2
end module comvar

module grvvar
  implicit none
  integer, parameter :: ndx2=66,ndy2=66,ndz2=66 !パラメータ属性必要
  DOUBLE PRECISION , dimension(-1:ndx2) :: x,y,z
  DOUBLE PRECISION , dimension(-1:ndx2,-1:ndy2,-1:ndz2) ::  Phicgm ,rho, Phi1step , Phi2step ,Phicgp
  DOUBLE PRECISION , dimension(-1:ndx2,-1:ndy2,-1:ndz2) :: Phidt,Phigrd,Phiexa
end module grvvar

program muscl1D
  !implicit none まちがった位置
  use comvar
  use grvvar
  implicit none
  DOUBLE PRECISION :: dt=0.0d0
  integer :: i,sv=0,iws,ws=2


  call INITIAL()
  call cllsub(1,dt)
  !call muslcslv1D(Phi,Phi1step,dt,13)
  ifEVO=1
  ifEVO2=1
  do i=1,laststep
     write(*,*) i,'timestep'
     call cllsub(2,dt)
     call cllsub(3,dt)
     call muslcslv1D(Phi1step,rho,dt*0.5d0,3,2)
     call muslcslv1D(Phi2step,rho,dt*0.5d0,3,2)
     call cllsub(3,dt)
     call cllsub(4,dt)
     call cllsub(3,dt)
     call cllsub(1,dt)
     call muslcslv1D(Phicgp,Phi1step,dt,4,2)
     call muslcslv1D(Phicgm,Phi2step,dt,4,2)
     call cllsub(1,dt)
     call cllsub(5,dt)
     call cllsub(3,dt)
     call muslcslv1D(Phi1step,rho,dt*0.5d0,3,2)
     call muslcslv1D(Phi2step,rho,0.5d0*dt,3,2)
     call cllsub(3,dt)
     call cllsub(4,dt)
     ifEVO=1+mod(i,6)
     ifEVO2=1+mod(i,6)
     if(mod(i,svv)==0) then
        call saveu(sv)
     end if
  end do
end program muscl1D


subroutine cnbn(phi1,phi2)
  use comvar
  !double precision :: ep , kappa
  DOUBLE PRECISION , dimension(-1:ndx) :: phin,phi1,phi2
  integer :: i
  do i = 1 , ndx-2
     phin(i)=(Phi1(i)+phi2(i))*0.5d0
  end do
  phi1(:)=phin(:)
  phi2(:)=phin(:)
end subroutine cnbn

subroutine INITIAL()
  use comvar
  use grvvar
  integer :: i,j,k
  double precision :: amp,pi=3.1415926535d0,haba

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

  !----------z--------------
  dz = Lbox/dble(ndz-2)
  z(1) = dz/2.0d0
  z(0) = z(1) - dz
  z(-1) = z(0) - dz
  do i=2,ndz
     z(i) = z(i-1) + dz
  end do
  !----------z--------------


  !---------Phi-------------
  Phicgp(:,:,:)=0.0d0
  Phicgm(:,:,:)=0.0d0
  !---------Phi-------------

  !-------Phi1step-----------
  Phi1step(:,:,:)=0.0d0
  Phi2step(:,:,:)=0.0d0
  !-------Phi1step-----------

  !-------Phidt-----------
  !Phidt(:)=0.0d0
  !-------Phdt-----------




  !---------rho-------------
  do k=-1,ndz
     do j=-1,ndy
        do i = -1,ndx
           if( dabs(x(i) - hcen) .le. h) then
              !if( dabs(y(i) - hcen) .le. h) then
              !if( dabs(z(i) - hcen) .le. h) then
              rho(i,j,k) = dinit1
              !rho(i) = 0.0d0
           else
              rho(i,j,k) = 0.0d0
              !rho(i) = dinit1*1.d-2
           end if
        end do
     end do
  end do
  !---------rho-------------



  !--------Phiexa-----------
  !goto 200
  open(142,file='/Users/maeda/Desktop/kaiseki/testcode5/phiexact.DAT')
  open(143,file='/Users/maeda/Desktop/kaiseki/testcode5/INIden.DAT')
  open(144,file='/Users/maeda/Desktop/kaiseki/testcode5/phigrd.DAT')

  do k=-1,ndz
     do j=-1,ndy
        do i= -1,ndx
           if( dabs(x(i) - hcen) .le. h ) then
              !if( dabs(y(i) - hcen) .le. h) then
              !if( dabs(z(i) - hcen) .le. h) then
              Phiexa(i,j,k) = G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2
              write(142,*) sngl(x(i)),sngl(y(j)),sngl(z(k)) ,  sngl(G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2)
           else
              Phiexa(i,j,k) = G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
              write(142,*) sngl(x(i)) ,sngl(y(j)),sngl(z(k)) ,&
                   sngl(G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2)
           end if
           write(143,*) sngl(rho(i,j,k))
        end do
     end do
  end do

  !--- x-exact -----
  do k=-1,ndz
     do j=-1,ndy
        do i=0,ndx-1
           Phigrd(i,j,k)=(-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx
           !write(144,*) sngl(x(i)) , Phigrd(i) , Phiexa(i-1),Phiexa(i+1)
        end do
        Phigrd(-1,j,k)=(-Phiexa(0,j,k)+Phiexa(1,j,k))/dx
        Phigrd(ndx,j,k)=(Phiexa(ndx-1,j,k)-Phiexa(ndx-2,j,k))/dx
     end do
  end do
        !Phigrd(-1)=(-Phiexa(0)+Phiexa(1))/dx
        !Phigrd(ndx)=(Phiexa(ndx-1)-Phiexa(ndx-2))/dx

  do k=1,ndz-2
      do j=1,ndy-2
         do i=1,ndx-2
            write(144,*) sngl(x(i)),sngl(y(j)),sngl(z(k)) , sngl(Phigrd(i,j,k)) !, Phiexa(i-1,j,k),Phiexa(i+1)
         end do
      end do
   end do

   do k=-1,ndz
      do j=-1,ndy

         bcphi1(1,j,k) = G4pi * dinit1 * h * dabs(x(1) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
         bcphi2(1,j,k) = G4pi * dinit1 * h * dabs(x(ndx-2) - hcen)  - G4pi/2.0d0 * dinit1 * h**2

         bcphi1(2,j,k) = G4pi * dinit1 * h * dabs(x(0) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
         bcphi2(2,j,k) = G4pi * dinit1 * h * dabs(x(ndx-1) - hcen)  - G4pi/2.0d0 * dinit1 * h**2

         bcphi1(3,j,k) = G4pi * dinit1 * h * dabs(x(-1) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
         bcphi2(3,j,k) = G4pi * dinit1 * h * dabs(x(ndx) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
      end do
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


  do i = -1, ndx
     amp = 1.d-3
     haba=10.0d0
     !Phi(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
     !Phi1step(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
  end do
  201 continue
  !---------wave--------


  !--------const------------
  !---------Phi-------------
  Phicgp(:,:,:)=bcphi1(1,1,1)
  Phicgm(:,:,:)=bcphi1(1,1,1)
  !---------Phi-------------

  !-------Phidt-----------
  !Phidt(:)=bcphi1(1)
  !-------Phdt-----------
  !-------Phi1step-----------
  !Phi1step(:)=bcphi1
  !-------Phi1step-----------
  !--------const------------
end subroutine INITIAL



subroutine BC(mode)
  use comvar
  use grvvar
  integer :: i,mode,j,k
  double precision , dimension(1:2) :: pl,pr

  if(mode==1) then
     !---------kotei-x-----------
     !goto 100
     !---------Phi-------------
     Phicgp(1,:,:)= bcphi1(1,:,:)
     Phicgp(0,:,:)= bcphi1(2,:,:)
     Phicgp(-1,:,:)= bcphi1(3,:,:)
     Phicgp(ndx-2,:,:)= bcphi2(1,:,:)
     Phicgp(ndx-1,:,:)= bcphi2(2,:,:)
     Phicgp(ndx,:,:)= bcphi2(3,:,:)

     Phicgm(1,:,:)= bcphi1(1,:,:)
     Phicgm(0,:,:)= bcphi1(2,:,:)
     Phicgm(-1,:,:)= bcphi1(3,:,:)
     Phicgm(ndx-2,:,:)= bcphi2(1,:,:)
     Phicgm(ndx-1,:,:)= bcphi2(2,:,:)
     Phicgm(ndx,:,:)= bcphi2(3,:,:)

     !---------Phi-------------

  end if

  if(mode==3)then
     !---------kotei-x-----------
     !-------Phi1step+cg-----------
     !goto 700
     Phi1step(1,:,:)= Phigrd(1,:,:)
     Phi1step(0,:,:)= Phigrd(0,:,:)
     Phi1step(-1,:,:)=Phigrd(-1,:,:)
     Phi1step(ndx-2,:,:)= Phigrd(ndx-2,:,:)
     Phi1step(ndx-1,:,:)= Phigrd(ndx-1,:,:)
     Phi1step(ndx,:,:)= Phigrd(ndx,:,:)
     !700 continue
     !-------Phi1step-----------
     !---------kotei-x-----------
  end if

  if(mode==4) then
     !---------kotei-x-----------
     !-------Phi1step-cg-----------
     !goto 701
     Phi2step(1,:,:)= -Phigrd(1,:,:)
     Phi2step(0,:,:)= -Phigrd(2,:,:)
     Phi2step(-1,:,:)=-Phigrd(3,:,:)
     Phi2step(ndx-2,:,:)= -Phigrd(ndx-2,:,:)
     Phi2step(ndx-1,:,:)= -Phigrd(ndx-1,:,:)
     Phi2step(ndx,:,:)= -Phigrd(ndx,:,:)
     !701 continue
     !-------Phi1step-----------
     !---------kotei-x-----------
  end if


  !----- y-periodic ------
  if(mode==8) then
     !-------period2-----------
     !     goto 102
     do k=-1,ndz
        do i=-1,ndx
           !---------Phi-------------
           pr(2)= Phicgp(i,ndx-2,k)
           pr(1)= Phicgp(i,ndx-3,k)
           pl(1)= Phicgp(i,1,k)
           pl(2)= Phicgp(i,2,k)
           Phicgp(i,-1,k)=pr(1)
           Phicgp(i,0,k)=pr(2)
           Phicgp(i,ndx-1,k)=pl(1)
           Phicgp(i,ndx,k)=pl(2)
           !---------Phi-------------
        end do
     end do
  end if
  if(mode==18) then
     do k=-1,ndz
        do i=-1,ndx
           !-------Phi1step-----------
           pr(2)= Phi1step(i,ndx-2,k)
           pr(1)= Phi1step(i,ndx-3,k)
           pl(1)= Phi1step(i,1,k)
           pl(2)= Phi1step(i,2,k)
           Phi1step(i,-1,k)=pr(1)
           Phi1step(i,0,k)=pr(2)
           Phi1step(i,ndx-1,k)=pl(1)
           Phi1step(i,ndx,k)=pl(2)
           !-------Phi1step-----------
!102  continue
           !-------period2-----------
        end do
     end do
  end if




  if(mode==28) then
     !-------period2-----------
     !     goto 102
     do k=-1,ndz
        do i=-1,ndx
           !---------Phi-------------
           pr(2)= Phicgm(i,ndx-2,k)
           pr(1)= Phicgm(i,ndx-3,k)
           pl(1)= Phicgm(i,1,k)
           pl(2)= Phicgm(i,2,k)
           Phicgm(i,-1,k)=pr(1)
           Phicgm(i,0,k)=pr(2)
           Phicgm(i,ndx-1,k)=pl(1)
           Phicgm(i,ndx,k)=pl(2)
           !---------Phi-------------
        end do
     end do
  end if
  if(mode==38) then
     do k=-1,ndz
        do i=-1,ndx
           !-------Phi1step-----------
           pr(2)= Phi2step(i,ndx-2,k)
           pr(1)= Phi2step(i,ndx-3,k)
           pl(1)= Phi2step(i,1,k)
           pl(2)= Phi2step(i,2,k)
           Phi2step(i,-1,k)=pr(1)
           Phi2step(i,0,k)=pr(2)
           Phi2step(i,ndx-1,k)=pl(1)
           Phi2step(i,ndx,k)=pl(2)
           !-------Phi1step-----------
!102  continue
           !-------period2-----------
        end do
     end do
  end if


  !----- y-periodic ------


  !----- z-periodic ------
  if(mode==9) then
     !-------period2-----------
     !     goto 102
     do j=-1,ndy
        do i=-1,ndx
           !---------Phi-------------
           pr(2)= Phicgp(i,j,ndx-2)
           pr(1)= Phicgp(i,j,ndx-3)
           pl(1)= Phicgp(i,j,1)
           pl(2)= Phicgp(i,j,2)
           Phicgp(i,j,-1)=pr(1)
           Phicgp(i,j,0)=pr(2)
           Phicgp(i,j,ndx-1)=pl(1)
           Phicgp(i,j,ndx)=pl(2)
           !---------Phi-------------
        end do
     end do
  end if
  if(mode==19) then
     do j=-1,ndy
        do i=-1,ndx
           !-------Phi1step-----------
           pr(2)= Phi1step(i,j,ndx-2)
           pr(1)= Phi1step(i,j,ndx-3)
           pl(1)= Phi1step(i,j,1)
           pl(2)= Phi1step(i,j,2)
           Phi1step(i,j,-1)=pr(1)
           Phi1step(i,j,0)=pr(2)
           Phi1step(i,j,ndx-1)=pl(1)
           Phi1step(i,j,ndx)=pl(2)
           !-------Phi1step-----------
!102  continue
           !-------period2-----------
        end do
     end do
  end if




  if(mode==29) then
     !-------period2-----------
     !     goto 102
     do j=-1,ndy
        do i=-1,ndx
           !---------Phi-------------
           pr(2)= Phicgm(i,j,ndx-2)
           pr(1)= Phicgm(i,j,ndx-3)
           pl(1)= Phicgm(i,j,1)
           pl(2)= Phicgm(i,j,2)
           Phicgm(i,j,-1)=pr(1)
           Phicgm(i,j,0)=pr(2)
           Phicgm(i,j,ndx-1)=pl(1)
           Phicgm(i,j,ndx)=pl(2)
           !---------Phi-------------
        end do
     end do
  end if
  if(mode==39) then
     do j=-1,ndy
        do i=-1,ndx
           !-------Phi1step-----------
           pr(2)= Phi2step(i,j,ndx-2)
           pr(1)= Phi2step(i,j,ndx-3)
           pl(1)= Phi2step(i,j,1)
           pl(2)= Phi2step(i,j,2)
           Phi2step(i,j,-1)=pr(1)
           Phi2step(i,j,0)=pr(2)
           Phi2step(i,j,ndx-1)=pl(1)
           Phi2step(i,j,ndx)=pl(2)
           !-------Phi1step-----------
!102  continue
           !-------period2-----------
        end do
     end do
  end if


  !----- z-periodic ------
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
  integer i,mode,j,k
  double precision :: dt,sdt,mindt,maxdt , epsl = 1.0d-4
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phiv,source

  !mindt=1000.0d0
  maxdt=0.0d0

  if(mode==1) then
     do k=1,ndx-2
        do j=1,ndx-2
           do i=1,ndx-2
              if((source(i,j,k) .ne. 0.0d0) .and. (Phiv(i,j,k) .ne. 0.0d0))then
                 sdt = 0.5d0*dabs(Phiv(i,j,k)) / (cg * G4pi * source(i,j,k) )
                 !sdt = 0.2d0*dabs(Phiv(i)) / (cg * G4pi * source(i) )
                 !mindt=dmin1(mindt,sdt)
                 maxdt=dmax1(maxdt,sdt)
              end if
           end do
        end do
     end do
     if( (maxdt < dt) .and. (maxdt .ne. 0.0d0)) then
        dt = sdt
     end if
  end if


  if(mode==2) then
     do k=1,ndx-2
        do j=1,ndx-2
           do i=1,ndx-2
              if((source(i,j,k) .ne. 0.0d0) .and. (Phiv(i,j,k) .ne. 0.0d0))then
                 sdt = 0.5d0*dabs(Phiv(i,j,k)) / ( cg * source(i,j,k) )
                 !sdt = 0.05d0*dabs(Phiv(i)) / ( cg * source(i) )
                 !mindt=dmin1(mindt,sdt)
                 maxdt=dmax1(maxdt,sdt)
              end if
           end do
        end do
     end do
     write(*,*) maxdt,'maxdt'
     if( (maxdt < dt) .and. (maxdt .ne. 0.0d0)) then
        dt = sdt
     end if
  end if


  write(*,*) 'time source' , dt
end subroutine timesource


subroutine timesource2(Phiv,Phidt,source,dt)
  use comvar
  !use grvver
  integer i,mode
  double precision :: dt,sdt,mindt,maxdt
  DOUBLE PRECISION, dimension(-1:ndx) :: Phiv,source,Phidt
  do i=1,ndx-2
     if((source(i) .ne. 0.0d0) .and. (Phiv(i) .ne. 0.0d0))then
        sdt = 0.5d0*dsqrt( dabs(2*Phiv(i)-Phidt(i)) / (cg * cg * G4pi * source(i) ) )
        !sdt = 0.2d0*dabs(Phiv(i)) / (cg * G4pi * source(i) )
        !mindt=dmin1(mindt,sdt)
        maxdt=dmax1(maxdt,sdt)
     end if
  end do
  if( (maxdt < dt) .and. (maxdt .ne. 0.0d0)) then
     dt = sdt
  end if


  write(*,*) 'time source2' , dt
end subroutine timesource2


subroutine muslcslv1D(Phiv,source,dt,mode,hazi)
  use comvar
  double precision :: nu2 , w=6.0d0 , dt2 , dt , deltap,deltam !kappa -> comver  better?
  integer :: direction , mode , invdt , loopmode , dloop,cnt=0
  !DOUBLE PRECISION :: fluxf(-1:ndx,-1:ndy,-1:ndz),fluxg(-1:ndx,-1:ndy,-1:ndz)
  !DOUBLE PRECISION, dimension(-1:ndx) :: Phigrad,Phipre,fluxphi,Phiv,source,Phi2dt,Phiu,sourcepre,sourcepri
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phigrad,Phipre,fluxphi&
       ,Phiv,source,Phi2dt,Phiu,sourcepre,sourcepri
  character(5) name
  integer Ncell,Ncm,Ncl,ix,jy,kz,Lnum,Mnum,hazi,is,ie


  if(iwx.eq.1) then; Ncell = ndx; Ncm = ndy; Ncl = ndz; endif!  BT1 = 2; BT2 = 3; VN = 2; end if
     if(iwy.eq.1) then; Ncell = ndy; Ncm = ndz; Ncl = ndx; endif! BT1 = 3; BT2 = 1; VN = 3; end if
        if(iwz.eq.1) then; Ncell = ndz; Ncm = ndx; Ncl = ndy; endif! BT1 = 1; BT2 = 2; VN = 4; end if


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
  Phipre(:,:,:) = Phiv(:,:,:)
  !write(name,'(i5.5)') cnt
  !------------ul.solver.+cg-------------
  if(mode==1) then
     call fluxcal(Phipre,Phipre,Phiu,0.0d0,1.d0/3.0d0,10,is,ie)
     !call fluxcal(Phipre,Phipre,Phiu,0.0d0,0.0d0,10)
     !------------calcurate dt/2------------
     DO Lnum = 1, Ncl
        DO Mnum = 1, Ncm
           do i = is-1,ie+1
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
     !do i=ist-1,ndx-ien+1 !一次なので大丈夫
              Phi2dt(ix,jy,kz) = Phipre(ix,jy,kz)- 0.5d0 * nu2 * ( Phiu(ix,jy,kz) - Phiu(ixm,jym,kzm))
           end do
        end DO
     end DO
     !------------calcurate dt/2------------
     call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,1.d0/3.0d0,1,is,ie)
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,1)
     !write(*,*) Phiu(127),'127-2'
     !do i = ist , ndx-ien
      DO Lnum = 1, Ncl
        DO Mnum = 1, Ncm
           do i = is,ie
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              Phiv(ix,jy,kz) = Phipre(ix,jy,kz) - nu2 * (Phiu(ix,jy,kz) - Phiu(ixm,jym,kzm))
           end do
        end DO
     end DO
  end if
  !------------ul.solver.+cg-------------



  !------------ul.solver.-cg-------------
  if(mode==2) then

     call fluxcal(Phipre,Phipre,Phiu,0.0d0,1.d0/3.0d0,11,is,ie)
     !call fluxcal(Phipre,Phipre,Phiu,0.0d0,0.0d0,11)
     !------------calcurate dt/2------------
     DO Lnum = 1, Ncl
        DO Mnum = 1, Ncm
           do i = is-1,ie+1
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
     !do i=ist-1,ndx-ien+1
              Phi2dt(ix,jy,kz) = Phipre(ix,jy,kz) + 0.5d0 * nu2 * ( Phiu(ixp,jyp,kzp) - Phiu(ix,jy,kz))
           end do
        end DO
     end DO
     !------------calcurate dt/2------------
     call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,1.d0/3.0d0,4,is,ie)
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,4)

     !do i = ist , ndx-ien
     DO Lnum = 1, Ncl
        DO Mnum = 1, Ncm
           do i = is,ie
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              Phiv(ix,jy,kz) = Phipre(ix,jy,kz) + nu2 * (Phiu(ixp,jyp,kzp) - Phiu(ix,jy,kz))
           end do
        end DO
     end DO

     !do i=-1,ndx
     !   write(202,*) i, Phiv(i)
     !end do

  end if
  !------------ul.solver.-cg-------------


  !--------------source------------------
  if(mode==3) then
     do k = 1,ndz-2
        do j = 1,ndy-2
           do i=is,ie
              Phiv(i,j,k) =  -cg * G4pi * source(i,j,k) * dt + Phipre(i,j,k)
           end do
        end do
     end do
  end if

  if(mode==4) then
     do k = 1,ndz-2
        do j = 1,ndy-2
           do i=is,ie
              Phiv(i,j,k) = cg * source(i,j,k) * dt + Phipre(i,j,k)
           end do
        end do
     end do
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
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phipre
  DOUBLE PRECISION, dimension(-1:dmein) :: Phigrad


  if(iwx.eq.1) Ncell = ndx
  if(iwy.eq.1) Ncell = ndy
  if(iwz.eq.1) Ncell = ndz

  do i = i_sta-1 , i_end+1
     ix  = iwx*i    + iwy*Lnum + iwz*Mnum
     jy  = iwx*Mnum + iwy*i    + iwz*Lnum
     kz  = iwx*Lnum + iwy*Mnum + iwz*i
     ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
     jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
     kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
     ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
     jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
     kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)

     delp = Phipre(ixp,jyp,kzp)-Phipre(ix,jy,kz)
     delm = Phipre(ix,jy,kz)-Phipre(ixm,jym,kzm)
     flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
     !Phigrad(ix,jy,kz) = flmt
     Phigrad(i) = flmt
  end do

end subroutine vanalbada

subroutine saveu(in1)
  use comvar
  use grvvar
  integer :: i,in1,j,k
  character(5) name

  write(name,'(i5.5)') in1
  open(21,file='/Users/maeda/Desktop/kaiseki/testcode5/phi'//name//'.dat')
  !open(23,file='/Users/maeda/Desktop/kaiseki/testcode5/phix'//name//'.dat')
  !open(24,file='/Users/maeda/Desktop/kaiseki/testcode5/phiy'//name//'.dat')
  !open(25,file='/Users/maeda/Desktop/kaiseki/testcode5/phiz'//name//'.dat')
  !goto 399
  do k=1,ndz-2
     do j=1,ndy-2
        do i=1,ndx-2
           write(21,*) sngl(x(i)),sngl(y(j)),sngl(z(k)),sngl(Phicgp(i,j,k)),sngl(Phi1step(i,j,k))&
                , sngl(Phicgm(i,j,k)),sngl(Phi2step(i,j,k)),sngl((Phicgp(i,j,k)+Phicgm(i,j,k))*0.5d0)
           !write(21,*) x(i), Phicgp(i),Phi1step(i) , Phicgm(i),Phi2step(i) ,&
           !(Phicgp(i)+Phicgm(i))*0.5d0,(Phi1step(i)+Phi2step(i))*0.5d0, (Phicgp(i)-Phicgm(i)),&
           !(Phicgp(i)+Phicgm(i))*0.5d0+dabs(Phicgp(i)-Phicgm(i))
        end do
     end do
  end do
  !399 continue
  close(21)

  !do i=1,ndx-2
  !  write(23,*) sngl(x(i))&!, sngl(Phicgp(i,30,30)),sngl(Phi1step(i,30,30))&
  !       , sngl(Phicgm(i,30,30)),sngl(Phi2step(i,30,30))!,sngl((Phicgp(i,30,30)+Phicgm(i,30,30))*0.5d0)
 !end do
 !close(23)
 !do i=1,ndx-2
 !   write(24,*) sngl(y(i))&!, sngl(Phicgp(30,i,30)),sngl(Phi1step(30,i,30))&
!         , sngl(Phicgm(30,i,30)),sngl(Phi2step(30,i,30))!,sngl((Phicgp(30,i,30)+Phicgm(30,i,30))*0.5d0)
! end do
! close(24)
! do i=1,ndx-2
!    write(25,*) sngl(z(i))&!, sngl(Phicgp(30,30,i)),sngl(Phi1step(30,30,i))&
!         , sngl(Phicgm(30,30,i)),sngl(Phi2step(30,30,i))!,sngl((Phicgp(30,30,i)+Phicgm(30,30,i))*0.5d0)
! end do
! close(25)
  in1=in1+1
end subroutine saveu



subroutine fluxcal(preuse,pre,u,ep,kappa,mode,is,ie)
  use comvar
  double precision :: ep , kappa
  DOUBLE PRECISION , dimension(-1:ndx,-1:ndy,-1:ndz) :: ul,ur,pre,preuse,u
  DOUBLE PRECISION , dimension(-1:ndx) :: slop
  integer :: i,mode,Ncell,Ncl,Ncm,j,k,Lnum,Mnum
  integer ix,jy,kz,ixp,jyp,kzp,ixm,jym,kzm,is,ie
  !u(:)=0.0d0
  if(iwx.eq.1) then; Ncell = ndx; Ncm = ndy; Ncl = ndz;  end if
     if(iwy.eq.1) then; Ncell = ndy; Ncm = ndz; Ncl = ndx;  end if
        if(iwz.eq.1) then; Ncell = ndz; Ncm = ndx; Ncl = ndy;  end if

           !call vanalbada(pre,slop)
           if(mode==1) then
              DO Lnum = 1, Ncl
              DO Mnum = 1, Ncm
              call vanalbada(Ncm,Ncl,pre,slop,is,ie,Ncell)
              do i = is-1,ie+1
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !call vanalbada(pre,slop)
              !do i = is,ie
              ul(ix,jy,kz) = preuse(ix,jy,kz) + 0.25d0 * ep * slop(i) &
                   * ((1.0d0-slop(i)*kappa)*(pre(ix,jy,kz)-pre(ixm,jym,kzm)) + &
                   (1.0d0+slop(i)*kappa)*(pre(ixp,jyp,kzp) - pre(ix,jy,kz))) !i+1/2
              u(ix,jy,kz)=ul(ix,jy,kz)
              end do
              end DO
              end DO
              !write(*,*) slop(127),'127slop'
              !u(:)=ul(:)
           end if


           if(mode==4) then
              DO Lnum = 1, Ncl
              DO Mnum = 1, Ncm
              call vanalbada(Ncm,Ncl,pre,slop,is,ie,Ncell)
              do i = is-1,ie+1
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i = ist-1,ndx-ien+1
              ur(ix,jy,kz) = preuse(ix,jy,kz) - 0.25d0 * ep * slop(i) &
                   * ((1.0d0+slop(i)*kappa)*(pre(ix,jy,kz)-pre(ixm,jym,kzm)) + &
                   (1.0d0-slop(i)*kappa)*(pre(ixp,jyp,kzp) - pre(ix,jy,kz))) !i-1/2
              u(ix,jy,kz)=ur(ix,jy,kz)
              end do
              end DO
              end DO
              !write(*,*) slop(127),'127slop'
              !write(*,*) slop(ndx-ien),ndx-ien,slop(ndx-ien+1)
              !write(*,*) u(2)
              !u(:)=ur(:)
           end if

           if(mode==10) then
              DO Lnum = 1, Ncl
              DO Mnum = 1, Ncm
              !call vanalbada(pre,slop)
              do i = is-2,ie+2
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i = ist-2,ndx-ien+2
              ul(ix,jy,kz) = preuse(ix,jy,kz)
              u(ix,jy,kz)=ul(ix,jy,kz)
              end do
              end DO
              end DO
           end if

           if(mode==11) then
              DO Lnum = 1, Ncl
              DO Mnum = 1, Ncm
              !call vanalbada(pre,slop)
              do i = is-2,ie+2
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i = is,ie
              ur(ix,jy,kz) = preuse(ix,jy,kz)
              u(ix,jy,kz)=ur(ix,jy,kz)
              end do
              end DO
              end DO
           end if

           !if(mode==12) then
           !   do i = ist-1,ndx-ien+1
           !      ur(i) = preuse(i+1)
           !      u(i)=ur(i)
           !   end do
           !end if

end subroutine fluxcal





subroutine cllsub(mode,dt)
  use comvar
  use grvvar
  integer :: mode !,ifEVO=1,ifEVO2=1
  double precision dt
  if(mode==1) then
     call BC(1)
     call BC(8)
     call BC(28)
     call BC(9)
     call BC(29)
  end if

  if(mode==2) then
     call time(dt)
     call timesource(Phicgp,Phi1step,dt,2)
     call timesource(Phi1step,rho,dt,1)
     call timesource(Phicgm,Phi2step,dt,2)
     call timesource(Phi2step,rho,dt,1)
  end if

  if(mode==3) then
     call BC(4)
     call BC(3)
     call BC(18)
     call BC(38)
     call BC(19)
     call BC(39)
  end if

  if(mode==4) then
     if(ifEVO.eq.1) then
        iwx=1; iwy=0; iwz=0
        call BC(3)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        call BC(4)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        iwx=0; iwy=1; iwz=0
        call BC(18)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BC(38)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=0; iwy=0; iwz=1
        call BC(19)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BC(39)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        !ifEVO = 2
        goto 1000
     end if
     if(ifEVO.eq.2) then
        iwx=0; iwy=1; iwz=0
        call BC(18)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BC(38)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=0; iwy=0; iwz=1
        call BC(19)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BC(39)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=1; iwy=0; iwz=0
        call BC(3)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        call BC(4)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        !ifEVO = 3
        goto 1000
     end if
     if(ifEVO.eq.3) then
        iwx=0; iwy=0; iwz=1
        call BC(19)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BC(39)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=1; iwy=0; iwz=0
        call BC(3)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        iwx=0; iwy=1; iwz=0
        call BC(18)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BC(38)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        !ifEVO = 4
        goto 1000
     end if
     if(ifEVO.eq.4) then
        iwx=1; iwy=0; iwz=0
        call BC(3)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        call BC(4)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        iwx=0; iwy=0; iwz=1
        call BC(19)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BC(39)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=0; iwy=1; iwz=0
        call BC(18)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BC(38)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        !ifEVO = 5
        goto 1000
     end if
     if(ifEVO.eq.5) then
        iwx=0; iwy=1; iwz=0
        call BC(18)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BC(38)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=1; iwy=0; iwz=0
        call BC(3)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        call BC(4)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        iwx=0; iwy=0; iwz=1
        call BC(19)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BC(39)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        !ifEVO = 6
        goto 1000
     end if
     if(ifEVO.eq.6) then
        iwx=0; iwy=0; iwz=1
        call BC(19)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BC(39)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=0; iwy=1; iwz=0
        call BC(18)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BC(38)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=1; iwy=0; iwz=0
        call BC(3)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        call BC(4)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        !ifEVO = 1
        goto 1000
     end if
1000 continue
  end if


  if(mode==5) then
     if(ifEVO2.eq.1) then
        iwx=1; iwy=0; iwz=0
        call BC(1)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call BC(1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        iwx=0; iwy=1; iwz=0
        call BC(8)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BC(28)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=0; iwy=0; iwz=1
        call BC(9)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BC(29)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        !ifEVO2 = 2
        goto 1100
     end if
     if(ifEVO2.eq.2) then
        iwx=0; iwy=1; iwz=0
        call BC(8)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BC(28)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=0; iwy=0; iwz=1
        call BC(9)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BC(29)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=1; iwy=0; iwz=0
        call BC(1)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call BC(1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        !ifEVO2 = 3
        goto 1100
     end if
     if(ifEVO2.eq.3) then
        iwx=0; iwy=0; iwz=1
        call BC(9)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BC(29)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=1; iwy=0; iwz=0
        call BC(1)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call BC(1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        iwx=0; iwy=1; iwz=0
        call BC(8)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BC(28)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        !ifEVO2 = 4
        goto 1100
     end if
     if(ifEVO2.eq.4) then
        iwx=1; iwy=0; iwz=0
        call BC(1)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call BC(1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        iwx=0; iwy=0; iwz=1
        call BC(9)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BC(29)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=0; iwy=1; iwz=0
        call BC(8)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BC(28)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        !ifEVO2 = 5
        goto 1100
     end if
     if(ifEVO2.eq.5) then
        iwx=0; iwy=1; iwz=0
        call BC(8)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BC(28)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=1; iwy=0; iwz=0
        call BC(1)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call BC(1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        iwx=0; iwy=0; iwz=1
        call BC(9)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BC(29)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        !ifEVO2 = 6
        goto 1100
     end if
     if(ifEVO2.eq.6) then
        iwx=0; iwy=0; iwz=1
        call BC(9)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BC(29)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=0; iwy=1; iwz=0
        call BC(8)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BC(28)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=1; iwy=0; iwz=0
        call BC(1)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call BC(1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        !ifEVO2 = 1
        goto 1100
     end if
1100 continue
  end if
end subroutine cllsub




SUBROUTINE PB()
USE comvar
USE mpivar
USE slfgrv
INCLUDE 'mpif.h'
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
DOUBLE PRECISION  :: VECU

double precision,  parameter :: pi = 3.14159265359d0
DOUBLE PRECISION, dimension(:,:,:), allocatable :: temp1,temp2

DOUBLE PRECISION, dimension(:,:,:),  allocatable :: data
complex*16, dimension(:,:), allocatable :: speq
DOUBLE PRECISION :: rho

DOUBLE PRECISION, dimension(:,:),  allocatable :: dat1,dat2
complex*16, dimension(:), allocatable :: spe1,spe2
double precision :: kap,temp1r,temp1i,temp2r,temp2i,facG,fac,dxx,dyy,dzz,zp1,zp2
double precision, dimension(:,:,:), allocatable :: fint0,fint1

character*4 fnum

iwx=1;iwy=1;iwz=1;N_MPI(20)=1;N_MPI(1)=1;CALL BC_MPI(2,1)

!*** Pointer for boundary ***!
!pointb1(1) = 1
!nl = 1
!nc=3
!2 continue
!pointb1(nl+1)=pointb1(nl)+nc**2
!nc=nc*2-1
!nl = nl+1
!if(nl.ne.NGcr+1) goto 2
!Needb = pointb1(NGcr+1)
!ALLOCATE(bphi1(Needb,2))

!nl = NGcr-1
!pointb2(nl) = 1
!nx=(2**NGcr)/NSPLTx+2; ny=(2**NGcr)/NSPLTy+2; nz=(2**NGcr)/NSPLTz+2
!3 continue
!pointb2(nl+1)=pointb2(nl)+max(nx*ny,ny*nz,nz*nx)
!nx=nx*2-2; ny=ny*2-2; nz=nz*2-2
!nl = nl+1
!if(nl.ne.NGL+1) goto 3
!Needb = pointb2(NGL+1)
!ALLOCATE(bphi2(Needb,2))


!MIYAMA method ---------------------------------------------------

ALLOCATE(data(Ncelly*NSPLTy,Ncellz*NSPLTz,Ncellx+1),speq(Ncellz*NSPLTz,Ncellx))
ALLOCATE(dat1(Ncelly*NSPLTy,Ncellz*NSPLTz),spe1(Ncellz*NSPLTz), &
         dat2(Ncelly*NSPLTy,Ncellz*NSPLTz),spe2(Ncellz*NSPLTz))

!nccy = Ncelly/NSPLTy; nccz = Ncellz/NSPLTz
nccy = Ncelly; nccz = Ncellz
do k=1,Ncellz; kz=KST*Ncellz+k
do j=1,Ncelly; jy=JST*Ncelly+j
do i=1,Ncellx
  data(jy,kz,i) = U(i,j,k,1)
end do;end do;end do

                    !count,blocklength,stride
!CALL MPI_TYPE_VECTOR(Ncellz,Ncelly,Ncelly*NSPLTy,MPI_REAL8,VECU,IERR)
!CALL MPI_TYPE_COMMIT(VECU,IERR)

!do Nlp = 1,NSPLTy*NSPLTz-1

!  isend = NRANK + NSPLTx*Nlp; if(isend.ge.NPE) isend = isend - NPE
!  KSs = isend/(NSPLTx*NSPLTy); JSs = isend/NSPLTx-NSPLTy*KSs
!  irecv = NRANK - NSPLTx*Nlp; if(irecv.lt.0  ) irecv = irecv + NPE
!  KSr = irecv/(NSPLTx*NSPLTy); JSr = irecv/NSPLTx-NSPLTy*KSr

!  Nis = JSs + NSPLTy*KSs
!  kls = Nis + 1
!  Nir = JST + NSPLTy*KST
!  klr = Nir + 1

!  if(kls.gt.Ncellx) then; isend = MPI_PROC_NULL; kls = Ncellx+1; end if
!  if(klr.gt.Ncellx) then; irecv = MPI_PROC_NULL; klr = Ncellx+1; end if
!  CALL MPI_SENDRECV(data(JST*Ncelly+1,KST*Ncellz+1,kls),1,VECU,isend,1, & !send
!                    data(JSr*Ncelly+1,KSr*Ncellz+1,klr),1,VECU,irecv,1, MPI_COMM_WORLD,MSTATUS,IERR) !recv
!end do

!CALL MPI_TYPE_FREE(VECU,IERR)


dxx = dy(1); dyy = dz(1); dzz = dx(1)
facG = -G4pi*dzz

dat1(:,:)=0.d0; dat2(:,:)=0.d0; spe1(:)=(0.d0,0.d0); spe2(:)=(0.d0,0.d0)

!Nir = JST + NSPLTy*KST
!klr = Nir + 1


nn1 = Ncelly*NSPLTy; nn2 = Ncellz*NSPLTz

if(klr.le.Ncellx) then

  call rlft3(data(1,1,klr),speq(1,klr),nn1,nn2,1)

  kz = klr
  zp1 = x(kz)-0.5d0*dzz
  zp2 = Lbox - zp1
  temp1r = dat1(1,1) - data(1,1,klr) * 0.5d0*zp1 * facG
  temp1i = dat1(2,1) - data(2,1,klr) * 0.5d0*zp1 * facG
  temp2r = dat2(1,1) - data(1,1,klr) * 0.5d0*zp2 * facG
  temp2i = dat2(2,1) - data(2,1,klr) * 0.5d0*zp2 * facG

  do m=1,nn2/2+1; do l=1,nn1/2
    kap = 4.d0*( sin(pi*(l-1)/nn1)**2/dxx**2 + sin(pi*(m-1)/nn2)**2/dyy**2 ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    kap = sqrt(kap)
    dat1(2*l-1,m) = dat1(2*l-1,m) + data(2*l-1,m,klr)* 0.5d0*exp(-zp1*kap)/kap *facG
    dat1(2*l  ,m) = dat1(2*l  ,m) + data(2*l  ,m,klr)* 0.5d0*exp(-zp1*kap)/kap *facG
    dat2(2*l-1,m) = dat2(2*l-1,m) + data(2*l-1,m,klr)* 0.5d0*exp(-zp2*kap)/kap *facG
    dat2(2*l  ,m) = dat2(2*l  ,m) + data(2*l  ,m,klr)* 0.5d0*exp(-zp2*kap)/kap *facG
  end do;end do

  l=nn1/2+1
  do m=1,nn2/2+1
    kap = 4.d0*( sin(pi*(l-1)/nn1)**2/dxx**2 + sin(pi*(m-1)/nn2)**2/dyy**2 )
    kap = sqrt(kap)
    spe1(m) = spe1(m) + speq(m,klr)* 0.5d0*exp(-zp1*kap)/kap *facG
    spe2(m) = spe2(m) + speq(m,klr)* 0.5d0*exp(-zp2*kap)/kap *facG
  end do

  do m2=nn2/2+2,nn2; m=nn2+2-m2; do l=1,nn1/2
    kap = 4.d0*( sin(pi*(l-1)/nn1)**2/dxx**2 + sin(pi*(m-1)/nn2)**2/dyy**2 )
    kap = sqrt(kap)
    dat1(2*l-1,m2) = dat1(2*l-1,m2) + data(2*l-1,m2,klr)* 0.5d0*exp(-zp1*kap)/kap *facG
    dat1(2*l  ,m2) = dat1(2*l  ,m2) + data(2*l  ,m2,klr)* 0.5d0*exp(-zp1*kap)/kap *facG
    dat2(2*l-1,m2) = dat2(2*l-1,m2) + data(2*l-1,m2,klr)* 0.5d0*exp(-zp2*kap)/kap *facG
    dat2(2*l  ,m2) = dat2(2*l  ,m2) + data(2*l  ,m2,klr)* 0.5d0*exp(-zp2*kap)/kap *facG
  end do;end do

  l=nn1/2+1
  do m2=nn2/2+2,nn2; m=nn2+2-m2
    kap = 4.d0*( sin(pi*(l-1)/nn1)**2/dxx**2 + sin(pi*(m-1)/nn2)**2/dyy**2 )
    kap = sqrt(kap)
    spe1(m2) = spe1(m2) + speq(m2,klr)* 0.5d0*exp(-zp1*kap)/kap *facG
    spe2(m2) = spe2(m2) + speq(m2,klr)* 0.5d0*exp(-zp2*kap)/kap *facG
  end do

  dat1(1,1) = temp1r
  dat1(2,1) = temp1i
  dat2(1,1) = temp2r
  dat2(2,1) = temp2i

end if

!CALL MPI_ALLREDUCE(dat1(1,1),data(1,1,1),Ncelly*NSPLTy*Ncellz*NSPLTz,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
!CALL MPI_ALLREDUCE(spe1(1)  ,speq(  1,1),Ncellz*NSPLTz,MPI_COMPLEX16,MPI_SUM,MPI_COMM_WORLD,IERR)
!CALL MPI_ALLREDUCE(dat2(1,1),data(1,1,2),Ncelly*NSPLTy*Ncellz*NSPLTz,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
!CALL MPI_ALLREDUCE(spe2(1)  ,speq(  1,2),Ncellz*NSPLTz,MPI_COMPLEX16,MPI_SUM,MPI_COMM_WORLD,IERR)


call rlft3(data(1,1,1),speq(1,1),nn1,nn2,-1)
call rlft3(data(1,1,2),speq(1,2),nn1,nn2,-1)

fac = 2*(nn1/2)**2
do j=1,nn2; do i=1,nn1
  data(i,j,1) = data(i,j,1)/fac
  data(i,j,2) = data(i,j,2)/fac
end do; end do

!If(NRANK.eq.0) then; open(3,file='Pot.dat')
!do j=1,nn2
!  do i=1,nn1
!    write(3,*) i,j,data(i,j,1),data(i,j,2)
!  end do; write(3,*) ' '
!end do
!close(3); END IF


ncx=Ncellx+1; ncy=Ncelly+1; ncz=Ncellz+1
do k=0,ncz; kk= (ncy+1)*k
do j=0,ncy; n = j+kk
  jb  = JST*Ncelly + j
  kbb = KST*Ncellz + k
  if((j.eq.ncy).and.(JST.eq.NSPLTy-1)) jb  = 1
  if((k.eq.ncz).and.(KST.eq.NSPLTz-1)) kbb = 1
  if((j.eq.0  ).and.(JST.eq.0       )) jb  = Ncelly*NSPLTy
  if((k.eq.0  ).and.(KST.eq.0       )) kbb = Ncellz*NSPLTz

  bphi2(pointb2(NGL)+n,1) = dble(data(jb,kbb,1))
  bphi2(pointb2(NGL)+n,2) = dble(data(jb,kbb,2))
end do; end do


!write(fnum,'(I3.3)') NRANK
!open(3,file='/work/inouety/SGte/bphi'//fnum//'.dat')
!ncx=Ncellx+1; ncy=Ncelly+1; ncz=Ncellz+1
!do k=0,ncz; kk= (ncy+1)*k
!do j=0,ncy; n = j+kk
!  write(3,*) JST*Ncelly+j,KST*Ncellz+k,bphi2(pointb2(NGL)+n,2)
!end do;write(3,*) ' '; end do
!close(3)


DEALLOCATE(data,speq)
DEALLOCATE(dat1,spe1,dat2,spe2)
!-----------------------------------------------------------------

END SUBROUTINE PB


SUBROUTINE rlft3(data,speq,nn1,nn2,isign)
INTEGER isign,nn1,nn2
COMPLEX*16 data(nn1/2,nn2),speq(nn2)
INTEGER i1,i2,j1,j2,nn(2)
DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
COMPLEX*16 c1,c2,h1,h2,w
c1=cmplx(0.5d0,0.0d0)
c2=cmplx(0.0d0,-0.5d0*isign) 
theta=6.28318530717959d0/dble(isign*nn1)
wpr=-2.0d0*sin(0.5d0*theta)**2
wpi=sin(theta)
nn(1)=nn1/2
nn(2)=nn2

if(isign.eq.1) then
  call fourn(data,nn,2,isign)
    do i2=1,nn2 
      speq(i2)=data(1,i2)
    enddo
endif
wr=1.0d0
wi=0.0d0
do i1=1,nn1/4+1 
  j1=nn1/2-i1+2 
  do i2=1,nn2 
    j2=1
    if(i2.ne.1) j2=nn2-i2+2 
    if(i1.eq.1) then
      h1=c1*(data(1,i2)+conjg(speq(j2))) 
      h2=c2*(data(1,i2)-conjg(speq(j2)))
      data(1,i2)=h1+h2 
      speq(j2)=conjg(h1-h2) 
    else
      h1=c1*(data(i1,i2)+conjg(data(j1,j2))) 
      h2=c2*(data(i1,i2)-conjg(data(j1,j2))) 
      data(i1,i2)=h1+w*h2 
      data(j1,j2)=conjg(h1-w*h2) 
    endif
  enddo
  wtemp=wr
  wr=wr*wpr-wi*wpi+wr
  wi=wi*wpr+wtemp*wpi+wi
  w=cmplx(wr,wi)
enddo

if(isign.eq.-1) call fourn(data,nn,2,isign)
END SUBROUTINE

SUBROUTINE fourn(data,nn,ndim,isign) 
INTEGER isign,ndim,nn(ndim) 
DOUBLE PRECISION data(*) 
INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,k2,n,nprev,nrem,ntot 
DOUBLE PRECISION tempi,tempr
DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
ntot=1 
do idim=1,ndim
  ntot=ntot*nn(idim)
enddo
nprev=1 
do idim=1,ndim
  n=nn(idim)
  nrem=ntot/(n*nprev) 
  ip1=2*nprev 
  ip2=ip1*n 
  ip3=ip2*nrem 
  i2rev=1 
  do i2=1,ip2,ip1
    if(i2.lt.i2rev)then 
      do i1=i2,i2+ip1-2,2 
        do i3=i1,ip3,ip2 
          i3rev=i2rev+i3-i2 
          tempr=data(i3) 
          tempi=data(i3+1) 
          data(i3)=data(i3rev) 
          data(i3+1)=data(i3rev+1) 
          data(i3rev)=tempr 
          data(i3rev+1)=tempi 
        enddo 
      enddo 
    endif 
    ibit=ip2/2 
1   if((ibit.ge.ip1).and.(i2rev.gt.ibit)) then 
      i2rev=i2rev-ibit 
      ibit=ibit/2 
      goto 1 
    endif 
    i2rev=i2rev+ibit 
  enddo
  ifp1=ip1
2 if(ifp1.lt.ip2)then 
    ifp2=2*ifp1 
    theta=isign*6.28318530717959d0/(ifp2/ip1)
    wpr=-2.d0*sin(0.5d0*theta)**2 
    wpi=sin(theta) 
    wr=1.d0 
    wi=0.d0 
    do i3=1,ifp1,ip1 
      do i1=i3,i3+ip1-2,2 
        do i2=i1,ip3,ifp2 
          k1=i2
          k2=k1+ifp1 
          tempr=wr*data(k2)-wi*data(k2+1) 
          tempi=wr*data(k2+1)+wi*data(k2) 
          data(k2)=data(k1)-tempr 
          data(k2+1)=data(k1+1)-tempi
          data(k1)=data(k1)+tempr 
          data(k1+1)=data(k1+1)+tempi 
        enddo
      enddo
      wtemp=wr
      wr=wr*wpr-wi*wpi+wr 
      wi=wi*wpr+wtemp*wpi+wi 
    enddo 
    ifp1=ifp2 
    goto 2 
  endif 
  nprev=n*nprev 
enddo 
END SUBROUTINE

