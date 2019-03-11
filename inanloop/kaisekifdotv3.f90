program main
  implicit none
  integer i,j,k,val,m,loop,n,n1,n2,n3,dm,time1,timeloop,jump,timestep,cnt10
  integer nx,ny,nz,rgn,zero,sloop,smesh
  integer,allocatable,dimension(:,:,:,:) ::tag
  real(4),allocatable,dimension(:,:,:,:) ::U
  character(3) num,timech
  character(1) num1
  double precision Lbox,Msun,dl,Mshd,ganmim1,Mtot
  integer,allocatable,dimension(:) :: cas,count1,count1S
  integer,allocatable,dimension(:) :: ig,jg,kg
  real(4) dm1
  double precision ,allocatable,dimension(:) :: MM,S,lgM,phi
  double precision ,allocatable,dimension(:) :: x,y,z,xg,yg,zg,vx,vy,vz,vxg,vyg,vzg
  double precision ,allocatable,dimension(:) :: emag,ekin,ethm,egrv,div
  double precision ,allocatable,dimension(:,:,:) :: jx,jy,jz,fx,fy,fz,fgx,fgy,fgz,fpx,fpy,fpz,fbx,fby,fbz,divv
  double precision ,allocatable,dimension(:,:,:) :: vdotf,vdotfg,vdotfp,vdotfb
  integer,allocatable,dimension(:):: lgMcnt
  double precision cs,nsisu,ff,vv,ffb,ffg,ffp

  !----parameter----
  nx=512
  ny=512
  nz=512
  val=9
  loop=5
  Lbox=100.d0
  Msun=1.473d-2 !1pc * 1pc * 1pc * 1m_p/cc
  rgn=1000
  Mshd=1.d2
  ganmim1=3.d0/2.d0 !1.d0/((5.d0/3.d0-1.d0)
  sloop=30
  smesh=5
  timeloop=7
  timestep=15
  jump=4
  !----parameter----

  dl=Lbox/dble(nx)

  allocate(U(-1:nx+2,-1:ny+2,-1:nz+2,val))
  !allocate(tag(1:nx,1:ny,1:nz,loop))
  allocate(tag(0:nx+1,0:ny+1,0:nz+1,loop))
  allocate(lgM(0:sloop),lgMcnt(sloop))
  allocate(x(nx),y(ny),z(nz),xg(rgn),yg(rgn),zg(rgn),MM(rgn),S(rgn))
  allocate(vx(rgn),vy(rgn),vz(rgn),vxg(rgn),vyg(rgn),vzg(rgn))
  allocate(emag(rgn),ekin(rgn),ethm(rgn),egrv(rgn))
  allocate(cas(0:rgn),count1(rgn),div(rgn),ig(rgn),jg(rgn),kg(rgn),count1S(rgn),phi(rgn))
  allocate(jx(1:nx,1:ny,1:nz),jy(1:nx,1:ny,1:nz),jz(1:nx,1:ny,1:nz),fx(1:nx,1:ny,1:nz),fy(1:nx,1:ny,1:nz),fz(1:nx,1:ny,1:nz)&
       ,fgx(1:nx,1:ny,1:nz),fgy(1:nx,1:ny,1:nz),fgz(1:nx,1:ny,1:nz),fpx(1:nx,1:ny,1:nz),fpy(1:nx,1:ny,1:nz),fpz(1:nx,1:ny,1:nz)&
       ,fbx(1:nx,1:ny,1:nz),fby(1:nx,1:ny,1:nz),fbz(1:nx,1:ny,1:nz),divv(1:nx,1:ny,1:nz))
  allocate(vdotf(1:nx,1:ny,1:nz),vdotfg(1:nx,1:ny,1:nz),vdotfp(1:nx,1:ny,1:nz),vdotfb(1:nx,1:ny,1:nz))

  !U(:,:,:,:)=0.e0
  tag(:,:,:,:)=0


  nsisu=2
  lgM(:)=0.d0
  lgM(0)=10.0d0**nsisu
  do j=1,sloop
     nsisu=2
     nsisu =nsisu + 1.d0/dble(smesh)*dble(j)
     lgM(j) = 10.0d0**nsisu
     write(*,*) lgM(j)
  end do


  do time1=1,timeloop
     open(unit=350,file='cnt3.dat')
     read(350,*) cnt10
     close(350)
     !write(cntc,'(i3.3)') cnt10
     write(timech,'(i3.3)') cnt10
     cnt10=cnt10+timestep
     open(unit=350,file='cnt3.dat')
     write(350,*) cnt10
     close(350)
     !write(timech,'(i3.3)') time1

  do m=1,loop
  !m=5
  !m=1
  write(num,'(i3.3)') m
  write(num1,'(i1.1)') m
  open(150,file='taglog'//timech//num//'.DAT',FORM='FORMATTED')
  write(*,*) num
  !m=1
  do k=1,nz
     do j=1,ny
        !write(*,*) i,j,k,m,'read-tag'
        do i=1,nx
           !write(*,*) i,j,k,m,'read-tag'
           read(150,*) tag(i,j,k,m)!,dm,dm1
           !write(*,*) i,j,k,tag(i,j,k,m)
        end do
     end do
     write(*,*) k,m,timech,'read-tag'
  end do
  close(150)
end do

  U(:,:,:,:)=0.e0
  open(110,file='Allnewbig'//timech//'.DAT',FORM='UNFORMATTED')
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
  U(:,:,-1,:)=U(:,:,nz-2,:)
  U(:,:, 0,:)=U(:,:,nz-1,:)
  U(:,:,nz+1,:)=U(:,:,1,:)
  U(:,:,nz+2,:)=U(:,:,2,:)
  U(:,-1,:,:)=U(:,ny-2,:,:)
  U(:, 0,:,:)=U(:,ny-1,:,:)
  U(:,ny+1,:,:)=U(:,1,:,:)
  U(:,ny+2,:,:)=U(:,2,:,:)
  !U(-1,:,:,:)=U(nx-2,:,:,:)
  !U( 0,:,:,:)=U(nx-1,:,:,:)
  !U(nx+1,:,:,:)=U(1,:,:,:)
  !U(nx+2,:,:,:)=U(2,:,:,:)
  U(-1,:,:,:)=U(1,:,:,:)
  U( 0,:,:,:)=U(1,:,:,:)
  U(nx+1,:,:,:)=U(nx,:,:,:)
  U(nx+2,:,:,:)=U(nx,:,:,:)

  jx(:,:,:)=0.d0
  jy(:,:,:)=0.d0
  jz(:,:,:)=0.d0
  fx(:,:,:)=0.d0
  fy(:,:,:)=0.d0
  fz(:,:,:)=0.d0
  fgx(:,:,:)=0.d0
  fgy(:,:,:)=0.d0
  fgz(:,:,:)=0.d0
  fpx(:,:,:)=0.d0
  fpy(:,:,:)=0.d0
  fpz(:,:,:)=0.d0
  fbx(:,:,:)=0.d0
  fby(:,:,:)=0.d0
  fbz(:,:,:)=0.d0
  divv(:,:,:)=0.d0

  do k=1,nz
     do j=1,ny
        do i=1,nx
           jx(i,j,k) = (dble(U(i,j,k+1,8))- dble(U(i,j,k-1,8))-( dble(U(i,j+1,k,7))- dble(U(i,j-1,k,7))))/dl
           jy(i,j,k) = (dble(U(i+1,j,k,6))- dble(U(i-1,j,k,6))-( dble(U(i,j,k+1,8))- dble(U(i,j,k-1,8))))/dl
           jz(i,j,k) = (dble(U(i,j+1,k,7))- dble(U(i,j-1,k,7))-( dble(U(i+1,j,k,6))- dble(U(i-1,j,k,6))))/dl

           fbx(i,j,k) = jy(i,j,k)*dble(U(i,j,k,8))-jz(i,j,k)*dble(U(i,j,k,7))
           fby(i,j,k) = jz(i,j,k)*dble(U(i,j,k,6))-jx(i,j,k)*dble(U(i,j,k,8))
           fbz(i,j,k) = jx(i,j,k)*dble(U(i,j,k,7))-jy(i,j,k)*dble(U(i,j,k,6))

           fpx(i,j,k)=-(dble(U(i+1,j,k,5))-dble(U(i-1,j,k,5)))/dl
           fpy(i,j,k)=-(dble(U(i,j+1,k,5))-dble(U(i,j-1,k,5)))/dl
           fpz(i,j,k)=-(dble(U(i,j,k+1,5))-dble(U(i,j,k-1,5)))/dl

           fgx(i,j,k)=-(dble(U(i+1,j,k,9))-dble(U(i-1,j,k,9)))/dl*dble(U(i,j,k,1))
           fgy(i,j,k)=-(dble(U(i,j+1,k,9))-dble(U(i,j-1,k,9)))/dl*dble(U(i,j,k,1))
           fgz(i,j,k)=-(dble(U(i,j,k+1,9))-dble(U(i,j,k-1,9)))/dl*dble(U(i,j,k,1))

           fx(i,j,k)=fbx(i,j,k)+fpx(i,j,k)+fgx(i,j,k)
           fy(i,j,k)=fby(i,j,k)+fpy(i,j,k)+fgy(i,j,k)
           fz(i,j,k)=fbz(i,j,k)+fpz(i,j,k)+fgz(i,j,k)

           ff=dsqrt(fx(i,j,k)**2+fy(i,j,k)**2+fz(i,j,k)**2)
           ffg=dsqrt(fgx(i,j,k)**2+fgy(i,j,k)**2+fgz(i,j,k)**2)
           ffb=dsqrt(fbx(i,j,k)**2+fby(i,j,k)**2+fbz(i,j,k)**2)
           ffp=dsqrt(fpx(i,j,k)**2+fpy(i,j,k)**2+fpz(i,j,k)**2)
           vv=dsqrt(dble(U(i,j,k,2))**2+dble(U(i,j,k,3))**2+dble(U(i,j,k,4))**2)

           vdotf(i,j,k) =(fx(i,j,k)*dble(U(i,j,k,2))+fy(i,j,k)*dble(U(i,j,k,3))+fz(i,j,k)*dble(U(i,j,k,4)))&
                /(ff*vv)
           !vdotfb(i,j,k)=fbx(i,j,k)*dble(U(i,j,k,2))+fby(i,j,k)*dble(U(i,j,k,3))+fbz(i,j,k)*dble(U(i,j,k,4))
           !vdotfg(i,j,k)=fgx(i,j,k)*dble(U(i,j,k,2))+fgy(i,j,k)*dble(U(i,j,k,3))+fgz(i,j,k)*dble(U(i,j,k,4))
           !vdotfp(i,j,k)=fpx(i,j,k)*dble(U(i,j,k,2))+fpy(i,j,k)*dble(U(i,j,k,3))+fpz(i,j,k)*dble(U(i,j,k,4))
           vdotfb(i,j,k)=(fbx(i,j,k)*dble(U(i,j,k,2))+fby(i,j,k)*dble(U(i,j,k,3))+fbz(i,j,k)*dble(U(i,j,k,4)))&
                /(ffb*vv)!vdotf(i,j,k)
           vdotfg(i,j,k)=(fgx(i,j,k)*dble(U(i,j,k,2))+fgy(i,j,k)*dble(U(i,j,k,3))+fgz(i,j,k)*dble(U(i,j,k,4)))&
                /(ffg*vv)!vdotf(i,j,k)
           vdotfp(i,j,k)=(fpx(i,j,k)*dble(U(i,j,k,2))+fpy(i,j,k)*dble(U(i,j,k,3))+fpz(i,j,k)*dble(U(i,j,k,4)))&
                /(ffp*vv)!vdotf(i,j,k)

!           cs = 5.d0/3.d0 * dble(U(i,j,k,5)) / dble(U(i,j,k,1))
!           divv(i,j,k) = (dble(U(i+1,j,k,2))-dble(U(i-1,j,k,2))+dble(U(i,j+1,k,3))&
!                -dble(U(i,j-1,k,3))+dble(U(i,j,k+1,4))-dble(U(i,j,k-1,4)))/dl/cs

        end do
     end do
     write(*,*) k, 'ok1'
  end do


  x(1)=dl*0.5d0
  y(1)=dl*0.5d0
  z(1)=dl*0.5d0
  do i=2,nx
     x(i)=x(i-1)+dl
     y(i)=y(i-1)+dl
     z(i)=z(i-1)+dl
  end do

  do m=1,loop
  write(num,'(i3.3)') m
  MM(:)=0.d0
  count1(:)=0
  count1S(:)=0
  S(:)=0.d0
  vx(:)=0.d0
  vy(:)=0.d0
  vz(:)=0.d0
  !vxg(:)=0.d0
  !vyg(:)=0.d0
  !vzg(:)=0.d0
  !xg(:)=0.d0
  !yg(:)=0.d0
  !zg(:)=0.d0
  ig(:)=0
  jg(:)=0
  kg(:)=0
  Mtot=0.d0
  tag(0,:,:,m)=tag(nx,:,:,m)
  tag(:,0,:,m)=tag(:,ny,:,m)
  tag(:,:,0,m)=tag(:,:,nz,m)
  tag(nx+1,:,:,m)=tag(1,:,:,m)
  tag(:,ny+1,:,m)=tag(:,1,:,m)
  tag(:,:,nz+1,m)=tag(:,:,1,m)


  !do n=1,rgn
     phi(:)=1.d5
     !Mtot=0.d0
     do k=1,nz
        do j=1,ny
           do i=1,nz
              n=tag(i,j,k,m)
              !if(n==1) then
                 Mtot=Mtot+dble(U(i,j,k,1))*(dl**3)
              !end if
                 !if(n==tag(i,j,k,m)) then
               if(n.ne.0) then
                 MM(n)=MM(n)+dble(U(i,j,k,1))*(dl**3)
                 vx(n)=dble(U(i,j,k,2))+vx(n)
                 vy(n)=dble(U(i,j,k,3))+vy(n)
                 vz(n)=dble(U(i,j,k,4))+vz(n)
                 !vxg(n)=dble(U(i,j,k,2))*dble(U(i,j,k,1))*(dl**3)+vxg(n)
                 !vyg(n)=dble(U(i,j,k,3))*dble(U(i,j,k,1))*(dl**3)+vyg(n)
                 !vzg(n)=dble(U(i,j,k,4))*dble(U(i,j,k,1))*(dl**3)+vzg(n)
                 count1(n)=count1(n)+1
                 count1S(n)=count1S(n)+(6-((tag(i,j,k,m)*tag(i-1,j,k,m)+tag(i,j,k,m)*tag(i+1,j,k,m)&
                      +tag(i,j,k,m)*tag(i,j-1,k,m)+tag(i,j,k,m)*tag(i,j+1,k,m)+tag(i,j,k,m)*tag(i,j,k-1,m)&
                      +tag(i,j,k,m)*tag(i,j,k+1,m)))/n/n)!S=dl*dl*count1S and if(count1S(i,j,k).ne.0)then surface
              !   write(*,*) ((tag(i,j,k,m)*tag(i-1,j,k,m)+tag(i,j,k,m)*tag(i+1,j,k,m)&
              !        +tag(i,j,k,m)*tag(i,j-1,k,m)+tag(i,j,k,m)*tag(i,j+1,k,m)+tag(i,j,k,m)*tag(i,j,k-1,m)&
              !        +tag(i,j,k,m)*tag(i,j,k+1,m)))/n/n ,'countS 6'
                 if(dble(U(i,j,k,9))<phi(n))then
                    phi(n)=dble(U(i,j,k,9))
                    ig(n)=i
                    jg(n)=j
                    kg(n)=k
                 end if
              end if
                 !xg(n)=dble(U(i,j,k,1))*x(i)*(dl**3)+xg(i)
                 !yg(n)=dble(U(i,j,k,1))*y(i)*(dl**3)+yg(i)
                 !zg(n)=dble(U(i,j,k,1))*z(i)*(dl**3)+zg(i)
              !end if
           end do
        end do
        write(*,*) k,n,'mass'
     end do

     do n=1,rgn
     if(count1(n).ne.0) then
     vx(n)=vx(n)/dble(count1(n))
     vy(n)=vy(n)/dble(count1(n))
     vz(n)=vz(n)/dble(count1(n))
     S(n)=dl*dl*dble(count1S(n))
     write(*,*) n, S(n),count1S(n) ,'S-s-S'
     !vxg(n)=vx(n)/MM(n)
     !vyg(n)=vy(n)/MM(n)
     !vzg(n)=vz(n)/MM(n)
     !xg(n)=xg(i)/MM(n)
     !yg(n)=yg(i)/MM(n)
     !zg(n)=zg(i)/MM(n)
     end if
     enddo
  !end do
  Mtot=Mtot*Msun
  write(*,*) Mtot,'Mtot'

  cas(:)=0
  do n=1,rgn
     MM(n)=MM(n)*Msun
     write(*,*) n,MM(n),'MM'
     if( MM(n) > Mshd) then
        cas(n) = 1
        write(*,*) n,cas(n),'ntag'
     else
        cas(n) = 0
     end if
  end do

  !-----nakatsu-disp------
!  for (i=il; i<=ir; i++)
!{
!	vdisp = vdisp + U[i]*U[i] + V[i]*V[i] + W[i]*W[i];
!	vxdis = vxdis + U[i]*U[i];
!	vydis = vydis + V[i]*V[i];
!	vzdis = vzdis + W[i]*W[i];
!	vxave = vxave + U[i];
!	vyave = vyave + V[i];
!	vzave = vzave + W[i];
!	cell++;
!}
!if (cell .ne. 0)
!{
!	vxave = vxave / cell;
!	vyave = vyave / cell;
!	vzave = vzave / cell;
!	vdisp = sqrt(vdisp / cell - vxave*vxave - vyave*vyave - vzave*vzave);
!	vxdis = sqrt(vxdis / cell - vxave*vxave);
!	vydis = sqrt(vydis / cell - vyave*vyave);
!	vzdis = sqrt(vzdis / cell - vzave*vzave);
! }
 !-----nakatsu-disp------


  div(:)=0.d0
  emag(:)=0.d0
  ethm(:)=0.d0
  ekin(:)=0.d0
  egrv(:)=0.d0
  !do n=1,rgn
     !if(cas(n)==1)then
        do k=1,nz
           do j=1,ny
              do i=1,nz
                 !if(n==tag(i,j,k,m)) then
                 n=tag(i,j,k,m)
                    div(n)=(dble(U(i+1,j,k,2))-dble(U(i-1,j,k,2))+dble(U(i,j+1,k,3))-dble(U(i,j-1,k,3))&
                         +dble(U(i,j,k+1,4))-dble(U(i,j,k-1,4)))*(dl**3)+div(n)
                    emag(n) = 0.5d0 * ( dble(U(i,j,k,6))**2+dble(U(i,j,k,7))**2+dble(U(i,j,k,8))**2)*(dl**3) + emag(n)
                    ethm(n) = ganmim1*dble(U(i,j,k,5))*(dl**3) + ethm(n)
                    egrv(n) = dble(U(i,j,k,1))*((x(i)-x(ig(n)))*(dble(U(i+1,j,k,9))&
                         -dble(U(i-1,j,k,9)))+(y(j)-y(jg(n)))*(dble(U(i,j+1,k,9))&
                         -dble(U(i,j-1,k,9)))+(z(k)-z(kg(n)))*(dble(U(i,j,k+1,9))&
                         -dble(U(i,j,k-1,9))))*(dl**3)+egrv(n)
                    ekin(n) = 0.5d0 * dble(U(i,j,k,1)) * ((dble(U(i,j,k,2)-vx(n)))**2 + (dble(U(i,j,k,3)-vy(n)))**2&
                         + (dble(U(i,j,k,4))-vz(n))**2)*(dl**3)+ekin(n)
                 !end if
              end do
           end do
           write(*,*) k,n,'div'
        end do
     !end if
  !end do

  zero=0
  cas(0)=0
  !open(510,file='/Users/maeda/Desktop/kaiseki/image/lb4.DAT',FORM='FORMATTED')
!  do n=1,rgn
!     if(cas(n)==1)then
  !write(num,'(i3.3)') 000
  !write(num,'(i3.3)') 000
    !    open(510,file='lb'//num1//'.DAT',FORM='FORMATTED')
        !if(cas(n)==1)then
    !    do k=1,nz
    !       do j=1,ny
    !          do i=1,nx
                 !if((cas(int(tag(i,j,k,m)))==1).and.(tag(i,j,k,m).ne.0))then
              !   if(tag(i,j,k,m)==n)then
    !                write(510,*) tag(i,j,k,m)!*cas(tag(i,j,k,m))
                    !write(*,*) tag(i,j,k,m)
              !   else
              !      write(510,*) zero
                    !write(*,*) zero
              !   end if
                 !write(510,*) tag(i,j,k)
    !          end do
    !       end do
    !       write(*,*) k,'tag'
    !    end do
        !end if
    !    close(510)
!     end if
!  end do
  !close(510)

  open(610,file='mdiv'//timech//num//'.DAT',FORM='FORMATTED')
  do n=1,rgn
     if(cas(n)==1)then
        write(610,*) n,MM(n),div(n),emag(n),ekin(n),ethm(n),-egrv(n),x(ig(n)),y(jg(n)),z(kg(n)),&
             ekin(n)+ekin(n)+ethm(n)-egrv(n),div(n)/S(n),S(n),dble(count1(n))*dl*dl*dl
     end if
     !write(510,*) tag(i,j,k)
  end do
  close(610)

  lgMcnt(:)=0
  do i=1,sloop
     do n=1,rgn
        if((lgM(i-1)<MM(n)).and.(lgM(i)>MM(n)))then
           lgMcnt(i)=lgMcnt(i)+1
        end if
     end do
  end do

  open(910,file='masfun'//timech//num//'.DAT',FORM='FORMATTED')
  do i=1,sloop
     write(910,*) 10.d0**((dlog10(lgM(i-1))+dlog10(lgM(i)))*0.5d0),lgMcnt(i)
  end do
  close(910)
end do
open(100,file='Allphy'//timech//'.DAT',FORM='FORMATTED')
do k=1,nz,jump
   do j=1,ny,jump
      do i=1,nx,jump
         write(100,*) U(i,j,k,1),U(i,j,k,2),U(i,j,k,3),U(i,j,k,4),U(i,j,k,5),U(i,j,k,6),U(i,j,k,7),&
              U(i,j,k,8),U(i,j,k,9),vdotf(i,j,k),vdotfb(i,j,k),vdotfg(i,j,k),vdotfp(i,j,k),real(tag(i,j,k,1))&
              ,real(tag(i,j,k,2)),real(tag(i,j,k,3)),real(tag(i,j,k,4)),real(tag(i,j,k,5))
      end do
   end do
   write(*,*) k
end do
close(100)

end do

  !ekin = 0.5d0 * (   U(i,j,k,2)**2 + U(i,j,k,3)**2 + U(i,j,k,4)**2 ) * U(i,j,k,1)
  !emag = 0.5d0 * ( Bcc(i,j,k,iBcc)**2+Blg(i,j,k,1)**2+Blg(i,j,k,2)**2 )
  !gammi1 =   3.d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))+5.d0*ndH2(i,j,k)
  !gammi1 = ( 2.d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))+2.d0*ndH2(i,j,k) )/gammi1
  !U(i,j,k,5) = (  U(i,j,k,5)-ekin-emag  )*gammi1

  deallocate(U,lgM,lgMcnt,phi)
  deallocate(tag,MM,S)
  deallocate(x,y,z,xg,yg,zg,vx,vy,vz,vxg,vyg,vzg)
  deallocate(emag,ekin,ethm,egrv)
  deallocate(cas,count1,div,ig,jg,kg,count1S)
  deallocate(jx,jy,jz,fx,fy,fz,fgx,fgy,fgz,fpx,fpy,fpz,fbx,fby,fbz,divv)
  deallocate(vdotf,vdotfg,vdotfp,vdotfb)

  !U(:,:,:,:)=0.e0
  !tag(:,:,:,:)=0


end program main



!subroutine calf


!jx(i,j,k) = (U(i,j,k+1,8)- U(i,j,k-1,8)-( U(i,j+1,k,7)- U(i,j-1,k,7)))/dl
!jy(i,j,k) = (U(i+1,j,k,6)- U(i-1,j,k,6)-( U(i,j,k+1,8)- U(i,j,k-1,8)))/dl
!jz(i,j,k) = (U(i,j+1,k,7)- U(i,j-1,k,7)-( U(i+1,j,k,6)- U(i-1,j,k,6)))/dl

!fbx(i,j,k) = jy(i,j,k)*U(i,j,k,8)-jz(i,j,k)*U(i,j,k,7)
!fby(i,j,k) = jz(i,j,k)*U(i,j,k,6)-jx(i,j,k)*U(i,j,k,8)
!fbz(i,j,k) = jx(i,j,k)*U(i,j,k,7)-jy(i,j,k)*U(i,j,k,6)

!fpx(i,j,k)=-(U(i+1,j,k,5)-U(i-1,j,k,5))/dl
!fpy(i,j,k)=-(U(i,j+1,k,5)-U(i,j-1,k,5))/dl
!fpz(i,j,k)=-(U(i,j,k+1,5)-U(i,j,k-1,5))/dl

!fgx(i,j,k)=-(U(i+1,j,k,17)-U(i-1,j,k,17))/dl
!fgy(i,j,k)=-(U(i,j+1,k,17)-U(i,j-1,k,17))/dl
!fgz(i,j,k)=-(U(i,j,k+1,17)-U(i,j,k-1,17))/dl

!fx(i,j,k)=fbx(i,j,k)+fpx(i,j,k)+fgx(i,j,k)
!fy(i,j,k)=fby(i,j,k)+fpy(i,j,k)+fgy(i,j,k)
!fz(i,j,k)=fbz(i,j,k)+fpz(i,j,k)+fgz(i,j,k)

!ff=dsqrt(fx(i,j,k)**2+fy(i,j,k)**2+fz(i,j,k)**2)

!vdotf(i,j,k) =fx(i,j,k)*U(i,j,k,2)+fy(i,j,k)*U(i,j,k,3)+fz(i,j,k)*U(i,j,k,4)
!vdotfb(i,j,k)=fbx(i,j,k)*U(i,j,k,2)+fby(i,j,k)*U(i,j,k,3)+fbz(i,j,k)*U(i,j,k,4)
!vdotfg(i,j,k)=fgx(i,j,k)*U(i,j,k,2)+fgy(i,j,k)*U(i,j,k,3)+fgz(i,j,k)*U(i,j,k,4)
!vdotfp(i,j,k)=fpx(i,j,k)*U(i,j,k,2)+fpy(i,j,k)*U(i,j,k,3)+fpz(i,j,k)*U(i,j,k,4)

!cs = 5.d0/3.d0 * U(i,j,k,5) / U(i,j,k,1)
!divv = (U(i+1,j,k,2)-U(i-1,j,k,2)+U(i,j+1,k,3)-U(i,j-1,k,3)+U(i,j,k+1,4)-U(i,j,k-1,4))/dl/cs

!end subroutine calf

!subroutine lag



!subroutine chng

!  do j=1,rgn
!     Mmax=0.d0
!     imax=0
!     do i=1,rgn
!        if(MM(i)>Mmax) then
!           Mmax=MM(i)
!           imax=i
!        end if
!     end do
!     Mnew(k)=Mmax
!     MM(imax)=0.d0
!  end do

!end subroutine chng
