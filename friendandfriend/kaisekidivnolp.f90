program main
  implicit none
  integer i,j,k,val,m,loop,n,n1,n2,n3,dm
  integer nx,ny,nz,rgn,zero
  integer,allocatable,dimension(:,:,:,:) ::tag
  real(4),allocatable,dimension(:,:,:,:) ::U
  character(3) num
  character(1) num1
  double precision Lbox,Msun,dl,Mshd,ganmim1,phi,Mtot
  integer,allocatable,dimension(:) :: cas,count1
  integer,allocatable,dimension(:) :: ig,jg,kg
  real(4) dm1
  double precision ,allocatable,dimension(:) :: MM
  double precision ,allocatable,dimension(:) :: x,y,z,xg,yg,zg,vx,vy,vz,vxg,vyg,vzg
  double precision ,allocatable,dimension(:) :: emag,ekin,ethm,egrv,div

  !----parameter----
  nx=512
  ny=512
  nz=512
  val=9
  loop=1
  Lbox=100.d0
  Msun=1.473d-2 !1pc * 1pc * 1pc * 1m_p/cc
  rgn=100
  Mshd=1.d4
  ganmim1=3.d0/2.d0 !1.d0/((5.d0/3.d0-1.d0)
  !----parameter----

  dl=Lbox/dble(nx)

  allocate(U(-1:nx+2,-1:ny+2,-1:nz+2,val))
  allocate(tag(1:nx,1:ny,1:nz,loop))
  allocate(x(nx),y(ny),z(nz),xg(rgn),yg(rgn),zg(rgn),MM(rgn))
  allocate(vx(rgn),vy(rgn),vz(rgn),vxg(rgn),vyg(rgn),vzg(rgn))
  allocate(emag(rgn),ekin(rgn),ethm(rgn),egrv(rgn))
  allocate(cas(0:rgn),count1(rgn),div(rgn),ig(rgn),jg(rgn),kg(rgn))

  U(:,:,:,:)=0.e0
  tag(:,:,:,:)=0

  ! do m=1,loop
  m=1
  !m=1
  write(num,'(i3.3)') m
  write(num1,'(i1.1)') m
  open(150,file='/Users/maeda/Desktop/kaiseki/image/tag'//num//'.DAT',FORM='FORMATTED')
  write(*,*) num
  m=1
  do k=1,nz
     do j=1,ny
        !write(*,*) i,j,k,m,'read-tag'
        do i=1,nx
           !write(*,*) i,j,k,m,'read-tag'
           read(150,*) tag(i,j,k,m)!,dm,dm1
           !write(*,*) i,j,k,tag(i,j,k,m)
        end do
     end do
     write(*,*) k,m,'read-tag'
  end do
  close(150)
  !end do

  open(110,file='/Users/maeda/Desktop/kaiseki/cnv100wbwg/all.DAT',FORM='FORMATTED')
  do k=1,nz
     do j=1,ny
        do i=1,nz
           read(110,*) U(i,j,k,1), U(i,j,k,2), U(i,j,k,3), U(i,j,k,4), U(i,j,k,5),&
                U(i,j,k,6), U(i,j,k,7), U(i,j,k,8), U(i,j,k,9)
        end do
     end do
     write(*,*) k, U(5,5,5,2),'read-phyval'
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

  x(1)=dl*0.5d0
  y(1)=dl*0.5d0
  z(1)=dl*0.5d0
  do i=2,nx
     x(i)=x(i-1)+dl
     y(i)=y(i-1)+dl
     z(i)=z(i-1)+dl
  end do

  MM(:)=0.d0
  count1(:)=0
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
  do n=1,rgn
     phi=1.d5
     !Mtot=0.d0
     do k=1,nz
        do j=1,ny
           do i=1,nz
              if(n==1) then
                 Mtot=Mtot+dble(U(i,j,k,1))*(dl**3)
              end if
              if(n==tag(i,j,k,m)) then
                 MM(n)=MM(n)+dble(U(i,j,k,1))*(dl**3)
                 vx(n)=dble(U(i,j,k,2))+vx(n)
                 vy(n)=dble(U(i,j,k,3))+vy(n)
                 vz(n)=dble(U(i,j,k,4))+vz(n)
                 !vxg(n)=dble(U(i,j,k,2))*dble(U(i,j,k,1))*(dl**3)+vxg(n)
                 !vyg(n)=dble(U(i,j,k,3))*dble(U(i,j,k,1))*(dl**3)+vyg(n)
                 !vzg(n)=dble(U(i,j,k,4))*dble(U(i,j,k,1))*(dl**3)+vzg(n)
                 count1(n)=count1(n)+1
                 if(dble(U(i,j,k,9))<phi)then
                    phi=dble(U(i,j,k,9))
                    ig(n)=i
                    jg(n)=j
                    kg(n)=k
                 end if
                 !xg(n)=dble(U(i,j,k,1))*x(i)*(dl**3)+xg(i)
                 !yg(n)=dble(U(i,j,k,1))*y(i)*(dl**3)+yg(i)
                 !zg(n)=dble(U(i,j,k,1))*z(i)*(dl**3)+zg(i)
              end if
           end do
        end do
        write(*,*) k,n,'mass'
     end do
     vx(n)=vx(n)/dble(count1(n))
     vy(n)=vy(n)/dble(count1(n))
     vz(n)=vz(n)/dble(count1(n))
     !vxg(n)=vx(n)/MM(n)
     !vyg(n)=vy(n)/MM(n)
     !vzg(n)=vz(n)/MM(n)
     !xg(n)=xg(i)/MM(n)
     !yg(n)=yg(i)/MM(n)
     !zg(n)=zg(i)/MM(n)
  end do
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
  do n=1,rgn
     if(cas(n)==1)then
        do k=1,nz
           do j=1,ny
              do i=1,nz
                 if(n==tag(i,j,k,m)) then
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
                 end if
              end do
           end do
           write(*,*) k,n,'div'
        end do
     end if
  end do

  zero=0
  cas(0)=0
  !open(510,file='/Users/maeda/Desktop/kaiseki/image/lb4.DAT',FORM='FORMATTED')
!  do n=1,rgn
!     if(cas(n)==1)then
  !write(num,'(i3.3)') 000
  !write(num,'(i3.3)') 000
        open(510,file='/Users/maeda/Desktop/kaiseki/image/lb'//num1//'.DAT',FORM='FORMATTED')
        !if(cas(n)==1)then
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 !if((cas(int(tag(i,j,k,m)))==1).and.(tag(i,j,k,m).ne.0))then
              !   if(tag(i,j,k,m)==n)then
                    write(510,*) tag(i,j,k,m)!*cas(tag(i,j,k,m))
                    !write(*,*) tag(i,j,k,m)
              !   else
              !      write(510,*) zero
                    !write(*,*) zero
              !   end if
                 !write(510,*) tag(i,j,k)
              end do
           end do
           write(*,*) k,'tag'
        end do
        !end if
        close(510)
!     end if
!  end do
  !close(510)

  open(610,file='/Users/maeda/Desktop/kaiseki/image/mdiv'//num1//'.DAT',FORM='FORMATTED')
  do n=1,rgn
     if(cas(n)==1)then
        write(610,*) n,MM(n),div(n),emag(n),ekin(n),ethm(n),-egrv(n),x(ig(n)),y(jg(n)),z(kg(n)),&
             ekin(n)+ekin(n)+ethm(n)-egrv(n)
     end if
     !write(510,*) tag(i,j,k)
  end do

  close(610)


  !ekin = 0.5d0 * (   U(i,j,k,2)**2 + U(i,j,k,3)**2 + U(i,j,k,4)**2 ) * U(i,j,k,1)
  !emag = 0.5d0 * ( Bcc(i,j,k,iBcc)**2+Blg(i,j,k,1)**2+Blg(i,j,k,2)**2 )
  !gammi1 =   3.d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))+5.d0*ndH2(i,j,k)
  !gammi1 = ( 2.d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))+2.d0*ndH2(i,j,k) )/gammi1
  !U(i,j,k,5) = (  U(i,j,k,5)-ekin-emag  )*gammi1

  deallocate(U)
  deallocate(tag,MM)
  deallocate(x,y,z,xg,yg,zg,vx,vy,vz,vxg,vyg,vzg)
  deallocate(emag,ekin,ethm,egrv)
  deallocate(cas,count1,div,ig,jg,kg)
end program main
