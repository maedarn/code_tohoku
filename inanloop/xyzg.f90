program main
  implicit none
  integer :: i,j,k,nx,ny,nz,lev,val,nx1,ny1,nz1,nx2,ny2,nz2,ilev,rgn,n,ii
  real(4),allocatable,dimension(:,:,:,:) :: U,Usave,Udummy
  integer,allocatable,dimension(:,:,:) :: tag,tagsave,tagdummy
  integer,allocatable,dimension(:,:,:,:,:) ::tagsph
  real(4) dl,Lbox,dl1,Msun,shd,Mshd,lx,ly,lz,r,rshd
  character(3) num,ch
  real(4),allocatable,dimension(:,:) :: xg,yg,zg,MM,mxrho1,Phi,x,y,z!,Ixx,Ixy,Ixz,Iyy,Iyz,Izz
  real(4),allocatable,dimension(:,:,:) :: Mss,Ixx,Ixy,Ixz,Iyy,Iyz,Izz
  integer,allocatable,dimension(:,:) :: count1,ig,jg,kg,igr,jgr,kgr!,cntnum

  !parameter
  nx=512
  ny=512
  nz=512
  Lbox=100.e0
  lev=2
  val=9
  rgn=1000
  Msun=1.473e-2 !1pc * 1pc * 1pc * 1m_p/cc
  shd=1.e4
  Mshd=1.e2
  rshd=5.e0
  !parameter
  dl=Lbox/real(nx)

  allocate(U(-1:nx+2,-1:ny+2,-1:nz+2,val))
  allocate(tag(0:nx+1,0:ny+1,0:nz+1),tagsph(1:nx,1:ny,1:nz,0:lev,3))
  allocate(x(nx,0:lev),y(ny,0:lev),z(nz,0:lev),xg(rgn,0:lev),yg(rgn,0:lev),zg(rgn,0:lev),MM(rgn,0:lev),mxrho1(rgn,0:lev)&
       ,Phi(rgn,0:lev),Ixx(3,rgn,0:lev),Ixy(3,rgn,0:lev),Ixz(3,rgn,0:lev)&
       ,Iyy(3,rgn,0:lev),Iyz(3,rgn,0:lev),Izz(3,rgn,0:lev),Mss(rgn,0:lev,3))
  allocate(count1(rgn,0:lev),ig(rgn,0:lev),jg(rgn,0:lev),kg(rgn,0:lev)&
       ,igr(rgn,0:lev),jgr(rgn,0:lev),kgr(rgn,0:lev))

  count1(:,:)=0
  MM(:,:)=0.e0
  Mss(:,:,:)=0.e0
  Phi(:,:)=1.e5
  mxrho1(:,:)=0.e0
  ig(:,:)=0
  jg(:,:)=0
  kg(:,:)=0
  igr(:,:)=0
  jgr(:,:)=0
  kgr(:,:)=0
  xg(:,:)=0.e0
  yg(:,:)=0.e0
  zg(:,:)=0.e0
  Ixx(:,:,:)=0
  Ixy(:,:,:)=0
  Ixz(:,:,:)=0
  Iyy(:,:,:)=0
  Iyz(:,:,:)=0
  Izz(:,:,:)=0
  tagsph(:,:,:,:,:)=0

  nx2=nx
  ny2=ny
  nz2=nz
  do ilev=0,lev,1
     !write(*,*) ilev,'ilev'
     dl1=Lbox/real(nx2)
     x(1,ilev)=dl1*0.5e0
     y(1,ilev)=dl1*0.5e0
     z(1,ilev)=dl1*0.5e0
     write(*,*) dl1,ilev,x(1,ilev),y(1,ilev),z(1,ilev),'ilev'
     do i=2,nx2
        x(i,ilev)=x(i-1,ilev)+dl1
        write(*,*) ilev,i,x(i,ilev),'x'
     end do
     do j=2,ny2
        y(j,ilev)=y(j-1,ilev)+dl1
        write(*,*) ilev,j,y(j,ilev),'y'
     end do
     do k=2,nz2
        z(k,ilev)=z(k-1,ilev)+dl1
        write(*,*) ilev,k,z(k,ilev),'z'
     end do
     nx2=nx2/2
     ny2=ny2/2
     nz2=nz2/2
  end do

  tag(:,:,:)=0
  open(150,file='tagloghighden098001.DAT',FORM='FORMATTED')
  do k=1,nz
     do j=1,ny
        do i=1,nx
           read(150,*) tag(i,j,k)
        end do
     end do
     write(*,*) k,'read-tag'
  end do
  close(150)
  tag(0,:,:)=tag(nx,:,:)
  tag(:,0,:)=tag(:,ny,:)
  tag(:,:,0)=tag(:,:,nz)
  tag(nx+1,:,:)=tag(1,:,:)
  tag(:,ny+1,:)=tag(:,1,:)
  tag(:,:,nz+1)=tag(:,:,1)


  U(:,:,:,:)=0.e0
  open(110,file='Allnewbigtime098.DAT',FORM='UNFORMATTED')
  do k=1,nz
     do j=1,ny
        do i=1,nz
           read(110) U(i,j,k,1), U(i,j,k,2), U(i,j,k,3), U(i,j,k,4), U(i,j,k,5),&
                U(i,j,k,6), U(i,j,k,7), U(i,j,k,8), U(i,j,k,9)
        end do
     end do
     write(*,*) k, U(5,5,k,2), U(256,256,k,9),'read-phyval'
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


  ilev=0
  do k=1,nz
     do j=1,ny
        do i=1,nz
           n=tag(i,j,k)
           if(n.ne.0) then
              count1(n,ilev)=count1(n,ilev)+1
              MM(n,ilev)=MM(n,ilev)+(U(i,j,k,1))*(dl**3)
              if(U(i,j,k,9)<Phi(n,ilev))then
                 Phi(n,ilev)=U(i,j,k,9)
                 ig(n,ilev)=i
                 jg(n,ilev)=j
                 kg(n,ilev)=k
                 write(*,*)'ii,jj,kkP',i,j,k
              end if
              !mxrho1(n,ilev)=amax1(mxrho1(n,ilev),U(i,j,k,1))
              !write(*,*) U(i,j,k,1),mxrho1(n,ilev),'U(i,j,k,1)'
              if((U(i,j,k,1))>mxrho1(n,ilev))then
                 mxrho1(n,ilev)=(U(i,j,k,1))
                 igr(n,ilev)=i
                 jgr(n,ilev)=j
                 kgr(n,ilev)=k
                 write(*,*)'ii,jj,kkr',i,j,k
              end if
              xg(n,ilev)=(U(i,j,k,1))*x(i,ilev)*(dl**3)+xg(n,ilev)
              yg(n,ilev)=(U(i,j,k,1))*y(j,ilev)*(dl**3)+yg(n,ilev)
              zg(n,ilev)=(U(i,j,k,1))*z(k,ilev)*(dl**3)+zg(n,ilev)
           end if
        end do
     end do
     write(*,*)k,'ilev=0'
  end do

  do i=1,rgn
     if(MM(i,0).ne.0.e0)then
        xg(i,ilev)=xg(i,ilev)/MM(i,ilev)
        yg(i,ilev)=yg(i,ilev)/MM(i,ilev)
        zg(i,ilev)=zg(i,ilev)/MM(i,ilev)
        write(*,*)'GGG',xg(i,ilev),yg(i,ilev),zg(i,ilev)
     end if
  end do
  write(*,*)'endrgn,ilev=0'


  nx1=nx
  ny1=ny
  nz1=nz
  do ilev=1,lev,1
     nx2=nx1/2
     ny2=ny1/2
     nz2=nz1/2
     allocate(Udummy(-1:nx1+2,-1:ny1+2,-1:nz1+2,val))
     allocate(tagdummy(0:nx1+1,0:ny1+1,0:nz1+1))
     Udummy(:,:,:,:)=0.e0
     tagdummy(:,:,:)=0
     write(*,*)ilev,'allocation'
     if(ilev==1)then
        Udummy(:,:,:,:)=U(:,:,:,:)
        tagdummy(:,:,:)=tag(:,:,:)
        allocate(Usave(-1:nx2+2,-1:ny2+2,-1:nz2+2,val))
        allocate(tagsave(0:nx2+1,0:ny2+1,0:nz2+1))
        Usave(:,:,:,:)=0.e0
        tagsave(:,:,:)=0
        call smooth(Udummy,Usave,nx1,ny1,nz1,nx2,ny2,nz2,val)
        call smoothtag(tagdummy,tagsave,Usave,nx1,ny1,nz1,nx2,ny2,nz2,shd)
        deallocate(Udummy)
        deallocate(tagdummy)
     else
        !allocate(Usave(-1:nx2+2,-1:ny2+2,-1:nz2+2,val))
        Udummy(:,:,:,:)=Usave(:,:,:,:)
        tagdummy(:,:,:)=tagsave(:,:,:)
        deallocate(Usave)
        deallocate(tagsave)
        allocate(Usave(-1:nx2+2,-1:ny2+2,-1:nz2+2,val))
        allocate(tagsave(0:nx2+1,0:ny2+1,0:nz2+1))
        call smooth(Udummy,Usave,nx1,ny1,nz1,nx2,ny2,nz2,val)
        call smoothtag(tagdummy,tagsave,Usave,nx1,ny1,nz1,nx2,ny2,nz2,shd)
        deallocate(Udummy)
        deallocate(tagdummy)
     end if
     write(*,*)ilev,'endsmooth'

     dl1=Lbox/real(nx2)
     do k=1,nz2
        do j=1,ny2
           do i=1,nz2
              n=tagsave(i,j,k)
              if(n.ne.0) then
                 count1(n,ilev)=count1(n,ilev)+1
                 MM(n,ilev)=MM(n,ilev)+(Usave(i,j,k,1))*(dl1**3)
                 if((Usave(i,j,k,9))<phi(n,ilev))then
                    phi(n,ilev)=(Usave(i,j,k,9))
                    ig(n,ilev)=i
                    jg(n,ilev)=j
                    kg(n,ilev)=k
                    write(*,*) i,j,k,'ijkP'
                 end if
                 !mxrho1(n,ilev)=amax1(mxrho1(n,ilev),Usave(i,j,k,1))
                 !write(*,*) Usave(i,j,k,1),mxrho1(n,ilev),'Usave(i,j,k,1)'
                 if((Usave(i,j,k,1))>mxrho1(n,ilev))then
                    mxrho1(n,ilev)=(Usave(i,j,k,1))
                    igr(n,ilev)=i
                    jgr(n,ilev)=j
                    kgr(n,ilev)=k
                    write(*,*) i,j,k,'ijkr'
                 end if
                 xg(n,ilev)=(Usave(i,j,k,1))*x(i,ilev)*(dl1**3)+xg(n,ilev)
                 yg(n,ilev)=(Usave(i,j,k,1))*y(j,ilev)*(dl1**3)+yg(n,ilev)
                 zg(n,ilev)=(Usave(i,j,k,1))*z(k,ilev)*(dl1**3)+zg(n,ilev)
              end if
           end do
        end do
        write(*,*) k,ilev,'lev'
     end do
     do i=1,rgn
        if(MM(i,ilev).ne.0.e0)then
           xg(i,ilev)=xg(i,ilev)/MM(i,ilev)
           yg(i,ilev)=yg(i,ilev)/MM(i,ilev)
           zg(i,ilev)=zg(i,ilev)/MM(i,ilev)
           write(*,*)'GGG2',xg(i,ilev),yg(i,ilev),zg(i,ilev)
        end if
     end do
     nx1=nx2
     ny1=ny2
     nz1=nz2
     write(*,*) ilev,'ilev'
  end do

  !Mass(10pc)
  !L(ian)

  do ilev=0,lev
     do k=1,nz
        do j=1,ny
           do i=1,nz
              !if()then
              n=tag(i,j,k)
              !Mss(n,ilev)=Mss(n,ilev)+(U(i,j,k,1))*(dl**3)
              !end if
              if(n.ne.0) then
              lx=abs(x(i,0)-x(igr(n,ilev),ilev))
              ly=abs(y(j,0)-y(jgr(n,ilev),ilev))
              lz=abs(z(k,0)-z(kgr(n,ilev),ilev))
              !r=sqrt(lx**2+ly**2+lz**2)
              write(*,*)x(igr(n,ilev),ilev),y(jgr(n,ilev),ilev),z(kgr(n,ilev),ilev),igr(n,ilev),jgr(n,ilev),kgr(n,ilev),'r'
              !if(r<rshd) then
              !   Mss(n,ilev,1)=Mss(n,ilev,1)+(U(i,j,k,1))*(dl**3)
              !end if
              Ixx(1,n,ilev)=(U(i,j,k,1)*(ly**2+lz**2))*dl**3+Ixx(1,n,ilev)
              Ixy(1,n,ilev)= -(U(i,j,k,1)*(lx*ly))*dl**3+Ixy(1,n,ilev)
              Ixz(1,n,ilev)= -(U(i,j,k,1)*(lx*lz))*dl**3+Ixz(1,n,ilev)
              !Iyx
              Iyy(1,n,ilev)=(U(i,j,k,1)*(lx**2+lz**2))*dl**3+Iyy(1,n,ilev)
              Iyz(1,n,ilev)= -(U(i,j,k,1)*(lz*ly))*dl**3+Iyz(1,n,ilev)
              !Izx
              !Izy
              Izz(1,n,ilev)=(U(i,j,k,1)*(ly**2+lx**2))*dl**3+Izz(1,n,ilev)

              lx=abs(x(i,0)-x(ig(n,ilev),ilev))
              ly=abs(y(j,0)-y(jg(n,ilev),ilev))
              lz=abs(z(k,0)-z(kg(n,ilev),ilev))
              write(*,*)x(ig(n,ilev),ilev),y(jg(n,ilev),ilev),z(kg(n,ilev),ilev),ig(n,ilev),jg(n,ilev),kg(n,ilev),'P'
              !r=sqrt(lx**2+ly**2+lz**2)
              !if(r<rshd) then
              !   Mss(n,ilev,2)=Mss(n,ilev,2)+(U(i,j,k,1))*(dl**3)
              !end if
              Ixx(2,n,ilev)=(U(i,j,k,1)*(ly**2+lz**2))*dl**3+Ixx(2,n,ilev)
              Ixy(2,n,ilev)= -(U(i,j,k,1)*(lx*ly))*dl**3+Ixy(2,n,ilev)
              Ixz(2,n,ilev)= -(U(i,j,k,1)*(lx*lz))*dl**3+Ixz(2,n,ilev)
              !Iyx
              Iyy(2,n,ilev)=(U(i,j,k,1)*(lx**2+lz**2))*dl**3+Iyy(2,n,ilev)
              Iyz(2,n,ilev)= -(U(i,j,k,1)*(lz*ly))*dl**3+Iyz(2,n,ilev)
              !Izx
              !Izy
              Izz(2,n,ilev)=(U(i,j,k,1)*(ly**2+lx**2))*dl**3+Izz(2,n,ilev)

              lx=abs(x(i,0)-xg(n,ilev))
              ly=abs(y(j,0)-yg(n,ilev))
              lz=abs(z(k,0)-zg(n,ilev))
              !r=sqrt(lx**2+ly**2+lz**2)
              write(*,*)xg(n,ilev),yg(n,ilev),zg(n,ilev),'G'
              !if(r<rshd) then
              !   Mss(n,ilev,3)=Mss(n,ilev,3)+(U(i,j,k,1))*(dl**3)
              !end if
              Ixx(3,n,ilev)=(U(i,j,k,1)*(ly**2+lz**2))*dl**3+Ixx(3,n,ilev)
              Ixy(3,n,ilev)= -(U(i,j,k,1)*(lx*ly))*dl**3+Ixy(3,n,ilev)
              Ixz(3,n,ilev)= -(U(i,j,k,1)*(lx*lz))*dl**3+Ixz(3,n,ilev)
              !Iyx
              Iyy(3,n,ilev)=(U(i,j,k,1)*(lx**2+lz**2))*dl**3+Iyy(3,n,ilev)
              Iyz(3,n,ilev)= -(U(i,j,k,1)*(lz*ly))*dl**3+Iyz(3,n,ilev)
              !Izx
              !Izy
              Izz(3,n,ilev)=(U(i,j,k,1)*(ly**2+lx**2))*dl**3+Izz(3,n,ilev)
              end if
           end do
        end do
        write(*,*)k, ilev,'ilev2'
     end do
  end do

  do ilev=0,lev
     !write(ch,'(i3.3)') ilev
     !open(710,file='tagshp-a'//ch//'.DAT',FORM='FORMATTED')
     !open(720,file='tagshp-b'//ch//'.DAT',FORM='FORMATTED')
     !open(730,file='tagshp-c'//ch//'.DAT',FORM='FORMATTED')
     do ii=1,rgn
        if(MM(ii,ilev)*Msun>Mshd)then
           do k=1,nz
              do j=1,ny
                 do i=1,nz
                    lx=abs(x(i,0)-x(igr(ii,ilev),ilev))
                    ly=abs(y(j,0)-y(jgr(ii,ilev),ilev))
                    lz=abs(z(k,0)-z(kgr(ii,ilev),ilev))
                    r=sqrt(lx**2+ly**2+lz**2)
                    if(r<rshd) then
                       Mss(ii,ilev,1)=Mss(ii,ilev,1)+(U(i,j,k,1))*(dl**3)
                       tagsph(i,j,k,ilev,1)= ii
                    end if

                    lx=abs(x(i,0)-x(ig(ii,ilev),ilev))
                    ly=abs(y(j,0)-y(jg(ii,ilev),ilev))
                    lz=abs(z(k,0)-z(kg(ii,ilev),ilev))
                    r=sqrt(lx**2+ly**2+lz**2)
                    if(r<rshd) then
                       Mss(ii,ilev,2)=Mss(ii,ilev,2)+(U(i,j,k,1))*(dl**3)
                       tagsph(i,j,k,ilev,2)= ii
                    end if

                    lx=abs(x(i,0)-xg(ii,ilev))
                    ly=abs(y(j,0)-yg(ii,ilev))
                    lz=abs(z(k,0)-zg(ii,ilev))
                    r=sqrt(lx**2+ly**2+lz**2)
                    if(r<rshd) then
                       Mss(ii,ilev,3)=Mss(ii,ilev,3)+(U(i,j,k,1))*(dl**3)
                       tagsph(i,j,k,ilev,3)= ii
                    end if
                 end do
              end do
              write(*,*)k, ilev,'ilev2'
           end do
        end if
     end do
     !close(710)
     !close(720)
     !close(730)
  end do
  Mss(:,:,:)=Mss(:,:,:)*Msun
  MM(:,:)=MM(:,:)*Msun


  open(810,file='xg.DAT',FORM='FORMATTED')
  do i=1,rgn
     if(MM(i,0)>Mshd) then
        Write(num,'(i3.3)') i
        open(910,file='Imoment'//num//'.DAT',FORM='FORMATTED')
        !write(810,*) i,MM(i),x(ig(i,ilev),ilev),y(jg(i,ilev),ilev),z(kg(i,ilev),ilev)&
        !     ,x(igr(i,ilev),ilev),y(igr(i,ilev),ilev),z(igr(i,ilev),ilev),xg(i,ilev),yg(i,ilev),zg(i,ilev)
        write(910,*) i,MM(i,0)
        do ilev=0,lev

           write(810,*) i,count1(i,ilev),MM(i,ilev),Mss(i,ilev,1),Mss(i,ilev,2),Mss(i,ilev,3)&
                ,x(ig(i,ilev),ilev),y(jg(i,ilev),ilev),z(kg(i,ilev),ilev)&
                ,x(igr(i,ilev),ilev),y(jgr(i,ilev),ilev),z(kgr(i,ilev),ilev),xg(i,ilev),yg(i,ilev),zg(i,ilev)


           write(910,*) 'ilev=',ilev
           write(910,*) Ixx(1,i,ilev)*Msun/MM(i,0),Ixy(1,i,ilev)*Msun/MM(i,0),Ixz(1,i,ilev)*Msun/MM(i,0)
           write(910,*) Ixy(1,i,ilev)*Msun/MM(i,0),Iyy(1,i,ilev)*Msun/MM(i,0),Iyz(1,i,ilev)*Msun/MM(i,0)
           write(910,*) Ixz(1,i,ilev)*Msun/MM(i,0),Iyz(1,i,ilev)*Msun/MM(i,0),Izz(1,i,ilev)*Msun/MM(i,0)
           write(910,*)
           write(910,*) Ixx(2,i,ilev)*Msun/MM(i,0),Ixy(2,i,ilev)*Msun/MM(i,0),Ixz(2,i,ilev)*Msun/MM(i,0)
           write(910,*) Ixy(2,i,ilev)*Msun/MM(i,0),Iyy(2,i,ilev)*Msun/MM(i,0),Iyz(2,i,ilev)*Msun/MM(i,0)
           write(910,*) Ixz(2,i,ilev)*Msun/MM(i,0),Iyz(2,i,ilev)*Msun/MM(i,0),Izz(2,i,ilev)*Msun/MM(i,0)
           write(910,*)
           write(910,*) Ixx(3,i,ilev)*Msun/MM(i,0),Ixy(3,i,ilev)*Msun/MM(i,0),Ixz(3,i,ilev)*Msun/MM(i,0)
           write(910,*) Ixy(3,i,ilev)*Msun/MM(i,0),Iyy(3,i,ilev)*Msun/MM(i,0),Iyz(3,i,ilev)*Msun/MM(i,0)
           write(910,*) Ixz(3,i,ilev)*Msun/MM(i,0),Iyz(3,i,ilev)*Msun/MM(i,0),Izz(3,i,ilev)*Msun/MM(i,0)
           write(910,*)
           !write(910,*) 'ilev=',ilev
           !write(910,*) Ixx(1,i,ilev)*Msun,Ixy(1,i,ilev)*Msun,Ixz(1,i,ilev)*Msun
           !write(910,*) Ixy(1,i,ilev)*Msun,Iyy(1,i,ilev)*Msun,Iyz(1,i,ilev)*Msun
           !write(910,*) Ixz(1,i,ilev)*Msun,Iyz(1,i,ilev)*Msun,Izz(1,i,ilev)*Msun
           !write(910,*)
           !write(910,*) Ixx(2,i,ilev)*Msun,Ixy(2,i,ilev)*Msun,Ixz(2,i,ilev)*Msun
           !write(910,*) Ixy(2,i,ilev)*Msun,Iyy(2,i,ilev)*Msun,Iyz(2,i,ilev)*Msun
           !write(910,*) Ixz(2,i,ilev)*Msun,Iyz(2,i,ilev)*Msun,Izz(2,i,ilev)*Msun
           !write(910,*)
           !write(910,*) Ixx(3,i,ilev)*Msun,Ixy(3,i,ilev)*Msun,Ixz(3,i,ilev)*Msun
           !write(910,*) Ixy(3,i,ilev)*Msun,Iyy(3,i,ilev)*Msun,Iyz(3,i,ilev)*Msun
           !write(910,*) Ixz(3,i,ilev)*Msun,Iyz(3,i,ilev)*Msun,Izz(3,i,ilev)*Msun
           !write(910,*)
           !write(910,*) 'ilev=',ilev
           !write(910,*) Ixx(1,ilev),Ixy(1,ilev),Ixz(1,ilev)
           !write(910,*) Ixy(1,ilev),Iyy(1,ilev),Iyz(1,ilev)
           !write(910,*) Ixz(1,ilev),Iyz(1,ilev),Izz(1,ilev)
           !write(910,*)
           !write(910,*) Ixx(2,ilev),Ixy(2,ilev),Ixz(2,ilev)
           !write(910,*) Ixy(2,ilev),Iyy(2,ilev),Iyz(2,ilev)
           !write(910,*) Ixz(2,ilev),Iyz(2,ilev),Izz(2,ilev)
           !write(910,*)
           !write(910,*) Ixx(3,ilev),Ixy(3,ilev),Ixz(2,ilev)
           !write(910,*) Ixy(3,ilev),Iyy(3,ilev),Iyz(2,ilev)
           !write(910,*) Ixz(3,ilev),Iyz(3,ilev),Izz(2,ilev)
           !write(910,*)
        end do
        close(910)
        write(*,*)i,'rgn'
     end if
  end do
  close(810)

  do ilev=0,lev
     write(ch,'(i3.3)') ilev
     open(710,file='tagshp-a'//ch//'.DAT',FORM='FORMATTED')
     open(720,file='tagshp-b'//ch//'.DAT',FORM='FORMATTED')
     open(730,file='tagshp-c'//ch//'.DAT',FORM='FORMATTED')
     do k=1,nz
        do j=1,ny
           do i=1,nz
              write(710,*) tagsph(i,j,k,ilev,1)
              write(720,*) tagsph(i,j,k,ilev,2)
              write(730,*) tagsph(i,j,k,ilev,3)
           end do
        end do
        write(*,*)k, ilev,'tag'
     end do
     close(710)
     close(720)
     close(730)
  end do


  write(*,*)'end'
  deallocate(Usave,U,tag,tagsave,Udummy,tagdummy)
  deallocate(xg,yg,zg,MM,mxrho1,Phi,Ixx,Ixy,Ixz,Iyy,Iyz,Izz)
  deallocate(count1,ig,jg,kg,igr,jgr,kgr,x,y,z,Mss)

end program main



subroutine smooth(Uin,Uout,nxin,nyin,nzin,nxout,nyout,nzout,val)
  implicit none
  real(4),dimension(-1:nxin+2,-1:nyin+2,-1:nzin+2,val) :: Uin
  real(4),dimension(-1:nxout+2,-1:nyout+2,-1:nzout+2,val) :: Uout
  integer nxin,nyin,nzin,nxout,nyout,nzout,i,j,k,ii,jj,kk,val

  do k=1,nzout
     do j=1,nyout
        do i=1,nxout
           Uout(i,j,k,:) = (Uin(i*2-1,j*2-1,k*2-1,:)+Uin(i*2,j*2,k*2,:)+Uin(i*2,j*2-1,k*2-1,:)+Uin(i*2-1,j*2,k*2-1,:)&
                +Uin(i*2-1,j*2-1,k*2,:)+Uin(i*2,j*2,k*2-1,:)+Uin(i*2,j*2-1,k*2,:)+Uin(i*2-1,j*2,k*2,:))/8.d0
        end do
     end do
     !write(*,*)'sm',Uout(nxout/2,nyout/2,k,1)
  end do
  Uout(-1,:,:,:)=Uout(nxout-1,:,:,:)
  Uout(0,:,:,:)=Uout(nxout,:,:,:)
  Uout(nxout+1,:,:,:)=Uout(1,:,:,:)
  Uout(nxout+2,:,:,:)=Uout(2,:,:,:)
  Uout(:,-1,:,:)=Uout(:,nyout-1,:,:)
  Uout(:,0,:,:)=Uout(:,nyout,:,:)
  Uout(:,nyout+1,:,:)=Uout(:,1,:,:)
  Uout(:,nyout+2,:,:)=Uout(:,2,:,:)
  Uout(:,:,-1,:)=Uout(:,:,nzout-1,:)
  Uout(:,:,0,:)=Uout(:,:,nzout,:)
  Uout(:,:,nzout+1,:)=Uout(:,:,1,:)
  Uout(:,:,nzout+2,:)=Uout(:,:,2,:)
end subroutine smooth

subroutine smoothtag(Uin,Uout,rho,nxin,nyin,nzin,nxout,nyout,nzout,shd)
  implicit none
  !double precision,allocatable,dimension(:,:,:) :: Uin,Uout
  integer,dimension(0:nxin+1,0:nyin+1,0:nzin+1) :: Uin
  integer,dimension(0:nxout+1,0:nyout+1,0:nzout+1) :: Uout
  real(4),dimension(-1:nxout+2,-1:nyout+2,-1:nzout+2) :: rho
  real(4) shd
  integer nxin,nyin,nzin,nxout,nyout,nzout,i,j,k,ii,jj,kk,val,maxint

  do k=1,nzout
     do j=1,nyout
        do i=1,nxout
           maxint=0
           maxint=max0(maxint,Uin(i*2-1,j*2-1,k*2-1),Uin(i*2,j*2,k*2),Uin(i*2,j*2-1,k*2-1),Uin(i*2-1,j*2,k*2-1)&
                ,Uin(i*2-1,j*2-1,k*2),Uin(i*2,j*2,k*2-1),Uin(i*2,j*2-1,k*2),Uin(i*2-1,j*2,k*2))
           Uout(i,j,k) = (Uin(i*2-1,j*2-1,k*2-1)+Uin(i*2,j*2,k*2)+Uin(i*2,j*2-1,k*2-1)+Uin(i*2-1,j*2,k*2-1)&
                +Uin(i*2-1,j*2-1,k*2)+Uin(i*2,j*2,k*2-1)+Uin(i*2,j*2-1,k*2)+Uin(i*2-1,j*2,k*2))/8
           if(maxint.ne.Uout(i,j,k)) then
              if(rho(i,j,k)>shd) then
                 Uout(i,j,k)=maxint
              else
                 Uout(i,j,k)=0
              end if
           end if
        end do
     end do
  end do
  Uout(0,:,:)=Uout(nxout,:,:)
  Uout(nxout+1,:,:)=Uout(1,:,:)
  Uout(:,0,:)=Uout(:,nyout,:)
  Uout(:,nyout+1,:)=Uout(:,1,:)
  Uout(:,:,0)=Uout(:,:,nzout)
  Uout(:,:,nzout+1)=Uout(:,:,1)
end subroutine smoothtag



