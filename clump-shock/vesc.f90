program main
  implicit none
  integer,allocatable,dimension(:,:,:) :: tag,tagtrace,tagtracesf,tagswp
  integer,allocatable,dimension(:) :: zeroone,ii,jj,kk
  real(4), allocatable :: u(:,:,:,:),v(:,:,:,:),sfgrv(:),vesc(:),Uswp(:,:,:)
  integer val,times,lasttime,ig,jg,kg,numtag,igg,jgg,kgg,countrsph
  integer :: i,j,k,kz=1,initime,loopnum1,loopnum,gg,iimax,iimin,jjmax,jjmin,kkmax,kkmin,ll,mm,nn,nx,ny,nz,m,n
  real(4) :: G,Lbox,mu,dxyz,Msun,vescmean,v0,vescmean2,Mass,Phimin,lx,ly,lz,xg,yg,zg,rsph,Massrsph,vescmeanrsph
  character(3) NPE,time,num
  character(55) :: dir='/glv0/maedarn/clst-form-HIcol/cnv100-wb-wg-sm-optthick/'

!----parameter----
nx=512
ny=512
nz=512
val=18
Lbox=100.e0
Msun=2.4d-2 !1pc * 1pc * 1pc * 1m_p/cc
G=0.00011142e0
!L=100.0e0
dxyz=100.e0/(512.e0)
m=98
v0=0.953e0
numtag=9
rsph=5.e0
!----parameter----

 allocate(U(-1:nx+2,-1:ny+2,-1:nz+2,val))
 allocate(tag(0:nx+1,0:ny+1,0:nz+1),tagtrace(0:nx+1,0:ny+1,0:nz+1),tagtracesf(0:nx+1,0:ny+1,0:nz+1))
 allocate(zeroone(0:20000))
 allocate(Uswp(-1:nx+2,-1:ny+2,-1:nz+2))
 allocate(tagswp(0:nx+1,0:ny+1,0:nz+1))


  zeroone(:)=0
  write(num,'(i3.3)') m
  open(150,file=dir//'tag/tagn'//num//'001.DAT',access='stream',FORM='UNFORMATTED')
  write(*,*) num
  !m=1
  do k=1,nz
     do j=1,ny
        !write(*,*) i,j,k,m,'read-tag'
        do i=1,nx
           !write(*,*) i,j,k,m,'read-tag'
           read(150) tag(i,j,k)!,dm,dm1
           !write(*,*) i,j,k,tag(i,j,k,m)
        end do
     end do
  end do
  close(150)

  U(:,:,:,:)=0.e0
  !open(110,file='/home/maedarn/cnv100wbg-cluster/Allhirshighden098.DAT',FORM='FORMATTED')
  open(110,file=dir//'All/All'//num//'.DAT',access='stream',FORM='UNFORMATTED')
  do k=1,nz
     do j=1,ny
        do i=1,nx
           read(110) (U(i,j,k,n),n=1,val)
        end do
     end do
     write(*,*) k, U(5,5,k,2),'read-phyval'
  end do

goto 290
do k=1,nz
   do j=1,ny
      do i=1,nx
         kz=k+(512-450)-512*((+512-450+k)/512)
         Uswp(i,j,kz)=U(i,j,k,1)
         tagswp(i,j,kz)=tag(i,j,k)
      end do
   end do
   write(*,*) k, U(5,5,k,2),'read-phyval'
end do
U(:,:,:,1)=Uswp(:,:,:)
tag(:,:,:)=tagswp(:,:,:)
290 continue

U(:,:,-1,:)=U(:,:,nz-2,:)
U(:,:, 0,:)=U(:,:,nz-1,:)
U(:,:,nz+1,:)=U(:,:,1,:)
U(:,:,nz+2,:)=U(:,:,2,:)
U(:,-1,:,:)=U(:,ny-2,:,:)
U(:, 0,:,:)=U(:,ny-1,:,:)
U(:,ny+1,:,:)=U(:,1,:,:)
U(:,ny+2,:,:)=U(:,2,:,:)
U(-1,:,:,:)=U(1,:,:,:)
U( 0,:,:,:)=U(1,:,:,:)
U(nx+1,:,:,:)=U(nx,:,:,:)
U(nx+2,:,:,:)=U(nx,:,:,:)


tag(0,:,:)=tag(nx,:,:)
tag(:,0,:)=tag(:,ny,:)
tag(:,:,0)=tag(:,:,nz)
tag(nx+1,:,:)=tag(1,:,:)
tag(:,ny+1,:)=tag(:,1,:)
tag(:,:,nz+1)=tag(:,:,1)


  zeroone(numtag)=1
  do k=0,nz+1
    do j=0,ny+1
      do i=0,nz+1
        tagtrace(i,j,k)=zeroone(tag(i,j,k))
      end do
    end do
  end do


zeroone(:)=0
zeroone(1)=1
zeroone(2)=1
zeroone(3)=1
zeroone(4)=1
zeroone(5)=1
zeroone(6)=1

tagtracesf(:,:,:)=0
do k=1,nz
  do j=1,ny
    do i=1,nz
      tagtracesf(i,j,k)=zeroone(tagtrace(i,j,k)+tagtrace(i+1,j,k)+tagtrace(i-1,j,k) &
             +tagtrace(i,j+1,k)+tagtrace(i,j-1,k)+tagtrace(i,j,k+1)+tagtrace(i,j,k-1))
      tagtracesf(i,j,k)=tagtracesf(i,j,k)*tagtrace(i,j,k)
    end do
  end do
end do
zeroone(:)=0

loopnum=sum(tagtracesf)

allocate(ii(1:loopnum),jj(1:loopnum),kk(1:loopnum))

loopnum1=1
Mass=0.e0
Phimin=1.e5
xg=0.e0
yg=0.e0
zg=0.e0
do k=1,nz
  do j=1,ny
    do i=1,nz
    if(tagtracesf(i,j,k)==1) then
    ii(loopnum1)=i
    jj(loopnum1)=j
    kk(loopnum1)=k
    loopnum1=loopnum1+1
    endif
    if(tagtrace(i,j,k)==1) then
    Mass=Mass+U(i,j,k,1)*dxyz*dxyz*dxyz*Msun
    xg=U(i,j,k,1)*dxyz*dxyz*dxyz*(dxyz*float(i))+xg
    yg=U(i,j,k,1)*dxyz*dxyz*dxyz*(dxyz*float(j))+yg
    zg=U(i,j,k,1)*dxyz*dxyz*dxyz*(dxyz*float(k))+zg
      if(U(i,j,k,17)<Phimin) then
      Phimin=(U(i,j,k,17))
      ig=i
      jg=j
      kg=k
      endif
    end if
    end do
  end do
end do

xg=xg/(Mass/Msun)
yg=yg/(Mass/Msun)
zg=zg/(Mass/Msun)
igg=nint(xg/dxyz)
jgg=nint(yg/dxyz)
kgg=nint(zg/dxyz)


iimax=maxval(ii)
iimin=minval(ii)
jjmax=maxval(jj)
jjmin=minval(jj)
kkmax=maxval(kk)
kkmin=minval(kk)

allocate(sfgrv(1:loopnum),vesc(1:loopnum))
do gg=1,loopnum
write(*,*) gg,'loopnum'
sfgrv(gg)=G*U(ii(gg),jj(gg),kk(gg),1)*float(tagtrace(ii(gg),jj(gg),kk(gg)))*dxyz*dxyz*dxyz/sqrt(1.e-10)
do nn=kkmin,kkmax
  do mm=jjmin,jjmax
    do ll=iimin,iimax
     lx=abs(float(ll)-float(ii(gg)))
     ly=abs(float(mm)-float(jj(gg)))
     lz=abs(float(nn)-float(kk(gg)))
!     lx=amin1(lx,abs(lx-float(nx)))
!     ly=amin1(ly,abs(ly-float(ny)))
!     lz=amin1(lz,abs(ly-float(nz)))
    !sfgrv(gg)=sfgrv(gg)-G*U(ll,mm,nn,1)*float(tagtrace(ll,mm,nn))*dxyz*dxyz*dxyz/sqrt(float((ll-ii(gg))**2+(mm-jj(gg))**2+(nn-kk(gg))**2)*dxyz**2+1.e-10)
    sfgrv(gg)=sfgrv(gg)-G*U(ll,mm,nn,1)*float(tagtrace(ll,mm,nn))*dxyz*dxyz*dxyz/sqrt(((lx)**2+(ly)**2+(lz)**2)*dxyz**2+1.e-10)
    end do
  end do
end do
enddo

vescmean=0.e0
vescmean2=0.e0
do gg=1,loopnum
vesc(gg)=sqrt(-2.0e0*sfgrv(gg))
vescmean=vescmean+sqrt(-2.0e0*sfgrv(gg))
vescmean2=vescmean2+vesc(gg)**2
enddo
vescmean=vescmean/float(loopnum)
vescmean2=sqrt(vescmean2/float(loopnum))

open(910,file=dir//'AllHDF/SF'//num//'.DAT',access='stream',FORM='UNFORMATTED')
do k=1,nz
  do j=1,ny
    do i=1,nz
      write(910) float(tagtracesf(i,j,k)),float(tagtrace(i,j,k)),U(i,j,k,1)
    end do
  end do
end do
write(*,*)k,'write'
close(910)

open(920,file=dir//'AllHDF/SFPhi-R'//num//'.DAT',access='stream',FORM='UNFORMATTED')
do gg=1,loopnum
lx=abs(float(ig)-float(ii(gg)))
ly=abs(float(jg)-float(jj(gg)))
lz=abs(float(kg)-float(kk(gg)))
!lx=amin1(lx,abs(lx-float(nx)))
!ly=amin1(ly,abs(ly-float(ny)))
!lz=amin1(lz,abs(ly-float(nz)))
write(920) sqrt(((lx)**2+(ly)**2+(lz)**2)*dxyz**2),vesc(gg)
write(*,*)gg,'write-loopnum'
end do
close(920)

write(*,*) vescmean*v0,vescmean2*v0,Mass,loopnum

Massrsph=0.e0
vescmeanrsph=0.e0
do k=1,nz
  do j=1,ny
    do i=1,nz
    lx=abs(float(ig)-float(i))
    ly=abs(float(jg)-float(j))
    lz=abs(float(kg)-float(k))
    if(sqrt((lx**2+ly**2+lz**2)*dxyz**2)<rsph) then
    Massrsph=Massrsph+U(i,j,k,1)*dxyz*dxyz*dxyz*Msun*float(tagtrace(i,j,k))
    end if
    end do
  end do
end do

countrsph=0
do gg=1,loopnum
lx=abs(float(ig)-float(ii(gg)))
ly=abs(float(jg)-float(jj(gg)))
lz=abs(float(kg)-float(kk(gg)))
if(sqrt((lx**2+ly**2+lz**2)*dxyz**2)<rsph) then
vescmeanrsph=vescmeanrsph+sqrt(-2.0e0*sfgrv(gg))
countrsph=countrsph+1
end if
enddo
vescmeanrsph=vescmeanrsph/float(countrsph)
    
write(*,*)'rsph=',rsph,Massrsph,vescmeanrsph


deallocate(U)
deallocate(tag,tagtrace,tagtracesf)
deallocate(zeroone)
deallocate(ii,jj,kk)
deallocate(sfgrv,vesc)
deallocate(Uswp,tagswp)

end program main

