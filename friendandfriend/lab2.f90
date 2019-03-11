module common
  implicit none
  integer,parameter :: nx=512,ny=512,nz=512
  integer,parameter :: rank=nx*ny*nz/4
  integer,dimension(0:rank) :: rl_table,next_table,tail_table,T
end module common

program main
  use common
  implicit none
  integer,allocatable,dimension(:,:,:) :: tag,den
  real(4),allocatable,dimension(:,:,:) :: rho
  integer lab,l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13
  integer i,j,k,Vb,r,m,mm,in
  real(4) rhoshd

  !allocate(tag(1:nx,1:ny,1:nz))
  !allocate(den(1:nx,1:ny,1:nz))
  !allocate(rho(1:nx,1:ny,1:nz))

  !allocate(tag(0:nx,0:ny,0:nz))
  !allocate(den(0:nx,0:ny,0:nz))
  !allocate(rho(0:nx,0:ny,0:nz))]


  allocate(tag(0:nx+1,0:ny+1,0:nz+1))
  allocate(den(0:nx+1,0:ny+1,0:nz+1))
  allocate(rho(0:nx+1,0:ny+1,0:nz+1))

  lab=1
  Vb=0
  rhoshd=5000.e0
  tag(:,:,:)=0
  rho(:,:,:)=0.e0
  rl_table(:)=0
  next_table(:)=0
  tail_table(:)=0
  T(:)=0

!  goto 1000
  open(150,file='/Users/maeda/Desktop/kaiseki/cnv100wbwg/all.DAT',FORM='FORMATTED')
  do k=1,nz
     do j=1,ny
        do i=1,nx
           read(150,*) rho(i,j,k)
        end do
     end do
     write(*,*) rho(256,256,k),k
  end do
   !rho(:,:,-1)=rho(:,:,nz-2)
   rho(:,:, 0)=rho(:,:,nz-1)
   rho(:,:,nz+1)=rho(:,:,1)
   !rho(:,:,nz+2)=rho(:,:,2)
   !rho(:,-1,:)=rho(:,ny-2,:)
   rho(:, 0,:)=rho(:,ny-1,:)
   rho(:,ny+1,:)=rho(:,1,:)
   !rho(:,ny+2,:)=rho(:,2,:)
   close(150)
!1000 continue

   goto 1001
   write(*,*)'qq'
   do k=1,nz
      do j=1,ny
         do i=1,nx
            r=i+j+k
            !write(*,*) j,k
            if(r==66) then
               rho(i,j,k)=10000.e0
            else
               rho=0.e0
            end if
            !write(*,*) j,k
         end do
         write(*,*) j,k
      end do
      !write(*,*) rho(256,256,k),k
   end do
1001 continue
   !rho(29,29,29)=10000e0
   !rho(29,29,28)=10000e0
   !rho(29,28,29)=10000e0
  !do k=1,nz
  !   do j=1,ny
  !      do i=1,nx
  do k=0,nz+1
     do j=0,ny+1
        do i=0,nx+1
           if(rho(i,j,k)>rhoshd) then
              !rho(i,j,k)=1.e0
              den(i,j,k)=1
           else
              !rho(i,j,k)=0.d0
              den(i,j,k)=0
           end if
        end do
     end do
     write(*,*) k
  end do


!  do m=1,2
  do k=1,nz
  do j=1,ny
  do i=1,nx
  if(den(i,j,k).ne.Vb)then

  l1 =tag(i-1,j  ,k  )
  l2 =tag(i-1,j-1,k  )
  l3 =tag(i  ,j-1,k  )
  l4 =tag(i+1,j-1,k  )
  l5 =tag(i-1,j-1,k-1)
  l6 =tag(i  ,j-1,k-1)
  l7 =tag(i+1,j-1,k-1)
  l8 =tag(i-1,j  ,k-1)
  l9 =tag(i  ,j  ,k-1)
  l10=tag(i+1,j  ,k-1)
  l11=tag(i-1,j+1,k-1)
  l12=tag(i  ,j+1,k-1)
  l13=tag(i+1,j+1,k-1)

  !write(*,*) l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,i,j,k,lab
  if(l9>0)then
     tag(i,j,k)=l9
  elseif(l3>0) then
     tag(i,j,k)=l3
     if((l12>0).and.(l8==0).and.(l10==0))then
        call resolve(tag(i  ,j-1,k  ),tag(i  ,j+1,k-1))
     else
        if((l11>0).and.(l8==0))then
           call resolve(tag(i  ,j-1,k  ),tag(i-1,j+1,k-1))
        end if
        if((l13>0).and.(l10==0))then
           call resolve(tag(i  ,j-1,k  ),tag(i+1,j+1,k-1))
        end if
     end if
  elseif(l6>0) then
     tag(i,j,k)=l6
     if((l12>0).and.(l8==0).and.(l10==0))then
        call resolve(tag(i  ,j-1,k-1),tag(i  ,j+1,k-1))
     else
        if((l11>0).and.(l8==0))then
           call resolve(tag(i  ,j-1,k-1),tag(i-1,j+1,k-1))
        end if
        if((l13>0).and.(l10==0))then
           call resolve(tag(i  ,j-1,k-1),tag(i+1,j+1,k-1))
        end if
     end if
  elseif(l1>0) then
     tag(i,j,k)=l1
     if((l10>0).and.(l12==0))then
        call resolve(tag(i-1,j  ,k  ),tag(i+1,j  ,k-1))
     elseif(l7>0) then
        call resolve(l1,tag(i+1,j-1,k-1))
        if(l13>0) then
           call resolve(tag(i-1,j  ,k  ),tag(i+1,j+1,k-1))
        end if
     elseif(l4>0) then
        !write(*,*)'bag'
        call resolve(tag(i-1,j  ,k  ),tag(i+1,j-1,k  ))
        if(l13>0) then
           !write(*,*)'bag'
           call resolve(tag(i-1,j  ,k  ),tag(i+1,j+1,k-1))
        end if
     elseif((l13>0).and.(l12==0))then
        call resolve(tag(i-1,j  ,k  ),tag(i+1,j+1,k-1))
     end if
  elseif(l8>0) then
     tag(i,j,k)=l8
     if((l10>0).and.(l12==0))then
        call resolve(tag(i-1,j  ,k-1),tag(i+1,j  ,k-1))
     elseif(l7>0) then
        call resolve(tag(i-1,j  ,k-1),tag(i+1,j-1,k-1))
        if(l13>0) then
           call resolve(tag(i-1,j  ,k-1),tag(i+1,j+1,k-1))
        end if
     elseif(l4>0) then
        call resolve(tag(i-1,j  ,k-1),tag(i+1,j-1,k  ))
        if(l13>0) then
           call resolve(tag(i-1,j  ,k-1),tag(i+1,j+1,k-1))
        end if
     elseif((l13>0).and.(l12==0))then
        call resolve(tag(i-1,j  ,k-1),tag(i+1,j+1,k-1))
     end if
  elseif(l10>0) then
     tag(i,j,k)=l10
     if((l11>0).and.(l12==0))then
        call resolve(tag(i+1,j  ,k-1),tag(i-1,j+1,k-1))
     end if
     if(l5>0) then
        call resolve(tag(i+1,j  ,k-1),tag(i-1,j-1,k-1))
     elseif(l2>0) then
        call resolve(tag(i+1,j  ,k-1),tag(i-1,j-1,k  ))
     end if
  elseif(l12>0) then !elseif(l2>0) then?
     tag(i,j,k)=l12
     if(l4>0) then
        call resolve(tag(i  ,j+1,k-1),tag(i+1,j-1,k  ))
     elseif(l7>0) then
        call resolve(tag(i  ,j+1,k-1),tag(i+1,j-1,k-1))
     end if
     if(l2>0) then
        call resolve(tag(i  ,j+1,k-1),tag(i-1,j-1,k  ))
     elseif(l5>0) then
        call resolve(tag(i  ,j+1,k-1),tag(i-1,j-1,k-1))
     end if
  elseif(l5>0) then
     tag(i,j,k)=l5
     if(l4>0) then
        call resolve(tag(i-1,j-1,k-1),tag(i+1,j-1,k  ))
     elseif(l7>0) then
        call resolve(tag(i-1,j-1,k-1),tag(i+1,j-1,k-1))
     end if
     if(l11>0) then
        call resolve(tag(i-1,j-1,k-1),tag(i-1,j+1,k-1))
     elseif(l13>0) then
        call resolve(tag(i-1,j-1,k-1),tag(i+1,j+1,k-1))
     end if
  elseif(l2>0) then
     tag(i,j,k)=l2
     if(l4>0) then
        call resolve(tag(i-1,j-1,k  ),tag(i+1,j-1,k  ))
     elseif(l7>0) then
        call resolve(tag(i-1,j-1,k  ),tag(i+1,j-1,k-1))
     end if
     if(l11>0) then
        call resolve(tag(i-1,j-1,k  ),tag(i-1,j+1,k-1))
     elseif(l13>0) then
        call resolve(tag(i-1,j-1,k  ),tag(i+1,j+1,k-1))
     end if
  elseif(l4>0) then
     tag(i,j,k)=l4
     if(l11>0) then
        call resolve(tag(i+1,j-1,k  ),tag(i-1,j+1,k-1))
     elseif(l13>0) then
        call resolve(tag(i+1,j-1,k  ),tag(i+1,j+1,k-1))
     end if
  elseif(l7>0) then
     tag(i,j,k)=l7
     if(l11>0) then
        call resolve(tag(i+1,j-1,k-1),tag(i-1,j+1,k-1))
     elseif(l13>0) then
        !call resolve(tag(i-1,j-1,k  ),tag(i+1,j+1,k-1))
        call resolve(tag(i+1,j-1,k-1),tag(i+1,j+1,k-1))
     end if
  elseif(l11>0) then
     tag(i,j,k)=l11
     if(l13>0) then
        !call resolve(tag(i-1,j-1,k  ),tag(i+1,j+1,k-1))
        call resolve(tag(i-1,j+1,k-1),tag(i+1,j+1,k-1))
     end if
  elseif(l13>0) then
     tag(i,j,k)=l13

  else
     !lab=lab+1
     write(*,*) lab,'lab'
     tag(i,j,k)=lab
     T(lab)=lab
     rl_table(lab)=lab
     next_table(lab)= -1
     tail_table(lab)= lab
     lab=lab+1
  end if

 else
   tag(i,j,k)=Vb
 endif
! write(*,*) i,j,k,'ijk'
 enddo
 write(*,*) i,j,k
 enddo
 !write(*,*) k
 enddo

 tag(:,:, 0)=tag(:,:,nz)
 tag(:,:,nz+1)=tag(:,:,1)
 tag(:, 0,:)=tag(:,ny,:)
 tag(:,ny+1,:)=tag(:,1,:)
 !enddo


 !変換
 T(:)=rl_table(:)
 T(Vb)=0
 mm=0
 !do i=1,rank
 do i=1,3000
    !if(rl_table(i)==i) then
    in=1
    do j=1,3000
       if(rl_table(j)==i) then
          if(in==1) then
             mm=mm+1
          end if
          T(j)=mm
          write(*,*) mm,in,'mm'
          in=in+1
          !write(*,*) mm
       end if
    end do
 !   mm=mm+1
    !end if
 end do
 write(*,*)'separate region number = ' , mm
 !変換




  !T(Vb)=0
  do k=1,nz
     do j=1,ny
        do i=1,nx
           tag(i,j,k)=T(tag(i,j,k))
           !tag(i,j,k)=rl_table(tag(i,j,k))
        end do
     end do
     write(*,*) k
  end do

  open(290,file='/Users/maeda/Desktop/kaiseki/image/kakuninrl.DAT',FORM='FORMATTED')
  do k=1,3000
     write(290,*) k,rl_table(k),T(k)
  end do
  close(290)

  open(110,file='/Users/maeda/Desktop/kaiseki/image/sa-3D.DAT',FORM='FORMATTED')
  do k=1,nz
     do j=1,ny
        do i=1,nx
           !write(110,*) u(i,j,k)-v(i,j,k)
           write(110,*) tag(i,j,k),den(i,j,k),rho(i,j,k)
        end do
     end do
  end do
  deallocate(tag,rho,den)
end program main

!subroutine resolve(u,v)
subroutine resolve(a,b)
  !implicit none
  use common
  implicit none
  integer :: u,v,i,ik,count=0,a,b
  character(3) chr
  write(chr,'(I3.3)') count
  count=count+1
  open(190,file='/Users/maeda/Desktop/kaiseki/image/kakunin'//chr//'.DAT',FORM='FORMATTED')
  do ik=1,150
     write(190,*) ik,rl_table(ik),next_table(ik),tail_table(ik),a,b
  end do
  close(190)
  !write(*,*)'aa',u,v,next_table(tail_table(u)),next_table(tail_table(v)),tail_table(u)
  v=max0(rl_table(a),rl_table(b))
  u=min0(rl_table(a),rl_table(b))

  if(v==u) then
     goto 600
  end if
  next_table(tail_table(u))=v
  tail_table(u)=tail_table(v)
  !write(*,*)tail_table(u),tail_table(v),next_table(v),next_table(next_table(v))
  i=v
  write(*,*)'aa',u,v,next_table(tail_table(u)),tail_table(u)
  do while(i.ne.-1)
     write(*,*) u,v,i,next_table(i),next_table(v)!,next_table(31)
     rl_table(i)=u
     i=next_table(i)
!     write(*,*) i
  end do
  600 continue
  !write(*,*)'bb'
end subroutine resolve
