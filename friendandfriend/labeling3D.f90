module common
  implicit none
  integer , parameter ::inin=100
  integer table(inin)
end module common

program main
  !implicit none
  use common
  implicit none
  real(4), allocatable:: rho(:,:,:)
  real(4) rhoshd
  integer,allocatable :: v(:,:,:),u(:,:,:)!,table(100,100)
  integer i,j,k,tag,a,b,l,m,n,c(13),r1,r2,sh,sh2(100),bc(9)
  integer mesh,lengthx,lengthy,lengthz

  !parameter
  mesh=512+2
  lengthx=512
  lengthy=512
  lengthz=512
  rhoshd=1000.e0
  !parameter

  allocate(v(-1:mesh,-1:mesh,-1:mesh), u(-1:mesh,-1:mesh,-1:mesh))
  allocate(rho(-1:mesh,-1:mesh,-1:mesh))
  rho(:,:,:)=1.e0
  write(*,*)rho(1,1,1)
  open(150,file='/Users/maeda/Desktop/kaiseki/cnv100wbwg/all.DAT',FORM='FORMATTED')
   do k=1,lengthz
      do j=1,lengthy
         do i=1,lengthx
            read(150,*) rho(i,j,k)
            !write(*,*) rho(i,j,k)
         end do
      end do
      write(*,*) rho(256,256,k),k
   end do
   rho(:,:,-1)=rho(:,:,lengthz-2)
   rho(:,:, 0)=rho(:,:,lengthz-1)
   rho(:,:,lengthz+1)=rho(:,:,1)
   rho(:,:,lengthz+2)=rho(:,:,2)
   rho(:,-1,:)=rho(:,lengthy-2,:)
   rho(:, 0,:)=rho(:,lengthy-1,:)
   rho(:,lengthy+1,:)=rho(:,1,:)
   rho(:,lengthy+2,:)=rho(:,2,:)
   close(150)

   v(:,:,:)=0
   do k=-1,lengthz+2
      do j=-1,lengthy+2
         do i=-1,lengthx+2
            if(rho(i,j,k)>rhoshd) then
               !rho(i,j,k)=1.e0
               v(i,j,k)=1
            else
               !rho(i,j,k)=0.d0
               v(i,j,k)=0
            end if
         end do
      end do
      write(*,*) k,'int'
   end do

  table(:)=0
  !do j=1,100
  do i=1,inin
     table(i)=i
  end do
  !end do

!  v(:,:,:)=0
  u(:,:,:)=0
  !initial
  goto 999
  do j=1,100
     do i=1,100
      !  l=i+j
      !  m=i-j
      !  if(l==20 .or. l==21) then
      !     v(i,j)=1
      !  end if
      !  if(m==10 .or. m==11) then
      !     v(i,j)=1
      !  end if
        r1=int(sqrt(real((i-50)*(i-50)+(j-50)*(j-50))))
        if(r1==20 .or. r1==21) then
           v(i,j,k)=2
        end if
        if(r1==30 .or. r1==31) then
           v(i,j,k)=1
        end if
     end do
  end do
  open(10,file='/Users/maeda/Desktop/kaiseki/image/int3D.DAT',FORM='FORMATTED')
  do j=1,100
     do i=1,100
        write(10,*) v(i,j,k)
     end do
     write(10,*)
  end do
  close(10)
  999 continue
  !initial


  tag=1
  do k=1,mesh-2
     do j=1,mesh-2
        do i=1,mesh-2
           if(v(i,j,k) .ne. 0) then
              c(1) =u(i-1,j-1,k-1)
              c(2) =u(i  ,j-1,k-1)
              c(3) =u(i+1,j-1,k-1)
              c(4) =u(i-1,j  ,k-1)
              c(5) =u(i  ,j  ,k-1)
              c(6) =u(i+1,j  ,k-1)
              c(7) =u(i-1,j+1,k-1)
              c(8) =u(i  ,j+1,k-1)
              c(9) =u(i+1,j+1,k-1)
              c(10)=u(i-1,j-1,k  )
              c(11)=u(i  ,j-1,k  )
              c(12)=u(i+1,j-1,k  )
              c(13)=u(i-1,j  ,k  )
              a=c(1)+c(1)+c(2)+c(3)+c(4)+c(5)+c(6)+c(7)+c(8)+c(9)+c(10)+c(11)+c(12)+c(13)
              b=(c(1)-a)+(c(1)-a)+(c(2)-a)+(c(3)-a)+(c(4)-a)+(c(5)-a)+(c(6)-a)+(c(7)-a)&
                   +(c(8)-a)+(c(9)-a)+(c(10)-a)+(c(11)-a)+(c(12)-a)+(c(13)-a)
              !b=(u(i-1,j-1)-a)*(u(i,j-1)-a)*(u(i+1,j-1)-a)*(u(i-1,j)-a)
              if(a.ne.0) then
                 !u(i,j)=min0(-u(i-1,j-1),-u(i,j-1),-u(i+1,j-1),-u(i-1,j))
                 !u(i,j)=-u(i,j)
                 !u(i,j)=min0(u(i-1,j-1),u(i,j-1))
                 !write(*,*)u(i,j) ,'inin'
                 u(i,j,k)=inin*10
                 do n=1,13
                    if(c(n).ne.0)then
                       u(i,j,k)=min0(c(n),u(i,j,k))
                       !write(*,*) u(i,j) ,'00'
                    end if
                 end do
                 if(b.ne.0) then
                    call tbl(c(1),u(i,j,k))
                    ! write(*,*)'in'
                 end if
              else
                 u(i,j,k)=tag
                 tag=tag+1
                 write(*,*) 'tag ' ,tag
              end if
           end if
        end do
     end do
     write(*,*) k,'aaa'
  end do

  !-----bc------
  goto 1001
  k=mesh-2
  j=mesh-2
  i=mesh-2
  if(v(i,j,k) .ne. 0) then
     c(1) =u(i-1,j-1,k-1)
     c(2) =u(i  ,j-1,k-1)
     c(3) =u(i+1,j-1,k-1)
     c(4) =u(i-1,j  ,k-1)
     c(5) =u(i  ,j  ,k-1)
     c(6) =u(i+1,j  ,k-1)
     c(7) =u(i-1,j+1,k-1)
     c(8) =u(i  ,j+1,k-1)
     c(9) =u(i+1,j+1,k-1)
     c(10)=u(i-1,j-1,k  )
     c(11)=u(i  ,j-1,k  )
     c(12)=u(i+1,j-1,k  )
     c(13)=u(i-1,j  ,k  )
     !bc(1)
     !bc(2)
     !bc(3)
     !bc(4)
     !bc(5)
     !bc(6)
     !bc(7)
     !bc(8)
     !bc(9)
     a=c(1)+c(1)+c(2)+c(3)+c(4)+c(5)+c(6)+c(7)+c(8)+c(9)+c(10)+c(11)+c(12)+c(13)
     b=(c(1)-a)+(c(1)-a)+(c(2)-a)+(c(3)-a)+(c(4)-a)+(c(5)-a)+(c(6)-a)+(c(7)-a)&
          +(c(8)-a)+(c(9)-a)+(c(10)-a)+(c(11)-a)+(c(12)-a)+(c(13)-a)
     !b=(u(i-1,j-1)-a)*(u(i,j-1)-a)*(u(i+1,j-1)-a)*(u(i-1,j)-a)
     if(a.ne.0) then
        !u(i,j)=min0(-u(i-1,j-1),-u(i,j-1),-u(i+1,j-1),-u(i-1,j))
        !u(i,j)=-u(i,j)
        !u(i,j)=min0(u(i-1,j-1),u(i,j-1))
        !write(*,*)u(i,j) ,'inin'
        u(i,j,k)=inin*10
        do n=1,13
           if(c(n).ne.0)then
              u(i,j,k)=min0(c(n),u(i,j,k))
              !write(*,*) u(i,j) ,'00'
           end if
        end do
        if(b.ne.0) then
           call tbl(c(1),u(i,j,k))
           ! write(*,*)'in'
        end if
     else
        u(i,j,k)=tag
        tag=tag+1
     end if
  end if
  1001 continue
  !-----bc------

!  sh=1
  do i=1,inin
     if(table(i)==i) then
        write(*,*) i,'table'
     end if
  end do

  do k=1,mesh-2
     do j=1,mesh-2
        do i=1,mesh-2
           u(i,j,k)=table(u(i,j,k))
        end do
     end do
     write(*,*) table(j)
  end do


  open(100,file='/Users/maeda/Desktop/kaiseki/image/rabel3D.DAT',FORM='FORMATTED')
  open(110,file='/Users/maeda/Desktop/kaiseki/image/sa3D.DAT',FORM='FORMATTED')
  do k=1,mesh-2
     do j=1,mesh-2
        do i=1,mesh-2
           write(110,*) u(i,j,k)-v(i,j,k)
           write(100,*) u(i,j,k)
        end do
     end do
     !write(100,*)
     !write(110,*)
  end do
  close(100)
  close(110)

  Deallocate(u,v,rho)
end program main

subroutine tbl(c,min)
  use common
  integer c(13)
  integer aa,min,i,k

  do k=1,13
  aa=(c(k)-min)*c(k)
  if(aa .ne. 0) then
     !do i=1,100
     !   if(c(k)==i) then
     !table(c(k))=min
     table(c(k))=table(min)
     !write(*,*) c(k) , min,'oo'
     !   end if
     !end do
  end if
  end do
end subroutine tbl

