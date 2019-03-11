module common
  implicit none
  integer table(100)
end module common

program main
  !implicit none
  use common
  implicit none
  !real(4) u(100,100)
  integer v(100,100),u(100,100)!,table(100,100)
  integer i,j,k,tag,a,b,l,m,n,c(4),r1,r2,sh,sh2(100)

  table(:)=0
  !do j=1,100
  do i=1,100
     table(i)=i
  end do
  !end do

  v(:,:)=0
  u(:,:)=0
  !initial
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
           v(i,j)=2
        end if
        if(r1==30 .or. r1==31) then
           v(i,j)=1
        end if
     end do
  end do
  open(10,file='/Users/maeda/Desktop/kaiseki/image/int.DAT',FORM='FORMATTED')
  do j=1,100
     do i=1,100
        write(10,*) v(i,j)
     end do
     write(10,*)
  end do
  close(10)
  !initial


  tag=1
  do j=2,100
     do i=2,100
        if(v(i,j) .ne. 0) then
           c(1)=u(i-1,j-1)
           c(2)=u(i,j-1)
           c(3)=u(i+1,j-1)
           c(4)=u(i-1,j)
           a=u(i-1,j-1)+u(i,j-1)+u(i+1,j-1)+u(i-1,j)
           b=(u(i-1,j-1)-a)*(u(i,j-1)-a)*(u(i+1,j-1)-a)*(u(i-1,j)-a)
           if(a.ne.0) then
              !u(i,j)=min0(-u(i-1,j-1),-u(i,j-1),-u(i+1,j-1),-u(i-1,j))
              !u(i,j)=-u(i,j)
              !u(i,j)=min0(u(i-1,j-1),u(i,j-1))
              !write(*,*)u(i,j) ,'inin'
              u(i,j)=1000
              do n=1,4
                 if(c(n).ne.0)then
                    u(i,j)=min0(c(n),u(i,j))
                    !write(*,*) u(i,j) ,'00'
                 end if
              end do
              if(b.ne.0) then
                 call tbl(u(i-1,j-1),u(i,j-1),u(i+1,j-1),u(i-1,j),u(i,j))
                ! write(*,*)'in'
              end if
           else
              u(i,j)=tag
              tag=tag+1
           end if
        end if
     end do
  end do

!  sh=1
  do i=1,100
     if(table(i)==i) then
        write(*,*) i,'table'
     end if
  end do

  do j=1,100
     do i=1,100
        u(i,j)=table(u(i,j))
     end do
     write(*,*) table(j)
  end do


  open(100,file='/Users/maeda/Desktop/kaiseki/image/rabel.DAT',FORM='FORMATTED')
  open(110,file='/Users/maeda/Desktop/kaiseki/image/sa.DAT',FORM='FORMATTED')
  do j=1,100
     do i=1,100
        write(110,*) u(i,j)-v(i,j)
        write(100,*) u(i,j)
     end do
     write(100,*)
     write(110,*)
  end do
  close(100)
  close(110)

end program main

subroutine tbl(c1,c2,c3,c4,min)
  use common
  integer c(4)
  integer aa,min,i,k,c1,c2,c3,c4
  c(1)=c1
  c(2)=c2
  c(3)=c3
  c(4)=c4

  do k=1,4
  aa=(c(k)-min)*c(k)
  if(aa .ne. 0) then
     !do i=1,100
     !   if(c(k)==i) then
     !table(c(k))=min
     table(c(k))=table(min)
     write(*,*) c(k) , min,'oo'
     !   end if
     !end do
  end if
  end do
end subroutine tbl

