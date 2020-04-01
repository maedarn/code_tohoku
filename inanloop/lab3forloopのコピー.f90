module common
  implicit none
  integer,parameter :: nx=512,ny=512,nz=512
  integer,parameter :: rank=nx*ny*nz/4
  integer,dimension(0:rank) :: rl_table,next_table,tail_table,T
  !integer,allocatable,dimension(:,:,:) :: tag
end module common

program main
  use common
  implicit none
  integer,allocatable,dimension(:,:,:) :: tag,den
  !integer,allocatable,dimension(:,:,:) :: den
  real(4),allocatable,dimension(:,:,:) :: rho
  integer lab,l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13
  integer l14,l15,l16,l17,l18,l19,l20,l21,l22
  integer i,j,k,Vb,r,m,mm,in,loop,ii,tmstp,timeloop,timestep,cnt10!,nsisu
  double precision nsisu
  integer itime,ftime
  real(4),allocatable,dimension(:) ::rhoshdlp
  integer,allocatable,dimension(:) :: rgn
  real(4) rhoshd,dmy1
  character(3) num,timech
  character(46):: dir='/glv0/maedarn/clst-form-HIcol/cnv100-wb-wg-sm/'

  !allocate(tag(1:nx,1:ny,1:nz))
  !allocate(den(1:nx,1:ny,1:nz))
  !allocate(rho(1:nx,1:ny,1:nz))

  !allocate(tag(0:nx,0:ny,0:nz))
  !allocate(den(0:nx,0:ny,0:nz))
  !allocate(rho(0:nx,0:ny,0:nz))]
  loop=1
  !timeloop=50
  !timestep=1
  timestep=1
  itime=98
  ftime=98

  allocate(tag(0:nx+1,0:ny+1,0:nz+1))
  allocate(den(0:nx+1,0:ny+1,0:nz+1))
  allocate(rho(0:nx+1,0:ny+1,0:nz+1))
  allocate(rhoshdlp(loop))
  allocate(rgn(loop))

  nsisu=1.5d0
  !nsisu=2.d0
  do j=1,loop
     nsisu=1.5d0
     nsisu =nsisu + 0.5d0*dble(j)
     rhoshdlp(j) = sngl(10.0d0**nsisu)
     !rhoshdlp(1)=50.e0
     !rhoshdlp(2)=500.e0
     !rhoshdlp(3)=1000.e0
     !rhoshdlp(4)=5000.e0
     !rhoshdlp(5)=10000.e0
     write(*,*) rhoshdlp(j)
  end do

  rhoshdlp(1) = 1.e4

  !do tmstp = 1,timeloop
   do tmstp = itime,ftime,timestep
   !  open(unit=350,file='cnt2.dat')
   !  read(350,*) cnt10
   !  close(350)
   !  write(timech,'(i3.3)') cnt10
   !  cnt10=cnt10-timestep !minusstep
   !  open(unit=350,file='cnt2.dat')
   !  write(350,*) cnt10
      !  close(350)
      write(timech,'(i3.3)') tmstp
!  lab=1
  Vb=0
  !rhoshd=5000.e0
  !rhoshd=23700.e0
!  tag(:,:,:)=0
  rho(:,:,:)=0.e0
  rgn(:)=0
!  rl_table(:)=0
!  next_table(:)=0
!  tail_table(:)=0
!  T(:)=0

!  goto 1000
  !open(150,file='Allnewbig'//timech//'.DAT',FORM='UNFORMATTED')
  open(150,file=dir//'All/All'//timech//'.DAT',access='stream',FORM='UNFORMATTED')
  do k=1,nz
     do j=1,ny
        do i=1,nx
           read(150) rho(i,j,k),dmy1,dmy1,dmy1,dmy1,dmy1,dmy1,dmy1,dmy1
        end do
     end do
     write(*,*) rho(256,256,k),k,timech
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
  !do k=1,nz
  !   do j=1,ny
   !      do i=1,nx
   do ii=1,loop
      Vb=0
      tag(:,:,:)=0
      den(:,:,:)=0
      rl_table(:)=0
      next_table(:)=0
      tail_table(:)=0
      T(:)=0
      lab=1
      write(*,*) 'ii',ii,rhoshdlp(ii)
  do k=0,nz+1
     do j=0,ny+1
        do i=0,nx+1
           if(rho(i,j,k)>rhoshdlp(ii)) then
              !rho(i,j,k)=1.e0
              den(i,j,k)=1
           else
              !rho(i,j,k)=0.d0
              den(i,j,k)=0
           end if
        end do
     end do
     write(*,*) k,rhoshdlp(ii),ii,timech
  end do


!  do m=1,2
  do k=1,nz
  do j=1,ny
  do i=1,nx

  if(k==nz) then
     tag(:,:,nz+1)=tag(:,:,1)
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

        l14=tag(i-1,j-1,k+1)
        l15=tag(i  ,j-1,k+1)
        l16=tag(i+1,j-1,k+1)
        l17=tag(i-1,j  ,k+1)
        l18=tag(i  ,j  ,k+1)
        l19=tag(i+1,j  ,k+1)
        l20=tag(i-1,j+1,k+1)
        l21=tag(i  ,j+1,k+1)
        l22=tag(i+1,j+1,k+1)
        call bcz(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,&
             l14,l15,l16,l17,l18,l19,l20,l21,l22,lab,tag(i,j,k))
     else
        tag(i,j,k)=Vb
     endif
     goto 400
  end if
  if(j==ny) then
     tag(:,ny+1,:)=tag(:,1,:)
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

        l14=tag(i-1,j+1,k  )
        l15=tag(i  ,j+1,k  )
        l16=tag(i+1,j+1,k  )
        call bcy(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,&
             l14,l15,l16,lab,tag(i,j,k))
     else
        tag(i,j,k)=Vb
     endif
     goto 400
  end if
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
 400 continue
 enddo
 tag(:, 0,:)=tag(:,ny,:)
 write(*,*) i,j,k
 enddo
 tag(:,:, 0)=tag(:,:,nz)
 !write(*,*) k
 enddo

 !tag(:,:, 0)=tag(:,:,nz)
 !tag(:,:,nz+1)=tag(:,:,1)
 !tag(:, 0,:)=tag(:,ny,:)
 !tag(:,ny+1,:)=tag(:,1,:)
 !enddo


 !変換
 T(:)=rl_table(:)
 T(Vb)=0
 mm=0
 !do i=1,rank
 do i=1,7000
    !if(rl_table(i)==i) then
    in=1
    !do j=1,rank
    do j=1,7000
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
 rgn(ii)=mm
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

!  open(290,file='/Users/maeda/Desktop/kaiseki/image/kakuninrl.DAT',FORM='FORMATTED')
!  do k=1,3000
!     write(290,*) k,rl_table(k),T(k)
!  end do
!  close(290)
  write(num,'(i3.3)') ii
  !open(110,file='taglog'//timech//num//'.DAT',FORM='FORMATTED')
   open(110,file=dir//'tag/tag'//timech//num//'.DAT',FORM='UNFORMATTED')
  do k=1,nz
     do j=1,ny
        do i=1,nx
           !write(110,*) u(i,j,k)-v(i,j,k)
           write(110) tag(i,j,k)!,den(i,j,k),rho(i,j,k)
        end do
     end do
     write(*,*) k,'savetag'
  end do
  close(110)

  !open(111,file='taglogkakunin'//timech//num//'.DAT',FORM='FORMATTED')
!  do k=1,nz,4
!     do j=1,ny,4
!        do i=1,nx,4
           !write(110,*) u(i,j,k)-v(i,j,k)
!           write(111,*) tag(i,j,k),den(i,j,k),rho(i,j,k)
!        end do
!     end do
!     write(*,*) k,'savetag'
!  end do
!  close(111)
enddo

!open(120,file='rgn'//timech//'.DAT',FORM='FORMATTED')
open(120,file=dir//'tag/rgn-low'//timech//'.DAT',FORM='FORMATTED')
  do i=1,loop
     write(120,*) i,rhoshdlp(i),rgn(i)
  end do
  close(120)

end do
  deallocate(tag,rho,den,rgn,rhoshdlp)
end program main

!subroutine resolve(u,v)
subroutine resolve(a,b)
  !implicit none
  use common
  implicit none
  integer :: u,v,i,ik,count=0,a,b
  character(3) chr
  !write(chr,'(I3.3)') count
  !count=count+1
  !open(190,file='/Users/maeda/Desktop/kaiseki/image/kakunin'//chr//'.DAT',FORM='FORMATTED')
  !do ik=1,150
  !   write(190,*) ik,rl_table(ik),next_table(ik),tail_table(ik),a,b
  !end do
  !close(190)
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

subroutine bcz(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,&
     l14,l15,l16,l17,l18,l19,l20,l21,l22,lab,tg)
  use common
  integer l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,&
       l14,l15,l16,l17,l18,l19,l20,l21,l22
  integer,dimension(22) :: l
  integer m,i,j,k,lab,tg

  l(1 )=l1
  l(2 )=l2
  l(3 )=l3
  l(4 )=l4
  l(5 )=l5
  l(6 )=l6
  l(7 )=l7
  l(8 )=l8
  l(9 )=l9
  l(10)=l10
  l(11)=l11
  l(12)=l12
  l(13)=l13
  l(14)=l14
  l(15)=l15
  l(16)=l16
  l(17)=l17
  l(18)=l18
  l(19)=l19
  l(20)=l20
  l(21)=l21
  l(22)=l22

  m=0
  do i=1,22
     m=max0(l(i),m)
  end do
  if(m.ne.0)then
     m=10000
     do i=1,22
        if(l(i).ne.0) then
           m=min0(l(i),m)
        end if
     end do
     tg=m
     do i=1,21
        k=i+1
        do j=k,22
           if((l(i)>0).and.(l(j)>0).and.(l(i).ne.l(j)))then
              call resolve(l(i),l(j))
           end if
        end do
     end do
  else
     tg=lab
     !T(lab)=lab
     rl_table(lab)=lab
     next_table(lab)= -1
     tail_table(lab)= lab
     lab=lab+1
  end if
end subroutine bcz

subroutine bcy(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,&
     l14,l15,l16,lab,tg)
  use common
  integer l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,&
       l14,l15,l16
  integer,dimension(22) :: l
  integer m,i,j,k,lab,tg

  l(1 )=l1
  l(2 )=l2
  l(3 )=l3
  l(4 )=l4
  l(5 )=l5
  l(6 )=l6
  l(7 )=l7
  l(8 )=l8
  l(9 )=l9
  l(10)=l10
  l(11)=l11
  l(12)=l12
  l(13)=l13
  l(14)=l14
  l(15)=l15
  l(16)=l16

  m=0
  do i=1,16
     m=max0(l(i),m)
  end do
  if(m.ne.0)then
     m=10000
     do i=1,16
        if(l(i).ne.0) then
           m=min0(l(i),m)
        end if
     end do
     tg=m
     do i=1,16-1
        k=i+1
        do j=k,16
           if((l(i)>0).and.((j)>0).and.(l(i).ne.l(j)))then
              call resolve(l(i),l(j))
           end if
        end do
     end do
  else
     tg=lab
     !T(lab)=lab
     rl_table(lab)=lab
     next_table(lab)= -1
     tail_table(lab)= lab
     lab=lab+1
  end if
end subroutine bcy
