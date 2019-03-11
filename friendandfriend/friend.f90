program computeraw
  implicit none
  real(4), allocatable :: u(:,:,:,:,:)
  integer i,j,k,mesh,tag,val

  !----parameter----
  rhoshd
  !----parameter----


  ALLOCATE(u(1:mesh,1:mesh,1:mesh,1:val,1))

  open(unit=150,file='/Users/maeda/Desktop/kaiseki/cnv100wbwg/135.dat',FORM='UNFORMATTED')

  do k = 1, Ncellz
     do j = 1, Ncelly
        do i= 1, Ncellx
           read(150) u(i,j,k,1)
        end do
     end do
  end do


  do k = 1, Ncellz
     do j = 1, Ncelly
        do i= 1, Ncellx
           if(u(i,j,k,1)<rhoshd) then
              u(i,j,k,1)=0.e0
           else
              u(i,j,k,1)=1.e0
           end if
        end do
     end do
  end do

  do k = 1, Ncellz
     do j = 1, Ncelly
        do i= 1, Ncellx
           if(u(i,j,k,1,tag)>rhoshd) then
              if
           end if
        end do
     end do
  end do
  DEALLOCATE(u)
end program computeraw
