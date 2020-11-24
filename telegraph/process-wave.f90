module combar
!implicit none
integer*8 :: prc=0,prcn=0
integer :: nx=1028,ny=1028,nz=1028,iwx,iwy,iwz
end module combar


program main
    use combar
!    implicit none
    integer Nlx,Nly,Nlz,Ncellx,Ncelly,Ncellz
    integer ncx,ncy,ncz,ngrid,NGL,NGcr,ncycle,NPRE,NPOST
    

!nx=64
!ny=64
!nz=64
!Nlx=8
!Nly=8
!Nlz=8
Nlx=int((1.d5)**(1.d0/3.d0))
Nly=int((1.d5)**(1.d0/3.d0))
Nlz=int((1.d5)**(1.d0/3.d0))
nx=int(1.d4/Nlx)
ny=int(1.d4/Nly)
nz=int(1.d4/Nlz)
!Nlx=int((2.d4)**(1.d0/3.d0))
!Nly=int((2.d4)**(1.d0/3.d0))
!Nlz=int((2.d4)**(1.d0/3.d0))
!nx=int(5.d3/Nlx)
!ny=int(5.d3/Nly)
!nz=int(5.d3/Nlz)

write(*,*) nx,ny,nz,Nlx,nx*Nlx,Nlx*Nly*Nlz



  prc=nx*ny*nz+prc
  !write(*,*) prc
  iwx=1;iwy=0;iwz=0
  call BCgrv()
  call muslcslv1D()
   
   iwx=0;iwy=1;iwz=0
   call BCgrv()
   call muslcslv1D()
   
   iwx=0;iwy=0;iwz=1
   call BCgrv()
   call muslcslv1D()
   
   prc=prc+nx*ny*nz*4
   !iwx=1;iwy=1;iwz=1
   !call BCgrv(100,1,8)
   iwx=0;iwy=0;iwz=1
   call BCgrv()
   call muslcslv1D()
   
   iwx=0;iwy=1;iwz=0
   call BCgrv()
   call muslcslv1D()
   
   iwx=1;iwy=0;iwz=0
   call BCgrv()
   call muslcslv1D()
   
   !%%%%%%%%%%%%%%%%%phi(t+0.5*dt)%%%%%%%%%%%%%%%%%%
!write(*,*) prc

   !%%%%%%%%%%%%%%%%%phigrd(t+0.5*dt)%%%%%%%%%%%%%%%%%%
   iwx=1;iwy=0;iwz=0
   call BCgrv()
   call muslcslv1D()
   
   iwx=0;iwy=1;iwz=0
   call BCgrv()
   call muslcslv1D()
   
   iwx=0;iwy=0;iwz=1
   call BCgrv()
   call muslcslv1D()
   
   !iwx=1;iwy=1;iwz=1
   prc=prc+nx*ny*nz*31
  

   iwx=0;iwy=0;iwz=1
   call BCgrv()
   call muslcslv1D()
   
   iwx=0;iwy=1;iwz=0
   call BCgrv()
   
   iwx=1;iwy=0;iwz=0
   call BCgrv()
   !%%%%%%%%%%%%%%%%%phigrd(t+0.5*dt)%%%%%%%%%%%%%%%%%%


  !%%%%%%%%%%%%%%%%%phi(t+0.5*dt)%%%%%%%%%%%%%%%%%%
  iwx=1;iwy=0;iwz=0
   call BCgrv()
   call muslcslv1D()
   
   iwx=0;iwy=1;iwz=0
   call BCgrv()
   call muslcslv1D()
   
   iwx=0;iwy=0;iwz=1
   call BCgrv()
   call muslcslv1D()
   
   prc=prc+nx*ny*nz*4
   !iwx=1;iwy=1;iwz=1
   !call BCgrv(100,1,8)
   iwx=0;iwy=0;iwz=1
   call BCgrv()
   call muslcslv1D()
   
   iwx=0;iwy=1;iwz=0
   call BCgrv()
   call muslcslv1D()
   
   iwx=1;iwy=0;iwz=0
   call BCgrv()
   call muslcslv1D()
   

write(*,*) prc,prcn,dble(prc+prcn)/dble(prcn+prc/Nlx/Nly/Nlz), &
    1.d0/dble(prc+prcn)*dble(prcn), 1-1.d0/dble(prc+prcn)*dble(prcn) , &  !,1.d0/dble(prc+prcn)*dble(prc/Nlx/Nly/Nlz)
    dble(prc+prcn)/dble(prcn+prc/Nlx/Nly/Nlz)/dble(Nlx*Nly*Nlz)*100
end program main


subroutine BCgrv()
use combar

prc=nx*ny*2*2+prc
if(iwx==1)then
!prcn=nx*ny*2*2+prcn
prcn=nx*ny*2+prcn
endif
end subroutine BCgrv

subroutine muslcslv1D()
  use combar


     call fluxcal()
     !call fluxcal(Phipre,Phipre,Phiu,0.0d0,0.0d0,10)
     !------------calcurate dt/2------------
     prc=prc + 49*nx*ny*nz
     !------------calcurate dt/2------------
     call fluxcal()
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,1)
     !write(*,*) Phiu(127),'127-2'
     !do i = ist , ndx-ien
     prc=prc + 55*nx*ny*nz

  !------------ul.solver.+cg-------------

end subroutine muslcslv1D

subroutine fluxcal()
  use combar
  integer Ncl,Ncm
  Ncl=nx
  Ncm=ny
  
  

              DO Lnum = 1, Ncl
              DO Mnum = 1, Ncm
              call vanalbada()
              prc=prc + 65*nx
              end do
              end DO

end subroutine fluxcal


subroutine vanalbada()
  use combar
  prc=prc+ndx*63
end subroutine vanalbada
