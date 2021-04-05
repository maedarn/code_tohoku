SUBROUTINE feedback(dt,mode)
USE comvar
USE mpivar
USE slfgrv
INCLUDE 'mpif.h'
DOUBLE PRECISION  :: dt,dxi,intt=0.d0
DOUBLE PRECISION  :: LSFE=0.02d0,mstarmn=1.d0,pstar
INTEGER :: LEFTt,RIGTt,TOPt,BOTMt,UPt,DOWNt,jtime=0
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
DOUBLE PRECISION  :: VECU
character(3) Nfinal,itime
integer :: i,j,k
double precision :: div(-1:ndx,-1:ndy,-1:ndz,1:4)
double precision :: rhoth=5.d0*1.d5*1.27d0,SFE,SFratio=0.5d0,nid=0.d0
DOUBLE PRECISION ran
double precision :: Ustar(1:9)


if(mode==1) then

intt=intt+dt

do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
div(i,j,k,1)=(U(i-1,j,k,2)-U(i+1,j,k,2))*0.5d0/dx1+(U(i,j-1,k,3)-U(i,j+1,k,3))*0.5d0/dy1+(U(i,j,k-1,4)-U(i,j,k+1,4))*0.5d0/dz1
div(i,j,k,2)=U(i-1,j,k,2)*U(i+1,j,k,2)
div(i,j,k,3)=U(i,j-1,k,3)*U(i,j+1,k,3)
div(i,j,k,4)=U(i,j,k-1,4)*U(i,j,k+1,4)
end do; end do; end do

do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
if((div(i,j,k,1)<0.d0).and.(div(i,j,k,2)<0.d0).and.(div(i,j,k,3)<0.d0).and.(div(i,j,k,4)<0.d0).and.(U(i,j,k,1)>rhoth)) then
SFE=LSFE*dsqrt(U(i,j,k,1)/(4.04d0*1.d3))*dt
SFE=LSFE*U(i,j,k,1)/dsqrt(G*U(i,j,k,1))*dx1*dy1*dz1*dt
pstar=(1-dexp(-SFE/mstarmn))
call ran0(ran,1)
if(ran<pstar) then
nid=nid+1.d0
Ustar(1)=i
Ustar(2)=j
Ustar(3)=k
Ustar(4)=U(i,j,k,2)
Ustar(5)=U(i,j,k,3)
Ustar(6)=U(i,j,k,4)
Ustar(7)=U(i,j,k,1)*dx1*dy1*dz1*SFratio
Ustar(8)=intt
Ustar(9)=nid
endif
endif
end do; end do; end do

endif




end SUBROUTINE feedback
