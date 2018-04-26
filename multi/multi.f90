SUBROUTINE mglin(u,n,ncycle)
INTEGER n,ncycle,NPRE,NPOST,NG,MEMLEN
DOUBLE PRECISION u(n,n)
PARAMETER (NG=5,MEMLEN=13*2**(2*NG)/3+14*2**NG+8*NG-100/3)
PARAMETER (NPRE=1,NPOST=1)
INTEGER j,jcycle,jj,jpost,jpre,mem,nf,ngrid,nn,ires(NG), irho(NG),irhs(NG),iu(NG),maloc
DOUBLE PRECISION z
COMMON /memory/ z(MEMLEN),mem  !Storage for grid functions is allocated by maloc
mem=0 !from array z.
nn=n/2+1     !nはx,yのメッシュ数,nnは半分の地点
ngrid=NG-1 !レベル2**Nみたいなもの（何回イタレーションするか）
irho(ngrid)=maloc(nn**2)   ! Allocate storage for r.h.s. on grid NG − 1,
call rstrct(z(irho(ngrid)),u,nn)   !and ll it by restricting from the ne grid.
1 if (nn.gt.3) then   !Similarly allocate storage and ll r.h.s. on all
nn=nn/2+1 !coarse grids.
ngrid=ngrid-1
irho(ngrid)=maloc(nn**2)
call rstrct(z(irho(ngrid)),z(irho(ngrid+1)),nn)
goto 1
endif
nn=3
iu(1)=maloc(nn**2)
irhs(1)=maloc(nn**2)
call slvsml(z(iu(1)),z(irho(1))) !Initial solution on coarsest grid.
ngrid=NG
do j=2,ngrid! Nested iteration loop.
nn=2*nn-1
iu(j)=maloc(nn**2)
irhs(j)=maloc(nn**2)
ires(j)=maloc(nn**2)
call interp(z(iu(j)),z(iu(j-1)),nn) !Interpolate from coarse grid to next ner grid.

if (j.ne.ngrid) then
call copy(z(irhs(j)),z(irho(j)),nn) !Set up r.h.s.
else
call copy(z(irhs(j)),u,nn)
endif

do jcycle=1,ncycle !V-cycle loop.
nf=nn
do jj=j,2,-1 !Downward stoke of the V.
   do jpre=1,NPRE !Pre-smoothing.
call relax(z(iu(jj)),z(irhs(jj)),nf)
enddo

call resid(z(ires(jj)),z(iu(jj)),z(irhs(jj)),nf)
nf=nf/2+1
call rstrct(z(irhs(jj-1)),z(ires(jj)),nf)

!Restriction of the residual is the next r.h.s.

call fill0(z(iu(jj-1)),nf) !Zero for initial guess in next relaxation.
enddo
call slvsml(z(iu(1)),z(irhs(1))) !Bottom of V: solve on coarsest grid.
nf=3
do jj=2,j !Upward stroke of V.
nf=2*nf-1
call addint(z(iu(jj)),z(iu(jj-1)),z(ires(jj)),nf)

!Use res for temporary storage inside addint.

do jpost=1,NPOST !Post-smoothing.
call relax(z(iu(jj)),z(irhs(jj)),nf)
enddo
enddo
enddo
enddo
call copy(u,z(iu(ngrid)),n) !Return solution in u.
return
END SUBROUTINE mglin

SUBROUTINE rstrct(uc,uf,nc)
INTEGER nc
DOUBLE PRECISION uc(nc,nc),uf(2*nc-1,2*nc-1)
!Half-weighting restriction. nc is the coarse-grid dimension. The ne-grid solution is input
!in uf(1:2*nc-1,1:2*nc-1), the coarse-grid solution is returned in uc(1:nc,1:nc).
INTEGER ic,iff,jc,jf
do  jc=2,nc-1 !Interior points.
jf=2*jc-1
do  ic=2,nc-1
iff=2*ic-1
uc(ic,jc)=0.5d0*uf(iff,jf)+0.125d0*(uf(iff+1,jf)+ uf(iff-1,jf)+uf(iff,jf+1)+uf(iff,jf-1)) !平均
enddo
enddo
do  ic=1,nc !Boundary points.
uc(ic,1)=uf(2*ic-1,1)
uc(ic,nc)=uf(2*ic-1,2*nc-1)
enddo
do  jc=1,nc !Boundary points.
uc(1,jc)=uf(1,2*jc-1)
uc(nc,jc)=uf(2*nc-1,2*jc-1)
enddo
return
END SUBROUTINE rstrct

SUBROUTINE interp(uf,uc,nf)
INTEGER nf
DOUBLE PRECISION uc(nf/2+1,nf/2+1),uf(nf,nf)
INTEGER ic,iff,jc,jf,nc
!Coarse-to-ne prolongation by bilinear interpolation. nf is the ne-grid dimension. The
!coarse-grid solution is input as uc(1:nc,1:nc), where nc = nf=2 + 1. The ne-grid
!solution is returned in uf(1:nf,1:nf).
nc=nf/2+1
do  jc=1,nc !Do elements that are copies. 増やしたグリッドに元の値（雑な）を代入
jf=2*jc-1
do  ic=1,nc
uf(2*ic-1,jf)=uc(ic,jc)
enddo
enddo

do jf=1,nf,2 !Do odd-numbered columns, interpolating verdo 内挿
do iff=2,nf-1,2 !tically.
uf(iff,jf)=0.5d0*(uf(if+1,jf)+uf(iff-1,jf))
enddo
enddo
do jf=2,nf-1,2 !Do even-numbered columns, interpolating hordo
do iff=1,nf !izontally.
uf(iff,jf)=0.5d0*(uf(iff,jf+1)+uf(iff,jf-1))
enddo
enddo
return
END SUBROUTINE interp

SUBROUTINE addint(uf,uc,res,nf)
INTEGER nf
DOUBLE PRECISION res(nf,nf),uc(nf/2+1,nf/2+1),uf(nf,nf)
!C USES interp
!Does coarse-to-ne interpolation and adds result to uf. nf is the ne-grid dimension. The
!coarse-grid solution is input as uc(1:nc,1:nc), where nc = nf=2 + 1. The ne-grid
!solution is returned in uf(1:nf,1:nf). res(1:nf,1:nf) is used for temporary storage.
INTEGER i,j
call interp(res,uc,nf)
do 12 j=1,nf
do 11 i=1,nf
uf(i,j)=uf(i,j)+res(i,j)
enddo 11
enddo 12
return
END SUBROUTINE addint

SUBROUTINE slvsml(u,rhs)
DOUBLE PRECISION rhs(3,3),u(3,3)
!C USES fill0
!Solution of the model problem on the coarsest grid, where h = 1
!2 . The right-hand side is
!input in rhs(1:3,1:3) and the solution is returned in u(1:3,1:3).
DOUBLE PRECISION h
call fill0(u,3)
h=0.5d0
u(2,2)= -h*h*rhs(2,2)/4.d0
return
END SUBROUTINE slvsml

SUBROUTINE relax(u,rhs,n) !ガウスサイデル
INTEGER n
DOUBLE PRECISION rhs(n,n),u(n,n)
!Red-black Gauss-Seidel relaxation for model problem. The current value of the solution
!u(1:n,1:n) is updated, using the right-hand side function rhs(1:n,1:n).
INTEGER i,ipass,isw,j,jsw
DOUBLE PRECISION h,h2
h=1.d0/(n-1) !メッシュ間隔
h2=h*h
jsw=1
do ipass=1,2 ! Red and black sweeps.
isw=jsw
do j=2,n-1
do i=isw+1,n-1,2 !Gauss-Seidel formula.
u(i,j) = 0.25d0*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-h2*rhs(i,j))
enddo
isw=3-isw
enddo
jsw=3-jsw
enddo
return
END SUBROUTINE relax

SUBROUTINE resid(res,u,rhs,n)
INTEGER n
DOUBLE PRECISION res(n,n),rhs(n,n),u(n,n)
!Returns minus the residual for the model problem. Input quantities are u(1:n,1:n) and
!rhs(1:n,1:n), while res(1:n,1:n) is returned.
INTEGER i,j
DOUBLE PRECISION h,h2i
h=1.d0/(n-1)
h2i=1.d0/(h*h)
do j=2,n-1 !Interior points.
do i=2,n-1
res(i,j)=-h2i*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4.d0*u(i,j))+rhs(i,j)
enddo
enddo
do i=1,n !Boundary points.
res(i,1)=0.d0
res(i,n)=0.d0
res(1,i)=0.d0
res(n,i)=0.d0
enddo
return
END SUBROUTINE resid

SUBROUTINE copy(aout,ain,n)
INTEGER n
DOUBLE PRECISION ain(n,n),aout(n,n)
!Copies ain(1:n,1:n) to aout(1:n,1:n).
INTEGER i,j
do i=1,n
do j=1,n
aout(j,i)=ain(j,i)
enddo
enddo
return
END SUBROUTINE copy

SUBROUTINE fill0(u,n)
INTEGER n
DOUBLE PRECISION u(n,n)
!Fills u(1:n,1:n) with zeros.
INTEGER i,j
do j=1,n
do i=1,n
u(i,j)=0.d0
enddo
enddo
return
END SUBROUTINE fill0

FUNCTION maloc(len) !len=格子の個数
INTEGER maloc,len,NG,MEMLEN
PARAMETER (NG=5,MEMLEN=13*2**(2*NG)/3+14*2**NG+8*NG-100/3) !for mglin
PARAMETER (NG=5,MEMLEN=17*2**(2*NG)/3+18*2**NG+10*NG-86/3) !for mgfas, N.B.!
INTEGER mem
DOUBLE PRECISION z
COMMON /memory/ z(MEMLEN),mem
!Dynamical storage allocation. Returns integer pointer to the starting position for len array
!elements in the array z. The preceding array element is lled with the value of len, and
!the variable mem is updated to point to the last element of z that has been used.
if (mem+len+1.gt.MEMLEN) pause 'insufficient memory in maloc'
z(mem+1)=len
maloc=mem+2
mem=mem+len+1
return
END FUNCTION maloc
