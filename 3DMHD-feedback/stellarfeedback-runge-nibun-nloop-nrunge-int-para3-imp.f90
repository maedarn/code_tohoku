SUBROUTINE feedback(dt,mode)
USE comvar
USE mpivar
USE slfgrv
USE chmvar
USE fedvar
INCLUDE 'mpif.h'
DOUBLE PRECISION  ::dt,dxi,intt=0.d0,ranout,Rst,rloop,Rst11,Rst12,Rstmin,Rstmax,Rsv,raninp,ransfe,denst,rdenst,drmv
DOUBLE PRECISION  :: emag1 ,ekin1, eth1, egr1, engy, grlocp
DOUBLE PRECISION  :: stx,sty,stz,volratio,meandrho
DOUBLE PRECISION  :: Rnbnmn,Rnbnmx,pred1
DOUBLE PRECISION  :: Q,lQ,denst1
DOUBLE PRECISION  :: mstarmn=1.d0,pstar,Mstar1,Mstarlp,T0=1.d3,pi=3.14159265359d0
INTEGER :: LEFTt,RIGTt,TOPt,BOTMt,UPt,DOWNt,jtime=0,nstdmy,nmsv1,nidpre,nidpost!,nid=0
!INTEGER, parameter :: nstar=1000000
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
DOUBLE PRECISION  :: VECU
character(3) Nfinal,itime
integer :: ISTini,JSTini,KSTini,idenref1,idenref2,nummsv
integer :: ist1,jst1,kst1,ist2,jst2,kst2,ist3,jst3,kst3,Iist1,Jjst1,Kkst1,Iist2,Jjst2,Kkst2,Iist3,Jjst3,Kkst3,nist1,njst1,nkst1
integer :: i,j,k,itrrst=10,rixmn,rjymn,rkzmn,rixmx,rjymx,rkzmx,nrixmn,nrjymn,nrkzmn,nrixmx,nrjymx,nrkzmx
integer :: rsix,rsjy,rskz,nstloop,NRANKdm,nstdm,lran1=0,lran2=0,idenref,nidck
integer :: nlpnm1,nlpnm6,nlpnm4,nmIST,nmJST,nmKST,nlpnm7,nlpnm3,nlpnm,nmIST1,nmJST1,nmKST1
DOUBLE PRECISION :: nlpnm2,nlpnm5,nlpnm8,nlpnm9,nlpnm10,nlpnm11,nwnid2,nlpnm21,nlpnm22,nlpnm23,nlpnm25,nlpnm26,nlpnm27,nlpnm88,nlpnm888
DOUBLE PRECISION  :: rsixc,rsjyc,rskzc,ndcore,nlpnmpre14,nlpnmpre16
double precision :: div(-1:ndx,-1:ndy,-1:ndz,1:4)
double precision :: SFE,SFratio=0.5d0,Mmassivetot=0.d0,rdumy,pre16,pre14
DOUBLE PRECISION ran
DOUBLE PRECISION :: rdumy1,rdumy2,rdumy3,nmeanrd,nmcellrad,nmeanrdden,nupdt,supdt,nwnid1=0.d0,submass1,submass2
DOUBLE PRECISION :: nrdumy1,nrdumy2,nrdumy3
DOUBLE PRECISION :: nmeanrdmax,nmcellradmax,nmeanrddenmax,nmeanrdmin,nmcellradmin,nmeanrddenmin,nrmv1,rdenst1,radint1,radint2,rhoimp,rhopre1
!double precision :: dnMPI(0:numstar,0:NPE-1),rdnumMPI(0:numstar,0:NPE-1)
!double precision :: dnMPI(0:NPE-1),rdnumMPI(0:NPE-1),rdcellMPI(0:NPE-1)
double precision :: dnMPI(0:NPE-1,3),rdnumMPI(0:NPE-1,3),rdcellMPI(0:NPE-1,3),nrmv(0:NPE-1),nwnid(0:NPE-1),nwnidex(0:NPE-1),stprc(0:NPE-1,3),submas(0:NPE-1,2)
!double precision :: Ustar(1:10,1:nstar)
double precision :: dmsiv(0:NPE-1),dmsiv_gt(0:NPE-1),submas1(0:NPE,2)
!DOUBLE PRECISION, dimension(:,:), allocatable :: submas1
double precision :: vv1,vv2,vv3,xx1,xx2,xx3,ffgas=1.d0,clsmass=10.d0,Rstdmy,masstoflux,SFE2
DOUBLE PRECISION, dimension(:,:)  , allocatable :: STF1
DOUBLE PRECISION, dimension(:,:,:)  , allocatable :: STF2


!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
if(mode==1) then

!Ustar(9,0)=0.d0
Nstinp(:)=0.d0
nidnw=0
nidpre=nid
nwnid1=0.d0
intt=intt+dt
!N_MPI(20)=5; N_MPI(1)=1; N_MPI(2)=2; N_MPI(3)=3; N_MPI(4)=4; N_MPI(5)=5; iwx = 1; iwy = 1; iwz = 1; CALL BC_MPI(2,1)
N_MPI(20)=2; N_MPI(1)=1; N_MPI(2)=5; iwx = 1; iwy = 1; iwz = 1; CALL BC_MPI(2,1)

!do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
!div(i,j,k,1)=(U(i+1,j,k,2)-U(i-1,j,k,2))*0.5d0/dx1+(U(i,j+1,k,3)-U(i,j-1,k,3))*0.5d0/dy1+(U(i,j,k+1,4)-U(i,j,k-1,4))*0.5d0/dz1
!div(i,j,k,2)=(U(i+1,j,k,2)-U(i-1,j,k,2))*0.5d0/dx1
!div(i,j,k,3)=(U(i,j+1,k,3)-U(i,j-1,k,3))*0.5d0/dy1
!div(i,j,k,4)=(U(i,j,k+1,4)-U(i,j,k-1,4))*0.5d0/dz1
!end do; end do; end do
MGtoS(NRANK)=0.d0
Fstar(NRANK)=0.d0
do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
!nlpnm2=(0.5d0-dsign(0.5d0,div(i,j,k,1)))*(0.5d0-dsign(0.5d0,div(i,j,k,2)))*(0.5d0-dsign(0.5d0,div(i,j,k,3)))*(0.5d0-dsign(0.5d0,div(i,j,k,4)))&
!*(0.5d0-dsign(0.5d0,-U(i,j,k,1)+rhoth))
!nlpnm2=(0.5d0-dsign(0.5d0,div(i,j,k,1)))*(0.5d0-dsign(0.5d0,div(i,j,k,2)))*(0.5d0-dsign(0.5d0,div(i,j,k,3)))*(0.5d0-dsign(0.5d0,div(i,j,k,4)))&
!*(0.5d0-dsign(0.5d0,-ndH2(i,j,k)+rhoth))

gammi1 =   3.d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))+5.d0*ndH2(i,j,k)
gammi1 = ( 2.d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))+2.d0*ndH2(i,j,k) )/gammi1

emag1 = 0.5d0 * ( U(i,j,k,6)**2 + U(i,j,k,7)**2 + U(i,j,k,8)**2 )
ekin1 = 0.5d0 * ( U(i,j,k,2)**2 + U(i,j,k,3)**2 + U(i,j,k,4)**2 ) * U(i,j,k,1)
eth1  = U(i,j,k,5)/(gammi1)
egr1  =Phi(i  ,j  ,k  )-&
       (Phi(i+1,j+1,k+1)+Phi(i+1,j+1,k  )+Phi(i+1,j+1,k-1)+ &
        Phi(i+1,j  ,k+1)+Phi(i+1,j  ,k  )+Phi(i+1,j  ,k-1)+ &
        Phi(i+1,j-1,k+1)+Phi(i+1,j-1,k  )+Phi(i+1,j-1,k-1)+ &
        Phi(i  ,j+1,k+1)+Phi(i  ,j+1,k  )+Phi(i  ,j+1,k-1)+ &
        Phi(i  ,j  ,k+1)                 +Phi(i  ,j  ,k-1)+ &
        Phi(i  ,j-1,k+1)+Phi(i  ,j-1,k  )+Phi(i  ,j-1,k-1)+ &
        Phi(i-1,j+1,k+1)+Phi(i-1,j+1,k  )+Phi(i-1,j+1,k-1)+ &
        Phi(i-1,j  ,k+1)+Phi(i-1,j  ,k  )+Phi(i-1,j  ,k-1)+ &
        Phi(i-1,j-1,k+1)+Phi(i-1,j-1,k  )+Phi(i-1,j-1,k-1)   )/26.0

engy=emag1+ekin1+eth1+U(i,j,k,1)*egr1
grlocp=(0.5d0-dsign(0.5d0,Phi(i,j,k)-Phi(i+1,j,k)))*(0.5d0-dsign(0.5d0,Phi(i,j,k)-Phi(i,j+1,k)))*(0.5d0-dsign(0.5d0,Phi(i,j,k)-Phi(i,j,k+1)))&
+(0.5d0-dsign(0.5d0,Phi(i,j,k)-Phi(i-1,j,k)))*(0.5d0-dsign(0.5d0,Phi(i,j,k)-Phi(i,j-1,k)))*(0.5d0-dsign(0.5d0,Phi(i,j,k)-Phi(i,j,k-1)))-0.5d0

nlpnm2=(0.5d0-dsign(0.5d0,-ndH2(i,j,k)+rhoth))*(0.5d0+dsign(0.5d0,-U(i,j,k,5)/( kb*(ndp(i,j,k)+ndH(i,j,k)+ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)))+Tth/T0)) &
    *(0.5d0-dsign(0.5d0,engy))*(0.5d0+dsign(0.5d0,grlocp))
!U(i,j,k,5)/( kb*(ndp(i,j,k)+ndH(i,j,k)+ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)) )
!write(*,*)NRANK,nlpnm2,-U(i,j,k,1)+rhoth,div(i,j,k,1),div(i,j,k,2),div(i,j,k,3),div(i,j,k,4),'div'
!nlpnm2=(0.5d0-dsign(0.5d0,-ndH2(i,j,k)+rhoth))*(0.5d0+dsign(0.5d0,-U(i,j,k,5)/( kb*(ndp(i,j,k)+ndH(i,j,k)+ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)))+Tth/T0))

!masstoflux=dsqrt(G*3.1415926535d0)*U(i,j,k,1)*dx1/dsqrt(U(i,j,k,6)**2.0+U(i,j,k,7)**2.0+U(i,j,k,8)**2.0)
!nlpnm2=nlpnm2*(0.5d0+dsign(0.5d0,masstoflux - 1.d0))
!******no dvi v*********
!nlpnm2=1.d0
!******no dvi v*********
SFE=LSFE*dsqrt(U(i,j,k,1)/(4.04d0*1.d3)) * dt * nlpnm2
SFE2=LSFE*dsqrt(U(i,j,k,1)/(4.04d0*1.d3)) * dt * 2.d0 * ndH2(i,j,k) / (2.d0 * ndH2(i,j,k)+ndH(i,j,k)+ndp(i,j,k))
SFE2=dmin1(SFE2,0.8d0)
!SFE=1.d0 * nlpnm2
Mstar1=U(i,j,k,1)*Msun*dx1*dy1*dz1*rmvratio
Mstarlp=1.d-15
Mmassivetot=0.d0

call ran0(ransfe,idum1)
nlpnm3=idint(0.5d0-dsign(0.5d0,ransfe-SFE))
nlpnm4=1
nmsv1=0.d0
rdenst1=0.d0
denst1=0.d0
nummsv=0

rhopre1=U(i,j,k,1)
do nstdm=1,(nlpmass)*idint(nlpnm2)*nlpnm3
nlpnm4=idint(0.5d0-dsign(0.5d0,Mstarlp-Mstar1))
call ran0(raninp,idum1)
call IMF(raninp,ranout)
Mstarlp=Mstarlp+ranout
Mtotrc1=Mtotrc1+ranout

nlpnm1=idint(0.5d0-dsign(0.5d0,Mmassive-ranout))
nummsv=int(nlpnm1)+nummsv !*int(nlpnm1)
nlpnm=nlpnm1*nlpnm3*nlpnm4
!write(*,*)NRANK,nlpnm,nlpnm1,nlpnm4,Mmassive,ranout,Mstarlp,Mstar1,'Massive'

Mmassivetot=Mmassivetot+ranout * dble(nlpnm)
!nid=nid+1
!nidnw=nidnw + 1 * nlpnm
nwnid1 = nwnid1 + 1.d0 * dble(nlpnm)

nidnw=idint(nwnid1)
Ustar(1,(nid+nidnw)*nlpnm)=dble(IST*Ncellx + i)*dx1!+dx1*0.5d0
Ustar(2,(nid+nidnw)*nlpnm)=dble(JST*Ncelly + j)*dy1!+dy1*0.5d0
Ustar(3,(nid+nidnw)*nlpnm)=dble(KST*Ncellz + k)*dz1!+dz1*0.5d0
Ustar(4,(nid+nidnw)*nlpnm)=U(i,j,k,2)
Ustar(5,(nid+nidnw)*nlpnm)=U(i,j,k,3)
Ustar(6,(nid+nidnw)*nlpnm)=U(i,j,k,4)
Ustar(7,(nid+nidnw)*nlpnm)=ranout !U(i,j,k,1)*dx1*dy1*dz1*SFratio
Ustar(8,(nid+nidnw)*nlpnm)=intt
Ustar(9,(nid+nidnw)*nlpnm)=dble(nid)
lQ=aq+bq*dlog10(ranout)+cq*dlog10(ranout)**2+dq*dlog10(ranout)**3+eq*dlog10(ranout)**4+fq*dlog10(ranout)**5
Ustar(10,(nid+nidnw)*nlpnm)=lQ
Ustar(11,(nid+nidnw)*nlpnm)=dble(IST*Ncellx + i)
Ustar(12,(nid+nidnw)*nlpnm)=dble(JST*Ncelly + j)
Ustar(13,(nid+nidnw)*nlpnm)=dble(KST*Ncellz + k)
Rst=10.d0*(10**(Ustar(10,(nid+nidnw)*nlpnm)-49.d0))**(1.d0/3.d0)*(U(i,j,k,1)/1.27d0/10.d0)**(-2.d0/3.d0) !loop
Ustar(14,(nid+nidnw)*nlpnm)=Rst
!Ustar(14,0)=1.d0
nstdmy=(nid+nidnw)*nlpnm
Ustar(15,(nid+nidnw)*nlpnm)=U(i,j,k,1)
!(U(i+1,j+1,k+1,1)+U(i+1,j+1,k  ,1)+U(i+1,j+1,k-1,1)&
!+U(i+1,j  ,k+1,1)+U(i+1,j  ,k  ,1)+U(i+1,j  ,k-1,1)&
!+U(i+1,j-1,k+1,1)+U(i+1,j-1,k  ,1)+U(i+1,j-1,k-1,1)&
!+U(i  ,j+1,k+1,1)+U(i  ,j+1,k  ,1)+U(i  ,j+1,k-1,1)&
!+U(i  ,j  ,k+1,1)+U(i  ,j  ,k  ,1)+U(i  ,j  ,k-1,1)&
!+U(i  ,j-1,k+1,1)+U(i  ,j-1,k  ,1)+U(i  ,j-1,k-1,1)&
!+U(i-1,j+1,k+1,1)+U(i-1,j+1,k  ,1)+U(i-1,j+1,k-1,1)&
!+U(i-1,j  ,k+1,1)+U(i-1,j  ,k  ,1)+U(i-1,j  ,k-1,1)&
!+U(i-1,j-1,k+1,1)+U(i-1,j-1,k  ,1)+U(i-1,j-1,k-1,1))/27.d0
Ustar(16,(nid+nidnw)*nlpnm)=2.2d0*(U(i,j,k,1)*ffgas)**(-1.d0/3.d0)*(ranout*clsmass)**(1.d0/3.d0) !10 -> cluster mass
!Ustar(16,(nid+nidnw)*nlpnm)=2.1d0*(U(i,j,k,1)*ffgas)**(-1.d0/3.d0)*(Ustar(7,(nid+nidnw)*nlpnm)*10.d0)**(1.d0/3.d0)
Ustar(14,(nid+nidnw)*nlpnm)=Rst+Ustar(16,(nid+nidnw)*nlpnm)
Ustar(17,(nid+nidnw)*nlpnm)=(U(i+1,j+1,k+1,1)+U(i+1,j+1,k  ,1)+U(i+1,j+1,k-1,1)&
                            +U(i+1,j  ,k+1,1)+U(i+1,j  ,k  ,1)+U(i+1,j  ,k-1,1)&
                            +U(i+1,j-1,k+1,1)+U(i+1,j-1,k  ,1)+U(i+1,j-1,k-1,1)&
                            +U(i  ,j+1,k+1,1)+U(i  ,j+1,k  ,1)+U(i  ,j+1,k-1,1)&
                            +U(i  ,j  ,k+1,1)+U(i  ,j  ,k  ,1)+U(i  ,j  ,k-1,1)&
                            +U(i  ,j-1,k+1,1)+U(i  ,j-1,k  ,1)+U(i  ,j-1,k-1,1)&
                            +U(i-1,j+1,k+1,1)+U(i-1,j+1,k  ,1)+U(i-1,j+1,k-1,1)&
                            +U(i-1,j  ,k+1,1)+U(i-1,j  ,k  ,1)+U(i-1,j  ,k-1,1)&
                            +U(i-1,j-1,k+1,1)+U(i-1,j-1,k  ,1)+U(i-1,j-1,k-1,1))*dx1*dy1*dz1*Msun
Ustar(18,(nid+nidnw)*nlpnm)=SFE
Ustar(19,(nid+nidnw)*nlpnm)=Mstarlp
Ustar(20,(nid+nidnw)*nlpnm)=dble(nummsv)


!nwnid2 = nwnid1*nlpnm
!Ustar(9,0)=dble(nidnw)
!--------check---------------
!nwnid(NRANK)=nwnid1*nlpnm
!call BC_ST_rad(nwnid1,nwnid,nlpnm)
!nidnw=idint(nwnid1)
!Ustar(9,0)=dble(nidnw)
Nstinp(NRANK)=dble(nidnw)
!write(*,*)NRANK,(nid+nidnw)*nlpnm,nid,nidnw,nwnid1,'N'
!--------check---------------
!U(i,j,k,1)=rhopre1*(1.d0-SFE)
!MGtoS(NRANK)=rhopre1*(SFE2)
Fstar(NRANK)=Mstarlp
enddo
U(i,j,k,1)   = rhopre1*(1.d0-SFE2)
ndH(i,j,k)   = ndH(i,j,k)*(1.d0-SFE2)
ndp(i,j,k)   = ndp(i,j,k)*(1.d0-SFE2)
ndH2(i,j,k)  = ndH2(i,j,k)*(1.d0-SFE2)
ndHe(i,j,k)  = ndHe(i,j,k)*(1.d0-SFE2)
ndHep(i,j,k) = ndHep(i,j,k)*(1.d0-SFE2)
ndC(i,j,k)   = ndC(i,j,k)*(1.d0-SFE2)
ndCO(i,j,k)  = ndCO(i,j,k)*(1.d0-SFE2)
ndCp(i,j,k)  = ndCp(i,j,k)*(1.d0-SFE2)
nde(i,j,k)   = nde(i,j,k)*(1.d0-SFE2)
ndtot(i,j,k) = ndtot(i,j,k)*(1.d0-SFE2)
MGtoS(NRANK) = MGtoS(NRANK)+rhopre1*(SFE2)
end do; end do; end do
call collectST()

!N_MPI(20)=1; N_MPI(1)=1; iwx = 1; iwy = 1; iwz = 1; CALL BC_MPI(2,1)
endif


if(mode==2) then
!write(*,*)NRANK,'mode2-en'
!call system_clock(time_end_c6)
!rdrgn(:,:,:,1)=1.d0
iwx=1
iwy=1
iwz=1
N_MPI(20) = 2
N_MPI(1)  = 1
N_MPI(2)  = 5
CALL BC_MPI(2,1)

allocate(STF1(1:nid,1:nid))
allocate(STF2(1:ndx-2,1:ndy-2,ndz-2))
nstloop2=0
STF2(:,:,:)=0.d0
do nstloop1=1,nid
nstloop2=1.d0+(0.5d0-dsign(0.5d0,-STF2(istf,jstf,kstf)+1.d0))
STF1(nstloop1,nstloop1)=dble(nstloop1)
STF1(nstloop1,nstloop1)=dble(nstloop1)*(0.5d0-dsign(0.5d0,-STF2(istf,jstf,kstf)+1.d0))
istf=idint((dmod(Ustar(1,nstloop1)-0.5d0*dx1,dble(Ncellx)*dx1)/dx1))+1
jstf=idint((dmod(Ustar(2,nstloop1)-0.5d0*dy1,dble(Ncelly)*dy1)/dy1))+1
kstf=idint((dmod(Ustar(3,nstloop1)-0.5d0*dz1,dble(Ncellz)*dz1)/dz1))+1

!if(rixcn.eq.-1 ) rixcn = NSPLTx-1
!if(rjycn.eq.-1 ) rjycn = NSPLTy-1
!if(rkzcn.eq.-1 ) rkzcn = NSPLTz-1
!if(rixcn.eq.NSPLTx) rixcn = 0
!if(rjycn.eq.NSPLTy) rjycn = 0
!if(rkzcn.eq.NSPLTz) rkzcn = 0

!ndcore=0.d0
!if((IST==rixcn).and.(JST==rjycn).and.(KST==rkzcn)) then
!ndcore=1.d0
!endif

STF2(istf,jstf,kstf)=10.d0**Ustar(10,nstloop)+STF2(istf,jstf,kstf)
enddo


!call system_clock(time_end_c7)

radint1=0.d0
radint2=0.d0
do nstloop=1,nid
!rdrgn(:,:,:,2)=1.d0

rixcn=idint((Ustar(1,nstloop)-0.5d0*dx1)/(dble(Ncellx)*dx1))
rjycn=idint((Ustar(2,nstloop)-0.5d0*dy1)/(dble(Ncelly)*dy1))
rkzcn=idint((Ustar(3,nstloop)-0.5d0*dz1)/(dble(Ncellz)*dz1))

if(rixcn.eq.-1 ) rixcn = NSPLTx-1
if(rjycn.eq.-1 ) rjycn = NSPLTy-1
if(rkzcn.eq.-1 ) rkzcn = NSPLTz-1
if(rixcn.eq.NSPLTx) rixcn = 0
if(rjycn.eq.NSPLTy) rjycn = 0
if(rkzcn.eq.NSPLTz) rkzcn = 0

ndcore=0.d0
if((IST==rixcn).and.(JST==rjycn).and.(KST==rkzcn)) then
ndcore=1.d0
endif

i=idint((dmod(Ustar(1,nstloop)-0.5d0*dx1,dble(Ncellx)*dx1)/dx1))+1
j=idint((dmod(Ustar(2,nstloop)-0.5d0*dy1,dble(Ncelly)*dy1)/dy1))+1
k=idint((dmod(Ustar(3,nstloop)-0.5d0*dz1,dble(Ncellz)*dz1)/dz1))+1

!Ustar(15,nstloop)=0.d0!Ustar(14,nstloop) * dble(nmIST1*nmJST1*nmKST1)
!Ustar(16,nstloop)=0.d0
!Rsv=Ustar(14,nstloop)

nlpnm26=(0.5d0-dsign(0.5d0,-U(i,j,k,1)+0.05d0))
nlpnm27=(0.5d0-dsign(0.5d0,radint(i,j,k,2)-dble(nstloop)-0.1d0))
radint1=radint(i,j,k,1)*nlpnm27+10.d48*(1.d0-nlpnm27)

meandrho=&!U(i  ,j  ,k  ,1)
         (U(i+1,j+1,k+1,1)+U(i+1,j+1,k  ,1)+U(i+1,j+1,k-1,1)&
         +U(i+1,j  ,k+1,1)+U(i+1,j  ,k  ,1)+U(i+1,j  ,k-1,1)&
         +U(i+1,j-1,k+1,1)+U(i+1,j-1,k  ,1)+U(i+1,j-1,k-1,1)&
         +U(i  ,j+1,k+1,1)+U(i  ,j+1,k  ,1)+U(i  ,j+1,k-1,1)&
         +U(i  ,j  ,k+1,1)+U(i  ,j  ,k  ,1)+U(i  ,j  ,k-1,1)&
         +U(i  ,j-1,k+1,1)+U(i  ,j-1,k  ,1)+U(i  ,j-1,k-1,1)&
         +U(i-1,j+1,k+1,1)+U(i-1,j+1,k  ,1)+U(i-1,j+1,k-1,1)&
         +U(i-1,j  ,k+1,1)+U(i-1,j  ,k  ,1)+U(i-1,j  ,k-1,1)&
         +U(i-1,j-1,k+1,1)+U(i-1,j-1,k  ,1)+U(i-1,j-1,k-1,1))/27.d0

pre16=Ustar(16,nstloop)
pre14=Ustar(14,nstloop)


 Ustar(16,nstloop)=2.1d0*(meandrho*ffgas)**(-1.d0/3.d0)*(Ustar(7,nstloop)*clsmass)**(1.d0/3.d0)*ndcore
!Ustar(16,nstloop)=2.1d0*(U(i,j,k,1)*ffgas)**(-1.d0/3.d0)*(Ustar(7,nstloop)*clsmass)**(1.d0/3.d0)*ndcore

call BC_ST_RS_one(nstloop,1,16)
Ustar(14,nstloop)=(Ustar(16,nstloop)+2.d0*dx1) * ndcore
call BC_ST_RS_one(nstloop,1,14)


do irdn= 1,iradloop

dnMPI(:,1)=0.d0
!rdnumMPI(:,:)=0.d0
!rdcellMPI(:,:)=0.d0
nmeanrd=0.d0
!nmeanrdden=0.d0
!nmcellrad=0.d0

do k = -1, Ncellz+2; do j = -1, Ncelly+2; do i = -1, Ncellx+2
!do k = 1*nkst1, Ncellz*nkst1; do j = 1*njst1, Ncelly*njst1; do i = 1*nist1, Ncellx*nist1
rdumy1=dsqrt((dble(IST*(ndx-2) + i)*dx1-Ustar(1,nstloop))**2.d0 )
!rdumy2=dble(JST*(ndy-2) + j)*dy1-Ustar(2,nstloop)
rdumy2=dsqrt((dble(JST*(ndy-2) + j)*dy1-Ustar(2,nstloop))**2.d0 )
nrdumy2=0.5d0+dsign(0.5d0,rdumy2-Lbox*0.5d0)  !dabs
!rdumy2=rdumy2-Lbox*0.5d0*nrdumy2
rdumy2=-rdumy2+Lbox*nrdumy2
!rdumy3=dble(KST*(ndz-2) + k)*dz1-Ustar(3,nstloop)
rdumy3=dsqrt((dble(KST*(ndz-2) + k)*dz1-Ustar(3,nstloop))**2.d0 )
nrdumy3=0.5d0+dsign(0.5d0,rdumy3-Lbox*0.5d0)
!rdumy3=rdumy3-Lbox*0.5d0*nrdumy3
rdumy3=-rdumy3+Lbox*nrdumy3


rdumy=dsqrt((rdumy1)**2.d0+(rdumy2)**2.d0+(rdumy3)**2.d0)
Rstdmy=10.d0*(10**(Ustar(10,nstloop)-49.d0))**(1.d0/3.d0)*(U(i,j,k,1)/1.27d0/10.d0)**(-2.d0/3.d0)
nlpnm8=(0.5d0-dsign(0.5d0,-Ustar(14,nstloop)+rdumy))
nlpnm88=(0.5d0+dsign(0.5d0,-Ustar(16,nstloop)+rdumy))
!nlpnm888=(0.5d0+dsign(0.5d0,2.d0-Ustar(14,nstloop)/Rstdmy))
nlpnm888=1.d0

!nmeanrd=nmeanrd+(U(i,j,k,1)/1.27d0)*(U(i,j,k,1)/1.27d0)*alpharcm*dx1*dy1*dz1*(3.d0*1.d18)**3.d0 *rdrgn(i,j,k,1) * nlpnm8
nmeanrd=nmeanrd+(U(i,j,k,1)/1.27d0)*(U(i,j,k,1)/1.27d0)*alpharcm*dx1*dy1*dz1*(3.d0*1.d18)**3.d0 * nlpnm8* nlpnm88 *nlpnm888
!nmeanrdden=U(i,j,k,1) * nlpnm8 * rdrgn(i,j,k,1)+nmeanrdden
!nmcellrad=nmcellrad+1.d0 *rdrgn(i,j,k,1) * nlpnm8
end do; end do; end do

dnMPI(NRANK,1)=nmeanrd
call BC_ST_rad(nmeanrd,dnMPI(0,1),1)

if(dabs(nmeanrd-10.d0**Ustar(10,nstloop)) < 10.d0**Ustar(10,nstloop) * 3.d-1) write(*,*)'conv',irdn; exit

if(nmeanrd < 10.d0**Ustar(10,nstloop)) then
!write(*,*)dabs(nmeanrd - 10.d0**Ustar(10,nstloop)), rdpcnt*10.d0**Ustar(10,nstloop),'shusoku'
nlpnm21=(0.5d0-dsign(0.5d0,-Ustar(14,nstloop)+dx1))
nlpnm22=(0.5d0-dsign(0.5d0,-dabs(nmeanrd - 10.d0**Ustar(10,nstloop)) + rdpcnt * 10.d0**Ustar(10,nstloop)))
Ustar(14,nstloop)=(Ustar(14,nstloop)+dx1*nlpnm21*nlpnm22) * ndcore
!goto 9192
endif
if(nmeanrd > 10.d0**Ustar(10,nstloop)) then
nlpnm21=(0.5d0-dsign(0.5d0,-Ustar(14,nstloop)+dx1))
 !write(*,*)dabs(nmeanrd - 10.d0**Ustar(10,nstloop)), rdpcnt*10.d0**Ustar(10,nstloop),'shusoku'
nlpnm22=(0.5d0-dsign(0.5d0,-dabs(nmeanrd - 10.d0**Ustar(10,nstloop)) + rdpcnt * 10.d0**Ustar(10,nstloop)))
Ustar(14,nstloop)=(Ustar(14,nstloop)-dx1*nlpnm22*nlpnm21) * ndcore
!goto 9192
endif
call BC_ST_RS_one(nstloop,1,14)

!if(nmeanrd > 10.d0**Ustar(10,nstloop)) then
enddo


!nlpnmpre14=(0.5d0-dsign(0.5d0,-Ustar(14,nstloop)+0.1d0*pre14))
!nlpnmpre16=(0.5d0+dsign(0.5d0,-Ustar(14,nstloop)+5.d0 *pre14))
nlpnmpre14=(0.5d0-dsign(0.5d0,-Ustar(14,nstloop)+dx1*2.d0))
!nlpnmpre16=(0.5d0+dsign(0.5d0,-Ustar(14,nstloop)+dx1*1.5d0))
!Ustar(14,nstloop)=(Ustar(14,nstloop)*nlpnmpre14*nlpnmpre16+pre14*(1.d0-nlpnmpre14*nlpnmpre16)) * ndcore
!Ustar(16,nstloop)=(Ustar(16,nstloop)*nlpnmpre14*nlpnmpre16+pre16*(1.d0-nlpnmpre14*nlpnmpre16)) * ndcore
Ustar(14,nstloop)=(Ustar(14,nstloop)*nlpnmpre14+dx1*2.d0*(1.d0-nlpnmpre14)) * ndcore
Ustar(16,nstloop)=(Ustar(16,nstloop)*nlpnmpre14+dx1*2.d0*(1.d0-nlpnmpre14)) * ndcore
call BC_ST_RS_one(nstloop,1,14)
call BC_ST_RS_one(nstloop,1,16)

!if((pre14>1.d0).and.(0.4d0*pre14>Ustar(14,nstloop)))then
!Ustar(14,nstloop)=pre14
!Ustar(16,nstloop)=pre16
!endif
!if((pre14>1.d0).and.(3.d0*pre14<Ustar(14,nstloop)))then
!Ustar(14,nstloop)=pre14
!Ustar(16,nstloop)=pre16
!endif


!nlpnm11=(0.5d0-dsign(0.5d0,+Ustar(14,nstloop)-dx1))
!*****cell heat up*****
 !volratio = (Ustar(14,nstloop)/dx1)**3
 !U(i,j,k,5)=(1.d4/T0)*kb*(ndp(i,j,k)+ndH(i,j,k)+ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))*ndcore * volratio + &
 !U(i,j,k,5)*ndcore *(1.d0-volratio)
!i=idint((dmod(Ustar(1,nstloop)-0.5d0*dx1,dble(Ncellx)*dx1)/dx1))+1
!j=idint((dmod(Ustar(2,nstloop)-0.5d0*dy1,dble(Ncelly)*dy1)/dy1))+1
!k=idint((dmod(Ustar(3,nstloop)-0.5d0*dz1,dble(Ncellz)*dz1)/dz1))+1
!U(i,j,k,5)=(1.d4/T0)*kb*(ndp(i,j,k)+ndH(i,j,k)+ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))*ndcore
!*****cell heat up*****

!write(*,*) Ustar(14,nstloop),Ustar(15,nstloop),Ustar(1,nstloop),Ustar(2,nstloop),Ustar(3,nstloop),nstloop,'-Rst-nibun-'

do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
!do k = 1*nkst1, Ncellz*nkst1; do j = 1*njst1, Ncelly*njst1; do i = 1*nist1, Ncellx*nist1
!rdumy1=dble(IST*(ndx-2) + i)*dx1-Ustar(1,nstloop)
rdumy1=dsqrt((dble(IST*(ndx-2) + i)*dx1-Ustar(1,nstloop))**2.d0 )
!rdumy2=dble(JST*(ndy-2) + j)*dy1-Ustar(2,nstloop)
rdumy2=dsqrt((dble(JST*(ndy-2) + j)*dy1-Ustar(2,nstloop))**2.d0 )
nrdumy2=0.5d0+dsign(0.5d0,rdumy2-Lbox*0.5d0)
!rdumy2=rdumy2-Lbox*0.5d0*nrdumy2
rdumy2=-rdumy2+Lbox*nrdumy2
!rdumy3=dble(KST*(ndz-2) + k)*dz1-Ustar(3,nstloop)
rdumy3=dsqrt((dble(KST*(ndz-2) + k)*dz1-Ustar(3,nstloop))**2.d0 )
nrdumy3=0.5d0+dsign(0.5d0,rdumy3-Lbox*0.5d0)
!rdumy3=rdumy3-Lbox*0.5d0*nrdumy3
rdumy3=-rdumy3+Lbox*nrdumy3
rdumy=dsqrt((rdumy1)**2.d0+(rdumy2)**2.d0+(rdumy3)**2.d0)

nlpnm11=(0.5d0-dsign(0.5d0,-Ustar(14,nstloop)+rdumy))
Rstdmy=10.d0*(10**(Ustar(10,nstloop)-49.d0))**(1.d0/3.d0)*(U(i,j,k,1)/1.27d0/10.d0)**(-2.d0/3.d0)
!nlpnm888=(0.5d0+dsign(0.5d0,2.d0-Ustar(14,nstloop)/Rstdmy))
nlpnm888=1.d0
rdrgn(i,j,k,2)=1.d0*(1.d0-nlpnm11)
U(i,j,k,5)=(1.d4/T0)*kb*(ndp(i,j,k)+ndH(i,j,k)+ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))*nlpnm11*nlpnm888+U(i,j,k,5)*(1.d0-nlpnm11*nlpnm888)
end do; end do; end do

iwx=1
iwy=1
iwz=1
N_MPI(20) = 1
N_MPI(1)  = 5
CALL BC_MPI(2,1)

enddo

!call system_clock(time_end_c8)
!write(*,*)NRANK,'mode2-en'
deallocate(STF1,STF2)
endif



if(mode==3) then
do nstloop=1,nid
rsix=idint((Ustar(1,nstloop))/(dble(Ncellx)*dx1))
rsjy=idint((Ustar(2,nstloop))/(dble(Ncelly)*dy1))
rskz=idint((Ustar(3,nstloop))/(dble(Ncellz)*dz1))

rsixc=dble(idint(dmod(Ustar(1,nstloop),dble(Ncellx)*dx1)/dx1)+1)
rsjyc=dble(idint(dmod(Ustar(2,nstloop),dble(Ncellx)*dy1)/dy1)+1)
rskzc=dble(idint(dmod(Ustar(3,nstloop),dble(Ncellx)*dz1)/dz1)+1)

xx1=0.d0
xx2=0.d0
xx3=0.d0
vv1=0.d0
vv2=0.d0
vv3=0.d0
if((IST==rsix).and.(JST==rsjy).and.(KST==rskz)) then
call runge(Ustar(1,nstloop),Ustar(2,nstloop),Ustar(3,nstloop),Ustar(4,nstloop),dt,nstloop,rsixc,rsjyc,rskzc,1,Ustar(7,nstloop)/Msun)
call runge(Ustar(2,nstloop),Ustar(3,nstloop),Ustar(1,nstloop),Ustar(5,nstloop),dt,nstloop,rsixc,rsjyc,rskzc,2,Ustar(7,nstloop)/Msun)
call runge(Ustar(3,nstloop),Ustar(1,nstloop),Ustar(2,nstloop),Ustar(6,nstloop),dt,nstloop,rsixc,rsjyc,rskzc,3,Ustar(7,nstloop)/Msun)
!Ustar(1,nstloop)=Ustar(1,nstloop)+Ustar(4,nstloop)*dt/dx
!Ustar(2,nstloop)=Ustar(2,nstloop)+Ustar(5,nstloop)*dt/dy
!Ustar(3,nstloop)=Ustar(3,nstloop)+Ustar(6,nstloop)*dt/dz
!call BC_ST(nstloop)
if(Ustar(1,nstloop)>Lbox) then
Ustar(1,nstloop)=Ustar(1,nstloop)-Lbox
endif
if(Ustar(2,nstloop)>Lbox) then
Ustar(2,nstloop)=Ustar(2,nstloop)-Lbox
endif
if(Ustar(3,nstloop)>Lbox) then
Ustar(3,nstloop)=Ustar(3,nstloop)-Lbox
endif

if(Ustar(1,nstloop)<0.d0) then
Ustar(1,nstloop)=Ustar(1,nstloop)+Lbox
endif
if(Ustar(2,nstloop)<0.d0) then
Ustar(2,nstloop)=Ustar(2,nstloop)+Lbox
endif
if(Ustar(3,nstloop)<0.d0) then
Ustar(3,nstloop)=Ustar(3,nstloop)+Lbox
endif

xx1=Ustar(1,nstloop)
xx2=Ustar(2,nstloop)
xx3=Ustar(3,nstloop)
vv1=Ustar(4,nstloop)
vv2=Ustar(5,nstloop)
vv3=Ustar(6,nstloop)
!write(*,*) xx1,xx2,xx3,vv1,vv2,vv3,'mode3-renew'
endif
Ustar(1,nstloop)=xx1
Ustar(2,nstloop)=xx2
Ustar(3,nstloop)=xx3
Ustar(4,nstloop)=vv1
Ustar(5,nstloop)=vv2
Ustar(6,nstloop)=vv3
enddo

call BC_ST_RS_tot(1,1)
call BC_ST_RS_tot(1,2)
call BC_ST_RS_tot(1,3)
call BC_ST_RS_tot(1,4)
call BC_ST_RS_tot(1,5)
call BC_ST_RS_tot(1,6)

!call BC_ST()

!call BC_ST_RS(nstloop,1,16)

endif

end SUBROUTINE feedback

!SUBROUTINE runge(xint,vint,tint,dt,nid)
SUBROUTINE runge(x1int,x2int,x3int,vint,dt,nid,rsixc,rsjyc,rskzc,mode,mst1)
  !implicit none
  double precision :: tstart,xstart,tend,vstart,x1int,x2int,x3int,vint
  double precision :: tnx,xnx,vnx,dt,rsixc,rsjyc,rskzc
  double precision :: tstp0,tstp1,tstp2,tstp3
  double precision :: kx0,kx1,kx2,kx3
  double precision :: kv0,kv1,kv2,kv3
  double precision :: f0,f1,f2,f3
  double precision :: v0,v1,v2,v3
  double precision :: funcf,funcv,mst1
  !double precision,allocatable,dimension(:) :: tint,xint,vint,xexa,vexa
  integer :: i,iend,nid,mode
  character(21) :: fname='/Users/maeda/Desktop/'

  !tstart=0.d0
  !tend  =100.d0
  !xstart=0.d0
  !vstart=1.d0
  !tstart=0.d0
  !dt=0.1d0
  !nid=100000
  
  !allocate(tint(1:nid),xint(1:nid),vint(1:nid),xexa(1:nid),vexa(1:nid))
  !allocate(xint(1:nid),vint(1:nid))
  !tint(:)=0.d0
  !xint(:)=0.d0
  !vint(:)=0.d0
 
 
  !t=tstart
  !x=xstart
  !v=vstart
  !do i=1,nid
    !write(10,'(2f13.5)') t,x
 
    !tstp0=tint(i)
    !i=int(xint/dx1)
    !j=int(yint/dy1)
    !k=int(zint/dz1)

    kv0=vint
    kx0=x1int
    !f0=funcf(int(rsixc),int(rsjyc),int(rskzc),mode,kx0,kv0,mst1)
    !v0=funcv(int(rsixc),int(rsjyc),int(rskzc),mode,kx0,kv0)
    f0=funcf(kx0,x2int,x3int,mode,kx0,kv0,mst1)
    v0=funcv(idint(rsixc),idint(rsjyc),idint(rskzc),mode,kx0,kv0)
 
    !tstp1=tint(i)+dt/2.0
    kv1=vint+f0*dt/2.0
    kx1=x1int+v0*dt/2.0
    !f1=funcf(int(rsixc),int(rsjyc),int(rskzc),mode,kx1,kv1,mst1)
    !v1=funcv(int(rsixc),int(rsjyc),int(rskzc),mode,kx1,kv1)
    f1=funcf(kx1,x2int,x3int,mode,kx1,kv1,mst1)
    v1=funcv(idint(rsixc),idint(rsjyc),idint(rskzc),mode,kx1,kv1)

    !tstp2=tint(i)+dt/2.0
    kv2=vint+f1*dt/2.0
    kx2=x1int+v1*dt/2.0
    !f2=funcf(int(rsixc),int(rsjyc),int(rskzc),mode,kx2,kv2,mst1)
    !v2=funcv(int(rsixc),int(rsjyc),int(rskzc),mode,kx2,kv2)
    f2=funcf(kx0,x2int,x3int,mode,kx2,kv2,mst1)
    v2=funcv(idint(rsixc),idint(rsjyc),idint(rskzc),mode,kx2,kv2)
 
    !tstp3=tint(i)+dt
    kv3=vint+f2*dt
    kx3=x1int+v2*dt
    !f3=funcf(int(rsixc),int(rsjyc),int(rskzc),mode,kx3,kv3,mst1)
    !v3=funcv(int(rsixc),int(rsjyc),int(rskzc),mode,kx3,kv3)
    f3=funcf(kx2,x2int,x3int,mode,kx3,kv3,mst1)
    v3=funcv(idint(rsixc),idint(rsjyc),idint(rskzc),mode,kx3,kv3)
 
    vnx=vint+(f0+f1*2.0+f2*2.0+f3)*dt/6.0
    xnx=x1int+(v0+v1*2.0+v2*2.0+v3)*dt/6.0

    !tnx=tint(i)+dt

    !tint(i)=tnx
    x1int=xnx
    vint=vnx

    !write(*,*) vnx,xnx,kv0,kv1,kv2,kv3,kx0,kx1,kx2,kx3,mode, 'runge'
    !xexa(i)=vstart*dsin(t)
    !vexa(i)=-vstart*dcos(t)
 
  !end do

  !open(10,file=fname//'runge4.dat')
  !do i=1,nid
  !write(10,*) tint(i),',',xint(i),',',vint(i),',',xexa(i),',',vexa(i)
  !end do
  !close(10)

  !deallocate(xint,vint)
  !deallocate(tint,xint,vint,xexa,vexa)
 
  !close(10)
 
end SUBROUTINE runge

 
!function funcf(i,j,k,mode,kx,kv,mst1)
function funcf(x1,x2,x3,mode,kx,kv,mst1)
  USE comvar
  USE slfgrv
  !double precision :: t,x,v
  double precision :: kx,kv,mst1,x1,x2,x3,xmod,ymod,zmod,delx1,delx2,delx
  double precision :: funcf
  integer :: i,j,k,mode,ii,jj,kk
 
  !funcf=(1.0-x**2)
  !funcf=-x!9.8d0!(1.0-x**2)

  !U(i,j,k,2) = U(i,j,k,2) - dt * ( -Phi(i+2,j,k)+8.d0*Phi(i+1,j,k)-8.d0*Phi(i-1,j,k)+Phi(i-2,j,k) ) * dxi *0.5d0
  !U(i,j,k,3) = U(i,j,k,3) - dt * ( -Phi(i,j+2,k)+8.d0*Phi(i,j+1,k)-8.d0*Phi(i,j-1,k)+Phi(i,j-2,k) ) * dxi *0.5d0
  !U(i,j,k,4) = U(i,j,k,4) - dt * ( -Phi(i,j,k+2)+8.d0*Phi(i,j,k+1)-8.d0*Phi(i,j,k-1)+Phi(i,j,k-2) ) * dxi *0.5d0

  !U(i,j,k,2) = -  ( -Phi(i+2,j,k)+8.d0*Phi(i+1,j,k)-8.d0*Phi(i-1,j,k)+Phi(i-2,j,k) ) * dxi *0.5d0
  !U(i,j,k,3) = -  ( -Phi(i,j+2,k)+8.d0*Phi(i,j+1,k)-8.d0*Phi(i,j-1,k)+Phi(i,j-2,k) ) * dyi *0.5d0
  !U(i,j,k,4) = -  ( -Phi(i,j,k+2)+8.d0*Phi(i,j,k+1)-8.d0*Phi(i,j,k-1)+Phi(i,j,k-2) ) * dzi *0.5d0
  if (mode==1)then
  i=mod(idint(x1/dx1),ndx-2)
  j=mod(idint(x2/dy1),ndy-2)
  k=mod(idint(x3/dz1),ndz-2)
  ii=idint(x1/dx1)
  jj=idint(x2/dy1)
  kk=idint(x3/dz1)
  xmod=x1-dble(ii)*dx1
  ymod=x2-dble(jj)*dy1
  zmod=x3-dble(kk)*dz1
  xmod=xmod/dx1
  ymod=ymod/dy1
  zmod=zmod/dz1
  !funcf = -  ( Phi(i+1,j,k)-2.d0*Phi(i,j,k)+Phi(i-1,j,k) ) / dx1 / dx1
  !funcf = -  ( Phi(i+1,j,k)-Phi(i-1,j,k) )*0.5d0 / dx1! *mst1!/Msun
  !funcf =( -( Phi(i+1,j+1,k)-Phi(i,j+1,k))/dx1 -( Phi(i+1,j,k)-Phi(i,j,k))/dx1)*0.5d0+( -( Phi(i+1,j,k+1)-Phi(i,j,k+1))/dx1 -( Phi(i+1,j,k)-Phi(i,j,k))/dx1)*0.5d0
  !funcf =( -( Phi(i+1,j+1,k)-Phi(i,j+1,k))/dx1 -( Phi(i+1,j,k)-Phi(i,j,k))/dx1)*0.5d0+( -( Phi(i+1,j,k+1)-Phi(i,j,k+1))/dx1 -( Phi(i+1,j,k)-Phi(i,j,k))/dx1)*0.5d0

delx1=( -( Phi(i+1,j+1,k)-Phi(i-1,j+1,k))*0.5d0/dx1 * (ymod) -( Phi(i+1,j,k)-Phi(i-1,j,k))*0.5d0/dx1 * (1.d0-ymod))*0.5d0 &
     +( -( Phi(i+1,j,k+1)-Phi(i-1,j,k+1))*0.5d0/dx1 * (zmod) -( Phi(i+1,j,k)-Phi(i-1,j,k))*0.5d0/dx1 * (1.d0-zmod))*0.5d0
delx2=( -( Phi(i+2,j+1,k)-Phi(i  ,j+1,k))*0.5d0/dx1 * (ymod) -( Phi(i+2,j,k)-Phi(i  ,j,k))*0.5d0/dx1 * (1.d0-ymod))*0.5d0 &
     +( -( Phi(i+2,j,k+1)-Phi(i  ,j,k+1))*0.5d0/dx1 * (zmod) -( Phi(i+2,j,k)-Phi(i  ,j,k))*0.5d0/dx1 * (1.d0-zmod))*0.5d0
  !delx1=-(Phi(i+1,j)-Phi(i-1,j))*0.5d0/dx * (1.d0-ymod) -(Phi(i+1,j+1)-Phi(i-1,j+1))*0.5d0/dx * (ymod)
  !delx2=-(Phi(i+2,j)-Phi(i  ,j))*0.5d0/dx * (1.d0-ymod) -(Phi(i+2,j+1)-Phi(i  ,j+1))*0.5d0/dx * (ymod)
  delx=delx1* (1.d0-xmod)+delx2* xmod
  funcf=delx
  !write(*,*) Phi(i+1,j,k),Phi(i-1,j,k),funcf,'phix'
  endif
  if (mode==2)then
  i=mod(idint(x3/dx1),ndx-2)
  j=mod(idint(x1/dy1),ndy-2)
  k=mod(idint(x2/dz1),ndz-2)
  ii=idint(x3/dx1)
  jj=idint(x1/dy1)
  kk=idint(x2/dz1)
  xmod=x3-dble(ii)*dx1
  ymod=x1-dble(jj)*dy1
  zmod=x2-dble(kk)*dz1
  xmod=xmod/dx1
  ymod=ymod/dy1
  zmod=zmod/dz1
  !funcf = -  ( Phi(i,j+1,k)-2.d0*Phi(i,j,k)+Phi(i,j-1,k) ) / dy1 / dy1
  !funcf = -  ( Phi(i,j+1,k)-Phi(i,j-1,k) )*0.5d0 / dy1! *mst1
  !funcf =( -( Phi(i+1,j+1,k)-Phi(i+1,j,k))/dy1 -( Phi(i,j+1,k)-Phi(i,j,k))/dy1)*0.5d0+( -( Phi(i,j+1,k+1)-Phi(i,j,k+1))/dy1 -( Phi(i,j+1,k)-Phi(i,j,k))/dy1)*0.5d0

delx1=( -( Phi(i+1,j+1,k)-Phi(i+1,j-1,k))*0.5d0/dy1 * (xmod) -( Phi(i,j+1,k)-Phi(i,j-1,k))*0.5d0/dy1 * (1.d0-xmod))*0.5d0 &
     +( -( Phi(i,j+1,k+1)-Phi(i,j-1,k+1))*0.5d0/dy1 * (zmod) -( Phi(i,j+1,k)-Phi(i,j-1,k))*0.5d0/dy1 * (1.d0-zmod))*0.5d0
delx2=( -( Phi(i+1,j+2,k)-Phi(i+1,j  ,k))*0.5d0/dy1 * (xmod) -( Phi(i,j+2,k)-Phi(i,j  ,k))*0.5d0/dy1 * (1.d0-xmod))*0.5d0 &
     +( -( Phi(i,j+2,k+1)-Phi(i,j  ,k+1))*0.5d0/dy1 * (zmod) -( Phi(i,j+2,k)-Phi(i,j  ,k))*0.5d0/dy1 * (1.d0-zmod))*0.5d0

  delx=delx1* (1.d0-ymod)+delx2* ymod
  funcf=delx
  !write(*,*) Phi(i,j+1,k),Phi(i,j-1,k),funcf,'phiy'
  endif
  if (mode==3)then
  i=mod(idint(x2/dx1),ndx-2)
  j=mod(idint(x3/dy1),ndy-2)
  k=mod(idint(x1/dz1),ndz-2)
  ii=idint(x2/dx1)
  jj=idint(x3/dy1)
  kk=idint(x1/dz1)
  xmod=x2-dble(ii)*dx1
  ymod=x3-dble(jj)*dy1
  zmod=x1-dble(kk)*dz1
  xmod=xmod/dx1
  ymod=ymod/dy1
  zmod=zmod/dz1
  !funcf = -  ( Phi(i,j,k+1)-2.d0*Phi(i,j,k)+Phi(i,j,k-1) ) / dz1 / dz1
  !funcf = -  ( Phi(i,j,k+1)-Phi(i,j,k-1) )*0.5d0 / dz1! *mst1
  !funcf =( -( Phi(i+1,j,k+1)-Phi(i+1,j,k))/dz1 -( Phi(i,j,k+1)-Phi(i,j,k))/dz1)*0.5d0+( -( Phi(i,j+1,k+1)-Phi(i,j+1,k))/dz1 -( Phi(i,j,k+1)-Phi(i,j,k))/dz1)*0.5d0

delx1=( -( Phi(i+1,j,k+1)-Phi(i+1,j,k-1))*0.5d0/dz1 * (xmod) -( Phi(i,j,k+1)-Phi(i,j,k-1))*0.5d0/dz1 * (1.d0-xmod))*0.5d0 &
     +( -( Phi(i,j+1,k+1)-Phi(i,j+1,k-1))*0.5d0/dz1 * (ymod) -( Phi(i,j,k+1)-Phi(i,j,k-1))*0.5d0/dz1 * (1.d0-ymod))*0.5d0
delx2=( -( Phi(i+1,j,k+2)-Phi(i+1,j,k  ))*0.5d0/dz1 * (xmod) -( Phi(i,j,k+2)-Phi(i,j,k  ))*0.5d0/dz1 * (1.d0-xmod))*0.5d0 &
     +( -( Phi(i,j+1,k+2)-Phi(i,j+1,k  ))*0.5d0/dz1 * (ymod) -( Phi(i,j,k+2)-Phi(i,j,k  ))*0.5d0/dz1 * (1.d0-ymod))*0.5d0

  delx=delx1* (1.d0-zmod)+delx2* zmod
  funcf=delx
  endif
 
  return
end function funcf


function funcv(i,j,k,mode,kx,kv)
  !double precision :: t,x,v
  double precision :: funcv
  integer :: i,j,k,mode
  double precision :: kx,kv
 

  funcv = kv
  return
end function funcv



SUBROUTINE IMF(raninp,ranout)
integer :: nmesh,cnc1,cnc2,nran!,idum
integer :: i
double precision :: raninp,ranout
double precision :: mmax,mmin,nsisu,scale,slide,dmy,Ncumcnc1,Ncumcnc2
double precision :: M1,M2,coefflow,coeffmid,coeffhig,pwllow,pwlmid,pwlhig,coefflow2,coeffmid2,coeffhig2
double precision, allocatable, dimension(:) :: M,lgM,dN,N,lgN,Ncum,Mran,ranini,ts1,ts2,ts3,ts4
double precision :: dm,lgmmax,lgmmin,ddN,Ntot
character(21) :: dir='/Users/maeda/Desktop/'

!-------paramerter------
!nmesh = 1000
!nran=10000000
mmin = 0.01d0
mmax = 150.d0
Ntot = 100.d0
!slide=5.d0
M1=0.08d0
M2=0.5d0
pwllow=-0.3d0
pwlmid=-1.3d0
pwlhig=-2.3d0
coeffhig=1.d0
coeffmid=(coeffhig*(M2)**pwlhig)/((M2)**pwlmid)
coefflow=(coeffmid*(M1)**pwlmid)/((M1)**pwllow)
coefflow2=1.d0
coeffmid2=(coefflow2*(M1)**(pwllow))/((M1)**(pwlmid))
coeffhig2=(coeffmid2*(M2)**(pwlmid))/((M2)**(pwlhig))
!-------paramerter------

!allocate(M(0:nmesh),lgM(0:nmesh),dN(0:nmesh),N(0:nmesh),lgN(0:nmesh),Ncum(0:nmesh),Mran(1:nran),ranini(1:nran))
!allocate(ts1(1:nran),ts2(1:nran),ts3(1:nran),ts4(1:nran))

!lgmmax = dlog10(mmax)
!lgmmin = dlog10(mmin)
!dm = (lgmmax-lgmmin)/dble(nmesh)
!cnc1=int((dlog10(M1)-dlog10(mmin))/dm)
!cnc2=int((dlog10(M2)-dlog10(mmin))/dm)



Ntot=coefflow*((M1  )**(pwllow+1.d0)/(pwllow+1.d0)-(mmin)**(pwllow+1.d0)/(pwllow+1.d0))+&
     coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0)-(M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
     coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2  )**(pwlhig+1.d0)/(pwlhig+1.d0))


Ncumcnc2=coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2     )**(pwlhig+1.d0)/(pwlhig+1.d0))
Ncumcnc1=coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0)-(M1     )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
         coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2     )**(pwlhig+1.d0)/(pwlhig+1.d0))


Ncumcnc2=Ncumcnc2/Ntot
Ncumcnc1=Ncumcnc1/Ntot



if (raninp>Ncumcnc1) then
ranout=((coefflow*((M1  )**(pwllow+1.d0)/(pwllow+1.d0))+&
          coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0))-coeffmid*((M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
          coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0))-coeffhig*((M2  )**(pwlhig+1.d0)/(pwlhig+1.d0)) &
          -Ntot*raninp)/coefflow*(pwllow+1.d0))**(1.d0/(pwllow+1.d0))

else if ((Ncumcnc1>raninp).and.(raninp>Ncumcnc2)) then
ranout=((coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
          coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2  )**(pwlhig+1.d0)/(pwlhig+1.d0)) &
          -Ntot*raninp)*(pwlmid+1.d0)/coeffmid)**(1.d0/(pwlmid+1.d0))

else if (Ncumcnc2>raninp) then
ranout=((coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0))-Ntot*raninp)*(pwlhig+1.d0)/coeffhig)**(1.d0/(pwlhig+1.d0))

endif


!deallocate(M,lgM,dN,N,lgN,Ncum,Mran,ranini)
!deallocate(ts1,ts2,ts3,ts4)
end SUBROUTINE IMF

!do i = 0, Ncell+1
!   N(i) = U(i,1)/1.27d0
!  Tn(i) = U(i,5)/( kb*(ndp(i)+ndH(i)+ndH2(i)+ndHe(i)+ndHep(i)) )
!  Pn(i) = U(i,5)
!end do
!do i = 0, Ncell
!  rtTx  = dsqrt(  ( dx(i)*Tn(i+1)+dx(i+1)*Tn(i) ) /( dx(i) + dx(i+1) )  )
!  Qx(i) = Kcond * rtTx * ( Tn(i+1)-Tn(i) )*2.d0 /( dx(i) + dx(i+1) )
!end do
!do i = 1, Ncell
!!----- Cooling ---------------------------------------------------------
!  Call Fcool( CooL,N(i),Tn(i))
!  gammi1 =   3.d0*(ndH(i)+ndp(i)+ndHe(i)+ndHep(i))+5.d0*ndH2(i)
!  gammi1 = ( 2.d0*(ndH(i)+ndp(i)+ndHe(i)+ndHep(i))+2.d0*ndH2(i) )/gammi1
!  if( dt .le. 0.2d0*Pn(i)/(gammi1*dabs(CooL)) ) then
!    U(i,5) = U(i,5) - gammi1*CooL*dt*0.5d0 !explicit
!  else
!    Call IMC( U(i,5),ndH(i)+ndp(i)+ndHe(i)+ndHep(i)+ndH2(i),dt*0.5d0,N(i) ) !implicit
!  end if
!!----- Conduction ------------------------------------------------------
!  U(i,5) = U(i,5) + gammi1*dt*0.5d0*(Qx(i)-Qx(i-1))/dx(i)
!  U(i,5) = dmax1(U(i,5),pmin)
!  U(i,5) = dmin1(U(i,5),pmax)
!end do


SUBROUTINE collectST()
USE comvar
USE mpivar
USE slfgrv
USE fedvar
INCLUDE 'mpif.h'
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
!double precision :: stMPI(1:valstar,0:numstar,0:NPE-1),u1(1:valstar,0:numstar)
integer :: num1,totst,stcnt,nwid,NRANKdm
integer :: rsix,rsjy,rskz
DOUBLE PRECISION, dimension(:,:,:), allocatable  :: stMPI
DOUBLE PRECISION, dimension(:,:), allocatable  :: u1

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

do Nroot=0,NPE-1
  CALL MPI_BCAST(Nstinp(Nroot),1,MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do

totst=0
do Nroot=0,NPE-1
totst=totst+idint(Nstinp(Nroot))
enddo

allocate(stMPI(1:valstar,0:totst,0:NPE-1),u1(1:valstar,0:totst))

stMPI(:,:,:)=0.d0
u1(:,:)=0.d0
do i=0,totst!nid+nidnw
do j=1,valstar
  stMPI(j,i,NRANK)=Ustar(j,i+nid)
  !write(*,*) NRANK,Ustar(j,i),'totst-1'
end do;end do

do Nroot=0,NPE-1
  CALL MPI_BCAST(stMPI(1,0,Nroot),(valstar)*(totst+1),MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do


stcnt=0
do Nroot=0,NPE-1
if (idint(Nstinp(Nroot)).ne.0) then
do i=1,idint(Nstinp(Nroot))
stcnt=stcnt+1
do iii=1,valstar
u1(iii,stcnt)=stMPI(iii,stcnt,Nroot)
!write(*,*) NRANK,u1(iii,stcnt),int(stMPI(9,0,Nroot)),iii,stcnt,nid+i,'totst12'
enddo
enddo
endif
enddo

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

!if(totst.ne.0) then
!do i=nid-totst+1,nid; do j=1,valstar!numstar
!  Ustar(j,i)=u1(j,i)
!  Ustar(9,i)=dble(i)
!write(*,*) NRANK,Ustar(j,i),u1(j,i),nid,j,i,'totst2'
!end do;end do
!endif

!if(totst.ne.0) then
do i=1,totst
do j=1,valstar!numstar
  Ustar(j,i+nid)=u1(j,i)
  !Ustar(9,i+nid)=dble(i+nid)
!write(*,*) NRANK,Ustar(j,i),u1(j,i),nid,j,i,'totst2'
end do
Ustar(9,i+nid)=dble(i+nid)
end do
nid=nid+totst
!endif
!nid=int(Ustar(9,0))+stcnt

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!nid=nid+totst
!stMPI(9,0,:)=0.d0

deallocate(stMPI,u1)
END SUBROUTINE collectST


SUBROUTINE BC_ST()
USE comvar
USE mpivar
USE slfgrv
USE fedvar
INCLUDE 'mpif.h'
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
double precision :: stMPI(1:valstar,0:numstar,0:NPE-1),u1(1:valstar,0:numstar)
integer :: num1,totst,stcnt,nwid,NRANKdm
integer :: rsix,rsjy,rskz,nmb

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

do j=1,valstar; do i=0,numstar
  stMPI(j,i,NRANK)=Ustar(j,i)
  u1(j,i)=Ustar(j,i)
end do;end do

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
do Nroot=0,NPE-1
  CALL MPI_BCAST(stMPI(1,0,Nroot),(valstar)*(numstar+1),MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do

! i=1,numstar
!rsix=int((Ustar(1,i))/(dble(Ncellx)*dx1))
!rsjy=int((Ustar(2,i))/(dble(Ncelly)*dy1))
!rskz=int((Ustar(3,i))/(dble(Ncellz)*dz1))
!NRANKdm = rsix + rsjy*NSPLTx + rskz*NSPLTx*NSPLTy
!do j=1,valstar
!u1(j,i)=stMPI(j,i,NRANKdm)
!end do
!end do

!do Nroot=0,NPE-1
!u1(iii,stcnt)=stMPI(iii,nid+i,Nroot)
!!write(*,*) NRANK,u1(iii,stcnt),int(stMPI(9,0,Nroot)),iii,stcnt,nid+i,'totst12'
!enddo

END SUBROUTINE BC_ST


SUBROUTINE BC_ST_RS(nmb,nlpnm,num)
USE comvar
USE mpivar
USE slfgrv
USE fedvar
INCLUDE 'mpif.h'
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
double precision :: stMPI(1:valstar,0:numstar,0:NPE-1),u1(1:valstar,0:numstar)
integer :: num1,totst,stcnt,nwid,NRANKdm
integer :: rsix,rsjy,rskz,nlpnm,num

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!write(*,*)'BC_ST_RS'
!if(nlpnm==1) then

j=num
do i=0,numstar
  stMPI(j,i,NRANK)=Ustar(j,i)
  u1(j,i)=Ustar(j,i)
end do

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
do Nroot=0,(NPE)*nlpnm-1
  CALL MPI_BCAST(stMPI(num,0,Nroot),(numstar+1),MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do


Ustar(num,nmb)=0.d0
do Nroot=0,NPE-1
  Ustar(num,nmb)=Ustar(num,nmb)+stMPI(num,nmb,Nroot)
  !write(*,*)num,nmb,Nroot,stMPI(num,nmb,Nroot)
end do
END SUBROUTINE BC_ST_RS

SUBROUTINE BC_ST_RS_one(nmb,nlpnm,num)
USE comvar
USE mpivar
USE slfgrv
USE fedvar
INCLUDE 'mpif.h'
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
double precision :: stMPI(0:NPE-1),u1(1:valstar,0:numstar)
integer :: num1,totst,stcnt,nwid,NRANKdm
integer :: rsix,rsjy,rskz,nlpnm,num

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!write(*,*)'BC_ST_RS'
!if(nlpnm==1) then

j=num
i=nmb
!do i=0,numstar
  stMPI(NRANK)=Ustar(j,i)
  u1(j,i)=Ustar(j,i)
!end do

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
do Nroot=0,(NPE)*nlpnm-1
  CALL MPI_BCAST(stMPI(Nroot),1,MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do


Ustar(num,nmb)=0.d0
do Nroot=0,NPE-1
  Ustar(num,nmb)=Ustar(num,nmb)+stMPI(Nroot)
  !write(*,*)num,nmb,Nroot,stMPI(num,nmb,Nroot)
end do
END SUBROUTINE BC_ST_RS_one

SUBROUTINE BC_ST_RS_tot(nlpnm,num)
USE comvar
USE mpivar
USE slfgrv
USE fedvar
INCLUDE 'mpif.h'
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
double precision :: u1(1:valstar,0:numstar)
integer :: num1,totst,stcnt,nwid,NRANKdm
integer :: rsix,rsjy,rskz,nlpnm,num
double precision,allocatable,dimension(:,:) :: stMPI!(1:valstar,0:numstar,0:NPE-1)

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!write(*,*)'BC_ST_RS'
!if(nlpnm==1) then

allocate(stMPI(0:nid,0:NPE-1))

j=num
do i=0,nid
  stMPI(i,NRANK)=Ustar(j,i)
  u1(j,i)=Ustar(j,i)
end do

do Nroot=0,(NPE)*nlpnm-1
  CALL MPI_BCAST(stMPI(0,Nroot),(nid+1),MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do


Ustar(num,:)=0.d0
do Nroot=0,NPE-1
do i=0,nid
  Ustar(num,i)=Ustar(num,i)+stMPI(i,Nroot)
  !write(*,*)num,nmb,Nroot,stMPI(num,nmb,Nroot)
enddo
end do

deallocate(stMPI)

END SUBROUTINE BC_ST_RS_tot


SUBROUTINE BC_ST_rad(Urad,UrdMPI,nlpnm)
USE comvar
USE mpivar
USE slfgrv
USE fedvar
INCLUDE 'mpif.h'
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
double precision  :: UrdMPI(0:NPE-1),Urad
integer Nroot,nlpnm

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!if(nlpnm==1) then
do Nroot=0,(NPE)*nlpnm-1
  CALL MPI_BCAST(UrdMPI(Nroot),1,MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do
!endif

Urad=0.d0
do Nroot=0,NPE-1
Urad=UrdMPI(Nroot)+Urad
end do
!endif

!do i=1,numstar
!  Ustar(14,i)=stMPI(14,i,NPE-1)
!end do

!do i=1,numstar
!rsix=int((Ustar(1,i))/(dble(Ncellx)*dx1))
!rsjy=int((Ustar(2,i))/(dble(Ncelly)*dy1))
!rskz=int((Ustar(3,i))/(dble(Ncellz)*dz1))
!NRANKdm = rsix + rsjy*NSPLTx + rskz*NSPLTx*NSPLTy
!do j=1,valstar
!u1(j,i)=stMPI(j,i,NRANKdm)
!end do
!end do

END SUBROUTINE BC_ST_rad

SUBROUTINE BC_ST_stprc(Urad1,Urad2,Urad3,UrdMPI,nlpnm)
USE comvar
USE mpivar
USE slfgrv
USE fedvar
INCLUDE 'mpif.h'
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
double precision  :: UrdMPI(0:NPE-1,3),Urad1,Urad2,Urad3
integer Nroot,nlpnm


!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

!if(nlpnm==1) then
do Nroot=0,(NPE)*nlpnm-1
  CALL MPI_BCAST(UrdMPI(Nroot,1),1,MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do
do Nroot=0,(NPE)*nlpnm-1
  CALL MPI_BCAST(UrdMPI(Nroot,2),1,MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do
do Nroot=0,(NPE)*nlpnm-1
  CALL MPI_BCAST(UrdMPI(Nroot,3),1,MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do
!endif

Urad1=0.d0
Urad2=0.d0
Urad3=0.d0
do Nroot=0,NPE-1
Urad1=UrdMPI(Nroot,1)+Urad1
Urad2=UrdMPI(Nroot,2)+Urad2
Urad3=UrdMPI(Nroot,3)+Urad3
end do
!endif

!do i=1,numstar
!  Ustar(14,i)=stMPI(14,i,NPE-1)
!end do

!do i=1,numstar
!rsix=int((Ustar(1,i))/(dble(Ncellx)*dx1))
!rsjy=int((Ustar(2,i))/(dble(Ncelly)*dy1))
!rskz=int((Ustar(3,i))/(dble(Ncellz)*dz1))
!NRANKdm = rsix + rsjy*NSPLTx + rskz*NSPLTx*NSPLTy
!do j=1,valstar
!u1(j,i)=stMPI(j,i,NRANKdm)
!end do
!end do

END SUBROUTINE BC_ST_stprc
