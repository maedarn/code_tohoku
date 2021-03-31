program computeraw
 implicit none
 !integer :: nx, ny, nz
 integer(kind=8) :: rec_len
 integer cnt10,jump,ts,timejump
 !integer :: i,j,k
 !parameter(nx = 100, ny = 20, nz = 50)
 !real, dimension(nx,ny,nz) :: array
 real(4), allocatable :: u(:,:,:,:) !, v(:,:,:,:)
 real(4), allocatable ::  v(:,:,:,:)
 real(8), allocatable ::  pbr(:,:),pbl(:,:)
 real(8), allocatable ::  pbrall(:,:),pblall(:,:)
 double precision , allocatable :: x(:),y(:),z(:),inputlength(:)
 real(4), allocatable ::Tn(:,:,:)
 real(4) mass,nth
 real(8) mass1
 integer(4), allocatable :: x1shock(:,:),x2shock(:,:)
 integer :: val,core=0,NSPLTx,NSPLTy,NSPLTz,times=0,mesh,IST,KST,JST,lengthx,&
      lengthy,lengthz,corexy,midle,corestart,lasttime,length
 integer :: i,j,k,kz=1,initime, Ncellx, Ncelly, Ncellz,usecore,dataname
 character(50) filename,data
 character(3) NPE,time,cntctj
 character(6) cntc
 !character(36) :: dir='/glv0/maedarn/test-grvwave/test-FFT/'
 !character(34) :: dir='/glv0/maedarn/test-grvwave/PHIINI/'
 !character(36) :: dir='/glv0/maedarn/test-telegraph/PHIINI/'
 !character(64) ::
 !dir='/glv0/maedarn/test-telegraph/cg1-T20-mash128-L100-dipole/PHIINI/'
character(63) :: dir='/glv0/maedarn/test-telegraph/cg1-T20-mash128-L100-eq-cn/PHIINI/' 
 !character(61) :: dir='/glv0/maedarn/test-telegraph/kp-1-10-cg1-mash256-L100/PHIINI/'
 !character(54) :: dir='/glv0/maedarn/test-grvwave/test-compre-exa-userhomean/'
double precision :: adiff=0.375d0,cg=1.d0,Tdiff=20.d0,dx1
DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G
real(4), allocatable ::  diff(:,:,:,:)

  !******parameter********
  core=64
  mesh=32+2
  !val=25
  !val=27
  val=8
  !val=28
  jump=1
  ts=1
  initime=1
  lasttime=101
  Ncellx=32
  Ncelly=32
  Ncellz=32
  timejump=1
  nth=10.e0
  dx1=100.d0/dble(128)
  !***********************

  mass=0.e0
  mass1=0.d0
if(core.eq.4)    then; NSPLTx = 2; NSPLTy = 2; NSPLTz = 1; end if
if(core.eq.8)    then; NSPLTx = 2; NSPLTy = 2; NSPLTz = 2; end if
if(core.eq.16)   then; NSPLTx = 4; NSPLTy = 2; NSPLTz = 2; end if
if(core.eq.32)   then; NSPLTx = 4; NSPLTy = 4; NSPLTz = 2; end if
if(core.eq.64)   then; NSPLTx = 4; NSPLTy = 4; NSPLTz = 4; end if
if(core.eq.128)  then; NSPLTx = 8; NSPLTy = 4; NSPLTz = 4; end if
if(core.eq.256)  then; NSPLTx = 8; NSPLTy = 8; NSPLTz = 4; end if
if(core.eq.512)  then; NSPLTx = 8; NSPLTy = 8; NSPLTz = 8; end if
if(core.eq.1024) then; NSPLTx = 8; NSPLTy = 8; NSPLTz =16; end if

   lengthx=(mesh-2)*NSPLTx
   lengthy=(mesh-2)*NSPLTy
   lengthz=(mesh-2)*NSPLTz

   ALLOCATE(v(1:lengthx,1:lengthy,1:lengthz,1:val))
   ALLOCATE(Tn(1:lengthx,1:lengthy,1:lengthz))
   ALLOCATE(u(-1:mesh,-1:mesh,-1:mesh,1:val))
   ALLOCATE(x1shock(1:lengthy,1:lengthz),x2shock(1:lengthy,1:lengthx))
   ALLOCATE(diff(1+1:lengthx-1,1+1:lengthy-1,1+1:lengthz-1,1:8))

  
write(*,*)'1',initime,lasttime,timejump
do times=initime,lasttime,timejump
     write(*,*)times
       write(cntc,'(i6.6)') times
      do usecore=0,core-1,1
         IST = mod(usecore,NSPLTx)
         KST = usecore/(NSPLTx*NSPLTy)
         JST = usecore/NSPLTx-NSPLTy*KST
         write(*,*) IST,JST,KST,time,usecore
         write(NPE,'(I3.3)') usecore
         open(unit=150,file=dir//'PHI'//cntc//NPE//'.DAT',FORM='UNFORMATTED')
         do k = -1, Ncellz+2 !*************************************************** +1 いる？ **********************************
            do j = -1, Ncelly+2
               do i= -1, Ncellx+2
                  read(150) u(i,j,k,1) ,u(i,j,k,2),u(i,j,k,3),u(i,j,k,4),u(i,j,k,5), &
                    u(i,j,k,6),u(i,j,k,7),u(i,j,k,8)!, &
                    !u(i,j,k,9),u(i,j,k,10)!,u(i,j,k,11)!,u(i,j,k,12), &
                    !u(i,j,k,13),u(i,j,k,14),u(i,j,k,15),u(i,j,k,16), u(i,j,k,17), &
                    !u(i,j,k,18),u(i,j,k,19),u(i,j,k,20),u(i,j,k,21), u(i,j,k,22), &
                    !u(i,j,k,23),u(i,j,k,24),u(i,j,k,25), & 
                    !u(i,j,k,26),u(i,j,k,27)!,u(i,j,k,28)
                end do
            end do
         end do
         close(150)
         do k = 1, Ncellz
            do j = 1, Ncelly
               do i = 1, Ncellx
                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),1)=u(i,j,k,1)
                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),2)=u(i,j,k,2)
                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),3)=u(i,j,k,3)
                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),4)=u(i,j,k,4)
                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),5)=u(i,j,k,5)
                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),6)=u(i,j,k,6)
                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),7)=u(i,j,k,7)
                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),8)=u(i,j,k,8)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),9)=u(i,j,k,9)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),10)=u(i,j,k,10)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),11)=u(i,j,k,11)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),12)=u(i,j,k,12)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),13)=u(i,j,k,13)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),14)=u(i,j,k,14)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),15)=u(i,j,k,15)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),16)=u(i,j,k,16)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),17)=u(i,j,k,17)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),18)=u(i,j,k,18)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),19)=u(i,j,k,19)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),20)=u(i,j,k,20)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),21)=u(i,j,k,21)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),22)=u(i,j,k,22)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),23)=u(i,j,k,23)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),24)=u(i,j,k,24)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),25)=u(i,j,k,25)
                  
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),26)=u(i,j,k,26)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),27)=u(i,j,k,27)
                  !v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),28)=u(i,j,k,28)
            end do
         end do
!write(*,*)'Cmetal',u(k,30,30,14),u(k,30,30,15)
      end do
   end do

   do k = 1+1, lengthz-1
        do j = 1+1, lengthy-1
            do i = 1+1, lengthx-1
                diff(i,j,k,1)=adiff*v(i+1,j+1,k,1)+adiff*v(i-1,j-1,k,1)+(adiff-0.5d0)*v(i+1,j-1,k,1)+(adiff-0.5d0)*v(i-1,j+1,k,1) &
                +(4.d0*adiff-1.d0)*v(i,j,k,1)+(-2.d0*adiff+0.5d0)*v(i+1,j,k,1)+(-2.d0*adiff+0.5d0)*v(i,j+1,k,1)+&
                (-2.d0*adiff+0.5d0)*v(i-1,j,k,1)+(-2.d0*adiff+0.5d0)*v(i,j-1,k,1)
                diff(i,j,k,2)=adiff*v(i,j+1,k+1,1)+adiff*v(i,j-1,k-1,1)+(adiff-0.5d0)*v(i,j+1,k-1,1)+(adiff-0.5d0)*v(i,j-1,k+1,1) &
                +(4.d0*adiff-1.d0)*v(i,j,k,1)+(-2.d0*adiff+0.5d0)*v(i,j+1,k,1)+(-2.d0*adiff+0.5d0)*v(i,j,k+1,1)+&
                (-2.d0*adiff+0.5d0)*v(i,j-1,k,1)+(-2.d0*adiff+0.5d0)*v(i,j,k-1,1)
                diff(i,j,k,3)=adiff*v(i+1,j,k+1,1)+adiff*v(i-1,j,k-1,1)+(adiff-0.5d0)*v(i-1,j,k+1,1)+(adiff-0.5d0)*v(i+1,j,k-1,1) &
                +(4.d0*adiff-1.d0)*v(i,j,k,1)+(-2.d0*adiff+0.5d0)*v(i,j,k+1,1)+(-2.d0*adiff+0.5d0)*v(i+1,j,k,1)+&
                (-2.d0*adiff+0.5d0)*v(i,j,k-1,1)+(-2.d0*adiff+0.5d0)*v(i-1,j,k,1)
                diff(i,j,k,4)=v(i,j,k,1) &
                -2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*diff(i,j,k,1)/dx1/dx1 &
                -2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*diff(i,j,k,2)/dx1/dx1 &
                -2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*diff(i,j,k,3)/dx1/dx1 &
                -cg*cg*4.d0*Tdiff*Tdiff*G4pi*v(i,j,k,8)
                diff(i,j,k,5)=adiff*v(i+1,j+1,k,6)+adiff*v(i-1,j-1,k,6)+(adiff-0.5d0)*v(i+1,j-1,k,6)+(adiff-0.5d0)*v(i-1,j+1,k,6) &
                +(4.d0*adiff-1.d0)*v(i,j,k,6)+(-2.d0*adiff+0.5d0)*v(i+1,j,k,6)+(-2.d0*adiff+0.5d0)*v(i,j+1,k,6)+&
                (-2.d0*adiff+0.5d0)*v(i-1,j,k,6)+(-2.d0*adiff+0.5d0)*v(i,j-1,k,6)
                diff(i,j,k,6)=adiff*v(i,j+1,k+1,6)+adiff*v(i,j-1,k-1,6)+(adiff-0.5d0)*v(i,j+1,k-1,6)+(adiff-0.5d0)*v(i,j-1,k+1,6) &
                +(4.d0*adiff-1.d0)*v(i,j,k,6)+(-2.d0*adiff+0.5d0)*v(i,j+1,k,6)+(-2.d0*adiff+0.5d0)*v(i,j,k+1,6)+&
                (-2.d0*adiff+0.5d0)*v(i,j-1,k,6)+(-2.d0*adiff+0.5d0)*v(i,j,k-1,6)
                diff(i,j,k,7)=adiff*v(i+1,j,k+1,6)+adiff*v(i-1,j,k-1,6)+(adiff-0.5d0)*v(i-1,j,k+1,6)+(adiff-0.5d0)*v(i+1,j,k-1,6) &
                +(4.d0*adiff-1.d0)*v(i,j,k,6)+(-2.d0*adiff+0.5d0)*v(i,j,k+1,6)+(-2.d0*adiff+0.5d0)*v(i+1,j,k,6)+&
                (-2.d0*adiff+0.5d0)*v(i,j,k-1,6)+(-2.d0*adiff+0.5d0)*v(i-1,j,k,6)
                diff(i,j,k,8)=v(i,j,k,6) &
                -2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*diff(i,j,k,5)/dx1/dx1 &
                -2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*diff(i,j,k,6)/dx1/dx1 &
                -2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*diff(i,j,k,7)/dx1/dx1 &
                -cg*cg*4.d0*Tdiff*Tdiff*G4pi*v(i,j,k,8)
            end do
        end do
    end do



   write(cntctj,'(i3.3)') times/timejump
   !open(100,file=dir//'All/All'//cntctj//'.DAT',access='stream',FORM='UNFORMATTED')
   !open(100,file=dir//'All/'//'All'//cntc//'.DAT',access='stream',FORM='UNFORMATTED')
   open(100,file=dir//'diff'//cntc//'.DAT',access='stream',FORM='UNFORMATTED')
   do k=1+1,lengthz-1,jump
      do j=1+1,lengthy-1,jump
         do i=1+1,lengthx-1,jump
         write(100) diff(i,j,k,1) ,diff(i,j,k,2),diff(i,j,k,3),diff(i,j,k,4),diff(i,j,k,5),diff(i,j,k,6),&
                    diff(i,j,k,7),diff(i,j,k,8)!,v(i,j,k,9),v(i,j,k,10)!,v(i,j,k,11)!,v(i,j,k,12),&
                    !v(i,j,k,13),v(i,j,k,14),v(i,j,k,15),v(i,j,k,16),v(i,j,k,17),v(i,j,k,18),&
                    !v(i,j,k,19),v(i,j,k,20),v(i,j,k,21),v(i,j,k,22),v(i,j,k,23),v(i,j,k,24),&
                    !v(i,j,k,25), &
                    !v(i,j,k,26),v(i,j,k,27)!,v(i,j,k,28)
         end do
      end do
      write(*,*) k,u(64,64,64,1),u(34,34,34,1)!,data
   end do
   close(100)
end do


DEALLOCATE(v)
DEALLOCATE(u)
DEALLOCATE(x1shock,x2shock)

end program computeraw
