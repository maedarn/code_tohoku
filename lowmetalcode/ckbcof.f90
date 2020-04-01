program main
  implicit none
  DOUBLE PRECISION, dimension(:,:,:,:), allocatable :: U
  DOUBLE PRECISION, dimension(:,:,:,:,:,:), allocatable :: Ux,Uy,Uz
  integer nx,ny,nz,i,j,k,NRANK,core,val,ival,NSPLTx,NSPLTy,NSPLTz,Ncellx,Ncelly,Ncellz
  integer :: inp=0,usecore,looplast,opnm
  character(5) sav
  character(3) frank
  character(50) name
  INTEGER :: IST,JST,KST,LEFT,RIGT,BOTM,TOP,UP,DOWN
  !character, dimension(:), allocatable :: stname


 !******parameter********
  core=512
  val=17
  Ncellx=64
  Ncelly=64
  Ncellz=64
  looplast=100
  name='BC_MPI'
  !***********************
  !ALLOCATE(U(-1:Ncellx,-1:Ncelly,-1:Ncellz,1:val))
  !ALLOCATE(Ux(1:2,-1:Ncelly,-1:Ncellz,1:val,1:4))
  !ALLOCATE(Uy(-1:Ncellx,1:2,-1:Ncellz,1:val,1:4))
  !ALLOCATE(Uz(-1:Ncellx,-1:Ncelly,1:2,1:val,1:4))
  !ALLOCATE(stname(len_trim(name)))
  !stname=trim(name)

  !Ux(:,:,:,:,:)=0.d0
  !Uy(:,:,:,:,:)=0.d0
  !Uz(:,:,:,:,:)=0.d0
  if(core.eq.4)    then
     NSPLTx = 2; NSPLTy = 2; NSPLTz = 1
  end if
  if(core.eq.8)    then
     NSPLTx = 2; NSPLTy = 2; NSPLTz = 2
  end if
  if(core.eq.16)   then
     NSPLTx = 4; NSPLTy = 2; NSPLTz = 2
  end if
  if(core.eq.32)   then
     NSPLTx = 4; NSPLTy = 4; NSPLTz = 2
  end if
  if(core.eq.64)   then
     NSPLTx = 4; NSPLTy = 4; NSPLTz = 4
  end if
  if(core.eq.128)  then
     NSPLTx = 8; NSPLTy = 4; NSPLTz = 4
  end if
  if(core.eq.256)  then
     NSPLTx = 8; NSPLTy = 8; NSPLTz = 4
  end if
  if(core.eq.512)  then
     NSPLTx = 8; NSPLTy = 8; NSPLTz = 8
  end if
  if(core.eq.1024) then
     NSPLTx = 8; NSPLTy = 8; NSPLTz =16
  end if

do opnm=1,3
   if(opnm==1) then
      name='BC_MPI'
      val=17
      ALLOCATE(Ux(1:2,-1:Ncelly+2,-1:Ncellz+2,1:val,1:4,0:core-1))
      ALLOCATE(Uy(-1:Ncellx+2,1:2,-1:Ncellz+2,1:val,1:4,0:core-1))
      ALLOCATE(Uz(-1:Ncellx+2,-1:Ncelly+2,1:2,1:val,1:4,0:core-1))
      !Ux(:,:,:,:,:)=0.d0
      !Uy(:,:,:,:,:)=0.d0
      !Uz(:,:,:,:,:)=0.d0
   elseif(opnm==2) then
      name='BC_MPI_OT'
      val=17
      ALLOCATE(Ux(1:2,-1:Ncelly+2,-1:Ncellz+2,1:val,1:4,0:core-1))
      ALLOCATE(Uy(-1:Ncellx+2,1:2,-1:Ncellz+2,1:val,1:4,0:core-1))
      ALLOCATE(Uz(-1:Ncellx+2,-1:Ncelly+2,1:2,1:val,1:4,0:core-1))
      !Ux(:,:,:,:,:)=0.d0
      !Uy(:,:,:,:,:)=0.d0
      !Uz(:,:,:,:,:)=0.d0
   elseif(opnm==3) then
      name='PHI_MPI'
      val=1
      ALLOCATE(Ux(1:2,-1:Ncelly+2,-1:Ncellz+2,1:val,1:4,0:core-1))
      ALLOCATE(Uy(-1:Ncellx+2,1:2,-1:Ncellz+2,1:val,1:4,0:core-1))
      ALLOCATE(Uz(-1:Ncellx+2,-1:Ncelly+2,1:2,1:val,1:4,0:core-1))
      !Ux(:,:,:,:,:)=0.d0
      !Uy(:,:,:,:,:)=0.d0
      !Uz(:,:,:,:,:)=0.d0
   end if


   do inp=0,looplast
      Ux(:,:,:,:,:,:)=0.d0
      Uy(:,:,:,:,:,:)=0.d0
      Uz(:,:,:,:,:,:)=0.d0
   do usecore=0,core-1
write(*,*) 'in',usecore,inp,opnm
write(sav,'(i5.5)') inp
write(frank,'(i3.3)') usecore
!open(12,file='inBC_MPIxr'//frank//sav//'.DAT',FORM='UNFORMATTED')
!open(17,file='inBC_MPIxl'//frank//sav//'.DAT',FORM='UNFORMATTED')
!open(13,file='inBC_MPIyr'//frank//sav//'.DAT',FORM='UNFORMATTED')
!open(14,file='inBC_MPIyl'//frank//sav//'.DAT',FORM='UNFORMATTED')
!open(15,file='inBC_MPIzr'//frank//sav//'.DAT',FORM='UNFORMATTED')
!open(16,file='inBC_MPIzl'//frank//sav//'.DAT',FORM='UNFORMATTED')

!open(22,file='bcBC_MPIxr'//frank//sav//'.DAT',FORM='UNFORMATTED')
!open(27,file='bcBC_MPIxl'//frank//sav//'.DAT',FORM='UNFORMATTED')
!open(23,file='bcBC_MPIyr'//frank//sav//'.DAT',FORM='UNFORMATTED')
!open(24,file='bcBC_MPIyl'//frank//sav//'.DAT',FORM='UNFORMATTED')
!open(25,file='bcBC_MPIzr'//frank//sav//'.DAT',FORM='UNFORMATTED')
!open(26,file='bcBC_MPIzl'//frank//sav//'.DAT',FORM='UNFORMATTED')

open(12,file='in'//trim(name)//'xr'//frank//sav//'.DAT',FORM='UNFORMATTED', status="old", err=100)
open(17,file='in'//trim(name)//'xl'//frank//sav//'.DAT',FORM='UNFORMATTED', status="old", err=100)
open(13,file='in'//trim(name)//'yr'//frank//sav//'.DAT',FORM='UNFORMATTED', status="old", err=100)
open(14,file='in'//trim(name)//'yl'//frank//sav//'.DAT',FORM='UNFORMATTED', status="old", err=100)
open(15,file='in'//trim(name)//'zr'//frank//sav//'.DAT',FORM='UNFORMATTED', status="old", err=100)
open(16,file='in'//trim(name)//'zl'//frank//sav//'.DAT',FORM='UNFORMATTED', status="old", err=100)

open(22,file='bc'//trim(name)//'xr'//frank//sav//'.DAT',FORM='UNFORMATTED', status="old", err=100)
open(27,file='bc'//trim(name)//'xl'//frank//sav//'.DAT',FORM='UNFORMATTED', status="old", err=100)
open(23,file='bc'//trim(name)//'yr'//frank//sav//'.DAT',FORM='UNFORMATTED', status="old", err=100)
open(24,file='bc'//trim(name)//'yl'//frank//sav//'.DAT',FORM='UNFORMATTED', status="old", err=100)
open(25,file='bc'//trim(name)//'zr'//frank//sav//'.DAT',FORM='UNFORMATTED', status="old", err=100)
open(26,file='bc'//trim(name)//'zl'//frank//sav//'.DAT',FORM='UNFORMATTED', status="old", err=100)

do k = -1,Ncellz+2
   do j = -1,Ncelly+2

      do i = 1,2
        read(17) (Ux(i,j,k,ival,2,usecore),ival=1,val)
     !end do

     !do i = Ncellx-2,Ncellx-1
        read(12) (Ux(i,j,k,ival,3,usecore),ival=1,val)
     !end do

     !do i = Ncellx+1,Ncellx+2
        read(22) (Ux(i,j,k,ival,4,usecore),ival=1,val)
     !end do

     !do i = -1,0
        read(27) (Ux(i,j,k,ival,1,usecore),ival=1,val)
     end do

  end do
end do

do i = -1,Ncellx+2
   do k = -1,Ncellz+2

      do j = 1,2
        read(14) (Uy(i,j,k,ival,2,usecore),ival=1,val)
     !end do

     !do j = Ncelly-2,Ncelly-1
        read(13) (Uy(i,j,k,ival,3,usecore),ival=1,val)
     !end do

     !do j = Ncelly+1,Ncelly+2
        read(23) (Uy(i,j,k,ival,4,usecore),ival=1,val)
     !end do

     !do j = -1,0
        read(24) (Uy(i,j,k,ival,1,usecore),ival=1,val)
     end do

  end do
end do

do j = -1,Ncelly+2
   do i = -1,Ncellx+2

      do k = 1,2
        read(16) (Uz(i,j,k,ival,2,usecore),ival=1,val)
     !end do

     !do k = Ncellz-2,Ncellz-1
        read(15) (Uz(i,j,k,ival,3,usecore),ival=1,val)
     !end do

     !do k = Ncellz+1,Ncellz+2
        read(25) (Uz(i,j,k,ival,4,usecore),ival=1,val)
     !end do

     !do k = -1,0
        read(26) (Uz(i,j,k,ival,1,usecore),ival=1,val)
     end do

  end do
end do

close(12)
close(13)
close(14)
close(15)
close(16)
close(17)

close(22)
close(23)
close(24)
close(25)
close(26)
close(27)
!open(12,file='BC_MPIxr'//NRANK//sav//'.DAT')
!inp=inp+1

!if(NSPLTx.ne.0) then
!100 continue
end do
do usecore=0,core-1
   IST = mod(usecore,NSPLTx); KST = usecore/(NSPLTx*NSPLTy); JST = usecore/NSPLTx-NSPLTy*KST
   LEFT = usecore - 1            ; if(IST.eq.0       ) LEFT = usecore + (NSPLTx-1)
   RIGT = usecore + 1            ; if(IST.eq.NSPLTx-1) RIGT = usecore - (NSPLTx-1)
   BOTM = usecore - NSPLTx       ; if(JST.eq.0       ) BOTM = usecore + NSPLTx*(NSPLTy-1)
   TOP  = usecore + NSPLTx       ; if(JST.eq.NSPLTy-1) TOP  = usecore - NSPLTx*(NSPLTy-1)
   DOWN = usecore - NSPLTx*NSPLTy; if(KST.eq.0       ) DOWN = usecore + NSPLTx*NSPLTy*(NSPLTz-1)
   UP   = usecore + NSPLTx*NSPLTy; if(KST.eq.NSPLTz-1) UP   = usecore - NSPLTx*NSPLTy*(NSPLTz-1)
   write(*,*) 'ck',usecore
   do k = -1,Ncellz+2
      do j = -1,Ncelly+2
         do i = 1,2
            do ival=1,val
               if((Ux(i,j,k,ival,1,IST).ne.Ux(i,j,k,ival,3,LEFT)).and.(IST.ne.0)) then
                  write(*,*) 'ERR1',i,j,k,ival,usecore,inp,opnm,Ux(i,j,k,ival,1,IST),Ux(i,j,k,ival,3,LEFT),&
                       Ux(i,j,k,ival,1,IST)-Ux(i,j,k,ival,3,LEFT)
               end if
               if((Ux(i,j,k,ival,4,IST).ne.Ux(i,j,k,ival,2,RIGT)).and.(IST.ne.NSPLTx-1)) then
                  write(*,*) 'ERR2',i,j,k,ival,usecore,inp,opnm,Ux(i,j,k,ival,4,IST),Ux(i,j,k,ival,2,RIGT),&
                       Ux(i,j,k,ival,1,IST)-Ux(i,j,k,ival,3,LEFT)
               end if
            end do
         end do
      end do
   end do

   do i = -1,Ncellx+2
      do k = -1,Ncellz+2
         do j = 1,2
            do ival=1,val
               if((Uy(i,j,k,ival,1,JST).ne.Uy(i,j,k,ival,3,BOTM))) then
                  write(*,*) 'ERR3',i,j,k,ival,usecore,inp,opnm,Uy(i,j,k,ival,1,JST),Uy(i,j,k,ival,3,BOTM),&
                       Ux(i,j,k,ival,1,IST)-Ux(i,j,k,ival,3,LEFT)
               end if
               if((Uy(i,j,k,ival,4,JST).ne.Uy(i,j,k,ival,2,TOP ))) then
                  write(*,*) 'ERR4',i,j,k,ival,usecore,inp,opnm,Uy(i,j,k,ival,4,JST),Uy(i,j,k,ival,2,TOP ),&
                       Ux(i,j,k,ival,1,IST)-Ux(i,j,k,ival,3,LEFT)
               end if
            end do
         end do
      end do
   end do

   do j = -1,Ncelly+2
      do i = -1,Ncellx+2
         do k = 1,2
            do ival=1,val
               if((Uz(i,j,k,ival,1,KST).ne.Uz(i,j,k,ival,3,DOWN))) then
                  write(*,*) 'ERR5',i,j,k,ival,usecore,inp,opnm,Uz(i,j,k,ival,1,KST),Uz(i,j,k,ival,3,DOWN),&
                       Ux(i,j,k,ival,1,IST)-Ux(i,j,k,ival,3,LEFT)
               end if
               if((Uz(i,j,k,ival,4,KST).ne.Uz(i,j,k,ival,2,UP  ))) then
                  write(*,*) 'ERR6',i,j,k,ival,usecore,inp,opnm,Uz(i,j,k,ival,4,KST),Uz(i,j,k,ival,2,UP  ),&
                       Ux(i,j,k,ival,1,IST)-Ux(i,j,k,ival,3,LEFT)
               end if
            end do
         end do
      end do
   end do

end do
!if(NSPLTx.ne.0) then
100 continue
end do
end do

end program main
