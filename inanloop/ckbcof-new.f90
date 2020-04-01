program main
  implicit none
  !DOUBLE PRECISION, dimension(:,:,:), allocatable :: Kyo
  !DOUBLE PRECISION, dimension(:,:,:,:), allocatable :: U,PHI,Kyo
  DOUBLE PRECISION, dimension(:), allocatable :: U,PHI,Kyo1,Kyo2
  DOUBLE PRECISION, dimension(:,:,:,:,:), allocatable :: Ux,Uy,Uz,V
  integer nx,ny,nz,i,j,k,NRANK,core,val,ival,NSPLTx,NSPLTy,NSPLTz,Ncellx,Ncelly,Ncellz
  integer cornum,rdmax,itr,itrmax
  integer :: inp=0,usecore,looplast,opnm,err,numloop
  character(5) sav
  character(4) chnum
  character(3) frank1,frank2
  character(2) n20_ch
  character(50) name
  INTEGER :: IST,JST,KST,LEFT,RIGT,BOTM,TOP,UP,DOWN
  integer cntnum,N20,cntnumend
  !character, dimension(:), allocatable :: stname


  !******parameter********
  cornum=292
  rdmax=4*10
  core=512
  val=17
  Ncellx=64
  Ncelly=64
  Ncellz=64
  looplast=100
  name='BC_MPI'
  cntnumend=47
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


  LEFT = cornum - 1            ; if(IST.eq.0       ) LEFT = cornum + (NSPLTx-1)
  RIGT = cornum + 1            ; if(IST.eq.NSPLTx-1) RIGT = cornum - (NSPLTx-1)
  BOTM = cornum - NSPLTx       ; if(JST.eq.0       ) BOTM = cornum + NSPLTx*(NSPLTy-1)
  TOP  = cornum + NSPLTx       ; if(JST.eq.NSPLTy-1) TOP  = cornum - NSPLTx*(NSPLTy-1)
  DOWN = cornum - NSPLTx*NSPLTy; if(KST.eq.0       ) DOWN = cornum + NSPLTx*NSPLTy*(NSPLTz-1)
  UP   = cornum + NSPLTx*NSPLTy; if(KST.eq.NSPLTz-1) UP   = cornum - NSPLTx*NSPLTy*(NSPLTz-1)

  allocate(Kyo1(rdmax),Kyo2(rdmax))

  write(frank1,'(i3.3)') cornum
  write(frank2,'(i3.3)') LEFT
  open(unit=210,file=trim(name)//frank1//'.dat',FORM='FORMATTED',iostat=err)
  open(unit=220,file=trim(name)//frank2//'.dat',FORM='FORMATTED',iostat=err)
  do itr=1,itrmax
     Kyo1(:)=0.d0
     Kyo2(:)=0.d0
     read(210,*,iostat=err) (Kyo1(numloop),numloop=1,rdmax)
     if(err>0)then
        goto 2990
     end if
2990 continue
     read(220,*,iostat=err) (Kyo2(numloop),numloop=1,rdmax)
     if(err>0)then
        goto 2999
     end if
2999 continue

     do i=1,rdmax
        if(Kyo1(i).ne.Kyo2(i))then
           write(*,*) 'err',Kyo1(i),Kyo2(i)
        end if
     end do
  end do

end program
