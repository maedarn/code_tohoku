program main
  implicit none
  !DOUBLE PRECISION, dimension(:,:,:), allocatable :: Kyo
  DOUBLE PRECISION, dimension(:,:,:,:), allocatable :: U,PHI,Kyo
  DOUBLE PRECISION, dimension(:,:,:,:,:), allocatable :: Ux,Uy,Uz,V
  integer nx,ny,nz,i,j,k,NRANK,core,val,ival,NSPLTx,NSPLTy,NSPLTz,Ncellx,Ncelly,Ncellz
  integer :: inp=0,usecore,looplast,opnm,err,numloop
  character(5) sav
  character(4) chnum
  character(3) frank
  character(2) n20_ch
  character(50) name
  INTEGER :: IST,JST,KST,LEFT,RIGT,BOTM,TOP,UP,DOWN
  integer cntnum,N20,cntnumend
  !character, dimension(:), allocatable :: stname


  !******parameter********
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


  !  LEFT = core - 1            ; if(IST.eq.0       ) LEFT = core + (NSPLTx-1)
  !  RIGT = core + 1            ; if(IST.eq.NSPLTx-1) RIGT = core - (NSPLTx-1)
  !  BOTM = core - NSPLTx       ; if(JST.eq.0       ) BOTM = core + NSPLTx*(NSPLTy-1)
  !  TOP  = core + NSPLTx       ; if(JST.eq.NSPLTy-1) TOP  = core - NSPLTx*(NSPLTy-1)
  !  DOWN = core - NSPLTx*NSPLTy; if(KST.eq.0       ) DOWN = core + NSPLTx*NSPLTy*(NSPLTz-1)
  !  UP   = core + NSPLTx*NSPLTy; if(KST.eq.NSPLTz-1) UP   = core - NSPLTx*NSPLTy*(NSPLTz-1)

  allocate(Kyo(1:2,0:Ncelly+1,0:Ncellz+1,0:core-1))
  do opnm=0,core-1
     write(frank,'(i3.3)') opnm
     open(unit=100,file='phibc'//frank//'.dat',FORM='UNFORMATTED')
     do k=0,Ncellz+1
        do j=0,Ncelly+1
           read(100) Kyo(1,j,k,opnm),Kyo(2,j,k,opnm)
        end do
     end do
     write(*,*)Kyo(1,1,1,opnm),Kyo(2,1,1,opnm),opnm
  end do

  do opnm=0,core-1
      write(*,*) '------------phibc------------ :',opnm,'/',511
     IST = mod(opnm,NSPLTx); KST = opnm/(NSPLTx*NSPLTy); JST = opnm/NSPLTx-NSPLTy*KST
     LEFT = opnm - 1            ; if(IST.eq.0       ) LEFT = opnm + (NSPLTx-1)
     RIGT = opnm + 1            ; if(IST.eq.NSPLTx-1) RIGT = opnm - (NSPLTx-1)
     BOTM = opnm - NSPLTx       ; if(JST.eq.0       ) BOTM = opnm + NSPLTx*(NSPLTy-1)
     TOP  = opnm + NSPLTx       ; if(JST.eq.NSPLTy-1) TOP  = opnm - NSPLTx*(NSPLTy-1)
     DOWN = opnm - NSPLTx*NSPLTy; if(KST.eq.0       ) DOWN = opnm + NSPLTx*NSPLTy*(NSPLTz-1)
     UP   = opnm + NSPLTx*NSPLTy; if(KST.eq.NSPLTz-1) UP   = opnm - NSPLTx*NSPLTy*(NSPLTz-1)

     do i=1,2
        do k = 1,Ncellz
           if((Kyo(i,1,k ,opnm).ne.Kyo(i,Ncelly+1,k ,BOTM))) then
              write(*,*) 'ERR-Kyo1',i,k ,opnm,Kyo(i,1,k ,opnm),Kyo(i,Ncelly+1,k ,BOTM),&
                   Kyo(i,1,k ,opnm)-Kyo(i,Ncelly+1,k ,BOTM)
           end if
           if((Kyo(i,Ncelly,k ,opnm).ne.Kyo(i,0,k ,TOP))) then
              write(*,*) 'ERR-Kyo2',i,k ,opnm,Kyo(i,Ncelly,k ,opnm),Kyo(i,0,k ,TOP),&
                   opnm,Kyo(i,Ncelly,k ,opnm)-Kyo(i,0,k ,TOP)
           end if
        end do

        do j = 1,Ncelly
           if((Kyo(i,j,1 ,opnm).ne.Kyo(i,j,Ncellz+1 ,DOWN))) then
              write(*,*) 'ERR-Kyo3',i,j ,opnm,Kyo(i,j,1 ,opnm),Kyo(i,j,Ncellz+1 ,DOWN),&
                   Kyo(i,j,1 ,opnm)-Kyo(i,j,Ncellz+1 ,DOWN)
           end if
           if((Kyo(i,j,Ncellz ,opnm).ne.Kyo(i,j,0 ,UP))) then
              write(*,*) 'ERR-Kyo4',i,j ,opnm,Kyo(i,j,Ncellz ,opnm),Kyo(i,j,0 ,UP),&
                   Kyo(i,j,Ncellz ,opnm)-Kyo(i,j,0 ,UP)
           end if
        end do
     end do
  end do
  deallocate(Kyo)

goto 5999
  !allocate(V(Ncellx+1,Ncelly+1,Ncellz+1,5,1:core))
  do cntnum = 0,cntnumend
     write(*,*) '-------------loop------------- :',cntnum,'/',cntnumend
     write(chnum,'(i4.4)') cntnum
     do N20=1,5
        write(n20_ch,'(i2.2)') N20
        allocate(V(-1:Ncellx+2,-1:Ncelly+2,-1:Ncellz+2,N20,0:core-1))
        do opnm=0,core-1
           write(frank,'(i3.3)') opnm
           open(unit=200,file='x'//chnum//n20_ch//frank//'.dat',FORM='UNFORMATTED',iostat=err)
           if(err>0)then
              goto 1990
           end if
           do k=-1,Ncellz+2
              do j=-1,Ncelly+2
                 do i=-1,Ncellx+2
                    write(200) (V(i,j,k,numloop,opnm),numloop=1,N20)
                 end do
              end do
           end do
           close(200)
        end do

        do opnm=0,core-1
           IST = mod(opnm,NSPLTx); KST = opnm/(NSPLTx*NSPLTy); JST = opnm/NSPLTx-NSPLTy*KST
           LEFT = opnm - 1            ; if(IST.eq.0       ) LEFT = opnm + (NSPLTx-1)
           RIGT = opnm + 1            ; if(IST.eq.NSPLTx-1) RIGT = opnm - (NSPLTx-1)
           BOTM = opnm - NSPLTx       ; if(JST.eq.0       ) BOTM = opnm + NSPLTx*(NSPLTy-1)
           TOP  = opnm + NSPLTx       ; if(JST.eq.NSPLTy-1) TOP  = opnm - NSPLTx*(NSPLTy-1)
           DOWN = opnm - NSPLTx*NSPLTy; if(KST.eq.0       ) DOWN = opnm + NSPLTx*NSPLTy*(NSPLTz-1)
           UP   = opnm + NSPLTx*NSPLTy; if(KST.eq.NSPLTz-1) UP   = opnm - NSPLTx*NSPLTy*(NSPLTz-1)
           if((IST.ne.0).and.(IST.ne.NSPLTx-1)) then
              !do k = -1,Ncellz+2
              !   do j = -1,Ncelly+2
              do k = 1,Ncellz
                 do j = 1,Ncelly
                    do i = 1,2
                       do ival=1,N20
                          if((V(i,j,k,ival,opnm).ne.V(Ncellx+i,j,k,ival,LEFT))) then
                             write(*,*) 'ERR1',i,j,k,ival,opnm,V(i,j,k,ival,opnm),V(Ncellx+i,j,k,ival,LEFT),&
                                  V(i,j,k,ival,opnm)-V(Ncellx+i,j,k,ival,LEFT)
                          end if
                          if((V(Ncellx+1-i,j,k,ival,opnm).ne.V(1-i,j,k,ival,RIGT))) then
                             write(*,*) 'ERR2',i,j,k,ival,opnm,V(Ncellx+1-i,j,k,ival,opnm),V(1-i,j,k,ival,RIGT),&
                                  V(Ncellx+1-i,j,k,ival,opnm)-V(1-i,j,k,ival,RIGT)
                          end if
                       end do
                    end do
                 end do
              end do
              !do k=1,Ncellz+1
              !   do j=1,Ncelly+1
              !      do i=1,Ncellx+1
              !         !if((V(i,j,k,N20,opnm).ne.V(i,j,k,N20,LEFT)).and.(V(i,j,k,N20,opnm).ne.V(i,j,k,N20,RIGHT))) then
              !      end do
              !   end do
              !end do
           end if
        end do
1990    continue

        !--------ypart-------
        do opnm=0,core-1
           open(unit=210,file='y'//chnum//n20_ch//frank//'.dat',FORM='UNFORMATTED',iostat=err)
           if(err>0)then
              goto 2990
           end if
           do k=-1,Ncellz+2
              do j=-1,Ncelly+2
                 do i=-1,Ncellx+2
                    write(210) (V(i,j,k,numloop,opnm),numloop=1,N20)
                 end do
              end do
           end do
           close(210)
        end do

        do opnm=0,core-1
           IST = mod(opnm,NSPLTx); KST = opnm/(NSPLTx*NSPLTy); JST = opnm/NSPLTx-NSPLTy*KST
           LEFT = opnm - 1            ; if(IST.eq.0       ) LEFT = opnm + (NSPLTx-1)
           RIGT = opnm + 1            ; if(IST.eq.NSPLTx-1) RIGT = opnm - (NSPLTx-1)
           BOTM = opnm - NSPLTx       ; if(JST.eq.0       ) BOTM = opnm + NSPLTx*(NSPLTy-1)
           TOP  = opnm + NSPLTx       ; if(JST.eq.NSPLTy-1) TOP  = opnm - NSPLTx*(NSPLTy-1)
           DOWN = opnm - NSPLTx*NSPLTy; if(KST.eq.0       ) DOWN = opnm + NSPLTx*NSPLTy*(NSPLTz-1)
           UP   = opnm + NSPLTx*NSPLTy; if(KST.eq.NSPLTz-1) UP   = opnm - NSPLTx*NSPLTy*(NSPLTz-1)
           !        if((JST.ne.0).and.(JST.ne.NSPLTy-1)) then
           !do k = -1,Ncellz+2
           !   do j = -1,Ncelly+2
           do k = 1,Ncellz
              do j = 1,2
                 do i = 1,Ncellx
                    do ival=1,N20
                       if((V(i,j,k,ival,opnm).ne.V(i,Ncelly+j,k,ival,BOTM))) then
                          write(*,*) 'ERR3',i,j,k,ival,opnm,V(i,j,k,ival,opnm),V(i,Ncelly+j,k,ival,BOTM),&
                               V(i,j,k,ival,opnm)-V(i,Ncelly+j,k,ival,BOTM)
                       end if
                       if((V(i,Ncelly+1-j,k,ival,opnm).ne.V(i,1-j,k,ival,TOP))) then
                          write(*,*) 'ERR4',i,j,k,ival,opnm,V(i,Ncelly+1-j,k,ival,opnm),V(i,1-j,k,ival,TOP),&
                               V(i,Ncelly+1-j,k,ival,opnm)-V(i,1-j,k,ival,TOP)
                       end if
                    end do
                 end do
              end do
           end do
           !        end if
        end do
2990    continue


        !--------zpart-------
        do opnm=0,core-1
           open(unit=220,file='z'//chnum//n20_ch//frank//'.dat',FORM='UNFORMATTED',iostat=err)
           if(err>0)then
              goto 3990
           end if
           do k=-1,Ncellz+2
              do j=-1,Ncelly+2
                 do i=-1,Ncellx+2
                    write(220) (V(i,j,k,numloop,opnm),numloop=1,N20)
                 end do
              end do
           end do
           close(220)
        end do

        do opnm=0,core-1
           IST = mod(opnm,NSPLTx); KST = opnm/(NSPLTx*NSPLTy); JST = opnm/NSPLTx-NSPLTy*KST
           LEFT = opnm - 1            ; if(IST.eq.0       ) LEFT = opnm + (NSPLTx-1)
           RIGT = opnm + 1            ; if(IST.eq.NSPLTx-1) RIGT = opnm - (NSPLTx-1)
           BOTM = opnm - NSPLTx       ; if(JST.eq.0       ) BOTM = opnm + NSPLTx*(NSPLTy-1)
           TOP  = opnm + NSPLTx       ; if(JST.eq.NSPLTy-1) TOP  = opnm - NSPLTx*(NSPLTy-1)
           DOWN = opnm - NSPLTx*NSPLTy; if(KST.eq.0       ) DOWN = opnm + NSPLTx*NSPLTy*(NSPLTz-1)
           UP   = opnm + NSPLTx*NSPLTy; if(KST.eq.NSPLTz-1) UP   = opnm - NSPLTx*NSPLTy*(NSPLTz-1)
           !     if((KST.ne.0).and.(KST.ne.NSPLTz-1)) then
           !do k = -1,Ncellz+2
           !   do j = -1,Ncelly+2
           do k = 1,2
              do j = 1,Ncelly
                 do i = 1,Ncellx
                    do ival=1,N20
                       if((V(i,j,k,ival,opnm).ne.V(i,j,Ncellz+k,ival,DOWN))) then
                          write(*,*) 'ERR5',i,j,k,ival,opnm,V(i,j,k,ival,opnm),V(i,j,Ncellz+k,ival,DOWN),&
                               V(i,j,k,ival,opnm)-V(i,j,Ncellz+k,ival,DOWN)
                       end if
                       if((V(i,j,Ncellz+1-k,ival,opnm).ne.V(i,j,1-k,ival,UP))) then
                          write(*,*) 'ERR6',i,j,k,ival,opnm,V(i,j,Ncellz+1-k,ival,opnm),V(i,j,1-k,ival,UP),&
                               V(i,j,Ncellz+1-k,ival,opnm)-V(i,j,1-k,ival,UP)
                       end if
                    end do
                 end do
              end do
           end do
           !     end if
        end do
3990    continue

        deallocate(V)
     end do
  end do

5999 continue

  do cntnum = 0,1
     allocate(PHI(-1:Ncellx+2,-1:Ncelly+2,-1:Ncellz+2,0:core-1))
     write(*,*) '-------------phi------------- :',cntnum,'/',1
     write(chnum,'(i4.4)') cntnum
     do opnm=0,core-1
        write(frank,'(i3.3)') opnm
        open(unit=300,file='ephi'//chnum//frank//'.dat',FORM='UNFORMATTED',iostat=err)
        !allocate(PHI(-1:Ncellx+2,-1:Ncelly+2,-1:Ncellz+2,0:core-1))
        do k=-1,Ncellz+2
           do j=-1,Ncelly+2
              do i=-1,Ncellx+2
                 write(300) PHI(i,j,k,opnm)
              end do
           end do
        end do
        close(300)
     end do
     do opnm=0,core-1
        IST = mod(opnm,NSPLTx); KST = opnm/(NSPLTx*NSPLTy); JST = opnm/NSPLTx-NSPLTy*KST
        LEFT = opnm - 1            ; if(IST.eq.0       ) LEFT = opnm + (NSPLTx-1)
        RIGT = opnm + 1            ; if(IST.eq.NSPLTx-1) RIGT = opnm - (NSPLTx-1)
        BOTM = opnm - NSPLTx       ; if(JST.eq.0       ) BOTM = opnm + NSPLTx*(NSPLTy-1)
        TOP  = opnm + NSPLTx       ; if(JST.eq.NSPLTy-1) TOP  = opnm - NSPLTx*(NSPLTy-1)
        DOWN = opnm - NSPLTx*NSPLTy; if(KST.eq.0       ) DOWN = opnm + NSPLTx*NSPLTy*(NSPLTz-1)
        UP   = opnm + NSPLTx*NSPLTy; if(KST.eq.NSPLTz-1) UP   = opnm - NSPLTx*NSPLTy*(NSPLTz-1)
        if((IST.ne.0).and.(IST.ne.NSPLTx-1)) then
           !do k = -1,Ncellz+2
           !   do j = -1,Ncelly+2
           do k = 1,Ncellz
              do j = 1,Ncelly
                 do i = 1,2
                    !do ival=1,N20
                    if((PHI(i,j,k,opnm).ne.PHI(Ncellx+i,j,k,LEFT))) then
                       write(*,*) 'ERR7',i,j,k,opnm,PHI(i,j,k ,opnm),PHI(Ncellx+i,j,k ,LEFT),&
                            PHI(i,j,k ,opnm)-PHI(Ncellx+i,j,k ,LEFT)
                    end if
                    if((PHI(Ncellx+1-i,j,k ,opnm).ne.PHI(1-i,j,k ,RIGT))) then
                       write(*,*) 'ERR8',i,j,k,opnm,PHI(Ncellx+1-i,j,k ,opnm),PHI(1-i,j,k ,RIGT),&
                            PHI(Ncellx+1-i,j,k ,opnm)-PHI(1-i,j,k ,RIGT)
                    end if
                    !end do
                 end do
              end do
           end do
        end if

        do k = 1,Ncellz
           do j = 1,2
              do i = 1,Ncellx
                 !do ival=1,N20
                 if((PHI(i,j,k ,opnm).ne.PHI(i,Ncelly+j,k ,BOTM))) then
                    write(*,*) 'ERR9',i,j,k ,opnm,PHI(i,j,k ,opnm),PHI(i,Ncelly+j,k ,BOTM),&
                         PHI(i,j,k ,opnm)-PHI(i,Ncelly+j,k ,BOTM)
                 end if
                 if((PHI(i,Ncelly+1-j,k ,opnm).ne.PHI(i,1-j,k ,TOP))) then
                    write(*,*) 'ERR10',i,j,k ,opnm,PHI(i,Ncelly+1-j,k ,opnm),PHI(i,1-j,k ,TOP),&
                         PHI(i,Ncelly+1-j,k ,opnm)-PHI(i,1-j,k ,TOP)
                 end if
                 !end do
              end do
           end do
        end do

        do k = 1,2
           do j = 1,Ncelly
              do i = 1,Ncellx
                 !do ival=1,N20
                 if((PHI(i,j,k ,opnm).ne.PHI(i,j,Ncellz+k ,DOWN))) then
                    write(*,*) 'ERR11',i,j,k ,opnm,PHI(i,j,k ,opnm),PHI(i,j,Ncellz+k ,DOWN),&
                         PHI(i,j,k ,opnm)-PHI(i,j,Ncellz+k ,DOWN)
                 end if
                 if((PHI(i,j,Ncellz+1-k ,opnm).ne.PHI(i,j,1-k ,UP))) then
                    write(*,*) 'ERR12',i,j,k ,opnm,PHI(i,j,Ncellz+1-k ,opnm),PHI(i,j,1-k ,UP),&
                         PHI(i,j,Ncellz+1-k ,opnm)-PHI(i,j,1-k ,UP)
                 end if
                 !end do
              end do
           end do
        end do

     end do
     deallocate(PHI)
  end do



end program
