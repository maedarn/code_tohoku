program main
  integer NRANK,NPE,NSPLTx,NSPLTy,NSPLTz,Ncellx
  NRANK=63
  NPE=8
  NSPLTx=2
  NSPLTy=2
  NSPLTz=2
  Ncellx=32

  do NN=0,NRANK

  KST = NN/(NSPLTx*NSPLTy); JST = NN/NSPLTx-NSPLTy*KST

  do Nlp = 1,NSPLTy*NSPLTz-1

     isend = NN + NSPLTx*Nlp; if(isend.ge.NPE) isend = isend - NPE
     KSs = isend/(NSPLTx*NSPLTy); JSs = isend/NSPLTx-NSPLTy*KSs
     irecv = NN - NSPLTx*Nlp; if(irecv.lt.0  ) irecv = irecv + NPE
     KSr = irecv/(NSPLTx*NSPLTy); JSr = irecv/NSPLTx-NSPLTy*KSr

     Nis = JSs + NSPLTy*KSs
     kls = Nis + 1
     Nir = JST + NSPLTy*KST
     klr = Nir + 1

     !***************fordebug*****************
     !CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
     !***************fordebug*****************

     if(kls.gt.Ncellx) then
        isend = -1
        kls = Ncellx+1
     end if
     if(klr.gt.Ncellx) then
        irecv = -1
        klr = Ncellx+1
     end if
     write(*,*)isend,irecv,kls,klr,Nis,Nir,JSr,KSr
  end do
  write(*,*)
end do

end program main
