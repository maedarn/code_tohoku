module combar
!implicit none
integer*8 :: prc=0,prcn=0
end module combar


program main
    use combar
!    implicit none
    integer nx,ny,nz,Nlx,Nly,Nlz,Ncellx,Ncelly,Ncellz
    integer ncx,ncy,ncz,ngrid,NGL,NGcr,ncycle,NPRE,NPOST
    

ncycle=2
NPRE=5
NPOST=5

!nx=64
!ny=64
!nz=64
!Ncellx=64
!Ncelly=64
!Ncellz=64
!Nlx=8
!Nly=8
!Nlz=8
!Nlx=int((1.d5)**(1.d0/3.d0))
!Nly=int((1.d5)**(1.d0/3.d0))
!Nlz=int((1.d5)**(1.d0/3.d0))
!nx=int(1.d4/Nlx)
!ny=int(1.d4/Nly)
!nz=int(1.d4/Nlz)
!Ncellx=int(1.d4/Nlx)
!Ncelly=int(1.d4/Nly)
!Ncellz=int(1.d4/Nlz)
Nlx=int((2.d4)**(1.d0/3.d0))
Nly=int((2.d4)**(1.d0/3.d0))
Nlz=int((2.d4)**(1.d0/3.d0))
nx=int(5.d3/Nlx)
ny=int(5.d3/Nly)
nz=int(5.d3/Nlz)
Ncellx=int(5.d3/Nlx)
Ncelly=int(5.d3/Nly)
Ncellz=int(5.d3/Nlz)


NGL = Ncellx*Nlx
NGL = int(dlog(dble(NGL))/dlog(2.d0)+1.d-3)

NGcr = Nlx
NGcr = int(dlog(dble(NGcr))/dlog(2.d0)+1.d-3) + 1


prc=nx*ny*nz+prc


ncx=(Ncellx+1)/2+1; ncy=(Ncelly+1)/2+1; ncz=(Ncellz+1)/2+1
ngrid=NGL-1

call rstrctMPI(ncx,ncy,ncz)

do while(ngrid.ne.NGcr)
  ncx=ncx/2+1; ncy=ncy/2+1; ncz=ncz/2+1
  ngrid=ngrid-1
  call rstrctMPI(ncx,ncy,ncz)
end do
ncxcr=ncx; ncycr=ncy; nczcr=ncz

ncx=2**NGcr+1;ncy=ncx;ncz=ncx
call collect(ncx,ncy,ncz,ncxcr,ncycr,nczcr)

do while(ngrid.ne.1)
  ncx=ncx/2+1; ncy=ncy/2+1; ncz=ncz/2+1
  ngrid=ngrid-1
  call rstrct(ncx,ncy,ncz)
end do

call slvsmlb() !BC set is necessary

!Here nc=3
ngrid = NGL
do j=2,ngrid
  
  IF(j.le.NGcr) THEN !*** generate candidate sol. from j-1 to j (upward) ***
    ncx=ncx*2-1; ncy=ncy*2-1; ncz=ncz*2-1
    call interp(ncx,ncy,ncz)  !BC set is necessary
    call copy(ncx,ncy,ncz)
    if(j.eq.NGcr) then
      ncx=ncxcr; ncy=ncycr; ncz=nczcr
      call divide(ncx,ncy,ncz,2**NGcr+1,2**NGcr+1,2**NGcr+1 )
      call copyMPI(ncx,ncy,ncz)
    end if
  ELSE
    ncx=ncx*2-1; ncy=ncy*2-1; ncz=ncz*2-1
    call interpMPI(ncx,ncy,ncz)  !BC set is necessary
    if(j.ne.ngrid) call copyMPI(ncx,ncy,ncz)
  END IF

  do jcycle=1,ncycle !V-cycle

    nfx=ncx; nfy=ncy; nfz=ncz
    do jj=j,2,-1          !*** DOWNWARD *****************************
                          !    phi + rhs --> res -->      (level N  )
                          !                          rhs  (level N-1)
      IF(jj.lt.NGcr) THEN !*** generate residual from jj to jj-1 ****
        do jpre=1,NPRE
          mode=2; if((jj.ne.j).and.(jpre.eq.1)) mode=1
          call relax(nfx,nfy,nfz)
        end do
        call resid(nfx,nfy,nfz)
        nfx=nfx/2+1; nfy=nfy/2+1; nfz=nfz/2+1
        call rstrct(nfx,nfy,nfz)  !fill0 at BC below this subroutine is necessary
        call  fill0(nfx,nfy,nfz)
      ELSE
        NPRE1 = NPRE; if(j.ge.NGL) NPRE1 = 2
        do jpre=1,NPRE1
          mode=2; if((jj.ne.j).and.(jpre.eq.1)) mode=1
          call relaxMPI(nfx,nfy,nfz)
        end do
        call residMPI(nfx,nfy,nfz)
        if(jj.eq.NGcr) then
          nfx=2**NGcr+1;nfy=nfx;nfz=nfx
          call collect(nfx,nfy,nfz,ncxcr,ncycr,nczcr ) !necessary at upward loop below
          call collect(nfx,nfy,nfz,ncxcr,ncycr,nczcr )
          nfx=nfx/2+1; nfy=nfy/2+1; nfz=nfz/2+1
          call rstrct(nfx,nfy,nfz)  !fill0 at BC below this subroutine is necessary
          call  fill0(nfx,nfy,nfz)
        else
          nfx=nfx/2+1; nfy=nfy/2+1; nfz=nfz/2+1
          call rstrctMPI(nfx,nfy,nfz)  !fill0 at BC below this subroutine is necessary
          call  fill0MPI(nfx,nfy,nfz)
        end if
      END IF
    end do

    call slvsml()  !BC set is unnecessary

    nfx=3; nfy=3; nfz=3
    do jj=2,j             !*** UPWARD **********************************
                          !    phi --> phi +rhs --> phi      (level N  )
                          !  + phi                           (level N-1)
      IF(jj.le.NGcr) THEN !*** generate new solution from jj-1 to jj ***
        nfx=2*nfx-1; nfy=2*nfy-1; nfz=2*nfz-1
        call addint(nfx,nfy,nfz)  !BC set is unnecessary
        if(jj.eq.NGcr) then
          nfx=ncxcr; nfy=ncycr; nfz=nczcr
          call divide(nfx,nfy,nfz,2**NGcr+1,2**NGcr+1,2**NGcr+1 )
          do jpost=1,NPOST
            call relaxMPI(nfx,nfy,nfz)
          end do
        else
          do jpost=1,NPOST
            call relax(nfx,nfy,nfz)
          end do
        end if
      ELSE
        nfx=2*nfx-1; nfy=2*nfy-1; nfz=2*nfz-1
        call addintMPI(nfx,nfy,nfz)  !BC set is unnecessary
        NPOST1 = NPOST; if(j.ge.NGL) NPOST1 = 2
        do jpost=1,NPOST1
          call relaxMPI(nfx,nfy,nfz)
        end do
      END IF
    end do

  end do
end do

!ncx = Ncellx+1
prc=prc+nx*ny*nz

write(*,*) prc,prcn,dble(prc+prcn)/dble(prcn+prc/Nlx/Nly/Nlz), &
    1.d0/dble(prc+prcn)*dble(prcn), 1-1.d0/dble(prc+prcn)*dble(prcn), &  !,1.d0/dble(prc+prcn)*dble(prc/Nlx/Nly/Nlz)
    dble(prc+prcn)/dble(prcn+prc/Nlx/Nly/Nlz)/dble(Nlx*Nly*Nlz)*100!,1.d0/dble(prc+prcn)*dble(prc/Nlx/Nly/Nlz)
end program main


subroutine collect(nx1,ny1,nz1,nx2,ny2,nz2)
USE combar
integer nx1,ny1,nz1,nx2,ny2,nz2

prc=prc+nz2*ny2*nx2 * 2
end subroutine collect

SUBROUTINE divide(nx2,ny2,nz2,nx1,ny1,nz1)
USE combar
integer nx1,ny1,nz1,nx2,ny2,nz2

prc=prc+nz2*ny2*nx2
END SUBROUTINE divide

SUBROUTINE rstrctMPI(nx,ny,nz)
USE combar
integer nx,ny,nz

prc=prc+nz*ny*nx
prcn=prcn+nz*ny*1 !6

END SUBROUTINE rstrctMPI

SUBROUTINE rstrct(nx,ny,nz)
USE combar
integer nx,ny,nz

prc=prc+nz*ny*nx
prc=prc+nz*ny*3

prcn=prcn+nz*ny*1

END SUBROUTINE rstrct

SUBROUTINE interpMPI(nx,ny,nz)
USE combar
integer nx,ny,nz

prc=prc+nz*ny*nx/2/2/2
prc=prc+nz*ny*nx/2/2/2
prc=prc+nz*ny*nx/2/2
prc=prc+nz*ny*nx/2

prcn=prcn+nz*ny*1
END SUBROUTINE interpMPI

SUBROUTINE interp(nx,ny,nz)
USE combar
integer nx,ny,nz

prc=prc+nz*ny*nx/2/2/2
prc=prc+nz*ny*nx/2/2/2
prc=prc+nz*ny*nx/2/2
prc=prc+nz*ny*nx/2

prcn=prcn+nz*ny*1
END SUBROUTINE interp

SUBROUTINE addintMPI(nx,ny,nz)
USE combar
integer nx,ny,nz

call interpMPI(nx,ny,nz)

prc=prc+nz*ny*nx
END SUBROUTINE addintMPI

SUBROUTINE addint(nx,ny,nz)
USE combar
integer nx,ny,nz

call interp(nx,ny,nz)

prc=prc+nz*ny*nx
END SUBROUTINE addint

SUBROUTINE slvsml()
USE combar

prc=prc+ 9*16
end SUBROUTINE slvsml

SUBROUTINE slvsmlb()
USE combar

prc=prc+ 9*60
end SUBROUTINE slvsmlb

SUBROUTINE relaxMPI(nx,ny,nz)
USE combar
integer nx,ny,nz

prc=prc+nz*ny*nx*12
END SUBROUTINE relaxMPI

SUBROUTINE relax(nx,ny,nz)
USE combar
integer nx,ny,nz

prc=prc+nz*ny*nx*6
END SUBROUTINE relax


SUBROUTINE residMPI(nx,ny,nz)
USE combar
integer nx,ny,nz

prc=prc+nx*ny*nz*10
prcn=prcn+nx*ny/2/2
END SUBROUTINE residMPI


SUBROUTINE resid(nx,ny,nz)
USE combar
integer nx,ny,nz

prc=prc+nx*ny*nz*10
!prcn=prcn+ndx*ndy/2/2
END SUBROUTINE resid


SUBROUTINE copy(nx,ny,nz)
integer nx,ny,nz

prc=prc+nx*ny*nz
END SUBROUTINE copy

SUBROUTINE copyMPI(nx,ny,nz)
USE combar
integer nx,ny,nz

prc=prc+nx*ny*nz
END SUBROUTINE copyMPI

SUBROUTINE fill0(nx,ny,nz)
USE combar
integer nx,ny,nz

prc=prc+nx*ny*6
END SUBROUTINE fill0

SUBROUTINE fill0MPI(nx,ny,nz)
USE combar
integer nx,ny,nz

prc=prc+nx*ny*9
END SUBROUTINE fill0MPI

