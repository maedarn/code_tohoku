program main
  implicit none
  integer kp
  integer i,M,t,tstep
  double precision :: ufl(0:1000)=0.0d0,ulw(0:1000)=0.0d0,ug(0:1000)=0.0d0,uan(0:1000)=0.0d0 !初期化
  double precision :: ufl2(0:1000)=0.0d0,ulw2(0:1000)=0.0d0,ug2(0:1000)=0.0d0 !初期化
  double precision pi,M2,Min,kappa,x(0:1000),L,dx,kkappa,timecount,kappasave
  double precision ct,ctoutput,tsave,tsavesub
  integer ctstep

  !pi=acos(-1.0d0)
  pi=3.14159265358979d0
  !*******************parameter*****************
  M=20 !wave length
  kappa=20.d0/25.d0 !courant condition
  kp=2021 !filename(kappa)
  L=1.0d0 !長さ
  ctoutput=L*0.5d0  !final step どこまで進めるか
  ctstep=5 !output (in final step)
  !*********************************************

  !*******************space*********************
  dx=L/1000.0d0
  x(0)=0.0d0
  do i=1,1000
     x(i)=x(i-1)+dx
  end do
  !*********************************************

  !*******************Initial*******************
  M2=dble(M)/2.0d0
  Min=1.0d0/dble(M)
  !do i=30-M/2,30+M/2
  do i=2,51
     ufl(i)= 10.d0 !* dsin(pi*dble(i-30+M2)*Min)
     ulw(i)=10.d0 !* dsin(pi*dble(i-30+M2)*Min)
     ug(i)=10.d0 !* dsin(pi*dble(i-30+M2)*Min)
     uan(i)=10.d0 !* dsin(pi*dble(i-30+M2)*Min)
     !ufl(i)=10.d0 * dsin(pi*dble(i-30+M2)*Min)
     !ulw(i)=10.d0 * dsin(pi*dble(i-30+M2)*Min)
     !ug(i)=10.d0 * dsin(pi*dble(i-30+M2)*Min)
     !uan(i)=10.d0 * dsin(pi*dble(i-30+M2)*Min)
  end do
  call save(ufl,ulw,ug,uan,x,kp,M)
  !*********************************************


  !*********************BC**********************
  !初期化で計算済み(Dirichlet)
  !*********************************************


  !********************evolve*******************
  kappasave=kappa
  kkappa=kappa*kappa
  ct=kappa*dx !c*dt
  timecount=0.0d0
  tsave=ctoutput/dble(ctstep)
  tsavesub=tsave
  tstep=int(ctoutput/ct)+ctstep+10 !coutput/(c*dt)
  !tsave=tstep/ctstep !save step
  !tsavesub=tsave
  write(*,*) tstep,tsave

  do t=1,tstep
     if(tsave<timecount) then
        kappa=(timecount-tsave)/dx
        write(*,*) t,tstep,kappa,kappasave
        call evolve(ufl,ulw,ug,ufl2,ulw2,ug2,kappa)
        kappa=kappasave
        kkappa=kappa*kappa
        uan(:)=0.0d0
        do i=30-M/2+int(650.d0*tsave/L+1.5),30+M/2+int(650.d0*tsave/L+1.5) !analytic
           uan(i)=10.d0 * dsin(pi*dble(i-int(650.d0*tsave/L+1.5)-30+M2)*Min)
        end do
        call save(ufl,ulw,ug,uan,x,kp,M)
        timecount=tsave
        tsave=tsave+tsavesub
     else
        call evolve(ufl,ulw,ug,ufl2,ulw2,ug2,kappa)
         timecount=timecount+ct
     end if
  end do
  !*********************************************
  !call save(ufl,ulw,ug,x,kp)
end program main



subroutine save(ufl,ulw,ug,uan,x,kp,M)
  character(28)  filenamefl,filenamelw,filenamegd,filenamean
  integer i,kp,M
  integer :: tnum=0
  double precision ufl(0:1000),ulw(0:1000),ug(0:1000),uan(0:1000),x(0:1000)

  write (filenamefl,'("hypbfl",i4.4,i2.2,i1.1,".dat")')kp, tnum,M/10
  write (filenamelw,'("hypblw",i4.4,i2.2,i1.1,".dat")')kp, tnum,M/10
  write (filenamegd,'("hypbgd",i4.4,i2.2,i1.1,".dat")')kp, tnum,M/10
  write (filenamean,'("hypban",i4.4,i2.2,i1.1,".dat")')kp, tnum,M/10
100 format(E19.10e3,E19.10e3)

  open (10, file=filenamefl)
  do i=0,1000
     write(10,100) x(i),ufl(i)
  end do
  close(10)
  open (11, file=filenamelw)
  do i=0,1000
     write(11,100) x(i),ulw(i)
  end do
  close(11)
  open (12, file=filenamegd)
  do i=0,1000
     write(12,100) x(i),ug(i)
  end do
  close(12)
  open (13, file=filenamean)
  do i=0,1000
     write(13,100) x(i),uan(i)
  end do
  close(13)

  tnum=tnum+1
end subroutine save

subroutine evolve(ufl,ulw,ug,ufl2,ulw2,ug2,kappa)
  double precision ufl(0:1000),ulw(0:1000),ug(0:1000)
  double precision ufl2(0:1000),ulw2(0:1000),ug2(0:1000)
  double precision kappa,kkappa
  kkappa=kappa*kappa
  !===============FL================
  do i=1,1000-1
     ufl2(i)=0.5d0*(1-kappa)*ufl(i+1)+0.5d0*(1+kappa)*ufl(i-1)
  end do
  do i=1,1000-1
     ufl(i)=ufl2(i)
  end do
  !===============LW================
  do i=1,1000-1
     ulw2(i)=0.5d0*(kkappa-kappa)*ulw(i+1)+(1.0d0-kkappa)*ulw(i)+0.5d0*(kkappa+kappa)*ulw(i-1)
  end do
  do i=1,1000-1
     ulw(i)=ulw2(i)
  end do
  !===============GD================
  do i=1,1000-1
     ug2(i)=(1-kappa)*ug(i)+kappa*ug(i-1)
  end do
  do i=1,1000-1
     ug(i)=ug2(i)
  end do
end subroutine evolve

subroutine evolve2(ufl,ulw,ug,ufl2,ulw2,ug2,kappa)
  double precision ufl(0:1000),ulw(0:1000),ug(0:1000)
  double precision ufl2(0:1000),ulw2(0:1000),ug2(0:1000)
  double precision kappa,kkappa,s0,s650
  kkappa=kappa*kappa
  !===============FL================
  do i=1,1000-1
     ufl2(i)=0.5d0*(1-kappa)*ufl(i+1)+0.5d0*(1+kappa)*ufl(i-1)
  end do
  do i=1,1000-1
     ufl(i)=ufl2(i)
  end do
  s0=ufl(1)
  s650=ufl(1000-1)
  ufl(0)=s650
  ufl(1000)=s0
  !===============LW================
  do i=1,1000-1
     ulw2(i)=0.5d0*(kkappa-kappa)*ulw(i+1)+(1.0d0-kkappa)*ulw(i)+0.5d0*(kkappa+kappa)*ulw(i-1)
  end do
  do i=1,1000-1
     ulw(i)=ulw2(i)
  end do
  s0=ulw(1)
  s650=ulw(1000-1)
  ulw(0)=s650
  ulw(1000)=s0
  !===============GD================
  do i=1,1000-1
     ug2(i)=(1-kappa)*ug(i)+kappa*ug(i-1)
  end do
  do i=1,1000-1
     ug(i)=ug2(i)
  end do
  s0=ug(1)
  s650=ug(1000-1)
  ug(0)=s650
  ug(1000)=s0
end subroutine evolve2
