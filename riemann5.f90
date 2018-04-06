program main
  implicit none
  character filename*128

  integer i , cells
  real(8) gamma , g1 , g2 , g3 , g4 , g5 , g6 , g7 , g8
  real(8) dl , ul , pl , cl , dr , ur , pr , cr
  real(8) diaph , domlen , ds , dx , pm , mpa , ps , s , timeout , um , us
  real(8) xpos

  integer  nriter
  real(8) change , fl , fld , fr , frd  , pold , pstart
  real(8) tolpre  , udiff

  real(8) cup , gel , ger  , pmax , pmin , ppv , pq
  real(8) ptl , ptr , qmax , quser , um2

  real(8) ak , bk  , prat , qrt


  real(8) c , cml , cmr  , pml , pmr
  real(8) shl , shr , sl , sr , stl , str

  integer dt , tstep
  real(8) a
  !initial
  domlen = 1.0d0
  diaph = 0.50d0
  cells = 1024
  gamma = 5.0d0 / 3.0d0

  timeout =  0.00000000000000000D+00
 !はじめのout put time











  
  dx = domlen / real(cells)
  a =  0.18730704064812873D+00 !時間刻み
  tstep = 1 !時間ステップ数
  nriter = 1000 !反復(ニュートンラフソン)
  tolpre = 1.0d-6 !ニュートンラフソン基準

  !allocate( dl(cells/2) , ul(cells/2) , pl(cells/2) , cl(cells/2) , dr(cells/2) , ur(cells/2) , pr(cells/2) , cr(cells/2))

  dl = 1.0d0
  ul = 0.0d0
  pl = 1.0d0
  dr = 0.1250d0
  ur = 0.0d0
  pr = 0.10d0

  mpa = 1.0d0 !圧力を見やすさのため規格化

  write (filename, '("riemanndata", i4.4, ".dat")') 0
  open (17, file=filename, status='replace')

  do  i = 1 , cells  !初期条件の記述
     if(i <= cells/2) then
        ds = dl
        us = ul
        ps = pl
     else
        ds = dr
        us = ur
        ps = pr
     endif

     xpos = (real(i) - 0.5d0) * dx
     write(17,*) xpos , ds , us , ps/mpa , ps/ds/g8/mpa

  end do

  close(17)

  !gamma

  g1 = (gamma - 1.0d0)/(2.0d0 * gamma)
  g2 = (gamma + 1.0d0)/(2.0d0 * gamma)
  g3 = 2.0d0 * gamma / (gamma - 1.0d0)
  g4 = 2.0d0 / (gamma - 1.0d0)
  g5 = 2.0d0 / (gamma + 1.0d0)
  g6 = (gamma - 1.0d0) / (gamma + 1.0d0)
  g7 = (gamma - 1.0d0) / 2.0d0
  g8 = gamma - 1.0d0

  !sound speed

  cl = sqrt(gamma * pl / dl)
  cr = sqrt(gamma * pr / dr)

  !真空になるときストップ

  if(g4 * (cl + cr) <= (ur - ul)) then

     write(*,*)
     write(*,*) '***Vacuum is generated by data***'
     write(*,*) '***Program stopped***'
     write(*,*)

     stop
  end if

  !clicurate p* for initial data (to start iteration : approxmation)

  quser = 2.0d0

  cup = 0.25d0 * (dl + dr) * (cl + cr)
  ppv = 0.5d0 * (pl + pr) + 0.5d0 * (ul - ur) * cup
  ppv = DMAX1(0.0d0 , ppv)   !https://jp.xlsoft.com/documents/intel/cvf/vf-html/az/az08_07.htm
  pmin = DMIN1(pl , pr) !小さい方 D : real(8)
  pmax = DMAX1(pl , pr)
  qmax = pmax/pmin

  if(qmax <= quser .AND. (pmin <= ppv .AND. ppv <= pmax)) then

     pstart = ppv

  else

     if(ppv <= pmin) then

        pq = (pl/pr) ** g1
        um = (pq*ul/cl + ur/cr + g4*(pq - 1.0d0))/(pq/cl + 1.0d0/cr)
        ptl = 1.0d0 + g7*(ul - um)/cl
        ptr = 1.0d0 + g7*(um - ur)/cr
        pstart = 0.5d0 * (pl*ptl**g3 + pr*ptr**g3)

     else

        gel = sqrt((g5/dl)/(g6*pl + ppv))
        ger = sqrt((g5/dr)/(g6*pr + ppv))
        pstart = (gel*pl + ger*pr - (ur - ul))/(gel + ger)

     endif

  endif

  pold = pstart
  udiff = ur - ul

  write(6,*)'----------------------------------'
  write(6,*)'iteration number   change'
  write(6,*)'----------------------------------'


  do  i = 1 , nriter

     !微分を出す ニュートンラフソン法
     if(pold <= pl) then

        prat = pold/pl
        fl = g4 * cl * (prat**g1 -1.0d0)
        fld = (1.0d0 / (dl * cl))*prat**(-g2)

     else

        ak = g5 / dl
        bk = g6 * pl
        qrt = sqrt(ak/(bk + pold))
        fl = (pold - pl)*qrt
        fld = (1.0d0 - 0.5d0*(pold - pl)/(bk + pold))

     end if

     if(pold <= pr) then

        prat = pold/pr
        fr = g4 * cr * (prat**g1 -1.0d0)
        frd = (1.0d0 / (dr * cr))*prat**(-g2)

     else

        ak = g5 / dr
        bk = g6 * pr
        qrt = sqrt(ak/(bk + pold))
        fr = (pold - pr)*qrt
        frd = (1.0d0 - 0.5d0*(pold - pr)/(bk + pold))

     end if


     pm = pold - (fl + fr + udiff)/(fld + frd)
     change = 2.0d0 * abs((pm - pold)/(pm + pold))
     write(*,*) i , change

     if(change <= tolpre) exit
     if(pm <= 0.0) pm = tolpre
     pold = pm

  end do


  write(6,*) 'Divergence in Newton-Raphson iteration'



  um = 0.5d0 * (ul + ur + fr - fl)

  write(*,*)'----------------------------------'
  write(*,*)'   pressure      velocity  '
  write(*,*)'----------------------------------'
  write(*,*) pm/mpa , um
  write(*,*)'----------------------------------'

  !30   format(5x , i5 , 15x , f12.7)
  !40   format(2(f14.6 , 5x))


  do dt = 1 , tstep

     timeout = timeout + a * dt
     !open(unit=2 , file = 'exact.dat' , STATUS = 'unknown') !unknown:ファイルがあればold、無ければnew
     write (filename, '("riemanndata", i4.4, ".dat")') dt
     open (17, file=filename, status='replace')

  do  i = 1 , cells

     xpos = (real(i) - 0.5d0) * dx

     s = (xpos - diaph) / timeout

     if(s <= um) then

        if(pm <= pl) then

           shl = ul - cl

           if(s <= shl) then

              ds = dl
              us = ul
              ps = pl

           else

              cml = cl*(pm/pl)**g1
              stl = um - cml

              if(s > stl) then

                 ds = dl * (pm/pl)**(1.0d0/gamma)
                 us = um
                 ps = pm
              else

                 us = g5*(cl + g7*ul + s)
                 c = g5*(cl + g7*(ul - s))
                 ds = dl*(c/cl)**g4
                 ps = pl*(c/cl)**g3

              endif
           endif
        else

           pml = pm/pl
           sl = ul - cl*sqrt(g2*pml + g1)

           if(s <= sl) then

              ds = dl
              us = ul
              ps = pl

           else

              ds = dl * (pml + g6) / (1.0d0 + pml*g6)
              us = um
              ps = pm

           endif
        endif
     else

        if(pm > pr) then

           pmr = pm/pr
           sr = ur + cr*sqrt(g2*pmr + g1)

           if(s >= sr) then

              ds = dr
              us = ur
              ps = pr

           else

              ds = dr*(pmr + g6)/(pmr*g6 + 1.0d0)
              us = um
              ps = pm

           endif
        else
           shr = ur + cr
           if(s >= shr) then

              ds = dr
              us = ur
              ps = pr

           else

              cmr = cr*(pm/pr)**g1
              str = um + cmr

              if(s <= str) then

                 ds = dr * (pm/pr)**(1.0d0/gamma)
                 us = um
                 ps = pm
              else

                 us = g5*(-cr + g7*ur + s)
                 c = g5*(cr - g7*(ur - s))
                 ds = dr*(c/cr)**g4
                 ps = pr*(c/cr)**g3

              endif
           endif
        endif
     endif

     write(17,*) xpos , ds , us , ps/mpa , ps/ds/g8/mpa , ps/(ds**(gamma)) , 0.5d0 *  (us**2)

  end do
end do

  close(17)
  !20      format(5(f14.6 , 2x))
end program main
