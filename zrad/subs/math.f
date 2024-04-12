      SUBROUTINE linear(xa,ya,m,x,y)
      IMPLICIT REAL*8(a-h,o-z)
      DIMENSION xa(m),ya(m)
      do 11 i=1,m
         if(x-xa(i).le.0.d0) then
            ms=i
            go to 12
         endif
 11   continue
      ms=m
 12   continue
      if(ms.eq.1) ms=2
      y1=ya(ms-1)
      y2=ya(ms)
      t=(x-xa(ms-1))/(xa(ms)-xa(ms-1))
      y=(1.d0-t)*y1+t*y2
      return
      END


      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
c      PARAMETER (NMAX=50,eps=1.d-13)
      PARAMETER (NMAX=50,eps=0.d0)
      DOUBLE PRECISION a(np,np),b(np,mp),b_max(np,mp)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      DOUBLE PRECISION big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue

      do ll=1,np
         do l=1,mp
            b_max(ll,l)=0.d0
         enddo
      enddo

      do 22 i=1,n
        big=0.d0
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                print *,'singular matrix in gaussj'
                go to 25 
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.d0) then
           print *, 'singular matrix in gaussj'
           go to 25
        endif
        pivinv=1.d0/a(icol,icol)
        a(icol,icol)=1.d0
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue 

        do 21 ll=1,n
           if(ll.ne.icol)then
              dum=a(ll,icol)
              a(ll,icol)=0.d0
              do 18 l=1,n
                 a(ll,l)=a(ll,l)-a(icol,l)*dum                 
                 if(a(icol,l)*dum+a(ll,l).ne.a(ll,l)) then
                    if(dabs(a(ll,l)/(a(icol,l)*dum)).lt.eps) then
                       a(ll,l)=0.d0
                    endif
                 endif
 18           continue
              do 19 l=1,m
                 b(ll,l)=b(ll,l)-b(icol,l)*dum
                 if(b(icol,l)*dum+b(ll,l).ne.b(ll,l)) then
                    if(dabs(b(ll,l)/(b(icol,l)*dum)).lt.eps) then
                       b(ll,l)=0.d0
                    endif
                 endif
                 b_max(ll,l)=max(b_max(ll,l),dabs(b(icol,l)*dum)) 
 19           continue
           endif
 21     continue
 22   continue

      do ll=1,np
         do l=1,mp
            if(b_max(ll,l).ne.0.d0) then
               if(dabs(b(ll,l)/b_max(ll,l)).le.eps) then
                  b(ll,l)=0.d0
               endif
            endif
         enddo
      enddo


      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
25    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<?4210(93Y"+91.d0

      SUBROUTINE bilinear(x1a,x2a,ya,m,n,x1,x2,y)
      IMPLICIT REAL*8(a-h,o-z)
      DIMENSION x1a(m),x2a(n),ya(m,n)
      do 11 i=1,m
         if(x1-x1a(i).le.0.d0) then
            ms=i
            go to 12  
         endif
 11   continue 
      ms=m
 12   continue
      do 13 i=1,n
         if(x2-x2a(i).le.0.d0) then
            ns=i
            go to 14  
         endif
 13   continue
      ns=n
 14   continue
      if(ms.eq.1) ms=2
      if(ns.eq.1) ns=2
      y1=ya(ms-1,ns-1)
      y2=ya(ms,ns-1)
      y3=ya(ms,ns)
      y4=ya(ms-1,ns)
      t=(x1-x1a(ms-1))/(x1a(ms)-x1a(ms-1))
      u=(x2-x2a(ns-1))/(x2a(ns)-x2a(ns-1))
      y=(1.d0-t)*(1.d0-u)*y1+t*(1.d0-u)*y2+t*u*y3+(1.d0-t)*u*y4
      return
      END
