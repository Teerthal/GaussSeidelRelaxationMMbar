c Here we define our differencing scheme (some differences are higher
c order than others -- better to make everything the same order).
c Note: this routine is for relaxation codes. E.g. derivatives of 
c fd(n,i,j,k) ($\dot f$ in evolution codes) are not calculated.
c
      subroutine derivatives(f,i,j,k,dfdx,dfdy,dfdz,
     1            d2fdx2,d2fdy2,d2fdz2)
      implicit none
      include 'parameters.inc'
      integer i,j,k,n,is
      real*8 f,fd
      real*8 dfdx,dfdy,dfdz,fs,d2fdx2,d2fdy2,d2fdz2
c
      dimension f(nf,-latx:latx,-laty:laty,-latz:latz)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)
c
      do 35 n=1,nf
c
c x-derivatives:
c
c 6th order 1st- and 2nd- derivatives:
       if(abs(i).le.latx-3) then
c central derivatives in lattice interior:
        dfdx(n)=(f(n,i+3,j,k)-9.*f(n,i+2,j,k)+45.*f(n,i+1,j,k)
     1   -45.*f(n,i-1,j,k)+9.*f(n,i-2,j,k)-f(n,i-3,j,k))/(60.*dx)
        d2fdx2(n)=(f(n,i+3,j,k)/90.-3.*f(n,i+2,j,k)/20.
     1   +3.*f(n,i+1,j,k)/2.-49.*f(n,i,j,k)/18.+3.*f(n,i-1,j,k)/2.
     2   -3.*f(n,i-2,j,k)/20.+f(n,i-3,j,k)/90.)/(dx**2)
       endif
c
c one-sided derivatives near boundaries:
        is=+1.
        if(i.lt.0) is=-1
c
        if(abs(i).eq.latx-2) then
        dfdx(n)=-float(is)*(f(n,i+is*2,j,k)/30.-2.*f(n,i+is*1,j,k)/5.
     1   -7.*f(n,i,j,k)/12.+4.*f(n,i-is*1,j,k)/3.-f(n,i-is*2,j,k)/2.
     2   +2.*f(n,i-is*3,j,k)/15.-f(n,i-is*4,j,k)/60.)/dx 
        d2fdx2(n)=(-13.*f(n,i+is*2,j,k)/180.+19.*f(n,i+is*1,j,k)/15.
     1   -7.*f(n,i,j,k)/3.+10.*f(n,i-is*1,j,k)/9.
     2   +f(n,i-is*2,j,k)/12.-f(n,i-is*3,j,k)/15.
     3   +f(n,i-is*4,j,k)/90.)/dx**2
        endif
c
        if(abs(i).eq.latx-1) then
        dfdx(n)=-float(is)*(-f(n,i+is*1,j,k)/6.-77.*f(n,i,j,k)/60.
     1   +5.*f(n,i-is*1,j,k)/2.-5.*f(n,i-is*2,j,k)/3.
     2   +5.*f(n,i-is*3,j,k)/6.-f(n,i-is*4,j,k)/4.
     3   +f(n,i-is*5,j,k)/30.)/dx 
        d2fdx2(n)=(137.*f(n,i+is*1,j,k)/180.-49.*f(n,i,j,k)/60.
     1   -17.*f(n,i-is*1,j,k)/12.+47.*f(n,i-is*2,j,k)/18.
     2   -19.*f(n,i-is*3,j,k)/12.+31.*f(n,i-is*4,j,k)/60.
     3   -13.*f(n,i-is*5,j,k)/180.)/dx**2
        endif
c
c        if(abs(i).eq.latx) then
cc 6th order:
c        dfdx(n)=-float(is)*(-49.*f(n,i,j,k)/20.+6.*f(n,i-is*1,j,k)
c     1   -15.*f(n,i-is*2,j,k)/2.+20.*f(n,i-is*3,j,k)/3.
c     2   -15.*f(n,i-is*4,j,k)/4.+6.*f(n,i-is*5,j,k)/5.
c     3   -f(n,i-is*6,j,k)/6.)/dx
c        dfddx(n)=-float(is)*(-49.*fd(n,i,j,k)/20.+6.*fd(n,i-is*1,j,k)
c     1   -15.*fd(n,i-is*2,j,k)/2.+20.*fd(n,i-is*3,j,k)/3.
c     2   -15.*fd(n,i-is*4,j,k)/4.+6.*fd(n,i-is*5,j,k)/5.
c     3   -fd(n,i-is*6,j,k)/6.)/dx
c        d2fdx2(n)=(203.*f(n,i,j,k)/45.-87.*f(n,i-is*1,j,k)/5.
c     1   +117.*f(n,i-is*2,j,k)/4.-254.*f(n,i-is*3,j,k)/9.
c     2   +33.*f(n,i-is*4,j,k)/2.-27.*f(n,i-is*5,j,k)/5.
c     3   +137.*f(n,i-is*6,j,k)/180.)/dx**2
c        endif
c
        if(abs(i).eq.latx) then
c 8th order:
        dfdx(n)=-float(is)*(-761.*f(n,i,j,k)/280.+8.*f(n,i-is*1,j,k)
     1   -14.*f(n,i-is*2,j,k)+56.*f(n,i-is*3,j,k)/3.
     2   -35.*f(n,i-is*4,j,k)/2.+56.*f(n,i-is*5,j,k)/5.
     3   -14.*f(n,i-is*6,j,k)/3.+8.*f(n,i-is*7,j,k)/7.
     4   -f(n,i-is*8,j,k)/8.)/dx
        d2fdx2(n)=(29531.*f(n,i,j,k)/5040.-962.*f(n,i-is*1,j,k)/35.
     1   +621.*f(n,i-is*2,j,k)/10.-4006.*f(n,i-is*3,j,k)/45.
     2   +691.*f(n,i-is*4,j,k)/8.-282.*f(n,i-is*5,j,k)/5.
     3   +2143.*f(n,i-is*6,j,k)/90.-206.*f(n,i-is*7,j,k)/35.
     4   +363.*f(n,i-is*8,j,k)/560.)/dx**2
        endif
c
c        if(abs(i).eq.latx) then
cc 4th order: (since one-sided higher order may be unstable -- see README)
c        dfdx(n)=-float(is)*(-25.*f(n,i,j,k)/12.+4.*f(n,i-is*1,j,k)
c     1   -3.*f(n,i-is*2,j,k)+4.*f(n,i-is*3,j,k)/3.
c     2   -f(n,i-is*4,j,k)/4.)/dx
c        dfddx(n)=-float(is)*(-25.*fd(n,i,j,k)/12.+4.*fd(n,i-is*1,j,k)
c     1   -3.*fd(n,i-is*2,j,k)+4.*fd(n,i-is*3,j,k)/3.
c     2   -fd(n,i-is*4,j,k)/4.)/dx
c        d2fdx2(n)=(35.*f(n,i,j,k)/12.-26.*f(n,i-is*1,j,k)/3.
c     1   +19.*f(n,i-is*2,j,k)/2.-14.*f(n,i-is*3,j,k)/3.
c     2   +11.*f(n,i-is*4,j,k)/12.)/dx**2
c        endif
c
c y-derivatives:
c
c 6th order 1st- and 2nd- derivatives:
       if(abs(j).le.laty-3) then
c central derivatives in lattice interior:
        dfdy(n)=(f(n,i,j+3,k)-9.*f(n,i,j+2,k)+45.*f(n,i,j+1,k)
     1   -45.*f(n,i,j-1,k)+9.*f(n,i,j-2,k)-f(n,i,j-3,k))/(60.*dx)
        d2fdy2(n)=(f(n,i,j+3,k)/90.-3.*f(n,i,j+2,k)/20.
     1   +3.*f(n,i,j+1,k)/2.-49.*f(n,i,j,k)/18.+3.*f(n,i,j-1,k)/2.
     2   -3.*f(n,i,j-2,k)/20.+f(n,i,j-3,k)/90.)/(dx**2)
       endif
c
c one-sided derivatives near boundaries:
        is=+1
        if(j.lt.0) is=-1
c
        if(abs(j).eq.laty-2) then
        dfdy(n)=-float(is)*(f(n,i,j+is*2,k)/30.-2.*f(n,i,j+is*1,k)/5.
     1   -7.*f(n,i,j,k)/12.+4.*f(n,i,j-is*1,k)/3.-f(n,i,j-is*2,k)/2.
     2   +2.*f(n,i,j-is*3,k)/15.-f(n,i,j-is*4,k)/60.)/dx 
        d2fdy2(n)=(-13.*f(n,i,j+is*2,k)/180.+19.*f(n,i,j+is*1,k)/15.
     1   -7.*f(n,i,j,k)/3.+10.*f(n,i,j-is*1,k)/9.
     2   +f(n,i,j-is*2,k)/12.-f(n,i,j-is*3,k)/15.
     3   +f(n,i,j-is*4,k)/90.)/dx**2
        endif
c
        if(abs(j).eq.laty-1) then
        dfdy(n)=-float(is)*(-f(n,i,j+is*1,k)/6.-77.*f(n,i,j,k)/60.
     1   +5.*f(n,i,j-is*1,k)/2.-5.*f(n,i,j-is*2,k)/3.
     2   +5.*f(n,i,j-is*3,k)/6.-f(n,i,j-is*4,k)/4.
     3   +f(n,i,j-is*5,k)/30.)/dx 
        d2fdy2(n)=(137.*f(n,i,j+is*1,k)/180.-49.*f(n,i,j,k)/60.
     1   -17.*f(n,i,j-is*1,k)/12.+47.*f(n,i,j-is*2,k)/18.
     2   -19.*f(n,i,j-is*3,k)/12.+31.*f(n,i,j-is*4,k)/60.
     3   -13.*f(n,i,j-is*5,k)/180.)/dx**2
        endif
c
        if(abs(j).eq.laty) then
c 8th order:
        dfdy(n)=-float(is)*(-761.*f(n,i,j,k)/280.+8.*f(n,i,j-is*1,k)
     1   -14.*f(n,i,j-is*2,k)+56.*f(n,i,j-is*3,k)/3.
     2   -35.*f(n,i,j-is*4,k)/2.+56.*f(n,i,j-is*5,k)/5.
     3   -14.*f(n,i,j-is*6,k)/3.+8.*f(n,i,j-is*7,k)/7.
     4   -f(n,i,j-is*8,k)/8.)/dx
        d2fdy2(n)=(29531.*f(n,i,j,k)/5040.-962.*f(n,i,j-is*1,k)/35.
     1   +621.*f(n,i,j-is*2,k)/10.-4006.*f(n,i,j-is*3,k)/45.
     2   +691.*f(n,i,j-is*4,k)/8.-282.*f(n,i,j-is*5,k)/5.
     3   +2143.*f(n,i,j-is*6,k)/90.-206.*f(n,i,j-is*7,k)/35.
     4   +363.*f(n,i,j-is*8,k)/560.)/dx**2
        endif
c
c        if(abs(j).eq.laty) then
cc 4th order: (since one-sided higher order may be unstable -- see README)
c        dfdy(n)=-float(is)*(-25.*f(n,i,j,k)/12.+4.*f(n,i,j-is*1,k)
c     1   -3.*f(n,i,j-is*2,k)+4.*f(n,i,j-is*3,k)/3.
c     2   -f(n,i,j-is*4,k)/4.)/dx
c        dfddy(n)=-float(is)*(-25.*fd(n,i,j,k)/12.+4.*fd(n,i,j-is*1,k)
c     1   -3.*fd(n,i,j-is*2,k)+4.*fd(n,i,j-is*3,k)/3.
c     2   -fd(n,i,j-is*4,k)/4.)/dx
c        d2fdy2(n)=(35.*f(n,i,j,k)/12.-26.*f(n,i,j-is*1,k)/3.
c     1   +19.*f(n,i,j-is*2,k)/2.-14.*f(n,i,j-is*3,k)/3.
c     2   +11.*f(n,i,j-is*4,k)/12.)/dx**2
c        endif
c
c
c z-derivatives:
c
c 6th order 1st- and 2nd- derivatives:
       if(abs(k).le.latz-3) then
c central derivatives in lattice interior:
        dfdz(n)=(f(n,i,j,k+3)-9.*f(n,i,j,k+2)+45.*f(n,i,j,k+1)
     1   -45.*f(n,i,j,k-1)+9.*f(n,i,j,k-2)-f(n,i,j,k-3))/(60.*dx)
        d2fdz2(n)=(f(n,i,j,k+3)/90.-3.*f(n,i,j,k+2)/20.
     1   +3.*f(n,i,j,k+1)/2.-49.*f(n,i,j,k)/18.+3.*f(n,i,j,k-1)/2.
     2   -3.*f(n,i,j,k-2)/20.+f(n,i,j,k-3)/90.)/(dx**2)
       endif
c
c one-sided derivatives near boundaries:
        is=+1
        if(k.lt.0) is=-1
c
        if(abs(k).eq.latz-2) then
        dfdz(n)=-float(is)*(f(n,i,j,k+is*2)/30.-2.*f(n,i,j,k+is*1)/5.
     1   -7.*f(n,i,j,k)/12.+4.*f(n,i,j,k-is*1)/3.-f(n,i,j,k-is*2)/2.
     2   +2.*f(n,i,j,k-is*3)/15.-f(n,i,j,k-is*4)/60.)/dx 
        d2fdz2(n)=(-13.*f(n,i,j,k+is*2)/180.+19.*f(n,i,j,k+is*1)/15.
     1   -7.*f(n,i,j,k)/3.+10.*f(n,i,j,k-is*1)/9.
     2   +f(n,i,j,k-is*2)/12.-f(n,i,j,k-is*3)/15.
     3   +f(n,i,j,k-is*4)/90.)/dx**2
        endif
c
        if(abs(k).eq.latz-1) then
        dfdz(n)=-float(is)*(-f(n,i,j,k+is*1)/6.-77.*f(n,i,j,k)/60.
     1   +5.*f(n,i,j,k-is*1)/2.-5.*f(n,i,j,k-is*2)/3.
     2   +5.*f(n,i,j,k-is*3)/6.-f(n,i,j,k-is*4)/4.
     3   +f(n,i,j,k-is*5)/30.)/dx 
        d2fdz2(n)=(137.*f(n,i,j,k+is*1)/180.-49.*f(n,i,j,k)/60.
     1   -17.*f(n,i,j,k-is*1)/12.+47.*f(n,i,j,k-is*2)/18.
     2   -19.*f(n,i,j,k-is*3)/12.+31.*f(n,i,j,k-is*4)/60.
     3   -13.*f(n,i,j,k-is*5)/180.)/dx**2
        endif
c
        if(abs(k).eq.latz) then
c 8th order:
        dfdz(n)=-float(is)*(-761.*f(n,i,j,k)/280.+8.*f(n,i,j,k-is*1)
     1   -14.*f(n,i,j,k-is*2)+56.*f(n,i,j,k-is*3)/3.
     2   -35.*f(n,i,j,k-is*4)/2.+56.*f(n,i,j,k-is*5)/5.
     3   -14.*f(n,i,j,k-is*6)/3.+8.*f(n,i,j,k-is*7)/7.
     4   -f(n,i,j,k-is*8)/8.)/dx
        d2fdz2(n)=(29531.*f(n,i,j,k)/5040.-962.*f(n,i,j,k-is*1)/35.
     1   +621.*f(n,i,j,k-is*2)/10.-4006.*f(n,i,j,k-is*3)/45.
     2   +691.*f(n,i,j,k-is*4)/8.-282.*f(n,i,j,k-is*5)/5.
     3   +2143.*f(n,i,j,k-is*6)/90.-206.*f(n,i,j,k-is*7)/35.
     4   +363.*f(n,i,j,k-is*8)/560.)/dx**2
        endif
c
c        if(abs(k).eq.latz) then
cc 4th order: (since one-sided higher order may be unstable -- see README)
c        dfdz(n)=-float(is)*(-25.*f(n,i,j,k)/12.+4.*f(n,i,j,k-is*1)
c     1   -3.*f(n,i,j,k-is*2)+4.*f(n,i,j,k-is*3)/3.
c     2   -f(n,i,j,k-is*4)/4.)/dx
c        dfddz(n)=-float(is)*(-25.*fd(n,i,j,k)/12.+4.*fd(n,i,j,k-is*1)
c     1   -3.*fd(n,i,j,k-is*2)+4.*fd(n,i,j,k-is*3)/3.
c     2   -fd(n,i,j,k-is*4)/4.)/dx
c        d2fdz2(n)=(35.*f(n,i,j,k)/12.-26.*f(n,i,j,k-is*1)/3.
c     1   +19.*f(n,i,j,k-is*2)/2.-14.*f(n,i,j,k-is*3)/3.
c     2   +11.*f(n,i,j,k-is*4)/12.)/dx**2
c        endif
c
35    continue
      return
      end
