c Just a monopole.
c
c Gauged SO(3) model. Initial conditions.
c
c 15 April 2015:
c 
c The initial BPS monopole profile is taken from Vilenkin & Shellard, 
c Sec. 14.1 (1994):
c phi=vev*h(r)*{\hat x}
c A_i^a = (1-K(r)) \epsilon^{aij} x^j/(er^2) (note: overall sign in V&S is incorrect)
c h = 1./tanh(\xi)- 1/\xi
c K = \xi /sinh(\xi)
c \xi = m_v r = gw*vev*r
c For the non-BPS case, this is modified to:
c h = 1./tanh(\xi)- (1+m \xi)e^(-m \xi)/\xi
c K = \xi /sinh(\xi)
c where $m^2 = 2\lambda \eta^2$ ($\lambda$ occurs in the Higgs potential).
c Note: $m$ is really the scalar mass $m_s$, here written as "ms".
c
      subroutine initialconditions(f,hatn,dxhatn,dyhatn,dzhatn)
      implicit none
      include 'parameters.inc'
      include 'initialParameters.inc'
      integer i,j,k,n
      integer itw,izm
c
c      real*8 zm (now declared in initialParameters.inc)
      real*8 x,y,z,xgm,ygm,zgm,rm,rhom
      real*8 sm,ms,mv
      real*8 pprofilem,wprofilem
c      real*8 pdashm,wdashm
      real*8 f
      real*8 hatn
      real*8 dxhatn,dyhatn,dzhatn
      real*8 phisq
      real*8 dfdx,dfdy,dfdz,d2fdx2,d2fdy2,d2fdz2
      real*8 twist
      real*8 correctionm
c 
      dimension f(nf,-latx:latx,-laty:laty,-latz:latz)
      dimension hatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension dxhatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension dyhatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension dzhatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)
c
c===================
c This section is for running the code on several nodes at the
c same time, each code with a different value of the parameters.
c read scan parameters and create files to write out successful runs:
      character*500 filename2
      character (len=200) winfile
      character (len=200) izmono
c
      call getarg(1,winfile)
      read(winfile,'(a500)')filename2
c
      call getarg(2,izmono)
      read(izmono,*) izm
      zm=(float(izm)+0.5)*dx
c     
      print*,' =========== '
      print*,' lambda, zm, dx, latx*dx ', lambda,zm,dx,float(latx)*dx
      print*,' =========== '
c
c initialize:
      do 35 i=-latx,latx
         do 25 j=-laty,laty
            do 15 k=-latz,latz
               do 50 n=1,nf
                  f(n,i,j,k)=0.
50            continue
            do 51 n=1,3
               hatn(n,i,j,k)=0.
               dxhatn(n,i,j,k)=0.
               dyhatn(n,i,j,k)=0.
               dzhatn(n,i,j,k)=0.
51          continue
15         continue
25      continue
35    continue
c
c =====================
c begin{dictionary}
c Gauged SO(3) model contains 3 real scalar fields:
c f(1)= phi(1), f(2)=phi(2), f(3)=phi(3), 
c and gauge fields: 
c f(4)=W^1_0, f(5)=W^1_1, f(6)=W^1_2, f(7)=W^1_3,
c f(8)=W^2_0, f(9)=W^2_1, f(10)=W^2_2, f(11)=W^2_3,
c f(12)=W^3_0, f(13)=W^3_1, f(14)=W^3_2, f(15)=W^3_3.
c Also: fd(n)=dot(f(n)). 
c end{dictionary}
c =====================
c scalar mass in model:
      ms=sqrt(2.*lambda)*vev
      mv=gw*vev
      print*, ' ms, mv, & ms/mv: ', ms,mv,ms/mv
c
c start: loop over larger lattice (+1) for scalar fields:
      do 30 i=-latx,latx
         do 20 j=-laty,laty
            do 10 k=-latz,latz
c
c -- begin: setup variables ------
c (x,y,z) for lattice site:
               x=float(i)*dx
               y=float(j)*dx
               z=float(k)*dx
       xgm=x-xm
       ygm=y-ym
       zgm=z-zm
c
c radial coordinate (with unboosted x) with origin on monopole:
       rm=sqrt(xgm**2+ygm**2+zgm**2)
c xy-radial unboosted coordinate with origin on monopole:
         rhom=sqrt(rm**2-zgm**2)
c dimensionless distances (in units of vector mass):
         sm=gw*vev*rm
c -- end: setup variables ------
c
c--begin scalar profile -----
c monopole and antimonopole scalar profile functions:
        pprofilem=vev*(1./tanh(sm)-(1.+ms*sm)*exp(-ms*sm)/sm)
c--end scalar profile -----
c
c \phi^a:
       if(rm.gt.0.001) then
      hatn(1,i,j,k)=xgm/rm
      hatn(2,i,j,k)=ygm/rm
      hatn(3,i,j,k)=zgm/rm
c
      dxhatn(1,i,j,k)=1./rm-xgm**2/rm**3
      dyhatn(1,i,j,k)=-xgm*ygm/rm**3
      dzhatn(1,i,j,k)=-xgm*zgm/rm**3
      dxhatn(2,i,j,k)=-ygm*xgm/rm**3
      dyhatn(2,i,j,k)=1./rm-ygm**2/rm**3
      dzhatn(2,i,j,k)=-ygm*zgm/rm**3
      dxhatn(3,i,j,k)=-zgm*xgm/rm**3
      dyhatn(3,i,j,k)=-zgm*ygm/rm**3
      dzhatn(3,i,j,k)=1./rm-zgm**2/rm**3
c
      f(1,i,j,k)=pprofilem*hatn(1,i,j,k)
      f(2,i,j,k)=pprofilem*hatn(2,i,j,k)
      f(3,i,j,k)=pprofilem*hatn(3,i,j,k)
c
      f(19,i,j,k)=pprofilem/vev
c
       else
       hatn(1,i,j,k)=0.
       hatn(2,i,j,k)=0.
       hatn(3,i,j,k)=1.
c
       f(1,i,j,k)=0.
       f(2,i,j,k)=0.
       f(3,i,j,k)=0.
       f(19,i,j,k)=0.
       endif
c
10         continue
20      continue
30    continue
c
c now for the initial gauge fields and the electric fields:
c
      do 32 i=-latx,latx
         do 22 j=-laty,laty
            do 12 k=-latz,latz
c 
c -- begin: setup variables ------
c (x,y,z) for lattice site:
               x=float(i)*dx
               y=float(j)*dx
               z=float(k)*dx
c
       xgm=(x-xm)
       ygm=(y-ym)
       zgm=(z-zm)
c 
c radial coordinate (with unboosted x) with origin on monopole:
       rm=sqrt(xgm**2+ygm**2+zgm**2)
c xy-radial unboosted coordinate with origin on monopole:
         rhom=sqrt(rm**2-zgm**2)
c dimensionless distances (in units of vector mass):
         sm=gw*vev*rm
c monopole gauge profile function (1-K(r)):
        correctionm=(1.+(1.-sqrt(lambda))*sm**2/4.+sm**4/16.)/
     1                (1.+sm**2/4.+sm**4/16.)
        wprofilem=(1.-(sm/sinh(sm))*correctionm)/gw
c
c W_\mu^1:
       f(4,i,j,k)=0.
       f(5,i,j,k)=-wprofilem
     1  *(hatn(2,i,j,k)*dxhatn(3,i,j,k)-hatn(3,i,j,k)*dxhatn(2,i,j,k))
       f(6,i,j,k)=-wprofilem
     1  *(hatn(2,i,j,k)*dyhatn(3,i,j,k)-hatn(3,i,j,k)*dyhatn(2,i,j,k))
       f(7,i,j,k)=-wprofilem
     1  *(hatn(2,i,j,k)*dzhatn(3,i,j,k)-hatn(3,i,j,k)*dzhatn(2,i,j,k))
c       f(5,i,j,k)=0.
c       f(6,i,j,k)=+gw*wprofilem*zgm/rm**2
c       f(7,i,j,k)=-gw*wprofilem*ygm/rm**2
c W_\mu^2:
       f(8,i,j,k)=0.
       f(9,i,j,k)=-wprofilem
     1  *(hatn(3,i,j,k)*dxhatn(1,i,j,k)-hatn(1,i,j,k)*dxhatn(3,i,j,k))
       f(10,i,j,k)=-wprofilem
     1  *(hatn(3,i,j,k)*dyhatn(1,i,j,k)-hatn(1,i,j,k)*dyhatn(3,i,j,k))
       f(11,i,j,k)=-wprofilem
     1  *(hatn(3,i,j,k)*dzhatn(1,i,j,k)-hatn(1,i,j,k)*dzhatn(3,i,j,k))
c       f(9,i,j,k)=-gw*wprofilem*zgm/rm**2
c       f(10,i,j,k)=0.
c       f(11,i,j,k)=+gw*wprofilem*xgm/rm**2
c W_\mu^3:
       f(12,i,j,k)=0.
       f(13,i,j,k)=-wprofilem
     1  *(hatn(1,i,j,k)*dxhatn(2,i,j,k)-hatn(2,i,j,k)*dxhatn(1,i,j,k))
       f(14,i,j,k)=-wprofilem
     1  *(hatn(1,i,j,k)*dyhatn(2,i,j,k)-hatn(2,i,j,k)*dyhatn(1,i,j,k))
       f(15,i,j,k)=-wprofilem
     1  *(hatn(1,i,j,k)*dzhatn(2,i,j,k)-hatn(2,i,j,k)*dzhatn(1,i,j,k))
c       f(13,i,j,k)=+gw*wprofilem*ygm/rm**2
c       f(14,i,j,k)=-gw*wprofilem*xgm/rm**2
c       f(15,i,j,k)=0.
c
12         continue
22      continue
32    continue
c
c gauge functions $\Gamma^a = \partial_i W^a_i$:
c This uses a call to derivatives. That is why it is necessary to
c evaluate the gauge functions in a separate do loop (after all the
c gauge fields have been specified).
c
      do 72 i=-latx,latx
         do 73 j=-laty,laty
            do 74 k=-latz,latz
c
      call derivatives(f,i,j,k,dfdx,dfdy,dfdz,d2fdx2,d2fdy2,d2fdz2)
c
       f(16,i,j,k)=dfdx(5)+dfdy(6)+dfdz(7)
       f(17,i,j,k)=dfdx(9)+dfdy(10)+dfdz(11)
       f(18,i,j,k)=dfdx(13)+dfdy(14)+dfdz(15)
c
74         continue
73      continue
72    continue
c
      return
      end
