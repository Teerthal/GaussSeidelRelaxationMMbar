c Gauged SO(3) model. Initial conditions.
c
c 15 April 2015:
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
c this declaration is now in initialParameters.inc:      
c     real*8 zm
c
      real*8 xa,ya,za
      real*8 x,y,z,xgm,ygm,zgm,xga,yga,zga,rm,ra,rhom,rhoa
      real*8 sm,sa,ms
      real*8 pprofilem,pprofilea,wprofilem,wprofilea
      real*8 f
      real*8 hatn
      real*8 phisq
      real*8 dfdx,dfdy,dfdz,d2fdx2,d2fdy2,d2fdz2
      real*8 twist,ctw,stw,ctw2,stw2
      real*8 correctionm,correctiona
      real*8 dxhatn,dyhatn,dzhatn
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
      character (len=200) itwist
      character (len=200) izmono
c
      call getarg(1,winfile)
      read(winfile,'(a500)')filename2
c
      call getarg(2,itwist)
      read(itwist,*) itw
c increase itw by 1 to increase twist by 10 degrees.
      twist=float(itw)*3.1416/18.
c
      call getarg(3,izmono)
      read(izmono,*) izm
      zm=(float(izm)+0.5)*dx
c      zm=float(izm)*dx
c     
      print*,' =========== '
      print*,' lambda, twist, zm ', lambda,twist,zm
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
c f(19)=|phi|
c end{dictionary}
c =====================
c scalar mass in model:
      ms=sqrt(2.*lambda)*vev
c
c Initial location of monopole is specified in initialParameters.inc
c (or passed on as an argument when executing) and is chosen so that 
c it is not on a lattice site.
c Antimonopole will be at (-xm,-ym,-zm):
      xa=-xm
      ya=-ym
      za=-zm
c cosine and sine of twist:
          ctw=cos(twist)
          stw=sin(twist)
          ctw2=cos(twist/2.)
          stw2=sin(twist/2.)
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
c
       xgm=(x-xm)
       ygm=(y-ym)
       zgm=(z-zm)
c 
       xga=(x-xa)
       yga=(y-ya)
       zga=(z-za)
c
c Notation is confusing -- xm,ym,zm are the location of the monopole
c (similarly antimonopole). xgm,ygm,zgm are the spatial point x,y,z
c as measured from the monopole (similarly antimonopole). But rm (ra)
c and rhom (rhoa) are the distances from the monopole and antimonopole.
c
c Radial coordinate (with boosted x) with origin on monopole:
       rm=sqrt(xgm**2+ygm**2+zgm**2)
c radial coordinate (with boosted x) with origin on antimonopole:
       ra=sqrt(xga**2+yga**2+zga**2)
c xy-radial boosted coordinate with origin on monopole:
         rhom=sqrt(rm**2-zgm**2)
c xy-radial boosted coordinate with origin on antimonopole:
         rhoa=sqrt(ra**2-zga**2)
c dimensionless distances (in units of vector mass):
         sm=gw*vev*rm
         sa=gw*vev*ra
c -- end: setup variables ------
c
c--begin scalar profiles and derivatives -----
c monopole and antimonopole scalar profile functions:
         if(sm.ne.0.) then
         pprofilem=vev*(1./tanh(sm)-(1.+ms*sm)*exp(-ms*sm)/sm)
         else
         pprofilem=0.
         endif
c
         if(sa.ne.0.) then
         pprofilea=vev*(1./tanh(sa)-(1.+ms*sa)*exp(-ms*sa)/sa)
         else
         pprofilea=0.
         endif
c
c--end scalar profiles and derivatives -----
c
c \phi^a:
c (note: monopole-antimonopole separation = 2.*z1)
c the direction of the Higgs at a point whose angular coordinates
c are (thetam,phi) and (thetaa,phia=phi) w.r.t monopole and antimonopole
c is: (sin(thetam-thetaa)*(cos(phi),sin(phi)),cos(thetam-thetaa)).
c Some algebra shows this to be:(xg*2*z1,yg*2*z1,zg*zga+rho*rhoa)/(r*ra).
c This implies that the antimonpole \phi^1 and \phi^2 are the same
c as those for the monopole in the untwisted case. However, \phi^3 
c for the antimonopole is -\phi^3 for the monopole, i.e. only the 
c third component is reversed. 
c Caution: this is not a simple "product ansatz".
c
       if(rm*ra.ne.0.) then
      hatn(1,i,j,k)=(zga*(x*ctw2+y*stw2)*ctw-zgm*(x*ctw2+y*stw2)
     1     -ra*(y*ctw2-x*stw2)*stw)/(rm*ra)
      f(1,i,j,k)=pprofilem*pprofilea*hatn(1,i,j,k)/vev
c
      hatn(2,i,j,k)=(zga*(y*ctw2-x*stw2)*ctw-zgm*(y*ctw2-x*stw2)
     1     +ra*(x*ctw2+y*stw2)*stw)/(rm*ra)
      f(2,i,j,k)=pprofilem*pprofilea*hatn(2,i,j,k)/vev
c
      hatn(3,i,j,k)=(zgm*zga+ctw*(x**2+y**2))/(rm*ra)
      f(3,i,j,k)=pprofilem*pprofilea*hatn(3,i,j,k)/vev
c
      f(19,i,j,k)=pprofilem*pprofilea/vev
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
c Following formulas assume there monopoles are not boosted (initial velocity=0).
c Also assumes monopole and antimonopole are initially on z-axis (xm=0=xa=ym=ya).
c
      if(rm*ra.ne.0.) then
         dxhatn(1,i,j,k)=(ctw2*(zga*ctw-zgm)+stw2*ra*stw
     1                 -(y*ctw2-x*stw2)*(x/ra)*stw)/(rm*ra)
     2                              -(x/rm**2+x/ra**2)*hatn(1,i,j,k)
         dyhatn(1,i,j,k)=(stw2*(zga*ctw-zgm)-ctw2*ra*stw
     1                 -(y*ctw2-x*stw2)*(y/ra)*stw)/(rm*ra)
     2                              -(y/rm**2+y/ra**2)*hatn(1,i,j,k)
         dzhatn(1,i,j,k)=((x*ctw2+y*stw2)*(ctw-1.)
     1                 -(y*ctw2-x*stw2)*(zga/ra)*stw)/(rm*ra)
     2                          -(zgm/rm**2+zga/ra**2)*hatn(1,i,j,k)
         dxhatn(2,i,j,k)=(-stw2*(zga*ctw-zgm)+ctw2*ra*stw
     1                 +(x*ctw2+y*stw2)*(x/ra)*stw)/(rm*ra)
     2                              -(x/rm**2+x/ra**2)*hatn(2,i,j,k)
         dyhatn(2,i,j,k)=(ctw2*(zga*ctw-zgm)+stw2*ra*stw
     1                 +(x*ctw2+y*stw2)*(y/ra)*stw)/(rm*ra)
     2                              -(y/rm**2+y/ra**2)*hatn(2,i,j,k)
         dzhatn(2,i,j,k)=((y*ctw2-x*stw2)*(ctw-1.)
     1                  +(x*ctw2+y*stw2)*(zga/ra)*stw)/(rm*ra)
     2                              -(zgm/rm**2+zga/ra**2)*hatn(2,i,j,k)
         dxhatn(3,i,j,k)=2.*x*ctw/(rm*ra)
     1                    -(x/rm**2+x/ra**2)*hatn(3,i,j,k)
         dyhatn(3,i,j,k)=2.*y*ctw/(rm*ra)
     1                    -(y/rm**2+y/ra**2)*hatn(3,i,j,k)
         dzhatn(3,i,j,k)=(zgm+zga)/(rm*ra)
     1                    -(zgm/rm**2+zga/ra**2)*hatn(3,i,j,k)
      else
         dxhatn(1,i,j,k)=0.
         dyhatn(1,i,j,k)=0.
         dzhatn(1,i,j,k)=0.
         dxhatn(2,i,j,k)=0.
         dyhatn(2,i,j,k)=0.
         dzhatn(2,i,j,k)=0.
         dxhatn(3,i,j,k)=0.
         dyhatn(3,i,j,k)=0.
         dzhatn(3,i,j,k)=0.
      endif
c
10         continue
20      continue
30    continue
c
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
       xga=(x-xa)
       yga=(y-ya)
       zga=(z-za)
c
c radial coordinate (with boosted x) with origin on monopole:
       rm=sqrt(xgm**2+ygm**2+zgm**2)
c radial coordinate (with boosted x) with origin on antimonopole:
       ra=sqrt(xga**2+yga**2+zga**2)
c xy-radial boosted coordinate with origin on monopole:
         rhom=sqrt(rm**2-zgm**2)
c xy-radial boosted coordinate with origin on antimonopole:
         rhoa=sqrt(ra**2-zga**2)
c dimensionless distances (in units of vector mass):
         sm=gw*vev*rm
         sa=gw*vev*ra
c
c monopole and antimonopole gauge profile function (1-K(r)):
c Note -- no factor of 1/r as in icMMbar.f and icMMbarTwisted.f
c Correction factors for profiles as evaluated by Ayush Saurabh:
        correctionm=(1.+(1.-sqrt(lambda))*sm**2/4.+sm**4/16.)/
     1                (1.+sm**2/4.+sm**4/16.)
        correctiona=(1.+(1.-sqrt(lambda))*sa**2/4.+sa**4/16.)/
     1                (1.+sa**2/4.+sa**4/16.)
        if(sm.ne.0.) then
        wprofilem=(1.-(sm/sinh(sm))*correctionm)/gw
        else
        wprofilem=0.
        endif
        if(sa.ne.0.) then
        wprofilea=(1.-(sa/sinh(sa))*correctiona)/gw
        else
        wprofilea=0.
        endif
c square of scalar field:
       phisq=f(1,i,j,k)**2+f(2,i,j,k)**2+f(3,i,j,k)**2
c
      call derivatives(f,i,j,k,dfdx,dfdy,dfdz,d2fdx2,d2fdy2,d2fdz2)
c
       if(phisq.ne.0.) then
c W_\mu^1:
       f(4,i,j,k)=0.
       f(5,i,j,k)=-gw*wprofilem*wprofilea
     1  *(hatn(2,i,j,k)*dxhatn(3,i,j,k)-hatn(3,i,j,k)*dxhatn(2,i,j,k))
       f(6,i,j,k)=-gw*wprofilem*wprofilea
     1  *(hatn(2,i,j,k)*dyhatn(3,i,j,k)-hatn(3,i,j,k)*dyhatn(2,i,j,k))
       f(7,i,j,k)=-gw*wprofilem*wprofilea
     1  *(hatn(2,i,j,k)*dzhatn(3,i,j,k)-hatn(3,i,j,k)*dzhatn(2,i,j,k))
c
c W_\mu^2:
       f(8,i,j,k)=0.
       f(9,i,j,k)=-gw*wprofilem*wprofilea
     1  *(hatn(3,i,j,k)*dxhatn(1,i,j,k)-hatn(1,i,j,k)*dxhatn(3,i,j,k))
       f(10,i,j,k)=-gw*wprofilem*wprofilea
     1  *(hatn(3,i,j,k)*dyhatn(1,i,j,k)-hatn(1,i,j,k)*dyhatn(3,i,j,k))
       f(11,i,j,k)=-gw*wprofilem*wprofilea
     1  *(hatn(3,i,j,k)*dzhatn(1,i,j,k)-hatn(1,i,j,k)*dzhatn(3,i,j,k))
c
c W_\mu^3:
       f(12,i,j,k)=0.
       f(13,i,j,k)=-gw*wprofilem*wprofilea
     1  *(hatn(1,i,j,k)*dxhatn(2,i,j,k)-hatn(2,i,j,k)*dxhatn(1,i,j,k))
       f(14,i,j,k)=-gw*wprofilem*wprofilea
     1  *(hatn(1,i,j,k)*dyhatn(2,i,j,k)-hatn(2,i,j,k)*dyhatn(1,i,j,k))
       f(15,i,j,k)=-gw*wprofilem*wprofilea
     1  *(hatn(1,i,j,k)*dzhatn(2,i,j,k)-hatn(2,i,j,k)*dzhatn(1,i,j,k))
c
       else
       f(4,i,j,k)=0.
       f(5,i,j,k)=0.
       f(6,i,j,k)=0.
       f(7,i,j,k)=0.
       f(8,i,j,k)=0.
       f(9,i,j,k)=0.
       f(10,i,j,k)=0.
       f(11,i,j,k)=0.
       f(12,i,j,k)=0.
       f(13,i,j,k)=0.
       f(14,i,j,k)=0.
       f(15,i,j,k)=0.
       endif
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
