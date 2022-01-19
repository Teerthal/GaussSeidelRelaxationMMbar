      subroutine evolveeuler(f,itime,hatn,dxhatn,dyhatn,dzhatn)
      implicit none
      include 'parameters.inc'
c need initialParameters to implement the *constrained* relaxation.
c we need to hold the position of the mmbar fixed and these are given
c in the initialParameters. Better would be to impose the constraints
c when calling the evolution subroutines and not in the evolution
c routine itself.
      include 'initialParameters.inc'
      integer n,i,j,k,is,js,ks,nn
      real*8 f,r,s
c      real*8 relaxparam
      real*8 error
      integer itime
      real*8 hatn,dxhatn,dyhatn,dzhatn
      real*8 one
      real*8 coeff
c
      dimension hatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension dxhatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension dyhatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension dzhatn(3,-latx:latx,-laty:laty,-latz:latz)
c
      dimension f(nf,-latx:latx,-laty:laty,-latz:latz)
      dimension r(nf),s(nf)
      dimension error(nf)
c
c e.g. eq.(107) in http://www.damtp.cam.ac.uk/lab/people/sd/lectures/nummeth98/pdes.htm
c
c
c evolve on lattice (one-sided differencing on boundaries as
c defined in subroutine derivatives):: 
c
           do 41 n=1,nf
           error(n)=0.
41         continue
c
c The Gauss-Seidel method is sensitive to how one loops over the
c lattice. I've tried looping from one end of the lattice to the
c other. But if I start from the center of the lattice, the code
c seems to do better (converge faster). Should explore this further.
c
c
      do 10 k=-latz,latz
         do 20 j=-laty,laty
            do 30 i=-latx,latx
c
      call fdflux(f,i,j,k,r,hatn,dxhatn,dyhatn,dzhatn)
c--------------------------------------------------------
c coeff is the coefficient of f(n,i,j,k) terms in the laplacian
c (see fluxesNumRelSO3.f) in the bulk of the lattice. 
       coeff=49./6.
c--------------------------------------------------------
c
c r(n) is the flux. It is also the error. Write out the
c error:
       if(j.eq.0.and.itime.eq.1) then
       write(51,*) i,k,r(1)
       write(52,*) i,k,r(2)
       write(53,*) i,k,r(3)
       write(54,*) i,k,r(5)
       write(55,*) i,k,r(6)
       write(56,*) i,k,r(7)
       endif
       if(j.eq.0.and.itime.eq.int(nt/2)) then
       write(61,*) i,k,r(1)
       write(62,*) i,k,r(2)
       write(63,*) i,k,r(3)
       write(64,*) i,k,r(5)
       write(65,*) i,k,r(6)
       write(66,*) i,k,r(7)
       endif
       if(j.eq.0.and.itime.eq.nt-1) then
       write(71,*) i,k,r(1)
       write(72,*) i,k,r(2)
       write(73,*) i,k,r(3)
       write(74,*) i,k,r(5)
       write(75,*) i,k,r(6)
       write(76,*) i,k,r(7)
       endif
c
      if(abs(i).ne.latx.and.abs(j).ne.laty.and.abs(k).ne.latz) then
c update scalar fields (note: r(1),r(2),r(3) not used):
      f(19,i,j,k)= f(19,i,j,k)+relaxparam*dx**2*r(19)/coeff
c
      f(1,i,j,k)=f(19,i,j,k)*hatn(1,i,j,k)
      f(2,i,j,k)=f(19,i,j,k)*hatn(2,i,j,k)
      f(3,i,j,k)=f(19,i,j,k)*hatn(3,i,j,k)
c
c update gauge fields:
        do 40 n=4,nf-1
c Gauss-Seidel method:
         f(n,i,j,k)= f(n,i,j,k)+relaxparam*dx**2*r(n)/coeff
c after relaxation, r(n) should be nearly zero:
c         error(n)=error(n)+r(n)**2
40      continue
c
         else
         endif
c
30         continue
7     continue
20      continue
8     continue
10    continue
9     continue
c
      return
      end
