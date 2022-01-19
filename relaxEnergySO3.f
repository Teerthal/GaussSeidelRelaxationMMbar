c SO(3) energy:
c
c 15 December 2015:
c Gives same values as energyElectroweak.f for electricEnergyW
c and magneticEnergyW. So staying with the Mathematica generated
c energyElectroweak.f. However, advantage here is that we can 
c also print out the field strengths.
c
      subroutine energy(f,totalenergy,itime)
      implicit none
      include 'parameters.inc'
      include 'initialParameters.inc'
      integer n,i,j,k,itime
      integer lx, ly, lz
      integer energyfile,potentialfile,rmagneticfile,fieldfile
      integer interval
      real*8 x, y, z, dV
      real*8 kineticEnergyPhi, gradientEnergyPhi, electricEnergyW
      real*8 magneticEnergyW, electricEnergyY, magneticEnergyY
      real*8 potentialEnergy, energyDensity
      real*8 totalKEPhi, totalGEPhi, totalEEW, totalMEW, totalEEY
      real*8 totalMEY, totalPE, totalEnergy
      real*8 f
      real*8 dfdx,dfdy,dfdz,fs,d2fdx2,d2fdy2,d2fdz2
      real*8 dfddx,dfddy,dfddz
      real*8 cd
      real*8 energyInSphere,r
c
      dimension f(nf,-latx:latx,-laty:laty,-latz:latz)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension dfddx(nf),dfddy(nf),dfddz(nf)
      dimension fs(3,0:3,0:3),d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)
c cd(0:3,number of scalar fields):
      dimension cd(0:3,3)
c
c total scalar KE, gradient energy; total gauge electric, magnetic; 
c total scalar potential; total energy.
        totalKEPhi=0.
        totalGEPhi=0.
        totalEEW=0.
        totalMEW=0.
        totalEEY=0.
        totalMEY=0.
        totalPE=0.
        totalEnergy=0.
        energyInSphere=0.
c
        lx=latx-1
        ly=laty-1
        lz=latz-1
c
c     volume element:
               dV=dx**3
c
c        print*, ' energy computation sub-box semi-size = ', lx 
c     
      do 30 k=-lz,lz
         do 20 j=-ly,ly
            do 10 i=-lx,lx
c     
c     if one wants to only count the energy within a small region:
               x=float(i)*dx
               y=float(j)*dx
               z=float(k)*dx
c     
      call derivatives(f,i,j,k,dfdx,dfdy,dfdz,d2fdx2,d2fdy2,d2fdz2)
c
c covariant derivative: cd(mu,a)=(D_\mu\Phi)^a (a=1 is the real 
c part of the upper component of the doublet; a=2 is the imag
c part; a=3 is real part of the lower component; a=4 is the
c imag part of the lower component.
      call covderiv(f,i,j,k,dfdx,dfdy,dfdz,cd)
c find gauge field strengths (uses derivatives):
      call fieldstrength(f,i,j,k,dfdx,dfdy,dfdz,fs)
c
      kineticEnergyPhi= 0.5*(cd(0,1)**2+cd(0,2)**2+cd(0,3)**2)
c
      gradientEnergyPhi= 0.5*(cd(1,1)**2+cd(2,1)**2+cd(3,1)**2  
     1                       +cd(1,2)**2+cd(2,2)**2+cd(3,2)**2  
     2                       +cd(1,3)**2+cd(2,3)**2+cd(3,3)**2)
c
      electricEnergyW=0.5*(fs(1,0,1)**2+fs(1,0,2)**2+fs(1,0,3)**2
     1                    +fs(2,0,1)**2+fs(2,0,2)**2+fs(2,0,3)**2
     1                    +fs(3,0,1)**2+fs(3,0,2)**2+fs(3,0,3)**2)
c
      magneticEnergyW=0.5*(fs(1,1,2)**2+fs(1,1,3)**2+fs(1,2,3)**2
     1                    +fs(2,1,2)**2+fs(2,1,3)**2+fs(2,2,3)**2
     2                    +fs(3,1,2)**2+fs(3,1,3)**2+fs(3,2,3)**2)
c
      electricEnergyY=0.
      magneticEnergyY=0.
c
      potentialEnergy=0.25*lambda*(f(1,i,j,k)**2+f(2,i,j,k)**2
     1                       +f(3,i,j,k)**2-vev**2)**2
c
c ==============
c     
       energyDensity=kineticEnergyPhi+gradientEnergyPhi+electricEnergyW+
     1   magneticEnergyW+electricEnergyY+magneticEnergyY+potentialEnergy
c     
        totalKEPhi=totalKEPhi+kineticEnergyPhi*dV
        totalGEPhi=totalGEPhi+gradientEnergyPhi*dV
        totalEEW=totalEEW+electricEnergyW*dV
        totalMEW=totalMEW+magneticEnergyW*dV
        totalEEY=totalEEY+electricEnergyY*dV
        totalMEY=totalMEY+magneticEnergyY*dV
        totalPE=totalPE+potentialEnergy*dV
        totalEnergy=totalEnergy+energyDensity*dV
c
c==================
c energy within sphere -- useful for single monopole:
c        r=sqrt((x-xm)**2+(y-ym)**2+(z-zm)**2)
c energy within sphere -- for mmbar the center is at origin:
        r=sqrt(x**2+y**2+z**2)
        if(r.lt.float(latx-1)*dx) then
        energyInSphere=energyInSphere+energyDensity*dV
        endif
c==================
c commented out write statements for now. uncomment when needed.
c
c take nsnaps (defined in parameters.inc) snapshots of energy distribution:
c (nsnaps shouldn't be too large otherwise memory will fill up.)
c         interval=int(nt/nsnaps)
c         if(int(nt/nsnaps).eq.0) interval=1
c        if(mod(itime,interval).eq.0) then
c
c         if(i.eq.0) then
c         write(51,*) j,k,energyDensity
c         endif
c         if(j.eq.0) then
c         write(52,*) k,i,energyDensity
c         endif
c         if(k.eq.0) then
c         write(53,*) i,j,energyDensity
c         endif
c
c         if(k.eq.0) then
c         write(40,*) i,j,energyDensity
c         write(41,*) i,j,kineticEnergyPhi
c         write(42,*) i,j,gradientEnergyPhi
c         write(43,*) i,j,electricEnergyW
c         write(44,*) i,j,magneticEnergyW
c         write(45,*) i,j,electricEnergyY
c         write(46,*) i,j,magneticEnergyY
c         write(47,*) i,j,potentialEnergy
c         else
c         endif
c        endif
5             continue
10         continue
20       continue
30     continue
c     
      print*, 'time = ',itime, '     KE,GE,EEW,MEW,EEY,MEY,PE,totalE' 
      print*, totalKEPhi,totalGEPhi,totalEEW,totalMEW,
     1        totalEEY,totalMEY,totalPE,totalEnergy
      print*, ' sphere radius ', float(latx-1)*dx, 
     1                            ' energyInSphere ', energyInSphere
      print*, ' =*=*=*= '
c
c      write(91,*) itime,totalKEPhi
c      write(92,*) itime,totalGEPhi
c      write(93,*) itime,totalEEW
c      write(94,*) itime,totalMEW
c      write(95,*) itime,totalEEY
c      write(96,*) itime,totalMEY
c      write(97,*) itime,totalPE
      write(98,*) itime,totalEnergy
      write(99,*) itime,energyInSphere
c
c      write(40,*) 
c      write(51,*) 
c      write(52,*) 
c      write(53,*) 
c      write(71,*) 
c      write(72,*) 
c     
      return
      end
