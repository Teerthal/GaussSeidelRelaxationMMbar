c Relax monopole-antimonopole field configuration using Gauss-Seidel method.
c
c ====================
c begin{declarations}
      implicit none
c ===================
c begin{input}
c model-dependent parameters:
      include 'parameters.inc'
c problem-dependent parameters:
      include 'initialParameters.inc'
c end{input}
c ===================
c
c     declare all fields and time derivatives -- nf is the number of
c     fields; f denotes field. 
c
      real*8 f 
      real*8 energyinitial, totalEnergy
      real*8 totalHelicity
      real*8 energymin
      real*8 hatn,dxhatn,dyhatn,dzhatn
      dimension hatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension dxhatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension dyhatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension dzhatn(3,-latx:latx,-laty:laty,-latz:latz)
      dimension f(nf,-latx:latx,-laty:laty,-latz:latz)
c
      integer it, itime, i, j, k 
c      integer initialenergyfile
c      integer relaxedenergyfile
c
c end{declarations}
c ==================
cc ouput files:
c      relaxedenergyfile=88
c      initialenergyfile=89
c ======================
      print*, ' no. of fields = ', nf
      print*, ' lattice size = ', latx 
      print*, ' lattice spacing = ', dx
      print*, ' Gauss-Seidel relaxparam = ', relaxparam
      print*, ' number of snap shots = ', nsnaps
      print*, ' gauge function parameter = ', gp2
      print*, ' theory parameters: gw,gy,lambda,vev = ', 
     1                                       gw, gy, lambda,vev
c
c ======================
c begin{initial conditions}
      call initialconditions(f,hatn,dxhatn,dyhatn,dzhatn)
c end{initial conditions}
c ======================
c begin{write out initial energies}
c 'it' in subroutine energy is the time step (zero right now):
      it=0
      call energy(f,totalEnergy,it)
         print*, ' Total Initial Energy ', totalEnergy
c remember initial energy in box:
       energyinitial=totalEnergy
c
c initalize energymin (used to stop after the energy has become minimum)
       energymin=energyinitial
c
c       open(unit=initialenergyfile, file='initialenergy.dat', 
c     1                                             access='append')
c       write(initialenergyfile,*) zm, energyinitial
c       close(initialenergyfile)
c end{write out initial energies}
cc =======================
c      print*, ' INITIAL ENERGY ONLY -- BYPASS EVOLUTION '
c      go to 101
c begin{evolution}
c
      do 10 itime=1,nt
      call evolveeuler(f,itime,hatn,dxhatn,dyhatn,dzhatn)
      call energy(f,totalEnergy,itime)
c Stop if energy minimum has been reached. Once the energy
c starts growing, it is seen to grow rapidly, signaling an
c instability. So best to stop if energy grows:
         if(totalEnergy.le.energymin) then
              energymin=totalEnergy
         else
c              open(unit=relaxedenergyfile, file='relaxedenergy.dat', 
c     1                                             access='append')
c              write(relaxedenergyfile,*) zm, energymin
c              close(relaxedenergyfile)
c              write(77,*) zm,energyinitial,energymin,itime
              print*, ' totalEnergy.gt.energymin =  STOP '
              go to 100
         endif
c
10      continue
100     continue
101    continue
c
      stop
      end
