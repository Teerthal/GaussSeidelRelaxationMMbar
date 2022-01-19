c For electroweak model:
c (following conventions in Achucarro and Vachaspati,
c  Phys. Rept. 327 (2000) 347-426.)
c
c nf=number of fields; nt=number of time steps;
c lat=semi-lattice size; dx=spatial spacing;
c alpha, beta etc. are parameters of the model
c
c Allows for different lattice sizes along different axes. 
c
c Fiducial lattice only goes up to lat-1. The last layer on the 
c lattice is a boundary layer where only some variables are calculated.
c
      integer nf,latx,laty,latz,nt,nsnaps,outerlayer
      real*8 dx,gw,gy,gp2,lambda,vev,relaxparam
c For abelian-Higgs model: gy is the gauge coupling.
c Note: conventions like those in electroweak but truncated to hypercharge.
c model parameters:
c the number of fields (nf) includes real scalar fields, gauge fields,
c and the gauge constraints. E.g. for electroweak there are 4 real
c scalars, 4x4=16 gauge fields, and 4 gauge constraints (Lorenz gauge
c conditions. (gauge constraint equations are first order.)
c For SO(3), there are 15+3 fields. 
c The 19th field for SO(3) is the absolute value of phi since now 
c working in the gauge where the direction of phi is fixed, only
c the magnitude varies.
c
      parameter(nf=19,gw=1.0,gy=0.0,lambda=1.0,vev=1.0)
ccc      parameter(nf=18,gw=1.0,gy=0.0,lambda=1.0,vev=1.0)
c
c simulation parameters next:
c gp2=gauge function parameter ($gp2=a^2$ of Num. Rel., Baumgarte & Shapiro,
c see Chapter 11). $0 < gp2 < 1$. 
c relaxparam = relaxation parameter in "successive over relaxation". (Similar
c to dt in the evolution code.)
c e.g. eq.(107) in http://www.damtp.cam.ac.uk/lab/people/sd/lectures/nummeth98/pdes.htm
c
      parameter(latx=32,laty=latx,latz=latx,nt=1000,dx=0.25,
     1                 gp2=0.75,relaxparam=1.0)
c
c Output parameter -- how many snapshots of the energy distribution?
      parameter(nsnaps=20)
c
c Calculate energy in a sub-box with dimension -lat+outerlayer,+lat-outerlayer.
c outerlayer must be greater than equal to 1 (see energy subroutine).
      parameter(outerlayer=2)
c
