c     --------------------------------
c     Control parameters and variables
c     --------------------------------
      

      parameter ( wsplit = 0.5                )  ! critical weight for split
      parameter ( wjoin  = 0.25                )  ! critical weight for join

c     ............numerical controls...........

      common /NUMC01/ niter,niter1,niter2     ! # of relax. iterations

c     ........interpolation coefficients.......

      common /INTER1/  wa, wbcd                  ! prolongation weights
      
c     ....run control parameters & variables...

      real AEXPN,AEXP0,AMPLT,ASTEP,PARTW,
     &               TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &               Om0,Oml0,hubble,Wp5

      common / RUN /    AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,Nspecs,Nseed,Om0,Oml0,hubble,Wp5,
     &                  Ocurv,extras(100)

      character*45      HEADER
      common / HEADDR/  HEADER

c     .................timing..................

      common / TIMECPU / CPU(10)

      parameter ( iostep = 1  )  ! output step
      parameter ( nsave  = 40 )

      dimension   asave (nsave)                  ! time moments to save
      character*6 aname (nsave)

c      parameter ( queue_time =  29.)     !  queue_time in minutes 
c     parameter ( queue_time =  1420.*6. )   !  queue_time 23.8 hrs
       parameter ( queue_time =  1.E8)    !  unlimited

      data asave /
     & 0.0250,0.0500,0.0750,0.1000,0.1250,0.1500,0.1750,
     & 0.2000,0.2250,0.2500,0.2750,0.3000,0.3250,0.3500,
     & 0.3750,0.4000,0.4250,0.4500,0.4750,0.5000,0.5250,
     & 0.5500,0.5750,0.6000,0.6250,0.6500,0.6750,0.7000,
     & 0.7250,0.7500,0.7750,0.8000,0.8250,0.8500,0.8750,
     & 0.9000,0.9250,0.9500,0.9750,1.0000/

      data aname /
     &'0.0250','0.0500','0.0750','0.1000','0.1250','0.1500','0.1750',
     &'0.2000','0.2250','0.2500','0.2750','0.3000','0.3250','0.3500',
     &'0.3750','0.4000','0.4250','0.4500','0.4750','0.5000','0.5250',
     &'0.5500','0.5750','0.6000','0.6250','0.6500','0.6750','0.7000',
     &'0.7250','0.7500','0.7750','0.8000','0.8250','0.8500','0.8750',
     &'0.9000','0.9250','0.9500','0.9750','1.0000'/
