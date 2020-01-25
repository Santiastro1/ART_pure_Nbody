      PARAMETER (NROW=128, NGRID =256, NPAGE=NROW**2, NMAX=NGRID/2)
      PARAMETER (NRECL= NPAGE*6, NARR =MAX(2*NROW+1,NGRID+1))
      PARAMETER (NF67=MAX(NROW,NGRID/2))
      PARAMETER (Nmaxpart = 2e6)     

c    General parameters

      INTEGER Nseed
      INTEGER NWANT,Nmo,Nd ,Nbin
      character*2 npath

      PARAMETER (npath ='./')
      PARAMETER ( G           = 6.67d-8                    )!G in cgs units
      PARAMETER (cKpcm        = 3.018d21                    ) !cm/kpc           
      PARAMETER (cKmcm        = 1.d-5                        ) !km/cm
      PARAMETER (SM           =2.d33          ) !solar mass in grams
      PARAMETER ( pi          = 3.1415926535d+0    )
      PARAMETER ( twopi       = 2.d0*3.14159265d+0)
      PARAMETER (Nmo          = Nmaxpart  ) ! Max. Number of particles
      PARAMETER (Nspmax       = 10)  !Max number of species


c    Parameters to move satellite from (0,0,0) to galactic coordinates

      PARAMETER (xsat = 0.0)! Kpc
      PARAMETER (ysat = 0.0)! Kpc
      PARAMETER (zsat = 0.0)! Kpc
      PARAMETER (vxsat = 0.0)! Km/s
      PARAMETER (vysat = 0.0)! Km/s
      PARAMETER (vzsat = 0.0)! Km/s


C-----------------------------------------------------------------


          Common /connection/RMvirMd,pMd,overdens,box,
     &  RMvirMb
C-----------------------------------------------------------------

C     Halo Quantities, System of Units: Mvir/FCON = Rs = G =1
        Common /halo/aMassh,Con,Xout

            Real*8 xc,yc,zc
          Common /param/ Rs
          Common /consts/ M_0,R_0,t_0,V_0,
     &                Fcon,
     &                Nrad
           Common /coords/ xc(Nmo),yc(Nmo),zc(Nmo),
     &              vxc(Nmo),vyc(Nmo),vzc(Nmo),pMass(Nmo),
     &            peso(Nmo),Ndr,xpt(Nmo),Xb(10)


                Real*8    M_0



C-----------------------------------------------------------------


c    Parameters for the input

      COMMON / CONTROL/ AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW, 
     +                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     +                  NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     +                  ,Ocurv,extras(100)
      COMMON / HEADDR/  HEADER
      CHARACTER*45      HEADER
      COMMON /FOURAR/Zf(NARR),Yf(NARR)
      COMMON /F67COM/   
     +                 IBC,      IP,       ISL,     L1,     N2,
     +                 N3,       N4,        N7,
     +                 SI(NF67),    INDEX(NF67)

c        Real*8 XPAR,YPAR,ZPAR,VX,VY,VZ,RECDAT    
      COMMON / ROW /	XPAR(NPAGE),YPAR(NPAGE),ZPAR(NPAGE),
     +			VX(NPAGE),VY(NPAGE),VZ(NPAGE)
      COMMON / BWEIG/ iWeight(Nmaxpart),RWeight(Nmaxpart)
      DIMENSION         RECDAT(NRECL)  ,wspecies(10),lspecies(10)
      DIMENSION      lspecies0(10)
      EQUIVALENCE    (RECDAT(1),XPAR(1))
     +                               ,(wspecies(1),extras(1)),
     +                               (lspecies(1),extras(11)),
     +                               (Id,extras(90)),
     +                               (Norb,extras(91)),  
     +                               (Rsnfw,extras(92)),
     +                               (diskmass,extras(93)),
     +                               (halomass,extras(94)),
     +                               (Rdisk,extras(95)),  
     +                               (Cnfw,extras(96)),
     +                               (Nbulbo,extras(97)), 
     +                               (Qt,extras(98)), 
     +                               (Rtrunc,extras(99)),   
     +                               (Caja,extras(100))               

      

  

 
