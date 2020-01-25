      PARAMETER (NROW=128, NGRID =256, NPAGE=NROW**2, NMAX=NGRID/2)
      PARAMETER (NRECL= NPAGE*6, NARR =MAX(2*NROW+1,NGRID+1))
      PARAMETER (NF67=MAX(NROW,NGRID/2))
      PARAMETER (Nmaxpart = 2e9)    
 
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

      COMMON / ROW /	XPAR(NPAGE),YPAR(NPAGE),ZPAR(NPAGE),
     +			VX(NPAGE),VY(NPAGE),VZ(NPAGE)
      COMMON / BWEIG/ iWeight(Nmaxpart),RWeight(Nmaxpart)
      DIMENSION         RECDAT(NRECL)  ,wspecies(10),lspecies(10)
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

      

  

 
