c     ---------------------------------------------
c     code setup : control constants and array sizes
c     ---------------------------------------------
c....                                                                global setup constants
      parameter ( ndim     = 3           ) ! # of spatial dimensions
      parameter ( mcell    = 60 000 000    ) ! # of refinement cells
      parameter ( nrow     = 128    ) ! # of particles in 1d
      parameter ( ng       = 128     ) ! # of 0 lv cells in 1d
      parameter ( MinLevel = 0          ) ! minimum allowed level
      parameter ( MaxLevel = 8      ) ! maximum allowed level <= 10
      parameter ( ioff     = 0          ) ! buffer size for arrays
      parameter ( nvar     = 2          ) ! # of physical variables
      parameter ( nchem    = 0          ) ! # of chemical species
      parameter ( nspecies = 3         ) ! # of particle species
      parameter ( ncell0 = ng**3        ) ! # of zero level cells
 
      parameter ( NCPUs    = 4    ) ! anticipated number of processors
      parameter ( Nproc    = NCPUs      )                                 
      parameter ( LBuffer  =  6000 000    ) ! Length of buffer for Dens
      parameter ( MLChunk  = 15000      ) ! Max Cells per processor for Dens
      parameter ( nbyteword = 4          )  ! defines length of direct-access record

c      parameter ( NBuffer = 4*(mcell -ncell0))
      parameter ( NBuffer = MAX(4*(mcell -ncell0),3*Nproc*LBuffer ) )

c.... secondary setup constants
     
c     ..........particles...............

      parameter ( np       = 2097152 )          ! # of particles
      parameter ( npage    = nrow**2        ) ! # of particles in a page
      parameter ( nrecl    = npage * 6      ) ! length of particle row in words

c     .........zero level...............

      parameter ( ng2    = ng**2        ) ! # of cells in a grid layer
      parameter ( narr   = ng + 1       ) ! FFT parameter
      parameter ( nf67   = ng/2         ) ! FFT parameter
      parameter ( xn     = float(ng)+1.-1.E-7 ) ! boundaries
      parameter ( yn     = float(ng)    ) 
      parameter ( nctot  = mcell ) ! total number of cells

c     ............tree..................

      parameter ( nneigh = 2*ndim       ) ! # of neighbors
      parameter ( neighb = 2*ndim       ) ! # of neighbors
      parameter ( nchild = 2**ndim      ) ! # of children
      parameter ( moct   = (mcell-ncell0)/nchild ) ! # of octs (old nky)








