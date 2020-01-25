c     --------------
c     Tree structure
c     --------------

      include 'a_numbers.h'
      include 'a_setup.h' 

c     ndim   - number of dimensions
c     nchild - number of children    ( nchild = 2**ndim        )
c     nneigh - number of neighbours  ( nneigh = 2*ndim         )
c     mcell  - number of tree cells  ( mcell  = ncell0 + nchild * moct  )
c     moct   - length of oct arrays  ( moct   = (mcell-ncell0) / nchild )
c     nctot  - total number of cells

c     iOctLv :    >0   - level of an oct
c     iOctPr :         - parent of an oct
c     iOctCh :    >0   - pointer to an oct of children
c                  0   - there are no children; the cell is a leaf
c     iOctNb :    >0   - pointers to neighbouring cells 
c     iOctPs :         - coordinates of Oct centers

      parameter ( iFreeLevel = -1000 )    ! fake level for free (unused) cells

c     ............positions.................

      parameter ( isize     = 2**MaxLevel            )		  
      parameter ( d_x       = 1.0 / 2**MaxLevel      ) ! differs from Lesha!
      parameter ( nshift    = nchild - 1             )
      parameter ( nbshift   = nshift - ncell0        ) ! big shift; see Tools
      dimension idelta(8,3)
      data      idelta / -1,  1, -1,  1, -1,  1, -1,  1,
     &                   -1, -1,  1,  1, -1, -1,  1,  1, 
     &                   -1, -1, -1, -1,  1,  1,  1,  1 /

c     ..........oct information..............

      common /TREE01/ iOctPs(0:moct,ndim)       ! coordinates of Oct centers 
      common /TREE02/ iOctNb(0:moct,nneigh)     ! neighbouring cells
      common /TREE03/ iOctPr(0:moct)            ! parents/linked list index
      common /TREE04/ iOctLv(0:moct)            ! Level 
      common /TREE05/ iOctCh(0:nctot)           ! children Octs
      common /TREE06/ iOctFree                  ! first free Oct
      common /TREE07/ nOct                      ! number of Octs in use
      common /TREE08/ iOctMax                   ! oct with max. index
c     ........linked list of octs............

      common /LIST01/ iOctLL1(0:moct)              ! doubly linked list of octs
      common /LIST02/ iOctLL2(0:moct)
      common /LIST03/ iHOLL(MinLevel:MaxLevel)  ! linked list header
      common /LIST04/ iNOLL(MinLevel:MaxLevel)  ! # of ll entries at a Level

c     ..........common arrays................

      common /COMM01/ iCL(nctot)                ! headers of particle LList
      common /COMM02/ iNp(nctot)                ! number of particles in cell

c     ........physical variables.............
c
c     var(1,*)    - density 
c     var(2,*)    - potential 
c
c     during refinement modifications var is used 
c     to store refinement/derefinement indicators
c
      common /VAR01/ var(2,nctot)               ! vectors of physical variables
      common /VAR02/ pot(nctot)
      common /REF01/ iBufferDummy(NBuffer)
      common /REF02/ ref(nctot)      ! density of only small particles 
c      common /REF03/ iBuffer(LBuffer,Nproc)
c      common /REF04/ wBuffer(LBuffer,Nproc)
c      common /REF05/ rBuffer(LBuffer,Nproc)

      Integer iLevNb(3,nctot-ncell0)
      Integer iBuffer(LBuffer,Nproc)
      Real*4  wBuffer(LBuffer,Nproc) 
      Real*4  rBuffer(LBuffer,Nproc)
      Equivalence (iBufferDummy(1),iLevNb(1,1))
      Equivalence (iBufferDummy(1),iBuffer(1,1))
      Equivalence (iBufferDummy(LBuffer*Nproc+1),wBuffer(1,1))
      Equivalence (iBufferDummy(2*LBuffer*Nproc+1),rBuffer(1,1))
c     ..........scratch arrays...............

      common /TEMP01/ iSelect(moct)             

c     ..........miscellaneous................

      common /MISC01/ CellSize(MinLevel:MaxLevel)                      

c     .........refinement controls...........

      common /RFMC01/ nsplit(MinLevel:MaxLevel)  ! number  thresholds 
      common /RFMC02/ trho  (MinLevel:MaxLevel)  ! density thresholds 

c     .........restricted refinement zone......

      common /RFMZONE/ fsmall, imin , imax , jmin , jmax , kmin , kmax

c              trho = wpar * nsplit

c     ..........tree related controls..........

      common /CELLS1/  ncell,                    ! number of cells in use
     &                 MaxLevelNow               ! current maximum level
      common /CUBE01/  xmin , ymin , zmin,       ! coordinates of cube corners
     &                 xmax , ymax , zmax


c     .............particle arrays.............
c
c     particles move with different time steps the current 
c     time moment, ptime, and time step, pdt, of a particle 
c     are stored in pt in the following manner
c               
c               ptime = int(pt)*dAmin + Aexpn
c               pdt   = pt - int(pt)
c     where
c               dAmin = astep / 2.0**iTimeBin(MaxLevel)
c
      real pt(np)
	integer ipt(np)
      double precision x(np), y(np), z(np), vx(np), vy(np), vz(np)	

      common / PART01 /      x
      common / PART02 /      y
      common / PART03 /      z
      common / PART04 /     vx
      common / PART05 /     vy
      common / PART06 /     vz
      common / PART07 / pdummy(np)                    ! scratch (for epot)
      common / PART08 /     pt,ipt                     ! particle time & step
      common / PART09 /  iLL(np,2)                     ! doubly linked list
      common / PART10 /  Wpart(np)                     ! weights of particles       

c      common / PART01 /      x(np)
c      common / PART02 /      y(np)
c      common / PART03 /      z(np)
c      common / PART04 /     vx(np)
c      common / PART05 /     vy(np)
c      common / PART06 /     vz(np)
c      common / PART07 / pdummy(np)                    ! scratch (for epot)
c      common / PART08 /     pt(np)                     ! particle time & step
c      common / PART09 /  iLL(np,2)                     ! doubly linked list

c     .........information about species........
      
      common / SPEC01 /  wpar(nspecies)                   ! particle weight
c      common / SPEC02 /    sw(nspecies)                   ! relative weight
      common /  SPEC03 /  nsp(nspecies,2), lsp(nspecies)
c      data nsp / 1       , ncold1 , 
c     &           ncold , np      / 
c      data lsp / ncold , np /
c      
c     ..................time...................

      common / TIME01 / iTimeBin(MinLevel:MaxLevel),     ! timebin 
     &                    Aexp(MinLevel:MaxLevel),     ! time moment 
     &                   dAexp(MinLevel:MaxLevel),     ! = 1/2**iTimeBin
     &                   dAmin                         ! = dAexp(MaxLevel)

c     ..........I/O auxiliary arrays..........

      real xpar(npage),ypar(npage),zpar(npage),
     &     vxx(npage),vyy(npage),vzz(npage)

      common / ROW /   xpar,
     &                 ypar,
     &                 zpar,
     &		        vxx, 
     &                  vyy, 
     &                  vzz

c      common / ROW /   xpar(npage),
c     &                 ypar(npage),
c     &                 zpar(npage),
c     &		        vxx(npage), 
c     &                  vyy(npage), 
c     &                  vzz(npage)

      real             recdat(nrecl)                   ! record storage
      equivalence      (recdat(1),xpar(1))              ! position at row



