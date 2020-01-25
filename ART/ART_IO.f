C     =============================================================
c
c      ART Version 3: ART_IO.f - routines dealing with Input/Output
c
c     =============================================================

c     ------------------------
      subroutine Run_Output_o ()
c     ------------------------
c
c     purpose: creates an output of current run information 
c 
      include 'a_control.h'
      
c      write (*,*) ' Close Log'
      call Close_Log ()
c      write (*,*) ' Open Log'
      call Open_Log  ()
      if ( mod(istep,iostep) .eq. 0 ) then 
         write (*,*) ' istep=',istep,' iostep=',iostep
         write (*,*) ' Write Particles'
            call Write_Particles_Binary1_o ( 'PMcrd.DAT' ,
     &                                 'PMcrs0.DAT')
         write (*,*) ' Write pt'
        call Write_Particle_Time     ( 'pt.dat' )
         write (*,*) ' Done writting'
      endif

      do imoment = 1 , nsave
         call Save_Moment_o ( imoment )
      enddo
c     call Save_Moment_o (ISTEP )  

      return
      end

c     ----------------------------------
      subroutine Save_Moment_o ( imoment )
c     ----------------------------------
c
      include 'a_control.h'
      integer imoment 
      character*16 SaveFile1
      character*17 SaveFile2
      character*13 SaveFile3

      SaveFile1 = 'PMcrda' // aname(imoment) // '.DAT'
      SaveFile2 = 'PMcrs0a'// aname(imoment) // '.DAT'
      SaveFile3 = 'pta'    // aname(imoment) // '.dat'

c     write(SaveFile1,'(a,i4.4,a)') 'PMcrda0.',ISTEP,'.DAT'
c     write(SaveFile2,'(a,i4.4,a)') 'PMcrs0a0.',ISTEP,'.DAT'
c     write(SaveFile3,'(a,i4.4,a)') 'pta0.',ISTEP,'.dat'

      ad = asave(imoment)
      if ((aexpn.ge.(ad-0.5*astep)).and.(aexpn.lt.(ad+0.5*astep))) then 
         write (*,*) ' Save Moment:'
         call Write_Particles_Binary1_o ( SaveFile1 , SaveFile2 )
         call Write_Particle_Time     ( SaveFile3 )
         write (*,*) ' Done Saving'
      endif
      
      return
      end
      
c     ----------------------
      subroutine Open_Log ()
c     ----------------------
c
c     purpose: open files used to log the current run
c    
c     the files are opened to be appended, not rewritten
c

c.... position='append' - for AIX; access = 'append' for others

      open (40 , file = 'timing.log', access = 'append')
      open (50 , file = 'run.log', access = 'append' )
c      open (40 , file = 'timing.log', position = 'append')
c      open (50 , file = 'run.log', position = 'append' )

      return
      end

c     -----------------------
      subroutine Close_Log ()
c     -----------------------
c
      close (40)
      close (50)
      return
      end


c     --------------------------------
      subroutine Read_Particle_Time ()
c     --------------------------------
c
c     purpose: when the code starts or restarts
c              this routine takes care of array pt
c 
      include 'a_tree.h'
      include 'a_control.h'
      DIMENSION         wspecies(nspecies),lspecies(nspecies)
      EQUIVALENCE    (wspecies(1),extras(1)),
     +                               (lspecies(1),extras(11))
      open ( 60 ,file = 'pt.dat',
     &       form = 'unformatted',status = 'old' )

c     read(60) ( pt(i), i=1,lspecies(Nspecies) )
      read(60) step
      Do i=1, lspecies(Nspecies)
         pt(i)= step
      End do
      close(60)

      return
      end

c     -------------------------------------------
      subroutine Write_Particle_Time ( FileName )
c     -------------------------------------------
c
c     purpose: write particle time moments in file FileName
c
      include 'a_tree.h'
      include 'a_control.h'
      DIMENSION         wspecies(nspecies),lspecies(nspecies)
      EQUIVALENCE    (wspecies(1),extras(1)),
     +                               (lspecies(1),extras(11))
      character FileName*(*)

      open ( 60 , file = FileName , form = 'unformatted' )

      write (60) ( pt(i), i=1,lspecies(Nspecies))
      close ( 60 )
      return
      end



c     -----------------------------------
      subroutine Read_Particles_Binary_o ()
c     -----------------------------------
c
c     purpose: opens files with control information and particle data
c              reads in particle coordinates and momenta
c
      include 'a_tree.h'
      include 'a_control.h'
      DIMENSION         wspecies(nspecies),lspecies(nspecies)
      EQUIVALENCE    (wspecies(1),extras(1)),
     +                               (lspecies(1),extras(11))
 
      ngrid = ng

      open ( 2,file = 'Result.DAT')
      open ( 3 ,file ='PMcrd.DAT', form = 'unformatted')

c.... read control information and check whether it has proper structure
      Ocurv = 0
c      lspecies(1) =NROW**3
c      wspecies(1) =NROW**3/FLOAT(ngrid**3)
      read      (3) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,Nspecs,Nseed,Om0,Oml0,hubble,Wp5
     &                   ,Ocurv,extras
      write (*,100) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  EKIN,EKIN1,EKIN2,
     &                  NROWC,NGRID,NRECL,Om0,Oml0,hubble
      write (2,100) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  EKIN,EKIN1,EKIN2,
     &                  NROWC,NGRID,NRECL,Om0,Oml0,hubble
 100  format (1X,'Header=>',A45,/
     &            1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     &            1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/
     &            1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I8,/
     &            1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,' Hubble=',f7.3)
      if(Nspecs.eq.0)Nspecs =1
      if ( Nspecs.ne. Nspecies ) then
        write (*,*)
     &      ' Nspecies in PARAMETER and in TAPE-FILE are different:'
         write (*,*) ' Nspecies=',Nspecies,' not equal ',Nspecs
         STOP
      endif
      If(np.lt. lspecies(Nspecies))Then
         write(*,*) ' Nspecies=',Nspecies,wspecies(1)
         write (*,*) ' Wrong number of particles !!! '
         write (*,*) ' It should be =',lspecies(Nspecies),' (lspecies)'
         write (*,*) ' but it was requested =',np,' (a_setup.h)'
         STOP
      Endif 
c.... open particle file on disk
c     for Exemplar, SP: nacces = nbyte; for PC: nacces = nrecl 

      nbyte  = nrecl * 4
      nacces = nbyte / nbyteword
 
      open ( 1 , file = 'PMcrs0.DAT', access = 'direct',
     &	         status = 'unknown', recl = nacces      )
 
      rewind 3

         N_particles =lspecies(Nspecies)   ! Total number of particles
         Npages = (N_particles -1)/NPAGE +1
         N_in_last=N_particles -NPAGE*(Npages-1)
      write (*,*) ' Pages=',Npages,' Species=',Nspecies
      write (*,*) ' Nparticles=',N_particles 

      DO  IROW=1, Npages         ! Loop over particle pages
            In_page =NPAGE
            If(IROW.eq.Npages)In_page =N_in_last
            write (*,*)' Read page=',IROW,' file=',ifile,' N=',In_page
            iL = NPAGE*(IROW-1)
         CALL GETROW(IROW,1) ! read in a page of particles
         DO  IN=1, In_page          ! Loop over particles
                ip =IN+iL                     ! current particle number
c                  WPAR =iWeight(ip)   ! particles weight
c check rounding errors (REAL*8 <=> REAL*4)
                x(ip) = xpar(in)
                IF(x(ip).LT.1.) x(ip) = x(ip) + yn
                IF(x(ip).GT.xn) x(ip) = x(ip) - yn
                y(ip) = ypar(in)
                IF(y(ip).LT.1.) y(ip) = y(ip) + yn
                IF(y(ip).GT.xn) y(ip) = y(ip) - yn
                z(ip) = zpar(in)
                IF(z(ip).LT.1.) z(ip) = z(ip) + yn
                IF(z(ip).GT.xn) z(ip) = z(ip) - yn
             If(  x(ip).ge.xn)Then
                x(ip) = x(ip) -1.e-4
                write (*,500) x(ip),ip,1
             endif
             If(  y(ip).ge.xn)Then
                y(ip) = y(ip) -1.e-4
                write (*,500) y(ip),ip,2
             endif
             If(  z(ip).ge.xn)Then
                z(ip) = z(ip) -1.e-4
                write (*,500) z(ip),ip,3
             endif
 500           format(' fixed boundary. It is now=',g14.7,i9,i2)

                vx(ip) = vxx(in)
                vy(ip) = vyy(in)
                vz(ip) = vzz(in)
         Enddo
      Enddo

      close (1)
      close (2)
      close (3)

      return
      end

c     ------------------------------------------------------------
      subroutine Write_Particles_Binary1_o ( FileName1 , FileName2 )
c     ------------------------------------------------------------
c
c     purpose: writes control information and particles to the specified 
c              files (this routine is to be used when nspec = 1)
c
c     input  : FileName1 - C3CRD*; FileName2 - C3crs0*
c

      include 'a_tree.h'
      include 'a_control.h'
      DIMENSION         wspecies(nspecies),lspecies(nspecies)
      EQUIVALENCE    (wspecies(1),extras(1)),
     +                               (lspecies(1),extras(11))

      character FileName1*(*) 
      character FileName2*(*)

c     for Exemplar, SP: nacces = nbyte; for PC: nacces = nrecl 

      nbyte  = nrecl * 4
      nacces = nbyte /nbyteword

      open (4 , file = FileName1,
     &           form = 'UNFORMATTED' , status = 'UNKNOWN')
 
      open (5 , file = FileName2 , access = 'DIRECT',
     &	         status = 'UNKNOWN', recl = NACCES)

c.... write header and control data

      write (4) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     &                   ,Ocurv,extras
      close (4)


         N_particles =lspecies(Nspecies)   ! Total number of particles
         Npages = (N_particles -1)/NPAGE +1
         N_in_last=N_particles -NPAGE*(Npages-1)
c      write (*,*) ' Pages=',Npages,' Species=',Nspecies
c      write (*,*) ' Nparticles=',N_particles

      DO  IROW=1, Npages         ! Loop over particle pages
            In_page =NPAGE
            If(IROW.eq.Npages)In_page =N_in_last
c            write (*,*)' Write page=',IROW,' file=',ifile,' N=',In_page
            iL = NPAGE*(IROW-1)
         DO  IN=1, In_page          ! Loop over particles
              ip =IN+iL                     ! current particle number
	        xpar(in) = x(ip) 
	        ypar(in) = y(ip) 
	        zpar(in) = z(ip) 
	         vxx(in)  =VX(ip)
	         vyy(IN)  =VY(ip)
	         vzz(IN)  =VZ(ip)
         EndDo
         CALL WRIROW(IROW,5)
      EndDo

      close (5)

      return
      end

c     ----------------------------------
      subroutine WriRow ( irow , ifile )
c     ----------------------------------
c
c     input:  irow - number of record
c             ifile - number of i/o unit (0-99)
c     nrecl - length of ROW block in words
c
      include 'a_tree.h'
      integer irow , ifile 
      write (ifile , rec = irow ) recdat
      return
      end
 
c     ----------------------------------
      subroutine GetRow ( irow , ifile )
c     ----------------------------------
c
c     input:  irow - number of record
c             ifile - number of i/o unit (0-99)
c     nrecl - length of ROW block in words
c
      include 'a_tree.h'
      read  ( ifile , rec = irow ) recdat
      return
      end

cc     
c     -----------------------------------------------------------
      subroutine Dump_Particles ( x1,x2,y1,y2,z1,z2 )
c     -----------------------------------------------------------
c
c     purpose: dumps particle coordinates 
c     uses   :              
c
      include 'a_tree.h'

      open( 8 , file = 'part_dump' , status = 'unknown' )
      do 159 ic1 = 1 , np
        if ( x(ic1) .lt. x1 ) goto 159
        if ( x(ic1) .ge. x2 ) goto 159
        if ( y(ic1) .lt. y1 ) goto 159
        if ( y(ic1) .ge. y2 ) goto 159
        if ( z(ic1) .lt. z1 ) goto 159
        if ( z(ic1) .ge. z2 ) goto 159
        icell = iFindCell (MaxLevelNow,x(ic1),y(ic1),z(ic1))
        write(8,160) ic1,icell, iLv(icell), var(1,icell), pot(icell),
     &                x(ic1) , y(ic1) , z(ic1) , 
     &                vx(ic1) , vy(ic1) , vz(ic1), pdummy(ic1)
 159  continue 

 160  format( i7,1x,i7,1x,i2,2(1x,f14.5),6(1x,f10.6),1x,f14.5)
      close (8)

      return
      end

c     ---------------------
      subroutine TreeDump()
c     ---------------------
c
c     purpose: dumps particle coordinates 
c     uses   :              
c
      include 'a_tree.h'

      open ( 32 , file = 'art.dat', form = 'unformatted' )
      ioct = 0
      do Level = 1 , MaxLevelNow
        nLevel = iNOLL(Level)         ! get level boundaries for index array 
        iO     = iHOLL(Level) 
        do j = 1 , nLevel
          ioct = ioct + 1
          iSelect(iOct) = iO
          iO            = iOctLL1(iO)
        enddo
      enddo
c      write(*,*) 'nOct = ',nOct, ' iOct =', iOct
      write(32) iOct, iOctFree, MinLevel, MaxLevel, MaxLevelNow
      write(32) ( iSelect(i), i=1,iOct )
      do i = 1 , iOct
        iO = iSelect(i)
        write(32) (iOctPs(iO,j),j=1,3)
        write(32) iOctLv(iO)
        write(32) (iOctNb(iO,j),j=1,6)
        write(32) iOctPr(iO)
        write(32) iOctLL1(iO)
        write(32) iOctLL2(iO)
c        write(32) ( (iOctPs(iO,j),j=1,3), i=1,iOct )
c        write(32) (  iOctLv(iO), i=1,iOct )
c        write(32) ( (iOctNb(iO,j),j=1,6), i=1,iOct )
c        write(32) (  iOctPr(iO), i=1,iOct )
c        write(32) (  iOctLL1(iO), i=1,iOct )
c        write(32) (  iOctLL2(iO), i=1,iOct )
      enddo
      write(32) (  iHOLL(i), i=MinLevel, MaxLevelNow )
      write(32) (  iNOLL(i), i=MinLevel, MaxLevelNow )
      write(32) ncell
      write(32) ( iOctCh(i) , i=1,ncell )
c
      close ( 32 )

      return
      end

c     -----------------------
      subroutine Init_Control_Files ()  
c     -----------------------
c
c     initializes control_files for queue management
c 
 
      include 'a_constant.h'
      OPEN(UNIT=19,file = 'ERROR_RUN', status = 'unknown')
      write(19,*) nil, 
     +    ' : NO errors, ART is running'
      close(19) 
      OPEN(UNIT=19,file = 'STOP_RUN', status = 'unknown')
      write(19,*) nil, ' : ART is running'
      close(19)

      return
      end
      

c     -----------------------
      subroutine check_remaining_time () 
c     -----------------------
c
c     check remaining time 
c 

      include 'a_control.h'
      Common /check/ time_start, time, continue_run, istepp
      Logical  continue_run

      delta_time= seconds() - time    ! time of last step 
      time = seconds() - time_start
      rest = 60.*queue_time - time

      open(777,file='queue_timing', access = 'append' )
c      open(777,file='queue_timing', position = 'append' )
      write(777,7777) istep,istepp,
     +             time/60.,delta_time/60.,rest/60.
 7777 FORMAT('istep = ',I4,' istepp = ',I4,' time = ',F8.2, 
     +        ' delta_time = ',F7.2,
     +        ' rest = ',F8.2)

      IF(rest.LT.delta_time*1.1) THEN
        Continue_Run= .false. 
        write (*,*) ' istep=',istep,' iostep=',iostep
        write (777,7778)  istep, iostep, delta_time, rest
 7778   FORMAT(' Continue_Run= .false.: istep =',I4,
     +                  ' iostep =', I3,
     +                  ' delta_time [s] = ', F10.1,
     +                  ' rest [s] = ', F10.1 )

c write step if not yet written
        if ( mod(istep,iostep) .ne. 0 ) then    !ELSE already written
          write (777,*) 'Write Particles and Particle_Time'
          write (*,*) ' Write Particles'
          call Write_Particles_Binary1_o ( 'PMcrd.DAT      ' ,
     &                                 'PMcrs0.DAT      ')
          write (*,*) ' Write pt'
          call Write_Particle_Time     ( 'pt.dat      ' )
        endif         
      ENDIF
      close(777)

      return
      end










