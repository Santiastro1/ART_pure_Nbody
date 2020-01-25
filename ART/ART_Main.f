C     ====================================================
c
c          Adaptive Refinement Tree (ART) N-body solver
c
c                 Version 3 - February 1997
c
c       Andrey Kravtsov, Anatoly Klypin, Alexei Khokhlov
c
c     ====================================================

      program ART

      include 'a_tree.h'
      include 'a_control.h'
      Common /AuxTime/t1,t2,t3,t4,t5
      Common /check/ time_start, time, continue_run, istepp

      logical Continue_Run

      integer nMove
      integer mtot
      iStepp =0
      write (*,*) ' Start '
c.....
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( ic1)
      do ic1 = 1 , mcell
         var(1,ic1) = -1.0
         var(2,ic1) = -1.0
         ref(ic1)   = zero
         pot(ic1)   = zero
      enddo
c     EndDo
c.... read in initial conditions, initialize variables and linked lists
      sss =seconds()
      write (*,*) ' Read Particles '
      call Read_Particles_Binary_o ()
      call Read_Particle_Time    ()
      call Init_Parameters       ()
      call Init_Control_Files    ()
      Continue_Run = .true. 

      time_start = seconds()
      time=0.

      DO WHILE ( Continue_Run )      ! main loop
        t1 =0.
        t2 =0.
        t3 =0.

        dummy =seconds()  

        call Init_Arrays           ()
        call Init_Tree             ()
        call LL_Construct          () ! construct linked list of particles
        call Open_Log              ()
c.... rebuild refinements if code is restarted     
        do ic1 = 1 , 10 
           CPU(ic1) = 0.0
        enddo
        dummy =seconds()
      
        mtot  = 1
        ncell = noct * nchild + ncell0
        if ( istep .ge. 0 ) then 
           do while ( (mtot .ne. nil) )
              call Get_MaxLevelNow ()
              ncell = noct * nchild + ncell0

              call Timing ( 1 , -1 )
              call Assign_Density ( MaxLevelNow , MaxLevelNow )
              call Timing ( 1 ,  1 )
              write(*,*) ' ncells=',ncell ,' MaxLevel=', MaxLevelNow
              call Timing ( 8 , -1 )

              call Modify ( mtot, MaxLevelNow, MaxLevelNow, MaxLevel, 1)
              call Timing ( 8 ,  1 )
              ncell = noct * nchild + ncell0

           enddo
        endif
c        CALL ViewDump (MaxLevelNow ,ww , 1 )
c        CALL ViewDump (MaxLevelNow ,ww , 2 )
c        CALL ViewDump (MaxLevelNow ,ww , 3 )
c        STOP
        t2 =t2 + seconds()
c     DO WHILE ( Continue_Run )

        call Get_MaxLevelNow ()
        do Level = MinLevel , MaxLevelNow 
          call LL_Update ( Level , MinModify , MaxModify )
        enddo
        t3 =t3 +seconds()

        call Move ()

c....   advance time variables, compute energies

        call Timing ( 10 ,  -1 )
        call Advance
        call Timing ( 10 ,  1 )

c....   write timing 
        Continue_Run = .true. 
        call WriteTiming ()
        call Run_Output_o  ()
        call check_remaining_time ()
        iStepp =iStepp +1
c
c      if ( aexpn .ge. 1.0001 .or.iStepp.ge.10000)Continue_Run= .false. 
c        Continue_Run= .false.

        if(aexpn .ge. 1.0001) then 
           OPEN(UNIT=19,file = 'STOP_RUN', status = 'unknown')
           write(19,'(A)') ' 1   STOP'
           Continue_Run= .false.           
           close(19)
        endif
      call Close_Log ()
      ENDDO                     ! main loop ends


      call Open_Log  ()
      call Write_Particles_Binary1_o ( 'PMcrd.DAT      ' ,
     &                                 'PMcrs0.DAT      ')
      call Write_Particle_Time     ( 'pt.dat      ' )
      
      STOP
      END

c     ---------------------------------------
      subroutine ViewDump ( Lev , zc , iDir )
c     ---------------------------------------
c         Dump data for visualization
c
c
      character*50 plotname,jobname
      character*5 fstep
      parameter ( nvars = 8 )
      character*10 varnames(nvars)
      real CellPos(3)
      integer iOD(3,2)
      data iOD / 2 , 1 , 1 , 
     &           3 , 3 , 2   / 
c
      include 'a_tree.h'
      include 'a_control.h'
c
      call Get_MaxLevelNow ()
       call Assign_Density ( MinLevel , MaxLevelNow )
c
      iProlongFlag = 1
      call Solve_Poisson  ( MinLevel , Lev , iProlongFlag )    
c
      hsmall = 0.7
      Box = 2.e3/hsmall*AEXPN
      write (*,*) ' hsmall =',hsmall,' a=',AEXPN 
      Xscale = Box/ng
      write (*,*) ' Box =',Box,' kpc',' Cell=',Xscale
      rho_dm  = -100000.
      phi_min = 100000.
      numcell = 0
          call Select_Cells ( Lev , nLevel ) 
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              ic = icell + ic2
                call Ps ( ic , CellPos(1) , CellPos(2) , CellPos(3) ) 
                  numcell = numcell + 1
                  If(rho_dm .lt.var(1,ic))Then
                     rho_dm = var(1,ic)
                     x1 = CellPos(1)
                     y1 = CellPos(2)
                     z1 = CellPos(3)
                  EndIf
                  If(phi_min .gt.pot(ic))Then
                     phi_min = pot(ic)
                     x2 = CellPos(1)
                     y2 = CellPos(2)
                     z2 = CellPos(3)
                  EndIf
            enddo ! ic2
          enddo ! ic1            
          write (*,*) ' Max density =',rho_dm,' x=',x1,y1,z1,Lev
          write (*,*) ' Min potentl =',phi_min,' x=',x2,y2,z2,Lev
c
      jobname = 'File'
      if ( iDir .eq. 1 ) then 
         plotname  = 'Dump_x.v'
         zc = x2
      endif
      if ( iDir .eq. 2 ) then 
         plotname  = 'Dump_y.v'
         zc = y2
      endif
      if ( iDir .eq. 3 ) then 
         plotname  = 'Dump_z.v'
         zc = z2
      endif
c
      nlayers = 0
      imoviestep =0
      open ( 25 , file = plotname )
      write (unit=25,fmt='(a10,i10,1pe12.4,i10)')
     &            jobname1,imoviestep,AEXPN, nlayers
      write (unit=25,fmt='(''zones         1 1 1'')')
c
      nfuncs=3
      varnames(1)='Rho_dm'
      varnames(2)='phi'
      varnames(3)='lev'
      write ( unit = 25 , fmt = '(i3,20(1x,a10))' )
     &               nfuncs,(varnames(i),i=1,nfuncs)
      ifuncs = 0 
      write ( unit = 25 , fmt = '(i3)' ) ifuncs
c

c
c.... figure out number of cells first 
c      
      numcell = 0 
      do Level = MinLevel , Lev
        cs = CellSize(Level)/2.          ! halfsize of a cell
        IF ( Level .eq. MinLevel ) THEN 
          do ic = 1 , ncell0
            if ( iOctCh(ic) .eq. 0 .or. Level. eq. Lev ) then 
              call Ps ( ic , CellPos(1) , CellPos(2) , CellPos(3) ) 
              if ( abs(CellPos(iDir) - zc) .le. cs ) then 
                numcell = numcell + 1
              endif
            endif
          enddo ! ic
        ELSE
          nLevel = iNOLL(Level)
          call Select_Cells ( Level , nLevel ) 
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              ic = icell + ic2
              if ( iOctCh(ic) .eq. 0 .or. Level. eq. Lev ) then     
                call Ps ( ic , CellPos(1) , CellPos(2) , CellPos(3) ) 
                if ( abs(CellPos(iDir) - zc) .le. cs ) then 
                  numcell = numcell + 1
                endif
              endif
            enddo ! ic2
          enddo ! ic1            
        ENDIF
      end do
c.... write cell corners
c
c      write(*,*) 'numcell =',numcell
      write (25,*)  numcell*4

      do Level = MinLevel , Lev
        cs = CellSize(Level)/2.          ! halfsize of a cell
        IF ( Level .eq. MinLevel ) THEN 
          do ic = 1 , ncell0
            if ( iOctCh(ic) .eq. 0 .or. Level. eq. Lev ) then 
              call Ps ( ic , CellPos(1) , CellPos(2) , CellPos(3) ) 
              if ( abs(CellPos(iDir) - zc) .le. cs ) then 
                xi = CellPos(iOD(iDir,1))  ! X-coord of cell center
                yi = CellPos(iOD(iDir,2))  ! Y-coord of cell center
             write(25,'(1p4e13.6)')Xscale*(xi-cs),Xscale*(yi-cs),0.,0.
             write(25,'(1p4e13.6)')Xscale*(xi-cs),Xscale*(yi+cs),0.,0.
             write(25,'(1p4e13.6)')Xscale*(xi+cs),Xscale*(yi+cs),0.,0.
             write(25,'(1p4e13.6)')Xscale*(xi+cs),Xscale*(yi-cs),0.,0.
              endif
            endif
          enddo ! ic
        ELSE
          nLevel = iNOLL(Level)
          call Select_Cells ( Level , nLevel ) 
          do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              ic = icell + ic2
              if ( iOctCh(ic) .eq. 0 .or. Level. eq. Lev )
     +             then     
              call Ps ( ic , CellPos(1) , CellPos(2) , CellPos(3) ) 
                if ( abs(CellPos(iDir) - zc) .le. cs ) then 
                  xi = CellPos(iOD(iDir,1))  ! X-coord of cell center
                  yi = CellPos(iOD(iDir,2))  ! Y-coord of cell center
               write(25,'(1p4e13.6)')Xscale*(xi-cs),Xscale*(yi-cs),0.,0.
               write(25,'(1p4e13.6)')Xscale*(xi-cs),Xscale*(yi+cs),0.,0.
               write(25,'(1p4e13.6)')Xscale*(xi+cs),Xscale*(yi+cs),0.,0.
               write(25,'(1p4e13.6)')Xscale*(xi+cs),Xscale*(yi-cs),0.,0.
                endif
              endif
            enddo ! ic2
          enddo ! ic1            
        ENDIF
      end do
c
c.... now write variables 
c      
      nzone = 1
      write (*,*)  numcell
      write (25,*)  numcell
      numcell = 0      
      a3 = aexpn**3
      a2 = aexpn**2

      do Level = MinLevel , Lev
        fjfact = 2.**Level*sqrt(2.0*pi**2/3.*2.5/aexpn)
        cs = CellSize(Level)/2.          ! halfsize of a cell
        cs3 = CellSize(Level)**3
        IF ( Level .eq. MinLevel ) THEN 
          do ic = 1 , ncell0
            if ( iOctCh(ic) .eq. 0 .or. Level. eq. Lev ) then 
              call Ps ( ic , CellPos(1) , CellPos(2) , CellPos(3) ) 
              if ( abs(CellPos(iDir) - zc) .le. cs ) then 
                numcell = numcell + 1
                np1 = numcell * 4 - 3
                np2 = np1 + 1
                np3 = np2 + 1
                np4 = np3 + 1
                rho_dm = var(1,ic) + 1.0

                write (unit=25,fmt='(4(i7,1x),i5,1x,8(e14.6,1x))')
     &                np1,np2,np3,np4,nzone,
     &                rho_dm,pot(ic),Float(Level)
              endif
            endif
          enddo ! ic
        ELSE
           nLevel = iNOLL(Level)
           call Select_Cells ( Level , nLevel ) 
           do ic1 = 1 , nLevel
             icell = iSelect(ic1)
             do ic2 = 0 , 7
               ic = icell + ic2
               if ( iOctCh(ic) .eq. 0 .or. Level. eq. Lev ) then     
                 call Ps ( ic , CellPos(1) , CellPos(2) , CellPos(3) ) 
                 if ( abs(CellPos(iDir) - zc) .le. cs ) then 
                   numcell = numcell + 1
                   np1 = numcell * 4 - 3
                   np2 = np1 + 1
                   np3 = np2 + 1
                   np4 = np3 + 1
                   rho_dm = (var(1,ic) / cs3 +1.) 
                 write (unit=25,fmt='(4(i7,1x),i5,1x,8(e14.6,1x))')
     &                np1,np2,np3,np4,nzone,
     &                rho_dm,pot(ic),Float(Level)
                 endif
               endif
             enddo ! ic2
           enddo ! ic1            
        ENDIF
      end do
      write (unit=25,fmt='(''interfaces   0'')')
      close ( 25 ) 
c

      write(*,*) 'Saved data for view...'
c

      return
      end
