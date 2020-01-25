C __________	Merge 2 (or more) RODIN ICs _____________
C    January 2020         S. Roca-Fabrega (sroca01@ucm.es) 
C		                     
C          Scales internal PM coordinates and velocities 
C          to  coordinates in Mpc/h and km/s
C          In PM the coordinates are in the range 1 - (NGRID+1)
C                     velocities are P = a_expansion*V_pec/(x_0H_0)
C                     where     x_0 = comoving cell_size=Box/Ngrid
C                                    H_0 = Hubble at z=0
C                   
C	     NROW = number of particles in 1D
C	     NGRID= number of cells        in 1D

C ------------

C ------------

       INCLUDE 'PMparART.h'
c       REAL     INPUT
       Character  FileASCII*50
       character sntime*6
c       character sntime*5  ! Ony for HART simulations (including hydro)
       character tail*4
       character tailasc*4
       character tailtip*4
       character head1*6
       character head2*7
       dimension xpart(10000000,11),ypart(10000000,11),zpart(10000000,11
     +),vxpart(10000000,11),vypart(10000000,11),vzpart(10000000,11),nums
     +(11),x(100000000),y(100000000),z(100000000)
      dimension vxx(100000000),vyy(100000000),vzz(100000000),
     &lspecies1(10)
      Box  =1.0 !(box size Mpc/h)                                                                                                               
      scrip=1 !Output step (allparticles/scrip will be writen in the output tipsy file
      BoxV =Box*100.    ! Box size in km/s                                                                                                 
      tail='.DAT'
      tailasc='.asc'
      tailtip='.tip'
      head1='PMcrda'
      head2='PMcrs0a'
      xmax =-1.e+9
      xmin = 1.e+9
      vmax =-1.e+9
      vmin = 1.e+9


C...................................................................
C			Read data and open RODIN files

C First we are going to read the satellite particles and set the 
C satellite position and velocities with respect to the mai galaxy
C The main galaxy will reside in the center of the simulated box (0,0,0)
C The main galaxy halo, bulge, disk will have (0,0,0) (0,0,0) and (0,Vcirc,0) mean velocities, respectively


c     READING RODIN OUTPUT HEADER FILE (SATELLITE)

      open(3,file='PMcrd.DAT',form='unformatted',status='unknown
     &')
      read      (3,err=10,end=10) HEADER,
     &              AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &              TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &              NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     &              ,Ocurv,extras
      ScaleV = BoxV/AEXPN/NGRID  ! scale factor for Velocities                                                                            \
                                                        
      ScaleC = Box*aexpn/NGRID         ! scale factor for Coordinates 
 
 100     FORMAT(1X,'Header=>',A45,/
     +      1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     +            1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',E12.3,/
     +            1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I6,/
     +            1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,' Hubble=',f7.3,/
     +            1x,' Omega_curvature=',F7.3)

 10   npaget    = nrow**2    ! # of particles in a page
      nrecli    = npaget * 6 ! length of particle row in words
      nbyte  = nrecli * 4
      nbyteword=1!  defines length of direct-access record, usually 4 for HART, 1 for the rest (ART)
      nacces = nbyte / nbyteword
      xn=float(ngrid)+1.-1.E-7
      yn=float(ngrid)

c READING RODIN OUTPUT MAIN FILE

       open ( 1 , file ='PMcrs0.DAT', access = 'direct',
     &           status = 'unknown', recl = nacces      )


         N_particles =lspecies(Nspecies)   ! Total number of particles                                                                    
         do i=1,nspecies
          write(*,*)'Particulas de la especie',i,' = ',lspecies(i),'weig
     &ht =', wspecies(i)
         enddo
         Npages = (N_particles -1)/NPAGE +1
         N_in_last=N_particles -NPAGE*(Npages-1)
      write (*,*) ' Pages=',Npages,' Species=',Nspecies
      write (*,*) ' Nparticles=',N_particles

      DO  IROW=1, Npages         ! Loop over particle pages                                                                                
            In_page =NPAGE
            If(IROW.eq.Npages)In_page =N_in_last
            iL = NPAGE*(IROW-1)
         CALL GETROW(IROW,1) ! read in a page of particles                                                                                 
         DO  IN=1, In_page          ! Loop over particles                                                                                  
                ip =IN+iL                     ! current particle number                                                                    
                  WPAR =iWeight(ip)   ! particles weight                                                                                  
c check rounding errors (REAL*8 <=> REAL*4)                                                                                                
                x(ip) = xpar(in)
                y(ip) = ypar(in)
                z(ip) = zpar(in)

             If(  x(ip).ge.xn)Then
                x(ip) = x(ip) -1.e-4
                write (*,*) x(ip),ip,1
             endif
             If(  y(ip).ge.xn)Then
                y(ip) = y(ip) -1.e-4
                write (*,*) y(ip),ip,2
             endif
             If(  z(ip).ge.xn)Then
                z(ip) = z(ip) -1.e-4
                write (*,*) z(ip),ip,3
             endif

                vxx(ip) = vx(in)
                vyy(ip) = vy(in)
                vzz(ip) = vz(in)
           X(ip)  =ScaleC* (X(Ip)-1.)*1000.0d0/hubble/(2.87e4)       
           Y(ip)  =ScaleC* (Y(Ip)-1.)*1000.0d0/hubble/(2.87e4)
           Z(ip)  =ScaleC* (Z(Ip)-1.)*1000.0d0/hubble/(2.87e4)
           Vxx(ip)=ScaleV* VXx(Ip)/691.0d0
           Vyy(ip)=ScaleV* VYy(Ip)/691.0d0
           Vzz(ip)=ScaleV* VZz(Ip)/691.0d0

c     Adding positions and velocities of the satellite with respect
c    to the center of the main galaxy

C  SATELLITE positions and velocities need to be indicated in the
C   "PMparART.h" file (lines 25 to 30):

           X(ip) = (X(ip)+xsat)/ScaleC/1000.0d0*hubble*2.87e4+1.
           Y(ip) = (Y(ip)+ysat)/ScaleC/1000.0d0*hubble*2.87e4+1.
           Z(ip) = (Z(ip)+zsat)/ScaleC/1000.0d0*hubble*2.87e4+1.
           Vxx(ip) = (Vxx(ip)+vxsat)/ScaleV*691.0d0
           Vyy(ip) = (Vyy(ip)+vysat)/ScaleV*691.0d0
           Vzz(ip) = (Vzz(ip)+vzsat)/ScaleV*691.0d0

         Enddo
      Enddo

      hfact=hubble
      nspecies0=nspecies
      do i=1,nspecies0
         lspecies0(i)=lspecies(i)
      enddo

      close(1)
      close(3)

c     Here we read the main simulation

      write(*,*)'Reading a snapshot of the main simulation'


C You can add the satellite not only in the ICs file of the main galaxy 
C but also to a time when the main galaxy system has already been self-relaxed

      write(*,*)'Write snapshot time in scale factor units in the form
     +: 0.6000'

C 0.6000 is the default ICs time (0Myrs evolution). BE CAREFUL no evolution
C means no relaxation of the system!!

      read(*,*)sntime
      
      open(3,file=head1//sntime//tail,form='unformatted',status='unknown
     &')
      read      (3,err=10,end=10) HEADER,
     &              AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &              TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &              NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     &              ,Ocurv,extras
      close(3)

      open ( 1 , file =head2//sntime//tail, access = 'direct',
     &           status = 'unknown', recl = nacces      )


         N_particles =lspecies(Nspecies)   ! Total number of particles                                                                    
         do i=1,nspecies
          write(*,*)'Particulas de la especie',i,' = ',lspecies(i),'weig
     &ht =', wspecies(i)
         enddo
         Npages = (N_particles -1)/NPAGE +1
         N_in_last=N_particles -NPAGE*(Npages-1)
      write (*,*) ' Pages=',Npages,' Species=',Nspecies
      write (*,*) ' Nparticles=',N_particles

      DO  IROW=1, Npages         ! Loop over particle pages                                                                                
            In_page =NPAGE
            If(IROW.eq.Npages)In_page =N_in_last
            iL = NPAGE*(IROW-1)
         CALL GETROW(IROW,1) ! read in a page of particles                                                                                 
         DO  IN=1, In_page          ! Loop over particles                                                                                  
                ip =IN+iL                     ! current particle number                                                                    
                  WPAR =iWeight(ip)   ! particles weight
c check rounding errors (REAL*8 <=> REAL*4)                                                                                                
                xc(ip) = xpar(in)
                yc(ip) = ypar(in)
                zc(ip) = zpar(in)

             If(  xc(ip).ge.xn)Then
                xc(ip) = xc(ip) -1.e-4
                write (*,*) xc(ip),ip,1
             endif
             If(  yc(ip).ge.xn)Then
                yc(ip) = yc(ip) -1.e-4
                write (*,*) yc(ip),ip,2
             endif
             If(  zc(ip).ge.xn)Then
                zc(ip) = zc(ip) -1.e-4
                write (*,*) zc(ip),ip,3
             endif

               vxc(ip) = vx(in)
                vyc(ip) = vy(in)
                vzc(ip) = vz(in)
         Enddo
      Enddo

      close(1)

c     We now add particles of the satellite to the main simulation
      do i=1,nspecies
      lspecies(i)=lspecies(i)+lspecies0(i)
      enddo
      write(*,*)'WARNING!!!! All incoming particles need to have the 
     +same mass as the ones of the same specie in the main simulation'

         N_particles =lspecies(Nspecies)   ! Total number of particles                                                                    
         ip=n_particles
         do i=1,nspecies
          write(*,*)'Particulas de la especie',i,' = ',lspecies(i)-
     +lspecies(i-1),'weight =', wspecies(i)
         enddo

      do i=lspecies(1)-lspecies0(1)+1,lspecies(1)
      xc(i)=x(i-lspecies(1)+lspecies0(1))
      yc(i)=y(i-lspecies(1)+lspecies0(1))
      zc(i)=z(i-lspecies(1)+lspecies0(1))
      vxc(i)=vxx(i-lspecies(1)+lspecies0(1))
      vyc(i)=vyy(i-lspecies(1)+lspecies0(1))
      vzc(i)=vzz(i-lspecies(1)+lspecies0(1))
      enddo

      do j=2,nspecies
      do i=lspecies(j)-(lspecies0(j)-lspecies0(j-1))+1,lspecies(j)
      xc(i)=x(i-lspecies(j)+lspecies0(j)-lspecies0(j-1))
      yc(i)=y(i-lspecies(j)+lspecies0(j)-lspecies0(j-1))
      zc(i)=z(i-lspecies(j)+lspecies0(j)-lspecies0(j-1))
      vxc(i)=vxx(i-lspecies(j)+lspecies0(j)-lspecies0(j-1))
      vyc(i)=vyy(i-lspecies(j)+lspecies0(j)-lspecies0(j-1))
      vzc(i)=vzz(i-lspecies(j)+lspecies0(j)-lspecies0(j-1))
      enddo
      enddo

c     Writting in new IC files

c     Write header PMcrd

      CALL RDTAPE
c      OPEN(UNIT=9,FILE=npath//'PMcrd_ini.DAT',form='unformatted')
c      write (9)             ! this clears old header file
      CALL WRTAPE
c      close (9)
c     write pt.dat file: tme-step for particles

      If(Nspecies.eq.0)Then
         Nparticles =lspecies(1)
      Else
         Nparticles =lspecies(Nspecies)
      EndIf
      Do ic1 =1,Nparticles
         xpt(ic1) =astep
      Enddo
      open ( 60 , file = npath//'pt_ini.dat' , form = 'unformatted' )
      write(60) (xpt(ic1),ic1=1,Nparticles) 

c     Write in main data file PMcrs

      CALL WriteData

      END


      subroutine GetRow ( irow , ifile )
c     ----------------------------------
c
c     input:  irow - number of record
c             ifile - number of i/o unit (0-99)
c     nrecl - length of ROW block in words
c
      include 'PMparART.h'
      read  ( ifile , rec = irow ) recdat
      return
      end

c     -----------------------------------

      SUBROUTINE WRTAPE
C----------------------------------------------
      INCLUDE 'PMparART.h'
C     write header and control data
      WRITE (9) HEADER,
     +           AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +           TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     +           NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5,
     +           Ocurv,extras
      REWIND 9
      RETURN
      END

C---------------------------------------------
C                       Write current data to  disk/tape
      SUBROUTINE WriteData
C----------------------------------------------
      INCLUDE 'PMparART.h'

      Ibuff =0
      KROW =0
      Do i=1,lspecies(Nspecies) !Over all the particles
             Ibuff              = Ibuff +1

             XPAR(Ibuff) = xc(i)
             YPAR(Ibuff) = yc(i)
             ZPAR(Ibuff) = zc(i)
             VX(Ibuff)   = vxc(i)
             VY(Ibuff)   = vyc(i)
             VZ(Ibuff)   = vzc(i)

             If(Ibuff.ge.NPAGE)Then
                KROW = KROW +1
                CALL WRIROW(KROW,1)
                Ibuff     =0
             EndIf
      EndDo

          If(Ibuff.ne.0)Then    ! write last incomplete page
                KROW = KROW +1
                write (*,*) ' Write page=',KROW,' i_buff=',Ibuff
                CALL WRIROW(KROW,1)
                Ibuff     =0
          EndIf
 38       format(i6,1x,7(1x,f8.4))

      RETURN
      END
C--------------------------------------------------
C                             Write current PAGE of particles (x,v) to
C                             disk
C                             NRECL - length of ROW block in words
      SUBROUTINE WRIROW(IROW,Ifile)
C--------------------------------------------------
      INCLUDE 'PMparART.h'
        WRITE (20+Ifile,REC=IROW) RECDAT
      RETURN
      END

C---------------------------------------------------
C                                  Read  current data from disk/tape,
C                                  Open files
C                                  Nrecl is the number of values in a
C                                  record
C                                  Npage is the number of particles in a
C                                  page
      SUBROUTINE RDTAPE
C---------------------------------------------------
      INCLUDE 'PMparART.h'

      Character  Hd*5,Tail*8,Nm*44
      Hd  ='PMcrs'
      Tail='_ini.DAT'
C                                     Open file on a tape
      OPEN(UNIT=9,FILE=npath//'PMcrd_ini.DAT',
     +                FORM='UNFORMATTED', STATUS='UNKNOWN')
C                                 Read control information
C                                 and see whether it has proper
C                                 structure
      READ  (9,err=10,end=10) HEADER,
     +                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     +                  NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     +                   ,Ocurv,extras
      Ocurv =0.
      WRITE (*,100) HEADER,
     +                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +                  EKIN,EKIN1,EKIN2,
     +                  NROWC,NGRID,NRECL,Om0,Oml0,hubble,
     +                  Ocurv,extras
      WRITE (16,100) HEADER,
     +                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +                  EKIN,EKIN1,EKIN2,
     +                  NROWC,NGRID,NRECL,Om0,Oml0,hubble,
     +                  Ocurv,extras
100   FORMAT(1X,'Header=>',A45,/
     +           1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.4,/
     +            1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/
     +            1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I6,/
     +            1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,' Hubble=',f7.3,/
     +            1x,' Omega_curvature=',F7.3)
      IF(NROW.NE.NROWC) THEN
         WRITE (*,*)
     +            ' NROW in PARAMETER and in TAPE-FILE are different' 
      ENDIF
      IF(NGRID.NE.NGRIDC) THEN
         WRITE (*,*)
     +           ' NGRID in PARAMETER and in TAPE-FILE are different:'
         write (*,*) ' Ngrid=',NGRID,' NgridC=',NGRIDC
      ENDIF
C                                         Open work file on disk 
 10   NBYTE = NRECL*4
      NACCES= NBYTE!/4           !the four is for athena    
                          
          DO ifile=0,0
         Nm =Hd//Char(ifile+48)//Tail
         iun=21+ifile
c         write (*,*) ' File>',ifile,' unit=',iun,' Name=',Nm
         OPEN(UNIT=iun,FILE=npath//Nm,ACCESS='DIRECT',
     +                 STATUS='UNKNOWN',RECL=NACCES)
      EndDo

      REWIND 9
      RETURN
      END

