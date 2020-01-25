c     =====================================================================
c                                                                         .
c       ART Version 3: Advancing control variables to the next time step  .
c                                                                         .
c                by Andrey Kravtsov and Anatoly Klypin (1997)             .
c                                                                         .
c     =====================================================================

c     ---------------------
      subroutine Advance ()
c     ---------------------
c
c     advances time variables, 
c     computes kinetic and potential energies
c     and energy conservation error
c
      include 'a_tree.h'
      include 'a_control.h'

      double precision epot_,dspl
      double precision ekin_(nspecies)                   ! kinetic energy 
      dimension Ndisplace(0:50)

      ahalf = aexpn + 0.5*astep
      dconst = astep/ sqrt (ahalf * (Om0+Oml0*ahalf**3)) /
     &     ahalf  ! constant to convert velocity to
                                ! displacement in units of 1 of zero-level cell
      dspl     = 0.
      ndspl   = 0
      dspl_max =0.
      do i=0,50
         Ndisplace(i) =0
      enddo 
      do ic1 = 1 , nspecies
        ekin_(ic1) = zero 
        epot_      = zero 
      enddo

      do ic0 = 1 , nspecies 
        do ic1 = nsp(ic0,1) , nsp(ic0,2)
           v2   = vx(ic1)**2 + vy(ic1)**2 + vz(ic1)**2
          ekin_(ic0) = ekin_(ic0) + wpar(ic0) * v2
          epot_ = epot_ + wpar(ic0)*pdummy(ic1)
          vv = dconst*sqrt(v2)
          iv = MIN(INT(10.*vv),50)            ! units of 1/10 of a cell
          Ndisplace(iv) =Ndisplace(iv) +1
          dspl = dspl + vv
          ndspl = ndspl +1
          dspl_max =max(vv,dspl_max)
        enddo
      enddo
      nn =0
      iv   =0
       do i=50,0,-1
         nn =nn +Ndisplace(i)
         if(nn.gt.100)Then
            iv =i
            goto 5
         endif 
      enddo 
 5    write (*,100) ndspl,dspl/ndspl,dspl_max,iv,nn
 100  format(i8,' mean displacement=',g11.3,' max=',g11.3,
     &         ' 100_fast=',i6,' nfast=',i5,' units: 1 cell')
      write (*,110) Ndisplace
 110  format(5x,' N with given displacement: 1/10 of cell',/6(10i6/))
      enkin = zero
      do ic0 = 1 , nspecies 
        ekin_(ic0) = ekin_(ic0) / 2.0 / (aexpn + 0.5*astep)**2 
        enkin      = enkin + ekin_(ic0)
      enddo
        ENpot = epot_ / 2.0
 
      istep = istep + 1
      aexpn = aexpn + astep
c
c.... energy conservation control
c
      IF ( istep .eq. 1 ) THEN
        Ekin1 = Ekin
        Ekin2 = 0.
        Ekin  = ENkin
        au0   = aexp0 * ENpot
        aeu0  = aexp0 * ENpot + aexp0 * (Ekin + Ekin1) / 2.0
        TINTG = 0.

        write(50,*) HEADER
        write(*,40) istep , aexpn , Enkin , ENpot , au0 , aeu0
        write(50, *) ' istep ' , ' aexpn ' , 
     &               ' Ekin ' , ' Ekin_s ' , ' Ekin_l ' , 'PE'
        write(50,40) istep , aexpn , Enkin , 
     &               (Ekin_(ic0),ic0=1,nspecies) , Enpot
        write(50, *) ' istep ' , ' aexpn ' , ' MaxLevel ' , ' ncell ',
     &               ' Ekin ' , ' Ekin1 ' , ' Ekin2 ' , 
     &               '  PE  ' , ' TINTG '

      ELSE
        Ekin2 = Ekin1
        Ekin1 = Ekin
        Ekin  = ENkin
        TINTG = TINTG + astep * (Ekin1 + (Ekin - 2.*Ekin1 + Ekin2)/24.)
        error = ((aexpn - astep) * ((Ekin  + Ekin1) / 2.0 + ENpot) - 
     &            aeu0  + TINTG) / ((aexpn - astep) * ENpot)

        ncell = noct * nchild + ncell0

        call Get_MaxLevelNow ()

        write  (*,50)  istep , aexpn , MaxLevelNow , ncell , 
     &                  Enkin ,  (Ekin_(ic0),ic0=1,nspecies) , 
     &                  Enpot , error , TINTG
        write (50,50)  istep , aexpn , MaxLevelNow , ncell , 
     &                  Enkin ,  (Ekin_(ic0),ic0=1,nspecies) ,
     &                  Enpot , error , TINTG

      ENDIF

c.... check if format conforms with the number of species
c     namely check number of e12.4

40    format (i5 , 1x , f10.4 , 1x , 4(1x,e12.4) )
50    format (i5 , 1x ,  f7.5 , 1x , i2 , 1x ,  i10 , 7(1x,E13.5) )


      return
      end

