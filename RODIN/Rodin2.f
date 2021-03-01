c New version by  S. Roca-Fàbrega
c modificado el 31 de mayo de 2011, corregimos el cálculo de kappa del halo, en la línea 2143
c tambien hemos modificado el calculo de la Vcircular (rr+delta/2) --> rr (linea 1725)
c modificado 16 de noviembre de 2011, corregimos la linea 2935, donde se asignava erroneamente
c el valor de iorb (que corresponde al número total de particulas) al numero de partíclas de la 
c primera especie cuando se le debe asignar el valor de partículas de la primera especie encontradas
c para el halo + para el disco y sumarle las del satelite, lspecies(1)+Nsat

C-----------------------------------------   
C    the disk functions are for a finite disk    
C---------------------- -------------------    
C     Initial Orbits for particles for Iconf= 
C     0:  NFW Halos 
C     1:  Exponential Disk     
C     2:  Exponential Disk + Hernquist Bulge  
C     3: Exponential Disk within a NFW halo  
C     4: Exponential Disk + Bulge inside a NFW halo  
C     Isat = 
C     0: Isolated Model  
C     1: A NFW Satellite   
C     Irigid: 1 satellites from cosmological simualtion will be rigid  ! If problem check Nspecies 
C Units:  mass =  Mvir/F(C)    
C               radius=           Rs    
C               time    = t_0   =1/sqrt(GMvir/(Rs^3F(C))   
C               velocity=V_0 = Rs/t_0  
C               density=rho_0 =Mvir/(F(C)Rs^3)   
cC              Distances in Mpc/h  
c              Masses in Solar masses/h   
c              Velocities in km/s          
c         
C    Hernquist  Units:  mass =  M = Mvir(1. + 1/C)**2    
C               radius=           a    
C               time    = t_0   =1/sqrt(GM/a^3)   
C               velocity=V_0 = a/t_0 = sqrt(GM/a)   
C               density=rho_0 =M/a^3    

 
      INCLUDE 'PMparameters.h'    
      INCLUDE 'Rodin2.h'   
      Character HeadTest 
      parameter(Isirtakis =0)! 0 -> Rodin, 1 Sirtakis    
   
      
      COMMON / KINENRG/ SKINE,SX,SY,SZ,SX2,SY2,SZ2     
      REAL*4        Mmx,AINTG,Mcl     
      External       Mmx    
       
             
     
      Nrad     =1000! 15000!1500!600! 10000 !700!600!1000          ! number of radii for initial conditions                 
c      write(*,*)'Nrad = ',Nrad 
      CALL InitValues(NBYTE,SCALEL)   ! cosmological model    
      
      Call Configuration      
      
      FCon  = F(Con)     
      aMass_vir = aMassh*coef(1) + aMassd*coef(2) + 
     *     aMassb*coef(3)+aMass_sat* coef(4)
      Rvir     = 442.97 * (aMass_vir/1.e+11/(Om0*Overdens))**0.3333 ! kpc/h     
      
      Rs   = Rvir/Con  
      Rsnfw = Rs      
      R_0       = Rs            ! in units of kpc/h   
      M_0       = aMassh/Fcon   ! in units of h^{-1}Msun, massscale   
      t_0       = 4.693d+11/sqrt(M_0/R_0**3) ! in   yrs     
      V_0       = 9.767d+8 *R_0/t_0 ! in units of km/sec  
      Vscale    = 2.08e-5*AEXPN*NGRID/Box*sqrt(M_0/R_0) ! ! velocity scale in km/s, 1e-5*sqrt(G*2.d33/3e21) =208.e-5 
C                                Set up radii       
      
         Vconv = 208.e-5*(sqrt(M_0/R_0)) 
c                write(*,*)'Vconv= ',Vconv,V_0,Vscale 
         xscale= 1.e-3 *R_0*NGRID/Box/AEXPN        
         Vcmax = 98.e-5*Sqrt(CON*M_0/Rvir)   
        
            write(*,*)'  Rvir     =',Rvir,' Vcmax = ',Vcmax  
            write(*,*)'  Rs(Kpc/h)=',Rs,  ' Vscale=',Vscale,
     &                 ' Xscale =',xscale   
            write(*,*)'  Mvir     =',aMass_vir,' Mhalo =',aMassh      
      
            Fhcon = FHernq(C_b )   ! Hernquist mass funtion  
            FHb = Fhcon  

C                  Define quantities to be stored in extras 
         diskmass = aMassd  
         halomass = aMassh  
         Rdisk = Rd  
         Cnfw = Con   
         Zsc = Zdkp   
         Qt = Q 
         Rtrunc = Xout     
         Caja = Box 
         Nbulbo = Nbulge     
c         write(*,*)Rdisk,Cnfw,Zsc,Qt,Rtrunc,Caja    
C       

         OPEN(UNIT=16,FILE=npath//'Results.DAT') 
         OPEN(UNIT=9,FILE=npath//'PMcrd.DAT',form='unformatted') 
      write (9)             ! this clears old header file
      CLOSE (9)
      CALL RDTAPE                  ! this opens files on disk

      Do i=1,Nspecies 
         lspecies(i)=0   
         wspecies(i)=0. 
      Enddo 
      write(*,*)'                RDTAPE done'  

      If(Iconf.gt.0)then ! no empty halo  
      CALL CDISKBIN      

      pMdi =aMassd/Ndr*FCon/aMassh !disk mass p/particle in halo units   
      Iorb =Ndr                  !Initialize counters  
      If(Iconf.le.2.)wspecies(1)=1. ! just 1 specie   
      
     
      If(Iconf.gt.2)write(*,*)'  Nparticles before halo generation=',Ndr !system with halo  
      If(Iconf.le.2)write(*,*)'  Number of disk particles= ',Ndr ! sustem with no halo  
      else   
         write(*,*)'  Nparticles before Halo generation=',Ndr   
         Id=Ndr  
         Id=Nd
c        pMdi =aMassd/Id*FCon/aMassh   
         pMdi = (aMassd/aMassh)/Id*Fcon     
c        Iorb =0  
         write(*,*)'  Inces of configuration Iconf= ',Iconf  
c          pause   
      Endif  
       
      write(*,'("  disk particle mass(Msun)=",g11.4,
     &          " in internal units= ",g11.4)')
     &          aMassd/Id,pMdi     

      If(Iconf.eq.2.or.Iconf.eq.4)then !Bulge       
           
c        Rs_b = xtrunc/C_b* R_0  ! Bulge Rs in halo units      
        
         xtrunc = Xout_b/Rs_b  ! in Rs_b units   
         Rslarge = R_0/Rs_b   ! conversion factor from  !xtrunc/C_b* R_0 !     
         Clarge = C_b    
         aMlarge =  aMassb     
         Nwant = (aMassb/aMassd)*Id ! same mass per particle        
c         write(*,*)' Nwant_bulge= ',Nwant,Id,aMassb,aMassd     
         wspecies(1) = 1.0       
         R0_b = (Rsnfw/ Rs_b)  
         V0_b = sqrt((aMlarge/aMassh)*Fcon*R0_b)      
c         write(*,*)' V0_b', V0_b  
c         V0_b = sqrt((aMlarge/aMassh)*Fcon/FHernq(C_b)*R0_b)  
         write(*,*)'                Bulge Initialization Starts' 
         CALL Initialb(Nwant,xtrunc) 
         write(*,*)'                Bulge generation is done'     
      endif  

      Ncount =0    
      Ibin =0    
      txcmh =0    
      tycmh =0    
      tzcmh =0     
      twh =0     
       
       If(Iconf.eq.0.or.Iconf.gt.2)then !single halo or system inside a halo     
          If(nspecies.gt.1)then  

             If(Iconf.eq.0)OPEN(68,file=npath//'trackgalaxy.dat',   
     &            status='unknown')      

             write(*,*)'              Halo Initialization Starts'   
             CALL Initialnew(pMdi) ! set  initial conditions for multiple masses   
             
          else 
             xtrunc = xout  
            Rslarge = Rsnfw  
             Nwant = (aMassh/aMassd)*Id  ! same mass per particle    
             wspecies(1) = 1.0  
             OPEN(68,file=npath//'trackgalaxy.dat',   
     &            status='unknown')      
             write(*,*)'              Halo Initialization Starts'                                              
             If(Isirt akis.eq.0)then  
 
              CALL Initialh(Nwant,xtrunc)         
                 
	     else if(Isirtakis.eq.1)then 
 
        OPEN(40,file=npath//'nfw_high.ascii',   
     &            status='unknown')  

c          OPEN(40,file='test.dat',   
c     &            status='unknown')  

c       Do il=1,10
c	READ(40,*)HeadTest 
c        WRITE(*,'(a)')HeadTest 
c       Enddo
        Nt = 0
 
	READ(40,*)idumb,am
11      READ(40,*,err=20,end=20)xj,yj,zj,vxj,vyj,vzj !Read data     

c11      READ(40,*,err=20,end=20)il,xj,yj,zj,vxj,vyj,vzj,am !Read data  
        Nt = Nt+1

        xc(Nt) = xj*.7/R_0 
        yc(Nt) = yj*.7/R_0     
        zc(Nt) = zj*.7/R_0 

        vxc(Nt) = vxj/V_0  
        vyc(Nt) = vyj/V_0
        vzc(Nt) = vzj/V_0
        peso(Nt) = am*.7
         goto 11 

20 	 CLOSE(40) 

         Norb = Nt
        endif
        Iorb = Norb
             lspecies(1) = Iorb     
          endif  
             write(*,*)'              Halo Initialization is done' 
       Endif     

       write(*,'("  Np after halo generation= ",3i8)')iorb,
     &        lspecies(nspecies),lspecies(Nspecies-3*Irigid)  
        close(68)  

c
c                                   Satellite Initialization   
        write(*,*)'Satellite?'  

        If(Isat.eq.1.and.ISSat.eq.0)then    
           write(*,*)'1 Satellite'
           If(iMM.eq.1)then  
              write(*,*)'Major Merger: Msat = Mhalo'
              xtrunc =Xout  
              Rslarge =R_0  
              aMlarge = aMassh   
              Nwant = lspecies(nspecies)    
           else    

              xtrunc = xout_sat    
              Rvir_sat     = 442.97 * 
     *             (aMass_sat/1.e+11/(Om0*Overdens))**0.3333 ! kpc/h    
              Rs_sat2=  Rvir_sat/Con_sat   
              
              Rslarge = R_0/Rs_sat2 ! This conversion is because we don't know Rs_sat     
              Write(*,*)' Rs_sat= ',Rslarge,R_0     

              aMlarge = aMass_sat  
              Nwant = (aMass_sat/aMassd)*Id ! same mass per particle   
              
              R0_sat = Rslarge  ! (Rsnfw/ Rslarge)    
              V0_sat = sqrt((aMlarge/aMassh*fcon)*R0_sat)         


              Dsatix = 0.  
              Dsatiy= 0.
              Dsatiz = 0.   

              

              Vx0_sati = Vx0_sat/V_0    
              Vy0_sati = Vy0_sat/V_0  
              Vz0_sati = Vz0_sat/V_0    

              write(*,*)' Satellite Initialization Starts'        
 
              CALL Initialsat(Nwant,xtrunc)               
              
           endif     
           wspecies(1) = 1.0       
           store = coef(2)        
           coef(2) =0  
           Dsatix = D0x_sat/R_0  
           Dsatiy= D0y_sat/R_0  
           Dsatiz = D0z_sat/R_0  
           write(*,*)'Coords-sat ',Dsatix,D0x_sat,Dsatiz,D0z_sat    
           Vx0_sati = Vx0_sat/V_0    
           Vy0_sati = Vy0_sat/V_0  
           Vz0_sati = Vz0_sat/V_0   
           R0_sat = (Rsnfw/ Rslarge) 
           V0_sat = sqrt((aMlarge/aMassh*fcon)*R0_sat)    
           

           coef(2) =store     
           
        write(*,*)' Separation along X-axys: ',Dsatix*R_0  
        Call Shift  
          
      Endif  

c*******************************************************      
      If(ISSat.eq.1.and.Isat.eq.0)then       
         write(*,*)'Satellite System' 

         Ialmacen =  0          !counting satellites particles   


         Call CM( )      

         Do icsat =1,Icount   
            xtrunc =   10. 
            Con_sat = Csat_i (icsat)  
            Rvir_sat     = 442.97 *   
     *           (.7*aMsat_i(icsat)/1.e+11/(Om0*Overdens))**0.3333 ! kpc/h      
            
            Rs_sat2=  Rvir_sat/Con_sat       
              
            Rslarge = R_0/Rs_sat2 ! This conversion is because we don't know Rs_sat      
            Write(*,*)' Rs_sat= ',Rslarge,R_0     
c              pause
            aMlarge =.7*aMsat_i(icsat)   
            Nwant = (.7*aMsat_i(icsat)/aMassd)*Id ! same mass per particle      
              
            R0_sat = Rslarge    ! (Rsnfw/ Rslarge)   
            V0_sat = sqrt((aMlarge/aMassh*fcon)*R0_sat)         
            write(*,*)' Individual Satellite Initialization Starts'       

               
            wspecies(1) = 1.0    
            store = coef(2)   
            coef(2) =0  
            Dsatix = x_i(icsat)/R_0  
            Dsatiy= y_i(icsat)/R_0   
            Dsatiz = z_i(icsat)/R_0   
            Vx0_sati = vx_i(icsat)/V_0    
            Vy0_sati = vy_i(icsat)/V_0   
            Vz0_sati = vz_i(icsat)/V_0    
           

            CALL Initialsat(Nwant,xtrunc)              
           
            write(*,*)' Coords ',Dsatix,D0x_sat,Dsatiz,D0z_sat        
c           Vx0_sati = vx_i(icsat)/V_0    
c           Vy0_sati = vy_i(icsat)/V_0   
c           Vz0_sati = vz_i(icsat)/V_0   
c           R0_sat = (Rsnfw/ Rslarge) 
c           V0_sat = sqrt((aMlarge/aMassh*fcon)*R0_sat)   

           coef(2) =store   
           write(*,*)'Separationalong X-axys: ',Dsatix*R_0,icsat,
     *          lspecies(Nspecies),iorb,lspecies(1)       
      Enddo 
      
      Call Shift_cent( )       

      Endif   
      goto 100
c*******************************************************   
c 9    Write(*,*)'The model will be read from a file'     
      
c      Call Readmod( )    
      
c      WRITE(*,*)'Model Reading  is done'            

c      Call Disk_File( )           

c      pause  


      
c*******************************************************  

 100  Norb = lspecies(Nspecies)  
        write(*,*)Norb,iorb,"flag"   

   
c        write(*,*)'Xcmh= ', txcmh/twh,  
c     *       'Ycmh= ', tycmh/twh, 
c     *       'Zcmh= ', tzcmh/twh  
  
       write(*,*)'  N particles=',Norb 
       pMass(1) = aMassd/Id!/FCon*aMassh           !   /FCon*F(xout)  go back to physical units
      write (*,10) Om0,Overdens,aMassh,Rs,Con,
     &     xout,Norb,pMass(1)
      write (16,10) Om0,Overdens,aMassh,Rs,Con, 
     &     xout,Norb,pMass(1) 
 10   format('=== Omega_0=',f6.3,' Overdensity=',F6.1,
     .               '    Virial Mass of halo(Msun/h)=',g10.3,
     .        ' R_s(kpc/h)=',f7.2,' c=',f7.2,/
     .       '    Outer radius/R_s=',f7.2,' N_part=',i9,      
     .     ' Mass-of-1part(Msun/h)=',G12.4)     
 

      If(Nspecies.ge.2)then  
      
       Do i =2,Nspecies    
         pMass(i) = pMass(i-1)*2.
c         write (*,'(" particles=",i8," particle mass=",g11.4)')
c     &               Norb,pMass(i) 
c         write (16,'(" particles=",i8," particle mass=",g11.4)')
c     &               Norb,pMass(i) 
      Enddo 

        Do i =Nspecies-Irigid*3 ,Nspecies
          write (*,10) Om0,Overdens,aMassh,Rs,Con,
     &         xout,Norb,pMass(i) 
          write (16,10) Om0,Overdens,aMassh,Rs,Con,
     &         xout,Norb,pMass(i)
    
        Enddo  
      Endif 

c      write(*,*)'flag1'
c      Call  Analyze(Nvir)  
c      write(*,*)'flag2'

c      write(*,*)'If you want a Tipsy file type 1'
c      read(*,*)Itipsy
      itipsy=0
      if(itipsy.eq.0)then


      Call  Rescale(xscale,vscale,Isirtakis)   
      

      Nefft = lspecies(1)  
      write(*,*)' lspecies(1) = ',lspecies(1)  
      Do i=2,nspecies
         Nefft = Nefft + (lspecies(i)- lspecies(i-1))*wspecies(i)
         write(*,*)'   # of particles in specie-i = ', 
     *        (lspecies(i)- lspecies(i-1)) 
      Enddo

c      write(*,*)'Neff= ',Nefft

       If(Isirtakis.eq.0)then


       Do i=1,nspecies  ! Masses in Msun/h, Box in Mpc/h 
          wspecies(i) =FLOAT(NGRID)**3/Nefft*
     &  ((F(xout)/FCon*aMassh*coef(1)+aMassd*coef(2)
     &   +    aMassb*coef(3)+aMass_sat* coef(4)  )*3.64e-12/Om0/Box**3)*
     &         wspecies(i)     ! weights of particles
c          write(*,*)' weights= ',wspecies(i) 
                                                                                                     ! if disk mass is included!+aMassd 
       Enddo

        else 

       wspecies(1) = (FLOAT(NGRID)**3)*peso(1)*3.64e-12/Om0/Box**3       


c       wspecies(1) =  FLOAT(NGRID)**3/Nefft* 
c     &  ((F(xout)/FCon*aMassh*coef(1)/.7+aMassd*coef(2) 
c     &   +    aMassb*coef(3)+aMass_sat* coef(4)  )*3.64e-12/Om0/Box**3) 
 

        endif  ! ask if sirtakis is active   


       write(*,*)'  Weight Ratios: 2/1=',wspecies(2)/wspecies(1),
     *      ' 3/1=',wspecies(3)/wspecies(1)
       write(*,*)'  first weight      = ',wspecies(1)
       write(*,*)'  mass1p(Msun/h)    = ',wspecies(1)*
     &                                    (Box/FLOAT(NGRID))**3/
     *      3.64e-12*Om0 

      Vscale = Box*100.*hubble/NGRID  
      AexpV  = AEXPN - ASTEP/2. ! leap frog 
      CALL WriteData(Vscale,AexpV,Wtotal)
            EKIN = 0.5*SKINE/AEXPV**2       
           WRITE (*,'('' Ekin='',E12.4,'' Weght per cell='',g12.5, 
     .                   '' (   )'')')   EKIN,Wtotal/NGRID**3 
c           WRITE (16,'('' Ekin='',E12.4,'' Weght per cell='',g12.5, 
c     .                   '' (must be 1)'')')   EKIN,Wtotal/NGRID**3
C			        
      CALL WRTAPE   ! write header and control data
C                      write pt.dat file: time-step for particles
      WRITE (*,'('' Ntotal='',i9,'' must be equal='',i9)')
     &      Norb,lspecies(Nspecies)
      If(Nspecies.eq.0)Then
         Nparticles =lspecies(1)
      Else
         Nparticles =lspecies(Nspecies) 
      EndIf  
      Do ic1 =1,Nparticles 
         xpt(ic1) =astep 
      Enddo 
      open ( 60 , file = npath//'pt.dat' , form = 'unformatted' ) 
      write(60) (xpt(ic1),ic1=1,Nparticles)
c      write (*,*) ' pt=', (xpt(ic1),ic1=1,10)  


      else
         call  wtipsy( ) 
      endif
      
      STOP
      END
     
C---------------------------------------
      REAL*4 Function Conc(m)   ! Concentration.  m is in  10^11h-1M_sun
C---------------------------------------
       IMPLICIT REAL*4(A-H,O-Z)
      REAL*4 m
         Conc =10.**(1.17-0.084*LOG10(m))
      Return
      End
C--------------------------------------- 
      REAL*4 Function F(x)   ! x =r/r_s,  F= mass(r)/4pi r_s^3
C---------------------------------------
      IMPLICIT REAL*4(A-H,O-Z)
         y =max(x,1.e-10) 
         F =Max(LOG(1.+y) -y/(1.+y),1.e-16)
      Return
      End  
C---------------------------------------
      REAL*4 Function Mmx(x)  ! Self contribution of halo to Vrms 
C---------------------------------------
      IMPLICIT REAL*4(A-H,O-Z)
c      REAL*4  x
         Mmx =(LOG(1.+x) -x/(1.+x))/x**3/(1.+x)**2
      Return
      End 
C---------------------------------------
      REAL*4 Function Mmhb(x)   ! Contribution of halo to Bulge Vrms 
C---------------------------------------
      INCLUDE 'PMparameters.h'   
      INCLUDE 'Rodin2.h'    


      y = x*Rs_b/Rsnfw 
      Mmhb = (LOG(1.+y) -y/(1.+y))/(x*(1.+x))**3     

      Return 
      End 
C---------------------------------------  
      SUBROUTINE Analyze(Nvir) !      
C---------------------------------------  
C
C
      INCLUDE 'PMparameters.h'   
      INCLUDE 'Rodin2.h'   
      PARAMETER( Lrad=30 )  
      REAL*4           Radd(Lrad),Vr(Lrad),Sr(Lrad),
     &                 Rdi(Lrad),V3rms(Lrad)
      Dimension        Nr(Lrad),wNr(Lrad),wwNr(Lrad)
      REAL*4           Mmx,AINTG
      External         Mmx,AINTG
      OPEN(66,FILE=npath//'hprof.dat',status='unknown', access='append') 
      
       

      vscale = 208.e-5*sqrt(M_0/Rs) ! velocity scale in km/s, 1e-5*sqrt(G*2.d33/3e21) =208.e-5
C                             Set up radii
      hr  =0.05    !0.5*Con/Lrad 
      Cons = xout/10.**((Lrad-1)*hr) 
      DO i=1,Lrad            
         Radd(i) = Cons*10.**((i-1)*hr)  
c         write(*,*)'nspecies ',nspecies
         Neff =0.
 555     Rdi(i)   = 0.
         Nr(i)   = 0 ! Cumulative number of particles
         wNr(i) = 0. ! Cumulative effective number of particles
       wwNr(i) = 0.! Cumulative mass weighted effective number of particles
         Vr(i)   = 0.
         Sr(i)   = 0.
         V3rms(i)= 0.
      ENDDO
      N200   =0 
      Nvir     =0
      istep =0 
      DO j=Id+1,Norb ! loop over all the halo particles  
            x =sqrt(xc(j)**2+yc(j)**2+zc(j)**2)
c            write(*,*)x ,j
            If(x.lt.CON)Then 
               Nvir =Nvir +1
            endif 
            Vradial =(vxc(j)*xc(j)+vyc(j)*yc(j)+vzc(j)*zc(j))/
     &                  max(x,1.e-5) 
            ir =1                                  ! find the bin
 10         If(x.lt.Radd(ir))Then      ! inside the bin
              Do jb =1,nspecies
                 If(x.le.Xb(jb))then
                     tw = wspecies(jb)
c                     write(*,*)x,Xb(jb),tw
                     goto 11 
                  Endif  
               Enddo
 11            wNr(ir) =wNr(ir) +tw !cumulative effective number of particles
               wwNr(ir) =wwNr(ir) +tw*aMassd/Id !cumulative mass
               Neff = Neff + tw
               Nr(ir) =Nr(ir) +1
               Rdi(ir) =Rdi(ir) +x
                 Vr(ir) =Vr(ir)  +Vradial
                 Sr(ir) =Sr(ir)  +Vradial**2
                 V3rms(ir)=V3rms(ir)+vxc(j)**2+vyc(j)**2+vzc(j)**2
c                 write(*,*) Nr(ir),Rd(ir),ir,x,Radd(ir) 
            Else  
c               Rdi(ir) =Rdi(ir) +x
                 ir =ir+1
                 if(ir.le.Lrad)goto 10
            Endif
      ENDDO

      V0 = 0. 
      Mr = 0
      wMr = 0.
c      dNorm = F(xout)/Norb
      dNorm = F(xout)/Neff
      write (*,200)
      write (66,200)
 200  format(' Bin Radius Particle',1x,
     &     ' N_part V_circ',4x,'Density',
     &     5x,'Expected  RadVelocity',
     &     ' RMS_V_r Expected  RMS_V-3D',/
     &     '  (kpc/h)   Rad(kpc/h)',8x,'(km/s)',
     &     30X,'(km/s)  (km/s)') 
      Do i=1,Lrad  
         Volume   = 4.188*Radd(i)**3!4Pi/3 = 4.188
         dVolume  = Volume -V0 ! V0 is the volume up to the last bin.

         Mr       = Mr + Nr(i) 
         wMr =wMr  + wwNr(i)

c 199     Vcirc    = 208.e-5*   
c     &     sqrt(max(Mr/(Rs*Radd(i)),1.e-12))
 199     Vcirc    = 208.e-5*   
     &        sqrt(max(wMr/(Rs*Radd(i)),1.e-12))

c     write (*,*) ' M=',Mr,' vscale=',vscale,Radd(i),aMass1 
c         dN       = Nr(i)/dVolume*dNorm
         dN =  wNr(i)/dVolume*dNorm !wNr(i) cumulative weighted numerber of particles
         Velrad   = Vr(i)/max(Nr(i),1)
         Raver    = Rdi(i)/max(Nr(i),1)
         If(Nr(i).eq.0)Raver =Radd(i) 
         x        = Raver
         xmax     =  max(x,5.e+2)
  
         Svel     =sqrt(max(Sr(i)/max(Nr(i),1)-Velrad**2,1.e-10))
         V3       =sqrt(max(V3rms(i),1.e-10)/max(Nr(i),1)) 
           Rho_ex   =1./(4.*3.14159*Raver*(1.+Raver)**2)
         Svel_ex  =x*(1.+x)**2*AINTG(Mmx,x,xmax)

         write (*,100) Radd(i)*Rs,Raver*Rs,        
     &                    Mr,Vcirc,dN,Rho_ex, 
     &                    Velrad*vscale,Svel*vscale,
     &                    sqrt(Svel_ex)*vscale,V3*vscale 
         write (66,100) Radd(i)*Rs,Raver*Rs,
     &        Mr,Vcirc,dN,Rho_ex,
     &        Velrad*vscale,Svel*vscale,
     &        sqrt(Svel_ex)*vscale,V3*vscale
         
 100     format(2g11.3,i7,3g11.4,g12.2,g11.3,3g11.3)
         V0=Volume
c         CLOSE(66)   

      EndDo 
      Return
      End 
   
C---------------------------------------   
      SUBROUTINE Initialnew(pMdi)  !   !  Set initial Conditions         
C---------------------------------------   
C                  xout = outer radius in units of Rs
C                  Xb - > boundaries btween species
C                  pMdi = mass p/particle in halo units
      INCLUDE 'PMparameters.h'
      INCLUDE 'Rodin2.h'

      PARAMETER(NradM=15000)  
      REAL*4    Radd(NradM),AINTG,Mmx 
c      Real*8 hr
      Dimension Nr(NradM)       !,Xb(Nspecies-3)  
      External  Mmx 
      hr  =xout/Nrad     
      write(*,*)'hr= ',hr,'xout= ',xout      
      k=0 
      write(*,*)'Nrad',Nrad 
      wspecies(1)=1.  !relative weights of particles in the 1st specie  
 200  format('   Bin  radius  Npart  Ntheo %error')! Specie N(r<))   sigr')  
      write (68,200)   
      write (*,200) 
      
      Do i=1,Nspecies-3*Irigid           !Loop over the species  
         Nc =0 
         aMlimit= 10.*pMdi*2**(i-1)              !10.*pMdi*2**(i-1)  !Ten particles of the corresponding specie    
c         write(*,*)'aMlimit= ',aMlimit 
         If(i.eq.1)then 
            correct = Xb(i)/hr 
            write(*,*)hr,i,Xb(i)  
            iYb = Int(correct) ! 
            correct = correct/float(iYb)
            iYbp=0
            aMi = F(Xb(i))     !mass inside the first big bin 
            write(*,*)'Xb(i), iYb,aMi ',Xb(i), iYb,aMi 
         Else
          write(*,*)'wspecies i, i-1',wspecies(i),wspecies(i-1) 
           wspecies(i)= wspecies(i-1)*2.d0    
          write(*,*)'Xb(i),hr',Xb(i),hr
           correct = Xb(i)/hr
          write(*,*)'correct',correct
           iYb = Int(correct)    ! Index for the outer boundary for the ith-specie  
          write(*,*)'iYb',iYb
           correct = correct/float(iYb)
          write(*,*)'correct',correct
           iYbp = Int(Xb(i-1)/hr) ! Index for the inner boundary for the ith-specie
          write(*,*)'iYbp, Xb(i-1), hr',iYbp,Xb(i-1),hr
           aMi = F(Xb(i)) - F(Xb(i-1)) ! mass inside the big bin
          write(*,*)'F(Xb(i)),F(Xb(i-1))',F(Xb(i)),F(Xb(i-1))   
           write(*,*)'Xb(i), iYb,aMi ',Xb(i), iYb,aMi 

         Endif    
            
         aMleft = aMi           !mass available for this specie   
         rmin =    Xb(i-1) 
          rbig = rmin     
          Nbig =0
          pMbig =0.
 
 887      k=k+1                 !scan over the small bins   
   
               
          If(k.eq.1)then  
             dMbig = F(hr*correct) !mass inside the 1st bin  
          Else
            dMbig  =   F(float(k)*hr*correct) - F(float(k-1)*hr*correct) !mass inside this small bin
         Endif
  
         pMbig = pMbig + dMbig  ! stores total mass inside  the  present bin 
         rbig = rbig+hr*correct    
         
c         if(k.gt.iYb.or.(rbig-Xb(i)).gt.0.)then  
         if(k.le.iYb)then      
c            write(*,*)'rbig > Xb(i)',k,iYb,i,rbig,Xb(i)
c            write(*,*)'****',k,iYb,i,rbig,Xb(i),iYb*hr*correct     
c            stop
         endif     
         

         rwrite = (rbig + rmin)/2.d0   
                       

         aMleft =  F(xb(i)) -  F(rbig) !aMleft -  pMbig ! leftover mass for that specie  
c         if(aMleft.lt.0.)then  
c            write(*,*)'Negative aMleft'
c         endif
         pass=abs(aMleft-F(xb(i))+F(rbig))/
     *        max((F(xb(i))-F(rbig)),1.e0)    
         if(pass.gt.5.e-2)then                
            write(*,*)'pMbig= ',pMbig,aMleft,i,k,  
     *     abs(aMleft - F(xb(i))+ F(rbig))/max((F(xb(i))-F(rbig)),1.e0),      
     *           dMbig,xb(i),rbig,rmin,  
     *           F(xb(i)),F(rbig)        
c     Stop       
         endif     
         If(aMleft.lt.aMlimit)then !The scan over this specie is over         
            If(k.eq. iYb)then  ! Both boundaries are equal  
               Nbig  =    Int( pMbig /(pMdi*2**(i-1))+.5) !Number in that bin   
               Call writec(Nbig,i,rbig,rmin,pMdi)   
               Nc = Nc + Nbig   
            Else   
c               Nbig = Int( (pMbig +aMleft)/(pMdi*2**(i-1))+.5) ! Include all the lefting particles 
   
               rbig = xb(i)    
               Nbig = Int( (F(rbig)-F(rmin))/(pMdi*2**(i-1))+.5) ! Include all the lefting particles   
    
               write(*,*)'Flag',rbig, xb(i),iYb*hr*correct,iYb,k,aMleft,   
     *              aMlimit,F(xb(i))-F(rbig),i      
                     
               
               Call writec (Nbig,i,rbig,rmin,pMdi) ! and go to the next specie.  
               Nc = Nc + Nbig   
               k = iYb   
            Endif   
      
         Else ! We are not done yet with the scan of that specie   
            If(pMbig.lt.aMlimit)goto 887 ! join this small bin to the following
            Nbig  =    Int( pMbig /(pMdi*2**(i-1))+.5)   
            Call writec(Nbig,i,rbig,rmin,pMdi)   
            Nc = Nc + Nbig  
            rmin =rbig   
            pMbig =0. ! it will star with a new small bin   
            GOTO  887  
         Endif      
    
         If(i.eq.1)then 
            lspecies(1)=Nc +  lspecies(1)  
         else
            lspecies(i)= Nc + lspecies(i-1) !total number of particles up to that specie
         endif
      Enddo                              ! end i  -- loop of species
      Return  
      End        
C--------------------------------------- 
      SUBROUTINE writec(Nbig,ik,rbig,rmin,Pmdi) 
                                !  Calculate Vhalo and Write in the arrays of coordinates and velocities 
C--------------------------------------- 
      INCLUDE 'PMparameters.h'  
      INCLUDE 'Rodin2.h'  
      REAL*4 AINTG 
      External Mmd,Mmx,F,AINTG,d, aMmxh         
      Open(70,FILE=npath//'sigmah.dat',status='unknown',access='append') 
c      write(*,*)'txcmt= ', txcmh,tycmh,tzcmh   
  
      ibin = ibin +1            !counting the small bins   
c      write(*,*)'rmin= ',rmin,'rbig= ',rbig,'Nbig= ',Nbig  
 
      Ncount = Ncount + Nbig  !counting the particles  
      x = rbig                     !(rbig + rmin)/2.d0  
      xmax   = max(x,Xout)  
      at = AINTG(Mmx,x,xmax)    
      bt = AINTG(Mmd,x,xmax)*Fcon/RMvirMd*coef(2) !disk contribution  
      ct = AINTG(aMmxh,x,xmax)*Fcon*aMassb/aMassh*coef(3) ! bulge contribution   

      sig_r2 = x*(1.+x)**2*(at + bt+ ct)              
c      write(*,*) ' r=',x*Rs,x,'  bt/at= ',bt/at
     
c      if(x.le.(1./Rs))then  
c        if(x.le.(0.1))then
c         sig_r  = sqrt(sig_r2)*(.5*Rs*x+ .5) ! rms of radial velocity !  (.8*Rs*x+ .2) (.6*Rs*x+.4)  
c     write(*,*)"flag",sig_r ,sqrt(sig_r2),x
c     Stop 
c      else      
         sig_r  = sqrt(sig_r2)  ! rms of radial velocity      
c      endif    
      write(70,71)x*Rs,sig_r*208.e-5*sqrt(aMassh/Fcon/Rs),          
     *     sqrt(x*(1.+x)**2*at)*208.e-5*sqrt(aMassh/Fcon/Rs),    
     *     sqrt(x*(1.+x)**2*bt)*208.e-5*sqrt(aMassh/Fcon/Rs)             

c      write(70,71)x*Rs,sig_r,        
c     *     sqrt(x*(1.+x)**2*at),     
c     *     sqrt(x*(1.+x)**2*bt)    
        
 222  Vescap2= 2.*(LOG(1.+x)/x +    
     *     (1. - Exp(-x*Rs/Rd) )/x*Fcon/RMvirMd*coef(2) + !disk contribution  
     *     1./(x*Rs_b/Rs+1.)*Fcon*aMassb/aMassh*coef(3)) ! bulge contribution          
                                ! square of escape velocity                
c     write(*,*)'pots ratio= ',(LOG(1.+x)/x)/dpot(x*Rs/Rd)      

      x=(rbig + rmin)/2.d0   
         hr  =xout/Nrad    
            theo = (F(rbig) - F(rmin))/wspecies(ik)/Pmdi 
            write(68,250)ibin,x,Nbig,theo,(Nbig-theo)/theo,   
     *           ik,Ncount,sig_r,rbig,rmin 
            If(abs((Nbig-theo)/theo).gt..05)then 
               write (*,200) 
               write(*,250)ibin,x,Nbig,theo,(Nbig-theo)/theo, 
     *              ik,Ncount,sig_r,rbig,rmin 
               write(*,*)'Stopped'    
                
 200           format('   Bin  radius  Npart Ntheo %error Specie N(r<)',
     *              'sigr  rbig rmin')   
             write(*,*)'Check Mass limit perhaps is too big(in writec)',  
     *              'Increase Nrad(Line 29)'   
            stop  
            Endif  
            
250       format(i5,f10.3,i8,f10.3,f10.3,i3,i8,3f10.3)  
c250       format(i7,1x,f10.3,1x,i7,f10.3,1x,f10.3)    
c          write(*,*)'iorb =',iorb
          sumv2 = 0.
         Do j=1,Nbig
            iorb =iorb +1       !starts where the disk stoped
            If(iorb.gt.Nmo)Then  
               write (*,*) ' Too many particles:',iorb,Nmo
               STOP  
            EndIf
 70         vxc(iorb) =sig_r*gauss(Nseed)
            vyc(iorb) =sig_r*gauss(Nseed)
            vzc(iorb) =sig_r*gauss(Nseed)
c            if(abs(vzc(iorb)).le.0.)then
c               write(*,*)'HERE',vxc(iorb),iorb
c               Stop  
c            Endif 
  
c            if( (Vescap2/(3.*sig_r**2)).gt.1.)then 
c               WRITE(*,*)'Vescape > sig ',x,Vescap2,(3.*sig_r**2)
c            endif   

            Vv       =vxc(iorb)**2 +
     &           vyc(iorb)**2 + 
     &           vzc(iorb)**2  
           
           

            If(Vv.ge.Vescap2)Then       
c     write(*,*)'unbound',Vv,Vescap2,x,Nbig,theo,ik,Ncount 
               goto 70  
            EndIf   
            sumv2 = sumv2 +  Vv       


            cc    = 2.*Randd(Nseed)-1.  
            costh = sqrt(abs(1.-cc**2))
            phi   = 2.*pi*Randd(Nseed)
            rsin  = x*cc        !           SIN(theta)
            rcos  = x*costh     ! COS(theta) 
            xc(iorb) =rcos*COS(phi) 
            yc(iorb) =rcos*SIN(phi) 
            zc(iorb) =rsin
           
            if(iorb.gt.id)then
               txcmh = txcmh + xc(iorb)* wspecies(ik)    
               tycmh = tycmh + yc(iorb)* wspecies(ik)        
               tzcmh = tzcmh + zc(iorb)* wspecies(ik)     
               
               twh = twh + wspecies(ik) 
               peso(iorb) = wspecies(ik) 
            endif
                   
         Enddo 
         close(70)
c         write(*,*)'twh',twh
c         pause 
c         write(*,*)sqrt(sumv2/(3.*Nbig))*208.e-5*sqrt(aMassh/Fcon/Rs),   
c     *        x*Rs  
c         write(*,*)'Halo cm= ',txcmh/twh,tycmh/twh,tzcmh/twh 
 71      format (4(1x,f10.4))   
         Return 
         End
C---------------------------------------
      SUBROUTINE Rescale(xscale,vscale ,Isirtakis)  !  Rescale NFW -> ART
C---------------------------------------
C                       
      INCLUDE 'PMparameters.h'      
      INCLUDE 'Rodin2.h'
       Real*8  scXcm,scYcm,scZcm,scXcmd,scYcmd,scZcmd,scWeight,
     *     xmin,ymin,zmin,xmax,ymax,zmax,Sv2

       svx=0.  
       svy=0.  
       svz=0.  

       Sv2 =0.
       scXcm =0.
       scYcm =0.
       scZcm =0.
       scXcmd =0.
       scYcmd =0. 
       scZcmd =0.
       scWeight =0.
      scWeightd =0
      xr =NGRID+1. 


       vscale2 = V_0/(Box*100/AEXPN/NGRID)  


      If(Iconf.eq.0)Ndr=0

      xmin =NGRID*2 
      ymin =NGRID*2
      zmin =NGRID*2  
      xmax =-xmin
      ymax =-ymin
      zmax =-zmin
      Vhubble = AEXPN*sqrt(Om0/AEXPN+Oml0*AEXPN**2)
      Vhh     = Vhubble*xscale  !aexp is cancelled with the 1/aexp from xscale
c      Vxcenter =Vcenter*AEXPN/(Box*100./NGRID) 
      

c      write(*,*)'Hubble',sqrt(Om0/AEXPN+Oml0*AEXPN**2),Oml0*AEXPN**2,
c    *     Om0/AEXPN,Vhh, xscale 
 
         Do iorb=1,Norb  
                          ! subtract hubble velocity
     

   
            If(iorb.ge.1000007)then       
               svx =  svx + vxc(iorb) 
               svy =  svy + vyc(iorb)    
               svz =  svz + vzc(iorb)   
c               write(*,*)'test-2',V_0*vxc(iorb),iorb  
c               stop  
            Endif      


            vpaso1 = vxc(iorb)  
            vpaso2 = vyc(iorb)     
            vpaso3 = vzc(iorb)   
            
c            write(*,*)'Vtest= ',vxc(iorb),vscale2*vxc(iorb),
c     &	    Vhh*(xc(iorb))
c        	pause        

	    If(Isirtakis.eq.0)then !Rodin


            vxc(iorb) =vscale*vxc(iorb)-0*Vhh*(xc(iorb))              
            vyc(iorb) =vscale*vyc(iorb)-0*Vhh*(yc(iorb))                
            vzc(iorb) =vscale*vzc(iorb)-0*Vhh*(zc(iorb))                   

            
	else 

 	  vxc(iorb) =vscale2*vxc(iorb)-0*Vhh*(xc(iorb))              
          vyc(iorb) =vscale2*vyc(iorb)-0*Vhh*(yc(iorb))                
          vzc(iorb) =vscale2*vzc(iorb)-0*Vhh*(zc(iorb))                   


        endif


c            if(iorb.le.Nd)then   
c               write(*,*)'Scaled difference V  disk',
c     *              vxc(iorb), vyc(iorb),vzc(iorb) ,iorb  
c            endif
    
            
c     *           (VXc(iorb)+Vhh*xc(iorb))/vscale-vpaso1,
c     *           (VYc(iorb)+Vhh*yc(iorb))/vscale-vpaso2,
c     *           (Vzc(iorb)+Vhh*zc(iorb))/vscale -vpaso3

c

            xc(iorb) = xscale*xc(iorb)+NGRID/2.
            yc(iorb) = xscale*yc(iorb)+NGRID/2.                        
            zc(iorb) = xscale*zc(iorb)+NGRID/2.

            If(xc(iorb).lt.1.)xc(iorb)=xc(iorb)+NGRID 
            If(xc(iorb).ge.xr)xc(iorb)=xc(iorb)-NGRID 
            If(yc(iorb).lt.1.)yc(iorb)=yc(iorb)+NGRID 
            If(yc(iorb).ge.xr)yc(iorb)=yc(iorb)-NGRID 
            If(zc(iorb).lt.1.)zc(iorb)=zc(iorb)+NGRID 
            If(zc(iorb).ge.xr)zc(iorb)=zc(iorb)-NGRID 
            xmin =min(xc(iorb),xmin)  
            ymin =min(yc(iorb),ymin)
            zmin =min(zc(iorb),zmin)
            xmax =max(xc(iorb),xmax)
            ymax =max(yc(iorb),ymax)    
            zmax =max(zc(iorb),zmax)

           

c            if(iorb.gt.Ndr)then
               scXcm =scXcm + xc(iorb)*peso(iorb) 
               scYcm =scYcm + yc(iorb)*peso(iorb)
               scZcm =scZcm + zc(iorb)*peso(iorb) 
               scWeight = scWeight + peso(iorb) 


c            else
c               scXcmd =scXcmd + xc(iorb)
c               scYcmd =scYcmd + yc(iorb)
c               scZcmd =scZcmd + zc(iorb)
             
c            endif 

            

               
            

  
         EndDo
         write (*,*) ' -- Done rescaling of coordinates --'


	write(*,*)'MM test',iorb,peso(1499999),peso(2999996) 

         If(Iconf.gt.2)then
            write(*,*)'Difference Xcmh -Xcmd'
            write(*,*)scXcm/scWeight-scXcmd/Ndr,
     *        scYcm/scWeight-scYcmd/Ndr,
     *           scZcm/scWeight-scZcmd/Ndr 
         else  If(Iconf.eq.0)then
            write(*,*)'Halo CM (scaled to the Box)'
            write(*,*)scXcm/scWeight,
     *        scYcm/scWeight, 
     *           scZcm/scWeight
	    write(*,*)'Total Mass',scWeight 
            else  If(Iconf.eq.1)then
               write(*,*)'Disk Xcmd'
               write(*,*)scXcmd/Ndr,
     *        scYcmd/Ndr,
     *         scZcmd/Ndr
         endif

c         write (*,*) '    Min=',xmin,ymin,zmin
c         write (*,*) '    Max=',xmax,ymax,zmax
         If(min(xmin,ymin,zmin).lt.1.0)Then
            write (*,*) ' wrong minimum of coordinates'
            stop
         EndIf
         If(max(xmax,ymax,zmax).lt.1.0)Then 
            write (*,*) ' wrong maximum of coordinates'
            stop
         EndIf 
         
         write(*,*)'This is written in the binary file'  
         write(*,*)'vcm-sat', svx/Nsat*V_0,
     &        svy/Nsat*V_0,svz/Nsat*V_0            

      Return 
      End
C----------------------------------------------------------------
      SUBROUTINE InitValues(NBYTE,SCALEL) 
C----------------------------------------------------------------
      INCLUDE 'PMparameters.h'
      INCLUDE 'Rodin2.h'  
      Character            Answer*1
      Real                 INPUT

C*************************************


      Call Readconf 
 
      If(Iconf.eq.1.or.Iconf.eq.2)then 
         Con = 15.
         aMassh = 2.e12
         Xout = Con
      Endif   

c      write(*,*)'Xb',Xb(1),Xb(2),Xb(3),
c     *     Xb(4),Xb(5),Nspecies-3   

      RMvirMd     = aMassh/aMassd 
      pMd             =      aMassd/Nd 
      Rref            =       Rrefkp/Rd  
      zd               =        Zdkp/Rd  
      aMd              =            1.e0  
      RT               =        Rtkp/Rd 
      zcoef          =         .995e0 
       a                  =     0.e0 
       delta            =    0.1e-2       
C*************************************




c       HEADER='N=128x256 L=20h-1CDMsig=0.239 z=30-----------'
c       write (*,*) '------ Enter Header for the run up to 45 characters'
c       write (*,*) '       Example:'
c       write (*,'(A)') HEADER                                   
c       read  (*,'(A)') HEADER
       write (*,'(A)') HEADER
c       AEXPN =INPUT(' Initial expansion parameter (0-1)=')
c       ASTEP =INPUT(' Step of the expansion parameter  =')
       AMPLT =0.
c       Ocurv =INPUT(' Omega_curvature   at present     =')
c       Nseed =INPUT(' Random seed (Integer  1-2^31)    =') 
 
c       hubble=INPUT(' The Hubble constant (h=H/100)    =')
c       Om0   =INPUT(' Omega_0 matter    at present     =')
c       Oml0  =INPUT(' Omega_lambda      at present     =')

       NBYTE = NPAGE*6*4
c       Nspecies = 1
c       Nspecies =INPUT(' Number of mass species (0-10) = ') 
       
             W     = (FLOAT(NGRID)/FLOAT(NROW))**3
            PARTW = W 
      If(Nspecies.eq.1)Then  ! constant mass resolution
         write(*,*) ' You will use ART code with one particle specie '
       Else    ! multiple mass resolution
              write(*,*) ' You will use multiple mass resolution'
              write(*,*) ' Number of mass species is equal ',Nspecies
c              Do i =1,nspecies
c                Xb(i)        = INPUT('Change radius in units of Rs => ') !At this radius 
c              Enddo             ! we change mass/p/partcle
       endif  
       ISTEP = 0
       TINTG = 0.
       AU0   = 0.
       AEU0  = 0. 
       EKIN  = 0.
       EKIN1 = 0.
       EKIN2 = 0.
       NROWC = NROW
       NGRIDC= NGRID
      RETURN
      END
C---------------------------------------------
C                       Write current data to  disk/tape
C                        for PMstartM - multiple masses
      SUBROUTINE WriteData(Vscale,AexpV,Wtotal) 
C----------------------------------------------
c     INCLUDE 'initial.h'
      INCLUDE 'PMparameters.h'
      INCLUDE 'Rodin2.h'      
       COMMON / KINENRG/ SKINE,SX,SY,SZ,SX2,SY2,SZ2
C
       Real*8 wcmx,wcmy,wcmz,wwr,wcmxd,wcmyd,wcmzd
c       Open(36,file='conf.dat',status='unknown') 
       

       wcmx =0.d0
       wcmy=0.d0  
       wcmz=0.d0
       wcmxd =0.d0
       wcmyd=0.d0  
       wcmzd=0.d0
       wwr=0.d0 
       sx2 = 0.d0
       sy2 = 0.d0
       sz2 = 0.d0
       sx = 0.d0
       sy = 0.d0
       sz= 0.d0
       

      XMAX  = FLOAT(NGRID) + 1.
      XSHF  = FLOAT(NGRID)
      SKINE =0.
      Wtotal =0.
       jstart =1
       
       Do j =1,Nspecies
          vvx =0.
          vvy =0.
          vvz =0.
          vx2 =0.
          vy2 =0.  
          vz2 =0.
          jend  =lspecies(j)   
          W     =wspecies(j)
c          write(*,*)'species= ',j,jend,jstart,w     

          If(jend.gt.(jstart-Irigid))Then   
          Do i=jstart,jend             !Total Kinetic Energy, and total weight
             vvx = vvx + vxc(i)     ! for that specie are calculated in this loop.
             vvy = vvy + vyc(i)
             vvz = vvz + vzc(i)
             vx2 = vx2 + vxc(i)**2
             vy2 = vy2 + vyc(i)**2
             vz2 = vz2 + vzc(i)**2
              SKINE    = SKINE + W*(vxc(i)**2 +vyc(i)**2 +vzc(i)**2) 
              
              Wtotal =Wtotal +w 
          EndDo
           npp =max(jend-jstart+1,1)
           v2 = sqrt((vx2 +vy2+vz2)/npp)/AexpV*Vscale !Scaling to ART units. 
          
           vvx =vvx/npp/AexpV*Vscale
           vvy =vvy/npp/AexpV*Vscale
           vvz =vvz/npp/AexpV*Vscale


           write (*,350) vvx,vvy,vvz,v2,jend-jstart+1,j,W 
 350       format('   V(km/s)=',4g11.3, ' Npart=',i8,' Species=',i3,
     &            ' weight =',g11.3)  
          jstart =  jend+1
          Endif 
       EndDo
      Ibuff =0
      KROW =0
      Do i=1,lspecies(Nspecies) !Over all the particles
C                                Periodical boundary conditions
c	       IF(XPt(i).GT.XMAX)	      XPt(i)=XPt(i)-XSHF
c	       IF(XPt(i).LE.1.)	            XPt(i)=XPt(i)+XSHF
c	       IF(YPt(i).GT.XMAX)	      YPt(i)=YPt(i)-XSHF 
c	       IF(YPt(i).LE.1.)	            YPt(i)=YPt(i)+XSHF
c	       IF(ZPt(i).GT.XMAX)	      ZPt(i)=ZPt(i)-XSHF
c	       IF(ZPt(i).LE.1.)	            ZPt(i)=ZPt(i)+XSHF
             Ibuff              = Ibuff +1
             XPAR(Ibuff) = xc(i)  
             YPAR(Ibuff) = yc(i)
             ZPAR(Ibuff) = zc(i)
             VX(Ibuff)   = vxc(i)
             VY(Ibuff)   = vyc(i)
             VZ(Ibuff)   = vzc(i)
             

c             write(36,38)i,xc(i)-63.01,yc(i)-63.01,zc(i)-63.01,  
c     *            vxc(i)/AexpV*Vscale,
c     *            vyc(i)/AexpV*Vscale,vzc(i)/AexpV*Vscale,
c     *            peso(i)    
c             write(*,*)'xyz',xc(i),yc(i),zc(i)
c     Center of mass calculation 
             if(i.gt.Ndr)then 
                wcmx = wcmx  + xc(i)*peso(i)
                wcmy = wcmy  + yc(i)*peso(i)
                wcmz = wcmz  + zc(i)*peso(i)  
                wwr = wwr  + peso(i) 
                sx = sx + vxc(i) 
                sy = sy + vyc(i)
                sz = sz + vzc(i)
                sx2 = sx2 + vxc(i)**2
                sy2 = sy2 + vyc(i)**2
                sz2 = sz2 + vzc(i)**2

     
c                write(36,38)i,xc(i)-63.01,yc(i)-63.01,zc(i)-63.01,  
c     *               vxc(i)/AexpV*Vscale,
c     *               vyc(i)/AexpV*Vscale,vzc(i)/AexpV*Vscale,
c     *               peso(i)  

             else
                wcmxd = wcmxd  + xc(i)
                wcmyd = wcmyd  + yc(i)
                wcmzd = wcmzd  + zc(i)
              
             endif 
             
             If(Ibuff.ge.NPAGE)Then
                KROW = KROW +1
c                write (*,*) ' Write page=',KROW,' i=',i,Ibuff 
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
c          Close(36)   
          write(*,*)'CM of written particles'
          If(Iconf.eq.0)then
             write(*,*)'Halo',wcmx/wwr,wcmy/wwr,wcmz/wwr   
          else    If(Iconf.eq.1)then
             write(*,*)'Disk',wcmxd/Ndr,wcmyd/Ndr,wcmzd/Ndr 
          else    If(Iconf.ge.3)then 
             write(*,*)'Halo',wcmx/wwr,wcmy/wwr,wcmz/wwr
             write(*,*)'Disk',wcmxd/Ndr,wcmyd/Ndr,wcmzd/Ndr 
          endif
      RETURN
      END  
C---------------------------------- Read in variables      
      REAL FUNCTION INPUT(text)
C--------------------------------------------------
      Character text*(*)
          write (*,'(A,$)')text
          read (*,*) x
          INPUT =x
      Return
      End
        
C----------------------------------- Simpson integration
      REAL*4 FUNCTION AINTG(FUNC,A,B)
C--------------------------------------------------
      IMPLICIT REAL*4 (A-H,O-Z)
      PARAMETER (EPS=2.0d-5, JMAX=22) 
      EXTERNAL FUNC
c      REAL*4 AINTG 
      OST=-1.d30 
      OS  = -1.d30
      ST   =0. 
      DO 11 J=1,JMAX
        CALL TRAPZD(FUNC,A,B,ST,J)
        AINTG=(4.0d0*ST-OST)/3.0d0
        IF (ABS(AINTG-OS).Le.EPS*ABS(OS).and.J.GT.6) RETURN
        OS=AINTG
        OST=ST
11    CONTINUE
      WRITE (16,*)'Integration did not converge'
      RETURN 
      END
C----------------------------------------------
      SUBROUTINE TRAPZD(FUNCC,A,B,S,N)
C--------------------------------------------------
      IMPLICIT REAL*4 (A-H,O-Z)
        SAVE IT
        EXTERNAL FUNCC
      IF (N.EQ.1) THEN
        S=0.5d0*(B-A)*(FUNCC(A)+FUNCC(B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5D0*DEL
        SUM=0.0D0
        DO 11 J=1,IT
          SUM=SUM+FUNCC(X)
          X=X+DEL
11      CONTINUE
        S=0.5D0*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END
  
C--------------------------------------------------
      REAL*4 function GAUSS(M)
C--------------------------------------------------
      IMPLICIT REAL*4 (A-H,O-Z)
      X=0.
      DO 1 I=1,5
	  X=X+RANDd(M)
1      CONTINUE
      X2=1.5491933*(X-2.5)
      GAUSS=X2*(1.-0.01*(3.-X2*X2))
      RETURN
      END
C--------------------------------------------------
      REAL*4 function RANDd(M)
C--------------------------------------------------
      IMPLICIT REAL*4 (A-H,O-Z)
       DATA LC,AM,KI,K1,K2,K3,K4,L1,L2,L3,L4
     +	/453815927,2147483648.,2147483647,536870912,131072,256,
     +	 16777216,4,16384,8388608,128/
      ML=M/K1*K1
      M1=(M-ML)*L1
      ML=M/K2*K2
      M2=(M-ML)*L2
      ML=M/K3*K3
      M3=(M-ML)*L3
      ML=M/K4*K4
      M4=(M-ML)*L4
      M5=KI-M
      IF(M1.GE.M5)M1=M1-KI-1
      ML=M+M1
      M5=KI-ML
      IF(M2.GE.M5)M2=M2-KI-1
      ML=ML+M2
      M5=KI-ML
      IF(M3.GE.M5)M3=M3-KI-1
      ML=ML+M3
      M5=KI-ML
      IF(M4.GE.M5)M4=M4-KI-1
      ML=ML+M4
      M5=KI-ML
      IF(LC.GE.M5)ML=ML-KI-1
      M=ML+LC
      RANDd=M/AM
      RETURN
      END
C---------------------------------------------------
C                                  Read  current data from disk/tape,
C                                  Open files
C                                  Nrecl is the number of values in a record
C                                  Npage is the number of particles in a page
      SUBROUTINE RDTAPE
C---------------------------------------------------
      INCLUDE 'PMparameters.h'
      INCLUDE 'Rodin2.h'

      Character  Hd*5,Tail*4,Nm*40
      Hd  ='PMcrs'
      Tail='.DAT' 
C                                     Open file on a tape
      OPEN(UNIT=9,FILE=npath//'PMcrd.DAT',
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
     +	               STATUS='UNKNOWN',RECL=NACCES)
      EndDo
 
      REWIND 9
      RETURN
      END
C---------------------------------------------
C                       Write current data to  disk/tape
C
      SUBROUTINE WRTAPE
C----------------------------------------------
      INCLUDE 'PMparameters.h'
      INCLUDE 'Rodin2.h' 
C                                       write header and control data
      WRITE (9) HEADER,
     +           AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +           TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     +           NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5,
     +           Ocurv,extras
      REWIND 9
      RETURN
      END
C--------------------------------------------------
C                             Write current PAGE of particles (x,v) to disk
C                             NRECL - length of ROW block in words
      SUBROUTINE WRIROW(IROW,Ifile)
C--------------------------------------------------
      INCLUDE 'PMparameters.h'
        WRITE (20+Ifile,REC=IROW) RECDAT
      RETURN
      END



C--------------------------------------------------------
      SUBROUTINE allhood(r1,r2,fnormalize)
C--------------------------------------------------------
     
      if(r2.le.1.)then
         fnormalize= r2*exp(-r2)
      else if(r1.ge.1.)then
         fnormalize = r1*exp(-r1)
      else 
          fnormalize = exp(-1.)
      endif

         return 
      end



C--------------------------------------------------------
C     This routine bins the disk density distribution in concentric rings
C     of size delta and height zd.
C     It also initializes the ellipsoid of velocities, according to Hernquist 1993.
C     Using the moments of boltzmann equation.
C------------------------------------------------------------  
      SUBROUTINE CDISKBIN       
C------------------------------------------------------------  
      INCLUDE 'PMparameters.h' 
      INCLUDE 'Rodin2.h' 
      Common/testF/  radi(Nmo) ,n(Nmo) 
      Common/zprof/zmin 
      External F,QSIGR,SIGPHI,d,AK2,OMEGAHD2 
      External cHsf,cBsf,DcHsf,DcBsf,aHD,aBInt 
      External alininterp,allhood  
      Real z_sigma(Nmo)
      Real*8  edelta,edelta1,dn0,d1,d2
c      real*8 dn0,radi

      write (*,*) '                    Disk Initialization Starts' 
      Id = 0  
      txcm =0. !center of mass components  
      tycm =0. 
      tzcm =0. 
      tz2 =0. 
      vrot =0.  
            
      Nbins = RT/delta           
      delta_z = 15.*delta*zd ! rd units
 
      Nbins_z = int(3.*zd/delta_z) + 1
      If(Nbins_z.le.1)then
       write(*,*)' Nbins_z < 1',Nbins_z
      Endif 

c      write(*,*)'Nbins_z =',Nbins_z,delta_z 


      Open(16,file=npath//'perfdiskbin.dat',status='unknown') 
      Open(17,file=npath//'diskin.dat',status='unknown') 
      Open(20,file=npath//'perfcildisp.dat',status='unknown')
      Open(11,file=npath//'table.dat',status='unknown')
     
c      OPEN(58,file='disktest.dat',status='unknown')    
      OPEN(68,file=npath//'trackgalaxy.dat',status='unknown')  
  
c      index = 0    

c
      Write(17,*)'  r (kpc/h)    Vrot   N(i)  Vrot(measured)   %error  
     *     Vcirc       Vc%error       vphirms/vphirm_th Vr_rms/Vr_rms_th  
     *   Vz_rms/Vz_rms_th  '   
c
      edelta = exp(-delta)
      edelta1 = 1.-edelta
      ind = 0 
      dn0 = 0.d0   
      Do i=1,Nbins   
         radi(i) = (i-0.5)*delta !The radial coordinate along the disk. 
         d1  = (i-1)*delta
         d2  = d1 + delta
c         dn0 = dn0 + radi(i)/EXP(radi(i)) !The normalization constant for n(i). 
         dn0 = dn0 + EXP(-d1)*(edelta1*(1.+d1)-delta*edelta) !The normalization constant for n(i). 
      Enddo   
      
    
       WRITE(68,*)'i r(i) N(i) N( i<)    N( i<)m      d(x)      %error'   
c      WRITE(*,*)'i r(i) N(i) N( i<)    N( i<)m      d(x)      %error'  
c
c      QSAV = 0.d0  
c
      Do i=1,Nbins    ! loop along the bins  starts 
         d1  = (i-1)*delta 
         d2  = d1 + delta 
         dummy =  EXP(-d1)*(edelta1*(1.+d1)-delta*edelta)/dn0 
c         dummy =radi(i)/EXP(radi(i))/dn0 
         n(i) = INT(Nd*dummy+0.5)  
c          write (*,'(2i7,8g11.4)')i,n(i),d1,d2,
c     &     radi(i)/EXP(radi(i))/dn0,dummy
        ind = ind + n(i)     
 
         radi(i) = radi(i) + 0.5*delta ! we go back to the edge of the bin  
                                !the term +delta/2.d0, is included because d describes the mass
                                !inside a given radius and  r(i) was defined at the bin center not at the edge.
         fr=d(radi(i))     
         aNm = (aMd/Nd)*ind  
         write (68,334)i,radi(i),n(i),ind,aNm,fr,       
     *        (aNm-fr)/fr 
c
c         IF(ABS(radi(i)-Rref).LT.delta)THEN     
c            write (68,*) ' Inside the reference bin'         
c            write (*,*) ' Inside the reference bin'            
c            write(*,*)radi(i),i ,RT   
c            QSAV = zcoef*3.36/TWOPI*QSIGR(radi(i))   
c            QSAV =  3.36*QSIGR(radi(i))      
c            write(*,*)QSAV,radi(i),i 
c         ENDIF  
c         write(*,*)'diskflag',i,Nbins 
c
      Enddo   ! The normalization constant for VR. ! loop along the bins  starts



 334  Format(i5,1x,f8.4,1x,i6,1x,i6,3(2x,1pe10.3))          
                     
c      CNORM =(QSAV*Q)*EXP(Rref/2.d0) ! The normalization constant for VR if Q = Q(r)    
                                
c        write(68,*)'CNORM=',CNORM,'sqrt(ek)',sqrt(AK2(Rref)),'Rref',Rref
c        write(*,*)'CNORM=',CNORM,'sqrt(ek)',sqrt(AK2(Rref)),'Rref',Rref
c
c      If(QSAV.le.1.e-10)Then  
c         write (68,*) ' Error in r and Rref. QSAV=0!' 
c         write(*,*)'diskflag2',radi(i),i,AK2(radi(i)) 
c
c         STOP
c      EndIf 
c
      write (11,*) ' rad  n_part V_c    Disp: rad   Z','         phi',  
     *'kappa Omega Q  <Vphi>  ASYMD  ' 
      write (68,*) ' rad  n_part V_c    Disp: rad     Z     phi' ,
     * '   kappa   Omega   Q   <Vphi>  ASYMD' 

               
      Do i=1,Nbins              !loop along the bins starts 
         rr = radi(i) 
         
         If(i.eq.1)then         !left edge of the bin   
            rr_p = 1.0e-7       ! delta/2.  
         else 
            rr_p = radi(i-1)
         endif  
          
         x = (rr+rr_p)/2.

c
	dzint = EXP(-rr)/(1.e0 - EXP(-RT)*(1.e0 + RT))
	dzint_p = EXP(-rr_p)/(1.e0 - EXP(-RT)*(1.e0 + RT))

 
        dzh = RMvirMd*(Rd/Rs)**2/Fcon/  
     *  (rr*Rd/Rs)/((rr*Rd/Rs)+1.0)**2*coef(1)*zd*Rd/Rs    
	dzb = aMassb/aMassd/Fhb*(Rd/Rs_b)**2/ 
     *  (rr*Rd/Rs_b)/((rr*Rd/Rs_b)+1.0)**3*coef(3)*zd*Rd/Rs_b 
c
c            dispVZ= SQRT(zd/2.e0*   
c     *           (1.0/EXP(rr)/(1.e0 - EXP(-RT)*(1.e0 + RT))+           
c     *        RMvirMd*(Rd/Rs)**2/Fcon*cHsf(rr)*coef(1)+   
c     *	      aMassb/aMassd/FHb*(Rd/Rs_b)**2*cBsf(rr))*coef(3))    

         dispVZ=  SQRT(zd/2.e0*dzint)

c-  
c     *   Rmvir*(Rd/Rs)/Fcon*cHdsigmaz(rr)*coef(1)-   
c     *   aMassb/aMassd/FHb*(Rd/Rs_b)**2*cBdsigmaz(rr)*coef(3))  

c	write(*,*)'test dispvz',dispVZ,zd/2.e0*dzint,	
c     *	Rmvir*(Rd/Rs)/Fcon*cHdsigmaz(rr) 

c	pause 

         dispVZ_p=  SQRT(zd/2.e0*dzint_p) 
c -    
c     *   Rmvir*(Rd/Rs)/Fcon*cHdsigmaz(rr_p)*coef(1) -   
c     *   aMassb/aMassd/FHb*(Rd/Rs_b)**2*cBdsigmaz(rr_p)*coef(3))        

c
c     Build and fill z-grid at radius rr  
c         Do iz=1,Nbins_z 
c             zmin = (iz-0.5)*delta_z  
c             z_sigma(iz) = SQRT(zd/2.e0*dzint-  
c     *   Rmvir*(Rd/Rs)/Fcon*cHdsigmaz(rr)*coef(1)-   
c     *   aMassb/aMassd/FHb*(Rd/Rs_b)**2*cBdsigmaz(rr)*coef(3))  
c
c 	write(*,*)'test z-grid',zmin,iz,zd,z_sigma(iz),Sqrt(zd/2.e0*dzint)  
c	pause  

c         Enddo    
c

c           write(*,*)'test_halo-Sigamz= ',        
c     *	RMvirMd*(Rd/Rs)/Fcon*cHdsigmaz(rr)*coef(1), 
c     *  aMassb/aMassd/FHb*(Rd/Rs_b)*cBdsigmaz(rr)*coef(3),  
c     *  zd/2.e0*EXP(-rr)/(1.e0 - EXP(-RT)*(1.e0 + RT)),rr,dispVZ         

c	  write(*,*)'test_Halo_BulgeSD= ',   
c     *    dzh,dzb,    
c     *    dzint,rr,(dzint + dzb+dzh)          
  
c           pause      
c         endif 
c
c     c     Sigma_R with a variable value  of  Q 

c            dispVR =     CNORM*dispVZ *sqrt(2.e0/zd)*   
c     *           (1.e0 - EXP(-RT)*(1.e0 + RT))             
c           CNORM*exp(-Sqrt(rr*rr + 2.d0*a*a)/2.d0)   
c
c     Sigma_R fixing the value of  Q  
         
            dispVR = 3.36*Q* QSIGR(rr) ! *(1.+1.55*exp(-x*Rd/0.200))  
            dispVR_p = 3.36*Q* QSIGR(rr_p) !*(1.+1.55*exp(-x*Rd/0.200))     

c

         dispVPR= dispVR*SIGPHI(rr)   ! phi dispersion      
         dispVPR_p= dispVR_p*SIGPHI(rr_p)   ! phi dispersion        

         Vcircular =sqrt(Vcd(rr+delta/2.)) 
c          Vcircular =sqrt(Vcd(rr))           

c         If(i.eq.1)then
            Vcircular_p =sqrt(Vcd(rr_p))      
c         else 
c            Vcircular_p =sqrt(Vcd(rr_p-delta/2.))     
c         endif
         


c	write(*,*)'Test smoothing',dispVPR/dispVPR_p,rr
c    	write(*,*)'Test smoothing',dispVR/dispVR_p,rr 
c	write(*,*)'Test smoothing',vcircular/vcircular_p,rr 
c	pause 

        Vcirc2 = sqrt(d(rr+delta/2.)+RMvirMd*f((rr+delta/2.)*Rd/Rs))   
         sigma_s= exp(-rr)/twopi/(1.e0 - EXP(-RT)*(1.e0 + RT))!*   
!     *        aMassd/(Rd*Rd)*(SM/cKpcm/cKpcm)   
C         sigma_s= exp(-rr)*zcoef/twop 
c         write(*,*)'sigma_s= ',sigma_s  
         Qq     = dispVR*sqrt(AK2(rr))/3.36/(sigma_s)!+
c     *     dzh +
c     *     dzb)  
 
         DVPHISTREAM = dispVR**2*  
     *        (1.d0 - SIGPHI(rr)**2   
     *        - 2.d0*rr) +  Vcircular**2  

         DVPHISTREAM_p = dispVR_p**2*  
     *        (1.d0 - SIGPHI(rr_p)**2   
     *        - 2.d0*rr_p) +  Vcircular_p**2      
  
c     *        - rr*DcHsf(rr)/max(cHsf(rr),1.e-7)*coef(1) 
c     *        - rr*DcBsf(rr)/max(cBsf(rr),1.e-7)*coef(3) ) 
c     *           +  Vcircular**2    

         If(DVPHISTREAM.LT.0.)write(68,*)' Negative DVPHISTREAM: ',
     *        DVPHISTREAM,rr,
     *       dispVPR* SIGPHI(rr), OMEGAHD2(rr)*rr**2, 
     *       dispVR**2
         DVPHISTREAM = Sqrt(Max(DVPHISTREAM,1.e-20))    
         DVPHISTREAM_p = Sqrt(Max(DVPHISTREAM_p,1.e-20))       

         Vdsc=sqrt(G*aMassd*SM/(Rd*ckpcm))*cKmcm!Scale to real Km/s
c                       write(*,*)'Vd scale',Vdsc  
         write (11,200)rr*Rd,n(i),Vcircular*Vdsc,  
     *              Vdsc*dispVR,Vdsc*dispVZ,Vdsc*dispVPR,sqrt(AK2(rr)),
     *               sqrt( OMEGAHD2(rr)),Qq,Vdsc*DVPHISTREAM,
     *          Vdsc*dispVR*sqrt(max( 2.*rr+Sigphi(rr)**2-1. ,1.e-10)),     
     *          sqrt(VC2(rr))*Vdsc
         factor = sqrt(aMassd/Rd*(G*SM)/cKpcm) 

         write (68,200)rr*Rd,n(i),Vcircular*Vdsc,
     *        Vdsc*dispVR,Vdsc*dispVZ,          
     *        Vdsc*dispVPR,sqrt(AK2(rr)),   
     *        sqrt( OMEGAHD2(rr)),Qq,Vdsc*DVPHISTREAM,
     *        Vdsc*sqrt(max(Vcircular**2  - DVPHISTREAM**2,1.e-10)),
     *        sqrt(VC2(rr))*Vdsc     

         diagnostic = DVPHISTREAM**2 -dispVR**2*(1.d0 - SIGPHI(rr)**2-
     *        2.d0*rr)  
c         if((diagnostic -Vcircular**2).gt.0.)then    
c            write(*,*)'Out of eq.',Vdsc*sqrt(diagnostic -Vcircular**2), 
c     *           radi(i)*rd,n(i)        
c            stop 
c         endif    
c         write(*,*)'Vcir-Vcirc2 = ',Vcircular- Vcirc2  
                      
c 200   format(f8.4,i5,1x,f7.2,3(1x,1pe10.3,1x),2(1x,f7.2,1x), 
c     *        4(1x,f9.2,1x)) 

         

 200   format(f8.4,i5,1x,f7.2,3(1x,f7.2,1x),2(1x,f7.2,1x),4(1x,f9.2,1x))
               r1 = radi(i) - delta
               r2 = radi(i)  
         CALL allhood(r1,r2,fnormalize)
 
         Do j=1,n(i) 
            Id =Id + 1
C
C     We generate the coordinates.   
C                           theta = phi = azimuthal angle              
            Theta = twopi*(randd(Nseed)) !Theta within the range from 0 : 2pi. 
c            costhetad = 2.*randd(Nseed)-1  
            

c            dx = radi(i)*cos(theta) !Cartesian coordinates on the disk. 
c            dy = radi(i)*sin(theta) 
             
c
 91         rad_new =  r1 + delta*randd(Nseed)
            Fran    =  randd(Nseed) 
            funcion =  rad_new*exp(-rad_new)/fnormalize
c            write(*,'(" test smoothing:",9G11.4)')funcion,Fran,rad_new,
c     &                                             r1,r2
            if(Fran.gt.funcion)goto 91 
c            write(*,'("         take  :",9G11.4)')rad_new
            dx = rad_new*cos(theta) !Cartesian coordinates on the disk.  
            dy = rad_new*sin(theta)   
            
c            dx = radi(i)*cos(theta) !Cartesian coordinates on the disk. 
c            dy = radi(i)*sin(theta) 
            

            xc(Id) = dx*Rd/Rs   
            yc(Id)= dy*Rd/Rs  
            peso(Id) = wspecies(1)     
c
 12         Fz = randd(Nseed)   ! The cumulative function for z. 
            If(Fz.ge.1.d0)goto 12 
            z = .5d0*dlog((1.d0+Fz)/(1.d0 - Fz)) !The vertical coordinate.
           If(z.gt.3.)goto12 
          sig1 = 2.d0*randd(Nseed)-1.d0 
           
          dz =sign(1.e0,sig1)*z*zd !z goes from -1 to 1 and 3*zd is the normalization.(See parameters.h,  zcoef)

c
            Index_z = int(abs(dz)/delta_z) + 1     
            dispVz2 = z_sigma(index_z)  

c	     write(*,*)'test vertical',dispVz2,Index_z,dz,rr
c             pause

c
            zc(Id) =dz*Rd/Rs  

            txcm = txcm + xc(Id) 
            tycm = tycm + yc(Id) 
            tzcm = tzcm + zc(Id) 
            tz2     = tz2 + zc(Id)**2  

  
C
C     Calculate The Velocity Ellipsoid

c         DVZ =  dispVZ *GAUSS(Nseed) 
          dispVZ2 = alininterp(dispVZ_p,dispVZ,rr_p,rr,rad_new) !neglecting changes in surf. density  
          dispVR2 = alininterp(dispVR_p,dispVR,rr_p,rr,rad_new)  
          dispVPR2 = alininterp(dispVPR_p,dispVPR,rr_p,rr,rad_new)   
          DVPHISTREAM2 = alininterp(DVPHISTREAM_p,DVPHISTREAM,rr_p, 
     *         rr,rad_new)      

    
          DVZ =  dispVZ2*GAUSS(Nseed)      ! if the Z dist is not isothermal

c          DVZ = dispVZ*GAUSS(Nseed)

         Vzc(Id) = DVZ/sqrt(RMvirMd/Fcon/Rs*Rd) 
         DVR =  dispVR2 *GAUSS(Nseed) 
c         DVR =  dispVR*GAUSS(Nseed) 


c         DVRt = CNORM*dispVZ*SQRT(2.d0/zd/zcoef)
c         If(abs(dispVR-DVRt)/dispVR.gt.0.01)Then
c          write (*,*) ' --Test DVR: ', dispVR,' =',
c     *             DVRt,abs(dispVR-DVRt)/dispVR
c          EndIf 

         DVPHIRAN = dispVPR2*GAUSS(Nseed)   

c         DVPHIRAN = dispVPR*GAUSS(Nseed)   

       
C
C     Now Convert to a Cartesian System of Coordiantes
C     First we calculate Vphi 

       DVPHI = DVPHIRAN + DVPHISTREAM2  
c       DVPHI = DVPHIRAN + DVPHISTREAM 

 70    dVX= DVR*cos(theta)-sin(theta)*DVPHI 
        VXc(Id) = dVX/sqrt(RMvirMd/Fcon/Rs*Rd) !in halo Units    
 
       dVY= DVR*sin(theta)+cos(theta)*DVPHI  
       VYc(Id)= dVY/sqrt(RMvirMd/Fcon/Rs*Rd) ! !in halo units 
 

   
C Now write useful quantities for analysis.

c         write(17,71)dispVz,dispVR,dispVPR,radi(i) ! The kinematics profile

         write(16,93)pMd,dX,dY,dZ, ! mi, x's, v's, sigma_t,omega
     *        dVX,dVY,dVZ,sqrt(DVPHI**2+DVR**2+DVZ**2),
     *        SQRT(OMEGAHD2( radi(i) )), !SQRT(OMEGA2( radi(i) ))- SQRT(AK2(radi(i)))/2.d0,
     *        SQRT(AK2(radi(i))),CNORM*dispVZ*SQRT(2.d0/zd/zcoef) 
     
C        Now write the positions and velocities in a file. 
c         write(58,73)pMd,dX,dY,dZ,dVX,dVY,dVZ 
        
          write(20,71)radi(i),DVR,DVPHIRAN,dVZ     
               
 93      FORMAT(11(1x,f15.6,1x))  
 73      FORMAT(73 (1x,f12.6,1x)) 
 71      FORMAT(4(1x,f10.6,1x))      
C
c         Measured Quantities from the realization for self-check           
c         vrot=vrot+ (dVy* dx -dVx* dy)/radi(i)      
          vrot2=vrot2 +  (Vyc(Id)*xc(Id) - Vxc(Id)* yc(Id))*Rs/Rd*     
     *        sqrt(RMvirMd/Fcon/Rs*Rd)/radi(i)                
           vrot = vrot +  DVPHISTREAM ! + DVPHIRAN           
          avphi =  avphi + sqrt(RMvirMd/Fcon/Rs*Rd)*
     *        ( vyc(Id)*xc(Id) -  vxc(Id)*yc(Id))*Rs/Rd/radi(i)
     

           avphi2 = avphi2 + (RMvirMd/Fcon/Rs*Rd)*
     *         (( vyc(Id)*xc(Id) -  vxc(Id)*yc(Id))*Rs/Rd/radi(i))**2 
          
           avpr =avpr  + sqrt(RMvirMd/Fcon/Rs*Rd)*
     *          (vxc(Id)*xc(Id) +  vyc(Id)*yc(Id))*Rs/Rd/radi(i)

           avpr2 = avpr2  +(RMvirMd/Fcon/Rs*Rd)*
     *          ((vxc(Id)*xc(Id) +  vyc(Id)*yc(Id))*Rs/Rd/radi(i))**2

            avz = avz + sqrt(RMvirMd/Fcon/Rs*Rd)*
     *          vzc(Id)

            avz2 = avz2 + (RMvirMd/Fcon/Rs*Rd)*
     *           vzc(Id)*vzc(Id)

c            if(id.eq.1)write(*,*)'vtest= ',vxc(id)*V_0,    
c     &           vyc(id)*V_0,  vzc(id)*V_0   
            

      enddo         

      sr = max((avpr2/max(n(i),1)- (avpr/max(n(i),1))**2),1.e-10)   
      sz = max((avz2/max(n(i),1)- (avz/max(n(i),1))**2),1.e-10)    
      sp = max((avphi2/max(n(i),1)- (avphi/max(n(i),1))**2),1.e-10)  
      if(sz.lt.0.or.sr.lt.0.or.sp.lt.0)then   
         write(*,*)'test detailed',avphi2,n(i),avphi,i
         write(*,*)'sz2 sr2 sp2 r n i =  ',sz,sr,sp,radi(i),n(i),i      
         stop  
      endif  


      write(17,105)radi(i)*rd,vrot*Vdsc/max(n(i),1),
     *     float(n(i)),vrot2*Vdsc/max(n(i),1),(vrot-vrot2)/vrot,    
     *     vcircular*Vdsc,
     *     (diagnostic-Vcircular**2)/(Vcircular**2)*100000,          
     *     sqrt(sp)/dispVPR , sqrt(sr)/dispVR ,sqrt(sz)/dispVZ    
c      write(*,*)'diskflag2',i
c      pause  
      vrot =0. 
      vrot2 =0.    
      avphi =0.
      avphi2 =0.  
      avpr =0. 
      avpr2 =0.
      avz =0.
      avz2 =0.
      enddo       
 105  format(11(1x,f10.5))   
      DMP = aMassd/Id 
      write(68,'(" Disk: Expected Mass-per-particle  Md/N=",2g11.3)')
     &           pMd,DMP 
      Ndr =Id
      write(*,*) ' Ndreal(# of disk particles)           = ',Id   
      lspecies(1) = Id          ! Include the disk particls in the 1st specie
      write(*,'(" Disk center of mass=",3g11.3," rms z_disk=",g11.3)')
     &    txcm/Id,tycm/Id,tzcm/Id,sqrt(max(tz2/Id-(tzcm/Id)**2,1.e-10))      
      close(16)
      close(17) 
      close(20) 
c      close(58)

      write (*,*) '                    Disk Initialization Done' 
      return
      end   
       

C----------------------------------------------------------
      REAL*4 Function d(x) ! This is the dimensionless disk mass. 
C---------------------------------------------------------- 
      INCLUDE 'PMparameters.h' 
      INCLUDE 'Rodin2.h'

       
c      write(*,*)'Flag d(x)',x,RT    
      If(x.le.RT)then  
         d = (1.0 - EXP(-X)*(1.0 + X))/        
     *        (1.0 - EXP(-RT)*(1.0 + RT)) 
      Else
         d =1.   
      Endif 
c     The normalization is made with respect the finite size of our disk
c     (1.d0 - 1.d0/EXP(RT)*(1.d0 + RT))* zcoef   not to the infinite disk 
c     This statement is wrong Zcoeff is cancelled out  because of the normalization factor (March 2 01)
c      d = (1.d0 - 1.d0/EXP(X)*(1.d0 + X))*zcoef   this the case of an infinite disk
       
      Return  
      End 
C-------------------------------------------------------------- 
      REAL*4 Function Mmd(x) !Disk contribution for the halo velocities.
C-----------------------------------------------------------
      INCLUDE 'PMparameters.h' 
      INCLUDE 'Rodin2.h'        
 
c      cMmd = (1.d0 - 1.d0*EXP(-X*Rs/Rd)*(1.d0 + X*Rs/Rd))*zcoef/
c     *     x**3/(1.d0+x)**2
      Mmd = d(x*Rs/(Rd))/x**3/(1.d0+x)**2        
         
c      write (*,*) '    mmd',Mmd-cMmd
      Return  
      End
C-----------------------------------------------------------  
      REAL*4 Function Mmdb(x) ! Disk contribution to the bulge velocities  
C-----------------------------------------------------------   
      INCLUDE 'PMparameters.h' 
      INCLUDE 'Rodin2.h'

      ww =  x*Rs_b/Rd  
      tcoef = (1.d0 - EXP(-RT)*(1.d0 + RT)) 

      Mmdb =  (1.d0 - EXP(-ww)*(1.d0 + ww))/  
     *     tcoef / (x*(x + 1.))**3              

      Return  
      End 
C-------------------------------------------------------------- 
      REAL*4 Function Dpot(x) ! Disk potential 
C-------------------------------------------------------------- 
      INCLUDE 'PMparameters.h' 
      INCLUDE 'Rodin2.h'
      External  bessi0,bessk0,bessi1,bessk1
      y=.5*x  
      dpot = y*(  bessi0(y) *bessk1(y)-bessi1(y) *bessk0(y) )    
      Return  
      End
C--------------------------------------------------------------
      REAL*4 function VC2(x)  !Disk Vcircular**2 
C--------------------------------------------------------------
      INCLUDE 'PMparameters.h' 
      INCLUDE 'Rodin2.h' 
      External  bessi0,bessk0,bessi1,bessk1 
      y=.5*x  
      VC2 =2.*y*y*( bessi0(y) *bessk0(y)-bessi1(y) *bessk1(y))   
   
      Return  
      End

C-------------------------------------------------------------- 
      REAL*4 function Vcd(x)!Total Vcircular**2  
C--------------------------------------------------------------
      INCLUDE 'PMparameters.h'  !Mvir/Md= 1.26d12/4.d10 =31.5
      INCLUDE 'Rodin2.h'

      External  bessi0,bessk0,bessi1,bessk1 

       y=.5*x 
       Vcd = 2.*y*y*( bessi0(y) *bessk0(y)-bessi1(y) *bessk1(y))  + 
     *     RMvirMd/Fcon*(F(x*Rd/Rs)/x)*coef(1)  +
     *     aMassb/aMassd*FHernq(x*Rd/Rs_b)/x*coef(3)   
 
      Return    
      End
C-------------------------------------------------------------- 
      REAL*4 function EK2(X)     !Square of the Epicyclic Frequency. 
C-------------------------------------------------------------- 
      INCLUDE 'PMparameters.h'  
      INCLUDE 'Rodin2.h' 

      EK2 = EXP(-X)/X/      
     *     (1.d0 - EXP(-RT)*(1.d0 + RT)) + 
     *     d(X)/X**3 + 
     *     RMvirMd/Fcon*((Rd/Rs)**2/X/
     *     (1.d0 + X*Rd/Rs)**2 +
     *     F(x*Rd/Rs)/X**3)*coef(1)  + 
     *     aMassb/aMassd/FHernq*((Rd/Rs_b)**2/X/
     *     (1.d0 + X*Rd/Rs_b)**2 +
     *     F(x*Rd/Rs_b)/X**3)*coef(3)       
       
      Return
      End
C---------------------------------------------------------------- 
      REAL*4 function OMEGA2(X)  !Square of the Circular Frequency.   
C----------------------------------------------------------------
      INCLUDE 'PMparameters.h'  !Mvir/Md= 1.26d12/4.d10 =31.5
      INCLUDE 'Rodin2.h'
         
       OMEGA2 =(d(X) + coef(1)*RMvirMd/Fcon*F(X*Rd/Rs))/X**3        
c       write(*,*)OMEGA2,x
      Return  
      End

C---------------------------------------------------------------- 
      REAL*4 function OMEGAHD2(X)  !Square of the Circular Frequency.   
C----------------------------------------------------------------  
      INCLUDE 'PMparameters.h'  !Mvir/Md= 1.26d12/4.d10 =31.5
      INCLUDE 'Rodin2.h'
      External F,FHernq, bessi0,bessk0,bessi1,bessk1 
 
      y=.5*x   
      
c     write(*,*)'flagOmHD2',x,bessk0(y)     
      OMEGAHD2=  .5*(bessi0(y) *bessk0(y)-bessi1(y) *bessk1(y)) +  
     *     RMvirMd/Fcon*(F(x*Rd/Rs)/x**3)* coef(1)   + 
     *     aMassb/aMassd*FHernq(x*Rd/Rs_b)/x**3* coef(3)                     

      Return   
      End     

C--------------------------------------------------------------  
      REAL*4 function AK2(X)     !Square of the Epicyclic Frequency.  
C--------------------------------------------------------------  
      INCLUDE 'PMparameters.h'  !Mvir/Md= 1.26d12/4.d10 =31.5 
      INCLUDE 'Rodin2.h'  
      External OMEGAHD2,F, bessk0,bessk1,bessi0,bessi1  

      y=.5*x  
      
      
c      write(*,*)'flagak2',F(x*Rd/Rs),OMEGAHD2(x) 
      

      aK2 = 4.*OMEGAHD2(x)  - 3.*RMvirMd/Fcon*(F(x*Rd/Rs)/x**3)* 
     *     coef(1) - 3.*aMassb/aMassd*(Fhernq(x*Rd/Rs_b)/x**3)*coef(3)+
     *     y*( bessi1(y)*bessk0(y)- bessi0(y)*bessk1(y) +bessi1(y)* 
     *     bessk1(y)/y) + 
     *     RMvirMd/Fcon*(Rd/Rs)**2*( 1.0/x/(1.0 + x*Rd/Rs)**2 )*coef(1)+ !anadimos un cuadrado a (Rd/Rs)        
     *     aMassb/aMassd*Fhernq(x*Rd/Rs_b)*2.*  
     *     ( 1.0/x*3/(1.0 + x*Rd/Rs) )*coef(3)                   

c     *     RMvirMd/Fcon*(Rd/Rs)**3*( 1.0/x/(1.0 + x*Rd/Rs)**2 )*coef(1)+   
C       (Rd/Rs)**3 was exchange by (Rd/Rs)   for halo contribution to aK2
      Return     
      End    
C--------------------------------------------------------------- 
      REAL*4 function SIGPHI(X) !Sqrt(sigmaphi**2/vr**2)eq. 2.26,Hernquist  
C--------------------------------------------------------------- 
      INCLUDE 'PMparameters.h' 
      INCLUDE 'Rodin2.h' 

      ee =AK2(X)
      eee=OMEGA2(X) 
c      write (*,*) ' x=',x,ee,eee

      SIGPHI = AK2(X)/4.0/OMEGAHD2(X)   
      SIGPHI = SQRT(SIGPHI)  
      Return 
      End 
 
C---------------------------------------------------------------- 
      REAL*4 function QSIGR(X) 
C---------------------------------------------------------------- 
      INCLUDE 'PMparameters.h' 
      INCLUDE 'Rodin2.h'    
      External AK2,cBsf,cHsf 
c      write(*,*)'flagqsigr',x,AK2(X)   
C     I corrected for the finite size of the disk:  (1.d0 - EXP(-RT)*(1.d0 + RT))   
      QSIGR = EXP(-X)/(1.d0 - EXP(-RT)*(1.d0 + RT)) /
     * SQRT(AK2(X))/6.2832  !+ 
c     *	aMassb/aMassd/Fhb*(Rd/Rs_b)**2/          !Bulge SF  
c     *  (x*Rd/Rs_b)/(( (x*Rd/Rs_b)+1.0)**3)*coef(3)*
c     *  zd*Rd/Rs_b*coef(3) +   
c     *  RMvirMd*(Rd/Rs)**2/Fcon/(x*Rd/Rs)/(((x*Rd/Rs)+1.0)**2)*
c     *  coef(1)*zd*Rd/Rs/   
c     *  SQRT(AK2(X))  
c     *  /6.2832        

c       write(*,*)'QSIGR',QSIGR  
      Return    
      End    
C---------------------------------------------------------------
      REAL*4 Function bessi0(x) ! Modified Bessel I_o Function  
C-----------------------------------------------------------  
      IMPLICIT REAL*4(A-H,O-Z) 
   
c      Real*4 x                  !bessi0, 
C     Returns the modified Bessel Function  I_o(x) for any real x
c      Real ax
      Real*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y ! Accumulate polinomials in double precision
      Save  p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      Data p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,
     *     1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      Data q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,
     *     0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,
     *     0.2635537d-1,-0.1647633d-1,0.392377d-2/ 
      
      if (abs(x).lt.3.75)then
         y=(x/3.75)**2
         bessi0 = p1 + y*(p2+y*(p3 + y*(p4 + y*(p5 + y*(p6 +y*p7)))))
      else
         ax = abs(x)
         y = 3.75/ax
         bessi0 = (exp(ax)/sqrt(ax))*(q1 + y*(q2 +y*(q3 + y*(q4 + 
     *        y*(q5 +y*(q6 + y*(q7 + y*(q8 + y*q9)))))))) 
      endif

      Return
      End 


C------------------------------------------------------------
      REAL*4 Function bessk0(x) ! Modified Bessel K_o Function 
C----------------------------------------------------------- 
 
C     Returns the modified Bessel Function  K_o(x) for positive real x  
      IMPLICIT REAL*4(A-H,O-Z)
c      Real*4 x                  
      external bessi0  
      Real*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      Save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7 
      Data p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,
     *     0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/ 
      Data q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,
     *     -0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/ 
      
      If(x.le.2.0)then  
         y = x*x/4.0
         bessk0 = (-log(x/2.0)* bessi0(x)) + (p1 + y*(p2 +y*(p3 + 
     *        y*(p4 + y*(p5 + y*(p6 + y*p7))))))
      Else 
         y = (2.0/x)
         bessk0 = ( exp(-x)/sqrt(x))*(q1 + y*(q2 + y*(q3 +
     *        y*(q4 + y*(q5 + y*(q6 + y*q7))))))
      Endif
      Return 
      End
C------------------------------------------------------------
      REAL*4 Function bessi1(x) ! Modified Bessel I_1 Function 
C----------------------------------------------------------- 
      IMPLICIT REAL*4(A-H,O-Z)
c      Real*4 x   
C     Returns the modified Bessel Function  I_o(x) for any real x 
c      Real ax
      Real*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y ! Accumulate polinomials in double precision
      Save  p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      Data p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,
     *     0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
      Data q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,
     *     -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,
     *     -0.2895312d-1,0.1787654d-1,-0.420059d-2/
      If(abs(x).lt.3.75)then
         y = (x/3.75)**2  
         bessi1 = x*(p1 + y*(p2 +y*(p3+ y*(p4+ y*(p5+ y*(p6 + y*p7)))))) 
      Else
         ax = abs(x) 
         y = 3.75/ax
         bessi1 = (exp(ax)/sqrt(ax))*(q1 + y*(q2 + y*(q3 + y*(q4 +
     *        y*(q5 + y*(q6 + y*(q7 + y*(q8 + y*q9))))))))
         if(x.lt.0.)bessi1= -bessi1
      Endif
      Return
      End

C------------------------------------------------------------
      REAL*4 Function bessk1(x) ! Modified Bessel K_1 Function 
C----------------------------------------------------------- 
      IMPLICIT REAL*4(A-H,O-Z)
c      Real*4 x   
c      Real*4 bessi1 
      external bessi1 
      Real*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      Save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7 
      Data p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,
     *     -0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
      Data q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1,
     *     0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/
      If(x.le.2.0)then
         y = x*x/4.0
         bessk1 = (log(x/2.0)*bessi1(x)) + (1.0/x)*(p1 + y*(p2+
     *        y*(p3+y*(p4+y*(p5 + y*(p6 + y*p7)))))) 
      Else
         y = 2.0/x  
         bessk1 = ( exp(-x)/sqrt(x))*(q1 + y*(q2+y*(q3 +
     *        y*(q4+y*(q5 + y*(q6 + y*q7))))))
      Endif
      Return
      End
  
      
 
c     -------------------------------------  
      subroutine Readconf ()  
c     -------------------------------------  
      Include 'PMparameters.h'   
      Include 'Rodin2.h'    
      write (*,'(/1x, ''****************************************'',
     &     /1x, ''* Rodin2 *'', 
     &     /1x, ''****************************************'')')
      
      open ( unit = 20 , file = npath//'Conf.R.dat' )   
       
      read  ( unit = 20 , fmt = '(10x,     a40)')Header  
      read  ( unit = 20 , fmt = '(70x,     i10)')Iconf
      read  ( unit = 20 , fmt = '(70x,     i10)')Isat 
      read  ( unit = 20 , fmt = '(70x,     i10)')IMM  
      read  ( unit = 20 , fmt = '(70x,     i10)')ISSat  
      read  ( unit = 20 , fmt = '(70x,     i10)')Nd
      read  ( unit = 20 , fmt = '(70x,   e10.3)')aMassd  
      read  ( unit = 20 , fmt = '(70x,   e10.3)')Rd  
      read  ( unit = 20 , fmt = '(70x,   e10.3)')Zdkp  
      read  ( unit = 20 , fmt = '(70x,   e10.3)')Q  
      read  ( unit = 20 , fmt = '(70x,   e10.3)')Rtkp
      read  ( unit = 20 , fmt = '(70x,   e10.3)')Rrefkp   
      read  ( unit = 20 , fmt = '(70x,   e10.3)')aMassb
      read  ( unit = 20 , fmt = '(70x,   e10.3)')C_b  
      read  ( unit = 20 , fmt = '(70x,   e10.3)')Rs_b   
      read  ( unit = 20 , fmt = '(70x,   e10.3)')Xout_b  
      read  ( unit = 20 , fmt = '(70x,   e10.3)')aMassh
      read  ( unit = 20 , fmt = '(70x,   e10.3)')Con   
      read  ( unit = 20 , fmt = '(70x,   e10.3)')Xout 
      read  ( unit = 20 , fmt = '(70x,   e10.3)')Box
      read  ( unit = 20 , fmt = '(70x,   e10.3)')aMass_sat 
      read  ( unit = 20 , fmt = '(70x,   e10.3)')Con_sat    
      read  ( unit = 20 , fmt = '(70x,   e10.3)')Xout_sat
      read  ( unit = 20 , fmt = '(70x,   e10.3)')D0x_sat 
      read  ( unit = 20 , fmt = '(70x,   e10.3)')D0y_sat 
      read  ( unit = 20 , fmt = '(70x,   e10.3)')D0z_sat  
      read  ( unit = 20 , fmt = '(70x,   e10.3)')Vx0_sat
      read  ( unit = 20 , fmt = '(70x,   e10.3)')Vy0_sat
      read  ( unit = 20 , fmt = '(70x,   e10.3)')Vz0_sat 
      read  ( unit = 20 , fmt = '(70x,   e10.3)')arot   
      read  ( unit = 20 , fmt = '(70x,   e10.3)')alpha1  
      read  ( unit = 20 , fmt = '(70x,   e10.3)')beta1 
      read  ( unit = 20 , fmt = '(70x,   e10.3)')alpha2 
      read  ( unit = 20 , fmt = '(70x,   e10.3)')beta2  
      read  ( unit = 20 , fmt = '(70x,   f10.6)') aexpn 
      read  ( unit = 20 , fmt = '(70x,   f10.6)') astep
      read  ( unit = 20 , fmt = '(70x,   e10.2)')  overdens
      read  ( unit = 20 , fmt = '(70x,   e10.2)') ocurv
      read  ( unit = 20 , fmt = '(70x,     i10)')  nseed
      read  ( unit = 20 , fmt = '(70x,   e10.2)') hubble
      read  ( unit = 20 , fmt = '(70x,   e10.2)')  Om0
      read  ( unit = 20 , fmt = '(70x,   e10.2)')  Oml0
      read  ( unit = 20 , fmt = '(70x,     i10)')  nspecies   
      read  ( unit = 20 , fmt = '(70x,   e10.2)')  Xb(1)
      read  ( unit = 20 , fmt = '(70x,   e10.2)')  Xb(2)
      read  ( unit = 20 , fmt = '(70x,   e10.2)')  Xb(3)
      read  ( unit = 20 , fmt = '(70x,   e10.2)')  Xb(4)
      read  ( unit = 20 , fmt = '(70x,   e10.2)')  Xb(5)
      read  ( unit = 20 , fmt = '(70x,   e10.2)')  Xb(6)
      read  ( unit = 20 , fmt = '(70x,   e10.2)')  Xb(7)
      
      Xb(Nspecies) = Xout
c      write(*,*)'Conf check',Rs_b 
      
      Return
      End 
C------------------------------------------
      Subroutine Configuration
C------------------------------------------
      INCLUDE 'PMparameters.h'     
      INCLUDE 'Rodin2.h'  

      
      If(Iconf.eq.0)then        !pure halo 
         coef(1)=1.
         do i=2,4 
            coef(i) = 0.
         enddo
         write(*,*)'Pure Halo'
      else if(Iconf.eq.1)then   !pure disk 
         coef(1)=0. 
         coef(2)=1.
         coef(3)=0.
         coef(4) = dble(Isat)
         write(*,*)'Pure Disk'
      else if(Iconf.eq.2)then   !disk + bulge  
         coef(1)=0.  
         coef(2)=1.  
         coef(3)=1.   
         coef(4) = dble(Isat)  
         write(*,*)'disk + bulge'      
      else if(Iconf.eq.3)then   !disk + halo  
         coef(1)=1.  
         coef(2)=1.
         coef(3)=0. 
         coef(4) = dble(Isat)  
         write(*,*)'halo + disk'  
      else if(Iconf.eq.4)then   !bulge + disk + halo
         coef(1)=1.
         coef(2)=1.
         coef(3)=1.  
         coef(4) = dble(Isat)
         write(*,*)'halo + disk + bulge'
      endif  
      Return
      End
         
   
C--------------------------------------------------------  
      SUBROUTINE Initialh(Nwant,xtrunc)  !   !  Set initial Conditions         
C--------------------------------------------------------   
C                  xtrunc = outer radius in units of Rs
C                  Indice1:  set equal 2 after bulge initialization 
C                  Indice2 : set equal 2 after main galaxy initialization    and before satellite
C                  Iconf:     0 pure halo, 1, pure disk, 2 disk-bulge, 3 disk-halo, 4 bulge + disk + halo
      INCLUDE 'PMparameters.h'  
      INCLUDE 'Rodin2.h' 

c      PARAMETER(pi =3.1415926535d+0)  
      DATA Indice1/1/   
      DATA Indice2/1/
      SAVE Indice1,Indice2  
      
      PARAMETER(NradM=15000)  
      REAL*4    Radd(NradM),AINTG,Mmx,Mmd,Mmhb,Mmdb,
     *     aMmxh 
      Dimension Nr(NradM)  
      External  Mmx,AINTG,Mmd,FHernq,aHdispr,  
     *     aMmxh,Mmhb,Mmdb 
 
C                             Set up radii and N orbits
  
      
      hr  =xtrunc/Nrad   ! based on the system size 
      n0  = 0 
c      write(*,*)'flag ',Indice1,Indice2,Iconf     
      If(Indice1.eq.1.and.coef(3).eq.1)then  ! If this is the first call and there is a bulge 
         write(*,*)'The profile is  Hernquist (Bulge)'     
         Fwant = Nwant/FHernq(xtrunc) ! this is the mas per particle in internal units   
c         write(*,*)'Fwant= ', Fwant,FHernq(xtrunc)    
      else   
         write(*,*)'The profile isNFW'  
            Fwant =Nwant/F(xtrunc)  
c            write(*,*)'Fwant= ', Fwant    
         endif   

      write (16,200)  
 200  format(' Bin  radius  Npart  Expected N')  
 250  format(i5,f10.3,i5,g11.3) 

      Do i=1,Nrad 
         Radd(i) = (i-0.5)*hr 
         If(Indice1.eq.1.and.coef(3).eq.1)then   ! If there is a bulge(coef(3)  and has not been built yet      
            fshould = Fwant*(FHernq(Radd(i)+0.5*hr)- 
     *           FHernq(Radd(i)-0.5*hr)) 
         else
            fshould = Fwant*(F(Radd(i)+0.5*hr)-F(Radd(i)-0.5*hr)) 
c            write(*,*)Fwant,fshould
         endif    

         Nr(i)   = INT(fshould +0.5)   
         n0      = n0 + Nr(i) 
         If(i.le.20.or.i.gt.Nrad-15)Then  
           write (68,250) i,Radd(i),Nr(i),fshould  
           write (16,250) i,Radd(i),Nr(i),fshould  
         EndIf 
      EndDo 
      Norb =0  
      write (*,*) ' Nwant=',Nwant,' Nreal=',n0 
      write (16,*) ' Nwant=',Nwant,' Nreal=',n0 
       iorb  =lspecies(1)    
c       pause 

      Do i=1,Nrad                             ! Loop over all bins             
        If(Nr(i).ne.0)Then 
          x      = Radd(i) 
          xmax   = max(x,xtrunc)          
 

c             write(*,*)'No bulge not satellite  or bulge was done'   
             
             at = AINTG(Mmx,x,xmax)  ! halo self contribution  
             bt = AINTG(Mmd,x,xmax)*Fcon/RMvirMd*coef(2) ! disk contribution 
             ct = AINTG(aMmxh,x,xmax)*Fcon*aMassb/aMassh*coef(3) ! bulge contribution  if it was done 
c             write(*,*)'flag2'  
             sig_r2 = x*(1.+x)**2*(at + bt + ct)                 
             Vescap2= 2.*(LOG(1.+x)/x +  
     *            (1. - Exp(-x*Rs/Rd) )/x*Fcon/RMvirMd*coef(2))!+  
c     *            2./(x*R0_b+1.)*Fcon*aMassb/aMassh*coef(3)  ) !square of escape velocity    
c*******************************************************************

          
c Si se quiere quitar el precooling, se tienen que comentar las siguientes 4 lineas dejar la 5 y comentartambién la 6  
             if(x.le.(1./Rslarge))then    
             sig_r  = sqrt(sig_r2)*(.6*Rslarge*x+ .4) ! rms of radial velocity   
c     write(*,*)"flag",sig_r ,sqrt(sig_r2),x 
          else  
             sig_r  = sqrt(sig_r2) ! rms of radial velocity  
          endif    
  
         Do j=1,Nr(i) 
            iorb =iorb +1  
c            write(*,*)'Flag1= ', iorb 
c            pause
            If(iorb.gt.Nmo)Then 
               write (*,*) ' Too many particles:',iorb,Nmo  
               STOP 
            EndIf  

 70         vxc(iorb) =sig_r*gauss(Nseed) !  No bulge, not satellite  or bulge was done  
               vyc(iorb) =sig_r*gauss(Nseed)  
               vzc(iorb) =sig_r*gauss(Nseed)    
            
               Vv       =vxc(iorb)**2 +
     &              vyc(iorb)**2 + 
     &              vzc(iorb)**2
               If(Vv.gt.Vescap2) goto 70
c     write (*,*) ' Escape: r=',x,
c     &             ' V=',sqrt(Vv),' Vesc=',sqrt(Vescap2)
               cc    = 2.*Randd(Nseed)-1.
               costh = sqrt(abs(1.-cc**2))
c     theta = asin(cc)
c     if(theta.lt.0.d+0)theta =theta +pi
c     theta = pi*Randd(Nseed)
               phi   = 2.*pi*Randd(Nseed) 
               rsin  = x*cc     !           SIN(theta)
               rcos  = x*costh  ! COS(theta) 
                                ! (Halo) 
                  xc(iorb) =rcos*COS(phi) 
                  yc(iorb) =rcos*SIN(phi)  
                  zc(iorb) =rsin   
               
               txcmh = txcmh + xc(iorb) 
               tycmh = tycmh + yc(iorb)
               tzcmh = tzcmh + zc(iorb)    
               peso(iorb) = 1.    
c     write(*,*)'flaginit',peso(iorb),iorb
c               write(*,*)'flag',xc(iorb), yc(iorb),zc(iorb),iorb   
            EndDo
         EndIf 
           
         !Main halo
            Norb = iorb   
         
      
      EndDo 

c      write(*,*)'Flag2= ', iorb 
              
      
   
      write (*,*) ' Number of Halo particles=',n0 
      write (16,*) ' Number of Halo particles=',Norb 
      
      twh = Norb  
     



      Return 
      End


C--------------------------------------------------------  
      SUBROUTINE Initialb(Nwant,xtrunc)  !   !  Set initial Conditions         
C--------------------------------------------------------   
C                  xtrunc = outer radius in units of Rs
C                  Indice1:  set equal 2 after bulge initialization 
C                  Indice2 : set equal 2 after main galaxy initialization    and before satellite
C                  Iconf:     0 pure halo, 1, pure disk, 2 disk-bulge, 3 disk-halo, 4 bulge + disk + halo
      INCLUDE 'PMparameters.h'    
      INCLUDE 'Rodin2.h'  

c      PARAMETER(pi =3.1415926535d+0) 
      DATA Indice1/1/ 
      DATA Indice2/1/
      SAVE Indice1,Indice2  
      
      PARAMETER(NradM=15000) 
      REAL*4    Radd(NradM),AINTG,Mmx,Mmd,Mmhb,Mmdb,
     *     aMmxh 
      Dimension Nr(NradM)  
      External  Mmx,AINTG,Mmd,FHernq,aHdispr,  
     *     aMmxh,Mmhb,Mmdb 
 
C                             Set up radii and N orbits 
  
      OPEN(67,FILE=npath//'bulge_control.dat',status='unknown')    
      OPEN(68,FILE=npath//'bulge_disp.dat',status='unknown')       


      hr  =xtrunc/Nrad   ! based on the system size    
      n0  = 0 
c      write(*,*)'flag ',Indice1,Indice2,Iconf     
      If(Indice1.eq.1.and.coef(3).eq.1)then ! If this is the first call and there is a bulge  
         write(*,*)'The profile is  Hernquist (Bulge)'      
         Fwant = Nwant/FHernq(xtrunc) ! this is the mas per particle in internal units    
         write(*,*)'Fwant= ', Fwant,FHernq(xtrunc)       
      else    
         write (*,*)'The profile isNFW'    
         Fwant =Nwant/F(xtrunc)  
         write(*,*)'Fwant= ', Fwant    
      endif     

      write (67,200)  
 200  format(' Bin  radius  Npart  Expected N')  
 250  format(i5,f10.3,i5,g11.3) 

      Do i=1,Nrad 
         Radd(i) = (i-0.5)*hr 
         If(Indice1.eq.1.and.coef(3).eq.1)then   ! If there is a bulge(coef(3)  and has not been built yet      
            fshould = Fwant*(FHernq(Radd(i)+0.5*hr)- 
     *           FHernq(Radd(i)-0.5*hr)) 
         else 
            fshould = Fwant*(F(Radd(i)+0.5*hr)-F(Radd(i)-0.5*hr)) 
c            write(*,*)Fwant,fshould
         endif    

         Nr(i)   = INT(fshould +0.5)   
         n0      = n0 + Nr(i) 
c         If(i.le.20.or.i.gt.Nrad-15)Then  
c           write (*,250) i,Radd(i),Nr(i),fshould
         write (67,250) i,Radd(i),Nr(i),fshould      
c         EndIf 
      EndDo 
      Norb =0  
      write (*,*) ' Nwant=',Nwant,' Nreal=',n0 
      write (67,*) ' Nwant=',Nwant,' Nreal=',n0 
       iorb  =lspecies(1)    
c       pause 

      Do i=1,Nrad                             ! Loop over all bins             
        If(Nr(i).ne.0)Then 
          x      = Radd(i) 
          xmax   = max(x,xtrunc)           


c          write(*,*)'for a bulge not initialized yet',i,radd(i),Nr(i)        
   
             at =AINTG(Mmhb,x,xmax)*aMassh/Fcon/aMassb*coef(1) !halo contribution  to bulge vrms       
             bt = AINTG(Mmdb,x,xmax)*Fcon/RMvirMd*coef(2) ! disk contribution to bulge vrms            
             ct = aHdispr(x)  !bulge self-contribution  to vrms                       
              

             sig_r2 = (at+bt )*x*(x+1.)**3 + ct ! Hernquist rms of radial velocity         
            Vescap2 = ( 2./(x+1.) + 2.*(LOG(1.+x*Rs_b/Rs)/(x*Rs_b/Rs))*  
     *            aMassh/aMassb/Fcon*coef(1)* Rs_b/Rs +  
     *            (1. - Exp(-x*Rs_b/Rd) )/x*aMassd/aMassb*Rs_b/Rd*
     *           coef(2) )*V0_b  !square of escape velocity                
              

             


  
c            if(x.le.(.3*Rs_b))then    ! if(x.le.(1./Rs_b))then 
c             sig_r  = sqrt(sig_r2)*(.6*Rs_b*x+ .4)!*V0_b ! rms of radial velocity     
c     write(*,*)"flag",sig_r ,sqrt(sig_r2),x  
c          else  
             sig_r  = sqrt(sig_r2) ! rms of radial velocity   
c          endif   
             
             write (67,711)x*Rs_b,sig_r*208.e-5*V0_b*    
     *            sqrt(aMassh/Fcon/Rs),sqrt(x*(1.+x)**2*at)*V0_b*  
     *            208.e-5*sqrt(aMassh/Fcon/Rs),    
     *            sqrt(x*(1.+x)**2*bt)*V0_b*208.e-5* 
     *            sqrt(aMassh/Fcon/Rs),sqrt(ct)*V0_b*
     *            208.e-5*sqrt(aMassh/Fcon/Rs) 


 711         format (5(1x,f10.4))   



c          endif     
             vdisp =0.
             vprom =0. 
  
         Do j=1,Nr(i) 
            iorb =iorb +1  
         
            If(iorb.gt.Nmo)Then 
               write (*,*) ' Too many particles:',iorb,Nmo  
               STOP 
            EndIf  

 
 
 70         vxc(iorb) =sig_r*gauss(Nseed)!*V0_b  
                  vyc(iorb) =sig_r*gauss(Nseed)!*V0_b  
                  vzc(iorb) =sig_r*gauss(Nseed)!*V0_b    

c                  write(*,*)vxc(iorb)*V0_b*208.e-5*sqrt(aMassh/Fcon/Rs),
c     *                 iorb, vxc(iorb)*208.e-5*sqrt(aMassh/Fcon/Rs)
c                  pause    


                  vxc(iorb) =vxc(iorb)*V0_b  
                  vyc(iorb) =vyc(iorb)*V0_b   
                  vzc(iorb) =vzc(iorb)*V0_b     


               Vv       =vxc(iorb)**2 +
     &              vyc(iorb)**2 + 
     &              vzc(iorb)**2
               If(Vv.gt.Vescap2) goto 70  
c               write (*,*) ' Escape: r=',x, 
c     &              ' V=',sqrt(Vv),' Vesc=',sqrt(Vescap2)   
               cc    = 2.*Randd(Nseed)-1.
               costh = sqrt(abs(1.-cc**2))
c     theta = asin(cc)
c     if(theta.lt.0.d+0)theta =theta +pi
c     theta = pi*Randd(Nseed)
               phi   = 2.*pi*Randd(Nseed) 
               rsin  = x*cc     !           SIN(theta)
               rcos  = x*costh  ! COS(theta) 



c                  write(*,*)vxc(iorb)*208.e-5*sqrt(aMassh/Fcon/Rs),
c     *                 iorb  
c                  pause       



                                ! Bulge 
               xc(iorb) =rcos*COS(phi)/Rslarge  
               yc(iorb) =rcos*SIN(phi)/Rslarge  
               zc(iorb) =rsin/Rslarge  
               
               txcmh = txcmh + xc(iorb)     
               tycmh = tycmh + yc(iorb) 
               tzcmh = tzcmh + zc(iorb)     
               peso(iorb) = 1.     
c     write(*,*)'flaginit',peso(iorb),iorb
c               write(*,*)'flag',xc(iorb), yc(iorb),zc(iorb),iorb   


               vdisp = vdisp + vxc(iorb)**2 

               
               vprom = vprom + vxc(iorb)
               

            EndDo  

            write(68,*)Sqrt(vdisp/Nr(i))*208.e-5*sqrt(aMassh/Fcon/Rs), 
     *           vprom/Nr(i),  
     *           sig_r*208.e-5*V0_b*sqrt(aMassh/Fcon/Rs), 
     *           x*Rs_b,Nr(i)  

         EndIf  
         
        ! Bulge 
         Nbulge = iorb -lspecies(1)    
         Nbulge2 = Nbulge      
      
         EndDo                  !large loop over the bins is done  
         close(68) 

      write (*,*) ' Bulge Number of particles=',Nbulge 
      write (67,*) ' Bulge Number of particles=',Nbulge 
      lspecies(1) =  iorb    
      twh = Norb  
       write (*,*) ' After bulge generation lspecies(1)= ',lspecies(1)    

c      write(*,*)'flag-test',xc(iorb),yc(iorb),zc(iorb)  


      Return 
      End 


C--------------------------------------------------------  
      SUBROUTINE Initialsat(Nwant,xtrunc)  !   !  Set initial Conditions            
C--------------------------------------------------------   
C                  xtrunc = outer radius in units of Rs 
C                  Indice1:  set equal 2 after bulge initialization  
C                  Indice2 : set equal 2 after main galaxy initialization    and before satellite  
C                  Iconf:     0 pure halo, 1, pure disk, 2 disk-bulge, 3 disk-halo, 4 bulge + disk + halo  
      INCLUDE 'PMparameters.h'     
      INCLUDE 'Rodin2.h'   

c      PARAMETER(pi =3.1415926535d+0) 
      DATA Indice1/1/ 
      DATA Indice2/1/ 
      SAVE Indice1,Indice2   
      
      PARAMETER(NradM=25000)
      REAL*4    Radd(NradM),AINTG,Mmx,Mmd,Mmhb,Mmdb, 
     *     aMmxh 
      Dimension Nr(NradM)  
      External  Mmx,AINTG,Mmd,FHernq,aHdispr,  
     *     aMmxh,Mmhb,Mmdb 
 
C                             Set up radii and N orbits
  
      
      hr  =xtrunc/Nrad   ! based on the system size  
      n0  = 0 
     
      If(Indice1.eq.1.and.coef(3).eq.1)then  ! If this is the first call and there is a bulge   
         write(*,*)'The profile is  Hernquist (Bulge)'     
         Fwant = Nwant/FHernq(xtrunc) ! this is the mas per particle in internal units   
c         write(*,*)'Fwant= ', Fwant,FHernq(xtrunc)      
      else   
         write(*,*)'The profile isNFW'   
            Fwant =Nwant/F(xtrunc)   
c            write(*,*)'Fwant= ', Fwant    
         endif    

      write (16,200)   
 200  format(' Bin  radius  Npart  Expected N')   
 250  format(i5,f10.3,i5,g11.3)  

      Do i=1,Nrad    
         Radd(i) = (i-0.5)*hr 
         If(Indice1.eq.1.and.coef(3).eq.1)then   ! If there is a bulge(coef(3)  and has not been built yet      
            fshould = Fwant*(FHernq(Radd(i)+0.5*hr)- 
     *           FHernq(Radd(i)-0.5*hr)) 
         else  
            fshould = Fwant*(F(Radd(i)+0.5*hr)-F(Radd(i)-0.5*hr))   
c             write(*,*)Fwant,fshould  
         endif    
    
         Nr(i)   = INT(fshould +0.5)   
         n0      = n0 + Nr(i) 
         If(i.le.20.or.i.gt.Nrad-15)Then  
c           write (*,250) i,Radd(i),Nr(i),fshould
           write (16,250) i,Radd(i),Nr(i),fshould  
         EndIf 
      EndDo 
      Norb =0  
      write (*,*) ' Nwant=',Nwant,' Nreal=',n0    
      write (16,*) ' Nwant=',Nwant,' Nreal=',n0   

     
      iorb  =lspecies(Nspecies)     
      
c      vavx = 0. 
c      vavy = 0.     
c      vavz = 0.   
  

      Do i=1,Nrad               ! Loop over all bins             
         If(Nr(i).ne.0)Then 
            x      = Radd(i) 
            xmax   = max(x,xtrunc)            


c******************************************************************* 
            ! for a satellite  not  yet  initialized  
             
c             write(*,*)'for a satellite  not  yet  initialized'

             at = AINTG(Mmx,x,xmax)      
             sig_r2 = x*(1.+x)**2*at 
             Vescap2= 2.*LOG(1.+x)/x*V0_sat!square of escape velocity       
 

             if(x.le.(1./Rslarge))then    
             sig_r  = sqrt(sig_r2)*(.6*Rslarge*x+ .4) ! rms of radial velocity   
c     write(*,*)"flag",sig_r ,sqrt(sig_r2),x 
          else  
             sig_r  = sqrt(sig_r2) ! rms of radial velocity  
          endif    
       
          

         Do j=1,Nr(i)  
            iorb =iorb +1    

            If(iorb.gt.Nmo)Then  
               write (*,*) ' Too many particles:',iorb,Nmo   
               STOP  
            EndIf    

                                ! for a satellite  not  yet  initialized  !Indice2=2  
 70         vxc(iorb) =sig_r*gauss(Nseed)*V0_sat    
            vyc(iorb) =sig_r*gauss(Nseed)*V0_sat     
            vzc(iorb) =sig_r*gauss(Nseed)*V0_sat   


       
               Vv       =vxc(iorb)**2 +     
     &              vyc(iorb)**2 + 
     &              vzc(iorb)**2 
               If(Vv.gt.Vescap2) goto 70 
c     write (*,*) ' Escape: r=',x,
c     &             ' V=',sqrt(Vv),' Vesc=',sqrt(Vescap2)

               vxc(iorb) =  vxc(iorb) !+Vx0_sati
               vyc(iorb) =  vyc(iorb) !+ Vy0_sati  
               vzc(iorb) =  vzc(iorb) !+ Vz0_sati  


               vavx =  vavx + vxc(iorb) 
               vavy =  vavy + vyc(iorb)  
               vavz =  vavz + vzc(iorb) 

             
               cc    = 2.*Randd(Nseed)-1.   
               costh = sqrt(abs(1.-cc**2))   
c     theta = asin(cc)   
c     if(theta.lt.0.d+0)theta =theta +pi
c     theta = pi*Randd(Nseed)
               phi   = 2.*pi*Randd(Nseed) 
               rsin  = x*cc     !           SIN(theta)  
               rcos  = x*costh  ! COS(theta)  
c                 satellite  
                  xc(iorb) =rcos*COS(phi)/Rslarge  ! - Dsatix
                  yc(iorb) =rcos*SIN(phi)/Rslarge  !  - Dsatiy
                  zc(iorb) =rsin/Rslarge   !- Dsatiz 

               
                  txcmh = txcmh + xc(iorb)      
                  tycmh = tycmh + yc(iorb) 
                  tzcmh = tzcmh + zc(iorb)     
                  peso(iorb) = 1.    

               EndDo  
            EndIf  
            
                                ! Satellite     
            Nsat = iorb -lspecies(Nspecies)        

         EndDo 

         
     
         write (*,*) ' Sat. Number of particles=',Nsat 
         write (16,*) ' Sat. Number of particles=',Nsat   
         

         ialmacen = ialmacen + Nsat 
         lspecies(1) =  lspecies(1)+Nsat   !Modified 16-11-2011, = iorb  
         write(*,*)lspecies(1)
         Do ij =2,Nspecies 
            lspecies(ij) = lspecies(ij) + Nsat  
            write(*,*)lspecies(ij)
         Enddo 
         

         twh = Norb  
         write(*,*)"Npart in sats",ialmacen 
        

c         write(*,*) 'vx= ',vavx/ Nsat *V_0          
c         write(*,*) 'vy= ',vavy/ Nsat *V_0            
c         write(*,*) 'vz= ',vavz/ Nsat *V_0          

      
         Return 
         End  

 
C--------------------------------------- 
      SUBROUTINE Shift( )                
C---------------------------------------  
C     
      INCLUDE 'PMparameters.h'  
      INCLUDE 	'Rodin2.h'    

      testx= 0.
      testy= 0.   
      testz= 0.
      tempx =0.  
      testpx= 0.
      testpy= 0.   
      testpz= 0.

      testgx= 0.
      testgy= 0.   
      testgz= 0.
      ilabg =0.  
      testgpx= 0.
      testgpy= 0.   
      testgpz= 0.
                                                                                 

      If(iMM.eq.0)then     
      Gal_shiftx = Dsatix*(aMass_sat/aMass_vir)
      Sat_shiftx = Dsatix*((aMass_vir-aMass_sat)/
     *     aMass_vir)    
       
      Gal_shifty = Dsatiy*(aMass_sat/aMass_vir) 
      Sat_shifty = Dsatiy*((aMass_vir-aMass_sat)/
     *     aMass_vir)   

      Gal_shiftz = Dsatiz*(aMass_sat/aMass_vir)
      Sat_shiftz = Dsatiz*((aMass_vir-aMass_sat)/
     *     aMass_vir)  


      Gal_shiftvx = Vx0_sati*(aMass_sat/aMass_vir) 
      Sat_shiftvx = Vx0_sati*((aMass_vir-aMass_sat)/
     *     aMass_vir)     
      
      Gal_shiftvy = Vy0_sati*(aMass_sat/aMass_vir) 
      Sat_shiftvy = Vy0_sati*((aMass_vir-aMass_sat)/
     *     aMass_vir)    

      Gal_shiftvz = Vz0_sati*(aMass_sat/aMass_vir) 
      Sat_shiftvz = Vz0_sati*((aMass_vir-aMass_sat)/ 
     *     aMass_vir)    
      

      else 
         Gal_shiftx = Dsatix*.5
         Sat_shiftx = Dsatix*.5
         
         Gal_shifty = Dsatiy*.5
         Sat_shifty = Dsatiy*.5
         
         Gal_shiftz = Dsatiz*.5
         Sat_shiftz = Dsatiz*.5


         Gal_shiftvx = Vx0_sati*.5
         Sat_shiftvx = Vx0_sati*.5
         
         Gal_shiftvy = Vy0_sati*.5
         Sat_shiftvy = Vy0_sati*.5   
         
         Gal_shiftvz = Vz0_sati*.5
         Sat_shiftvz = Vz0_sati*.5
      

         Do i=1,lspecies(Nspecies) ! duplicating the model  
            xc(i+lspecies(Nspecies)) =  xc(i) 
            yc(i+lspecies(Nspecies)) =  yc(i) 
            zc(i+lspecies(Nspecies)) =  zc(i) 

            vxc(i+lspecies(Nspecies)) =  vxc(i)  
            vyc(i+lspecies(Nspecies)) =  vyc(i)  
            vzc(i+lspecies(Nspecies)) =  vzc(i)  

            peso(i+lspecies(Nspecies)) = peso(i) 
            
            Nsat = lspecies(Nspecies)   
         enddo  

      endif

      write(*,*)'Shift xamplitud',Gal_shiftx,Sat_shiftx,Dsatix,D0x_sat  
      write(*,*)'Shift yamplitud',Gal_shifty,Sat_shifty,Dsatiy,D0y_sat  
      write(*,*)'Shift zamplitud',Gal_shiftz,Sat_shiftz,Dsatiz,D0z_sat   
      write(*,*)'Vx_Gal,Vx_Sat',Gal_shiftvx,Sat_shiftvx, V_0   
      write(*,*)'Vy_Gal,Vy_Sat',Gal_shiftvy,Sat_shiftvy ,V_0   
      write(*,*)'Vz_Gal,Vz_Sat',Gal_shiftvz,Sat_shiftvz,  V_0  
      
      If(arot.gt.0.)then 
      Call Rotx(0,alpha1) 
      Call Roty(0,beta1) 
      Call Rotx(1,alpha2)         
      Call Roty(1,beta2)       
      Endif  
      itemp = lspecies(Nspecies)*(1 + iMM)     
      lspecies(Nspecies) =  lspecies(Nspecies)*(1 + iMM)       
      icounter =0     

        



      Do i=1,itemp   

         If(i.gt.(itemp-Nsat))then       !displace the satellite/companion  
            xc(i) =  xc(i) + Sat_shiftx  
            yc(i) =  yc(i) + Sat_shifty    
            zc(i) =  zc(i) + Sat_shiftz  

            vxc(i) =  vxc(i) + Sat_shiftvx  
            vyc(i) =  vyc(i) + Sat_shiftvy     
            vzc(i) =  vzc(i) + Sat_shiftvz    

            testx = testx +  vxc(i)
            testy = testy +  vyc(i)    
            testz = testz +  vzc(i)
            tempx =  tempx + vxc(i)  

            testpx = testpx +  xc(i)
            testpy = testpy +  yc(i)    
            testpz = testpz +  zc(i)

            icounter =  icounter + 1 

         else   ! the primary system 
            xc(i) =  xc(i) - Gal_shiftx 
            yc(i) =  yc(i) - Gal_shifty 
            zc(i) =  zc(i) - Gal_shiftz 

            vxc(i) =  vxc(i) - Gal_shiftvx   
            vyc(i) =  vyc(i) - Gal_shiftvy
            vzc(i) =  vzc(i) - Gal_shiftvz   

            testgx = testx +  vxc(i)
            testgy = testy +  vyc(i)    
            testgz = testz +  vzc(i)
            ilabg =  ilabg + 1 

            testgpx = testgpx +  xc(i)
            testgpy = testgpy +  yc(i)    
            testgpz = testgpz +  zc(i)
         endif  

          

      Enddo 
 
     
     
      write(*,*)'vtest-sat',testx /Nsat*V_0 ,Sat_shiftvx ,Sat_shiftvy, 
     &     Sat_shiftvz               

       write(*,*)'rtest-sat',testpx/Nsat*R_0,testpz/Nsat*R_0,
     & testpz/Nsat*R_0

	write(*,*)'vtest-g',testy/ilabg*V_0,Gal_shiftvy,Gal_shiftvx,
     &	  Gal_shiftvz 				
      

	write(*,*)'rtest-g',testgpx/ilabg*R_0,testgpy/ilabg*R_0,
     &	testgpz/ilabg*R_0
      return 
      end


C---------------------------------------  
      SUBROUTINE Shift_cent( )                
C--------------------------------------- 
      INCLUDE 'PMparameters.h'   
      INCLUDE 'Rodin2.h' 

      Do i=1,lspecies(Nspecies)  

          xc(i) = xc(i) -Xcmt/R_0 
          yc(i) = yc(i) -Ycmt/R_0  
          zc(i) = zc(i) -Zcmt /R_0     

          vxc(i) = vxc(i) -VXcmt/R_0 
          vyc(i) = vyc(i) -VYcmt/R_0  
          vzc(i) = vzc(i) -VZcmt /R_0     

      Enddo

      return
      end     
C--------------------------------------- 
      REAL*4 Function FHernq(x) ! x =r/r_s,  F= mass(r)/M  for a Hernquist profile 
C--------------------------------------- 
      INCLUDE 'PMparameters.h'  
      INCLUDE 'Rodin2.h' 

      FHernq =  (1.e0 + 1.e0/C_b)**2*(x/(x+1.e0))**2     

      Return  
      End
c------------------------------------ 
      REAL*4 Function aMmxh(x)  ! Contribution of hernquist bulge  to Halo vrms   
c------------------------------------  
      INCLUDE 'PMparameters.h' 
      INCLUDE 'Rodin2.h' 

c     aMmxh =  (1.e0 + 1.e0/Clarge)**2*(x/(x+1.e0))**2/ x**3/(1.+x)**3 
c      aMmxh =  (1.e0 + 1.e0/Clarge)**2/x/(x+1)**5    
      y = x*Rs/Rs_b 


      aMmxh=  (1.e0 + 1.e0/C_b)**2*(y/(y+1.e0))**2/x**3/(1.+x)**2       
      
      Return 
      End 
C--------------------------------------    
      REAL*4 Function aHdispr(x) ! hernquist rms velocities   
C---------------------------------------    
      INCLUDE 'PMparameters.h'  
      INCLUDE 'Rodin2.h'   
      
      aHdispr=(1.e0+ 1.e0/C_b)**2*(12.*x*(x+1.)**3*log((x+1.)/x) -
     *     x/(x+1.e0)*(25. + 52.*x +42.*x**2 + 12.*x**3))/12.e0    

      Return     
      End     
C---------------------------------------    
       SUBROUTINE  CM( )                
C---------------------------------------      
C     
      INCLUDE 'PMparameters.h'    
      INCLUDE 'Rodin2.h'      

      CHARACTER*100 Encab 
      Open (69,FILE=npath//'sats_params.dat',status='unknown')         

      Icount = 0  
      Xcmt = 0.  
      Ycmt = 0. 
      Zcmt = 0.  
      VXcmt = 0.   
      VYcmt = 0.    
      VZcmt = 0.   
      aMassc = 0.   

      READ(69,*)Encab  
      
 10   READ(69,*,err=20,end=20)Csat_i(Icount+1),aMsat_i(Icount+1), 
     *     Rv_i(Icount+1),Rs_i(Icount + 1),Vcm_i(Icount+1), !Read data   
     *     x_i(Icount+1),y_i(Icount+1),z_i(Icount+1), 
     *     vx_i(Icount+1),vy_i(Icount+1),vz_i(Icount+1)            
      
      
      x_i(Icount+1) = x_i(Icount+1)*1000.
      y_i(Icount+1) =  y_i(Icount+1)*1000.
      z_i(Icount+1) = z_i(Icount+1)*1000.

      
      Xcmt = Xcmt + aMsat_i(Icount + 1)*x_i(Icount + 1)   
      Ycmt = Ycmt + aMsat_i(Icount + 1)*y_i(Icount + 1)   
      Zcmt = Zcmt + aMsat_i(Icount + 1)*z_i(Icount + 1)   

      VXcmt = VXcmt + aMsat_i(Icount + 1)*vx_i(Icount + 1)  
      VYcmt = VYcmt + aMsat_i(Icount + 1)*vy_i(Icount + 1)   
      VZcmt = VZcmt + aMsat_i(Icount + 1)*vz_i(Icount + 1)   


      aMassc  = aMassc + aMsat_i(Icount + 1)
      Icount = Icount + 1      

      goto 10  

 20     Close(69)  
      
        Xcmt = Xcmt/ aMassc   
        Ycmt = Ycmt/aMassc   
        Zcmt = Zcmt/ aMassc   
        
        VXcmt = VXcmt/ aMassc 
        VYcmt = VYcmt/aMassc   
        VZcmt = VZcmt/ aMassc 

      Return   
      End   

C---------------------------------------  
      SUBROUTINE Rotx(ifact,alphad)                !Rot around x 
C---------------------------------------  
C     
      INCLUDE 'PMparameters.h'  
      INCLUDE 'Rodin2.h'  

      ilower = 1 + lspecies(Nspecies)*ifact
      iupper = lspecies(Nspecies)*(1+ifact) 


      alphad = alphad*pi/180.d0

      Do i=ilower,iupper 
         yc(i) = yc(i)*cos(alphad) - zc(i)*sin(alphad)    
         zc(i) = zc(i)*cos(alphad) + yc(i)*sin(alphad)     
      Enddo
      Return 
      End  
      
C---------------------------------------
      SUBROUTINE Roty(ifact,betad)     !Rot around y 
C---------------------------------------
C     
      INCLUDE 'PMparameters.h' 
      INCLUDE 'Rodin2.h' 

      ilower = 1 + lspecies(Nspecies)*ifact 
      iupper = lspecies(Nspecies)*(1+ifact) 


      betad= betad*pi/180.d0


      Do i=ilower,iupper 
         xc(i) = xc(i)*cos(betad) - zc(i)*sin(betad)  
         zc(i) = zc(i)*cos(betad) + xc(i)*sin(betad)    
      Enddo   

      Return  
      End  


C--------------------------------------- 
      real*4 Function cHdsigmaz(yy)     !Halo contribution to disk_sigam_z in the  the disk plane  
C--------------------------------------- 
C     
      INCLUDE 'PMparameters.h' 
      INCLUDE 'Rodin2.h'  
       common/sph/RRc 
       Common/zprof/zmin 

       External AINTG,aHInt     
       REAL*4 AINTG,aHInt  

       RRc = yy*Rd         
       w = RRc/Rs   
	zmin= zmin*Rd/Rs !sigma_z in the plane    
	zmax = 3.*zd*Rd/Rs   
        cHdsigmaz =  Aintg(aHInt,zmin,zmax)      
c	write(*,*)'test zmin',zmin,zmax,cHdsigmaz   


       return  
       end   
C---------------------------------------  
      real*4 Function aHInt(zz)     !Integrand Halo contribution  
C---------------------------------------  
C     
      INCLUDE 'PMparameters.h' 
      INCLUDE 'Rodin2.h' 
      common/sph/RRc 
      Common/zprof/zmin 


       w = RRc/Rs  

       rsph = Sqrt(w**2 + zz**2)  !w=Rc/Rs, zz = z/Rs    

       apass =  rsph**3*(1. + rsph)    
       aHInt = zz*(rsph-log(1. + rsph)*(1.+ rsph))/ 
     * max(apass,1.e-17)*4./(exp(zz) + exp(-zz))**2        

       return	 
       end     

C--------------------------------------- 
      real*4 Function cBdsigmaz(yy)     !Bulge contribution to sigma_z in the disk plane 
C--------------------------------------- 
C     
      INCLUDE 'PMparameters.h' 
      INCLUDE 'Rodin2.h'  
       common/sph/RRc 
       Common/zprof/zmin 

       External AINTG,aBint    
       REAL*4 AINTG,aBint 

       RRc = yy*Rd        
       w = RRc/Rs_b  
	zmin= zmin*Rd/Rs_b !sigma_z in the plane     
	zmax = 3.*zd*Rd/Rs_b  
	
        cBdsigmaz =  Aintg(aBInt,zmin,zmax)      
       return 
       end 
C---------------------------------------   
      real*4 Function aBInt(zz)     !Integrand Bulge  contribution  
C---------------------------------------  
C     
      INCLUDE 'PMparameters.h'  
      INCLUDE 'Rodin2.h'  

       common/sph/RRc 
       Common/zprof/zmin  
 

       w = RRc/Rs_b   

       rsph = Sqrt(w**2 + zz**2)  !r/Rs, zz = z/Rs_b  

       apass =  rsph*(1. + rsph)**2   
       aBInt = zz/ 
     * max(apass,1.e-17)*4./(exp(zz) + exp(-zz))**2       

       return
       end    


C--------------------------------------- 
      real*4 Function cHsf(yy)     !Halo contribution to Sfdensity along the disk plane 
C--------------------------------------- 
C     
      INCLUDE 'PMparameters.h' 
      INCLUDE 'Rodin2.h'  
       common/sph/RRc  
       External AINTG,aHd    
       REAL*4 AINTG,aHd  

       RRc = yy*Rd       
       w = RRc/Rs  
	zmin= 0.0  
	zmax = 3.*zd*Rd/Rs     

c        cHsf =  Aintg(aHD,zmin,zmax)      
        cHsf = RMvirMd*(Rd/Rs)**2/Fcon/ 
     *  (yy*Rd/Rs)/((yy*Rd/Rs)+1.)**2*  
     * 	coef(1)*zd*Rd/Rs  		    


       return	
       end  
C---------------------------------------   
      real*4 Function aHD(zz)     !Halo density 
C---------------------------------------  
C     
      INCLUDE 'PMparameters.h'   
      INCLUDE 'Rodin2.h' 
      common/sph/RRc 

       
       w = RRc/Rs  

       rsph = Sqrt(w**2 + zz**2)  !r/Rs, zz = z/Rs 

       apass =  rsph*(1. + rsph)**2   
       aHD = 1.0/max(apass,1.e-17)     

       return	 
       end   

C--------------------------------------- 
      real*4 Function cBsf(yy)     !Bulge contribution to Sfdensity 
C--------------------------------------- 
C     
      INCLUDE 'PMparameters.h' 
      INCLUDE 'Rodin2.h' 
       common/sph/RRc   
       External AINTG,aBd    
       REAL*4 AINTG 

       RRc = yy*Rd       
       w = RRc/Rs_b 
	zmin= 0.0  
	zmax = 3.*zd*Rd/Rs_b   

c       cBsf =  Aintg(aBD,zmin,zmax)     
        cBsf  = aMassb/aMassd/Fhb*(Rd/Rs_b)**2/          !Bulge SF  
     *  (yy*Rd/Rs_b)/((yy*Rd/Rs_b)+1.0)**3*coef(3)*zd*Rd/Rs_b  


       return	
       end  
C---------------------------------------  
      real*4 Function aBD(zz)     !Bulge Density  
C---------------------------------------    
C     
      INCLUDE 'PMparameters.h'  
      INCLUDE 'Rodin2.h' 
      common/sph/RRc  

       w = RRc/Rs_b    

       rsph = Sqrt(w**2 + (zz)**2)  

       apass =  rsph*(1. + rsph)**3  
       aBD = 1.0/max(apass,1.e-17)           

       return	 
       end       

  
C--------------------------------------- 
      real*4 Function DcHsf(yy)     !Halo contribution to Derivative of Sfdensity 
C--------------------------------------- 
C     
      INCLUDE 'PMparameters.h' 
      INCLUDE 'Rodin2.h' 
       common/sph/RRc 
       External AINTG,aHDd    
       REAL*4 AINTG 

       RRc = yy*Rd        
       w = RRc/Rs  
	zmin= 0.0  
	zmax = 3.*zd*Rd/Rs  

       DcHsf =  Aintg(aHDd,zmin,zmax)    


       return	
       end  
C---------------------------------------     
       real*4 Function aHDd(zz)     !Halo Density Cylindrical Derivative  
C---------------------------------------     
C     
      INCLUDE 'PMparameters.h'  
      INCLUDE 'Rodin2.h'  
      common/sph/RRc   

       w = RRc/Rs   

       rsph = Sqrt(w**2 + (zz)**2)  !zz = z/Rs  

       apass = (rsph*(1. + rsph))**3     
       aHDd = w*(1.0 + 3.*rsph)/max(apass,1.e-17)            

       return	 
       end      
C---------------------------------------   
      real*4 Function DcBsf(yy) !Bulge contribution to Derivative of Sfdensity  
C---------------------------------------  
C     
      INCLUDE 'PMparameters.h' 
      INCLUDE 'Rodin2.h' 
       common/sph/RRc 
       External AINTG,aBDd    
       REAL*4 AINTG 

       RRc = yy*Rd       
       w = RRc/Rs_b  
	zmin= 0.0  
	zmax = 3.*zd*Rd/Rs_b  

       DcBsf =  Aintg(aBDd,zmin,zmax)     


       return	
       end  
C---------------------------------------  
       real*4 Function aBDd(zz)     !Bulge Density Cylindrical Derivative    
C---------------------------------------    
C     
      INCLUDE 'PMparameters.h'  
      INCLUDE 'Rodin2.h' 
      common/sph/RRc  

       w = RRc/Rs_b  

       rsph = Sqrt(w**2 + zz**2)  !zz = z/Rs_b  

       aBDd =  (rsph*(1. + rsph))**3     
       aBDd = w*(1. + 4.*rsph)/max(aBDd,1.e-7)             

       return	 
       end     
 
c************************************************************************      
       real*4 function alininterp(u1,u2,r1,r2,x0)   
c************************************************************************ 
       include 'PMparameters.h'   
c       integer n   
       real x0   
c       call locate(u1,ngrid,x0,i)   
c       i=max(1,i)   
c       i=min(n-1,i)  
       alininterp=u1+((u2-u1)/(r2-r1))*(x0-r1)    
       return 
       end      

c************************************************************************       


C---------------------------------------   
      SUBROUTINE wtipsy( )                ! write tipsy ascii  
C---------------------------------------    
      INCLUDE 'PMparameters.h'    
      INCLUDE 'Rodin2.h'  

      Character*40 OUTFILE   
      real  tform(Nmaxpart),zmetals(Nmaxpart),pote(Nmaxpart),
     *     amg(Nmaxpart),eps(Nmaxpart),epss(Nmaxpart) 


      rscale = 2.87e4
      vscale = 691.0
      amscale = 3.17e15  

      Ngas = 0 
      Nt = lspecies(nspecies) 
      If(Iconf.eq.0)then
         Nstar = 0
      else
         Nstar = Id 
      endif

      Ndim = 3 
      time = 0.0  
      Ndark = Nt-Ngas -Nstar 

      pote_def = -1.e3 
      tform_def = -1.0 
      zmetals_def = -100.0 
c      eps_def = 0.056/rscale    
c       eps_def = 0.09/rscale ! khtk
        eps_def = 0.07/rscale !pseud  
       


      Do i=1,Nt  
       tform(i) =  tform_def  
       zmetals(i) = zmetals_def  
       pote(i) = pote_def    
       eps(i) = eps_def  
       epss(i) = eps_def 

       k=1  
 13    if(i.le.lspecies(k))then  
          amg(i) = pmass(1)*wspecies(k)/wspecies(1)/hubble/amscale   
       else   
          k = k+1  
          if(k.gt.nspecies)then  
             write(*,*)'error k > nspecies'  
             stop  
          endif   
          go to 13  
       endif  
        
       xc(i) = xc(i)*Rs/hubble/rscale 
       yc(i) = yc(i)*Rs/hubble/rscale       
       zc(i) = zc(i)*Rs/hubble/rscale 
          
       vxc(i) = vxc(i)*208.e-5*sqrt(aMassh/Fcon/Rs)/vscale 
       vyc(i) = vyc(i)*208.e-5*sqrt(aMassh/Fcon/Rs)/vscale        
       vzc(i) = vzc(i)*208.e-5*sqrt(aMassh/Fcon/Rs)/vscale 
          

      enddo     


      OUTFILE = 'model.tipsy.dat'  
      
      OPEN(10,FILE=npath//OUTFILE,STATUS='UNKNOWN')    
      
      write(10,*)Nt,Ngas,Nstar          
      write(10,*)Ndim              
      write(10,*)time   

      write(*,*)'Nt Ngas Nstar ',Nt,Ngas,Nstar 
      
      
         

      Do i=Nt,1,-1 
        
         if(i.le.Nstar)then  
            write(10,*)amg(i)  
         else 
            write(10,*)amg(i)  
         endif
         write(*,*)i 
      enddo 
       write(*,*)'-----------------------------'
       write(*,*)'tipsy mass ok'
      Do i=Nt,1,-1 
        if(i.le.Nstar)then 
         write(10,*)xc(i)  
         else  
            write(10,*)xc(i)   
         endif  
      enddo 
      write(*,*)'tipsy x ok'
      Do i=Nt,1,-1 
         if(i.le.Nstar)then 
         write(10,*)yc(i)  
         else  
            write(10,*)yc(i)   
         endif  
      enddo 
         
      Do i=Nt,1,-1  
         if(i.le.Nstar)then 
            write(10,*)zc(i)  
         else  
            write(10,*)zc(i)   
         endif  
      enddo 
         
    
      Do i=Nt,1,-1 
         if(i.le.Nstar)then 
            write(10,*)vxc(i)  
         else  
            write(10,*)vxc(i)   
         endif  
      enddo 

      Do i=Nt,1,-1 
         if(i.le.Nstar)then 
            write(10,*)vyc(i)  
         else  
            write(10,*)vyc(i)   
         endif  
      enddo 


       Do i=Nt,1,-1 
         if(i.le.Nstar)then 
            write(10,*)vzc(i)  
         else  
            write(10,*)vzc(i)   
         endif  
      enddo  
     

     
        


         
        

        

      if(Ndark.gt.0)then 
         Do i=1,Ndark  
            write(10,*)eps(i+ Ngas)  
         enddo 
         
      endif


      if(Nstar.gt.0)then 
         Do i=1,Nstar 
            write(10,*)epss(i)   
         enddo  
      endif 
 
c      if(Ngas.gt.0)then 
c         write(10,*)(densg(i),i=1,Ngas)  
c         write(10,*)(tempg(i),i=1,Ngas)  
c         write(10,*)(ahsmooth(i),i=1,Ngas)   
c         write(10,*)(zmetalg(i),i=1,Ngas)    
c      endif 
      if(Nstar.gt.0)then 
         Do i=1,Nstar 
            write(10,*)zmetals(i) 
         enddo 
         Do i=1,Nstar 
            write(10,*)tform(i) 
         enddo 
            
            
      endif 
      
      Do i=Nt,1,-1  
         write(10,*)pote(i)  
      enddo 
            

      return    
      end  

