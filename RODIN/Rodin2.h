C The general system of units used is  Rs=Mh/F(c)=G=1(Halo  System), 
C but physical quantities are asked for the input   
C in physical units(M-> Solarmass/h  R->kpc/h) 
C  afterwards they are scaled  
C in the corresponding local system of units for each component   
C and at the end all the quantities will be   
C expressed back in physical  units.  
C The Q value at the solar neighborhood,Md,Rd and Zd were  
C selected according to the table 1-1 from Binney&Tremaine book, 
C see also the discussion at pages 363-364.  
C The extension of the disk is fixed by 97% of the mass. 
        
          
c          IMPLICIT REAL*4(A-H,O-Z)    

C     Commom Variables  and Parameters
      INTEGER Nseed     
      INTEGER NWANT,Nmo,Nd ,Nbin   
c       Real*8 aMlimit,pMd,pMdi,aMleft,dMbig,pMbig         
      character*2 npath

      PARAMETER (npath ='./')
      PARAMETER ( G           = 6.67d-8                    )!G in cgs units  
      PARAMETER (cKpcm        = 3.018d21                    ) !cm/kpc               
      PARAMETER (cKmcm        = 1.d-5                        ) !km/cm  
      PARAMETER (SM           =2.d33                        ) !solar mass in grams  
      PARAMETER ( pi          = 3.1415926535d+0    )       
      PARAMETER ( twopi       = 2.d0*3.14159265d+0)      
      PARAMETER (Nmo          = Nmaxpart  ) ! Max. Number of particles    
      PARAMETER (Nspmax       = 10)  !Max number of species 
       Common /pass/ Rslarge ,aMlarge, aMass_vir,Clarge  

C------------------------------------------------------------------------  
C     Disk Quantities           
C     Rref=8.5 kpc  is  assumed to be the solar radius where     
C     measurements are more accurate.      
C     The units are translated to the local system of units 
C     rd=Md=G=1.  
C     The assumed values are: 
C     For the Dwarf model we assumed Rd = 1.4 kpc, Md=8.4d8 solar masses. 
C     Mh= 1.d10   

        Common /disk/Nd,aMassd,Rd,Zdkp,Q,Rtkp,Rrefkp  
        Common /disk2/zcoef,zd,Rref,aMd,RT,a,delta 
 
        Common /bulge/aMassb,C_b,Xout_b,Nbulge,R0_b,V0_b,Rs_b,FHb       
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

            Common /satellites/aMass_sat,Con_sat ,Xout_sat, 
     &    D0_sat,Vx0_sat,Vy0_sat,Vz0_sat,Nsat,Dsati,R0_sat, 
     &    V0_sat,Vx0_sati,Vy0_sati,Vz0_sati,D0x_sat,D0y_sat, 
     &    D0z_sat,Dsatix ,Dsatiy ,Dsatiz, IMM     
 
            parameter(Nssat=50)  

          Common/ssystem/aMsat_i(Nssat),x_i (Nssat),  
     &  y_i (Nssat),z_i (Nssat),vx_i (Nssat),vy_i (Nssat),vz_i (Nssat), 
     & Rv_i(Nssat),Rs_i(Nssat),Vcm_i(Nssat),Csat_i(Nssat),
     &  VXcmt,VYcmt, VZcmt,Xcmt,Ycmt,Zcmt 

             Common /stat/Ibin,Ncount,Iorb, txcmh, tycmh, tzcmh,twh, 
     &   Iconf,coef(6),Isat,ISSat,Icount,Ialmacen,alpha1,beta1,          
     &   alpha2,beta2 

            Parameter(Irigid =0)     

            Dimension nspc(3)   
  

  
  


