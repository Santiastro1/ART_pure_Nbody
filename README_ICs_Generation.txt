This is a short description of the process to follow in order to obtain a simulation with multiple interacting galactic systems.

Please chance it at your convenience and add the changes you make (and its purpose) at the end of this document. If you have any question regarding the process of making new ART ICs, do not hesitate on contacting S. Roca-Fabrega (sroca01@ucm.es).

---------------------

The first step you need to do is to choose a snapshot of the simulation of your main galactic system. You can also generate new ICs for your main galactic system and run it up to 500Myr (recommended: a=0.6500 to 0.6750), in order to relax th ICs.

-------------  Choosing your main galactic system ----------------------------

---- Generation of new ICs for your main galactic system ---------------------
---- Programs you need: Rodin2.f + Rodin2.h + PMparameters.h + Conf.R.dat ----
------------------------------------------------------------------------------
ICs have to be generated with the RODIN software (see the RODIN folder):

1- You first need to create your own Conf.R.dat file. In this file you will choose the main properties of your galactic system. The system can be a simple NFW dark matter spheroid, a NFW DM spheroid + stellar bulge (Hernquist), a NFW DM spheroid + stellar bulge (Hernquist) + stellar disk (Myamoto Nagai), or a combination of these three components. 

!!! all units need to be '/h' where h=0.7 !!!

	a) Stellar disk: 
		-Nd: Number of particles (will set your mass resolution)
		-aMassd: Disk mass (combined with the Nd will set your mass resolution) in Msun/h
		-Rd: Disk exponential scale length (in kpc/h)
		-Zdkp: Vertical scale length (kpc/h)
		-Q: Toomre parameter (~1 Disk unstable to gravitational perturbations, >>1 disk is stable)
		-Rtkp: Disk truncation radius (kpc/h)
		-Rrefkp: Disk normalization radius (kpc/h)

	b) Stellar bulge:
		- aMassb: Bulge mass (Msun/h) - mass/particle will be the same as the on for disk particles (all star particles will have the same mass).
		- C_b: Bulge Hernquist concentration 
		- Rs_b: Bulge .25 mass radius (kpc/h)
		- Xout_b: Bulge truncation radius (kpc/h)

	c) Dark matter halo (Rs needs to be computed by Rvir=c*Rs):
		- aMassh: Total halo mass (Msun/h)
		- Con: Navarro Frenk and White concentration parameter
		- Xout: Halo truncation radius in Rs units (usually Xout = Con, so Rvir = Truncation radius)
		-**nspecies**!!!: Very important parameter. This parameter sets the number of DM species in the simulation. Each new specie has double mass than the previous one. This approach allows to create models with high mass/space resolution in the disk, with using a reasonable amount of CPU time. The first specie will have the same mass as "stellar disk/bulge" particles. It is highly recommended to use 5-6 species in a MW-mass galaxy. BE CAREFUL!!!! The more massive species should never be placed closer to the disk than 5Rs!!! The second DM mass specie should be placed at least at 2Rs.
		-Xb(i): Outer boundary for each mass specie (inner boundary is set by the previous specie outer boundary. WARNING! The last specie outer boundary needs to coincide with NFW virial radius.

#######  It is mandatory to test, after getting the final result of your simulation, that the disk contamination by less massive dark matter species is lower tha 1/1000. The most massive DM particles should have a much lower contamination rate!!! Each massive DM particle will have a satellite-like impact to the disk of the galaxy #####  THE LACK OF CONTAMINATION BY MASSIVE DM PARTICLES HAVE BEEN WELL TESTED IN PREVIOUS SIMULATIONS (when placing massive dark matter species > 5Rs)

	d) General parameters: 
		-Box: Boxsize for the simulation (Mpc/h) - periodical boundary conditions are applied, so it is recommended to use a 5*Rvir box ~ 1Mpc/h for a MW mass galaxy.
		-Iconf: It defines the system you are willing to generate
      			0:  NFW Halos 
      			1:  Exponential Disk     
      			2:  Exponential Disk + Hernquist Bulge  
      			3: Exponential Disk within a NFW halo  
      			4: Exponential Disk + Bulge inside a NFW halo  
		-aexpn: Initial scale factor, recommended 0.6000
		-astep: minimum timestep, recommended 0.0005 (time resolution)
		-overdens: critical density of the universe (default 340.)
		-ocurv: Universe curvature parameter (=0 in non-cosmological)
		-nseed: random seed (do not touch it! this will make the result difficult to reproduce)
		-hubble: hubble constant/100 (km/s/Mpc), default=0.7
		-Om0: Omega Matter (default = 0.3, unused in non-cosmological)
		-Oml0: Omega Lambda (default = 0.7, usused in non-cosmological)

	e) Unused parameters (not well tested yet):
		-Isat: Recommended =0
		-IMM: Recommended =0
		-ISSat: Recommended =0
		-*_sat, arot, alpha1, beta1, alpha2, beta2: Recommended = 0


2- Generating the initial conditions:

	a) Check the PMparameters.h file (Nmaxpart, Nrow and Ngrid, default values should be ok)
	b) Check the Rodin2.F file (initialization density-velocities Widrow and Dubinski approach):
		-Nrad need to be changed to higher values when a large number of particles is generated, this variable is related with the generation of a distribution of particles / radii, following density profiles and avoiding Poisson noise.

	c)Make Rodin2
	d)Execute Rodin2 (./Rodin2) - you can run it also in OpenMP
	e)Revise all the diagnostic outputs generated (.dat) - density and velocity profiles, ... etc.
	f)Change filenames (X.XXXX = aexpn in Conf.R.dat):
		PMcrs0.DAT --> PMcrs0aX.XXXX.DAT
		PMcrd.DAT --> PMcrdaX.XXXX.DAT
		pt.dat --> ptaX.XXXX.dat


---- Selection of a snapshot of an existing simulation ------

The only thing you need to check if you plan to use an already relaxed simulation is the number of dark matter species and the minimum particle mass. All galactic systems you add to the main simulation need to have the same mass / specie.

The files you need are the PMcrs0aX.XXXX.DAT, PMcrdaX.XXXX.DAT, pta0.XXXX.dat. 

------------------------------------------------------------------------------
------------------------------------------------------------------------------

You now have the PMcrs0aX.XXXX.DAT, PMcrdaX.XXXX.DAT, pta0.XXXX.dat files, these are the main galaxy ICs.

------------------------------------------------------------------------------
------------------------------------------------------------------------------

-------  Generating the satellite ICs ----------------------------------------
---- Programs you need: Rodin2.f + Rodin2.h + PMparameters.h + Conf.R.dat ----

You need to follow the same instructions for the main galaxy ICs  generation (see instructons above), BUT, you need to ensure the minimum particle mass is the same in the satellite than in the main galactic system.


KEEP THE ORIGINAL IC FILENAMES FOR THE SATELLITE:

PMcrs0.DAT, PMcrd.DAT, pt.dat


The minimum particle mass is set in Conf.R.dat. 
	In all cases, also when you are generated a NFW dark matter only, the minimum mass is set by aMassd/Nd.

------------------------------------------------------------------------------
------------------------------------------------------------------------------

------- Generating the Main Galaxy + satellite ICs ---------------------------
------ Programs you need: Add_to_ART.f + PMparART.h --------------------------


Now you have the IC files of the main galactic system (PMcrs0aX.XXXX.DAT, PMcrdaX.XXXX.DAT, pta0.XXXX.dat), and the ones of the satellite (PMcrs0.DAT, PMcrd.DAT, pt.dat). The following steps will merge the two set of files to a single one with information of both systems (PMcrs0_ini.DAT, PMcrd_ini.DAT, pt_ini.dat).

1- Check that parameters in PMparART.h are the same as the ones in the PMparameters.h and Rodin2.h (RODIN parameter files you used in the previous steps).
2- Choose the satellite position and velocities. (x,y,z)=(0,0,0) is the location of the main galactic system center. Velocities will also be in the Main galaxy (non-rotating) reference frame. -- you need to add satellite positions and velocities in PMparART.h (xsat,ysat,zsat in Kpc, vxsat, vysat, vzsat in Km/s/)
3- Compile: gfortran -O2 Add_to_ART.f -o Add2ART.exe
4- Move Add2ART.exe to the same folder as the Main galaxy and satellite IC files.
5- ./Add2ART.exe

------------------------------------------------------------------------------
----------  Analysis ---------------------------------------------------------
------------------------------------------------------------------------------
To make a first visual inspection of the generated ICs (and also of the final run snapshots), it is recommended to use the Tipsy software: https://github.com/N-BodyShop/tipsy
Tipsy is easy to install and has a README file with a full description of how to read and visualize data.

Tipsy authomatically generates a binary tipsy file from the provided ascii tipsy file, this binary file can be easily read by "yt": https://yt-project.org/doc/examining/loading_data.html

The first commands to use when running tipsy are:

openascii NAME_ASCII_FILE
readascii NAME_OUTPUT_BINARY(different from NAME_ASCII_FILE, if not, it will be overwriten)
loadb 0
"visualization commands", see: http://faculty.washington.edu/trq/hpcc/tipsy/man/

-------Generating the tipsy ascii file from the crude snapshots---------------
------------------PM_to_ascii_tipsy folder------------------------------------


We developed a "translator" from ART to ascii_tipsy format:
PM_to_ascii_tipsy

You just need to add the exact number of stellar particles in your simulation in the "PMparART.h" parameter file (diskpart variable). You also need to add the box size in the "PM_to_ASCIIoTIPSY_ART.f" (line 29).

Compile it by: gfortran -O2 PM_to_ASCIIoTIPSY_ART.f -o PM_to_ASCIIoTIPSY_ART.exe

Run it: ./PM_to_ASCIIoTIPSY_ART.exe

You need to have the snapshots (PMcrdaX.XXXX.DAT and PMcrsaX.XXXX.DAT) in the same folder as the .exe
You also can convert to tipsy the ICs snapshots (PMcrs0_ini.DAT, PMcrd_ini.DAT) by slightly modifying the PM_to_ASCIIoTIPSY_ART.f program (lines 56 to 58, and, 89 to 92)


------------------------------------------------------------------------------
------------------------------------------------------------------------------
Now you are ready to start the simulation using the pure N-body (no hydro) ART code.
------------------------------------------------------------------------------
------------------------------------------------------------------------------
---------------Folder ART ----------------------------------------------------

The N-body code you will use is the ART code (Adaptative Refinement Tree, By Klypin and Kravtsov).

-You just need to take a look on the "a_setup.h" file and change the following parameters according to the ones you used in the ICs generation:
	1- nrow
	2- ngrid
	3- nspecies
	4- nbyteword (1 or 4, depending on the length of the binary file).
	5- np (not critical)

-You need to decide the max spatial resolution you want to archeive:
	-MaxLevel (max spatial resolution = box[kpc]/ngrid/2**MaxLevel) (*2, if you want to compare with softening parameter of SPH simulation)

------------------------------------------------------------------------------
Compile the code by "make" (change the compiler options in the makefile, at your convenience).
THERE IS A HIGH CHANCE YOU WILL NEED TO CHANGE THE COMPILATION OPTIONS ACCORDING TO THE COMPILERs VERSION YOU ARE USING, if you have any problem with the compilation please do not hesitate on contacting sroca01@ucm.es

Now that you got the "art" executable you are almost ready to run your simulation.

Create a new folder "RUN", and copy the "art" file and  the IC files "PMcrd_ini.DAT, PMcrs0_ini.DAT and pt_ini.dat", changing their name to "PMcrd.DAT, PMcrs0.DAT and pt.dat".

You can initiate the simulation by launching the "art" job (./art , mpirun art, ...).

Remmember this is an OPENMP job, you need to set the OPENMP number of threads by:


set environment OMP_NUM_THREADS= number of threads (csh)

or

export OMP_NUM_THREADS= number of threads (bash)

------------------------------------------------------------------------------
--------------- Analysis -----------------------------------------------------

You can use the same procedure as for the ICs (conversor ART --> Tipsy + Tipsy + yt).

------------------------------------------------------------------------------
If you have any problem, please do not hesitate on contacting sroca01@ucm.es





