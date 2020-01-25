#FFLAGS  = -$(MODE) -model=F1 \
#         -nosave -excnvf -precexp=4 -noprecexpio \
#          -opt=ss -nopredicate \
#          -pvfunc=2 \
#          -procnum=8 -omp \
#          -verbose_loglist \
#          -uinline=2:ext
#LDFLAGS = -$(MODE) -parallel -pmfunc

#LDFLAGS = -q64  -qsmp=omp -O3 -qarch=pwr4 -qipa=level=2 -qfixed=72 -qtune=pwr4
#FFLAGS  = -q64  -qsmp=omp -O3 -qipa=level=2 -qarch=pwr4  -qfixed=72 -qtune=pwr4
#FC      = xlf90_r 
#LDFLAGS =  -O3 -g -mp  -byteswapio -fastsse -tp k8-64 -Knoieee -mcmodel=medium -Mlarge_arrays -Mbounds
#FFLAGS  =  -O3 -g -mp  -byteswapio -fastsse -tp k8-64 -Knoieee  -mcmodel=medium -Mlarge_arrays  -Mbounds
#LDFLAGS =  -O3 -g -mp  -byteswapio -fastsse -tp k8-64 -Knoieee -mcmodel=medium -Mlarge_arrays -Mbounds
#FFLAGS  =  -O3  -mp  -byteswapio -fastsse -tp k8-64 -Knoieee  -mcmodel=medium -Mlarge_arrays
#LDFLAGS =  -O3  -mp  -byteswapio -fastsse -tp k8-64 -Knoieee -mcmodel=medium -Mlarge_arrays
#FFLAGS = -fastsse -Munroll=n:4  -tp k8-64 -mcmodel=medium -byteswapio -mp -O3 
#LDFLAGS = -fastsse -Munroll=n:4 -tp k8-64 -mcmodel=medium -byteswapio -mp -O3
FFLAGS = -O3 -openmp -ipo -no-prec-div -i-dynamic -mcmodel=medium -convert big_endian 
LDFLAGS = -O3  -openmp -ipo -no-prec-div -i-dynamic -mcmodel=medium -convert big_endian

FC      = ifort
#FC      = /opt/appl/intel/fce/9.1.036/bin/ifort

MODE    =
LIB     = 


BIN     = art
OBJS    = ART_Advance.o ART_Auxiliary.o ART_Density1.o \
          ART_FFT.o ART_Init.o ART_IO.o \
          ART_Gravity.o ART_Main.o ART_Modify.o \
          ART_Particles.o ART_Tools.o
HEADERS = a_tree.h a_control.h a_setup.h a_constant.h a_numbers.h

$(BIN) : $(OBJS)
	        $(FC) $(LDFLAGS) $(OBJS) -o $(BIN) $(LIB)

$(OBJS): $(HEADERS)

ART_Tools.o: ART_Tools.f
	        $(FC) $(FFLAGS) -c $(LIB) $<


.f.o:
	        $(FC) $(FFLAGS) -c $(LIB) $<
#       $(FC) -c $(FFLAGS) -pmfunc $<

clean:
	        rm -f $(BIN) $(OBJS) 
