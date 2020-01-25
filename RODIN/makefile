
#opt =  -byteswapio -O3  -mcmodel=medium -Mlarge_arrays
#opt =  -O3 -tpp6 -xK -w -Vaxlib
#opt = -O3 -xW -openmp -ipo -no-prec-div -i-dynamic -mcmodel=medium -convert big_endian
#opt = -O3 -openmp -i-dynamic -mcmodel=medium -convert big_endian
opt=-O2  -fopenmp -mcmodel=medium -fconvert=big-endian

lib =      
#fort = ifort                    
fort = gfortran
#fort = pgf77               
#fort = ifc
load = $(fort) $(opt)  $(olist) $(lib) -o Rodin1.x  

comp = $(name).o : $(name).f $(dep);          \
                    $(fort) $(opt) -c  $(name).f $(lib) 

$(this)/Rodin1 : $(olist) ; $(load)

statistics.b :  parameters.h PMparameters.h statistics.h statistics.b.o ; \
            $(fort) $(opt) statistics.b.o -o statistics.b

statistics.b.o : parameters.h PMparameters.h statistics.h statistics.b.f; \
            $(fort) $(opt) -c statistics.b.f

statistics.t :  parameters.h PMparameters.h statistics.h statistics.t.o ; \
            $(fort) $(opt) statistics.t.o -o statistics.t

statistics.t.o : parameters.h PMparameters.h statistics.h statistics.t.f; \
            $(fort) $(opt) -c statistics.t.f 


PM_to_ASCII :PMparameters.h PM_to_ASCII.o; \
            $(fort) $(opt)  PM_to_ASCII.o PM_to_ASCII

PM_to_ASCII.o :PMparameters.h PM_to_ASCII.f PMauxiliary.f; \
            $(fort) $(opt)  PM_to_ASCII.f PMauxiliary.f

Rodin2 : PMparameters.h  Rodin2.h Rodin2.o; \
          $(fort) $(opt) Rodin2.o -o Rodin2 

Rodin1 : PMparameters.h  Rodin1.h Rodin1.o; \
          $(fort) $(opt) Rodin1.o -o Rodin1 

Rodin1.o : Rodin1.f PMparameters.h  Rodin1.h; \
            $(fort) $(opt) -c Rodin1.f    

Rodin2.o : Rodin2.f PMparameters.h  Rodin2.h; \
            $(fort) $(opt) -c Rodin2.f 




