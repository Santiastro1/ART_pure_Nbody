opt = -O3 -i-dynamic -mcmodel=medium -convert big_endian

lib =      
fort = ifort                    

Add2ART :PMparART.h Add_to_ART.o; \
	$(fort) $(opt) Add_to_ART.o -o Add2ART

Add_to_ART.o ::PMparART.h Add_to_ART.f
	$(fort) $(opt) -c Add_to_ART.f

clean :
	rm *.o Add2ART



