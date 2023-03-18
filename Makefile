CC = g++-12
CFLAGS = -Wall -g -fopenmp
 
# ****************************************************
# Targets needed to bring the executable up to date
 
abp_2D_confine: abp_2D_confine.o initialization.o
	$(CC) $(CFLAGS) -o abp_2D_confine abp_2D_confine.o initialization.o 
# The main.o target can be written more simply
 
abp_2D_confine.o: abp_2D_confine.cpp initialization.h
	$(CC) $(CFLAGS) -c abp_2D_confine.cpp

initialization.o: initialization.h

clean:
	rm *.o