CC = g++-12
CFLAGS = -Wall -g -fopenmp

abp_2D_confine: abp_2D_confine.o
	$(CC) $(CFLAGS) -o abp_2D_confine.out abp_2D_confine.o  

abp_2D_confine.o: abp_2D_confine.cpp
	$(CC) -c abp_2D_confine.cpp

clean:
	rm *.o