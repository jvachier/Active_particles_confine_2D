CC = g++-12
CFLAGS = -Wall -g -fopenmp

abp_2D_confine: abp_2D_confine.o print_file.o reflective_boundary_conditions.o circular_reflective_boundary_conditions.o
	$(CC) $(CFLAGS) -o abp_2D_confine.out abp_2D_confine.o print_file.o reflective_boundary_conditions.o circular_reflective_boundary_conditions.o

abp_2D_confine.o: abp_2D_confine.cpp
	$(CC) -c abp_2D_confine.cpp

print_file.o: print_file.cpp
	$(CC) -c print_file.cpp

reflective_boundary_conditions.o: reflective_boundary_conditions.cpp
	$(CC) -c reflective_boundary_conditions.cpp

circular_reflective_boundary_conditions.o: circular_reflective_boundary_conditions.cpp
	$(CC) -c circular_reflective_boundary_conditions.cpp

clean:
	rm *.o