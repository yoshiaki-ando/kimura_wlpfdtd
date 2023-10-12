OBJS = \
main.o \
calc_conductivity.o \
jz.o \
../assign_node_number.o \
../calc_const_pml.o \
../calc_wlp.o \
../compose_coef_matrix.o \
../compose_right_side.o \
../display_progress.o \
../expand_field.o \
../get_date_str.o \
../output_field.o \
../update_field.o \
../allocate_memory.o \
../csparse.o \
../get_time_sec.o \
../output_conductivity.o \
../make_profile.o \
../calc_d.o \
../get_coordinate.o \
../calc_const_param.o \
../calc_fourier.o \
../calc_coef_jz.o

OPTS = -O3 -Wall -std=c++1y -std=c++17
FDTD_INC = -I./
FFTW_LIB = -lfftw3 
FFTW_INC = 
GAUSS_INT_LIB = -L/home/ando/lib -lAndoLab -Wl,-R/home/ando/lib
GAUSS_INT_INC = -I/home/ando/include

HEADERS = wlp_fdtd2d.h 

.PHONY: all clean plot

all: main

main: $(OBJS)
	g++ -o $@ $(OBJS) $(OPTS) $(FFTW_LIB) $(GAUSS_INT_LIB)

%.o: %.cpp $(HEADERS)
	g++ -c $< -o $@ $(OPTS) $(FFTW_INC) $(GAUSS_INT_INC) $(FDTD_INC)

clean:
	rm -rf main *.o
	rm -rf ../*.o
