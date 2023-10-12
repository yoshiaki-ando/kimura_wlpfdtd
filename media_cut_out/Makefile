VER_SUFFIX = _18

OBJS = main.o func.o
LIBS = -L/usr/lib/x86_64-linux-gnu -lmatio -lz -lhdf5_cpp \
	-L$(HOME)/lib -Wl,-R$(HOME)/lib -lAndoLab$(VER_SUFFIX)

INC = -I$(HOME)/include
OPTS = $(INC) -O3 -Wall

main: $(OBJS)
	g++ -o $@ $(OBJS) $(LIBS)

%.o: %.cpp cut_out.h params.h
	g++ -c $<

clean:
	rm -f main *.o
