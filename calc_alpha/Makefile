main: main.o
	g++ -o $@ $< 

%.o: %.cpp params.h
	g++ -c $<

clean:
	rm main *.o