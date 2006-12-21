all: libspimage.so
clean: 
	rm *.o && rm libspimage.so

CFLAGS += -O3 -Wall -I/home/filipe/c/fftw-3.1.1/build/include -I/home/filipe/c/hdf5-1.6.4/AMD/include -fPIC
LDFLAGS += -O3 -L/home/filipe/c/fftw-3.1.1/build/lib -L/home/filipe/c/hdf5-1.6.4/AMD/lib

libspimage.so: image_util.o fft.o 
	$(CC) -shared -o libspimage.so image_util.o fft.o 

image_util.o: ../image_util.c ../image_util.h
	$(CC) $(CFLAGS) -c -o image_util.o ../image_util.c

fft.o: ../fft.c ../fft.h
	$(CC) $(CFLAGS) -c -o fft.o ../fft.c

