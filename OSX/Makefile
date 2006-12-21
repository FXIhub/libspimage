all: libspimage.dynlib
clean: 
	rm *.o && rm libspimage.so

CFLAGS += -fno-common -O3 -Wall -I/sw/include -fPIC -I/Users/filipe/home/c/fftw/OSX/include
LDFLAGS += -O3 -L/sw/lib -L/Users/filipe/home/c/fftw/OSX/include 
LOADLIBES += -lfftw3f -lfftw3f_threads -lhdf5 -ltiff -lpng

libspimage.dynlib: image_util.o fft.o 
	$(CC) -dynamiclib ${LDFLAGS} ${LOADLIBES} -o libspimage.dynlib image_util.o fft.o 

image_util.o: ../image_util.c ../image_util.h
	$(CC) $(CFLAGS) -c -o image_util.o ../image_util.c

fft.o: ../fft.c ../fft.h
	$(CC) $(CFLAGS) -c -o fft.o ../fft.c

