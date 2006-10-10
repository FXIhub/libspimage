all: libspimage.so
clean: 
	rm *.o && rm libspimage.so

CFLAGS += -g3 -Wall
LDFLAGS += -g3

libspimage.so: image_util.o fft.o 
	$(CC) -shared -o libspimage.so image_util.o fft.o /usr/lib/libfftw3f.so /usr/lib/libhdf5.so /usr/lib/libtiff.so /usr/lib/libpng.so
