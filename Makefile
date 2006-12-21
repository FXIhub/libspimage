all: libspimage.so
clean: 
	rm *.o && rm libspimage.so

CFLAGS += -O3 -Wall 
LDFLAGS += -O3

libspimage.so: image_util.o fft.o 
	$(CC) -shared -o libspimage.so image_util.o fft.o 
