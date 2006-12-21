SOURCES += fft.c image_util.c
HEADERS += fft.h image_util.h image.h
TEMPLATE = lib
CONFIG -= qt

macx: INCLUDEPATH += /sw/include
macx: LIBS += -L/sw/lib


LIBS += -lpng -ltiff -lhdf5 -lfftw3f -lfftw3f_threads
