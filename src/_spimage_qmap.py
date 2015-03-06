import numpy

def generate_absqmap(X,Y,pixel_size,detector_distance):
    qmap = generate_qmap(X,Y,pixel_size,detector_distance)
    qmap = numpy.sqrt(qmap[:,:,0]**2+qmap[:,:,1]**2+qmap[:,:,2]**2)
    return qmap

def generate_qmap(X,Y,pixel_size,detector_distance):
    qx = x_to_qx(X,pixel_size,detector_distance)
    qx = y_to_qy(Y,pixel_size,detector_distance)
    qz = xy_to_qz(X,Y,pixel_size,detector_distance)
    qmap = numpy.zeros(shape=(X.shape[0],Y.shape[1],3))
    qmap[:,:,0] = qz[:,:]
    qmap[:,:,1] = qy[:,:]
    qmap[:,:,2] = qx[:,:]
    qmap = qmap
    return qmap

_xy_to_qxy = lambda X_or_Y,pixel_size,detector_distance: 2*numpy.sin(numpy.arctan2(pixel_size*X_or_Y,detector_distance)/2.) / (pixel_size/detector_distance)
x_to_qx = lambda x,pixel_size,detector_distance: _xy_to_qxy(x,pixel_size,detector_distance)
y_to_qy = lambda y,pixel_size,detector_distance: _xy_to_qxy(y,pixel_size,detector_distance)
xy_to_qz = lambda x,y,pixel_size,detector_distance: 1 - numpy.cos( numpy.arctan2(pixel_size*numpy.sqrt(x**2+y**2),detector_distance) ) / (pixel_size/detector_distance)

    
