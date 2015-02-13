import numpy as np

# ====================== #
# Conversion conventions #
# ====================== #

# These conversions are compatible with numpy's fftshift procedure:
# Center = 0 at the middle of array (index = size / 2)
# Position = 0 at the beginning of array (index = 0)
center_to_pos  = lambda center,size: int(size)/2 + center
pos_to_center  = lambda pos,size: pos - int(size)/2

# Position conversion for a downsampled / upsampled array:
downsample_pos = lambda pos,size,binsize: (pos-(binsize-1)/2.)*(size/(1.*binsize)-1)/(1.*(size-binsize))
upsample_pos   = lambda pos,size,binsize: pos*(size*binsize-binsize)/(1.*(size-1))+(binsize-1)/2.
