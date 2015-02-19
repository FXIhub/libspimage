import numpy as np

# ====================== #
# Conversion conventions #
# ====================== #

# The middle of an array is defined as its center:
# x = (size-1)/2.
# Examples: - array of length 2 has it center at position 0.5.
#           - array of length 3 has it center at position 1.
# The absolute position (counted from the zeroth element) is determined from the center deviation (often referred to by cx and cy):
# x = cx - (size-1)/2.
center_to_pos = lambda center,size: (size-1)/2. + center
pos_to_center = lambda pos,size: pos - (size-1)/2.

# Position conversion for a downsampled / upsampled array:
downsample_pos = lambda pos,size,binsize: (pos-(binsize-1)/2.)*(size/(1.*binsize)-1)/(1.*(size-binsize))
upsample_pos   = lambda pos,size,binsize: pos*(size*binsize-binsize)/(1.*(size-1))+(binsize-1)/2.
