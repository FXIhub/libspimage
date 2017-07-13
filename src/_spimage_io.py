"""
A python module for reading/writing.

"""
import h5py



class CXILoader(object):
    """
    
    """
    def __init__(self, filename):
        self.filename = filename
        try:
            self.file = h5py.File(filename, "r")
        except IOError: print("Could not open file: %s" %filename)
    def close(self):
        self.file.close()
    def load_dataset(self, key):
        try:
            ds = self.file[key][:] 
        except: print("Could not find key: %s" %key)
        return ds
