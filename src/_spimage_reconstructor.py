import numpy as np
import spimage

# Set up logger
import logging
logger = logging.getLogger('RECONSTRUCTOR')
logger.setLevel("WARNING")

class Reconstructor:
    def __init__(self, use_gpu=True):
        self._mask = None
        self._sp_amplitudes = None
        self._sp_amplitudes_dirty = True
        self._sp_initial_support = None
        self._sp_initial_support_sh = None
        self._sp_initial_support_dirty = True
        self._sp_phaser = None
        self._sp_phaser_dirty = True
        self._use_gpu = use_gpu
        self.clear()

    def clear(self):
        """
        Clears the instance reconstructor and deletes all configuration and data associated to it.
        """
        self._clear_amplitudes()
        self._clear_initial_support()
        self._clear_iterations()
        self._clear_support_algorithms()
        self._clear_phasing_algorithms()
        self._clear_phaser()
        self._reconstruction = None
        self._iteration = None
        self._log("Reconstructor initialized.","DEBUG")

    def _clear_amplitudes(self):
        for img in [self._sp_amplitudes]:
            if img is not None:
                spimage.sp_image_free(img)
        self._amplitudes = None
        self._sp_amplitudes = None
        self._amplitudes_dirty = True

    def _clear_iterations(self):
        self._number_of_iterations = None
        self._number_of_outputs_images = None
        self._number_of_outputs_scores = None
        self._iterations_dirty = True       

    def _clear_initial_support(self):
        for img in [self._sp_initial_support,self._sp_initial_support_sh]:
            if img is not None:
                spimage.sp_image_free(img)
        self._sp_initial_support = None
        self._sp_initial_support_sh = None
        self._initial_support_config = None
        self._initial_support_dirty = True

    def _clear_support_algorithms(self):
        self._support_algorithms = []
        self._support_algorithms_configs = []
        self._i_support_algorithms = None
        self._support_algorithms_dirty = True
        
    def _clear_phasing_algorithms(self):
        self._phasing_algorithms = []
        self._phasing_algorithms_configs = []
        self._i_phasing_algorithms = None
        self._phasing_algorithms_dirty = True

    def _clear_phaser(self):
        if self._sp_phaser is not None:
            spimage.sp_phaser_free(self._sp_phaser)          
        self._sp_phaser = None
        self._phaser_dirty = True

    def _get_scores(self):
        out = {}
        FcFo = np.zeros((1), dtype=np.float32)
        out["fourier_error"] = spimage.sp_phaser_efourier(self._sp_phaser,FcFo)
        out["FcFo"] = FcFo[0]
        out["real_error"] = spimage.sp_phaser_ereal(self._sp_phaser)
        out["support_size"] = spimage.sp_phaser_support_fraction(self._sp_phaser)
        return out
    
    def _get_scores_slow(self):
        out = {}
        [I,M] = self._get_curr_model(shifted=False)
        [fI,fM] = self._get_curr_fmodel(shifted=False)
        out["real_error"] = np.sqrt((abs(I[M==0])**2).sum()/( (abs(I[M==1])**2).sum() + np.finfo("float32").eps ))
        out["fourier_error"] = np.sqrt((abs(abs(fI[fM==1])-abs(np.fft.fftshift(self._sp_amplitudes.image)[fM==1]))**2).sum()/( (abs(np.fft.fftshift(self._sp_amplitudes.image)[fM==1])**2).sum() + np.finfo("float32").eps ))
        out["FcFo"] = abs(fI[fM==1]).sum()/(abs(self._sp_amplitudes.image[fM==1]).sum()+np.finfo("float32").eps)
        out["support_size"] = M.sum()/float(M.size)
        return out    

    def _log(self,s,mode="INFO"):
        if logger is not None:
            prefix = []
            ps = ""
            if len(prefix) > 0:
                ps += "("
                for p in prefix:
                    ps += p + " "
                ps = ps[:-1] + ")"
                ps += " "
            if mode == "INFO":
                f = logger.info
            elif mode == "DEBUG":
                f = logger.debug
            elif mode == "WARNING":
                f = logger.warning
            elif mode == "ERROR":
                f = logger.error
            f(" %s %s" % (ps,s))

    def __del__(self):
        self.clear() 

    def set_intensities(self,intensities,shifted=False):
        """
        Sets the intensity pattern that shall be phased. By default it is expected that the provided image is the non-shifted version of the diffraction pattern given as a 2D numpy array (diffraction pattern representation). If desired change to \"shifted=False\". If \"shifted=True\" pixel (0,0) is located in the corner of the physical diffraction pattern.
        """
        self._amplitudes_dirty = True
        self._initial_support_dirty = True
        A = np.array(intensities,dtype="float")
        if (A<0).sum() > 0:
            self._log("Intensities found that have negative values. These values are set to zero.")
            A[A<0] = 0.
        A = np.sqrt(A)
        if not shifted:
            self._amplitudes = np.fft.fftshift(A)
        else:
            self._amplitudes = A
        self._log("Intensities set.","DEBUG")

    def set_mask(self,mask,shifted=False):
        """
        The mask is a 2D boolean numpy array with the same shape as the provided intensities. Values that equal True indicate valid pixels and values that equal False indicate unknown intensity values at the respective pixel location. By default it is expected that the provided image is the non-shifted version of the diffraction pattern given as a 2D numpy array (diffraction pattern representation). If desired change to \"shifted=False\". If \"shifted=True\" pixel (0,0) is located in the corner of the physical diffraction pattern.
        """
        self._amplitudes_dirty = True
        self._initial_support_dirty = True
        M = np.array(mask,dtype="bool")
        if not shifted:
            self._mask = np.fft.fftshift(M)
        else:
            self._mask = M
        self._log("Mask set.","DEBUG")

    def set_number_of_iterations(self,number_of_iterations):
        """
        Set the total number of iterations.
        """
        self._number_of_iterations = number_of_iterations
        self._iterations_dirty = True
        self._log("Number of iterations set.","DEBUG")

    def set_number_of_outputs_images(self,number_of_outputs_images):
        """
        Specify the number of times that images (such as real space and fourier space and support images etc.) shall be written to the output. The output interval will be calculated from the number of iterations specified by set_number_of_iterations.
        """
        self._number_of_outputs_images = number_of_outputs_images
        self._iterations_dirty = True      
        self._log("Number of output images set.","DEBUG")

    def set_number_of_outputs_scores(self,number_of_outputs_scores):
        """
        Specify the number of times that scores (such as real space error and fourier space error and support size etc.) shall be written to the output. The output interval will be calculated from the number of iterations specified by set_number_of_iterations.
        """
        self._number_of_outputs_scores = number_of_outputs_scores
        self._iterations_dirty = True
        self._log("Number of output scores set.","DEBUG")

    def set_number_of_outputs(self,number_of_outputs):
        """
        Specify the number of times that scores (such as real space error and fourier space error and support size etc.) and images (such as real space error and fourier space error and support size etc.) shall be written to the output. The output interval will be calculated from the number of iterations specified by set_number_of_iterations. If you like to specify the number of outputs of images and scores separately call the functions set_number_of_outputs_scores and set_number_of_outputs_images.
        """
        self._number_of_outputs_images = number_of_outputs
        self._number_of_outputs_scores = number_of_outputs
        self._iterations_dirty = True
        self._log("Number of outputs set.","DEBUG")

    def set_initial_support(self,**kwargs):
        """
        Set the initial support by either giving:
        - the \"radius=initial_support_radius_in_pixels\" of the circular support mask
        - the \"area=initial_support_area_fraction_of_whole_image\" of the circular support mask
        - an explicit support mask \"support_mask=2D_boolean_numpy_array\", non-shifted and of the same shape as the intensity pattern
        """
        if not ("radius" in kwargs or "support_mask" in kwargs or "area" in kwargs):
            self._log("set_initial_support requires one of the following key word arguments: radius, area, support_mask","ERROR")
            return 
        self._support_dirty = True
        self._initial_support_config = {}
        if "radius" in kwargs:
            self._initial_support_config["radius"] = kwargs["radius"]
        elif "area" in kwargs:
            self._initial_support_config["area"] = kwargs["area"]
        else:
            self._initial_support_config["support_mask"] = kwargs["support_mask"]
        self._log("Initial support configuration set.","DEBUG")

    def set_support_algorithm(self, type, number_of_iterations=None, center_image=False, **kwargs):
        """
        Set one of the following support-update algorithms with the respective keyord arguments:
        - area: update_period, blur_init*, blur_final*, area_init**, area_final**
        - threshold: update_period, blur_init*, blur_final*, threshold_init, threshold_final
        - static (no keyword arguments for this algorithm)
        * a gaussian blur is applied with a standard deviation of the kernel with the given blur value in pixel
        ** area as a fraction of the entire image [0.0...1.0]
        NOTE: If you like to use a series of support-update algorithms please use instead the function append_support_algorithm.
        """
        self._clear_support_algorithms()
        self.append_support_algorithm(type, number_of_iterations=number_of_iterations, center_image=center_image, **kwargs)

    def append_support_algorithm(self, type, number_of_iterations, center_image=False, **kwargs):
        """
        Append one of the following support-update algorithms with the respective keyord arguments:
        - area: number_of_iterations, update_period, blur_init*, blur_final*, area_init**, area_final**
        - threshold: number_of_iterations, update_period, blur_init*, blur_final*, threshold_init, threshold_final
        - static: number_of_iterations
        * a gaussian blur is applied with a standard deviation of the kernel with the given blur value in pixel
        ** area as a fraction of the entire image [0.0...1.0]
        NOTE: If you like to use only a single support-update algorithm you might want to use instead the function set_support_algorithm.
        """
        alg_conf = {"type":type}
        alg_conf["center_image"] = center_image
        alg_conf["number_of_iterations"] = number_of_iterations
        if alg_conf["number_of_iterations"] is None and len(self._support_algorithms_configs) > 0:
            self._log("You can not have more than one support algorithm of unspecified number of iterations if the total number of iterations is set to None. Set the total number of iterations by calling set_number_of_iterations and try again.","ERROR")
            return
        necessary_kwargs = {"area":["update_period","blur_init","blur_final","area_init","area_final"],
                            "threshold":["update_period","blur_init","blur_final","threshold_init","threshold_final"],
                            "static":[]}
        if type not in necessary_kwargs:
            self._log("append_support_algorithm accepts algorithms of the following types: " + str(necessary_kwargs.keys()),"ERROR")
            return
        for k in necessary_kwargs[type]:
            if k not in kwargs:
                self._log("append_support_algorithm with the support algorithm \"area\" requires the following argument: " + k, "ERROR")
                return
            else:
                alg_conf[k] = kwargs[k]
        if alg_conf["number_of_iterations"] is not None:
            if len(self._support_algorithms_configs) == 0:
                self._i_support_algorithms = 0
            self._i_support_algorithms += alg_conf["number_of_iterations"]
        self._support_algorithms_configs.append(alg_conf)
        self._support_algorithms_dirty = True
        self._log("Support algorithm set/appended.","DEBUG")

    def set_phasing_algorithm(self, type, **kwargs):
        """
        Set one of the following phasing algorithms with the respective keyord arguments:
        - er (no keyword arguments for this algorithm)
        - hio: beta_init, beta_final
        - raar: beta_init, beta_final
        - diffmap: beta_init, beta_final, gamma1, gamma2
        If desired specify additional constraints in a list of strings:
        - enforce_positivity
        - enforce_real
        - enforce_centrosymmetry
        NOTE: If you like to use a series of phasing algorithms please use instead the function append_phasing_algorithm.
        """
        self._clear_phasing_algorithms()
        kwargs1 = dict(kwargs)
        if "number_of_iterations" not in kwargs:
            kwargs1["number_of_iterations"] = None
        self.append_phasing_algorithm(type, **kwargs1)

    def append_phasing_algorithm(self, type,**kwargs):
        """
        Append one of the following phasing algorithms with the respective keyord arguments:
        - er: number_of_iterations
        - hio: number_of_iterations, beta_init, beta_final
        - raar: number_of_itearations, beta_init, beta_final
        - diffmap: number_of_iterations, beta_init, beta_final, gamma1, gamma2
        NOTE: If you like to use only a single phasing algorithms you might want to use instead the function set_phasing_algorithm.
        """
        alg_conf = {"type":type}
        # check input
        if "number_of_iterations" not in kwargs:
            self._log("append_phasing_algorithm requires the keyword argument \'number_of_iterations\'.","ERROR")
            return
        alg_conf["number_of_iterations"] = kwargs["number_of_iterations"]
        if alg_conf["number_of_iterations"] is None and len(self._phasing_algorithms_configs) > 0:
            self._log("You can not have more than one phasing algorithm of unspecified number of iterations if the total number of iterations is set to None. Set the total number of iterations by calling set_number_of_iterations and try again.","ERROR")
            return
        alg_conf["constraints"] = kwargs.get("constraints","")
        # phasing algorithm
        necessary_kwargs = {"raar":["beta_init","beta_final"],
                            "hio":["beta_init","beta_final"],
                            "er":[],
                            "diffmap":["beta_init","beta_final","gamma1","gamma2"]}
        if type not in necessary_kwargs.keys():
            self._log("append_phasing_algorithm accepts algorithms of the following types: %s" % str(necessary_kwargs.keys()),"ERROR")
            return
        for k in necessary_kwargs[type]:
            if k not in kwargs:
                self._log("append_phasing_algorithm with the phasing algorithm " + type + " requires the following argument: " + k, "ERROR")
                return

        for k in kwargs.keys():
            alg_conf[k] = kwargs[k]

        if alg_conf["number_of_iterations"] is not None:
            if len(self._phasing_algorithms_configs) == 0:
                self._i_phasing_algorithms = 0
            self._i_phasing_algorithms += alg_conf["number_of_iterations"]
        self._phasing_algorithms_configs.append(alg_conf)
        self._phasing_algorithms_dirty = True
        self._log("Phasing algorithm set/appended.","DEBUG")

    def _prepare_reconstruction(self):
        self._ready = True
        if self._amplitudes is None:
            self._log("Reconstruction cannot start! You need to set the intensities.","ERROR")
            self._ready = False
        if self._mask is None:
            self._log("You did not set a mask, therefore initializing without masking any pixels of intensity pattern.","INFO")
            self._ready = True
        if self._initial_support_config is None:
            self._log("Reconstruction cannot start! You need to set the initial support.","ERROR")
            self._ready = False
        if self._support_algorithms_configs == []:
            self._log("Reconstruction cannot start! You need to set the support algorithm.","ERROR")
            self._ready = False           
        if self._phasing_algorithms_configs == []:
            self._log("Reconstruction cannot start! You need to set the phasing algorithm.","ERROR")
            self._ready = False           
        if self._number_of_outputs_images is None or self._number_of_outputs_scores is None:
            self._log("Connot prepare reconstruction. Number of outputs need to be set.","ERROR")
            self._ready = False           
        if self._i_support_algorithms is not None:
            if self._i_support_algorithms != self._number_of_iterations:
                self._log("Connot prepare reconstruction. Support algorithms initialised are not in line with the set number of iterations.","ERROR")
                self._ready = False
        if self._i_phasing_algorithms is not None:
            if self._i_phasing_algorithms != self._number_of_iterations:
                self._log("Connot prepare reconstruction. Phasing algorithms initialised are not in line with the set number of iterations.","ERROR")
                self._ready = False
        if not self._ready:
            return
        self._init_amplitudes()
        self._init_initial_support()
        self._init_iterations()
        self._init_support_algorithms()
        self._init_phasing_algorithms()
        self._init_phaser()

    def _init_amplitudes(self):
        if not self._amplitudes_dirty:
            self._log("Amplitudes already initialised.","DEBUG")
            return
        A = self._amplitudes
        if len(A.shape) == 2:
            self._Nx = A.shape[1]
            self._Ny = A.shape[0]
            if self._mask is not None:
                M = self._mask.copy()
            else:
                M = np.ones(shape=A.shape,dtype="bool")
            for img in [self._sp_amplitudes]:
                if img is not None:
                    spimage.sp_image_free(img)
            self._sp_amplitudes = spimage.sp_image_alloc(A.shape[0],A.shape[1],1)
            self._sp_amplitudes.image[:,:] = np.float32(A[:,:])
            self._sp_amplitudes.mask[:,:] = np.int32(M[:,:])
            self._sp_amplitudes.scaled = 1
            self._sp_amplitudes.phased = 0
            self._amplitudes_dirty = False
            self._log("Amplitudes initialised.","DEBUG")
        if len(A.shape) == 3:
            self._Nx = A.shape[2]
            self._Ny = A.shape[1]
            self._Nz = A.shape[0]
            if self._mask is not None:
                M = self._mask.copy()
            else:
                M = np.ones(shape=A.shape,dtype="bool")
            for img in [self._sp_amplitudes]:
                if img is not None:
                    spimage.sp_image_free(img)
            self._sp_amplitudes = spimage.sp_image_alloc(A.shape[0],A.shape[1],A.shape[2])
            self._sp_amplitudes.image[:,:,:] = np.float32(A[:,:,:])
            self._sp_amplitudes.mask[:,:,:] = np.int32(M[:,:,:])
            self._sp_amplitudes.scaled = 1
            self._sp_amplitudes.phased = 0
            self._amplitudes_dirty = False
            self._log("Amplitudes initialised.","DEBUG")
    # The array for the initial support is created just before reconstruction as we might not know the dimensions of the image earlier
    def _init_initial_support(self):
        if not self._initial_support_dirty:
            self._log("Initial support already initialised.","DEBUG")
            return
        for img in [self._sp_initial_support,self._sp_initial_support_sh]:
            if img is not None:
                spimage.sp_image_free(img)
        A = self._amplitudes
        if len(A.shape) == 2:
            self._sp_initial_support = spimage.sp_image_alloc(self._Ny,self._Nx,1)
            if "radius" in self._initial_support_config or "area" in self._initial_support_config:
                X,Y = np.meshgrid(np.arange(self._Nx),np.arange(self._Ny))
                X = X-(self._Nx-1)/2.
                Y = Y-(self._Ny-1)/2.
                R = np.sqrt(X**2 + Y**2)
                self._phaser_dirty = True
                if "radius" in self._initial_support_config:
                    r = self._initial_support_config["radius"]
                else:
                    r = np.sqrt( self._initial_support_config["area"] * self._Nx * self._Ny / np.pi )
                self._sp_initial_support.image[:] = np.float32(np.fft.fftshift(R) < r)
            else:
                self._phaser_dirty = True
                S = self._initial_support_config["support_mask"]
                self._sp_initial_support.image[:] = np.float32(S)
                self._initial_support_dirty = False
                self._log("Initial support initialised.","DEBUG")
        if len(A.shape) == 3:
            self._sp_initial_support = spimage.sp_image_alloc(self._Nz, self._Ny,self._Nx)
            if "radius" in self._initial_support_config or "area" in self._initial_support_config:
                X,Y,Z = np.meshgrid(np.arange(self._Nx),np.arange(self._Ny),np.arange(self._Nz))
                X = X-(self._Nx-1)/2.
                Y = Y-(self._Ny-1)/2.
                Z = Z-(self._Nz-1)/2.
                R = np.sqrt(X**2 + Y**2 + Z**2)
                self._phaser_dirty = True
                if "radius" in self._initial_support_config:
                    r = self._initial_support_config["radius"]
                else:
                    r = np.sqrt( self._initial_support_config["area"] * self._Nx * self._Ny *  self._Nz / np.pi )
                self._sp_initial_support.image[:] = np.float32(np.fft.fftshift(R) < r)
            else:
                self._phaser_dirty = True
                S = self._initial_support_config["support_mask"]
                self._sp_initial_support.image[:] = np.float32(S)
                self._initial_support_dirty = False
                self._log("Initial support initialised.","DEBUG")

    def _init_iterations(self):
        if not self._iterations_dirty:
            self._log("Iterations already initialised.","DEBUG")
            return
        self._out_iterations_images = np.int64((np.linspace(0,self._number_of_iterations,self._number_of_outputs_images)).round())
        self._out_iterations_scores = np.int64((np.linspace(0,self._number_of_iterations,self._number_of_outputs_scores)).round())
        self._iterations_dirty = False
        self._support_algorithms_dirty = True
        self._phasing_algorithms_dirty = True
        self._log("Iterations initialised.","DEBUG")

    def _init_support_algorithms(self):
        if not self._support_algorithms_dirty:
            self._log("Support algorithms already initialised.","DEBUG")
            return       
        i = 0
        for alg_conf in self._support_algorithms_configs:
            alg = dict(alg_conf)
            if alg_conf["number_of_iterations"] is None:
                if len(self._support_algorithms_configs) == 1:
                    alg["number_of_iterations"] = self._number_of_iterations
                else:
                    self._log("Number of iterations can not be None if many support algorithms are specified. Please report this error.","ERROR")
                    return
            if alg["type"] == "area" or alg["type"] == "threshold":            
                blur_radius = spimage.sp_smap_alloc(2)
                spimage.sp_smap_insert(blur_radius, i, alg["blur_init"])
                spimage.sp_smap_insert(blur_radius, i + alg["number_of_iterations"], alg["blur_final"])
                support_area = spimage.sp_smap_alloc(2)
                spimage.sp_smap_insert(support_area, i,alg["area_init"])
                spimage.sp_smap_insert(support_area, i + alg["number_of_iterations"],alg["area_final"])            
            if alg["type"] == "area":
                alg["spimage_support_array"] = spimage.sp_support_array_init(spimage.sp_support_area_alloc(blur_radius, support_area),alg["update_period"])
            elif alg["type"] == "threshold":
                alg["spimage_support_array"] = spimage.sp_support_array_init(spimage.sp_support_threshold_alloc(self._blur_radius, self._threshold),alg["update_period"])
            elif alg["type"] == "static":
                alg["update_period"] = alg["number_of_iterations"]
                alg["spimage_support_array"] = spimage.sp_support_array_init(spimage.sp_support_static_alloc(),alg["update_period"])
            else:
                self._log("No valid support algorithm set. This error should be reported!","ERROR")
                return
            if alg["center_image"]:
                spimage.sp_support_array_append(alg["spimage_support_array"],spimage.sp_support_centre_image_alloc())
                self._log("Enabling center image in real space, it will be done with every support update.")
            else:
                self._log("Center real space image disabled.")
            i += alg["number_of_iterations"]
            self._support_algorithms.append(alg)
        self._support_algorithms_dirty = False
        self._log("Support algorithms initialised.","DEBUG")

    def _init_phasing_algorithms(self):
        if not self._phasing_algorithms_dirty:
            self._log("Phasing algorithms already initialised.","DEBUG")
            return       
        i = 0
        for alg_conf in self._phasing_algorithms_configs:
            alg = dict(alg_conf)
            if alg_conf["number_of_iterations"] is None:
                if len(self._phasing_algorithms_configs) == 1:
                    alg["number_of_iterations"] = self._number_of_iterations
                else:
                    self._log("Number of iterations can not be None if many phasing algorithms are specified. Please report this error.","ERROR")
                    return
            # constraints
            constraints = spimage.SpNoConstraints
            if "constraints" in alg:
                if "enforce_positivity" in alg["constraints"] and "enforce_real" in alg["constraints"]:
                    constraints |= spimage.SpPositiveRealObject
                elif "enforce_real" in alg["constraints"]:
                    constraints |= spimage.SpRealObject
                elif "enforce_positivity" in alg["constraints"]:
                    constraints |= spimage.SpPositiveComplexObject
                if "enforce_centrosymmetry" in alg["constraints"]:
                    constraints |= spimage.SpCentrosymmetricObject
            # phasing algorithm
            if alg["type"] in ["raar","hio","diffmap"]:
                alg["beta"] = spimage.sp_smap_alloc(1)
                spimage.sp_smap_insert(alg["beta"], i,alg["beta_init"])
                spimage.sp_smap_insert(alg["beta"], i+alg["number_of_iterations"], alg["beta_final"])
            if alg["type"] == "raar":
                alg["spimage_phasing"] = spimage.sp_phasing_raar_alloc(alg["beta"], constraints)
            elif alg["type"] == "hio":
                alg["spimage_phasing"] = spimage.sp_phasing_hio_alloc(alg["beta"], constraints)
            elif alg["type"] == "diffmap":
                alg["spimage_phasing"] = spimage.sp_phasing_diff_map_alloc(alg["beta"], alg["gamma1"], alg["gamma2"], constraints)
            elif alg["type"] == "er":
                alg["spimage_phasing"] = spimage.sp_phasing_er_alloc(constraints)
            else:
                self._log("Phasing algorithm %s not implemented. Please report this error!" % type,"ERROR")
                return
            i += alg["number_of_iterations"]
            self._phasing_algorithms.append(alg)
        self._phasing_algorithms_dirty = False
        self._log("Phasing algorithm initialised.","DEBUG")

    def _init_phaser(self):
        #if not self._phaser_dirty:
        #    self._log("Phaser already initialised.","DEBUG")
        #    return
        self._clear_phaser()
        self._sp_phaser = spimage.sp_phaser_alloc()
        if self._use_gpu:
            pe = spimage.SpEngineCUDA
        else:
            pe = spimage.SpEngineCPU
        self._log("Initialising phaser with the phasing algorithm %s and the support algorithm %s." % (self._phasing_algorithms[0]["type"],self._support_algorithms[0]["type"]))
        spimage.sp_phaser_init(self._sp_phaser, self._phasing_algorithms[0]["spimage_phasing"], self._support_algorithms[0]["spimage_support_array"], pe)
        spimage.sp_phaser_set_amplitudes(self._sp_phaser, self._sp_amplitudes)
        spimage.sp_phaser_init_model(self._sp_phaser, None, spimage.SpModelRandomPhases)
        spimage.sp_phaser_init_support(self._sp_phaser, self._sp_initial_support, 0, 0)
        self._phaser_dirty = False
        self._log("Phaser initialized.","DEBUG")

    def reconstruct(self):
        """
        Performs a single reconstruction with the specified configuration and outputs a dictionary with the results.
        """
        self._prepare_reconstruction()
        if not self._ready:
            return
        out = {}
        #self._init_new_phaser()
        i_alg = 0
        i_sup = 0
        i_out_images = 0
        i_out_scores = 0
        iteration0_alg = 0
        iteration0_sup = 0
        iterations_propagated = 0
        fourier_error = np.zeros(self._number_of_outputs_scores)
        real_error = np.zeros(self._number_of_outputs_scores)
        support_size = np.zeros(self._number_of_outputs_scores)
        A = self._amplitudes
        if len(A.shape) == 2:
            real_space = np.zeros(shape=(self._number_of_outputs_images,self._Ny,self._Nx),dtype="complex128")
            support = np.zeros(shape=(self._number_of_outputs_images,self._Ny,self._Nx),dtype="bool")
            fourier_space = np.zeros(shape=(self._number_of_outputs_images,self._Ny,self._Nx),dtype="complex128")
            mask = np.zeros(shape=(self._number_of_outputs_images,self._Ny,self._Nx),dtype="bool")
            for self._iteration in range(self._number_of_iterations+1):
                change_algorithm = (self._iteration-iteration0_alg == self._phasing_algorithms[i_alg]["number_of_iterations"]) and (self._iteration != self._number_of_iterations)
                change_support = (self._iteration-iteration0_sup == self._support_algorithms[i_sup]["number_of_iterations"]) and (self._iteration != self._number_of_iterations)
                if change_algorithm or change_support or self._iteration in self._out_iterations_images or self._iteration in self._out_iterations_scores:
                    if self._reconstruction is None:
                        self._log("Iteration %i" % (self._iteration),"INFO")
                    else:
                        self._log("Reconstruction %i - Iteration %i" % (self._reconstruction,self._iteration),"INFO")
                    # iterate to self._iteration
                    if self._iteration != 0:
                        spimage.sp_phaser_iterate(self._sp_phaser,self._iteration-iterations_propagated)
                        iterations_propagated = self._iteration
                    if change_algorithm:
                        i_alg += 1
                        iteration0_alg = self._iteration
                        self._sp_phaser.algorithm = self._phasing_algorithms[i_alg]["spimage_phasing"]
                        self._log("Change of phasing algorithm to %s." % self._phasing_algorithms[i_alg]["type"],"INFO")
                    if change_support:
                        i_sup += 1
                        iteration0_sup = self._iteration
                        self._sp_phaser.sup_algorithm = self._support_algorithms[i_sup]["spimage_support_array"]
                        self._log("Change of support algorithm to %s." % self._support_algorithms[i_sup]["type"],"INFO")
                    if self._iteration in self._out_iterations_images:
                        [real_space[i_out_images,:,:],support[i_out_images,:,:]] = self._get_curr_model(shifted=False)
                        [fourier_space[i_out_images,:,:],mask[i_out_images,:,:]] = self._get_curr_fmodel(shifted=False)
                        i_out_images += 1
                        self._log("Outputting images.","DEBUG")
                    if self._iteration in self._out_iterations_scores:
                        scores = self._get_scores()
                        fourier_error[i_out_scores] = scores["fourier_error"]
                        real_error[i_out_scores] = scores["real_error"]
                        support_size[i_out_scores] = scores["support_size"]
                        i_out_scores += 1
                        self._log("Outputting scores.","DEBUG")
            self._iteration = None
            out = {"iteration_index_images":self._out_iterations_images,
                   "iteration_index_scores":self._out_iterations_scores,
                   "real_space":real_space,
                   "support":support,
                   "fourier_space":fourier_space,
                   "mask":mask,
                   "real_error":real_error,
                   "fourier_error":fourier_error,
                   "support_size":support_size}

            return out

        if len(A.shape) == 3:
            real_space = np.zeros(shape=(self._number_of_outputs_images,self._Nz,self._Ny,self._Nx),dtype="complex128")
            support = np.zeros(shape=(self._number_of_outputs_images,self._Nz,self._Ny,self._Nx),dtype="bool")
            fourier_space = np.zeros(shape=(self._number_of_outputs_images,self._Nz,self._Ny,self._Nx),dtype="complex128")
            mask = np.zeros(shape=(self._number_of_outputs_images,self._Nz,self._Ny,self._Nx),dtype="bool")
            for self._iteration in range(self._number_of_iterations+1):
                change_algorithm = (self._iteration-iteration0_alg == self._phasing_algorithms[i_alg]["number_of_iterations"]) and (self._iteration != self._number_of_iterations)
                change_support = (self._iteration-iteration0_sup == self._support_algorithms[i_sup]["number_of_iterations"]) and (self._iteration != self._number_of_iterations)
                if change_algorithm or change_support or self._iteration in self._out_iterations_images or self._iteration in self._out_iterations_scores:
                    if self._reconstruction is None:
                        self._log("Iteration %i" % (self._iteration),"INFO")
                    else:
                        self._log("Reconstruction %i - Iteration %i" % (self._reconstruction,self._iteration),"INFO")
                    # iterate to self._iteration
                    if self._iteration != 0:
                        spimage.sp_phaser_iterate(self._sp_phaser,self._iteration-iterations_propagated)
                        iterations_propagated = self._iteration
                    if change_algorithm:
                        i_alg += 1
                        iteration0_alg = self._iteration
                        self._sp_phaser.algorithm = self._phasing_algorithms[i_alg]["spimage_phasing"]
                        self._log("Change of phasing algorithm to %s." % self._phasing_algorithms[i_alg]["type"],"INFO")
                    if change_support:
                        i_sup += 1
                        iteration0_sup = self._iteration
                        self._sp_phaser.sup_algorithm = self._support_algorithms[i_sup]["spimage_support_array"]
                        self._log("Change of support algorithm to %s." % self._support_algorithms[i_sup]["type"],"INFO")
                    if self._iteration in self._out_iterations_images:
                        [real_space[i_out_images,:,:,:],support[i_out_images,:,:,:]] = self._get_curr_model(shifted=False)
                        [fourier_space[i_out_images,:,:,:],mask[i_out_images,:,:,:]] = self._get_curr_fmodel(shifted=False)
                        i_out_images += 1
                        self._log("Outputting images.","DEBUG")
                    if self._iteration in self._out_iterations_scores:
                        scores = self._get_scores()
                        fourier_error[i_out_scores] = scores["fourier_error"]
                        real_error[i_out_scores] = scores["real_error"]
                        support_size[i_out_scores] = scores["support_size"]
                        i_out_scores += 1
                        self._log("Outputting scores.","DEBUG")
            self._iteration = None
            out = {"iteration_index_images":self._out_iterations_images,
                   "iteration_index_scores":self._out_iterations_scores,
                   "real_space":real_space,
                   "support":support,
                   "fourier_space":fourier_space,
                   "mask":mask,
                   "real_error":real_error,
                   "fourier_error":fourier_error,
                   "support_size":support_size}

            return out

    def reconstruct_loop(self,Nrepeats,full_output=False):
        """
        Repeats the reconstruction Nrepeats times with the specified parameters from random seeds of starting phases and outputs the results as a dictionary.
        """
        self._prepare_reconstruction()
        if not self._ready:
            return

        A = self._amplitudes
        if len(A.shape) == 2:
            single_outputs = []
            fourier_space_final = np.zeros(shape=(Nrepeats,self._Nx,self._Ny),dtype="complex128")
            real_space_final = np.zeros(shape=(Nrepeats,self._Nx,self._Ny),dtype="complex128")
            support_final = np.zeros(shape=(Nrepeats,self._Nx,self._Ny),dtype="bool")
            scores_final = {}
            scores_final["real_error"] = np.zeros(Nrepeats)
            scores_final["fourier_error"] = np.zeros(Nrepeats)
            scores_final["support_size"] = np.zeros(Nrepeats)
            for self._reconstruction in range(Nrepeats):
                self._log("Reconstruction %i started" % (self._reconstruction),"INFO")
                output = self.reconstruct()
                real_space_final[self._reconstruction,:,:] = output["real_space"][-1,:,:]
                fourier_space_final[self._reconstruction,:,:] = output["fourier_space"][-1,:,:]
                support_final[self._reconstruction,:,:] = output["support"][-1,:,:]
                scores_final["real_error"][self._reconstruction] = output["real_error"][-1]
                scores_final["fourier_error"][self._reconstruction] = output["fourier_error"][-1]
                scores_final["support_size"][self._reconstruction] = output["support_size"][-1]
                single_outputs.append(output)
                self._log("Reconstruction %i exited" % (self._reconstruction),"INFO")
            self._reconstruction = None
            out =  {"real_space_final":real_space_final,
                    "fourier_space_final":fourier_space_final,
                    "support_final":support_final,
                    "scores_final":scores_final}

        if len(A.shape) == 3:
            single_outputs = []
            fourier_space_final = np.zeros(shape=(Nrepeats,self._Nx,self._Ny,self._Nz),dtype="complex128")
            real_space_final = np.zeros(shape=(Nrepeats,self._Nx,self._Ny,self._Nz),dtype="complex128")
            support_final = np.zeros(shape=(Nrepeats,self._Nx,self._Ny,self._Nz),dtype="bool")
            scores_final = {}
            scores_final["real_error"] = np.zeros(Nrepeats)
            scores_final["fourier_error"] = np.zeros(Nrepeats)
            scores_final["support_size"] = np.zeros(Nrepeats)
            for self._reconstruction in range(Nrepeats):
                self._log("Reconstruction %i started" % (self._reconstruction),"INFO")
                output = self.reconstruct()
                real_space_final[self._reconstruction,:,:,:] = output["real_space"][-1,:,:,:]
                fourier_space_final[self._reconstruction,:,:,:] = output["fourier_space"][-1,:,:,:]
                support_final[self._reconstruction,:,:,:] = output["support"][-1,:,:,:]
                scores_final["real_error"][self._reconstruction] = output["real_error"][-1]
                scores_final["fourier_error"][self._reconstruction] = output["fourier_error"][-1]
                scores_final["support_size"][self._reconstruction] = output["support_size"][-1]
                single_outputs.append(output)
                self._log("Reconstruction %i exited" % (self._reconstruction),"INFO")
            self._reconstruction = None
            out =  {"real_space_final":real_space_final,
                    "fourier_space_final":fourier_space_final,
                    "support_final":support_final,
                    "scores_final":scores_final}

        if full_output:
            out["single_outputs"] = single_outputs
        return out

    def _get_curr_fmodel(self, shifted=False):
        tmp = spimage.sp_phaser_fmodel(self._sp_phaser)
        fimg = tmp.image.copy()
        tmp = spimage.sp_phaser_amplitudes(self._sp_phaser)
        fmsk = tmp.mask.copy()
        if not shifted:
            return [np.fft.fftshift(fimg),np.fft.fftshift(fmsk)]
        else:
            return [fimg,fmsk]
            
    def _get_curr_model(self,shifted=False, state='before_projection'):
        if self._iteration > 0:
            if state == "before_projection":
                tmp1 = spimage.sp_phaser_model_before_projection(self._sp_phaser)
                tmp2 = spimage.sp_phaser_support(self._sp_phaser)
                img = tmp1.image.copy()
                msk = np.int32(tmp2.image.real)
                #spimage.sp_image_free(tmp1)
                #spimage.sp_image_free(tmp2)
            elif state == "after_projection":
                tmp = spimage.sp_phaser_model_with_support(self._sp_phaser)
                img = tmp.image.copy()
                msk = tmp.mask.copy()
                #spimage.sp_image_free(tmp)
        else:
            img = self._sp_phaser.model.image.copy()
            msk = np.array(self._sp_initial_support.image.real,dtype="bool")
        if not shifted:
            return [np.fft.fftshift(img),np.fft.fftshift(msk)]
        else:
            return [img,msk]


def get_test_data():
    import scipy.misc
    real = scipy.misc.lena()
    fourier = np.fft.fftshift(np.fft.ifftn(real))
    o = {}
    o["intensities"] = np.float64(abs(fourier)**2)
    o["mask"] = np.ones(shape=fourier.shape,dtype="bool")
    return o
    
