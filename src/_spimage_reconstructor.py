from pylab import *
import spimage
from _spimage_log import logger

class Reconstructor:
    def __init__(self):
        # mask and intensities
        self._intensities = None
        self._mask = None
        # wrapped C-object instances
        self._sp_intensities = None
        self._sp_intensities_sh = None
        self._sp_amplitudes = None
        self._sp_amplitudes_sh = None
        self._sp_amplitudes_dirty = True
        self._sp_initial_support = None
        self._sp_initial_support_sh = None
        self._sp_initial_support_dirty = True
        self._sp_phaser = None
        self._sp_phaser_dirty = True
        # other stuff
        self.clear()
        self._log("Reconstructor initialized.","DEBUG")

    def clear(self):
        self._clear_intensities()
        self._clear_amplitudes()
        self._clear_initial_support()
        self._clear_iterations()
        self._clear_support_algorithms()
        self._clear_phasing_algorithms()
        self._clear_phaser()
        self._log("Reconstructor initialized.","DEBUG")

    def _clear_intensities(self):
        for img in [self._sp_intensities,self._sp_intensities_sh]:
            if img != None:
                spimage.sp_image_free(img)
        self._sp_intensities = None
        self._sp_intensities_sh = None
        self._amplitudes_dirty = True

    def _clear_amplitudes(self):
        for img in [self._sp_amplitudes,self._sp_amplitudes_sh]:
            if img != None:
                spimage.sp_image_free(img)
        self._sp_amplitudes = None
        self._sp_amplitudes_sh = None
        self._amplitudes_dirty = True

    def _clear_iterations(self):
        self._number_of_iterations = None
        self._number_of_outputs_images = None
        self._number_of_outputs_scores = None
        self._iterations_dirty = True       

    def _clear_initial_support(self):
        for img in [self._sp_initial_support,self._sp_initial_support_sh]:
            if img != None:
                spimage.sp_image_free(img)
        self._sp_initial_support = None
        self._sp_initial_support_sh = None
        self._initial_support_config = None
        self._initial_support_dirty = True

    def _clear_support_algorithms(self):
        self._support_algorithms = []
        self._support_algorithms_configs = []
        self._i_support_algorithms = 0
        self._support_algorithms_dirty = True

    def _clear_phasing_algorithms(self):
        self._phasing_algorithms = []
        self._phasing_algorithms_configs = []
        self._i_phasing_algorithms = 0
        self._phasing_algorithms_dirty = True

    def _clear_phaser(self):
        if self._sp_phaser != None:
            spimage.sp_phaser_free(self._sp_phaser)          
        self._sp_phaser = None
        self._phaser_dirty = True

    def _get_scores(self):
        out = {}
        [I,M] = self._get_curr_model(shifted=True)
        [fI,fM] = self._get_curr_fmodel(shifted=True)
        out["real_error"] = sqrt((abs(I[M==0])**2).sum()/( (abs(I[M==1])**2).sum() + finfo("float32").eps ))
        out["fourier_error"] = sqrt((abs(abs(fI[fM==1])-abs(fftshift(self._sp_amplitudes.image)[fM==1]))**2).sum()/( (abs(fftshift(self._sp_amplitudes.image)[fM==1])**2).sum() + finfo("float32").eps ))
        out["FcFo"] = sqrt(abs(fI[fM==1]).sum()/(abs(self._sp_amplitudes.image[fM==0]).sum()+finfo("float32").eps))
        out["support_size"] = M.sum()
        return out

    def _log(self,s,mode="INFO"):
        if logger != None:
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

    def set_intensities(self,intensities,shifted=True):
        self._intensities_dirty = True
        self._initial_support_dirty = True
        if shifted:
            self._intensities = fftshift(intensities)
        else:
            self._intensities = intensities.copy()
        self._log("Intensities set.","DEBUG")

    def set_mask(self,mask,shifted=True):
        self._intensities_dirty = True
        self._initial_support_dirty = True
        if shifted:
            self._mask = fftshift(mask)
        else:
            self._mask = mask.copy()
        self._log("Mask set.","DEBUG")

    def set_number_of_iterations(self,number_of_iterations):
        self._number_of_iterations = number_of_iterations
        self._iterations_dirty = True
        self._log("Number of iterations set.","DEBUG")

    def set_number_of_outputs_images(self,number_of_outputs_images):
        self._number_of_outputs_images = number_of_outputs_images
        self._iterations_dirty = True      
        self._log("Number of output images set.","DEBUG")

    def set_number_of_outputs_scores(self,number_of_outputs_scores):
        self._number_of_outputs_scores = number_of_outputs_scores
        self._iterations_dirty = True
        self._log("Number of output scores set.","DEBUG")

    def set_number_of_outputs(self,number_of_outputs):
        self._number_of_outputs_images = number_of_outputs
        self._number_of_outputs_scores = number_of_outputs
        self._iterations_dirty = True
        self._log("Number of outputs set.","DEBUG")

    def set_initial_support(self,**kwargs):
        if not ("radius" in kwargs or "support_mask" in kwargs):
            self._log("set_initial_support requires one of the following key word arguments: radius, support_mask","ERROR")
            return 
        self._initial_support_config = {}
        if "radius" in kwargs:
            self._initial_support_config["radius"] = kwargs["radius"]
        else:
            self._initial_support_config["support_mask"] = kwargs["support_mask"]
            self._initlal_support_congig["support_mask_shifted"] = kwargs.get("support_mask_shifted",True)
        self._log("Initial support configuration set.","DEBUG")

    def set_support_algorithm(self, type, number_of_iterations, update_period, **kwargs):
        self._clear_support_algorithms()
        self.append_support_algorithm(type,number_of_iterations,update_period,**kwargs)

    def append_support_algorithm(self, type, number_of_iterations, update_period, **kwargs):
        alg_conf = {"type":type,"number_of_iterations":number_of_iterations,"update_period":update_period}
        # check input
        if type not in ["area","threshold"]:
            self._log("append_support_algorithm accepts algorithms of the following types: area, threshold","ERROR")
            return
        if type == "area":
            for k in ["blur_init","blur_final","area_init","area_final"]:
                if k not in kwargs:
                    self._log("append_support_algorithm with the support algorithm \"area\" requires the following argument: " + k, "ERROR")
                    return
                else:
                    alg_conf[k] = kwargs[k]
        elif type == "threshold":
            for k in ["blur_radius_init","blur_radius_final","threshold_init","threshold_final"]:
                if k not in kwargs:
                    self._log("append_support_algorithm with the support algorithm \"threshold\" requires the following argument: " + k, "ERROR")
                    return
                else:
                    alg_conf[k] = kwargs[k]
        self._support_algorithms_configs.append(alg_conf)
        self._i_support_algorithms += number_of_iterations
        self._support_algorithms_dirty = True
        self._log("Support algorithms appended.","DEBUG")

    def set_phasing_algorithm(self, type, number_of_iterations, **kwargs):
        self._clear_phasing_algorithms()
        self.append_phasing_algorithm(type, number_of_iterations,**kwargs)

    def append_phasing_algorithm(self, type, number_of_iterations,**kwargs):
        alg_conf = {"type":type,"number_of_iterations":number_of_iterations}
        # check input
        alg_conf["constraints"] = kwargs.get("constratints","")
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
            else:
                alg_conf[k] = kwargs[k]
        self._phasing_algorithms_configs.append(alg_conf)
        self._i_phasing_algorithms += number_of_iterations
        self._phasing_algorithms_dirty = True
        self._log("Phasing algorithm appended.","DEBUG")

    def _prepare_reconstruction(self):
        self._ready = True
        if self._intensities == None:
            self._log("Reconstruction cannot start! You need to set the intensities.","ERROR")
            self._ready = False
        if self._mask == None:
            #self._log("Reconstruction cannot start! You need to set the mask.","ERROR")
            self._log("You did not set a mask, therefore initializing without any missing intensity values.","WARNING")
            self._ready = False
        if self._initial_support_config == None:
            self._log("Reconstruction cannot start! You need to set the initial support.","ERROR")
            self._ready = False
        if self._support_algorithms_configs == []:
            self._log("Reconstruction cannot start! You need to set the support algorithm.","ERROR")
            self._ready = False           
        if self._phasing_algorithms_configs == []:
            self._log("Reconstruction cannot start! You need to set the phasing algorithm.","ERROR")
            self._ready = False           
        if self._i_support_algorithms != self._number_of_iterations:
            self._log("Connot prepare reconstruction. Support algorithms initialised are not in line with the set number of iterations.","ERROR")
            self._ready = False
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
        I = self._intensities.copy()
        self._Nx = I.shape[1]
        self._Ny = I.shape[0]
        I[I<0] = 0
        if self._mask == None:
            M = self._mask.copy()
        else:
            M = ones(shape=I.shape,dtype="bool")
        self._clear_intensities()
        self._sp_intensities_sh = spimage.sp_image_alloc(I.shape[1],I.shape[0],1)
        self._sp_intensities_sh.image.real[:,:] = float32(I[:,:])
        self._sp_intensities_sh.mask[:,:] = int32(M[:,:])
        self._clear_amplitudes()
        self._sp_amplitudes = spimage.sp_image_shift(self._sp_intensities_sh)
        self._sp_amplitudes.image[:] = sqrt(abs(self._sp_amplitudes.image.real))
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
            if img != None:
                spimage.sp_image_free(img)
        self._sp_initial_support = spimage.sp_image_alloc(self._Ny,self._Nx,1)
        if "radius" in self._initial_support_config:
            X,Y = meshgrid(arange(self._Nx),arange(self._Ny))
            X = X-(self._Nx-1)/2.
            Y = Y-(self._Ny-1)/2.
            R = sqrt(X**2 + Y**2)
            self._phaser_dirty = True
            self._sp_initial_support.image[:] = float32(fftshift(R) < self._initial_support_config)
        else:
            self._phaser_dirty = True
            S = self._initial_support_config["support_mask"]
            if self._initial_support_config["support_mask_shifted"]:
                S = fftshift(S)
            self.sep_initial_support.image[:] = float32(S)
        self._initial_support_dirty = False
        self._log("Initial support initialised.","DEBUG")

    def _init_iterations(self):
        if not self._iterations_dirty:
            self._log("Iterations already initialised.","DEBUG")
            return
        self._out_iterations_images = int64((linspace(0,self._number_of_iterations,self._number_of_outputs_images)).round())
        self._out_iterations_scores = int64((linspace(0,self._number_of_iterations,self._number_of_outputs_scores)).round())
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
            if alg_conf["type"] == "area":
                blur_radius = spimage.sp_smap_alloc(2)
                spimage.sp_smap_insert(blur_radius, i, alg_conf["blur_init"])
                spimage.sp_smap_insert(blur_radius, i+alg_conf["number_of_iterations"], alg_conf["blur_final"])
                support_area = spimage.sp_smap_alloc(2)
                spimage.sp_smap_insert(support_area, i,alg_conf["area_init"])
                spimage.sp_smap_insert(support_area, i+alg_conf["number_of_iterations"],alg_conf["area_final"])
                alg["spimage_support_array"] = spimage.sp_support_array_init(spimage.sp_support_area_alloc(blur_radius, support_area),alg_conf["update_period"])
            elif alg_cong["type"] == "threshold":
                blur_radius = spimage.sp_smap_alloc(2)
                spimage.sp_smap_insert(blur_radius, i, alg_conf["blur_radius_init"])
                spimage.sp_smap_insert(blur_radius, i + alg_conf["number_of_iterations"], alg_conf["blur_radius_final"])
                threshold = spimage.sp_smap_alloc(2)
                spimage.sp_smap_insert(threshold, i, alg_conf["threshold_init"])
                spimage.sp_smap_insert(threshold, i + alg_conf["number_of_iterations"], alg_conf["threshold_final"])
                alg["spimage_support_array"] = spimage.sp_support_array_init(spimage.sp_support_threshold_alloc(blur_radius, threshold),alg_conf["update_period"])
            i += alg_conf["number_of_iterations"]
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
            # constraints
            constraints = spimage.SpNoConstraints
            if "constraints" in alg_conf:
                if "enforce_positivity" in alg_conf["constraints"] and "enforce_real" in alg_conf["constraints"]:
                    constraints |= spimage.SpPositiveRealObject
                elif "enforce_real" in alg_conf["constraints"]:
                    constraints |= spimage.SpRealObject
                elif "enforce_positivity" in alg_conf["constraints"]:
                    constraints |= spimage.SpPositiveComplexObject
                if "enforce_centrosymmetry" in alg_conf["constraints"]:
                    constraints |= spimage.SpCentrosymmetricObject
            # phasing algorithm
            if alg_conf["type"] in ["raar","hio","diffmap"]:
                alg["beta"] = spimage.sp_smap_alloc(1)
                spimage.sp_smap_insert(alg["beta"], i,alg_conf["beta_init"])
                spimage.sp_smap_insert(alg["beta"], i+alg_conf["number_of_iterations"], alg_conf["beta_final"])
            if alg_conf["type"] == "raar":
                alg["spimage_phasing"] = spimage.sp_phasing_raar_alloc(alg["beta"], constraints)
            elif alg_conf["type"] == "hio":
                alg["spimage_phasing"] = spimage.sp_phasing_hio_alloc(alg["beta"], constraints)
            elif alg_conf["type"] == "diffmap":
                alg["spimage_phasing"] = spimage.sp_phasing_diffmap_alloc(alg["beta"], alg_conf["gamma1"], alg_conf["gamma2"], constraints)
            elif alg_conf["type"] == "er":
                alg["spimage_phasing"] = spimage.sp_phasing_er_alloc(constraints)
            else:
                self._log("Phasing algorithm %s not implemented. Please report this error!" % type,"ERROR")
                return
            i += alg_conf["number_of_iterations"]
            self._phasing_algorithms.append(alg)
        self._phasing_algorithms_dirty = False
        self._log("Phasing algorithm initialised.","DEBUG")

    def _init_phaser(self):
        if not self._phaser_dirty:
            self._log("Phaser already initialised.","DEBUG")
            return
        self._clear_phaser()
        self._phaser = spimage.sp_phaser_alloc()
        pe = spimage.SpEngineCUDA
        spimage.sp_phaser_init(self._phaser, self._phasing_algorithms[0]["spimage_phasing"], self._support_algorithms[0]["spimage_support_array"], pe)
        spimage.sp_phaser_set_amplitudes(self._phaser, self._sp_amplitudes)
        spimage.sp_phaser_init_model(self._phaser, None, spimage.SpModelRandomPhases)
        spimage.sp_phaser_init_support(self._phaser, self._sp_initial_support, 0, 0)
        self._phaser_dirty = False
        self._log("Phaser initialized.","DEBUG")

    def reconstruct(self,i_rec=0):
        self._prepare_reconstruction()
        if not self._ready:
            return
        Out = {}
        #self._init_new_phaser()
        i_alg = 0
        i_sup = 0
        i_out_images = 0
        i_out_scores = 0
        iteration0_alg = 0
        iteration0_sup = 0
        iterations_propagated = 0
        fourier_error = zeros(self._number_of_outputs_scores)
        real_error = zeros(self._number_of_outputs_scores)
        support_size = zeros(self._number_of_outputs_scores)
        real_space = zeros(shape=(self._number_of_outputs_images,self._Ny,self._Nx),dtype="complex128")
        support = zeros(shape=(self._number_of_outputs_images,self._Ny,self._Nx),dtype="bool")
        fourier_space = zeros(shape=(self._number_of_outputs_images,self._Ny,self._Nx),dtype="complex128")
        mask = zeros(shape=(self._number_of_outputs_images,self._Ny,self._Nx),dtype="bool")
        for self._iteration in range(self._number_of_iterations+1):
            change_algorithm = (self._iteration-iteration0_alg == self._phasing_algorithms[i_alg]["number_of_iterations"]) and (self._iteration != self._number_of_iterations)
            change_support = (self._iteration-iteration0_sup == self._support_algorithms[i_sup]["number_of_iterations"]) and (self._iteration != self._number_of_iterations)
            if change_algorithm or change_support or self._iteration in self._out_iterations_images or self._iteration in self._out_iterations_scores:
                self._log("Reconstruction %i - Iteration %i" % (i_rec,self._iteration),"INFO")
                # iterate to self._iteration
                if self._iteration != 0:
                    spimage.sp_phaser_iterate(self._phaser,self._iteration-iterations_propagated)
                    iterations_propagated = self._iteration
                if change_algorithm:
                    i_alg += 1
                    iteration0_alg = self._iteration
                    self._phaser.algorithm = self._phasing_algorithms[i_alg]["spimage_phasing"]
                if change_support:
                    i_sup += 1
                    iteration0_sup = self._iteration
                    self._phaser.sup_algorithm = self._support_algortihms[i_sup]["spimage_support_array"]
                if self._iteration in self._out_iterations_images:
                    [real_space[i_out_images,:,:],support[i_out_images,:,:]] = self._get_curr_model(shifted=True)
                    [fourier_space[i_out_images,:,:],mask[i_out_images,:,:]] = self._get_curr_fmodel(shifted=True)
                    i_out_images += 1
                if self._iteration in self._out_iterations_scores:
                    scores = self._get_scores()
                    fourier_error[i_out_scores] = scores["fourier_error"]
                    real_error[i_out_scores] = scores["real_error"]
                    support_size[i_out_scores] = scores["support_size"]
                    i_out_scores += 1
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
    
    def reconstruct_loop(self,Nrepeats):
        self._prepare_reconstruction()
        if not self._ready:
            return
        outputs = []
        real_space_final = zeros(shape=(Nrepeats,self._Nx,self._Ny),dtype="complex128")
        support_final = zeros(shape=(Nrepeats,self._Nx,self._Ny),dtype="bool")
        for i_rec in range(Nrepeats):
            self._log("Reconstruction %i started" % (i_rec),"INFO")
            outputs = self.reconstruct(i_rec)
            real_space_final[i_rec,:,:] = outputs["real_space"][-1,:,:]
            support_final[i_rec,:,:] = outputs["support"][-1,:,:]
            self._log("Reconstruction %i exited" % (i_rec),"INFO")
        out =  {"single_outputs":outputs,
                "real_space_final":real_space_final,
                "support_final":support_final}
        return out

    def _get_curr_fmodel(self,**kwargs):
        shifted = kwargs.get("shifted",False)
        #normalize = kwargs.get("normalize",False)
        if self._iteration > 0:
            #fimg = spimage.sp_phaser_fmodel(self._phaser).image.copy()
            tmp = spimage.sp_phaser_fmodel_with_mask(self._phaser)
            fimg = tmp.image.copy()
            fmsk = tmp.mask.copy()
        else:
            fimg = self._phaser.model.image.size*fftn(self._phaser.model.image)
            fmsk = self._sp_amplitudes.mask.copy()           
        #if normalize:
        #    fimg = fimg/abs(fimg).sum()
        if shifted:
            return [fftshift(fimg),fftshift(fmsk)]
        else:
            return [fimg,fmsk]
            
    def _get_curr_model(self,**kwargs):
        shifted = kwargs.get("shifted",False)
        state = kwargs.get("state","before_projection")
        if self._iteration > 0:
            if state == "before_projection":
                tmp1 = spimage.sp_phaser_model_before_projection(self._phaser)
                tmp2 = spimage.sp_phaser_support(self._phaser)
                img = tmp1.image.copy()
                msk = int32(tmp2.image.real)
                #spimage.sp_image_free(tmp1)
                #spimage.sp_image_free(tmp2)
            elif state == "after_projection":
                tmp = spimage.sp_phaser_model_with_support(self._phaser)
                img = tmp.image.copy()
                msk = tmp.mask.copy()
                #spimage.sp_image_free(tmp)
        else:
            img = self._phaser.model.image.copy()
            msk = array(self._sp_initial_support.image.real,dtype="bool")
        if shifted:
            return [fftshift(img),fftshift(msk)]
        else:
            return [img,msk]


def get_test_data():
    import scipy.misc
    real = scipy.misc.lena()
    fourier = fftshift(ifftn(real))
    o = {}
    o["intensities"] = float64(abs(fourier)**2)
    o["mask"] = ones(shape=fourier.shape,dtype="bool")
    return o
    