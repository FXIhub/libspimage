import os
import numpy
import pickle
import spimage
import sys

class Material:
    def __init__(self,**kwargs):
        if "massdensity" in kwargs:
            self.materialtype = "custom"
            for key in kwargs:
                if key[0] == 'c' or key == 'massdensity':
                    exec("self." + key + " = kwargs[key]")
                else:
                    logger.error("%s is no valid argument for custom initialization of Material." % key)
                    return
        elif "material_type" in kwargs:
            self.material_type = kwargs['material_type']
            self.massdensity = DICT_massdensity[self.material_type]
            self.cH = DICT_atomic_composition[self.material_type][0]
            self.cC = DICT_atomic_composition[self.material_type][1]
            self.cN = DICT_atomic_composition[self.material_type][2]
            self.cO = DICT_atomic_composition[self.material_type][3]
            self.cP = DICT_atomic_composition[self.material_type][4]
            self.cS = DICT_atomic_composition[self.material_type][5]
            self.cAu = DICT_atomic_composition[self.material_type][6]
        else:
            logger.error("No valid arguments for Material initialization.")
            return

    def get_fX(self,element,photon_energy_eV):
        """
        Get the scattering factor for an element through linear interpolation.
        """
        SF_X = DICT_scattering_factors[element]
        e = spimage.DICT_physical_constants['e']
        c = spimage.DICT_physical_constants['c']
        h = spimage.DICT_physical_constants['h']
        f1 = numpy.interp(photon_energy_eV,SF_X[:,0],SF_X[:,1])
        f2 = numpy.interp(photon_energy_eV,SF_X[:,0],SF_X[:,2])
        return complex(f1,f2) 
 
    def get_n(self,photon_energy_eV):
        """
        Obtains complex refractive index.
        Henke (1994): n = 1 - r_0/(2pi) lambda^2 sum_q rho_q f_q(0)
        r_0: classical electron radius
        rho_q: atomic number density of atom species q
        f_q(0): atomic scattering factor (forward scattering) of atom species q
        """

        re = spimage.DICT_physical_constants['re']
        h = spimage.DICT_physical_constants['h']
        c = spimage.DICT_physical_constants['c']
        qe = spimage.DICT_physical_constants['e']

        photon_wavelength = h*c/photon_energy_eV/qe

        f = self.get_f(photon_energy_eV)
        atom_density = self.get_atom_density()
        
        n = 1 - re/2/numpy.pi * photon_wavelength**2 * f * atom_density

        return n

    def get_dn(self,photon_energy_eV):
        return (1-self.get_n(photon_energy_eV))

    # convenience functions
    # n = 1 - delta - i beta
    def get_delta(self,photon_energy_eV):
        return (1-self.get_n(photon_energy_eV).real)
    def get_beta(self,photon_energy_eV):
        return (-self.get_n(photon_energy_eV).imag)
    def get_photoabsorption_cross_section(self,photon_energy_eV):
        re = spimage.DICT_physical_constants['re']
        h = spimage.DICT_physical_constants['h']
        c = spimage.DICT_physical_constants['c']
        qe = spimage.DICT_physical_constants['e']
        photon_wavelength = h*c/photon_energy_eV/qe
        mu = 2*re*photon_wavelength*self.get_f(photon_energy_eV).imag
        return mu
    def get_transmission(self,thickness,photon_energy_eV):
        n = self.get_n(photon_energy_eV)
        mu = self.get_photoabsorption_cross_section(photon_energy_eV)
        rho = self.get_atom_density()
        return numpy.exp(-rho*mu*thickness)

    def get_f(self,photon_energy_eV):

        h = spimage.DICT_physical_constants['h']
        c = spimage.DICT_physical_constants['c']
        qe = spimage.DICT_physical_constants['e']

        photon_wavelength = h*c/photon_energy_eV/qe

        atomic_composition = self.get_atomic_composition_dict()

        f_sum = 0
        for element in atomic_composition.keys():
            # sum up average atom factor
            f = self.get_fX(element,photon_energy_eV)
            f_sum += atomic_composition[element] * f
        
        return f_sum


    def get_atom_density(self):
                
        u = spimage.DICT_physical_constants['u']

        atomic_composition = self.get_atomic_composition_dict()

        M = 0
        for element in atomic_composition.keys():
            # sum up mass
            M += atomic_composition[element]*DICT_atomic_mass[element]*u

        number_density = self.massdensity/M
        
        return number_density


    def get_electron_density(self):

        u = spimage.DICT_physical_constants['u']

        atomic_composition = self.get_atomic_composition_dict()

        M = 0
        Q = 0
        for element in atomic_composition.keys():
            # sum up electrons
            M += atomic_composition[element]*DICT_atomic_mass[element]*u
            Q += atomic_composition[element]*DICT_atomic_number[element]

        electron_density = Q*self.massdensity/M
        
        return electron_density
        
        
    def get_atomic_composition_dict(self):

        atomic_composition = {}
        
        for key in self.__dict__.keys():
            if key[0] == 'c':
                atomic_composition[key[1:]] = getattr(self, key)
 
        tmp_sum = float(sum(atomic_composition.values()))
        for element in atomic_composition.keys():
            atomic_composition[element] /= tmp_sum 
        
        return atomic_composition

# Realative atomic compositions of certain material types (order: H,C,N,O,P,S,Au)
DICT_atomic_composition = {'protein':[86,52,13,15,0,3,0],
                           'cell':[23,3,1,10,0,1,0], # Bergh et al. 2008
                           'latex':[1,1,0,0,0,0,0], 
                           'water':[2,0,0,1,0,0,0], 
                           'dna':[11,10,4,6,1,0,0],
                           'lipid':[69,36,0,6,1,0,0],
                           'genophore':[205,134,38,48,3,6,0],
                           'virus':[72.43,49.85,16.32,24.49,2.57,1.39,0],
                           'mimivirus':[23,3,1,10,0,1,0],
                           'carboxysome':[0.51,0.30,0.07,0.10,0.,0.02,0],
                           'sucrose':[22,12,0,11,0,0,0],
                           'gold':[0,0,0,0,0,0,1]}

# Estimated mass densities of certain material types
DICT_massdensity = {'protein':1350,
                    'cell':1000,
                    'latex':1050,
                    'water':998,
                    'gold':19300,
                    'dna':1700,
                    'lipid':1000,
                    'genophore':1560,
                    'virus':1381,
                    'mimivirus':1100,
                    'carboxysome':1250,
                    'sucrose':1587}
DICT_atomic_number = {'H':1,
                      'He':2,
                      'Li':3,
                      'Be':4,
                      'B':5,
                      'C':6,
                      'N':7,
                      'O':8,
                      'F':9,
                      'Ne':10,
                      'Na':11,
                      'Mg':12,
                      'Al':13,
                      'Si':14,
                      'P':15,
                      'S':16,
                      'Cl':17,
                      'Ar':18,
                      'K':19,
                      'Ca':20,
                      'Sc':21,
                      'Ti':22,
                      'V':23,
                      'Cr':24,
                      'Mn':25,
                      'Fe':26,
                      'Co':27,
                      'Ni':28,
                      'Cu':29,
                      'Zn':30,
                      'Ga':31,
                      'Ge':32,
                      'As':33,
                      'Se':34,
                      'Br':35,
                      'Kr':36,
                      'Rb':37,
                      'Sr':38,
                      'Y':39,
                      'Zr':40,
                      'Nb':41,
                      'Mo':42,
                      'Tc':43,
                      'Ru':44,
                      'Rh':45,
                      'Pd':46,
                      'Ag':47,
                      'Cd':48,
                      'In':49,
                      'Sn':50,
                      'Sb':51,
                      'Te':52,
                      'I':53,
                      'Xe':54,
                      'Cs':55,
                      'Ba':56,
                      'La':57,
                      'Ce':58,
                      'Pr':59,
                      'Nd':60,
                      'Pm':61,
                      'Sm':62,
                      'Eu':63,
                      'Gd':64,
                      'Tb':65,
                      'Dy':66,
                      'Ho':67,
                      'Er':68,
                      'Tm':69,
                      'Yb':70,
                      'Lu':71,
                      'Hf':72,
                      'Ta':73,
                      'W':74,
                      'Re':75,
                      'Os':76,
                      'Ir':77,
                      'Pt':78,
                      'Au':79,
                      'Hg':80,
                      'Tl':81,
                      'Pb':82,
                      'Bi':83,
                      'Po':84,
                      'At':85,
                      'Rn':86,
                      'Fr':87,
                      'Ra':88,
                      'Ac':89,
                      'Th':90,
                      'Pa':91,
                      'U':92,
                      'Np':93,
                      'Pu':94,
                      'Am':95,
                      'Cm':96,
                      'Bk':97,
                      'Cf':98,
                      'Es':99,
                      'Fm':100,
                      'Md':101,
                      'No':102,
                      'Lr':103,
                      'Rf':104,
                      'Db':105,
                      'Sg':106,
                      'Bh':107,
                      'Hs':108,
                      'Mt':109,
                      'Ds':110,
                      'Rg':111,
                      'Cp':112,
                      'Uut':113,
                      'Uuq':114,
                      'Uup':115,
                      'Uuh':116,
                      'Uus':117,
                      'Uuo':118}

def unpickle_scattering_factors():
    global DICT_atomic_mass
    DICT_atomic_mass = {}
    global DICT_scattering_factors
    DICT_scattering_factors = {}
    this_dir = os.path.dirname(os.path.realpath(__file__))
    with open('%s/elements.dat' % this_dir,'rb') as ELEMENTS_FILE:
        if sys.version_info >= (3, 0):
            DICT_atomic_mass,DICT_scattering_factors = pickle.load(ELEMENTS_FILE, encoding="latin1")
        else:
            DICT_atomic_mass,DICT_scattering_factors = pickle.load(ELEMENTS_FILE)
    F_MIN_ENERGY_EV = 0
    F_MAX_ENERGY_EV = 0
    for var in DICT_scattering_factors.values():
        if F_MIN_ENERGY_EV < var[0,0] or F_MIN_ENERGY_EV == 0: F_MIN_ENERGY_EV = var[0,0]
        if F_MAX_ENERGY_EV > var[-1,0] or F_MAX_ENERGY_EV == 0: F_MAX_ENERGY_EV = var[-1,0]

unpickle_scattering_factors()
