"""
.. module:: particle
   :synopsis: Defines the particle class and particle list class, their methods and related functions

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

#from smodels.tools.smodelsLogging import logger
#from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.physicsUnits import GeV
from copy import copy


class Particle(object):
	
    """ An instance of this class represents a single particle. 
        The properties are: label, pdg, mass, electric charge, color charge, width 
    """

    def __init__(self, **kwargs):
	    
        """ 
        Initializes the particle.
        Possible properties: 
        Z2parity: str, 'even' or 'odd'
        label: str, e.g. 'e-'
        pdg: number in pdg
        mass: mass of the particle
        echarge: electric charge as multiples of the unit charge
        colordim: color dimension of the particle 
        spin: spin of the particle
        width: total width
        decays: possible decays in pyslha format e.g. [ 0.5[1000022, -1], 0.5[1000022, -2] ]         
        """  
        
        #if kwargs is None or not 'Z2parity' in kwargs:
        #    logger.error("The Z2-Parity (zParity) of a particle must always be defined")
        #    raise SModelSError

        for key,value in kwargs.items():
            self.addProperty(key,value)

    def __cmp__(self,other): 
        """
        Compares particle with other.
        The comparison is done based on the label.
        :param other:  particle to be compared (Particle object)
        :return: -1 if self < other, 0 if self == other, +1, if self > other.

        Comparison to '*' always returns True. (= Particle Wild Card)
        """    
        
        if not isinstance(other,ParticleList) and not isinstance(other,Particle):
            return +1
        
        if isinstance(other,ParticleList):
            return other.__cmp__(self)        
                
        if hasattr(self, 'label') and hasattr(other, 'label'):
            if self.label != other.label:
            
                if ( self.label == '*' ) or ( other.label == '*' ): return 0                
            
                comp = self.label > other.label
                if comp: return 1
                else: return -1
       
            else:
                return 0  #Particles are equal
                 

    def __lt__( self, p2 ):
        return self.__cmp__ ( p2 ) == -1

    def __eq__( self, p2 ):
        return self.__cmp__ ( p2 ) == 0
        

    def __str__(self): 
        if hasattr(self, 'label'):
            return self.label
        else: return ''
    
    def describe(self):
        return str(self.__dict__)
    
    
    def addProperty(self,label,value,overwrite=False):
        """
        Add property with label and value.
        If overwrite = False and property already exists, it will not be overwritten
        :parameter label: property label (e.g. "mass", "label",...)
        :parameter value: property value (e.g. 171*GeV, "t+",...)
        """
        
        if overwrite:
            setattr(self, label, value)
        elif not hasattr(self, label) or getattr(self, label) is None:
            setattr(self, label, value)    




    def chargeConjugate(self):
        """
        Returns the charge conjugate particle (flips the sign of eCharge).        
        If it has a pdg property also flips its sign.
        The charge conjugate name is defined as the original name plus "~" or
        if the original name ends in "+" ("-"), it is replaced by "-" ("+")

        :return: the charge conjugate particle (Particle object)
        """

        pConjugate = copy(self)
        if hasattr(pConjugate, 'pdg'):
            pConjugate.pdg *= -1
        if hasattr(pConjugate, 'eCharge') and not pConjugate.eCharge: return pConjugate    # for neutral particle only conjugate pdg      
        
        if hasattr(pConjugate, 'eCharge'):
            pConjugate.eCharge *= -1    
        if hasattr(pConjugate, 'label'):    
            if pConjugate.label[-1] == "+":
                pConjugate.label = pConjugate.label[:-1] + "-"
            elif pConjugate.label[-1] == "-":
                pConjugate.label = pConjugate.label[:-1] + "+"
            #elif pConjugate.label[-1] == "~":
            #    pConjugate.label = pConjugate.label[:-1]
            #else: pConjugate.label = pConjugate.label + "~"            
                
        return pConjugate



    def isStable(self):
        """
        Return True if particle is stable, false otherwise.
        A particle is considered stable if it has zero width or if
        the width has not been defined
        :return: True/False
        """
        if self.width is None:
            return True
        elif type(self.width) == type(1.*GeV) and self.width > 0.*GeV:
            return False
        else:
            return True 


class ParticleList(object):

    """ An instance of this class represents a list of particle object to allow for inclusive expresions such as jet. 
        The properties are: label, pdg, mass, electric charge, color charge, width 
    """
    
    def __init__(self, label, particles):

        """ 
        Initializes the particle list.

        :parameter label: string, e.g. 'jet'
        :parameter particles: particle objects to be added to list
        """        
        self.label = label
        self.particles = particles
        
     
    def __cmp__(self,other):
        """
        Compares the list with another list or a single particle.
        :param other: particle list (ParticleList) or single particle (Particle) to be compared
        :return: If other is a single particle, returns:
                    -1 if all particles in self < other, 
                     0 if any particle in self == other,
                    +1, if any particle in self > other.
                If other is a particle list, returns 
                -1 if self < other, 0 if self == other, +1, if self > other.     
        """        
        
        if isinstance(other,ParticleList):
            if self.particles != other.particles:
                comp = self.particles > other.particles
                if comp: return 1
                else: return -1 
            else:
                return 0
        
        if isinstance(other,Particle):
            if other in self.particles:
                return 0
            else:
                for p in self.particles:
                    if p > other: return +1
                return -1


    def __str__(self): 
        return self.label
    
    def describe(self):
        return str(self.__dict__)  
        
        
    def getPdgs(self):
        """
        pdgs appearing in ParticleList
        :return: list of pgds of particles in the ParticleList
        """        
        pdgs = [particle.pdg for particle in self.particles]
        return pdgs
        
        
    def getLabels(self):
        """
        labels appearing in ParticleList
        :return: list of labels of particles in the ParticleList
        """
        labels = [particle.label for particle in self.particles]
        return labels
        