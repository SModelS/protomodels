from smodels.base.physicsUnits import UnitTime, UnitEnergy, eV, ns, GeV

hbar = 6.582119569*10**-16 * eV * ( 10**9 * ns )

def lifetimeToWidth ( lifetime : UnitTime ) -> UnitEnergy:
    """ convert a lifetime to a width """
    width = hbar / lifetime
    return width

def widthToLifetime ( width : UnitEnergy ) -> UnitTime:
    """ convert a lifetime to a width """
    lifetime = hbar / width
    return lifetime
