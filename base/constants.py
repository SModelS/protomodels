""" physics constants, FIXME clean up
"""

## vector boson masses and widths
#mass_W = 80.377
#mwidth_W = 2.14
#mass_Z = 91.1876
#mwidth_Z = 2.5

particleNames = { "W": 24, "Z": 23, "t": 6, "h": 25, "tau": 15,
                  "c": 4, "b": 5 }

smMasses = {  "t": 173., "W": 80.377, "Z": 91.1876, "h": 125.,
              "tau": 1.77, "c": 1.2, "b": 5.0 }

smWidths = { "W": 2.14, "Z": 2.5 }

for name, pid in particleNames.items():
    if name in smMasses:
        smMasses[pid]=smMasses[name]
    if name in smWidths:
        smWidths[pid]=smWidths[name]
