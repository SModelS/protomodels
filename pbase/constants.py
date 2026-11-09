""" physics constants, FIXME clean up
"""

## vector boson masses and widths
#mass_W = 80.377
#mwidth_W = 2.14
#mass_Z = 91.1876
#mwidth_Z = 2.5

# for widths we will consider only values of 0 (stable),
# 1 (promptly decaying) and values in the range given below
# lmin = log_10 ( min )

# electroweakinos
LSP = 1000022  # N1
X1Z = 1000022  # N1
X2Z = 1000023  # N2
X1W = 1000024  # C1
X3Z = 1000025  # N3
X4Z = 1000035  # N4
X2W = 1000037  # C2

# others
Xd = 1000001 # sdown
Xu = 1000002 # sup
Xs = 1000003 # sstrange
Xc = 1000004 # scharm
Xb = 1000005 # sbottom
Xt = 1000006 # sstop
Xb2 = 2000005 # sbottom
Xt2 = 2000006 # sstop
Xg = 1000021 # gluino

llpWidthRange = { "min": 1e-23, "max": 1e-13, "lmin": -23, "lmax": -13 }

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
