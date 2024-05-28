#!/usr/bin/env python3

from combinationsmatrix import getMatrix
import yaml

matrix = getMatrix()

ATLAS8, ATLAS13, CMS8, CMS13 = {}, {}, {}, {}

for key, value in matrix.items():
    if key.startswith('ATLAS'):
        if 12 <= int(key.split('-')[2][-2:]) <= 14:
            ATLAS8[key] = value
        elif 15 <= int(key.split('-')[2][-2:]) <= 24:
            ATLAS13[key] = value
        else:
            print('Which Run is this:',key,'?!')
    elif key.startswith('CMS'):
        if 12 <= int(key.replace('-PAS','').split('-')[2][-2:]) <= 14:
            CMS8[key] = value
        elif 15 <= int(key.replace('-PAS','').split('-')[2][-2:]) <= 24:
            CMS13[key] = value
        else:
            print('Which Run is this:',key,'?!')
    else:
        print('Not ATLAS nor CMS:',key,'?!')

with open('ATLAS8.yaml','w') as outputFile:
    outputFile.write('---\n')
    outputFile.write(yaml.dump(ATLAS8))
    outputFile.write('...')

with open('ATLAS13.yaml','w') as outputFile:
    outputFile.write('---\n')
    outputFile.write(yaml.dump(ATLAS13))
    outputFile.write('...')

with open('CMS8.yaml','w') as outputFile:
    outputFile.write('---\n')
    outputFile.write(yaml.dump(CMS8))
    outputFile.write('...')

with open('CMS13.yaml','w') as outputFile:
    outputFile.write('---\n')
    outputFile.write(yaml.dump(CMS13))
    outputFile.write('...')
