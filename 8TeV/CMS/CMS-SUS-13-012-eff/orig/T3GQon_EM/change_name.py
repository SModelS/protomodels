import os,sys

base = '/scratch/federico/SModelS_GitHub/smodels-database-TGQ/8TeV/CMS/CMS-SUS-13-012-eff/orig/T3GQon'

dirs = [ name for name in os.listdir(base) if 'TGQ' in name ]


for name in dirs:
    path=  base + '/' + name

    os.chdir(path)
    for name in os.listdir('.'):                         
        NAME = name.replace('MA5_EM_TGQon_MAPS_','MA5_EM_TGQon_MAPS_SR')
        os.system('mv ' + name + ' ' + NAME )
