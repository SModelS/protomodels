#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import scipy
import pandas as pd
import glob
import os
from smodels.base.physicsUnits import pb,GeV

class XSecFitter():
    
    def __init__(self, pid_pair, sqrt):

        self.pid = pid_pair
        self.sqrt = sqrt

        #self.computeRelativeDifference(xsec_dict,xsec_dict8)
        #self.data_13 = pd.DataFrame(xsec_dict13)

        
    def computeRelativeDifference(self,xs13, xs8):
        #print(xs13.keys())
        #print(xs8.keys())
        
        xsec_8 = [val for key,val in xs8.items() if key in xs13.keys() ]
        xsec_13 = [val for val in xs13.values()]
        
        diff_8 = np.array([xsec_8[i-1] - xsec_8[i] for i in range(1,len(xsec_8))])
        diff_13 = np.array([xsec_13[i-1] - xsec_13[i] for i in range(1,len(xsec_13))])
        
        diff_factors = (diff_13 - diff_8)/diff_13
        sum_diff_factors = np.sum(diff_factors)/len(xsec_13)
        print("\n8 Tev: ", diff_8)
        print("\n13 Tev: ", diff_13)
        print("\n Relative Difference: ", diff_factors)
        print("\nAverage of Difference:", sum_diff_factors)
    
    
    
    def getData(self, filename):
        
        fname = filename.split('/')[-1].replace("csv", "txt")
        files = glob.glob(f"xsecTables/{fname}")
        
        file = open(files[0], 'r')
        lines = file.readlines()

        out = open(filename,"w")
        
        for line in lines:
            if '#' in line:
                if 'Mass' in line: out.write("Mass(GeV),Xsec(pb)\n")
                else: continue
            else:
                line = line.replace("\t","").replace("GeV",",")
                if '±' in line: index = line.rfind('±')
                elif '+' in line: index = line.rfind('+')
                line = line[:index-1].replace(" ","") + "\n"
                out.write(line)
                

        out.close()
        
        self.data = pd.read_csv(filename)
        
    def plotfig(self, x, y, bs, inverse):
          
        '''
        xsec_pid = {"stoppb" : [1000006], "gluino":[1000021], "squark": [1000001, 1000002, 1000003, 1000004, 1000005], "slep": [1000011, 1000013, 1000015],
                    "N2C1p": [1000023, 1000024], "N2C1m": [1000023, -1000024], "C1C1": [1000024, -1000024] }
        
        xsec_unit = {"fb":[],"pb":["stop",]}
        prod = [pname for key,value in xsec_pid.items() if particle_id in value]
        fname13 = f"xsecTables/plots/xsec{prod[0]}13.csv"
        fname8 = f"xsecTables/plots/xsec{prod[0]}8.csv"

        if os.path.exists(fname13): self.data = pd.read_csv(fname13)
        else: self.getData(fname13)
        
        if os.path.exists(fname8): self.data = pd.read_csv(fname8)
        else: self.getData(fname8)
        '''
        prod_name = {"stop" : [1000006], "gluino":[1000021], "squark-antisquark": [1000001, 1000002, 1000003, 1000004, 1000005], "slepton-pair": [1000011, 1000013, 1000015],
                    "chargino-neutralino": [1000023], "chargino-chargino": [1000024]}
        
        pname = [key for key,value in prod_name.items() if self.pid[1] in value]
        plot_title = pname[0]+"_"+str(self.sqrt)
        #if self.pid == 1000006: plot_title = "stop"
        #elif self.pid == 1000021: plot_title = "gluino"
        
        x_lab = "Mass (GeV)"
        y_lab = "Xsec (pb)"
        fig_title = "fit_function"
        if inverse:
            x_lab = "Xsec (pb)"
            y_lab = "Mass (GeV)"
            fig_title = "inv_fit_function"
        
        plt.clf()
        plt.scatter(x, y, label="Data")
        plt.plot(x, bs(x), color='r', label="Fit using B-Spline")
        plt.xlabel(x_lab)
        #plt.xlim(120,2000)
        #plt.ylim(10**(-6), 10**(4))
        plt.ylabel(y_lab)
        plt.title(plot_title)
        plt.legend()
        #plt.show()
        plt.savefig(f"xsecTables/plots/{plot_title}_{fig_title}.png")
        plt.show()


    def getValueFromFit(self, input_value, inverse, show_plot=False):
        
        from refxsecComputer import RefXSecComputer
        xsec_obj = RefXSecComputer()
        
        xsec_dict = xsec_obj.getXSecsFor(pid1 = self.pid[0], pid2 = self.pid[1], sqrts=self.sqrt, ewk=None, masses=[200.0, 300.0])[0]
        if xsec_dict is None: return None
        #xsec_dict8  = xsec_obj.getXSecsFor(pid1 = pid_pair[0], pid2 = pid_pair[1], sqrts=8, ewk=None, masses=[200.0, 300.0])[0]
        xsec_mass_dict = {'Mass(GeV)':xsec_dict.keys(), 'Xsec(pb)': xsec_dict.values()}
        self.data  = pd.DataFrame(xsec_mass_dict)
        
        data_x = np.array(self.data["Mass(GeV)"])
        data_y = np.array(self.data["Xsec(pb)"])
        
        if inverse:
            ind = np.argsort(data_y)
            plot_x = np.sort(data_y)
            plot_y = np.array([data_x[i] for i in ind])

            from scipy.interpolate import make_interp_spline
            bspl = make_interp_spline(plot_x, plot_y, k=3)
            xsec = input_value.asNumber(pb)
            mass = bspl(xsec).item()*GeV
            if show_plot: self.plotfig(plot_x, plot_y, bspl, inverse)
            return mass

        plot_x = data_x
        plot_y = data_y
        
        from scipy.interpolate import make_interp_spline
        
        bspl = make_interp_spline(plot_x, plot_y, k=3)
        mass = input_value.asNumber(GeV)
        xsec = bspl(mass).item()*pb
        if show_plot: self.plotfig(data_x, data_y, bspl, inverse)
       
        return xsec
    
if __name__ == '__main__':
    
    func = findXsecFit(pid_pair=(-1000024,1000023), sqrt=13)
    #func.computeRelativeDifference()
    '''
    xsec = 10.0*pb
    #m = 500*GeV
    mass = func.getValueFromFit(xsec, inverse= True, show_plot=True)
    #xs = func.getValueFromFit(m, inverse= False, show_plot=True)
    print("Mass for cross-sec 10.0 pb is %s"%(mass._value))
    #print(f"Xsec for mass {m} GeV is %s"%(xs._value))
    '''
