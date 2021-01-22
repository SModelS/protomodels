#!/usr/bin/env python3

"""

.. module:: approximate likelihood analysis
   :synopsis: code to study the approximate analysis, come up with a conservative
              procedure

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com> 

"""

# coding: utf-8

import sys,os,copy,glob,time
import numpy as np
sys.path.append(os.path.abspath('../smodels'))
from smodels.tools.physicsUnits import fb
from smodels.tools import runtime
from smodels.tools import statistics
from smodels.tools.simplifiedLikelihoods import UpperLimitComputer, Data
runtime._experimental = True
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats,integrate,optimize,special
sns.set_style('ticks',{'font.family':'Times New Roman', 'font.serif':'Times New Roman'})
sns.set_context('paper', font_scale=2.0)
sns.set_palette(sns.color_palette("Paired"))

# ### Define number of observed, expected and observed events:
# ### Compute exact  likelihood

def likelihood(mu,nsig,nExp,nExpErr,nobs):
    """Marginalized likelihood for mu"""
    
    def integrand(nbg):
        nTot = nbg + mu*nsig
        nErr = nExpErr
        p = stats.poisson.pmf(k=nobs,mu=nTot)*stats.norm.pdf(x=nbg,loc=nExp,scale=nErr)
        return p
        
    #Marginalize over the background uncertainty
    result = integrate.quad(integrand, 0, max(2*(nobs-nExp),nExp+5*nExpErr))
    
    return result[0]

def p(mu,nsig,nExp,nExpErr,nobs):
    """Integral of the likelihood from zero to mu"""
    
    result = integrate.quad(likelihood, 0, mu,args=(nsig,nExp,nExpErr,nobs,))
    return result[0]

# ### Compute observed and expected upper limits

def run ( nobs, nExp, nExpErr, nsig ):
    """ run the procedure, with:
    :param nobs: number of observed
    :param nExp: number expected
    :param nExpErr: error on number expected
    """
    print ( "Starting run with", nobs, nExp )

    norm = p((nobs-nExp)+5*nExpErr,nsig,nExp,nExpErr,nobs)
    t0=time.time()
    ULobs = optimize.brentq(lambda x: p(x,nsig,nExp,nExpErr,nobs)/norm - 0.95,
                            0,2*(nobs-nExp+3*nExpErr))
    t1=time.time()
    print ( "ULobs", ULobs, "dt=",t1-t0 )
    data = Data ( nobs, nExp, nExpErr**2, nsignal = nsig )
    computer = UpperLimitComputer ( 10000 )
    marginalize = True
    ULobs2 = computer.ulSigma ( data, marginalize=marginalize )
    t2=time.time()
    print ( "ULobs2", ULobs2, "dt=", t2-t1 )

    normExp = p(5*nExp,nsig,nExp,nExpErr,nExp)
    ULexp = optimize.brentq(lambda x: p(x,nsig,nExp,nExpErr,nExp)/normExp - 0.95,
                            0,2*(3*nExpErr))
    t3=time.time()
    print ( "ULexp", ULexp, "dt=",t3-t2 )
    ULexp2 = computer.ulSigma ( data, expected=True, marginalize=marginalize )
    t4=time.time()
    print ( "ULexp2", ULexp2, "dt=",t4-t3 )


    print(r'Nobs = %1.2f, Nbg = %1.2f +- %1.2f, Nsig < %1.2f, Nsig (expected) < %1.2f'
         %(nobs,nExp,nExpErr,ULobs,ULexp))

    sigma = ULexp/1.96
    sigma2 = ULexp2/1.96
    mu0 = optimize.brentq(lambda x: 0.95- (special.erf((ULobs-x)/(np.sqrt(2)*sigma)) 
                     +  special.erf(x/(np.sqrt(2)*sigma)))/(1+special.erf(x/(np.sqrt(2)*sigma))),
                        0,2*(nobs-nExp+3*nExpErr))
    mu0B = optimize.brentq(lambda x: 0.95- (special.erf((ULobs2-x)/(np.sqrt(2)*sigma2)) 
                     +  special.erf(x/(np.sqrt(2)*sigma2)))/(1+special.erf(x/(np.sqrt(2)*sigma2))),
                        0,2*(nobs-nExp+3*nExpErr))
    print ( "mu0", mu0, "mu0B", mu0B )

    def llhdFromLimits(mu,mu0,sigma):
        return stats.norm.pdf(x=mu,loc=mu0,scale=sigma)

    normLim = 1 - stats.norm.cdf(0,loc=mu0,scale=sigma)
    normLim2 = 1 - stats.norm.cdf(0,loc=mu0B,scale=sigma2)

    nsteps = 100 # 100
    muvals = np.linspace(0,nobs-nExp+4*nExpErr, nsteps )

    llhds = np.array([[mu,likelihood(mu,nsig,nExp,nExpErr,nobs)/norm] for mu in muvals])
    llhdsApp = np.array([[mu,llhdFromLimits(mu,mu0,sigma)/normLim] for mu in muvals])
    llhdsFunc = np.array([[mu,llhdFromLimits(mu,mu0B,sigma2)/normLim2] for mu in muvals])
    print ( "llhds", llhds )
    print ( "llhdsApp", llhdsApp )
    print ( "llhdsFUnc", llhdsFunc )

    mumax = llhds[np.argmax(llhds[:,1])][0]
    mumaxApp = llhdsApp[np.argmax(llhdsApp[:,1])][0]

    f = plt.figure(figsize=(8,4))
    plt.plot(llhds[:,0],llhds[:,1],label='Full',linewidth=3,color="black")
    plt.plot(llhdsApp[:,0],llhdsApp[:,1],label='From Limits',linestyle='--',linewidth=3)
    plt.plot(llhdsFunc[:,0],llhdsFunc[:,1],label='From Real Limits', color="green", linestyle='-.',linewidth=3)
    plt.legend()
    plt.xlabel(r'$\mu$')
    plt.ylabel(r'$L$')
    plt.text(0,llhdsApp[:,1].min(),'$\mu_{max} = %1.2f$\n$\mu_{max}^{app} = %1.2f$' %(mumax,mumaxApp),fontsize=14)
    plt.title(r'$N_{obs} = %1.2f, N_{bg} = %1.2f \pm %1.2f, \sigma_{UL}^{obs} = %1.2f, \sigma_{UL}^{exp} = %1.2f, \Delta = %1.2f$' 
              %(nobs,nExp,nExpErr,ULobs,ULexp,abs(ULobs-ULexp)/(ULobs+ULexp)),fontsize=14)
    plt.savefig("llhd.png")
    # plt.show()


if __name__ == "__main__":
    nobs = 110
    nExp = 50.0
    nExpErr = 20.0
    nsig = 1 #Keep it at one in order to compare the likelihoods

    run ( nobs, nExp, nExpErr, nsig )
