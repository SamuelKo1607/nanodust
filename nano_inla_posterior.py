import numpy as np
import pyreadr
import os
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import interpolate

from paths import inla_results

import figure_standards as figstd
axes_size = figstd.set_rcparams_dynamo(mpl.rcParams, num_cols=1, ls='thin')
mpl.rcParams['figure.dpi']= 600

class InlaResult:
    """
    The container for INLA results coming in the form of .RData file 
    produced by with "inla_fitting.R" containing
    sampels and prior, posterior functions.

    The usual structure of such .RData file:
    
    odict_keys(['sample_l_bg', 'sample_l_isd', 'sample_l_b', 'sample_v_b_r',
                'sample_e_v', 'sample_e_r', 'fx_l_bg', 'fy_l_bg', 'fx_l_isd',
                'fy_l_isd', 'fx_l_b', 'fy_l_b', 'fx_v_b_r', 'fy_v_b_r',
                'fx_e_v', 'fy_e_v', 'fx_e_r', 'fy_e_r', 'px_l_bg', 'py_l_bg',
                'px_l_isd', 'py_l_isd', 'px_l_b', 'py_l_b', 'px_v_b_r',
                'py_v_b_r', 'px_e_v', 'py_e_v', 'px_e_r', 'py_e_r',
                'model_definition'])

    """
    def __init__(self,datafile):
        self.contents = pyreadr.read_r(datafile)
        self.sample_size = len(self.contents["sample_l_bg"])
        self.atts = ["l_bg","l_isd","l_b","v_b_r","e_v","e_r"]

    def inla_function(self):
        """
        Prints the function definition if the RGeneric model used,
        this was obtained as deparse(three_component_model).

        Returns
        -------
        None.

        """
        print(self.contents["model_definition"])

    def summary_posterior(self):
        """
        Prints a short summary (means and variances) of the posteriors
        for all the attributes.

        Returns
        -------
        None.

        """
        for att in ["l_bg","l_isd","l_b","v_b_r","e_v","e_r"]:
            mean = np.mean(self.contents[f"sample_{att}"].to_numpy())
            stdev = np.std(self.contents[f"sample_{att}"].to_numpy())
            print(f"{att}:\t mean = {mean:.3}\t +- {stdev:.2}")

    def sample(self,atts=None,sample_size=None):
        """
        Provides a sample from the posterior of the attribute of choice.

        Parameters
        ----------
        atts : str or list of str or None
            The attrbiute of interest, one or several of self.atts. 
            In None, then the multivariate sample of all 
            the attributes is returned. 
        sample_size : int, optional
            The size of the sample. If None, then the full
            sample as provided by INLA is returned,
            i.e. (sample_size = self.sample_size). 
            The default is None, hence self.sample_size.

        Returns
        -------
        samples : np.ndarray
            If a single attribute is requested, then the sample 
                is 1D array of len sample_size.
            If several (n) attributes are requested, then the sample
                is 2D array of shape (n,sample_size), i.e. the first
                dimension is the attribute order.
            If atts == None, then all the attributes as in self.atts
                are requested and the shape of the 2D array is 
                therefore (len(self.atts),sample_size).

        """
        if sample_size is None:
            sample_size = self.sample_size
            replace = False
        else:
            replace = True

        if atts is None:
            atts = self.atts

        if isinstance(atts, str):
            samples = np.random.choice(
                        self.contents[f"sample_{atts}"].to_numpy().flatten(),
                        size=sample_size,
                        replace=replace)
        else:
            indices = np.random.choice(np.arange(self.sample_size),
                                       size = sample_size,
                                       replace = replace)

            samples = np.zeros((0,sample_size))
            for att in atts:
                row = self.contents[f"sample_{att}"].to_numpy(
                                                     ).flatten()[indices]
                samples = np.vstack((samples,row))

        return samples

    def pdf(self,att,prior=False):
        """
        Returns the PDF, either posterior or prior, for the attribute 
        of choice.

        Parameters
        ----------
        att : str
            The attrbiute of interest, one of self.atts.
        prior : bool, optional
            Whther to return prior. If False, then posterior is returned.
            The default is False.

        Returns
        -------
        interpolated : scipy.interpolate._interpolate.interp1d 
            The pdf function, callable.

        """
        if prior:
            x = self.contents[f"px_{att}"].to_numpy().flatten()
            y = self.contents[f"py_{att}"].to_numpy().flatten()
        else:
            x = self.contents[f"fx_{att}"].to_numpy().flatten()
            y = self.contents[f"fy_{att}"].to_numpy().flatten()
        x = x[y>max(y)/10**4]
        y = y[y>max(y)/10**4]
        interpolated = interpolate.interp1d(x, y,
                                            kind=3,
                                            bounds_error=False,
                                            fill_value=0)
        return interpolated

    def plot_prior_posterior(self,
                             atts=None,
                             xrange = [[0,5e-4],
                                       [0,5e-4],
                                       [0,5e-4],
                                       [0,100],
                                       [1,3],
                                       [-2,-1.5]]):

        if atts is None:
            atts = self.atts
        elif isinstance(atts, str):
            atts = [atts]

        fig,ax = plt.subplots(ncols=1,nrows=len(atts),
                              figsize=(3,0.666*len(atts)))

        for i,att in enumerate(atts):
            x = np.linspace(xrange[i][0],xrange[i][1],num=1000)
            prior = self.pdf(att,prior=True)(x)
            posterior = self.pdf(att,prior=False)(x)
            ax[i].plot(x,prior/max(prior),label="prior")
            ax[i].plot(x,posterior/max(posterior),label="posterior")
            ax[i].set_xlim(xrange[i])
            ax[i].set_ylim(bottom=0)
            ax[i].set_ylabel(att)
        ax[0].legend(fontsize="x-small")

        fig.show()











def main():
    results = glob.glob(inla_results+"*.RData")

    result = InlaResult(results[0])

    result.summary_posterior()
    result.plot_prior_posterior()

    x=np.linspace(0,5e-4,num=1000)
    plt.plot(x,result.pdf("l_b",prior=0)(x))
    plt.show()



main()