import numpy as np
import pyreadr
import os
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy import interpolate
from scipy import stats
import pandas as pd

from conversions import jd2date

from paths import inla_results
from paths import readable_data
from nano_plot_flux import apo_jd
from nano_plot_flux import peri_jd

import figure_standards as figstd
axes_size = figstd.set_rcparams_dynamo(mpl.rcParams, num_cols=1, ls='thin')
mpl.rcParams['figure.dpi'] = 600

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
    def __init__(self,
                 datafile,
                 input_csv_readable=readable_data):
        """
        The creator for the InlaResult. Loads the RData output and 
        the csv input. 

        Parameters
        ----------
        datafile : str
            The path to the RData file with the INLA resutl.
        input_csv_readable : str, optional
            The path to the readable CSV data file (input). 
            The default is paths.readable_data. 

        Returns
        -------
        None.

        """
        self.contents = pyreadr.read_r(datafile)
        self.sample_size = len(self.contents["sample_l_bg"])
        self.atts = ["l_bg","l_isd","l_b","v_b_r","e_v","e_r"]
        self.input_df = pd.read_csv(input_csv_readable)

    def inla_function(self):
        """
        Prints the function definition if the RGeneric model used,
        this was obtained as deparse(three_component_model).

        Returns
        -------
        contents["model_definition"] : str
            The original INLA function.

        """
        return self.contents["model_definition"]

    def rate_function(self,
                      careful=True):
        """
        The function that checks what was the function used in R. If it is 
        the usual rate function, that this tranlated is returned. 

        Parameters
        ----------
        careful : bool, optional
            Whether to stop if the R rate function is not the usual function.
            If true, then Exception is raised if they do nt match. If false,
            then Warning is raised if they do not match. The default is True.

        Raises
        ------
        Warning
            If the R function does not match the expectation and the 
            careful is set to False.

        Exception
            If the R function does not match the expectation and the 
            careful is set to True.

        Returns
        -------
        usual_rate : function
            The rate function that is usually used.

        """

        first_row = np.argmax(self.inla_function()[
                                                    "model_definition"
                                                  ].str.find(
                                                        "rate <- function"
                                                            )>0)

        last_row = np.argmax(self.inla_function()[
                                                    "model_definition"
                                                  ].str.find(
                                                        "return(hourly_rate)"
                                                            )>0)

        r_function = np.array(self.inla_function()[first_row:last_row])
        r_function_list = [row[0] for row in r_function]

        usual_r_function_list = [
            '    rate <- function(v_sc_r, v_sc_t, r_sc, v_sc_x, v_sc_y, v_sc_z, ',
            '        l_bg, l_isd, l_b, v_b_r, e_v, e_r, S = 8, v_isd = 26, ',
            '        phi_isd = 259, th_isd = 8, v_b_a = 9, v_earth_a = 29.8) {',
            '        deg2rad <- function(deg) {', '            rad = deg/180 * pi',
            '            return(rad)', '        }', '        L_bg = l_bg',
            '        v_isd_x = -v_isd * sin(deg2rad(90 - th_isd)) * cos(deg2rad(phi_isd))',
            '        v_isd_y = -v_isd * sin(deg2rad(90 - th_isd)) * sin(deg2rad(phi_isd))',
            '        v_isd_z = -v_isd * cos(deg2rad(90 - th_isd))',
            '        L_isd = l_isd * (((v_isd_x - v_sc_x)^2 + (v_isd_y - v_sc_y)^2 + ',
            '            (v_isd_z - v_sc_z)^2)^0.5)/(v_isd)',
            '        ksi = -2 - e_r', '        r_factor = r_sc/1',
            '        v_factor = (((v_sc_r - (v_b_r * (r_factor^ksi)))^2 + ',
            '            (v_sc_t - (v_b_a * (r_factor^(-1))))^2)^0.5)/(((v_b_r)^2 + ',
            '            (v_earth_a - v_b_a)^2)^0.5)',
            '        L_b = l_b * (v_factor)^e_v * (r_factor)^e_r',
            '        hourly_rate = 3600 * S * (L_bg + L_isd + L_b)']

        if r_function_list == usual_r_function_list:
            pass
        elif not careful:
            raise Warning("""the rate function in R is not the usual one,
                          proceed with caution""")
        else:
            raise Exception("""the rate function in R is not the usual one,
                            proceed with caution""")

        def usual_rate(v_sc_r, v_sc_t, r_sc, v_sc_x, v_sc_y, v_sc_z,
                       l_bg, l_isd, l_b, v_b_r, e_v, e_r, S=8, v_isd=26,
                       phi_isd=259, th_isd=8, v_b_a=9, v_earth_a=29.8):
    
            def deg2rad(deg):
                return deg / 180 * np.pi

            L_bg = l_bg
            v_isd_x = -v_isd * (np.sin(deg2rad(90 - th_isd))
                                * np.cos(deg2rad(phi_isd)))

            v_isd_y = -v_isd * (np.sin(deg2rad(90 - th_isd))
                                * np.sin(deg2rad(phi_isd)))

            v_isd_z = -v_isd * (np.cos(deg2rad(90 - th_isd)))

            L_isd = l_isd * np.sqrt((v_isd_x - v_sc_x)**2 +
                                    (v_isd_y - v_sc_y)**2 +
                                    (v_isd_z - v_sc_z)**2) / v_isd
            ksi = -2 - e_r

            r_factor = r_sc / 1

            v_factor = (np.sqrt((v_sc_r - (v_b_r * (r_factor**ksi)))**2 +
                                 (v_sc_t - (v_b_a * (r_factor**(-1))))**2)
                        / np.sqrt((v_b_r**2 + (v_earth_a - v_b_a)**2)))

            L_b = l_b * (v_factor)**e_v * (r_factor)**e_r
            hourly_rate = 3600 * S * (L_bg + L_isd + L_b)

            return hourly_rate

        return usual_rate

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
        """
        The procedure that tshows the priors and the posteriors as per
        the InlaResult. 

        Parameters
        ----------
        atts : list of str, optional
            The list of attributes to show. If None, then 
            self.atts = ["l_bg",
                         "l_isd",
                         "l_b",
                         "v_b_r",
                         "e_v",
                         "e_r"] is used.
            The default is None.
        xrange : list of lists of float, optional
            The x axis limits. The default is [[0,5e-4],
                                               [0,5e-4],
                                               [0,5e-4],
                                               [0,100],
                                               [1,3],
                                               [-2,-1.5]].

        Returns
        -------
        None.

        """

        if atts is None:
            atts = self.atts
        elif isinstance(atts, str):
            atts = [atts]

        fig = plt.figure(figsize=(3,0.666*len(atts)))
        gs = fig.add_gridspec(len(atts), hspace=.6)
        ax = gs.subplots()

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

    def get_errors(self,
                   counts,
                   duty_hours,
                   prob_coverage = 0.9):
        """
        The function to calculate the errorbars for flux 
        assuming Poisson distribution and taking into account
        the number of detection.

        Parameters
        ----------
        counts : array of float
            Counts per day.
        duty_hours : array of float
            Duty hours per day.
        prob_coverage : float, optional
            The coverage of the errobar interval. The default is 0.9.

        Returns
        -------
        err_plusminus_flux : np.array of float
            The errorbars, lower and upper bound, shape (2, n).

        """

        counts_err_minus = -stats.poisson.ppf(0.5-prob_coverage/2,
                                              mu=counts)+counts
        counts_err_plus  = +stats.poisson.ppf(0.5+prob_coverage/2,
                                              mu=counts)-counts
        err_plusminus_flux = np.array([counts_err_minus,
                                       counts_err_plus]) / (duty_hours/(24))

        return err_plusminus_flux

    def plot_fit(self,
                 samples=1000,
                 overplot=None,
                 styles=None,
                 aspect=1.5,
                 apo_jd=apo_jd,
                 peri_jd=peri_jd):
        """
        A plot of daily flux is made with the data from 
        the provided input dataframe and overplotted 
        with the rate based on the posterior.

        Parameters
        ----------
        samples : int
            The number of samples from the posterior that 
            make up the overplotted rate.

        overplot : list of functions: np.array of dt.datetime -> np.array 
                   of float, optional
            A list of functons that will be used to overplot over the data. 
            Default is None, in which case nothing is overplotted.

        styles : list of str, optional
            list of fmt string, such as "g:2" or something, this assigns 
            a style to the overplot lines. Default is None, in which 
            case "b-" is used.

        aspect : float, optional
            The aspect ratio of the resulting plot. The default is 1.2666.

        apo_jd : list of float, optional
            The list of aphelia julian dates of Solar Orbiter.
            The default is nano_plot_flux.apo_jd.

        peri_jd : list of float, optional
            The list of perihelia julian dates of Solar Orbiter.
            The default is nano_plot_flux.peri_jd.
            
        Returns
        -------
        None.
        """

        rate = self.rate_function() # This is correct,
                                    # the rate_function
                                    # returns a function.

        if styles == None and overplot != None:
            styles = ["b-"]*len(overplot)

        # Extract the data from the dataframe.
        jds = np.array(self.input_df["Julian date"])
        dates = np.array([jd2date(jd) for jd in jds])
        counts = np.array(self.input_df["Fluxes [/day]"])
        duty_hours = np.array(self.input_df["Detection time [hours]"])

        # Calculate the errors.
        err_plusminus_flux = self.get_errors(counts, duty_hours)

        # Create the plot.
        fig, ax = plt.subplots(figsize=(3*aspect, 3))
        ax.set_ylabel("Impact rate (duty-cycle corrected) [$day^{-1}$]"
                      , fontsize="medium")
        ax.set_title('CNN Dust Impacts: '+str(np.round(sum(counts)))
                     , fontsize="medium", fontweight="bold")
        ax.tick_params(axis='x',labeltop=False,labelbottom=True)
        ax.tick_params(axis='y',labelleft=True,labelright=False)
        ax.tick_params(axis='y',which="minor",left=True,right=False)
        if overplot!= None:
            for i, line in enumerate(overplot):
                ax.plot(dates,line(dates),styles[i],lw=1,ms=0)
        ax.scatter(dates,counts/duty_hours*24,
                   c="black", s=1,zorder=100)
        ax.errorbar(dates, counts/duty_hours*24,
                    err_plusminus_flux,
                    c="gray", lw=0, elinewidth=0.4,alpha=0.35,zorder=100)
        ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 3, 5,
                                                                7, 9, 11)))
        ax.xaxis.set_minor_locator(mdates.MonthLocator())
        ax.set_xlim(left = min(dates), right = max(dates))
        ax.tick_params(axis='x',which="minor",bottom=True,top=True)
        ax.tick_params(axis='x',labelrotation=60)
        ax.tick_params(labelsize="medium")
        ax.set_ylim(0,1400)
        ax.tick_params(labelleft=True,
                       labelright=True,
                       labelbottom = True,
                       labeltop = False)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        for jd in apo_jd:
            ax.axvline(x=jd2date(jd),color="darkgray",lw=0.7)
        for jd in peri_jd:
            ax.axvline(x=jd2date(jd),color="darkgray",ls="--",lw=0.7)
        ax.text(.82, .96, 'Aphelion', ha='left', va='top', rotation=90,
                color="gray", fontsize="small", transform=ax.transAxes)
        ax.text(.74, .96, 'Perihelion', ha='left', va='top', rotation=90,
                color="gray", fontsize="small", transform=ax.transAxes)

        # Overplot the possibe rates.
        sampled_posterior = self.sample(sample_size=samples)
        v_sc_r = np.array(self.input_df['Radial velocity [km/s]'])
        v_sc_t = np.array(self.input_df['Tangential velocity [km/s]'])
        r_sc   = np.array(self.input_df['Radial distance [au]'])
        v_sc_x = np.array(self.input_df['V_X (HAE) [km/s]'])
        v_sc_y = np.array(self.input_df['V_Y (HAE) [km/s]'])
        v_sc_z = np.array(self.input_df['V_Z (HAE) [km/s]'])

        possible_rates = [rate(v_sc_r,v_sc_t,r_sc,v_sc_x,v_sc_y,v_sc_z,
                               sampled_posterior[0,i],
                               sampled_posterior[1,i],
                               sampled_posterior[2,i],
                               sampled_posterior[3,i],
                               sampled_posterior[4,i],
                               sampled_posterior[5,i])
                          for i in range(samples)]

        for rates in possible_rates:
            ax.plot(dates,24*rates,c="red",alpha=0.02,lw=1)

        fig.tight_layout()
        fig.show()


def list_datafiles(location=inla_results):
    datafiles = glob.glob(inla_results+"*.RData")
    return datafiles


def main(file):
    result = InlaResult(file)
    result.summary_posterior()
    result.plot_prior_posterior()
    result.plot_fit()

    return result

#%%
if __name__ == "__main__":

    files = ["sample_20231124093221_beta_only.RData",
             "sample_20231124093450_beta_bg.RData",
             "sample_20231124093709_beta_isd.RData",
             "sample_20231124094052_beta_bg_isd.RData",
             "sample_20231128114149_inner_fit_only.RData"]

    for file in files:
        main('998_generated\\inla\\'+file)
