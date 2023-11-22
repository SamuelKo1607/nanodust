import numpy as np
import pyreadr
import os
import glob

from paths import inla_results


results = glob.glob(inla_results+"*.RData")


def chew_inla_results(file):

    # odict_keys(['sample_l_bg', 'sample_l_isd', 'sample_l_b', 'sample_v_b_r',
    #             'sample_e_v', 'sample_e_r', 'fx_l_bg', 'fy_l_bg', 'fx_l_isd',
    #             'fy_l_isd', 'fx_l_b', 'fy_l_b', 'fx_v_b_r', 'fy_v_b_r',
    #             'fx_e_v', 'fy_e_v', 'fx_e_r', 'fy_e_r', 'px_l_bg', 'py_l_bg',
    #             'px_l_isd', 'py_l_isd', 'px_l_b', 'py_l_b', 'px_v_b_r',
    #             'py_v_b_r', 'px_e_v', 'py_e_v', 'px_e_r', 'py_e_r',
    #             'model_definition'])


    #TBD: return interpolated functions for priors, posteriors

    samples = pyreadr.read_r(file)


    b1 = np.array(samples["b1"])
    b2 = np.array(samples["b2"])
    c1 = np.array(samples["c1"])
    c2 = np.array(samples["c2"])
    v1 = np.array(samples["v1"])


    pass


chew_inla_results(results[0])