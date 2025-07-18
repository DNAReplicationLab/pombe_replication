import numpy as np
import pandas as pd
import itertools
import uuid
import random
import re
import os

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # stop unnecessary TensorFlow messages

import tensorflow.compat.v2 as tf
import tensorflow_probability as tfp
from scipy.special import expit
from scipy.stats import bernoulli

tf.enable_v2_behavior()
tfd = tfp.distributions


class BrdUVsPosModel:
    """Bare bones theoretical mean BrdU vs ref genome position on a fork.
    We use a linear model y = p1 * x for illustration.

    Attributes:
        params (list): List of just one parameter
        win_size (int): size of BrdU windows
        noiseStrength (float): noise amplitude
    """

    def __init__(self, win_size):
        """ Bare-bones initialization function. Simple linear model y = p1 * x.

        Args:
            win_size (int): size of BrdU windows
        """
        self.params = [1]
        self.win_size = win_size
        self.noiseStrength = 0

    def get_model(self, x):
        """Return model values at given locations"""

        return [self.params[0] * k for k in x]

    def get_model_grad_wrt_params(self, x):
        """Return list w entries = model grad at given locations w.r.t. 
           corresponding parameter
        """
        return [x]

    def __str__(self):
        """ Returns a string representation of the class """
        return f"win_size = {self.win_size}; params = {self.params}"

    def set_noise_strength(self, noise_strength):
        """ Set noise strength """
        self.noiseStrength = noise_strength

    def get_model_with_noise(self, x):
        """Return model values at given locations with noise"""

        return [m[0] + m[1] for m in zip(self.get_model(x),
                                         np.random.normal(0, self.noiseStrength, len(x)))]

    def generate_data(self, fork_len_in_wins, set_noise=True,
                      ori_loc_in_wins=0,
                      pause_site_in_bases=np.nan, pause_time_in_wins=0):
        """ Generate windowed BrdU data from the model.

        NOTE: (1) x coordinates in comments below refer to the x coordinate
        used by the model, and not any real x coordinate like genomic
        coordinates. 
        (2) It may not be the case that the model begins at x = 0. Meaning
        of x values depend on the model definition.

        Args:
            fork_len_in_wins (int): fork length in multiples of window size
            set_noise (bool): default True, use noisy version of model or not
            ori_loc_in_wins (int): origin location at x = num
            pause_site_in_bases (int): def. nan, fork pause at x = origin location + num
            pause_time_in_wins (int): def. 0 (no pause), fork pause time, in win num

        Returns:
            a tuple of (window start, window end, mean brdu density)
        """

        # set up output variables
        win_start = []
        win_end = []
        mean_brdu = []

        # index of current window
        curr_win_indx = 0

        # get model values
        if set_noise:
            brdu_model_vals = self.get_model_with_noise(
                [ori_loc_in_wins + k for k in
                 range(fork_len_in_wins + pause_time_in_wins)])
        else:
            brdu_model_vals = self.get_model(
                [ori_loc_in_wins + k for k in
                 range(fork_len_in_wins + pause_time_in_wins)])

        for k in range(fork_len_in_wins):

            # record window positions
            win_start.append(k * self.win_size + 1)
            win_end.append((k + 1) * self.win_size)

            if (pause_time_in_wins > 0 and
                    win_start[-1] <= pause_site_in_bases <= win_end[-1]):

                # if pause site is within current window,
                # model interpolates between two different
                # brdU densities

                frac_win = (pause_site_in_bases - win_start[-1]) / self.win_size
                mean_brdu.append(frac_win * brdu_model_vals[curr_win_indx] +
                                 (1 - frac_win) * brdu_model_vals[curr_win_indx +
                                                                  pause_time_in_wins])
                curr_win_indx += pause_time_in_wins + 1

            else:

                # if no pause site in current window,
                # then just record mean brdu values from the model
                mean_brdu.append(brdu_model_vals[curr_win_indx])
                curr_win_indx += 1

        return win_start, win_end, mean_brdu


class BrdUVsPosModelSimpleLinear(BrdUVsPosModel):
    """Linear mean BrdU vs window coordinate. *for testing only*

    Note: a linear brdu vs time model makes no sense as brdu densities
        cannot go negative! Only use this model for testing

    Attributes:
        params (list): List of two model parameters
        win_size (int): size of BrdU windows
    """

    def __init__(self, win_size, params):
        """Initialize linear function
        f(x) = p0 + p1*x
        
        Args:
            win_size (int): size of BrdU windows
            params (list): two params, all float
        """

        super().__init__(win_size)

        if len(params) != 2:
            raise ValueError('Need 2 parameters!')

        self.params = params

    def get_model(self, x):
        """Return model values at given locations"""

        return [self.params[0] + self.params[1] * k for k in x]

    def get_model_grad_wrt_params(self, x):
        """Return list w entries = model grad at given locations w.r.t. 
           corresponding parameter
        """

        return [[1 for _ in x],
                [_ for _ in x]]


class BrdUVsPosModelSimpleSigmoidal(BrdUVsPosModel):
    """Sigmoidal mean BrdU vs window coordinate

    Attributes:
        params (list): List of one model parameter
        win_size (int): size of BrdU windows
    """

    def __init__(self, win_size, params=None):
        """see base class"""

        super().__init__(win_size)

        if params is None:
            params = [1]

        self.params = params

    def get_model(self, x):
        """Return model values at given locations"""

        return [1 / (1 + np.exp(k / self.params[0]))
                for k in x]

    def get_model_grad_wrt_params(self, x):
        """Return list w entries = model grad at given locations w.r.t. 
           corresponding parameter
        """

        step1 = [np.exp(k / self.params[0]) for k in x]
        return [[m[0] / self.params[0] ** 2 * m[1] / (1 + m[1]) ** 2
                 for m in zip(x, step1)]]


class BrdUVsPosModelSigmoidal(BrdUVsPosModel):
    """Sigmoidal mean BrdU vs window coordinate

    Attributes:
        params (list): List of model parameters
        win_size (int): size of BrdU windows
    """

    def __init__(self, win_size, params):
        """Initialize sigmoidal function
        f(x) = p0 + p1/(1+e^(x/p2))
        
        Args:
            win_size (int): size of BrdU windows
            params (list): four params, all float
        """

        super().__init__(win_size)

        if len(params) != 3:
            raise ValueError('Need 3 parameters!')

        self.params = params

    def get_model(self, x):
        """Return model values at given locations x
        
        Args:
            x (list): list of floats of window locations
        """

        return [self.params[0] +
                self.params[1] / (1 + np.exp(k / self.params[2]))
                for k in x]

    def get_model_grad_wrt_params(self, x):
        """Return list w entries = model grad at given locations w.r.t. 
           corresponding parameter
        """

        nx = len(x)
        exp1 = [np.exp(k / self.params[2]) for k in x]
        exp2 = [self.params[1] / self.params[2] * y / (1 + y) ** 2 for y in exp1]
        return [[1 for _ in range(nx)],
                [1 / (1 + y) for y in exp1],
                [m[0] * m[1] / self.params[2]
                 for m in zip(x, exp2)]]


class BrdUVsPosModelSigmoidalClipped(BrdUVsPosModelSigmoidal):
    """Sigmoidal mean BrdU vs window coordinate clipped to lie b/w 0 and 1

    Attributes:
        params (list): List of model parameters
        win_size (int): size of BrdU windows
    """

    def __init__(self, win_size, params):
        """Initialize sigmoidal function.
        f(x) = p0 + p1/(1+e^(x/p2)). We require
        f to be a decreasing function over x (p0 + p1 > p0)
        and params to be such that 0 <= f(x) <= 1.
        
        Args:
            win_size (int): size of BrdU windows
            params (list): four params, all float
        """

        super().__init__(win_size, params)

        if (params[0] + params[1] > 1 or
                params[0] + params[1] < 0 or
                params[0] > 1 or
                params[0] < 0 or
                params[2] <= 0):
            raise ValueError("Invalid parameters!")

        if params[1] < 0:
            raise NotImplementedError("Do not use this parameter set!")

    def get_model_with_noise(self, x):
        """Return model w noise & clipped at given locations x
        
        Args:
            x (list): list of floats of window locations
        """

        # only let through values b/w 0 and 1, else
        # convert negative vals to 0 and vals > 1 to 1.

        return list(np.clip(super().get_model_with_noise(x), a_min=0, a_max=1))

    def get_model_grad_wrt_params(self, x):
        """Do not use as gradients of a clipped fn are not defined"""

        raise NotImplementedError("Do not use this function!")

    def invert_model(self, y, tol=None):
        """Invert the model

        As a sigmoid has two almost flat sections, it cannot be
        inverted very well except in the middle portion where there
        is a strong decrease. Thus, we invert only from y = p0 + p1*tol[1]
        to y = p0 + p1*tol[0]. tol[0], tol[1] are elements of the tol list.
        For definitions of p0, p1 see __init__.
        
        Args:
            y (list of floats): list of model values
            tol (list of floats): only invert over tol interval (see desc). df is
                [0.05, 0.95]

        Returns:
            list of floats, the corresponding x values
        
        """

        if tol is None:
            tol = [0.05, 0.95]

        if (tol[0] < 0 or tol[1] < 0 or tol[0] > 1 or tol[1] > 1
                or tol[0] >= tol[1]):
            raise ValueError("Invalid tolerances specified!")

        p0 = self.params[0]
        p1 = self.params[1]
        p2 = self.params[2]

        def fltTol(x):
            return (np.nan, x)[p0 + p1 * tol[1] >= x >= p0 + p1 * tol[0]]

        return [p2 * np.log(p1 / (f - p0) - 1) for f in
                map(fltTol, y)]


class BernoulliNot0to1(tfp.glm.ExponentialFamily):
    """Sigmoid that doesn't go from 0 to 1."""

    def __init__(self, p_low: float = 0, p_high: float = 1):
        self._pLow = p_low
        self._pHigh = p_high
        self._pDelta = p_high - p_low
        super(BernoulliNot0to1, self).__init__(name=None)

    def _call(self, r):
        sig = tf.math.sigmoid(r)
        mean = self._pLow + self._pDelta * sig
        variance = mean * (1 - mean)
        grad_mean = self._pDelta * sig * (1 - sig)
        return mean, variance, grad_mean

    def _as_distribution(self, r):
        return tfd.Bernoulli(probs=self._pLow + self._pDelta * tf.math.sigmoid(r))


class ErrorInfo:
    """ Store errors counts and msgs that are alphanumeric strs """

    def __init__(self, error_list):
        """ Set some basic parameters 
        
        Args:
            error_list (list of str): possible error names, each nm alphanumeric
        """

        self._N = dict()  # a count associated with each error

        for k in error_list:
            if not k.isalnum():
                raise ValueError('Only alphanumeric strings!')
            self._N[k] = 0

    def msg(self, msg=""):
        """ Set/get the error message """
        if msg == "":
            return self.tolist()

        if msg in self._N:
            self._N[msg] += 1
        else:
            raise ValueError('Only alphanumeric strings!')

    def incr(self, msg):
        """ Increment the counter """
        self.msg(msg)

    def count(self):
        """ report the total count """
        return sum(self._N.values())

    def isSet(self):
        """ Return T if error class under use """
        return any(self._N.values())

    def tolist(self):
        """ Return a list of errors recorded """
        return {k: str(v) for k, v in self._N.items()}


def fit_model_grad_descent(dataset_brdu_pos, offset, model, rel_param_tol,
                           win_pos, step_size):
    """ Fits given model w grad descent to mn brdU density vs ref. gnm. pos data

    Args:
        dataset_brdu_pos (list): a list of lists of mean brdU amounts from successive
            windows of equal size.
        offset (list): list of offsets
        model (BrdUVsPosModel): model initialized to some param values
        rel_param_tol (float): stop iterations when param changes are this small
            i.e. when abs(param_n - param_(n-1))/param_(n-1) < relParamTol.
            Assumption here is that all params converge to non-zero values.
        win_pos (list): only consider model in window numbers x = win_pos[0]
            to x = win_pos[1] for fitting
        step_size (float): param that sets size of each parameter update step

    Returns:
        list w best-fit model (BrdUVsPosModel), offset (list), lsq err (float)
    """

    num_params = len(model.params)
    delta_params = [rel_param_tol * k for k in model.params]
    windows = range(win_pos[0], win_pos[1] + 1)

    while (max(np.abs([delta_params[k] / model.params[k] for k in range(num_params)]))
           >= rel_param_tol):
        model_curve = model.get_model(windows)
        grad_model = model.get_model_grad_wrt_params(windows)
        goodness_of_fit = sum(calc_sq_diff_w_offset(model_curve, dataset_brdu_pos,
                                                    offset))
        grad_goodness_of_fit = [2 * sum(calc_diff_dot_f_w_offset(model_curve,
                                                                 x, dataset_brdu_pos, offset)) for x in grad_model]

        delta_params = [step_size * grad_goodness_of_fit[k] for k in range(num_params)]
        model.params = [model.params[k] + delta_params[k] for k in range(num_params)]
        offset = calc_offset(model.get_model(windows), dataset_brdu_pos)

    return [model, offset, goodness_of_fit]


def lsq_error_model_brute_force(dataset_brdu_pos, model, cand_params, win_pos,
                                bin=False, lowest_n_per_bin=3):
    """ Fits given model w brute force to mn brdU density v ref. gnm. pos data,
    if multiple global minima are seen, only one of them is reported, and
    logic for which one is reported is not meant to be reliable.

    Args:
        dataset_brdu_pos (list): a list of lists of mean brdU amounts from successive
            windows of equal size.
        model (BrdUVsPosModel): model whose best-fit param values are sought
        cand_params (list): list of lists, each item a point in Ndim param spce
            where N = number of params in the model.
        win_pos (list): only consider model in window numbers x = win_pos[0]
            to x = win_pos[1] for fitting
        bin (binary): default False. If True, then bin ip data along x
        lowest_n_per_bin (int): default 3. mark as nan if fewer pts in bin.

    Returns:
        list w best-fit model (BrdUVsPosModel), best-fit offset (list),
            best-fit model's goodness of fit (float),
            lsq errors of each set of params (list of floats)
    """

    windows = range(win_pos[0], win_pos[1] + 1)

    gof_list = []
    gof_min = np.inf
    best_fit_params = []
    best_fit_offset = []

    for k in cand_params:

        # set params and get model curve
        model.params = k
        model_curve = model.get_model(windows)

        # get offset and calculate goodness of fit
        offset = calc_offset(model_curve, dataset_brdu_pos)

        if not bin:
            gof = sum(calc_sq_diff_w_offset(model_curve, dataset_brdu_pos,
                                            offset))
            gof_list.append(gof)
        else:
            brdu_pos_mean, brdu_pos_err = bin_data_w_offset(dataset_brdu_pos, offset,
                                                            lowest_n_per_bin)
            gof = calc_sq_diff_w_error(model_curve, brdu_pos_mean,
                                       brdu_pos_err)
            gof_list.append(gof)

        # save if gof is lower than the min thus far
        if gof < gof_min:
            gof_min = gof
            best_fit_params = [l for l in model.params]
            best_fit_offset = [l for l in offset]

    best_fit_model = model
    best_fit_model.params = best_fit_params
    return [best_fit_model, best_fit_offset, gof_min, gof_list]


def calc_offset(cndt_brdu_pos, dataset_brdu_pos):
    """Align mean brdU density vs ref. genome pos curves to candidate curve.

    Example:
        >>> cndt_brdu_pos1 = [0.9,0.8,0.7,0.6]
        >>> dataset_brdu_pos1 = [[0.9,0.8],[0.82,0.71,0.62],[0.8,np.nan,0.6]]
        >>> print(calc_offset(cndt_brdu_pos1, dataset_brdu_pos1))
        [0,1,1]

    NOTE: we assume that for both the candidate curve and the dataset curves,
        mean brdU amounts were calculated using non-overlapping windows of
        same size w no gaps b/w them. While candidate curve is likely a
        theoretical curve, the dataset curves come from real data and there
        could be a distribution in window sizes and some gaps b/w windows.
        Function designed to ignore such variations.

    NOTE: If there are multiple alignments that are equally good, only one
        offset is reported. Program is not designed to be reliable in which
        offset is reported.

    Args:
        cndt_brdu_pos (list): mean brdU amts msr from (theoretical) successive
            windows of equal size.
        dataset_brdu_pos (list): a list of lists of mean brdU amts from successive
            windows of equal size.

    Returns:
        list of offsets
    """

    # pretty straightforward -  minimize least squares difference
    # between reference curve and data curve. only free parameter
    # is the horizontal offset of the data curve.

    offset_list = []

    N = len(cndt_brdu_pos)

    for curve in dataset_brdu_pos:

        min_diff = np.inf
        curr_off = 0

        M = len(curve)

        for j in range(N - M + 1):

            curr_diff = np.nansum([(curve[i] - cndt_brdu_pos[i + j]) ** 2
                                   for i in range(M)])

            if curr_diff < min_diff:
                curr_off = j
                min_diff = curr_diff

        offset_list.append(curr_off)

    return offset_list


def bin_data_w_offset(dataset_brdu_pos, offset, lowestNperBin):
    """Assign each data point to a window and get mean, s.d. per window.

    Example:
        >>> dataset1 = [[0.9,0.8,0.7,0.6,0.5], [1.0,0.9,0.8,0.7,0.6], [0.8,0.7,0.6,0.5,0.4], [0.2,0.1], [0.3,0.2]]
        >>> offset1  = [0,0,0,1,3]
        >>> # this means first,second,third curve x vals are = 0,...,4
        >>> # fourth x vals x=1,2 fifth x vals = 3,4
        >>> t1 = bin_data_w_offset(dataset1, offset1, 4)
        >>> for m in range(1,len(t1[1])):
        ...     t1[1][m] = round(t1[1][m],2)
        ...     t1[0][m] = round(t1[0][m],3)
        >>> print(t1)
        ([nan, 0.65, 0.55, 0.525, 0.425], [nan, 0.31, 0.31, 0.17, 0.17])

    Args:
        dataset_brdu_pos (list): list of lists of windowed brdU data
        offset (list): horizontal offset of each brdU curve frm x = 0
        lowestNperBin (int): if fewer data pts per bin, then assign NaNs.

    Returns:
        tuple of two lists, mean & s.d. per bin
    """

    windowed_data = dict()
    windowed_data[0] = list()
    max_indx = 0

    op_mean = []
    op_sd = []

    for (curve, currOff) in zip(dataset_brdu_pos, offset):

        # extend the windowed_data dictionary if need be
        last_indx = currOff + len(curve) - 1

        if last_indx > max_indx:
            for k in range(max_indx + 1, last_indx + 1):
                windowed_data[k] = list()
            max_indx = last_indx

        # collect data into the appropriate bin
        for k in zip(curve, range(currOff, last_indx + 1)):
            windowed_data[k[1]].append(k[0])

    # perform binning if requisite min num of pts are available
    for key, value in windowed_data.items():
        if len(value) >= lowestNperBin:
            op_mean.append(np.nanmean(value))
            op_sd.append(np.nanstd(value, ddof=1))
        else:
            op_mean.append(np.nan)
            op_sd.append(np.nan)

    return op_mean, op_sd


def calc_sq_diff_w_error(candidate_brdu_v_pos, bindataset_brdu_pos, errdataset_brdu_pos):
    """Sq diff calcn b/w candidate curve and data w mean and sd per x coord.

    Example: 
        >>> import numpy
        >>> model = [0.9,0.8,0.7,0.6]
        >>> binDt = [[0.9,0.8],[np.nan,0.82,0.71,0.62]]
        >>> errDt  = [[0.1,0.2],[np.nan,0.02,0.02,0.02]]
        >>> print(numpy.round(calc_sq_diff_w_error(model,binDt[0],errDt[0])))
        >>> 0.0
        >>> print(numpy.round(calc_sq_diff_w_error(model,binDt[1],errDt[1]),2))
        >>> 2.25

    Args:
        candidate_brdu_v_pos (list): mean brdU amts msr from (theoretical) successive
            windows of equal size.
        bindataset_brdu_pos (list): each entry one mean brdU amt at that window
        errdataset_brdu_pos (list): s.d. of brdU amt at that window

    Returns:
        sq differences w errors as weights
    """

    return np.nansum([((bindataset_brdu_pos[k] - candidate_brdu_v_pos[k]) / errdataset_brdu_pos[k]) ** 2
                      for k in range(len(bindataset_brdu_pos))])


def calc_sq_diff_w_offset(candidate_brdu_v_pos, dataset_brdu_pos, offset):
    """Sq diff calcn b/w candidate curve and data given an offset per curve.

    Example:
        >>> candidate_brdu_v_pos_test = [0.9,0.8,0.7,0.6]
        >>> dataset_brdu_pos_test = [[0.9,0.8],[0.82,0.71,0.62],[0.8,np.nan,0.6]]
        >>> offset_test = [0,1,1]
        >>> print(list(map(lambda x: round(x,4), calc_sq_diff_w_offset(candidate_brdu_v_pos_test, dataset_brdu_pos_test,
        ...              offset_test))))
        [0.0, 0.0009, 0.0]

    Args:
        candidate_brdu_v_pos (list): mean brdU amts msr from (theoretical) successive
            windows of equal size.
        dataset_brdu_pos (list): each entry a list of mean brdU amts from successive
            windows of equal size.
        offset (list): horizontal offset of each data curve

    Returns:
        list of sq differences
    """

    sq_diff_list = []
    model_curve_length = len(candidate_brdu_v_pos)

    for (curve, currOff) in zip(dataset_brdu_pos, offset):
        n = len(curve)

        if n + currOff > model_curve_length:
            raise IndexError(f"Current offset {currOff}, data length {n}, "
                             f"but model curve length is {model_curve_length} -> insufficient!")

        curr_diff = np.nansum([(curve[i] - candidate_brdu_v_pos[i + currOff]) ** 2
                               for i in range(n)])

        sq_diff_list.append(curr_diff)

    return sq_diff_list


def calc_diff_dot_f_w_offset(candidate_brdu_v_pos, fn_brdu_pos, dataset_brdu_pos, offset):
    """Dot product of a given function and diff b/w candidate curve and
    data w an offset per curve. Function useful when calculating gradients
    of model fits w.r.t. parameters.

    Example:
        >>> model = [0.9,0.8,0.7,0.6]
        >>> data = [[0.9,0.8],[0.82,0.71,0.62],[0.8,np.nan,0.6]]
        >>> fn = [1,2,3,4]
        >>> data_offset = [0,1,1]
        >>> print(list(map(lambda x: round(x,2), calc_diff_dot_f_w_offset(model, fn, data, data_offset))))
        [0.0, 0.15, 0.0]

    Args:
        candidate_brdu_v_pos (list): mean brdU amts msr from (theoretical) successive
            windows of equal size.
        fn_brdu_pos (list): some function from (theoretical) successive
            windows of equal size.
        dataset_brdu_pos (list): each entry a list of mean brdU amts from successive
            windows of equal size.
        offset (list): horizontal offset of each data curve

    Returns:
        list of dot products
    """

    diff_list = []

    model_curve_length = len(candidate_brdu_v_pos)
    assert(model_curve_length == len(fn_brdu_pos))

    for (curve, currOff) in zip(dataset_brdu_pos, offset):
        M = len(curve)

        if M + currOff > model_curve_length:
            raise IndexError(f"Current offset {currOff}, data length {M}, "
                             f"but model curve length is {model_curve_length} -> insufficient!")

        curr_diff = np.nansum([(curve[i] - candidate_brdu_v_pos[i + currOff]) *
                               fn_brdu_pos[i + currOff] for i in range(M)])

        diff_list.append(curr_diff)

    return diff_list


def calc_offset_bw_models(model1, model2, x_vals, cnd_offset_list):
    """ Align one model curve to another by displacement along x-axis.

    Args:
        model1 (BrdUVsPosModel): model of BrdU data
        model2 (BrdUVsPosModel): model of BrdU data
        x_vals (list of floats): x coordinates for model curve
        cnd_offset_list (list of floats): candidate offset values

    Returns:
        best-fit offset of modl2 w.r.t modl1 i.e. m2(x + offset) fits m1(x) best
    """

    m1 = model1.get_model(x_vals)
    lowest_diff_bw_models = len(x_vals)
    best_offset = None

    for k in cnd_offset_list:
        m2 = model2.get_model([b + k for b in x_vals])
        diff_bw_models = sum([(n[1] - n[0]) ** 2 for n in zip(m1, m2)])
        if diff_bw_models < lowest_diff_bw_models:
            lowest_diff_bw_models = diff_bw_models
            best_offset = k

    if not best_offset:
        raise NotImplementedError('Unknown error!')

    return best_offset


def process_for_model_fit(df_val, window, model, param_list, x_lim, win_len_tol):
    """ Process input data frame and fit to model
    
    Args:
        df_val (pandas df): must have cols detectIndex, mean_brdU, start, end
        window (int): window size in bases on reference
        model (BrdUVsPosModel): model of BrdU data
        param_list (list of lists): each list is a point in param space
        x_lim (list of two ints): x coordinate limits (both incl.) for model
        win_len_tol (list 2 nums asc. sorted): reject windows w lengths outside limits

    Returns:
        tuple of pandas df, dict w fit params plus other info
    """

    # part 1
    # ======
    # * extract the brdU traces for model fitting
    # * assign window indices 0,1,2,.. to the data
    # * also exclude windows due to undesirable lengths

    window_indices = []
    window_rejected = []

    dataset_brdu_pos = []
    nforks = 0
    fork_table_row_count = []

    # check each window is within size tolerance before
    # adding it to the data for fitting
    for name, group in df_val.groupby('detectIndex', sort=False):

        N = len(group['mean_brdU'])
        brdu_data = []

        first_window_mid = (group['end'].iloc[0] + group['start'].iloc[0]) / 2

        for cnt in range(N):

            # get the current window's index (i.e. first window = 0,
            #     2nd window = 1 and so on)
            curr_window_mid = (group['end'].iloc[cnt] +
                               group['start'].iloc[cnt]) / 2

            curr_window_num = (curr_window_mid - first_window_mid) / window

            if curr_window_num < 0:
                curr_window_num = -curr_window_num

            curr_window_num = int(np.round(curr_window_num))

            # we want window index to have incremented by one. so:
            # (1) if window numbers haven't incremented for whatever reason,
            #     reject the window
            # (2) if they've incremented by one, all's well.
            # (3) if they've incremented by more than one, we need to add
            #     nans for the missing windows, and reject the missing windows.
            if cnt == 0:
                window_indices.append(curr_window_num)
            elif curr_window_num <= window_indices[-1]:
                window_rejected.append(1)
                brdu_data.append(np.nan)
                window_indices.append(window_indices[-1] + 1)
                continue
            else:
                for winIndx in range(window_indices[-1], curr_window_num - 1):
                    brdu_data.append(np.nan)
                    # note: as these windows do not exist in the input data,
                    # we are not recording their acceptance status
                    # or the window index.
                    # window_rejected.append(1)
                    # window_indices.append(winIndx)
                window_indices.append(curr_window_num)

            # get window size and see if it falls within our tolerance
            # accept or reject suitably.
            win_size = group['end'].iloc[cnt] - group['start'].iloc[cnt]
            if win_size < win_len_tol[0] or win_size > win_len_tol[1]:
                brdu_data.append(np.nan)
                window_rejected.append(1)
            else:
                brdu_data.append(group['mean_brdU'].iloc[cnt])
                window_rejected.append(0)

        dataset_brdu_pos.append(brdu_data)
        fork_table_row_count.append(N)
        nforks += 1

    # part 2
    # ======
    # add new columns to input data frame,
    # obtain the best-fit model,
    # write model data to dataframe,
    # write binned data to dataframe

    df_val['windowIndex'] = window_indices
    df_val['windowRejected'] = window_rejected
    df_val['gof'] = np.nan
    df_val['deltaLower'] = np.nan
    df_val['deltaUpper'] = np.nan

    # fit model
    op_list = lsq_error_model_brute_force(dataset_brdu_pos, model,
                                          param_list, x_lim)

    # use best-fit offsets to move each fork horizontally
    # to the correct location to align with the model prediction
    cnt = 0
    for k in range(nforks):
        for m in range(fork_table_row_count[k]):
            df_val.at[cnt, 'windowIndex'] += op_list[1][k]
            if df_val.at[cnt, 'windowRejected'] == 0:
                df_val.at[cnt, 'gof'] = (df_val.at[cnt, 'mean_brdU'] -
                                         op_list[0].get_model([df_val.at[cnt, 'windowIndex'] +
                                                               x_lim[0]])[0]) ** 2
            else:
                df_val.at[cnt, 'gof'] = np.nan
            cnt += 1

    # best fit goodness of fit
    gof_model_min = op_list[2]

    # best fit params
    best_fit_params = op_list[0].params

    # to obtain a distribution of model values for every x,
    # we need to set a ceiling for goodness of fit and see
    # which parameter values lead to models that lie between
    # minimum goodness of fit and the ceiling

    # 95% conf int for goodness of fit (w=1.96 sigma is approx 95% of gaussian)
    Ndata = sum(1 - k for k in list(df_val.windowRejected))
    Nparams = len(op_list[0].params)
    gof_limit = gof_model_min * (1 + 1.96 * 1.96 / (Ndata - Nparams))

    # get params that lie in the 95% conf int
    ciParamList = [param for param, gofParam in zip(param_list, op_list[3])
                   if gofParam <= gof_limit]

    # add curves of the best-fit model value at each x,
    # and the lowest and highest model values at each x which
    # form the limits of the 95% confidence interval
    new_rows = []
    x_coords = list(range(x_lim[0], x_lim[1] + 1))
    window_coords = list(range(x_lim[1] - x_lim[0] + 1))

    best_fit_model_data = op_list[0].get_model(x_coords)
    best_fit_model_delta_lower_data = []
    best_fit_model_delta_upper_data = []
    other_model_data = []

    for k in ciParamList:
        model.params = k
        other_model_data.append(model.get_model(x_coords))

    for k in range(len(best_fit_model_data)):
        n1 = min(other_model_data[m][k] for m in range(len(ciParamList)))
        n2 = max(other_model_data[m][k] for m in range(len(ciParamList)))
        best_fit_model_delta_lower_data.append(best_fit_model_data[k] - n1)
        best_fit_model_delta_upper_data.append(n2 - best_fit_model_data[k])

    for k in zip(window_coords, best_fit_model_data,
                 best_fit_model_delta_lower_data, best_fit_model_delta_upper_data):
        new_rows.append(['bestFitModel', k[1], k[0] * window + 1,
                         (k[0] + 1) * window, k[0], 0, gof_model_min,
                         k[2], k[3]])

    # bin the data
    mnBrdUList, sdBrdUList = bin_data_w_offset(dataset_brdu_pos, op_list[1], 3)

    for modl, mn, sd, cnt in zip(best_fit_model_data, mnBrdUList, sdBrdUList,
                                 itertools.count()):
        new_rows.append(['aggData', mn, cnt * window + 1,
                         (cnt + 1) * window, cnt, 0, ((modl - mn) / sd) ** 2, sd, sd])

    df_new = pd.DataFrame(new_rows, columns=['detectIndex', 'mean_brdU',
                                             'start', 'end', 'windowIndex', 'windowRejected', 'gof',
                                             'deltaLower', 'deltaUpper'])

    # stitch data frames together and return
    df_val = pd.concat([df_val, df_new], sort=False)

    return df_val, {'bestFitParamList': best_fit_params,
                    'ciParamList': ciParamList}


def mark_outliers95_using_model(gof):
    """ Using model, flag points that lie outside 95% CI

    Args:
        gof (list): ith element = (data_i - model_prediction)^2 or np.nan

    Returns:
        list w 1s where element not within CI or 0s otherwise

    """

    # calculate variance in the model
    sd_sq = np.nanmean(gof)

    # now, we calculate that x^2 where exp(-x^2/2*sigma^2) = 5/100
    # or x^2 = 2*sigma^2*log(20)

    x_crit_sq = 2 * sd_sq * np.log(20)

    return [1 if k > x_crit_sq else 0 for k in gof]


def devel_fix_gof(f_name):
    """ *dev fn, do not use!* Fix gof errors

    Args:
        f_name (str): fit file with wrong gof columns

    Returns:
        None
    """

    col_names = ["dummyLabel", "forkid", "mean_brdU", "winStart", "winEnd",
                 "modelAlignWinIndx", "winRejected", "goodnessOfFit", "limLow", "limUpr", "label"]

    df = pd.read_csv(f_name, sep=" ", header=None,
                     names=col_names)

    brdUList = [k for k in df[df.forkid == "bestFitModel"]["mean_brdU"]]

    for k in range(len(df)):
        if df.at[k, "dummyLabel"] == "winDetect":
            df.at[k, "goodnessOfFit"] = (df.at[k, "mean_brdU"] -
                                         brdUList[df.at[k, "modelAlignWinIndx"]]) ** 2

    with open(f_name, 'w') as fp:
        print(df.to_csv(index=False, header=False, float_format='%.10f',
                        sep=" ", na_rep="NA"), file=fp, end="")


def make_dummy_windowed_fork_data(model, contig, originTime,
                                  model_x_for_t0, pauseSite, pauseTime, forkStop,
                                  forkSpeed):
    """ Make dummy, right-moving fork data with given origin, pause site

    At given contig, there's a dummy origin at the first base.
    It fires at a given time, where t = 0 is the start of S phase
    and 1 rightward moving fork at uniform speed emerges.
    Forks pause at the pauseSite for the given time (which can be 0),
    and continue till they are terminated at forkStop.

    Args:
        model (BrdUVsPosModel): model for B/B+T v reference coordinate
        contig (str): contig where forks operate
        originTime  (float): start time
        model_x_for_t0 (float): t = 0 corresponds to x = num on model coordinate
        pauseSite (int): fork pauses at given site on contig
        pauseTime (float): pause time
        forkStop (int): forks stop at given site on contig, has to be a
            multiple of winSize specified in model
        forkSpeed (float): in b/min, speed of forks

    Returns:
        Pandas dframe w cols index, label, detectIndex, mean_brdu, start, end
        Example:
        ...
        0 winDetect 64756d6d-79fe-4094-9fae-d60a178dcc8f_chrVI_1_20000_fwd_R_1_20000 0.747 194336 195238
        1 winDetect 64756d6d-79fe-4094-9fae-d60a178dcc8f_chrVI_1_20000_fwd_R_1_20000 0.690 195249 196425
        ...
        NOTE: the string "64756d6d79" is "dummy" if considered as a string
        of five ascii characters represented using hexadecimal numbers.
        The other hex characters "fe-4094-9fae-d60a178dcc8f" are generated
        at a random.
    """

    # ascertain fork length in multiples of window size of model
    if forkStop % model.win_size != 0:
        raise ValueError('Invalid input for fork length!')

    # ascertain each value of time in the origin time distribution
    # and pause time distributions correspond to a length
    # an integer multiple of window size
    if (originTime * forkSpeed) % model.win_size != 0:
        raise ValueError('Invalid input in origin time!')

    if (pauseTime * forkSpeed) % model.win_size != 0:
        raise ValueError('Invalid input in pause time!')

    fork_len_in_wins = int(forkStop / model.win_size)

    fork_name_prefix = "64756d6d-79"

    # set origin location
    origin_loc = round(originTime * forkSpeed / model.win_size)

    # set fork name
    rand_uuid = str(uuid.uuid4())
    fork_name = (fork_name_prefix + rand_uuid[11:] +
                 f"_{contig}_{1}_{forkStop}_fwd_R_{1}_{forkStop}")

    # get brdu data
    if pauseTime > 0:
        pause_time_in_wins = round(pauseTime * forkSpeed / model.win_size)
        win_start, win_end, mean_brdu = model.generate_data(fork_len_in_wins,
                                                            True, model_x_for_t0 + origin_loc, pauseSite,
                                                            pause_time_in_wins)
    else:
        win_start, win_end, mean_brdu = model.generate_data(fork_len_in_wins,
                                                            True, model_x_for_t0 + origin_loc)

    # make a pandas dataframe
    data_for_pandas = {'label': ['winDetect' for _ in mean_brdu],
                       'detectIndex': [fork_name for _ in range(len(win_start))],
                       'mean_brdU': mean_brdu,
                       'start': win_start,
                       'end': win_end}

    return pd.DataFrame(data=data_for_pandas,
                        index=range(len(mean_brdu)))


def get_pause_raw_simpSgm(x, y, w, L, analytical=True, cnv_thres=0.1):
    """ Legacy function. Use get_pause_raw_sgm instead. """
    pause, gof, pos, _, _ = get_pause_raw_sgm(x, y, w, L, analytical, cnv_thres)
    return pause, gof, pos


def get_pause_raw_sgm(x, y, w, L, analytical=True, cnv_thres=0.1,
                      p_low=0, p_high=1, max_iter=None):
    """Calc pause duration, model fit using sgm from threshold raw detect

    Cut raw brdU data into two, fit each part separately to the 
    sigmoid p_low + (p_high - p_low)*sigmoid((x-x0)/w) where
    sigmoid(z) = 1/(1+exp(z)) with two offsets (x0) as a fitting parameter,
    calculate log loss and return duration, goodness of fit (gof) vs
    each candidate pause site.

    NOTE: duration of fit can be negative. Must be discarded by the calling
        function.

    NOTE: an (approximate) analytical solution exists only for the simple
        sigmoid i.e. for p_low = 0, p_high = 1. Do not use analytical = True
        with any other p_low, p_high.

    NOTE: in gof vs x, let's say a = gof obtained by aligning all the
        fork data to the reference sigmoid without any cuts. 
        if function returns gof = a at all sites, then there
        were errors or the instances of non-convergence exceeded
        the set threshold.

    Args:
        x (list of floats): x coordinate, either ref genome pos or equivalent
        y (list of int): thresholded BrdU probability, can be 0 or 1
        w (float): width of sigmoid
        L (float): data length, x ranges from 0 to L
        analytical (bool): (default T) use analytical formulae to fit vs numerical
        cnv_thres (float): if numerical, max fraction of total runs that can fail
        p_low (float): sigmoid probabilities go from p_low to p_high
        p_high (float): sigmoid probabilities go from p_low to p_high
        max_iter (int): default None (inf), max number iterations for each fit in numerical approach

    Returns:
        tuple of pause durations (list, same units as w), gofs (list), pause position (list), err (dict),
            left_cut_pt (list). err has keys 'NoPauseFitNoConv', 'TfInvalidArgError', 'SomeNonConvergence' and
            integer values. If 'NoPauseFitNoConv' > 0, then the routine has failed. The other two integers
            are counts of how many times the respective error occurred.
    """

    # convert to numpy arrays for convenience
    d_type = np.float32
    y = np.array(y, d_type)
    x = np.array(x, d_type)

    # two parts are labelled a and s respectively
    N = len(y)
    s = sum(y)
    a = 0
    cnt = 0

    # ensure x starts from 0
    x0 = x[0]
    x = x - x0

    # do some preliminary checks
    if analytical and not (p_low == 0 and p_high == 1):
        raise ValueError("If using analytical, ensure p_low = 0, p_high = 1!")

    # standard p_low, p_high checks
    if not (0 <= p_low < p_high <= 1):
        raise ValueError('Invalid probability values!')

    # set up an error message structure
    error_msg = ErrorInfo(['NoPauseFitNoConv', 'TfInvalidArgError',
                           'SomeNonConvergence'])

    # initialize gof and pause arrays. to do that
    # we have to first calculate goodness of fit without pauses
    if analytical:

        t1 = np.exp(L * np.mean(y) / w)
        t2 = np.exp(-L / w)

        if t1 == 1 or t1 * t2 == 1:
            raise NotImplementedError('Inputs do not appear to be normal')

        offset = w * np.log((t1 - 1) / (1 - t1 * t2))

        # the standard tensorflow sigmoid is 1/(1+exp(-x)),
        # whereas we have used the definition 1/(1+exp(x)).
        # so we need an extra minus sign in front of x below.

        # although name says Not0to1, will go from 0 to 1
        # if levels are not specified
        gof_no_pause = -np.mean(BernoulliNot0to1().log_prob(y,
                                                            -(x - offset) / w))

    else:
        # we have to use -w instead of w as tf sigmoid and
        # our sigmoid follow opposite sign conventions.
        # i.e. we use 1/(1+exp(x)) whereas tf uses 1/(1+exp(-x))
        try:
            [_, _, is_converged,
             _, log_loss] = [t.numpy() for t in
                             tf_fit_sgm(
                                 np.ones(N, d_type).reshape((N, 1)),
                                 y,
                                 -x / w,
                                 p_low,
                                 p_high,
                                 max_iter)]
        except tf.errors.InvalidArgumentError:
            is_converged = False
            log_loss = [0]

        if not is_converged:
            error_msg.incr('NoPauseFitNoConv')
            log_loss = [0]

        gof_no_pause = np.mean(log_loss)

    gof = np.full((N,), gof_no_pause)
    pause = np.zeros((N,))
    left_cut_pt = np.zeros((N,))

    if error_msg.isSet():
        return [k.tolist() for k in (pause, gof, x + x0, error_msg, left_cut_pt)]

    # cut curve into two parts, see which cut offers
    # best fit and identify that as pause and the difference
    # between offsets as the pause duration

    # analytical approach (works only for simple sigmoid p_low = 0, p_high = 1)
    # ======================================================================
    # logic based on mathematical formulae for
    # best fit offsets of a simple sigmoid, assuming
    # the y data points are evenly spaced along x,
    # i.e. after the cut is performed, we assume
    # y values are evenly spaced over the length of each cut.
    # Assumption shouldn't make a big difference.

    # basically, following logistic regression theory,
    # best-fit offsets to dataset {x_i,y_i} happen when 
    # \sum_i{y_i} = \sum_i{\sigma(z_i)} where
    # z_i = (x_i - x0)/w where x0 is the offset we
    # are looking for and w is known already.
    # assuming uniform spacing in x, we can simplify this to
    # 1/N \sum_i{y_i} = 1/L \int_0^L {\sigma((x-x0)/w) dx}
    # and we can evaluate the integral analytically, and get formula for x0.
    # here \sigma(z) = 1/(1+\exp(z))

    # numerical approach
    # ==================
    # fit both parts of the curve numerically to a simple sigmoid

    for i in range(N - 1):
        a = a + y[i]
        s = s - y[i]
        cnt = cnt + 1

        a_mean = a / cnt
        s_mean = s / (N - cnt)

        if a_mean <= p_low or a_mean >= p_high or s_mean <= p_low or s_mean >= p_high:
            continue

        if analytical:
            l_a = cnt / N * L
            l_s = L - l_a

            ta1 = np.exp(l_a * a_mean / w)
            ta2 = np.exp(-l_a / w)
            ts1 = np.exp(l_s * s_mean / w)
            ts2 = np.exp(-l_s / w)

            off_a = w * np.log((ta1 - 1) / (1 - ta1 * ta2))
            off_s = w * np.log((ts1 - 1) / (1 - ts1 * ts2))

            # the standard tensorflow sigmoid is 1/(1+exp(-x)),
            # whereas we have used the definition 1/(1+exp(x)).
            # so we need an extra minus sign in front of x below.

            # we also need x coordinate to always start from zero
            # for formulae above to be correct

            linear_response = np.concatenate(
                (-(x[0:cnt] - off_a) / w, -(x[cnt:N] - off_s - x[cnt]) / w)
            )

            # although name says Not0to1, will go from 0 to 1
            # if levels are not specified
            log_loss = -1 * BernoulliNot0to1().log_prob(y, linear_response)
            gof[i] = np.mean(log_loss)

            pause[i] = - x[cnt - 1] + off_a - off_s
            left_cut_pt[i] = x[cnt - 1] - off_a

        else:

            # fit each cut separately to a sigmoid
            # NOTE: sigmoid in our rightward fork looks like
            #   p_low + (p_high - p_low) * 1/(1 + exp((x-x0)/w))
            #   where x0 is the offset and w the width
            #   in tensorflow notation, eta = beta_1 * x + beta_0, and sigmoid
            #   is 1/(1 + exp(-eta)). 
            #   Thus, matching coefficients we get beta_0 = x0/w and
            #   beta_1 = -1/w. Also note that only beta_0 is obtained
            #   as a fit parameter, beta_1 is not.

            beta_1 = -1 / w

            try:

                [beta_a, _, is_converged_a,
                 _, log_loss_a] = [t.numpy() for t in
                                   tf_fit_sgm(
                                       np.ones(cnt, d_type).reshape((cnt, 1)),
                                       y[0:cnt],
                                       beta_1 * x[0:cnt],
                                       p_low,
                                       p_high,
                                       max_iter)]

                [beta_s, _, is_converged_s,
                 _, log_loss_s] = [t.numpy() for t in
                                   tf_fit_sgm(
                                       np.ones(N - cnt, d_type).reshape((N - cnt, 1)),
                                       y[cnt:N],
                                       beta_1 * (x[cnt:N] - x[cnt]),
                                       p_low,
                                       p_high,
                                       max_iter)]

            except tf.errors.InvalidArgumentError:
                # unclear why this error is occurring.
                # if it happens, just abort, and act like
                # no pauses were found. most conservative approach.
                error_msg.incr('TfInvalidArgError')
                continue

            # don't bother continuing if routines haven't converged
            if not (is_converged_a and is_converged_s):
                error_msg.incr('SomeNonConvergence')
                continue

            gof[i] = (np.sum(log_loss_a) + np.sum(log_loss_s)) / N

            off_a = w * beta_a[0]
            off_s = w * beta_s[0]

            pause[i] = - x[cnt - 1] + off_a - off_s
            left_cut_pt[i] = x[cnt - 1] - off_a

    # if we get many instances of non-convergence,
    # reset gof and pause arrays.
    if error_msg.count() > cnv_thres * N:
        gof = np.full((N,), gof_no_pause)
        pause = np.zeros((N,))
        left_cut_pt = np.zeros((N,))
        error_msg.incr('NoPauseFitNoConv')

    return [k.tolist() for k in (pause, gof, x + x0, error_msg, left_cut_pt)]


@tf.function(autograph=False, experimental_relax_shapes=True)
def tf_fit_sgm(x, y, offset=None, p_low=0, p_high=1, max_iter=None):
    """Fit sigmoid to given data and return fit params

    Say p = (probability of y = 1).
    Fit p = p_low + (p_high - p_low)*sgm(linear_response + offset)
        where sgm = 1/(1+exp(-x))
        and linear_response = model_coefficients[0] * x

    Example 1:
        >>> x_test = np.linspace(-3, 3, 100, dtype = np.float32).reshape((100,1))
        >>> y_test = (1/(1 + np.exp(-x_test * 0.6))).reshape((100,))
        >>> coeff_test, resp_test, conv_test, n_iter_test, log_loss_test = tf_fit_sgm(x_test, y_test)
        >>> print(list(map(lambda z: round(z, 3), coeff_test.numpy())))
        [0.6]

    Example 2:
        >>> x_test = np.concatenate((np.linspace(-3, 3, 100, dtype = np.float32).reshape((100,1)),
        ...                          np.ones(100, dtype = np.float32).reshape((100,1))), axis = 1)
        >>> y_test = (0.1 + 0.5/(1 + np.exp(-np.dot(x_test,[0.6, -0.18]), dtype = np.float32))).reshape((100,))
        >>> coeff_test, resp_test, conv_test, n_iter_test, log_loss_test = tf_fit_sgm(x_test, y_test,
        ...     p_low=0.1, p_high=0.6)
        >>> print(list(map(lambda z: round(z, 3), coeff_test.numpy())))
        [0.6, -0.18]

    Args:
        x (numpy ndarray shape (n,1)): x coordinate
        y (numpy ndarray shape (n,)): thresholded or raw BrdU probability
        offset (numpy ndarray shape (n,)): (default None) term added to linear response
        p_low (float): default 0, lowest probability level
        p_high (float): default 1, highest probability level
        max_iter (int): default None (means inf), maximum number of iterations

    Returns:
        tuple w 5 items: model_coefficients, linear_response, is_converged, 
            num_iter, log_loss
    """
    model_coefficients, linear_response, is_converged, num_iter = tfp.glm.fit(
        model_matrix=x, response=y, offset=offset,
        maximum_iterations=max_iter,
        model=BernoulliNot0to1(p_low, p_high))
    log_loss = -1 * BernoulliNot0to1(p_low, p_high).log_prob(y,
                                                             linear_response)
    return (model_coefficients, linear_response, is_converged, num_iter,
            log_loss)


def map_known_pause_to_msr_pause_sgm(fork_len, pauseSite,
                                     pauseDuration, width, p_low=0, p_high=1, numeric=False,
                                     max_iter=None):
    """ Set up a right fork w known pause and perform pause measurements

    Goal is to figure out how pause sites and durations are distorted
    by our method. g.o.f means goodness of fit.

    Args:
        fork_len: length of fork in bases
        pauseSite: site of pause in bases
        pauseDuration: duration of pause in bases
        width: width of the simple sigmoid pT->B(x) ref curve in bases
        p_low: (default 0) lowest possible T->B probability
        p_high: (default 1) highest possible T->B probability
        numeric (bool): (default F) T/F = numeric/analytical method
        max_iter (int): default None (=inf), max iterations for numerical method

    Returns:
        fork name, pause site, pause duration, g.o.f min, g.o.f list, seq, err, left cut pt at min,
        pause site list (corresponding to positions of g.o.f list)
    """

    contig = 'chrI'  # fake contig

    fork_name_prefix = "64756d6d-79"  # 'dummy' in ascii

    # set fork name
    rand_uuid = str(uuid.uuid4())
    fork_name = (fork_name_prefix + rand_uuid[11:] +
                 f"_{contig}_{1}_{fork_len}_fwd_R_{1}_{fork_len}")

    # get a dummy sequence and coords of each thymidine
    seq = [random.choice(['A', 'T', 'G', 'C']) for _ in range(fork_len)]
    x = [k for k in range(fork_len) if seq[k] == 'T']  # thymidines

    # some checks on inputs
    if p_low < 0 or p_high > 1 or p_low > p_high or p_low + p_high > 1:
        raise ValueError('Invalid inputs!')

    # get BrdU calls
    # expit calculates 1/(1+exp(-x))
    mid = fork_len / 2

    def step(arg): return 1 if arg > pauseSite else 0

    p = expit([-(k - mid + pauseDuration * step(k)) / width for k in x])
    y = bernoulli.rvs([p_low + (p_high - p_low) * k for k in p])

    # get pause measurements
    pause, gof, pos, err, left_cut_pt = get_pause_raw_sgm(x, y, width, fork_len,
                                                          not numeric, 3 / 100, p_low, p_high, max_iter)

    gof_pos_pause_min = min(zip(gof, pos, pause, left_cut_pt), key=lambda z: z[0] if z[2] >= 0 else np.inf)
    gof_min = gof_pos_pause_min[0]
    pause_at_min = gof_pos_pause_min[2]
    pos_at_min = int(gof_pos_pause_min[1])
    left_cut_pt_min = gof_pos_pause_min[3]

    return fork_name, pos_at_min, pause_at_min, gof_min, gof, seq, err, left_cut_pt_min, pos


def convert_false_fork_and_fit_data_to_fastq(fork_name, seq,
                                             gof_pause_at_thym_list, other_info):
    """ Convert sequence, fit landscpe frm cut-and-align method on frk to fastq

    Our cut-and-align procedure on fake forks with a fake sequence
    generates a pause duration and a goodness of fit per thymidine.
    Storing such data poses a challenge when number of fake forks is large.
    So, we are going to use the fastq format. 

    All non-T bases are mapped to the lowest ascii symbol allowed by fastq
    '!'. At T, the best pause site as fit by the model is marked with
    the next ascii symbol '"', the worst pause site with the last
    ascii symbol "~", and a linear interpolation (not log as is usually done
    by basecallers!) gives symbols for every other site.
    The header contains numerical values of goodness of
    fit for best and worst pause sites and some other information.

    NOTE: the other_info dict must only have alphanumeric characters or . or _
        in keys and values. Items with unsuitable keys/values will be ignored.

    Args:
        fork_name (str): name of fork
        seq (str): DNA sequence
        gof_pause_at_thym_list (list of floats): one non-inf val for each T in seq
        other_info (dict): all keys & vals must be strs w chars alphanumr/./_/-

    Returns:
        str of data in fastq format

    Example: 
    >>> print(convert_false_fork_and_fit_data_to_fastq("name",'ATAGTAT',[0.4,0.5,0.45],{"info1":"4","info2":"5"}),\
        end="")
    @name warning=not_real_genomic_data max=0.50000 min=0.40000 info1=4 info2=5
    ATAGTAT
    +
    !"!!~!P
    """

    ascii_vals = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
    gof_min = min(gof_pause_at_thym_list)
    gof_max = max(gof_pause_at_thym_list)

    if np.isinf(gof_max):
        raise ValueError('Infinite values not allowed!')

    if gof_max != gof_min:
        def map_val_to_char(x):
            return ascii_vals[round(
                (x - gof_min) / (gof_max - gof_min) * (len(ascii_vals) - 2) + 1
            )]
    else:
        def map_val_to_char(x):
            return '"'

    # initialize counter
    cnt = itertools.count()

    # regex criterion
    re_crt = re.compile('^[a-zA-Z\d._-]*$')

    # make fastQ string
    fast_q_str = (f'@{fork_name} warning=not_real_genomic_data '
                  f'max={gof_max:.5f} min={gof_min:.5f} ')

    # statement below converts a dict like {'a':'b','c':'d'}
    # to 'a=b c=d'. Also rejects keys and values
    # of dictionary that have bad characters.
    fast_q_str += ' '.join(filter(None, map(lambda x: x[0] + '=' + x[1] if (bool(re_crt.match(x[0])) and
                                                                            bool(re_crt.match(x[1]))) else '',
                                            zip(other_info.keys(), other_info.values()))))

    fast_q_str += '\n' + seq + '\n' + '+' + '\n'
    fast_q_str += ''.join(
        map(lambda x: map_val_to_char(gof_pause_at_thym_list[next(cnt)]) if x == 'T' else '!', seq)) + '\n'

    return fast_q_str


def get_fit_pause_theor_simpSgm_cntr_fork_width_mismatch(L, w_, w, P):
    """ If ref and fork have diff widths, obtain thr fit pause vs actual pause

    If the reference is a simple sigmoid of some width, and the theoretical
    fork is a simple sigmoid with a central pause but a different width,
    then the fit and the actual pause will presumably not agree, and we
    want to characterize this.

    NOTE: the theoretical fork is not composed of binary data here,
    but we are using the actual underlying sigmoid, as we want to
    perform an ensemble measurement.

    We are basically going to obtain the solution to the equation below
    which gives x0. It can be solved analytically.

    \integral_{P/2}^{L/2 + P/2} \sigma(x/w_) dx = 
    \integral_{P/2}^{L/2 + P/2} \sigma((x-x0)/w) dx

    where sigma(x) = 1/(1 + exp(-x))

    Args:
        L (float): length of the theoretical fork
        w_ (float): width of the theoretical fork
        w (float): width of the reference fork
        P (float): pause duration (in bases) in the theoretical fork

    Returns:
        (float) The fit pause duration

    """

    t1 = np.exp((L + P) / (2 * w_))
    t2 = np.exp(P / (2 * w_))

    t3 = np.exp((L + P) / (2 * w))
    t4 = np.exp(P / (2 * w))

    lhs = w_ * np.log((1 + t1) / (1 + t2))
    t5 = np.exp(lhs / w)

    t6 = (t5 - 1) / (t3 - t4 * t5)

    x0 = -w * np.log(t6)

    return 2 * (P / 2 - x0)


def get_analytical_pause_site_duration_err(ref, refprime, dx, L, xp, P, four_end_points=None):
    """ Analytically get uncertainty in pause duration, site of known pause

    We are given a reference curve yref of a left-moving fork.
    A pause of duration P (in bases) is inserted at the point xp.
    Then, we assume pause site and duration are what we put in,
    and estimate parametric noise associated with the detection i.e.
    uncertainty in pause site and duration.

    NOTE: About outputs. Along the pause site axis, the minimum is
    sharp. So, we assume the underlying parameter distribution is
    exponential. Along the duration axis, the minimum is smooth,
    so we assume the underlying parameter distribution is gaussian.
    (1) So, the pdf for pause sites, f ~ exp(-|xp - xp*|/sigma), so log loss
    chi = -log f = A + |xp-xp*|/sigma. So 
    \lim_{xp \to xp*+} \frac{\partial chi}{\partial xp} = 1/sigma
    \lim_{xp \to xp*-} \frac{\partial chi}{\partial xp} = -1/sigma
    In general, the left and right slopes can be different in magnitude as well,
    so we report both separately.
    (2) Pdf for offsets, f ~ exp(-(xo - xo*)^2/2 * sigma^2), so log loss
    chi = -log f = A + (xo - xo*)^2/2 * sigma^2. So, we report sigma.
    In general, the two offsets can have different errors, so we report
    both. The difference b/w two best-fit offsets xo* is the pause duration.

    Args:
        ref (lambda function): standard fork vs genome coordinate x
        refprime (lambda function): derivative of ref w.r.t. x
        dx (float): grid spacing, use 4b for standard thymidine spacing
        L (float): fork length, fork runs from -L/2 to L/2 by default (see param four_end_points)
        xp (float): pause site
        P (float): pause duration
        four_end_points (list): (default None) alternate way to specify pause location and duration; details follow.
            four entries are left-hand-side and right-hand-side limits of left and right pieces respectively.
            In the default description, fork moves from x=L/2 to x=xp where a pause
            of duration P is inserted, and the rest of the fork moves from x=xp to x=-L/2. An alternate description
            is that there are two pieces of the fork: the right-hand-side before the pause that runs from x=L/2
            to x=xp, and the left-hand-side that runs from x=xp-P to x=-L/2-P. So the four end points are
            [-L/2-P, xp-P, xp, L/2]. NOTE: in this description, ensure that the other parameters of fork length
            and pause duration are consistent.

    Returns:
        (tuple of floats) errors in pause site, offset (see note above)
    """

    if four_end_points is None:
        four_end_points = [-L / 2 - P, xp - P, xp, L / 2]

    if not (abs(four_end_points[3] - four_end_points[0] - L - P) <= dx and
            abs(four_end_points[2] - four_end_points[1] - P) <= dx and
            abs(four_end_points[2] - xp) <= dx):
        raise ValueError("Fork end points disagree with fork length/pause duration")

    # fail if cut portions are too short
    len_left = four_end_points[1] - four_end_points[0]
    len_right = four_end_points[3] - four_end_points[2]
    if len_left <= dx or len_right <= dx or P <= dx:
        return np.nan, np.nan, np.nan, np.nan

    def xLogX(x):
        return x[0] * np.log(x[1])

    def I(x):
        return -xLogX((x[0], x[1])) - xLogX((1 - x[0], 1 - x[1]))

    val_to_left_of_pause = ref(xp - P)
    val_to_right_of_pause = ref(xp)

    s1 = 1 / dx * (
            I((val_to_left_of_pause, val_to_left_of_pause)) -
            I((val_to_left_of_pause, val_to_right_of_pause))
    )
    s2 = 1 / dx * (
            I((val_to_right_of_pause, val_to_left_of_pause)) -
            I((val_to_right_of_pause, val_to_right_of_pause))
    )

    def fn(x):
        return (refprime(x) ** 2) / (ref(x) * (1 - ref(x)))

    def integrate(x):
        return sum(
            x[0](k) * x[3]
            for k in np.arange(start=x[1],
                               stop=x[2], step=x[3])
        )

    j_r = integrate((fn, four_end_points[2], four_end_points[3], dx))
    j_l = integrate((fn, four_end_points[0], four_end_points[1], dx))
    p_e_r = np.sqrt(dx / j_r)
    p_e_l = np.sqrt(dx / j_l)

    return -1 / s1, 1 / s2, p_e_l, p_e_r


def get_loss_function_profile_around_offset(f, x1, x2, dx, off_low, off_high):
    """ Get profile of log loss fn arnd best-fit offset in given interval

    Args:
        f (lambda function): takes one argument, fn whose log loss is needed
        x1 (float): lower limit of interval
        x2 (float): upper limit of interval
        dx (float): grid size
        off_low (float): lower limit of offset
        off_high (float): upper limit of offset

    Returns:
        list of tuples, first entry = offset, second entry = gof val
    """

    # log loss function
    def xLogX(x): return x[0] * np.log(x[1])

    def loss(x): return -xLogX((x[0], x[1])) - xLogX((1 - x[0], 1 - x[1]))

    # standard integration function
    def integrate(x): return sum(
        x[0](m) * x[3]
        for m in np.arange(start=x[1],
                           stop=x[2], step=x[3])
    )

    # integrate and return result
    k = np.arange(start=off_low, stop=off_high, step=dx)

    def gof(x_off): return integrate((
        lambda x: loss((f(x), f(x - x_off))),
        x1, x2, dx)) * 1 / dx

    return list(zip(k, map(gof, k)))


def generate_fitness_landscape_LUT_along_duration_axis(f, log_f, log_one_minus_f, L, dx,
                                                       xm_grid_flt=lambda x: True,
                                                       x0_grid_flt=lambda x: True):
    """ Msr fitness landscape of fitting a fn to itself using log-loss error.

    Preamble
        Consider some mathematical function f(x) from x = a to x = b.
        If f satisfies some conditions, it can be shown that x0 = 0
        minimizes \int_a^b \! g[x0,f(x)] \, \mathrm{d}x where
        g[x0,f(x)] = -f(x) * log(f(x - x0)) - (1 - f(x)) * log(1 - f(x - x0)).
        In other words, minimizing log loss associated with aligning
        a function will end up aligning the function to itself.
    Our problem
        We know best-fit x0 = 0, but what does the fitness landscape
        around x0 = 0 look like?
    Approach
        We consider f(x) from x = -L/2 to x = L/2. Divide x into
        intervals of length dx. In each interval, compute g[x0,f(xm)]
        where xm is the midpoint of the interval and x0 varies from -L/2
        to L/2 in units of dx. So, we end up with a lookup table with
        each row = xm, offset, gof.
    
    Args:
        f (lambda function): takes 1 arg, function whose landscape is needed
        log_f (lambda function): takes 1 arg, log of f
        log_one_minus_f (lambda function): takes 1 arg, log of (1-f)
        L (float): measurements carried out from x = -L/2 to x = L/2
        dx (float): grid size
        xm_grid_flt (lambda fn): default always T, 1 arg -> bool, filter xm grid
        x0_grid_flt (lambda fn): default always T, 1 arg -> bool, filter x0 grid

    Returns:
        pandas df w index xm, columns offset, vals gof. For xm meaning see above

    """

    x0_list = list(itertools.takewhile(lambda x: x <= L / 2,
                                       (-L / 2 + k * dx for k in itertools.count())))

    xm_list = list(itertools.takewhile(lambda x: x <= L / 2,
                                       (-L / 2 + dx / 2 + k * dx for k in itertools.count())))

    return pd.DataFrame(
        [
            (
                xm, x0,
                -f(xm) * log_f(xm - x0) -
                (1 - f(xm)) * log_one_minus_f(xm - x0)
            )
            for xm in filter(xm_grid_flt, xm_list)
            for x0 in filter(x0_grid_flt, x0_list)
        ], columns=['xm', 'x0', 'gof']).pivot(index='x0', columns='xm',
                                              values='gof')


def generate_conf_int_along_duration_axis(lut_ref, dx, cut1_lt_pt_on_ref,
                                          cut1_rt_pt_on_ref, cut2_lt_pt_on_ref, cut2_rt_pt_on_ref, confInt=0.99,
                                          spacBwThyms=4 / 1000, durn1Minus2=True, opPdfs=False):
    """ Obtains conf int in durn given endpoints on ref of two cut segments

    NOTE: (1) we assume coordinate of ref. curve increases from left to right.
          (2) Because of (1), cut1LtPtOnRef < cut1RtPtOnRef and same for 2
          (3) duration = cut1Offset - cut2Offset or cut2Offset - cut1Offset,
              where offset means offset parameter of the reference curve.
              User makes the choice through the parameter durn1Minus2.
              CAUTION: think through this choice carefully.

    Args:
        lut_ref (pandas dataframe): from op of 
            generate_fitness_landscape_LUT_along_duration_axis
        dx (float): grid spacing of offset & xm in lutRef
        cut1_lt_pt_on_ref (float): left-hand side coord on ref of cut segment 1
        cut1_rt_pt_on_ref (float): right-hand side coord on ref of cut segment 1
        cut2_lt_pt_on_ref (float): left-hand side coord on ref of cut segment 2
        cut2_rt_pt_on_ref (float): right-hand side coord on ref of cut segment 2
        confInt (float): default 0.99, b/w 0 and 1, conf. interval desired
        spacBwThyms (float): default 0.004, avg dist bw Ts, set to dx's units
        durn1Minus2 (bool): default T, see note in docstring
        opPdfs (bool): default F, output offset pdfs, durn pdf

    Returns:
        default = 3-float tuple, duration, lower & upper uncertainty in duration;
        if opPdfs is set, then 6-element tuple including above three plus
            pdf of offset1, pdf of offset2, pdf of duration
    """

    if not (0 < confInt < 1 and
            cut1_lt_pt_on_ref <= cut1_rt_pt_on_ref and
            cut2_lt_pt_on_ref <= cut2_rt_pt_on_ref):
        raise ValueError('Suitable params needed')

    # set up filters to extract relevant portions of look up table
    interval1 = filter(lambda x: cut1_lt_pt_on_ref <= x <= cut1_rt_pt_on_ref, lut_ref.columns)
    interval2 = filter(lambda x: cut2_lt_pt_on_ref <= x <= cut2_rt_pt_on_ref, lut_ref.columns)

    # flip one of the vars as we are interested in the pdf of
    # the difference of two variables
    def flipIf1Minus2(x):
        return np.flip(x) if durn1Minus2 else x

    def flipIf2Minus1(x):
        return np.flip(x) if not durn1Minus2 else x

    b1 = flipIf2Minus1(lut_ref[interval1].sum(axis=1))
    b2 = flipIf1Minus2(lut_ref[interval2].sum(axis=1))

    # generate pdf, convolve, get cdf, and then get confidence interval
    pdf1 = np.exp(-dx / spacBwThyms * (b1 - min(b1)))
    pdf2 = np.exp(-dx / spacBwThyms * (b2 - min(b2)))
    pdf1 /= sum(pdf1)
    pdf2 /= sum(pdf2)

    pdf = np.convolve(pdf1, pdf2, mode='full')
    cdf = np.cumsum(pdf, dtype=np.float32)

    def confIntToCdfLevels(x):
        return [(1 - x) / 2, (1 + x) / 2]

    pts_left_right = np.searchsorted(cdf, confIntToCdfLevels(confInt))

    # raise error if pdfs dont have the same length
    assert (len(pdf1) == len(pdf2))

    # find most probable value
    max_pt = np.argmax(pdf)

    op_array = [
        (max_pt - len(pdf1) + 1) * dx,
        (max_pt - pts_left_right[0]) * dx,
        (pts_left_right[1] - max_pt) * dx,
        pdf1.tolist(),
        pdf2.tolist(),
        pdf.tolist()
    ]

    if not opPdfs:
        return op_array[0], op_array[1], op_array[2]
    else:
        return op_array[0], op_array[1], op_array[2], op_array[3], op_array[4], \
               op_array[5]


def fit_quality_to_numeric_vals(qual_str, low_val, high_val, ignore='!'):
    """ Converts ascii string to numeric values

    Args:
        qual_str (string): quality string in format used by basecallers
        low_val (float): lowest quality value, used in mapping
        high_val (float): highest quality value, used in mapping
        ignore (char): default is '!', character to be ignored 

    Returns:
        list of two-element tuples (ref coord, val); ref coord runs from 1-N
    """

    if not ignore == '!':
        raise NotImplementedError('Cannot use any other ignore character')

    ascii_vals = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
    mapped_vals = list(np.linspace(low_val, high_val, len(ascii_vals) - 1))

    return list(
        map(lambda x: (x[0], mapped_vals[ord(x[1]) - ord('\"')]),
            itertools.filterfalse(
                lambda x: x[1] == '!',
                zip(itertools.count(start=1), iter(qual_str))
            )
            )
    )
