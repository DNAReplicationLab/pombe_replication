import unittest
import numpy as np
import itertools
import pandas as pd
import uuid
import random
import math
from io import StringIO
from scipy.special import expit, log_expit
from scipy.stats import bernoulli, norm, gamma
from numpy.testing import assert_almost_equal
import os

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # stop unnecessary TensorFlow messages

import tensorflow.compat.v2 as tf
import tensorflow_probability as tfp
from fork_arrest_model_tools import calc_offset, calc_sq_diff_w_offset, \
    calc_diff_dot_f_w_offset, BrdUVsPosModel, BrdUVsPosModelSimpleSigmoidal, \
    fit_model_grad_descent, BrdUVsPosModelSigmoidal, \
    lsq_error_model_brute_force, bin_data_w_offset, calc_sq_diff_w_error, \
    process_for_model_fit, calc_offset_bw_models, \
    BrdUVsPosModelSimpleLinear, make_dummy_windowed_fork_data, \
    BrdUVsPosModelSigmoidalClipped, get_pause_raw_simpSgm, tf_fit_sgm, \
    convert_false_fork_and_fit_data_to_fastq, \
    get_fit_pause_theor_simpSgm_cntr_fork_width_mismatch, \
    get_analytical_pause_site_duration_err, \
    get_loss_function_profile_around_offset, \
    fit_quality_to_numeric_vals, \
    generate_fitness_landscape_LUT_along_duration_axis, \
    generate_conf_int_along_duration_axis

tf.enable_v2_behavior()
tfd = tfp.distributions


class TestModelFunctions(unittest.TestCase):
    """ Test whether model classes work correctly """

    def test_BrdUVsPosModel_basic(self):
        """Initialize and test the basic (linear) BrdU model class"""

        model = BrdUVsPosModel(1000)
        model.params = [2]

        self.assertEqual(model.get_model([1, 2, 3, 4]), [2, 4, 6, 8])
        self.assertEqual(model.get_model_grad_wrt_params([1, 2, 3]), [[1, 2, 3]])

        model.set_noise_strength(0.1)
        model_data = model.get_model_with_noise(list(range(1000)))
        residual_data = [m[0] - 2 * m[1] for m in zip(model_data, itertools.count())]

        self.assertAlmostEqual(float(np.mean(residual_data)), 0, places=1)
        self.assertAlmostEqual(float(np.std(residual_data)), 0.1, places=1)

    def test_BrdUVsPosModel_generate_data(self):
        """Test if generating windowed fork data in basic model works """

        model = BrdUVsPosModel(1000)
        model.params = [2]
        model.set_noise_strength(0.1)

        # first, check function with no noise and no pauses
        win_start, win_end, mean_brdu = model.generate_data(4, False, 10)

        self.assertEqual(win_start, [1, 1001, 2001, 3001])
        self.assertEqual(win_end, [1000, 2000, 3000, 4000])
        self.assertEqual(mean_brdu, [20, 22, 24, 26])

        # next, check the function with some noise
        win_start, win_end, mean_brdu = model.generate_data(4, True, 10)

        self.assertEqual(win_start, [1, 1001, 2001, 3001])
        self.assertEqual(win_end, [1000, 2000, 3000, 4000])

        for m in zip(mean_brdu, [20, 22, 24, 26]):
            self.assertNotEqual(m[0], m[1])
            self.assertAlmostEqual(m[0], m[1], delta=0.3)

        # next, check the function with a pause and noise
        model.set_noise_strength(0.01)
        win_start, win_end, mean_brdu = model.generate_data(5, True, 10, 1800, 6)

        self.assertEqual(win_start, [1, 1001, 2001, 3001, 4001])
        self.assertEqual(win_end, [1000, 2000, 3000, 4000, 5000])

        for m in zip(mean_brdu, [20, 24.4, 36, 38, 40]):
            self.assertNotEqual(m[0], m[1])
            self.assertAlmostEqual(m[0], m[1], delta=0.1)

        # perform a rigorous check w noise now
        model.set_noise_strength(0.2)
        win_start, win_end, mean_brdu = model.generate_data(1_000_000, True, 0)
        res = [m[1] - m[0] for m in
               zip(mean_brdu, (2 * k for k in itertools.count()))]

        self.assertAlmostEqual(float(np.mean(res)), 0, delta=0.01)
        self.assertAlmostEqual(float(np.std(res)), 0.2, delta=0.01)

    def test_BrdUVsPosModelSimpleLinear(self):
        """Initialize and test the linear BrdU model class"""

        model = BrdUVsPosModelSimpleLinear(1000, [1, 2])

        ip = [10, 20]

        for n1, n2 in zip(model.get_model(ip), [21, 41]):
            self.assertAlmostEqual(n1, n2)

        op = model.get_model_grad_wrt_params(ip)

        for n1, n2 in zip(op[0], [1, 1]):
            self.assertAlmostEqual(n1, n2)

        for n1, n2 in zip(op[1], ip):
            self.assertAlmostEqual(n1, n2)

        self.assertEqual(len(op), 2)

    def test_BrdUVsPosModelSimpleSigmoidal(self):
        """Initialize and test the sigmoidal BrdU model class"""

        model = BrdUVsPosModelSimpleSigmoidal(1000)
        model.params = [2]

        ip = [np.log(4), np.log(9)]

        # 1/(1+exp(log(4)/2)) = 1/3
        # 1/(1+exp(log(9)/2)) = 1/4

        for n1, n2 in zip(model.get_model(ip), [1 / 3, 1 / 4]):
            self.assertAlmostEqual(n1, n2)

        # we use an identity to verify sigmoidal differentiation
        # if f = 1/(1+e^(x/a)), f' = f*(1-f)*x/a^2, where prime
        # is differentiation w.r.t. a 

        op = model.get_model_grad_wrt_params(ip)
        for n1, n2 in zip(op[0], [1 / 3 * 2 / 3 * ip[0] / 4, 1 / 4 * 3 / 4 * ip[1] / 4]):
            self.assertAlmostEqual(n1, n2)

        self.assertEqual(len(op), 1)

    def test_BrdUVsPosModelSigmoidal(self):
        """Initialize and test the full sigmoidal BrdU model class"""

        model = BrdUVsPosModelSigmoidal(1000, [1, 2, 4])

        ip = [- 4 * 100, - 4 * np.log(2), 0, 4 * np.log(2), 4 * 100]

        expected_op = [1 + 2, 1 + 2 / 3 * 2, 1 + 2 / 2, 1 + 1 / 3 * 2, 1]

        for n1, n2 in zip(model.get_model(ip), expected_op):
            self.assertAlmostEqual(n1, n2)

        # we will now perform a check of derivatives of the function
        # by comparing numerical derivatives with analytical derivatives
        # specified in the model's definition

        op_grad = model.get_model_grad_wrt_params(ip)

        expected_op_grad = []
        dp = 1 * 10 ** (-7)

        for k in range(3):
            model.params[k] = model.params[k] + dp
            expected_op_grad.append([(m[1] - m[0]) / dp for m in zip(expected_op,
                                                                     model.get_model(ip))])
            model.params[k] = model.params[k] - dp

        for k in range(3):
            for m in range(len(ip)):
                self.assertAlmostEqual(expected_op_grad[k][m], op_grad[k][m])

    def test_BrdUVsPosModelSigmoidalClipped(self):
        """ Initialize and test the sigmoidal BrdU-clipped-invert class"""

        # check that weird parameter values are flagged
        self.assertRaises(NotImplementedError,
                          lambda: BrdUVsPosModelSigmoidalClipped(1000, [0.1, -0.1, 0.2]))

        self.assertRaises(ValueError,
                          lambda: BrdUVsPosModelSigmoidalClipped(1000, [0.8, 0.8, 0.2]))

        self.assertRaises(ValueError,
                          lambda: BrdUVsPosModelSigmoidalClipped(1000, [-0.8, 0.7, 0.2]))

        self.assertRaises(ValueError,
                          lambda: BrdUVsPosModelSigmoidalClipped(1000, [1.8, 0.7, 0.2]))

        self.assertRaises(ValueError,
                          lambda: BrdUVsPosModelSigmoidalClipped(1000, [0.8, 0.8, -0.2]))

        # set up the model
        model = BrdUVsPosModelSigmoidalClipped(1000, [0.2, 0.8, 2])

        # check that gradient function does not work
        self.assertRaises(NotImplementedError,
                          lambda: model.get_model_grad_wrt_params([1]))

        # check if model inversion works
        y = [k / 100 for k in [97, 93, 27, 20, 60, 80, 40]]
        expected_ans = [np.nan, np.nan, np.nan, np.nan,
                        0, -2 * np.log(3), 2 * np.log(3)]

        for m in zip(expected_ans,
                     model.invert_model(y, [0.1, 0.9])):
            if np.isnan(m[0]) or np.isnan(m[1]):
                self.assertTrue(np.isnan(m[0]) and np.isnan(m[1]))
            else:
                self.assertAlmostEqual(m[0], m[1])

    def test_calc_sq_diff_w_error(self):
        """Test calculating sq difference where data is mean +- s.d."""

        cndt_brdu_pos = [0.9, 0.8, 0.7, 0.6]
        bin_dt = [[0.9, 0.8], [np.nan, 0.82, 0.71, 0.62]]
        err_dt = [[0.1, 0.2], [np.nan, 0.02, 0.02, 0.02]]

        self.assertEqual(calc_sq_diff_w_error(cndt_brdu_pos, bin_dt[0], err_dt[0]),
                         0)
        self.assertAlmostEqual(
            calc_sq_diff_w_error(cndt_brdu_pos, bin_dt[1], err_dt[1]), 2.25)


class TestOffsetRelatedFunctions(unittest.TestCase):
    """Test functions related to offset in BrdU data w.r.t. candidate curve"""

    def test_calc_offset(self):
        """Align linear/almost linear data w a linear reference."""

        cndt_brdu_pos = [0.9, 0.8, 0.7, 0.6]
        dataset_brdu_pos = [[0.9, 0.8], [0.82, 0.71, 0.62], [0.8, np.nan, 0.6]]
        expected_op = [0, 1, 1]
        op = calc_offset(cndt_brdu_pos, dataset_brdu_pos)
        self.assertEqual(op, expected_op)

    def test_calc_sq_diff_w_offset(self):
        """Calc sq diff of almost linear/linear data w offset."""

        cndt_brdu_pos = [0.9, 0.8, 0.7, 0.6]
        dataset_brdu_pos = [[0.9, 0.8], [0.82, 0.71, 0.62], [0.8, np.nan, 0.6]]
        offset = [0, 1, 1]
        expected_op = [0, 0.0009, 0]
        op = calc_sq_diff_w_offset(cndt_brdu_pos, dataset_brdu_pos, offset)

        # due to floating point errors, have to do almost equal
        for n1, n2 in zip(op, expected_op):
            self.assertAlmostEqual(n1, n2)

    def test_calc_diff_dot_f_w_offset(self):
        """Calc diff dot f of almost linear/linear data w offset."""

        cndt_brdu_pos = [0.9, 0.8, 0.7, 0.6]
        dataset_brdu_pos = [[0.9, 0.8], [0.82, 0.71, 0.62], [0.8, np.nan, 0.6]]
        fn_brdu_pos = [1, 2, 3, 4]
        offset = [0, 1, 1]
        expected_op = [0, 0.15, 0]
        op = calc_diff_dot_f_w_offset(cndt_brdu_pos, fn_brdu_pos, dataset_brdu_pos, offset)

        # due to floating point errors, have to do almost equal
        for n1, n2 in zip(op, expected_op):
            self.assertAlmostEqual(n1, n2)

    def test_bin_data_w_offset(self):
        """Tests if data binning is working correctly"""

        # set some dummy data and bin it.
        dataset_brdu_pos = [[0.9, 0.8, 0.7, 0.6, 0.5], [1.0, 0.9, 0.8, 0.7, 0.6], [0.8, 0.7, 0.6, 0.5, 0.4], [0.2, 0.1],
                            [0.3, 0.2]]

        offset = [0, 0, 0, 1, 3]

        t1 = bin_data_w_offset(dataset_brdu_pos, offset, 4)

        # check to see if means and s.d.s are calc correctly.
        # we round as it is easier to compare w expected nums.
        for k in range(1, len(t1[1])):
            t1[1][k] = round(t1[1][k], 2)
            t1[0][k] = round(t1[0][k], 3)

        self.assertEqual(t1, ([np.nan, 0.65, 0.55, 0.525, 0.425],
                              [np.nan, 0.31, 0.31, 0.17, 0.17]))

    def test_calc_offset_bw_models(self):
        """ Test calculation of offset between models """

        m1 = BrdUVsPosModelSimpleLinear(1000, [-1, 2])
        m2 = BrdUVsPosModelSimpleLinear(1000, [-3, 2])
        x_range = [-30 + 2 * 30 * k / 100 for k in range(101)]
        cnd_offsets = [-3 + k * 0.2 for k in range(31)]
        self.assertEqual(calc_offset_bw_models(m1, m2, x_range, cnd_offsets), 1)


class TestModelFittingLinearGradDescent(unittest.TestCase):
    """ Test that model fitting routines work w a simple lin model """

    def test_fit_model_w_already_fit_lin_model(self):
        """ Test fit model w linear model that already fits the data """

        # set up simple linear model w slope 2
        model = BrdUVsPosModel(1000)
        model.params = [2]

        # give dummy data that is a linear model w slope 2
        dataset_brdu_pos = [[0, 2, 4, 6]]
        offset = [0]
        rel_param_tol = 0.01
        end_win_pos = 20

        n_brdu_wins_data = sum(len(k) for k in dataset_brdu_pos)
        step_size = 0.1 / n_brdu_wins_data

        op_list = fit_model_grad_descent(dataset_brdu_pos, offset, model,
                                         rel_param_tol, [0, end_win_pos], step_size)

        # input model is already the best-fit model, so should
        # come out unchanged.
        self.assertEqual(len(op_list[0].params), 1)
        self.assertAlmostEqual(op_list[0].params[0], 2)
        self.assertEqual(op_list[1], [0])

    def test_fit_model_w_one_array_linear_data(self):
        """ Test fit model w one array of simple linear data """

        # set up simple linear model w slope close to 2 but not 2
        model = BrdUVsPosModel(1000)
        model.params = [1.7]

        # give dummy data that is a linear model w slope 2
        dataset_brdu_pos = [[0, 2, 4, 6]]
        offset = [0]
        rel_param_tol = 0.0001
        end_win_pos = 20

        n_brdu_wins_data = sum(len(k) for k in dataset_brdu_pos)
        step_size = 0.01 / n_brdu_wins_data

        op_list = fit_model_grad_descent(dataset_brdu_pos, offset, model,
                                         rel_param_tol, [0, end_win_pos], step_size)

        # model must have a slope close to 2 now
        self.assertEqual(len(op_list[0].params), 1)
        self.assertAlmostEqual(op_list[0].params[0], 2, delta=0.01)
        self.assertEqual(op_list[1], [0])

    def test_fit_model_w_mult_array_linear_data_w_noise(self):
        """ Test fit model w multiple arrays of linear data + noise """

        # set up simple linear model w slope close to 2 but not 2
        # we need to be very close to 2 in our initial guess as we
        # have a small dataset which produces a landscape with
        # many minima
        model = BrdUVsPosModel(1000)
        model.params = [1.95]

        # make dummy data that is close to a linear model w slope 2
        dataset_brdu_pos = [[9.76, 11.77, 14.02, 16.15, 17.94, 20.1],
                            [5.99, 7.91, 9.94, 11.97, 13.92, 16.15, 17.97],
                            [9.85, 12.13, 14.01, 15.99, 17.85, 19.88, 22.13, 23.91, 25.93,
                             27.92, 30.05, 32.07, 34.08, 36.19],
                            [15.89, 18.15, 20.01, 22.0, 24.08]]
        offset = [0, 0, 0, 0]
        rel_param_tol = 0.0000002
        end_win_pos = 20

        n_brdu_wins_data = sum(len(k) for k in dataset_brdu_pos)
        step_size = 0.0001 / n_brdu_wins_data

        op_list = fit_model_grad_descent(dataset_brdu_pos, offset, model,
                                         rel_param_tol, [0, end_win_pos], step_size)

        # model must have a slope close to 2 now
        self.assertEqual(len(op_list[0].params), 1)
        self.assertAlmostEqual(op_list[0].params[0], 2, delta=0.01)
        self.assertEqual(op_list[1], [5, 3, 5, 8])

    def test_fit_model_w_big_array_linear_data_w_noise(self):
        """ Test fit model w multiple arrays of linear data + noise """

        model = BrdUVsPosModel(1000)
        model.params = [1.4]

        # make dummy data that is close to a linear model w slope 2
        n_data = 200
        dataset_brdu_pos = []
        slp = 2 + 0.1 * np.random.randn(n_data)
        """slopes have a mean of 2 and s.d. of 0.1 """

        for k in range(n_data):
            b = [0, 0]

            while b[1] - b[0] <= 3:
                b = np.sort(np.random.randint(0, 20 + 1, 2))

            dataset_brdu_pos.append([(b[0] + m) * slp[k] + 0.1 * np.random.randn()
                                     for m in range(b[1] - b[0])])

        # obtain best-fit model
        offset = [0 for _ in range(n_data)]
        rel_param_tol = 0.0000002
        end_win_pos = 20

        n_brdu_wins_data = sum(len(k) for k in dataset_brdu_pos)
        step_size = 0.0001 / n_brdu_wins_data

        op_list = fit_model_grad_descent(dataset_brdu_pos, offset, model,
                                         rel_param_tol, [0, end_win_pos], step_size)

        # model must have a slope close to 2 now, within an error.
        # error = 0.1 seems reasonable, as that is the s.d. of slopes
        # used to construct dummy data
        self.assertEqual(len(op_list[0].params), 1)
        self.assertAlmostEqual(op_list[0].params[0], 2, delta=0.1)


class TestModelFittingLinearBruteForce(unittest.TestCase):
    """ Test model fitting brute force routines work w a simp lin modl """

    def test_fit_model_w_mult_array_linear_data_w_noise(self):
        """ Test fit model w multiple arrays of linear data + noise """

        # set up simple linear model w slope close to 2 but not 2
        model = BrdUVsPosModel(1000)

        dataset_brdu_pos = [[9.76, 11.77, 14.02, 16.15, 17.94, 20.1],
                            [5.99, 7.91, 9.94, 11.97, 13.92, 16.15, 17.97],
                            [9.85, 12.13, 14.01, 15.99, 17.85, 19.88, 22.13, 23.91, 25.93,
                             27.92, 30.05, 32.07, 34.08, 36.19],
                            [15.89, 18.15, 20.01, 22.0, 24.08]]
        end_win_pos = 20

        op_list = lsq_error_model_brute_force(dataset_brdu_pos, model,
                                              [[k] for k in [0, 0.5, 1, 1.5, 2, 2.5, 3]], [0, end_win_pos])

        # model must have a slope equal to 2
        self.assertEqual(len(op_list[0].params), 1)
        self.assertEqual(op_list[0].params[0], 2)
        self.assertEqual(op_list[1], [5, 3, 5, 8])

    def test_gof_w_mult_array_linear_data_no_noise(self):
        """ Test goodness of fit calc w mult arrays linear data, no noise """

        # set up simple linear model w slope equal to 2 and no offset
        model = BrdUVsPosModel(1000)

        dataset_brdu_pos = [[0, 2, 4, 6, 8]]
        end_win_pos = 20

        # the param range has to be very small so that the offsets
        # are all the same so that we can test the least squares
        # goodness of fit values around slope = 2
        op_list = lsq_error_model_brute_force(dataset_brdu_pos, model,
                                              [[k] for k in [1.8, 1.9, 2.0, 2.1, 2.2]], [0, end_win_pos])

        # model must have a slope equal to 2 and a perfect fit
        self.assertEqual(len(op_list[0].params), 1)
        self.assertEqual(op_list[0].params[0], 2)
        self.assertEqual(op_list[1], [0])
        self.assertEqual(op_list[2], 0)

        # make sure other goodness of fit nums are correct.
        # goodness of fit at slope m = sum (2*i - m*i)^2
        # with i running from 0 to 4.
        for k in zip(op_list[3],
                     [30 * (2 - m) ** 2 for m in [1.8, 1.9, 2.0, 2.1, 2.2]]):
            self.assertAlmostEqual(k[0], k[1])

        # now, we measure goodness of fit for candidate slope = 1
        # which has a different offset for the data than
        # the best fit slope = 2.
        # NOTE: it's a small calculation to show that the goodness
        # of fit equals 10 and offset equals 2 for slope = 1
        op_list = lsq_error_model_brute_force(dataset_brdu_pos, model,
                                              [[1]], [0, end_win_pos])

        self.assertEqual(len(op_list[0].params), 1)
        self.assertEqual(op_list[0].params[0], 1)
        self.assertEqual(op_list[1], [2])
        self.assertEqual(op_list[2], 10)


class TestModelFittingSigmoidal(unittest.TestCase):
    """ Test that model fitting routines work w a sigmoidal model """

    def __init__(self, *args, **kwargs):
        """ Set up a dummy fork dataset """

        super(TestModelFittingSigmoidal, self).__init__(
            *args, **kwargs)

        self.dataset_brdu_pos = []
        """List of lists with dummy fork data """
        self.xMin = -20
        self.xMax = 20
        """ Maximum, minimum x coordinate possible """

        params = [0.05, 0.75, 2]
        param_noise = [0.01, 0.1, 0.1]  # strengths of parametric noise
        model = BrdUVsPosModelSigmoidal(1000, params)

        # make dummy data that is close to sigmoidal model
        # w parametric noise and measurement noise
        n_data = 200
        min_frk_size = 10

        for k in range(n_data):

            b = [0, 0]
            while b[1] - b[0] <= min_frk_size:
                b = np.sort(np.random.randint(self.xMin, self.xMax + 1, 2))

            model.params = [m1 + m2 * np.random.randn() for m1, m2 in
                            zip(params, param_noise)]

            self.dataset_brdu_pos.append([m + 0.03 * np.random.randn() for m in
                                          model.get_model(list(range(b[0], b[1])))])

    def test_fit_model_w_big_array_data_w_noise_grad_descent(self):
        """ Test fit model w multiple arrays of data + noise w grad descent """

        # obtain best-fit model. of the three params used for initialization,
        # only two are different from the values used to generate the dummy
        # data. If all three are changed, convergence time might be larger
        # than that would be convenient for a test.

        model = BrdUVsPosModelSigmoidal(1000, [0.15, 0.55, 2])
        offset = [0 for _ in range(len(self.dataset_brdu_pos))]
        rel_param_tol = 0.0002

        n_brdu_wins_data = sum(len(k) for k in self.dataset_brdu_pos)
        step_size = 1 / n_brdu_wins_data

        op_list = fit_model_grad_descent(self.dataset_brdu_pos, offset, model,
                                         rel_param_tol, [self.xMin, self.xMax], step_size)

        # ensure requisite num of params are found
        self.assertEqual(len(op_list[0].params), 3)

        # print params and goodness of fit
        print(op_list[0].params)
        print(op_list[2])

        # model must have params close to values used to generate dummy data
        # we give a little room, up to 4x the s.d. param used in dummy data
        self.assertAlmostEqual(op_list[0].params[0], 0.05, delta=0.03)
        self.assertAlmostEqual(op_list[0].params[1], 0.75, delta=0.3)
        self.assertAlmostEqual(op_list[0].params[2], 2, delta=0.4)

    def test_fit_model_w_big_array_data_w_noise_brute_force(self):
        """ Test fit model w multiple arrays of data + noise w brute force """

        # obtain best-fit model. 
        model = BrdUVsPosModelSigmoidal(1000, [0, 0, 0])
        cand_params = list(itertools.product([k / 100 for k in range(10)],
                                             [k / 100 for k in range(50, 90, 2)],
                                             [k for k in np.arange(1, 3, 0.05)]))

        op_list = lsq_error_model_brute_force(self.dataset_brdu_pos, model,
                                              cand_params, [self.xMin, self.xMax])

        # ensure requisite num of params are found
        self.assertEqual(len(op_list[0].params), 3)

        # print params and goodness of fit
        print(op_list[0].params)
        print(op_list[2])

        # model must have params close to values used to generate dummy data
        # we give a little room, up to 4x the s.d. param used in dummy data
        self.assertAlmostEqual(op_list[0].params[0], 0.05, delta=0.03)
        self.assertAlmostEqual(op_list[0].params[1], 0.75, delta=0.3)
        self.assertAlmostEqual(op_list[0].params[2], 2, delta=0.4)


class TestModelFittingSimpleSigmoidal(unittest.TestCase):
    """ Test that model fitting routines work w a simple sigmoidal model """

    def __init__(self, *args, **kwargs):
        """ Set up a dummy fork dataset """

        super(TestModelFittingSimpleSigmoidal, self).__init__(
            *args, **kwargs)

        self.xMin = -20  # smallest window coordinate
        self.xMax = 20  # largest window coordinate
        self.dataset_brdu_pos = []  # list of lists w dummy fork data

        model = BrdUVsPosModelSimpleSigmoidal(1000)
        model.params = [2]

        # make dummy data that is close to sigmoidal model
        # we use some parametric noise to make goodness of
        # fit landscape smoother, as one would expect in real data.

        n_data = 200
        min_frk_size = 10

        for k in range(n_data):

            b = [0, 0]
            while b[1] - b[0] <= min_frk_size:
                b = np.sort(np.random.randint(self.xMin, self.xMax + 1, 2))

            model.params = [2 + 0.2 * np.random.randn()]
            self.dataset_brdu_pos.append([m + 0.1 * np.random.randn() for m in
                                          model.get_model(range(b[0], b[1]))])

    def test_fit_model_w_big_array_data_w_noise_grad_descent(self):
        """ Test fit model w mult arrays of data + noise w grad descent """

        # obtain best-fit model
        model = BrdUVsPosModelSimpleSigmoidal(1000)
        model.params = [4]
        offset = [0 for _ in range(len(self.dataset_brdu_pos))]
        rel_param_tol = 0.0002

        n_brdu_wins_data = sum(len(k) for k in self.dataset_brdu_pos)
        step_size = 100 / n_brdu_wins_data

        op_list = fit_model_grad_descent(self.dataset_brdu_pos, offset, model,
                                         rel_param_tol, [self.xMin, self.xMax], step_size)

        # model must have params close to actual values
        # the estimated parameter could be a little off from 2
        self.assertEqual(len(op_list[0].params), 1)
        self.assertTrue(1.5 < op_list[0].params[0] < 2.5)

    def test_fit_model_w_big_array_data_w_noise_brute_force(self):
        """ Test fit model w mult arrays of data + noise w brute force """

        model = BrdUVsPosModelSimpleSigmoidal(1000)

        op_list = lsq_error_model_brute_force(self.dataset_brdu_pos, model,
                                              [[k] for k in [0.5, 1, 1.5, 2, 2.5, 3]], [self.xMin, self.xMax])

        # model must have params close to actual values
        self.assertEqual(len(op_list[0].params), 1)
        self.assertEqual(op_list[0].params[0], 2)

    def test_fit_model_w_big_array_data_w_noise_brute_force_binned(self):
        """ Test fit model w mult arrays data + noise w brt force w binning """
        model = BrdUVsPosModelSimpleSigmoidal(1000)

        op_list = lsq_error_model_brute_force(self.dataset_brdu_pos, model,
                                              [[k] for k in [0.5, 1, 1.5, 2, 2.5, 3]], [self.xMin, self.xMax],
                                              True, 3)

        # model must have params close to actual values
        self.assertEqual(len(op_list[0].params), 1)
        self.assertEqual(op_list[0].params[0], 2)


class TestProcessForModelFit(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """ Make a file w fake fork data

        Returns:
            None
        """

        dataset_brdu_pos = [
            [
                (9.76, 1, 1000),
                (11.77, 1001, 2000),
                (14.02, 2001, 3000),
                (16.15, 3001, 4000),
                (17.94, 4001, 5000),
                (20.1, 5001, 6000)
            ],
            [
                (5.99, 1, 1000),
                (7.91, 1001, 2000),
                (9.94, 2001, 3000),
                (11.97, 3001, 4000),
                (13.92, 4001, 5000),
                (16.15, 5001, 6000),
                (17.97, 6001, 7000)
            ],
            [
                (9.85, 1, 1000),
                (12.13, 1001, 2000),
                (14.01, 2001, 3000),
                (15.99, 3001, 4000),
                (17.85, 4001, 5000),
                (19.88, 5001, 6000),
                (22.13, 6001, 7000),
                (23.91, 7001, 8000),
                (25.93, 8001, 9000),
                (27.92, 9001, 10000),
                (30.05, 10001, 11000),
                (32.07, 11001, 12000),
                (34.08, 12001, 13000),
                (36.19, 13001, 14000)
            ],
            [
                (15.89, 1, 1000),
                (18.15, 1001, 2000),
                (20.01, 2001, 3000),
                (22.0, 3001, 4000),
                (24.08, 4001, 5000),
                (26.08, 10001, 11999)
            ]
        ]

        dataset_indx = ['fork1', 'fork2', 'fork3', 'fork4']

        # make a fake fork dataframe
        panda_str = "detectIndex mean_brdU start end\n"

        for k in zip(dataset_indx, dataset_brdu_pos):
            for m in k[1]:
                panda_str += f"{k[0]} {m[0]} {m[1]} {m[2]}\n"

        cls.df = pd.read_csv(StringIO(panda_str), sep='\s+',
                             comment="#", index_col=None)

    def test_get_model_fit(self):
        """ Test that model fitting works """

        model = BrdUVsPosModel(1000)

        # IMPORTANT NOTE: best-fit slope is actually 2.001,
        # but we have not included it in the candidate param
        # list. This is because calculation of values
        # output by the routine is much easier if we pretend
        # the best-fit is actually 2.
        df, param = process_for_model_fit(self.df, 1000,
                                          model,
                                          [[k] for k in [0, 0.5, 1, 1.5, 1.996, 1.997,
                                                         1.998, 1.999, 2, 2.002, 2.003,
                                                         2.004, 2.005, 2.5, 3]],
                                          [0, 20],
                                          [500, 1500])

        # prepare our expected-result dataframe
        # to compare to the output
        self.df['windowIndex'] = [5, 6, 7, 8, 9, 10,
                                  3, 4, 5, 6, 7, 8, 9,
                                  5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                                  8, 9, 10, 11, 12, 18]

        self.df['windowRejected'] = ([0 for _ in range(len(self.df.end) - 1)]
                                     + [1])

        self.df['gof'] = [k / 10000 for k in [576, 529, 4, 225, 36, 100,
                                              1, 81, 36, 9, 64, 225, 9,
                                              225, 169, 1, 1, 225, 144, 169, 81, 49, 64, 25, 49, 64, 361,
                                              121, 225, 1, 0, 64, np.nan]]

        self.df['deltaLower'] = np.nan
        self.df['deltaUpper'] = np.nan

        list_df_expected = [
            ('bestFitModel', 0, 1, 1000, 0, 0, 0.3933, 0, 0),
            ('bestFitModel', 2, 1001, 2000, 1, 0, 0.3933, 3 / 1000, 4 / 1000),
            ('bestFitModel', 4, 2001, 3000, 2, 0, 0.3933, 6 / 1000, 8 / 1000),
            ('bestFitModel', 6, 3001, 4000, 3, 0, 0.3933, 9 / 1000, 12 / 1000),
            ('bestFitModel', 8, 4001, 5000, 4, 0, 0.3933, 12 / 1000, 16 / 1000),
            ('bestFitModel', 10, 5001, 6000, 5, 0, 0.3933, 15 / 1000, 20 / 1000),
            ('bestFitModel', 12, 6001, 7000, 6, 0, 0.3933, 18 / 1000, 24 / 1000),
            ('bestFitModel', 14, 7001, 8000, 7, 0, 0.3933, 21 / 1000, 28 / 1000),
            ('bestFitModel', 16, 8001, 9000, 8, 0, 0.3933, 24 / 1000, 32 / 1000),
            ('bestFitModel', 18, 9001, 10000, 9, 0, 0.3933, 27 / 1000, 36 / 1000),
            ('bestFitModel', 20, 10001, 11000, 10, 0, 0.3933, 30 / 1000, 40 / 1000),
            ('bestFitModel', 22, 11001, 12000, 11, 0, 0.3933, 33 / 1000, 44 / 1000),
            ('bestFitModel', 24, 12001, 13000, 12, 0, 0.3933, 36 / 1000, 48 / 1000),
            ('bestFitModel', 26, 13001, 14000, 13, 0, 0.3933, 39 / 1000, 52 / 1000),
            ('bestFitModel', 28, 14001, 15000, 14, 0, 0.3933, 42 / 1000, 56 / 1000),
            ('bestFitModel', 30, 15001, 16000, 15, 0, 0.3933, 45 / 1000, 60 / 1000),
            ('bestFitModel', 32, 16001, 17000, 16, 0, 0.3933, 48 / 1000, 64 / 1000),
            ('bestFitModel', 34, 17001, 18000, 17, 0, 0.3933, 51 / 1000, 68 / 1000),
            ('bestFitModel', 36, 18001, 19000, 18, 0, 0.3933, 54 / 1000, 72 / 1000),
            ('bestFitModel', 38, 19001, 20000, 19, 0, 0.3933, 57 / 1000, 76 / 1000),
            ('bestFitModel', 40, 20001, 21000, 20, 0, 0.3933, 60 / 1000, 80 / 1000),
            ('aggData', np.nan, 1, 1000, 0, 0, np.nan, np.nan, np.nan),
            ('aggData', np.nan, 1001, 2000, 1, 0, np.nan, np.nan, np.nan),
            ('aggData', np.nan, 2001, 3000, 2, 0, np.nan, np.nan, np.nan),
            ('aggData', np.nan, 3001, 4000, 3, 0, np.nan, np.nan, np.nan),
            ('aggData', np.nan, 4001, 5000, 4, 0, np.nan, np.nan, np.nan),
            ('aggData', 9.85, 5001, 6000, 5, 0, 2.777778, 0.09, 0.09),
            ('aggData', 11.95667, 6001, 7000, 6, 0, 0.057719, 0.18037, 0.18037),
            ('aggData', 13.98333, 7001, 8000, 7, 0, 0.091575, 0.055076, 0.055076),
            ('aggData', 16.045, 8001, 9000, 8, 0, 0.123727, 0.127932, 0.127932),
            ('aggData', 17.9775, 9001, 10000, 9, 0, 0.031991, 0.125797, 0.125797),
            ('aggData', 19.99667, 10001, 11000, 10, 0, 0.000908, 0.110604, 0.110604),
            ('aggData', np.nan, 11001, 12000, 11, 0, np.nan, np.nan, np.nan),
            ('aggData', np.nan, 12001, 13000, 12, 0, np.nan, np.nan, np.nan),
            ('aggData', np.nan, 13001, 14000, 13, 0, np.nan, np.nan, np.nan),
            ('aggData', np.nan, 14001, 15000, 14, 0, np.nan, np.nan, np.nan),
            ('aggData', np.nan, 15001, 16000, 15, 0, np.nan, np.nan, np.nan),
            ('aggData', np.nan, 16001, 17000, 16, 0, np.nan, np.nan, np.nan),
            ('aggData', np.nan, 17001, 18000, 17, 0, np.nan, np.nan, np.nan),
            ('aggData', np.nan, 18001, 19000, 18, 0, np.nan, np.nan, np.nan),
        ]

        df_expected = pd.DataFrame(list_df_expected, columns=['detectIndex',
                                                              'mean_brdU', 'start', 'end', 'windowIndex',
                                                              'windowRejected',
                                                              'gof', 'deltaLower', 'deltaUpper'])

        self.df = pd.concat([self.df, df_expected], sort=False)

        # check the output
        self.assertEqual(param['bestFitParamList'], [2])
        self.assertEqual(param['ciParamList'],
                         [[1.997], [1.998], [1.999], [2], [2.002], [2.003], [2.004]])
        pd.testing.assert_frame_equal(df, self.df, rtol=0.01)


class TestMakeDummyWindowedForkData(unittest.TestCase):

    def test_make_dummy_data(self):
        """ Test that making dummy fork data works """

        model = BrdUVsPosModel(1000)
        model.params = [2]

        contig = 'chrI'

        origin_time = 2
        model_x_for_t0 = 10

        pause_site = 8_700
        pause_time = 3.5

        fork_stop = 10_000
        fork_speed = 2_000

        df = make_dummy_windowed_fork_data(model, contig,
                                           origin_time, model_x_for_t0, pause_site, pause_time,
                                           fork_stop, fork_speed)

        sample_uuid = ""

        # first run checks w a pause on
        expected_mean_brdu = [28, 30, 32, 34, 36, 38, 40, 42, 48.214, 60]

        self.assertTrue(len(df.columns) == 5)
        self.assertTrue(all(k == 'winDetect' for k in df.label))
        self.assertTrue(list(df.index) ==
                        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        self.assertTrue(list(df.start) ==
                        [1, 1001, 2001, 3001, 4001, 5001, 6001, 7001, 8001, 9001])
        self.assertTrue(list(df.end) ==
                        [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000])
        self.assertTrue(list(df.mean_brdU) ==
                        expected_mean_brdu)

        for k in df.detectIndex:
            uuid.UUID(k[0:36])
            sample_uuid = k[0:36]
            self.assertTrue(k[0:11] == "64756d6d-79")
            self.assertTrue(k[36:] == "_chrI_1_10000_fwd_R_1_10000")

        self.assertTrue(all(k == df.detectIndex[0] for k in df.detectIndex))

        # now, turn noise on and perform checks
        model.set_noise_strength(0.2)
        res = []

        for k in range(10_000):
            df = make_dummy_windowed_fork_data(model, contig,
                                               origin_time, model_x_for_t0, pause_site, pause_time,
                                               fork_stop, fork_speed)
            res.extend([m[1] - m[0] for m in
                        zip(expected_mean_brdu, list(df.mean_brdU))])

        # ensure the same uuid is not output again
        for k in df.detectIndex:
            self.assertNotEqual(sample_uuid, k[0:36])

        self.assertAlmostEqual(float(np.mean(res)), 0, delta=0.01)
        self.assertAlmostEqual(float(np.std(res)), 0.2, delta=0.01)

        # remove the pause and perform checks
        pause_time = 0
        expected_mean_brdu = [28, 30, 32, 34, 36, 38, 40, 42, 44, 46]

        res = []

        for k in range(10_000):
            df = make_dummy_windowed_fork_data(model, contig,
                                               origin_time, model_x_for_t0, pause_site, pause_time,
                                               fork_stop, fork_speed)
            res.extend([m[1] - m[0] for m in
                        zip(expected_mean_brdu, list(df.mean_brdU))])

        self.assertAlmostEqual(float(np.mean(res)), 0, delta=0.01)
        self.assertAlmostEqual(float(np.std(res)), 0.2, delta=0.01)


class TestGetPauseRawSimpSgm(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        """ Set up a dummy fork dataset """

        super(TestGetPauseRawSimpSgm, self).__init__(
            *args, **kwargs)

        self.N = 30_000  # number of bases in read
        self.width = 5_000  # width of sigmoid
        mid = 15_000  # location of midpoint of sigmoid

        self.pauseSite = 10_000  # site of pause
        self.pauseDuration = 2 * self.width  # duration of pause

        # dna mol
        seq = [random.choice(['A', 'T', 'G', 'C']) for _ in range(self.N)]
        self.x = [k for k in range(self.N) if seq[k] == 'T']  # thymidines

        # expit calculates 1/(1+exp(-x))
        def step(x): return 1 if x > self.pauseSite else 0

        p = expit([-(k - mid + self.pauseDuration * step(k)) / self.width
                   for k in self.x])
        self.y = bernoulli.rvs(p)

    def test_get_pause_raw_simpSgm_rightward(self, analytical=True):
        """ Test that getting pause data works for a rightward fork """

        self.test_get_pause_raw_simpSgm(analytical, direction="right")

    def test_get_pause_raw_simpSgm(self, analytical=True, direction="left"):
        """ Test that getting pause data works for a leftward/rightward fork
        
        We are only going to do a simple test here, as calculations
        are complex and doing a full test will take time.

        Basically, I will generate a fork with a pause in it, and 
        see if the program catches it i.e. see if the detected
        pause location, pause duration are within some tolerance
        """

        if direction == "left":
            pause, gof, pos = get_pause_raw_simpSgm(
                self.x, self.y[::-1], -self.width, self.N, analytical)
        elif direction == "right":
            pause, gof, pos = get_pause_raw_simpSgm(
                self.x, self.y, self.width, self.N, analytical)
        else:
            raise NotImplementedError("direction must be left or right")

        gof_min = min(gof)
        pause_at_min = pause[gof.index(gof_min)]
        pos_at_min = pos[gof.index(gof_min)]

        # we have used quite a bit of tolerance here,
        # that does not mean we expect the variation
        # in fit parameters to be this wide. this is just
        # so that we can test that the function works.
        self.assertEqual(pos, self.x)
        print(f'pos_at_min {pos_at_min}')
        print(f'pause_at_min {pause_at_min}')
        self.assertTrue(self.pauseDuration * 0.8 < pause_at_min < self.pauseDuration * 1.2)

        if direction == "left":
            self.assertTrue((self.N - self.pauseSite) - 1000 < pos_at_min < (self.N - self.pauseSite) + 1000)
        elif direction == "right":
            self.assertTrue(self.pauseSite - 1000 < pos_at_min < self.pauseSite + 1000)
        else:
            raise NotImplementedError("direction must be left or right")

    def test_get_pause_raw_simpSgm_rightward_numerical(self):
        """ Test that getting pause data works if done numerically """
        self.test_get_pause_raw_simpSgm_rightward(analytical=False)

    def test_get_pause_raw_simpSgm_leftward_numerical(self):
        """ Test that getting pause data works if done numerically """
        self.test_get_pause_raw_simpSgm(analytical=False)


class TestTfFitSgmAndBernoulliNot0to1(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        """ Set up a dummy sigmoid dataset 
        
        Dummy dataset is
        (prob y = 1) = 0.1 + 0.5 * sigmoid(z)
            where z = 0.5*x + 0.4.
                sigmoid(z) = 1/(1+exp(-z))
        and x = n points sampled uniformly between -3 and 3.
        """

        super(TestTfFitSgmAndBernoulliNot0to1, self).__init__(
            *args, **kwargs)

        dtype = np.float32  # data type
        n = int(1e4)  # number of data points
        self.beta0 = 0.5  # model slope
        self.beta1 = 0.4  # model offset
        model_coefficients = tf.constant([self.beta0, self.beta1])

        ones = tf.ones([n, 1], dtype)
        x = tfd.Uniform(
            low=-3, high=np.array(3, dtype)).sample([n, 1], seed=45)
        model_matrix = tf.concat([x, ones], 1)

        # multiplies matrices
        linear_response = tf.linalg.matvec(model_matrix, model_coefficients)

        # get response variable
        response = tf.cast(
            tfd.Bernoulli(probs=0.1 + 0.5 * tf.math.sigmoid(linear_response)).
            sample(seed=45),
            dtype)

        self.n = n
        self.ones = ones.numpy()
        self.x = x.numpy()
        self.response = response.numpy()

    def test_tf_fit_sgm(self):
        """ Obtain offset of sigmoid as fitting parameter and check if correct
        
        As shown in the init function we are generating data using a sigmoid
        with p_high = 0.6, p_low = 0.1, and sigmoid argument = 0.5 * x + 0.4.
        In this function, we know the slope 0.5 x, and act like we don't
        know the intercept 0.4 and obtain it as a fitting parameter.
        Then, we check if our answer is close to 0.4.
        """

        [model_coefficients, _, _, _,
         _] = [t.numpy() for t in tf_fit_sgm(
            self.ones,
            self.response,
            (self.beta0 * self.x).reshape(self.n, ),
            0.1,
            0.6)]

        # check outputs are correct
        print(model_coefficients[0])
        self.assertTrue(0.8 * self.beta1 < model_coefficients[0] < 1.2 * self.beta1)
        # just because we use +- 20% tolerance does not mean the real
        # tolerance is that large. it's just what we use for testing.


class TestConvertFalseForkAndFitDataToFastq(unittest.TestCase):

    def test_convert_false_fork_and_fit_data_to_fastq(self):
        """ Test conversion of fork and fit data to fastq """

        gen_str1 = convert_false_fork_and_fit_data_to_fastq("name", 'ATAGTAT',
                                                            [0.4, 0.5, 0.45], {'a': '0.6rs', 'b': 'bb ', 'c=': 'dd'})
        exp_str1 = ("@name warning=not_real_genomic_data max=0.50000 "
                    "min=0.40000 a=0.6rs") + "\n"
        exp_str1 += "ATAGTAT" + "\n"
        exp_str1 += "+" + "\n"
        exp_str1 += "!\"!!~!P" + "\n"
        self.assertTrue(gen_str1 == exp_str1)


class TestGetFitPauseTheorSimpSgmCntrForkWidthMismatch(unittest.TestCase):

    def test_get_fit_pause_theor_simpSgm_cntr_fork_width_mismatch(self):
        """ Test getting a best fit pause using a ref curve  """

        p_actual = 10
        fork_len = 30
        w = 10
        w_ = 10

        p_mes = get_fit_pause_theor_simpSgm_cntr_fork_width_mismatch(
            fork_len, w_, w, p_actual
        )

        print('Measured and actual pause')
        print(p_mes)
        print(p_actual)

        self.assertAlmostEqual(p_mes, p_actual)


class TestGetAnalyticalPauseSiteDurationErr(unittest.TestCase):

    def test_get_analytical_pause_site_duration_err(self):
        """ Test getting pause errors from a ref curve  """

        w = 5

        # expit calculates 1/(1+exp(-x))
        def ref(x): return expit(x / w)

        def ref_prime(x): return -1 / w * expit(x / w) * (1 - expit(x / w))

        dx = 4 / 1000
        fork_len = 30

        xp_p_pairs = [(10, 3), (0, 5), (0, 2), (-5, 4), (-10, 4), (-12, 2)]

        for xp, P in xp_p_pairs:
            xp_e_1, xp_e_2, p_e_1, p_e_2 = get_analytical_pause_site_duration_err(
                ref, ref_prime, dx, fork_len, xp, P)
            p_e = np.sqrt(p_e_1 ** 2 + p_e_2 ** 2)
            xp_e = 0.5 * (xp_e_1 + xp_e_2)

            # now, we use analytical formulae obtained otherwise
            # to check the numerical calculations
            f3 = (np.log(ref(P - xp) / ref(-xp))) / (P / w)
            xp_e_expected = dx * w / (2 * P) * (
                    1 / (ref(xp) - f3) - 1 / (ref(xp - P) - f3)
            )

            p_e_expected = np.sqrt(w * dx) * np.sqrt(
                1 / (ref(fork_len / 2) - ref(xp)) + 1 / (ref(xp - P) - ref(-fork_len / 2 - P))
            )

            print("obtained errors, pause and duration")
            print(xp_e)
            print(p_e)
            print("expected errors, pause and duration")
            print(xp_e_expected)
            print(p_e_expected)

            self.assertAlmostEqual(xp_e, xp_e_expected,
                                   delta=0.1 * max(xp_e, xp_e_expected))
            self.assertAlmostEqual(p_e, p_e_expected,
                                   delta=0.1 * max(p_e, p_e_expected))


class TestGetLossFunctionProfileAroundOffset(unittest.TestCase):

    def test_get_loss_function_profile_around_offset(self):
        """ Test getting offset uncertainties works """

        w = 5

        # expit calculates 1/(1+exp(-x))
        def ref(x):
            return 0.1 + 0.6 * expit(x / w)

        def ref_prime(x):
            return 0.6 * (-1 / w) * expit(x / w) * (1 - expit(x / w))

        # below, we are going to set up a fork with a pause in it,
        # use our analytical formulae to calculate uncertainty associated
        # with pause duration, then measure goodness of fit function
        # shape around the pause, and see if uncertainty calculated
        # using this shape matches our analytical value.

        dx = 4 / 1000
        fork_len = 2 * w

        xp_e_1, xp_e_2, p_e_1, p_e_2 = get_analytical_pause_site_duration_err(
            ref, ref_prime, dx, fork_len, 0, w)

        l1 = get_loss_function_profile_around_offset(ref, 0, w, dx,
                                                     0, w)

        l2 = get_loss_function_profile_around_offset(ref, -2 * w, -w, dx,
                                                     0, w)

        l3 = get_loss_function_profile_around_offset(ref, 0, w, dx,
                                                     -w, 0)

        l4 = get_loss_function_profile_around_offset(ref, -2 * w, -w, dx,
                                                     -w, 0)

        # find a valley where gof increases from the minimum by 1/2.
        # this corresponds to a 68% CI, assuming a gaussian error profile.

        def find_width(l0, level=1 / 2):
            min_gof = next(l0)[1]
            for k in l0:
                if k[1] >= min_gof + level:
                    return k[0]

        width1 = (find_width(iter(l1)) - find_width(reversed(l3))) / 2
        width2 = (find_width(iter(l2)) - find_width(reversed(l4))) / 2

        p_e_num = np.sqrt(width1 ** 2 + width2 ** 2)
        p_e = np.sqrt(p_e_1 ** 2 + p_e_2 ** 2)

        # test if the two calculations match
        print('Analytical pause uncertainty')
        print(p_e)
        print('Numerical pause uncertainty')
        print(p_e_num)

        self.assertAlmostEqual(p_e, p_e_num,
                               delta=0.1 * max(p_e, p_e_num))


class TestFitQualityToNumericVals(unittest.TestCase):

    def test_fit_quality_to_numeric_vals(self):
        """ Test conversion of fit quality string to numbers """

        qual_str = '~!!|!!!!!!!!q!!!!!!c!!Tdo'

        ll = fit_quality_to_numeric_vals(qual_str, 0, 92)

        self.assertEqual(list(m[0] for m in ll),
                         [1, 4, 13, 20, 23, 24, 25])
        self.assertEqual(list(m[1] for m in ll),
                         [92, 90, 79, 65, 50, 66, 77])


class TestGenerateFitnessLandscapeLUTAlongDurationAxis(unittest.TestCase):

    def test_generate_fitness_landscape_LUT_along_duration_axis(self):
        """ Test generation of fitness landscape """

        # for a simple function like the sigmoid, parts of the log loss
        # function can be calculated analytically, and can be used
        # as an independent check to make sure our routine is working.

        fork_len = 10
        dx = 0.1

        b = generate_fitness_landscape_LUT_along_duration_axis(expit,
                                                               log_expit, lambda x: log_expit(-x), fork_len, dx)

        x0 = b.index
        xm = b.columns

        for k in zip(itertools.count(), b.values):
            assert_almost_equal(
                -expit(xm) * (xm - x0[k[0]]) - log_expit(-(xm - x0[k[0]])),
                k[1],
                decimal=6
            )

        # also ensure x0 is symmetric
        assert_almost_equal(list(x0), list(np.flip(-x0)))


class TestGenerateConfIntAlongDurationAxis(unittest.TestCase):

    def test_gaussian(self):
        """ Test conf int calc using fitness landscape and gaussian functions

        We are going to use some mathematical tricks to verify the function
        is working. 

        Preamble
        ==========
        The log-loss functional for a function and its version shifted by
        x0 is chi = integrate {-f(x) * log(f(x - x0)) - 
            (1 - f(x)) * log(1 - f(x-x0))}.
        We generate probability distribution functions from it using
        p(x0) ~ exp(-chi) and convolve the prob distribution with another
        and then use the cdf to ascertain confidence intervals.
        Due to this complicated process, it's probably impossible
        to check the procedure analytically for even a simple f
        such as a sigmoid.

        Overall idea: Cheat and use incorrect log(f) and log(1-f) functions
        ===================================================================
        For sake of avoiding numerical problems, in defining
        generate_fitness_landscape_LUT_along_duration_axis, we
        do not evaluate log(f) and log(1-f), but instead ask the user
        to supply it. let's say g = log(f(x)) and h = log(1 - f(x)).
        Then, log loss functional is -f(x) * g(x - x0) - (1 - f(x)) * h(x - x0).
        Now, cheat and use g(x) = h(x) = -x^2, and any choice of f.
        Ensuing functional is - h(x-x0) = (x-x0)^2 i.e. the log loss functional
        has acquired a squared-difference form. Hence, the associated pdf
        is p ~ exp(-chi) is gaussian! So, we can use analytical results
        associated with a gaussian. For example, the convolution of two
        gaussians is a gaussian with comb_variance = sqrt(sum of variances).

        Approach 1
        ============
        To check aligning a function to itself in the interval [-l,l]
        chi = \int_{-l}^l (x_m - x_0)^2 dx_m
                 = 2/3 * l * (3 x_0^2 + l^2)

        To convert to pdf, do p ~ exp(-(chi - min(chi))/b) where 
        b = mean spacing b/w thymidines. Hence,
        p ~ exp(-2 x_0^2 l / b) i.e. a gaussian with variance b / 4l.
        Convolution with itself makes the variance b / 2l.
        Thus 68.27% CI size = +- sqrt(b/2l)
             95.45% CI size = +- 2 * sqrt(b/2l)
             99.73% CI size = +- 3 * sqrt(b/2l)

        Approach 2
        ============
        To check aligning a function in [0,l_1] to [0,l_2]
        chi = \int_{0}^l_1 (x_m - x_0)^2 dx_m
                 = const + l_1 * (x_0 - l_1/2)^2

        To convert to pdf, do p ~ exp(-(chi - min(chi))/b) where 
        b = mean spacing b/w thymidines. Hence
        p ~ exp(-(x_0 - l_1/2)^2 l_1 / b) i.e. a gaussian with var b / 2l_1
        and mean l_1/2.
        Convolution with pdf of [0, l_2] (remembering that we are subtracting
        2 from 1, not adding them), the new mean is (l_1 - l_2)/2
        and variance is b / 2l_1 + b / 2l_2

        Thus 68.27% CI size = +- sqrt( (b/2) * (1/l_1 + 1/l_2) )
             95.45% CI size = +- 2 * sqrt( (b/2) * (1/l_1 + 1/l_2) )
             99.73% CI size = +- 3 * sqrt( (b/2) * (1/l_1 + 1/l_2) )
        """

        x_lim: float = 10
        dx = 10 / 1000
        l: float = 0.5

        # gaussian function test
        # ======================

        lut_ref = generate_fitness_landscape_LUT_along_duration_axis(expit,
                                                                     lambda x: -x ** 2, lambda x: -x ** 2, x_lim, dx)

        lt68 = generate_conf_int_along_duration_axis(lut_ref, dx,
                                                     -l, l, -l, l, confInt=0.6827, spacBwThyms=dx)
        lt95 = generate_conf_int_along_duration_axis(lut_ref, dx,
                                                     -l, l, -l, l, confInt=0.9545, spacBwThyms=dx)
        lt99 = generate_conf_int_along_duration_axis(lut_ref, dx,
                                                     -l, l, -l, l, confInt=0.9973, spacBwThyms=dx,
                                                     opPdfs=True)

        one_std = np.sqrt(dx / (2 * l))
        obtained_uncertainty = [lt68[1], lt68[2], lt95[1], lt95[2],
                                lt99[1], lt99[2]]
        expected_uncertainty = [one_std, one_std, 2 * one_std, 2 * one_std,
                                3 * one_std, 3 * one_std]
        obtained_durations = [lt68[0], lt95[0], lt99[0]]

        x_vals_for_pdf_1 = np.arange(start=-x_lim / 2, stop=x_lim / 2 + dx,
                                     step=dx, dtype=np.float64)
        first_pdf = norm.pdf(x_vals_for_pdf_1, loc=0, scale=one_std / np.sqrt(2))
        x_vals_for_pdf_final = np.arange(start=-x_lim, stop=x_lim + dx,
                                         step=dx, dtype=np.float64)
        final_pdf = norm.pdf(x_vals_for_pdf_final, loc=0, scale=one_std)
        first_pdf /= sum(first_pdf)
        final_pdf /= sum(final_pdf)

        print('--- First Gaussian test ---')
        print('Obtained & expected uncertainty, obtained & expected durations')
        print(obtained_uncertainty)
        print(expected_uncertainty)
        print(obtained_durations)
        print([0, 0, 0])

        for k in zip(obtained_uncertainty, expected_uncertainty):
            self.assertTrue(math.isclose(
                k[0], k[1], rel_tol=0.1))

        for k in zip(obtained_durations, [0, 0, 0]):
            self.assertTrue(math.isclose(
                k[0], k[1], rel_tol=0.1))

        print('--- Check pdfs ---')
        assert_almost_equal(first_pdf, lt99[3])
        assert_almost_equal(first_pdf, lt99[4])
        assert_almost_equal(final_pdf, lt99[5])
        print('done')

        print('--- Check lengths of some output arrays ---')
        self.assertEqual(len(lt68), 3)
        self.assertEqual(len(lt95), 3)
        self.assertEqual(len(lt99), 6)
        print('done')

        # second gaussian function test
        # ====================================

        lt68 = generate_conf_int_along_duration_axis(lut_ref, dx,
                                                     0, l, 0, 1.5 * l, confInt=0.6827, spacBwThyms=dx)
        lt95 = generate_conf_int_along_duration_axis(lut_ref, dx,
                                                     0, l, 0, 1.5 * l, confInt=0.9545, spacBwThyms=dx)
        lt99 = generate_conf_int_along_duration_axis(lut_ref, dx,
                                                     0, l, 0, 1.5 * l, confInt=0.9973, spacBwThyms=dx,
                                                     durn1Minus2=False)

        one_std = np.sqrt(dx / 2 * (1 / l + 1 / (1.5 * l)))
        obtained_uncertainty = [lt68[1], lt68[2], lt95[1], lt95[2],
                                lt99[1], lt99[2]]
        expected_uncertainty = [one_std, one_std, 2 * one_std, 2 * one_std,
                                3 * one_std, 3 * one_std]
        obtained_durations = [lt68[0], lt95[0], lt99[0]]
        expected_durations = [-0.5 * l / 2, -0.5 * l / 2, 0.5 * l / 2]

        print('--- Second gaussian test ---')
        print('Obtained & expected uncertainty, obtained & expected durations')
        print(obtained_uncertainty)
        print(expected_uncertainty)
        print(obtained_durations)
        print(expected_durations)

        for k in zip(obtained_uncertainty, expected_uncertainty):
            self.assertTrue(math.isclose(
                k[0], k[1], rel_tol=0.1))

        for k in zip(obtained_durations, expected_durations):
            self.assertTrue(math.isclose(
                k[0], k[1], rel_tol=0.1))

    def test_exp(self):
        """ Test conf int calc using fitness landscape and exp functions 

        Approach same as test_gaussian above, but using exp functions.
        h(x) = x

        To check aligning a function to itself in the interval [0,s]
        chi = \int_{0}^s -(x_m - x_0) dx_m
                 = const + x_0 * s

        To convert to pdf, do p ~ exp(-(chi - min(chi))/b) where 
        b = mean spacing b/w thymidines. Hence
        p ~ exp(- x_0 * s / b).
        NOTE: above is a gamma function with parameter alpha = 1.
        
        Now, we subtract two random variables generated
        using the intervals [0,s_1] and [0,s_2].
        The resulting pdf looks like

        p_{1-2}(x_0) ~ exp(-x_0*s_1/b) for x_0 >= 0 and
        p_{1-2}(x_0) ~ exp(x_0*s_2/b) for x_0 <= 0

        For such a function, peak is at x_0 = 0
            and 99% CI size can be calculated
        """

        fork_len = 30
        dx = 10 / 1000
        s = 0.5
        b = 1

        lut_ref = generate_fitness_landscape_LUT_along_duration_axis(expit,
                                                                     lambda x: x, lambda x: x, fork_len, dx,
                                                                     x0_grid_flt=lambda x: round(x, 8) >= 0)

        lt99 = generate_conf_int_along_duration_axis(lut_ref, dx,
                                                     0, s, 0, 2 * s, confInt=0.99, spacBwThyms=b, opPdfs=True)

        # prepare obtained and expected uncertainties, durations
        obtained_durn = lt99[0]
        expected_durn = 0
        obtained_uncertainty = [lt99[1], lt99[2]]
        expected_uncertainty = [b / (2 * s) * np.log(200 / 3), b / s * np.log(400 / 3)]

        # calculate pdfs for x0
        x0_iter = itertools.takewhile(lambda x: x <= fork_len / 2,
                                      (-fork_len / 2 + k * dx for k in itertools.count()))
        x_vals_for_pdf_1 = list(filter(lambda x: round(x, 5) >= 0, x0_iter))

        first_pdf = gamma.pdf(
            x_vals_for_pdf_1,
            a=1, scale=b / s)

        second_pdf = gamma.pdf(
            x_vals_for_pdf_1,
            a=1, scale=b / (2 * s))

        final_pdf = np.concatenate((
            np.flip(second_pdf * first_pdf[0] / second_pdf[0]),
            first_pdf[1:]
        ))

        first_pdf /= sum(first_pdf)
        second_pdf /= sum(second_pdf)
        final_pdf /= sum(final_pdf)

        # print some measurements and perform tests
        print('--- Gamma test ---')
        print('Obtained & expected uncertainty, obtained & expected durations')
        print(obtained_uncertainty)
        print(expected_uncertainty)
        print(obtained_durn)
        print(expected_durn)

        print('--- Check pdf ---')

        def assert_almost_equal_mod(a1, a2):
            """ Check max diff b/w two arrays is approx 1% of max value 
                of one of the arrays """
            assert_almost_equal(a1, a2,
                                decimal=int(np.floor(-np.log10(max(a1) * 0.01))))

        assert_almost_equal_mod(first_pdf, lt99[3])
        assert_almost_equal_mod(np.flip(second_pdf), lt99[4])
        assert_almost_equal_mod(final_pdf, lt99[5])
        print('done')

        for k in zip(obtained_uncertainty, expected_uncertainty):
            self.assertTrue(math.isclose(
                k[0], k[1], rel_tol=0.1))

        self.assertTrue(math.isclose(obtained_durn, expected_durn,
                                     rel_tol=0.1))


if __name__ == '__main__':
    unittest.main()
