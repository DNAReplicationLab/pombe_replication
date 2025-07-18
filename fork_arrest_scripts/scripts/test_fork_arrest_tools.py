import unittest
import pandas as pd
from io import StringIO
from fork_arrest_tools import window_data_per_index, extract_reads_pauses
from DNAscentTools.modBAM_tools_additional import cigar_to_ref_to_query_tbl


class TestLibToolsWindowing(unittest.TestCase):

    def test_without_and_with_windowing(self):
        """
        Test that we can perform averaging with and without windows.
        """
        data = pd.DataFrame({'val': [1, 2, 3, 4, 5, 6, 7],
                             'detectIndex': [1, 1, 1, 1, 1, 2, 2],
                             'posOnRef': [1, 2, 3, 4, 5, 10, 11]}, index=None)

        # average over whole index blocks
        result = window_data_per_index(data)
        expected_result = pd.DataFrame({'mean_val': [3, 6.5], 'detectIndex': [1, 2],
                                        'start': [1, 10], 'end': [5, 11]})

        for k in ['mean_val', 'detectIndex', 'start', 'end']:
            self.assertEqual(list(result[k].values),
                             list(expected_result[k].values))

        # average in windows of size = 2
        result = window_data_per_index(data, win_size=2)
        expected_result = pd.DataFrame({'mean_val': [1.5, 3.5, 6.5], 'detectIndex': [1, 1, 2],
                                        'start': [1, 3, 10], 'end': [2, 4, 11]})

        for k in ['mean_val', 'detectIndex', 'start', 'end']:
            self.assertEqual(list(result[k].values),
                             list(expected_result[k].values))

        # use right aligned windows now
        result = window_data_per_index(data, win_size=2, dirn='left')
        expected_result = pd.DataFrame({'mean_val': [4.5, 2.5, 6.5], 'detectIndex': [1, 1, 2],
                                        'start': [4, 2, 10], 'end': [5, 3, 11]})

        for k in ['mean_val', 'detectIndex', 'start', 'end']:
            self.assertEqual(list(result[k].values),
                             list(expected_result[k].values))

    def test_sliding_windowing(self):
        """
        Test that we can perform sliding windowing.
        """

        data = pd.DataFrame({'val': [1, 2, 3, 4, 5, 6, 7],
                             'detectIndex': [1, 1, 1, 1, 1, 2, 2],
                             'posOnRef': [1, 2, 3, 4, 5, 10, 11]}, index=None)

        # average in windows of size = 2
        result = window_data_per_index(data, win_size=2, sliding=True)
        expected_result = pd.DataFrame({'mean_val': [1.5, 2.5, 3.5, 4.5, 6.5], 'detectIndex': [1, 1, 1, 1, 2],
                                        'start': [1, 2, 3, 4, 10], 'end': [2, 3, 4, 5, 11]})

        for k in ['mean_val', 'detectIndex', 'start', 'end']:
            self.assertEqual(list(result[k].values),
                             list(expected_result[k].values))

        # average in windows of size = 2 but use direction 'left'
        result = window_data_per_index(data, win_size=2, sliding=True, dirn='left')
        expected_result = pd.DataFrame({'mean_val': [4.5, 3.5, 2.5, 1.5, 6.5], 'detectIndex': [1, 1, 1, 1, 2],
                                        'start': [4, 3, 2, 1, 10], 'end': [5, 4, 3, 2, 11]})

        for k in ['mean_val', 'detectIndex', 'start', 'end']:
            self.assertEqual(list(result[k].values),
                             list(expected_result[k].values))

        # average in windows of size = 2 but use direction 'infer'
        data = pd.DataFrame({'val': [1, 2, 3, 4, 5, 6, 7],
                             'detectIndex': ["_L_", "_L_", "_L_", "_L_", "_L_", "_L_2", "_L_2"],
                             'posOnRef': [1, 2, 3, 4, 5, 10, 11]}, index=None)

        result = window_data_per_index(data, win_size=2, sliding=True, dirn='infer')
        expected_result = pd.DataFrame({'mean_val': [4.5, 3.5, 2.5, 1.5, 6.5],
                                        'detectIndex': ["_L_", "_L_", "_L_", "_L_", "_L_2"],
                                        'start': [4, 3, 2, 1, 10], 'end': [5, 4, 3, 2, 11]})

        for k in ['mean_val', 'detectIndex', 'start', 'end']:
            self.assertEqual(list(result[k].values),
                             list(expected_result[k].values))

        # average in windows of size = 2 but use direction 'infer'
        data = pd.DataFrame({'val': [1, 2, 3, 4, 5, 6, 7],
                             'detectIndex': ["_R_", "_R_", "_R_", "_R_", "_R_", "_R_2", "_R_2"],
                             'posOnRef': [1, 2, 3, 4, 5, 10, 11]}, index=None)

        result = window_data_per_index(data, win_size=2, sliding=True, dirn='infer')
        expected_result = pd.DataFrame({'mean_val': [1.5, 2.5, 3.5, 4.5, 6.5],
                                        'detectIndex': ["_R_", "_R_", "_R_", "_R_", "_R_2"],
                                        'start': [1, 2, 3, 4, 10], 'end': [2, 3, 4, 5, 11]})

        for k in ['mean_val', 'detectIndex', 'start', 'end']:
            self.assertEqual(list(result[k].values),
                             list(expected_result[k].values))

    def test_without_windowing_with_threshold(self):
        """
        Test that we can use thresholds
        """
        data = pd.DataFrame({'val': [1, 2, 3], 'detectIndex': [1, 1, 1],
                             'posOnRef': [1, 2, 3]}, index=None)
        result = window_data_per_index(data, thres=2)
        expected_result = pd.DataFrame({'mean_val': [2 / 3], 'detectIndex': [1],
                                        'start': [1], 'end': [3]}, index=[0])

        for k in ['mean_val', 'detectIndex', 'start', 'end']:
            self.assertEqual(result[k].values,
                             expected_result[k].values)

    def test_duplicate_ignorance(self):
        """
        Test that duplicates are ignored when flag is set,
        and not ignored when the flag is not set.
        """
        data = pd.DataFrame({'val': [1, 2, 3, 4, 5, 6, 7],
                             'detectIndex': [1, 1, 1, 1, 1, 2, 2],
                             'posOnRef': [1, 2, 4, 4, 4, 10, 11]}, index=None)

        # test if duplicates are ignored
        result = window_data_per_index(data, win_size=2)
        expected_result = pd.DataFrame({'mean_val': [6.5], 'detectIndex': [2],
                                        'start': [10], 'end': [11]})

        for k in ['mean_val', 'detectIndex', 'start', 'end']:
            self.assertEqual(list(result[k].values),
                             list(expected_result[k].values))

        # test if duplicates are not ignored
        result = window_data_per_index(data, win_size=2, not_rm_dupls=True)
        expected_result = pd.DataFrame({'mean_val': [1.5, 3.5, 6.5], 'detectIndex': [1, 1, 2],
                                        'start': [1, 4, 10], 'end': [2, 4, 11]})

        for k in ['mean_val', 'detectIndex', 'start', 'end']:
            self.assertEqual(list(result[k].values),
                             list(expected_result[k].values))

    def test_infer_averaging_direction(self):
        """
        Test that we can perform averaging detecting direction automatically.
        """
        data = pd.DataFrame({'val': [0.26, 0.76, 0.25, 0.85, 0.55, 0.52, 0.73, 0.91, 0.07, 0.45],
                             'detectIndex': ["_L_", "_L_", "_L_", "_L_", "_L_",
                                             "_R_", "_R_", "_R_", "_R_", "_R_"],
                             'posOnRef': [1, 2, 3, 4, 5, 100, 101, 102, 103, 104]}, index=None)

        # test if direction is detected automatically
        result = window_data_per_index(data, win_size=2, dirn="infer")
        expected_result = pd.DataFrame({'mean_val': [0.7, 0.505, 0.625, 0.49],
                                        'detectIndex': ["_L_", "_L_", "_R_", "_R_"],
                                        'start': [4, 2, 100, 102],
                                        'end': [5, 3, 101, 103]})

        for k in ['mean_val', 'detectIndex', 'start', 'end']:
            self.assertEqual(list(result[k].values),
                             list(expected_result[k].values))


class TestExtractReadsPauses(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """ Make a fake probability dataframe

        Returns:
            None
        """

        fake_data = ("#dummy fork data\n"
                     "posOnRef\tprobBrdU\tdetectIndex\n"
                     "1\t0.9\treadID1\n"
                     "4\t0.9\treadID1\n"
                     "7\t0.9\treadID1\n"
                     "10\t0.9\treadID1\n"
                     "13\t0.5\treadID1\n"
                     "16\t0.5\treadID1\n"
                     "19\t0.5\treadID1\n"
                     "22\t0.5\treadID1\n"
                     "1\t0.9\treadID2\n"
                     "4\t0.9\treadID2\n"
                     "7\t0.9\treadID2\n"
                     "10\t0.0\treadID2\n"
                     "13\t0.9\treadID2\n"
                     "16\t0.9\treadID2\n"
                     "19\t0.9\treadID2\n"
                     "22\t0.9\treadID2\n")

        print(fake_data)

        cls.df = pd.read_csv(StringIO(fake_data), sep='\s+',
                             comment="#", index_col=None)

    def test_program(self):
        """ Test the program works """

        self.assertEqual(extract_reads_pauses(self.df, 4, 0.3), ['readID1'])
        self.assertEqual(extract_reads_pauses(self.df, 4, 0.4), ['readID1'])
        self.assertEqual(extract_reads_pauses(self.df, 4, 0.5), [])
        self.assertEqual(extract_reads_pauses(self.df, 1, 0.9), ['readID2'])
        self.assertEqual(extract_reads_pauses(self.df, 2, 0.9), [])


if __name__ == '__main__':
    unittest.main()
