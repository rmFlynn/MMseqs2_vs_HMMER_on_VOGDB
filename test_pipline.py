"""Do a fast unittest"""
import os
import unittest
from datetime import datetime
from mmseqs_from_setup_to_sweep import download_all_vog_data, \
    process_vogdb_profile, process_vogdb_target, mmseqs_search, WORKING_DIR

TESTING_DIR = ("/home/projects/DRAM/hmmer_mmseqs2_testing_take_2/mmseqs/"
               "test_" + datetime.now().strftime("%Y_%m_%d_%H_%M"))

class TestAll(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestAll, self).__init__(*args, **kwargs)
        self.sensitivity = 1
        self.evalue = 0.1
        self.threads = 64

    def test_download_all_vog_data(self):
        download_all_vog_data()
        self.assertTrue(os.path.exists(WORKING_DIR + "/VOGDB_data" + \
              "/vog.members.tsv.gz"))
        self.assertTrue(os.path.exists(WORKING_DIR + "/VOGDB_data" + \
              "/vog.proteins.all.fa.gz"))
        self.assertTrue(os.path.exists(WORKING_DIR + "/VOGDB_data" + \
              "/vog.raw_algs.tar.gz"))

    def test_process_vogdb_profile(self):
        process_vogdb_profile(self.threads)
        profile_dir = WORKING_DIR + "/mmseqs_profile"
        self.assertTrue(os.path.exists(profile_dir + '/profile'))

    def test_process_vogdb_target(self):
        process_vogdb_target(self.threads)
        profile_dir ="/mmseqs_profile"
        self.assertTrue(
            os.path.exists(WORKING_DIR + '/to_annotate/vog.proteins.all'))

    def test_search(self):
        mmseqs_search(self.sensitivity, self.evalue, self.threads)
        self.assertTrue(WORKING_DIR + \
                        '/sweep_output/evalue_%s_sens_%f' \
                        %(str(self.evalue), self.sensitivity))

    @classmethod
    def setUpClass(cls):
        pass

    @staticmethod
    def disconnect():
        os.rename(WORKING_DIR, TESTING_DIR)

    @classmethod
    def tearDownClass(cls):
        cls.disconnect()


if __name__ == '__main__':
    unittest.main()


