"""Download, format data and make a sweep of mmseqs prams"""
import os
from subprocess import Popen
from datetime import datetime

BOUTFMT6_COLUMNS = ['qId', 'tId', 'seqIdentity', 'alnLen', 'mismatchCnt',
                    'gapOpenCnt', 'qStart', 'qEnd', 'tStart', 'tEnd', 'eVal',
                    'bitScore']

# put the file into a unique dir
WORKING_DIR = ("/home/projects/DRAM/hmmer_mmseqs2_testing_take_2/mmseqs/"
               "results_"+ datetime.now().strftime("%Y_%m_%d_%H"))

THREADS = 64

###############
# Vogdb setup #
###############

def download_all_vog_data():
    """Download all the vog db data"""
    assert not os.path.exists(WORKING_DIR), "Results already exist"
    os.mkdir(WORKING_DIR)
    os.chdir(WORKING_DIR)
    os.mkdir('VOGDB_data')
    os.chdir('VOGDB_data')
    os.system("wget http://fileshare.csb.univie.ac.at/vog/latest/"
              "vog.members.tsv.gz")
    os.system("wget http://fileshare.csb.univie.ac.at/vog/latest/"
              "vog.proteins.all.fa.gz")
    os.system("wget http://fileshare.csb.univie.ac.at/vog/latest/"
              "vog.raw_algs.tar.gz")
    os.system("mkdir vog")
    os.chdir(WORKING_DIR)

def process_vogdb_profile(threads):
    os.chdir(WORKING_DIR)
    os.mkdir("mmseqs_profile")
    os.chdir("mmseqs_profile")
    Popen(['mmseqs', 'tar2db',
           WORKING_DIR + "/VOGDB_data/vog.raw_algs.tar.gz", # Input
           "vog_msa",                                      # Output
           "--output-dbtype", "11",
           # "--tar-include", "'.*msa$'",
           "--threads",  str(threads)]).wait()
    os.system("mmseqs msa2profile vog_msa profile --match-mode 1")
    os.chdir(WORKING_DIR)

def process_vogdb_target(threads, sensitivity=None, index=False):
    """process to proteins to be annotated"""
    os.chdir(WORKING_DIR)
    os.system("mkdir to_annotate")
    Popen(['mmseqs', 'createdb',
          WORKING_DIR + '/VOGDB_data/vog.proteins.all.fa.gz',
          'to_annotate/vog.proteins.all']).wait()
    # indexes target DB
    if index:
        assert sensitivity is not None, \
            "you must provide sensitivity to use in the index"
        os.system("mkdir to_annotate/temp")
        Popen(['mmseqs', 'createindex', 'to_annotate/vog.proteins.all',
               'to_annotate/temp', '-k', '6', '-s', str(sensitivity),
               '--threads', str(threads)]).wait()

def mmseqs_search(sensitivity, evalue, threads):
    os.chdir(WORKING_DIR)
    os.system("mkdir raw_sweep_output")
    Popen(['mmseqs', 'search',
           'to_annotate/vog.proteins.all',
           'mmseqs_profile/profile',
           'raw_sweep_output/evalue_%s_sens_%f' %(str(evalue), sensitivity),
           'temp', '--threads', str(threads), '-e', str(evalue), '-s',
           str(sensitivity), '-k', '6']).wait()
    os.system("mkdir sweep_output")
    Popen(['mmseqs', 'convertalis',
           'to_annotate/vog.proteins.all',
           'mmseqs_profile/profile',
           'raw_sweep_output/evalue_%s_sens_%f' %(str(evalue), sensitivity),
           'sweep_output/evalue_%s_sens_%f' %(str(evalue), sensitivity),
           ]).wait()

def main():
    download_all_vog_data()
    process_vogdb_profile(THREADS)
    process_vogdb_target(THREADS)
    for evalue in range(-20, 0, 2):
        for sens in range(2, 16, 1):
            mmseqs_search(sens/2, float("1e%i"%evalue), THREADS)


if __name__ == '__main__':
    main()
