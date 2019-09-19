import argparse
import staNMF as st
import numpy as np
import os

parser = argparse.ArgumentParser()
parser.add_argument('--dataset_dir',
                    default='../data/clean_x.csv', type=str)
parser.add_argument('--folderID',
                    default='', type=str)
parser.add_argument('--k1', default = 15, type = int)
parser.add_argument('--k2', default = 30, type = int)
parser.add_argument('--replicates', default = 100, type = int)
parser.add_argument('--seed', type = int)


args = parser.parse_args()

assert(os.path.isfile(args.dataset_dir))

exampleNMF = st.staNMF(filename = args.dataset_dir,
                       folderID = args.folderID,
                       K1 = args.k1,
                       K2 = args.k2,
                       replicates = args.replicates,
                       parallel = True,
                       seed = args.seed)

nthreads = 8 # int(os.environ['SLURM_CPUS_PER_TASK'])

exampleNMF.runNMF(numThreads = nthreads)
exampleNMF.instability()

# Calculate instability for each K
# np.savetxt("./staNMFDicts/instabilities.csv", \
#             exampleNMF.instabilityarray, \
#             delimiter=",")
