import multiprocessing as mp
from subprocess import Popen, PIPE
import itertools
import numpy as np
import random
# AN_Ne = [10000]
# FW_Ne = [1000]
# M_AN_TO_FW = np.round(10 ** np.arange(-4, -0.99, 0.5), 5)
#M_FW_TO_AN = 10*M_AN_TO_FW

AN_Ne = [10000]
FW_Ne = [1000]
M_AN_TO_FW = [0.001]
M_FW_TO_AN = [0.01]

#FW_founding = [1000,2000,3000]
FW_founding = [250,500,750,1000,2000,3000]

param_matrix = list(itertools.product(AN_Ne, FW_Ne, M_AN_TO_FW, FW_founding))

seeds = random.sample(range(1, 1000000000), 100)

jobs = [(i, seed) for i in range(len(param_matrix)) for seed in seeds] 

nbthreads = 100
batches = [jobs[i:i + nbthreads] for i in range(0, len(jobs), nbthreads)]

def process_slim(param_index, seed, output):
    an_ne, fw_ne, m, FW_found = param_matrix[param_index]

    slim_command = (f"/gpfs/software/Anaconda/envs/slim/bin/slim -d AN_Ne={an_ne} "f"-d FW_Ne={fw_ne} "f"-d MIG_AN_TO_FW={m} "f"-d MIG_FW_TO_AN={10*m} " f"-d New_Lake_size={FW_found} " f"-d SEED={seed} 341_independent_loci_simulations.slim")

    try:
        process = Popen(slim_command, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()

        if process.returncode == 0:
            output.put(f"Finished seed {seed} with param set {param_index}")
        else:
            output.put(f"Error for seed {seed}, param set {param_index}:\n{stderr.decode()}")

    except Exception as e:
        output.put(f"Exception for seed {seed}, param set {param_index}: {str(e)}")

for batch_index, batch in enumerate(batches):
    print(f"\nLaunching batch {batch_index + 1}/{len(batches)}")
    
    output = mp.Queue()
    processes = []

    for param_index, seed in batch:
        p = mp.Process(target=process_slim, args=(param_index, seed, output))
        processes.append(p)

    for p in processes:
        p.start()

    for p in processes:
        p.join()

    results = [output.get() for _ in processes]
    print(f"\n Batch {batch_index + 1} Results:")
    for r in results:
        print(r)
