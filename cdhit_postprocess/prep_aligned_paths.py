#!/usr/bin/python3

import os
import numpy as np

# List all files in "/mfs/gdouglas/projects/ocean_mags/clusters/cluster_seqs_aligned/"
files = np.array(os.listdir("/mfs/gdouglas/projects/ocean_mags/clusters/cluster_seqs_aligned/"))

chunks = np.array_split(files, 100)

for i, chunk in enumerate(chunks):
    file_name = "/mfs/gdouglas/projects/ocean_mags/clusters/aligned_paths/" + "chunk_" + str(i) +".txt"
    with open(file_name, 'w') as f:
        for item in chunk:
            f.write("/mfs/gdouglas/projects/ocean_mags/clusters/cluster_seqs_aligned/" + item + "\n")
