#!/usr/bin/python3

import os
import numpy as np

aligned_paths_dir = "/mfs/gdouglas/projects/ATRAPP/cdhit_postprocess/aligned_paths/"
cluster_seqs_aligned_dir = "/mfs/gdouglas/projects/ATRAPP/cdhit_postprocess/cluster_seqs_aligned/"

files = np.array(os.listdir(cluster_seqs_aligned_dir))

chunks = np.array_split(files, 100)

os.makedirs(aligned_paths_dir, exist_ok=True)

for i, chunk in enumerate(chunks):
    file_name = aligned_paths_dir + "chunk_" + str(i) +".txt"
    with open(file_name, 'w') as f:
        for item in chunk:
            f.write(cluster_seqs_aligned_dir + item + "\n")
