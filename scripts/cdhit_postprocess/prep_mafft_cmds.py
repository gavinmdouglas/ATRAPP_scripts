import os

input_dir = 'cluster_seqs'
output_dir = '/mfs/gdouglas/projects/ATRAPP/cdhit_postprocess/cluster_seqs_aligned'
cmds_file = 'mafft_cmds.sh'

with open(cmds_file, 'w') as f:
	for fasta_file in os.listdir(input_dir):
		if fasta_file.endswith('.fa'):
			fasta_path = os.path.join(input_dir, fasta_file)
			cmd = f"mafft --retree 2 --maxiterate 2 --thread 1 --preservecase {fasta_path} > {output_dir}/{fasta_file}"
			f.write(cmd + '\n')
