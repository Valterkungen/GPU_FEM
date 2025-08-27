#!/bin/bash -l

#SBATCH -A uppmax2022-2-19  # project name
#SBATCH -M snowy            # name of system
#SBATCH -p node             # request a full node
#SBATCH -N 1                # request 1 node
#SBATCH -t 1:00:00          # job takes at most 1 hour
#SBATCH -J stream_triad_cpu # name of the job
#SBATCH --res=XXX			# enter reservation label or delete
#SBATCH -D ./               # stay in current working directory

./stream_triad -min 8 -max 1e8
