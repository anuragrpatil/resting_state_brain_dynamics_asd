#!/bin/bash
#SBATCH --partition=long
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -n 20
#SBATCH --gres=gpu:0
#SBATCH --mem-per-cpu=20000M
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END

module add matlab/R2017b
matlab -nodisplay -r "DMF_main_parallel; quit"
