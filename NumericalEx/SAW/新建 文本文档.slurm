#!/bin/bash
#!/bin/bash
#SBATCH -J SAW
#SBATCH -p kshctest
#SBATCH -N 1
#SBATCH --ntasks-per-node=32

module purge
module load compiler/intel/2017.5.239 mpi/hpcx/2.4.1/intel-2017.5.239
module load apps/abinit/8.10.3/hpcx-2.4.1-intel2017
srun  --mpi=pmix_v3 abinit < example.in