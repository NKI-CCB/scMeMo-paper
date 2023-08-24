#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --partition=thin
#SBATCH --time=72:00:00

module load 2022 
module load netCDF/4.9.0-gompi-2022a 
module load CMake/3.23.1-GCCcore-11.3.0 
module load Boost/1.79.0-GCC-11.3.0 
export LIBRARY_PATH=$LIBRARY_PATH:/home/bthijssen/libs/lib64 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/bthijssen/libs/lib64 
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/home/bthijssen/libs/include 
export BCM3_ROOT=/home/bthijssen/bcm3 

cd $TMPDIR
rm -rf *
cp -rf /home/bthijssen/finalized-analyses/1-cell-population/5-cell-cycle-model/data .
cp /home/bthijssen/finalized-analyses/1-cell-population/5-cell-cycle-model/prior_eward.xml .
cp /home/bthijssen/finalized-analyses/1-cell-population/5-cell-cycle-model/likelihood_eward.xml .
cp /home/bthijssen/finalized-analyses/1-cell-population/5-cell-cycle-model/model.xml .
cp /home/bthijssen/finalized-analyses/1-cell-population/5-cell-cycle-model/config20.txt .

/home/bthijssen/bcm3/bin/bcminf -c config20.txt --prior=prior_eward.xml --likelihood=likelihood_eward.xml -j 12 -k 16  --output.folder=output_eward_t24n20e2_j12k16 --progress_update_time=30

cp -rf output_* /home/bthijssen/finalized-analyses/1-cell-population/5-cell-cycle-model/
