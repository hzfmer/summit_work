#!/use/bin/sh

module swap xl gcc
module load fftw/3.3.5
export LD_LIBRARY_PATH=/gpfs/alpine/proj-shared/geo112/CyberShake/utils/libmemcached_1.0.18/lib:$LD_LIBRARY_PATH
