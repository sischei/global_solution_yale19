# Set up the environment
module load python/2.7-2015q2

export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH
export IPOPT_DIR="`pwd`/Ipopt-3.12.5/build"
export LD_LIBRARY_PATH=$IPOPT_DIR/lib:$LD_LIBRARY_PATH

