
# setup environment
source setup_env.sh


#Install IPOPT
echo "installing IPOPT now"

cd Ipopt-3.12.5

mkdir -p build

cd build

../configure

make -j12

make test

make install

IPOPT_DIR=`pwd`

echo $IPOPT_DIR

export LD_LIBRARY_PATH=IPOPT_DIR/lib:$LD_LIBRARY_PATH

echo $LD_LIBRARY_PATH

cd ../../
echo " IPOPT is installed"

#install PYIPOPT
tar xfv pyipopt.tar

cd pyipopt

python setup.py build

python setup.py install --user

cd examples

python hs071.py

echo " PYIPOPT is tested and installed"





