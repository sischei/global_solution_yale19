# Install entire sparse grid libraries, IPOPT and PYIPOT at once

# make matlab available 
module load matlab

# unpack Matlab sparse grid library
unzip spinterp_v5.1.1.zip
echo " SPINTERP is unpacked "

# unpack and compile Tasmanian Sparse grid library, and test python example
unzip TasmanianSparseGrids.zip
cd TasmanianSparseGrids
make
# cd InterfacePython
# python example.py  ##test
echo " Tasmanian library is installed "

# Install IPOPT and PYIPOPT
cd ../
cd pyipopt_grace
./install_grace.sh
echo " IPOPT and PYIPOPT is installed "


