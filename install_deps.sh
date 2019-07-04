#!/bin/sh
echo "***********************************************************"
echo "Starting installer for VarSpeedDubins"
echo "***********************************************************"

echo "==========================================================="
echo " Updating apt-get..."
echo "==========================================================="
sudo apt-get update

echo "==========================================================="
echo " Installing g++, cmake, make..."
echo "==========================================================="
sudo apt-get install g++ cmake build-essential -y


echo "================================================================="
echo "Installing math libraries (LAPACK,BLAS,ARMADILLO)..."
echo "================================================================="
sudo apt-get install libopenblas-dev liblapack-dev libblas-dev libarpack++2-dev -y
cd external
tar xvzf armadillo-6.600.5.tar.gz
cd armadillo-6.600.5
cmake .
make 
sudo make install	
cd ..
cd ..
echo "================================================================="
echo "Installing Vertex Enumeration Library CDDLIB..."
echo "================================================================="
sudo apt-get install  libgmp3-dev -y
cd external
tar xvzf cddlib-094h.tar.gz
cd cddlib-094h
./configure
make 
sudo make install	
cd ..
cd ..    #
echo "==========================================================="
echo " Installing plotting tools..."
echo "==========================================================="
sudo apt-get install octave -y
echo "==========================================================="
echo " Installing IPOPT, HSL ..."
echo "==========================================================="
sudo apt-get install gfortran gcc libopenblas-dev liblapack-dev libblas-dev -y
cd external
tar xvzf CoinIpopt_SvnRev3584.tar.gz
tar xvzf coinhsl-linux-x86_64-2014.01.10.tar.gz
mv coinhsl-linux-x86_64-2014.01.10 coinhsl
sudo cp -a coinhsl/lib/. /usr/local/lib/
sudo cp -a coinhsl/include/. /usr/local/include/
sudo ldconfig
mkdir ./CoinIpopt/build
cd ./CoinIpopt/build
../configure
make 
make test
sudo make install	
cd ../../..
 


echo "==========================================================="
echo " Install script complete."
echo "==========================================================="



