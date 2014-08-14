#!/bin/bash -e

# Install dependencies needed for vergini and count. This script is
# not exhaustive, you may still be missing dependencies. This script
# is not very smart and doesn't check for dependencies before
# attempting to install.
#
# This script was tested on 64-bit Ubuntu 14.04.

TMP_INSTALL_DIR='/tmp/count_install'

confirm_install() {
    read -p "Install $1? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]  # Regex match
    then
        return 0
    else
        return 1
    fi
}

mkdir $TMP_INSTALL_DIR
pushd $TMP_INSTALL_DIR

# GSL
if confirm_install GSL
then
    wget http://gnu.askapache.com/gsl/gsl-1.16.tar.gz
    tar -xzf gsl-1.16.tar.gz
    pushd gsl-1.16
    sudo ./configure && sudo make && sudo make install
    popd
fi

# gfortran
sudo apt-get install gfortran # 4.8

# LAPACK (depends on gfortran)
sudo apt-get install liblapack-dev

# ATLAS
sudo apt-get install libatlas3gf-base
sudo apt-get install libatlas-dev
sudo mkdir /usr/lib/atlas
sudo ln -s /usr/lib/atlas-base/libatlas.so.3 /usr/lib/atlas/libatlas.so
sudo ln -s /usr/lib/atlas-base/libcblas.so.3 /usr/lib/atlas/libcblas.so
sudo ln -s /usr/lib/atlas-base/libf77blas.so.3 /usr/lib/atlas/libf77blas.so
sudo ln -s /usr/lib/atlas-base/liblapack.so.3 /usr/lib/atlas/liblapack.so
