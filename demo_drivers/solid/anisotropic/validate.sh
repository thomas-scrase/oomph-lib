#! /bin/sh

# Get the OOMPH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1

#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

# Validation for Anisotropic - for now just runs the executable
#-------------------------------

cd Validation
mkdir RESLT

echo "Running Anisotropic validation "
../square > OUTPUT

echo "done"
echo " " >> validation.log
echo "Running Anisotropic validation " >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

echo "Running Anisotropic validation "
../square > OUTPUT

if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  echo "dummy"
fi