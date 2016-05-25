#!/bin/bash

source `dirname $0`/build-ifem-module.sh

# Create symlink so build_module can find the test result converter
mkdir deps
ln -sf $WORKSPACE deps/IFEM

# Default to a serial build if no types are given
if test -z "$BTYPES"
then
  BTYPES="serial"
fi

# Convert to arrays for easy looping
BTYPES_ARRAY=($BTYPES)
TOOLCHAINS=($CMAKE_TOOLCHAIN_FILES)

for BTYPE in "${!BTYPES_ARRAY[@]}"
do
  pushd .
  mkdir -p ${BTYPES_ARRAY[$BTYPE]}/build-IFEM
  cd ${BTYPES_ARRAY[$BTYPE]}/build-IFEM
  build_module "-DCMAKE_INSTALL_PREFIX=$WORKSPACE/${BTYPES_ARRAY[$BTYPE]}/install -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAINS[$BTYPE]}" 1 $WORKSPACE
  test $? -eq 0 || exit 1
  popd
  cp $WORKSPACE/${BTYPES_ARRAY[$BTYPE]}/build-IFEM/testoutput.xml $WORKSPACE/${BTYPES_ARRAY[$BTYPE]}
done

# If no downstream builds we are done
if ! grep -q "with downstreams" <<< $ghprbCommentBody
then
  # Add testsuite names
  for BTYPE in "${!BTYPES_ARRAY[@]}"
  do
    sed -e "s/classname=\"TestSuite\"/classname=\"${BTYPES_ARRAY[$BTYPE]}\"/g" ${WORKSPACE}/${BTYPES_ARRAY[$BTYPE]}/build-IFEM/testoutput.xml > ${WORKSPACE}/${BTYPES_ARRAY[$BTYPE]}/testoutput.xml
  done
  exit 0
fi

# remove cmake rule so apps do not get confused
mv $WORKSPACE/cmake/Modules/FindIFEM.cmake $WORKSPACE

# Downstream revisions
declare -a downstreams
downstreams=(IFEM-Stokes
             IFEM-AdvectionDiffusion
             IFEM-NavierStokes
             IFEM-Elasticity
             IFEM-BeamEx
#             IFEM-FiniteDeformation
             IFEM-ThermoElasticity
#             IFEM-PoroElasticity
             IFEM-OpenFrac)

declare -A downstreamRev
downstreamRev[IFEM-AdvectionDiffusion]=master
downstreamRev[IFEM-BeamEx]=master
downstreamRev[IFEM-Elasticity]=master
downstreamRev[IFEM-FiniteDeformation]=master
downstreamRev[IFEM-NavierStokes]=master
downstreamRev[IFEM-OpenFrac]=master
downstreamRev[IFEM-PoroElasticity]=thermal
downstreamRev[IFEM-Stokes]=master
downstreamRev[IFEM-ThermoElasticity]=master

build_downstreams IFEM

# move cmake rule back in place
mv $WORKSPACE/FindIFEM.cmake $WORKSPACE/cmake/Modules

test $? -eq 0 || exit 1
