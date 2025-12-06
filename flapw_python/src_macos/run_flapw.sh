#!/bin/sh

export DYLD_LIBRARY_PATH=$CONDA_PREFIX/lib:$DYLD_LIBRARY_PATH

# Create result directory
# mkdir -p ./flapw_out 

# Copy a new lapwin to the result directory
# cp ./infiles/lapwin ./flapw_out

# Run subprocess to run the flapw file in the result directory
# (cd ./flapw_out && ../bin/FLclean)
# (cd ./flapw_out && ../src_macos/flapw)
# (cd ./flapw_out && ../bin/FLcopy SCF && ../bin/FLclean)

# ./src_macos/flapw
# ./bin/FLcopy SCF
# ./bin/FLclean
# cp lapwinSCF lapwin

../flapw_python/bin/FLclean
../flapw_python/src_macos/flapw
../flapw_python/bin/FLcopy SCF
../flapw_python/bin/FLclean
cp lapwinSCF lapwin