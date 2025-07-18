export DYLD_LIBRARY_PATH=$CONDA_PREFIX/lib:$DYLD_LIBRARY_PATH
./src_macos/flapw
./bin/FLcopy SCF
./bin/FLclean
cp lapwinSCF lapwin

