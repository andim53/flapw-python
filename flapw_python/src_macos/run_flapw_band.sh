#!/bin/sh

export DYLD_LIBRARY_PATH=$CONDA_PREFIX/lib:$DYLD_LIBRARY_PATH
# ./bin/FLrst SCF
# sed -i '' 's/ 12  12  12 / 33  33  33 /' lapwin
# ./src_macos/flapw
# ./bin/FLcopy SCF
# ./bin/FLclean
# cp lapwinSCF lapwin
# sed -i '' 's/Density of states:F/Density of states:T/' lapwin
# cp ./infiles/dosin .
# ./bin/FLrst SCF
# ./src_macos/flapw
# ./src_macos/xdos
# mkdir BAND/
# mv $(find . -maxdepth 1 -type f -name '*in') BAND/
# mv $(find . -maxdepth 1 -type f -name '*out') BAND/
# mv $(find . -maxdepth 1 -type f -name '*xy') BAND/
# ./bin/FLclean

../flapw_python/bin/FLclean
../flapw_python/bin/FLrst SCF
cp lapwinSCF lapwin
sed -i '' 's/Band structure:T/Band structure:F/' lapwin
# sed -i '' 's/Density plot:T/Density plot:F/' lapwin
# sed -i '' 's/Density of states:F/Density of states:T/' lapwin
cp ../flapw_python/infiles/bandin .
cp ../flapw_python/infiles/spkpt.init .
../flapw_python/src_macos/flapw 
../flapw_python/src_macos/xband

mkdir -p BAND/
# mv $(find . -maxdepth 1 -type f -name '*in') BAND/
# mv $(find . -maxdepth 1 -type f -name '*out') BAND/
# mv $(find . -maxdepth 1 -type f -name '*xy') BAND/
../flapw_python/bin/FLclean