export DYLD_LIBRARY_PATH=$CONDA_PREFIX/lib:$DYLD_LIBRARY_PATH
./bin/FLrst SCF
sed -i '' 's/ 12  12  12 / 12  12  12 /' lapwin
./src_macos/flapw
./bin/FLcopy SCF
./bin/FLclean
cp lapwinSCF lapwin
sed -i '' 's/Density of states:F/Density of states:T/' lapwin
cp ./infiles/dosin .
./bin/FLrst SCF
./src_macos/flapw
./src_macos/xdos
mkdir DOS/
mv $(find . -maxdepth 1 -type f -name '*in') DOS/
mv $(find . -maxdepth 1 -type f -name '*out') DOS/
mv $(find . -maxdepth 1 -type f -name '*xy') DOS/
./bin/FLclean