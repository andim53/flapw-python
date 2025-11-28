@echo off

mkdir .\flapw_out
pushd .\flapw_out
bash ..\bin\FLclean
..\src_windows\flapw.exe
bash ..\bin\FLcopy SCF
bash ..\bin\FLclean