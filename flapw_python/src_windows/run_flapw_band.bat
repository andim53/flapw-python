@echo off

mkdir .\flapw_out
pushd .\flapw_out
cp ..\infiles\bandin .
cp ..\infiles\spkpt.init .
bash ..\bin\FLrst SCF 
..\src\flapw.exe
..\src_windows\xband.exe