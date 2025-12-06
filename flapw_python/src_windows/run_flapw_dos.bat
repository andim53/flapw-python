@echo off

mkdir .\flapw_out
pushd .\flapw_out
cp ..\infiles\dosin .
bash ..\bin\FLrst SCF 
..\src\flapw.exe 
..\src\xdos