@echo off

@REM initial time
set /A _tic=%time:~0,2%*3600 +%time:~3,1%*10*60 +%time:~4,1%*60 +%time:~6,1%*10 +%time:~7,1% >nul



set includes= -I "./blast"
set output= -o "./cuda/build/blast_cuda.exe"
set options= --std c++17 -m64 -G -noeh --gpu-architecture=compute_75 -D_DEBUG -Xcompiler "/nologo /Z7"
nvcc cuda/main.cu %options% %includes% %output%



@REM final time computation
set /A _toc=%time:~0,2%*3600 +%time:~3,1%*10*60 +%time:~4,1%*60 +%time:~6,1%*10 +%time:~7,1% >nul
set /A _elapsed=%_toc%-%_tic%
echo build time: %_elapsed% seconds.