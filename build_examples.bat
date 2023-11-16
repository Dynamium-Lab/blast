@echo off
echo 'building examples'

cd examples/build

set includes=-I "../../blast/" -I "../../blast/optional" -I "../nlopt/include/api" -I "../nlopt/include/algs" -I "multiple_compilation_units/"
@REM set options= --std c++17 -m64 -O2 -D BLAST_DEBUG
set options= --std c++17 -Xcompiler "-std:c++17 -O2 /EHsc /MT /Z7 /arch:AVX2 /fp:fast"
@REM set options= --std c++17 -Xcompiler "-std:c++17 -Od -D BLAST_DEBUG /EHsc /MT /Z7 /arch:AVX2 /fp:fast"
set libr= -L "../nlopt/lib"
@REM nvcc ../eval.cu nlopt.lib -o eval.exe %includes% %options% %libr%
@REM nvcc ../main.cu -o main.exe %includes% %options% %libr%
nvcc ../main_sqp.cu nlopt.lib -o main_sqp.exe %includes% %options% %libr%
cl ../multiple_compilation_units/main.cpp ../multiple_compilation_units/second.cpp /EHsc %includes% -std:c++17 -O2 /nologo /link /out:example.exe

cd ../