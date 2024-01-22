@echo off
echo 'building examples'

cd examples/build


set includes=-I "../../blast/" -I "../../extern/nlopt/include/api" -I "../../extern/nlopt/include/algs" -I "multiple_compilation_units/" -I "../../extern/" -I "../../extern/eigen"
set options= --std c++17 -m64 -O2 --expt-relaxed-constexpr -Xcompiler "-std:c++17 /MT /Z7 /arch:AVX2"
@REM set options= --std c++17 -m64 -O0 -g --expt-relaxed-constexpr -D BLAST_DEBUG -Xcompiler "-std:c++17 /MT /Zi /arch:AVX2"
set libr= -L "../../extern/nlopt/lib"
@REM nvcc ../eval.cu nlopt.lib -o eval.exe %includes% %options% %libr%
@REM nvcc ../eval_collisions.cu nlopt.lib -o eval_collisions.exe %includes% %options% %libr%
nvcc ../main_collisions.cu nlopt.lib -o main_collisions.exe %includes% %options% %libr%
@REM nvcc ../main.cu -o main.exe %includes% %options% %libr%
@REM cl ../multiple_compilation_units/main.cpp ../multiple_compilation_units/second.cpp /EHsc %includes% -std:c++17 -O2 /nologo /link /out:example.exe
@REM cl ../print_basis_functions.cpp /EHsc %includes% -std:c++17 -O2 /nologo /link /out:print_basis_functions.exe


cd ../