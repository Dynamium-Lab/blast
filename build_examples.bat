@echo off
echo 'building examples'

cd examples/build

set includes=-I "../../blast/" -I "../../blast/optional" -I "../nlopt/include/api" -I "../nlopt/include/algs" -I "multiple_compilation_units/"
set options= --std c++17 -m64 -O2 -D BLAST_DEBUG
set libr= -L "../nlopt/lib"
nvcc ../eval.cu nlopt.lib -o eval.exe %includes% %options% %libr%
nvcc ../main.cu -o main.exe %includes% %options% %libr%
cl ../multiple_compilation_units/main.cpp ../multiple_compilation_units/second.cpp /EHsc %includes% -std:c++17 -O2 /nologo /link /out:example.exe

cd ../