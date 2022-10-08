@echo off
echo 'building examples'

cd examples/build

set includes=-I "../../blast/" -I "../../blast/optional" -I "../nlopt/include/api" -I "../nlopt/include/algs"
set options= --std c++17 -m64 -g
set libr= -L "../nlopt/lib"
nvcc ../eval.cu nlopt.lib -o eval.exe %includes% %options% %libr%

cd ../