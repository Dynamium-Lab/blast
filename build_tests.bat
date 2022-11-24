@echo off
echo 'building tests and benchmarks'

pushd "examples/build"

set includes=-I "../../blast/" -I "../../blast/optional"
set options= --std c++17 -m64 -O0 -g -D BLAST_DEBUG
nvcc ../tests.cu -o tests.exe %includes% %options% %libr%

popd
