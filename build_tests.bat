@echo off
echo 'building tests and benchmarks'

pushd "examples/build"

set includes=-I "../../blast/" -I "../../blast/optional"
set options= --std c++17 -m64 -O0 -g -D BLAST_DEBUG -Xcompiler "-std:c++17 -Od /EHsc /MT /Zi /arch:AVX2"
@REM set options= --std c++17 -m64 -O2 -Xcompiler "-std:c++17 -O2 /EHsc /MT /Z7 /arch:AVX2"
nvcc ../tests.cu -o tests.exe %includes% %options% %libr%

popd
