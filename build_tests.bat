@echo off
echo 'building tests and benchmarks'

pushd "examples/build"

set includes=-I "../../blast/" -I "../../blast/optional"
@REM set options= --std c++17 -m64 -O0 -g -D BLAST_DEBUG
@REM set options= --std c++17 -m64 -O2 -Xcompiler /arch:AVX512
set options= --std c++17 -Xcompiler "-std:c++17 -Od /EHsc /MT /Z7 /arch:AVX2 /fp:fast"
nvcc ../tests.cu -o tests.exe %includes% %options% %libr%

popd
