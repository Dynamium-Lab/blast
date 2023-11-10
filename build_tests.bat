@echo off
echo 'building tests and benchmarks'

pushd "examples/build"

set includes=-I "../../blast/" -I "../../blast/optional" -I "../../extern" -I "../../extern/eigen"
set options= --std c++17 --expt-relaxed-constexpr -m64 -O0 -g -D BLAST_DEBUG -Xcompiler "/EHsc /MT /Zi /arch:AVX2 /fp:fast"
@REM set options= --std c++17 -m64 -O2 -Xcompiler /arch:AVX512
@REM set options= --std c++17 --expt-relaxed-constexpr -Xcompiler "-std:c++17 -O2 /EHsc /MT /Z7 /arch:AVX2 /fp:fast"
nvcc ../tests.cu -o tests.exe %includes% %options% %libr%

popd
