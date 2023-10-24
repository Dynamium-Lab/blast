@echo off
echo 'building tests and benchmarks'

pushd "examples/build"

set includes=-I "../../blast/" -I "../../blast/optional"
@REM set options= /DEBUG -D BLAST_DEBUG /EHsc /MDd /std:c++17 /Z7
set options= -std:c++17 -Od /EHsc /MT /Z7 /arch:AVX2 /fp:fast
cl ../tests.cpp -o tests.exe %includes% %options%

popd