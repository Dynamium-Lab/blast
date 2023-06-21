@echo off
echo 'building tests and benchmarks'

pushd "examples/build"

set includes=-I "../../blast/" -I "../../blast/optional"
set options= /DEBUG -D BLAST_DEBUG /EHsc /MDd /std:c++17 /Z7
@REM set options= -std:c++17 -O2
cl ../tests.cpp -o tests.exe %includes% %options%

popd