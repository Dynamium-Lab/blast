

start=`date +%s`
pushd "examples/build"

includes="-I../../blast/ -I../../blast/optional"
options="--std c++17 -m64 -O0 -g -G -D BLAST_DEBUG -Xcompiler -march=native" 
# options="--std c++17 -m64 -O3 -Xcompiler -march=native" 
nvcc ../tests.cu -o tests $includes $options

# options="-std=c++17 -O3 -march=native"
# g++ ../tests.cpp -o tests $includes $options

popd


end=`date +%s`
echo Built in `expr $end - $start` seconds.