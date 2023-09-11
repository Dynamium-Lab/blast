
pushd "examples/build"

includes="-I../../blast/ -I../../blast/optional"
options="--std c++17 -m64 -O0 -g -D BLAST_DEBUG -march=native" 
# set options= --std c++17 -m64 -O3 -march=native
nvcc ../tests.cu -o tests $includes $options

popd