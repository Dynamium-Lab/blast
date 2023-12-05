

start=`date +%s`
pushd "examples/build"

includes="-I../../blast/ -I../../blast/optional -I../../extern -I../../extern/eigen"
# options="--std c++17 -m64 -O0 -g -G -D BLAST_DEBUG -Xcompiler -march=native" 
# options="--std c++17 -m64 -O3 -arch=native --expt-relaxed-constexpr -Xcudafe --diag_suppress=unrecognized_gcc_pragma -Xcudafe --diag_suppress=68" 
# nvcc ../tests.cu -o tests $includes $options -Xcompiler "-std=c++17 -Ofast -march=native -DEIGEN_NO_CUDA"

options="-std=c++17 -O3 -march=native"
g++ ../tests.cpp -o tests $includes $options

# Build with Tracy profiler enabled
# options="-std=c++17 -O3 -march=native"
# g++ ../tests.cpp ../../extern/tracy/public/TracyClient.cpp -o tests $includes $options -DTRACY_ENABLE

popd


end=`date +%s`
echo Built in `expr $end - $start` seconds.
