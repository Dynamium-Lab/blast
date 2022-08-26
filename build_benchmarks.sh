

echo 'building benchmarks using Google Benchmark'
cd benchmarks
cmake . -B "build" && cmake --build ./build --config RelWithDebInfo
cd ../