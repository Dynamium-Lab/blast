

echo 'building benchmarks'
cd benchmarks
cmake . -B "build" && cmake --build ./build --config RelWithDebInfo
cd ../