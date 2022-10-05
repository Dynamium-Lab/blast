
echo 'building examples'
cd examples
cmake . -B "build" && cmake --build ./build --config Debug
cd ../
