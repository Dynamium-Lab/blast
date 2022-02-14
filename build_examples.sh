
echo 'building tests'
cd examples
cmake . -B "build"
cmake --build ./build --config Debug
cd ../
