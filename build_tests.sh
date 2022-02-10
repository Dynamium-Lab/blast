

echo 'building tests using Google Test'
cd tests
cmake . -B "build" && cmake --build ./build --config Release
cd ../