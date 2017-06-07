rm run
rm reachdists.bmp

export CC=/usr/bin/gcc-6
export CXX=/usr/bin/g++-6

g++ -std=c++11 -g test_main.cpp -o run -lX11 -pthread
./run
