rm -f run
rm -f reachdists.bmp

export CC=/usr/bin/gcc-6
export CXX=/usr/bin/g++-6
gcc -v && g++ -v && cmake --version
 
g++ -std=c++11 -g test_main.cpp -o run -lX11 -pthread
./run
