rm -f run
rm -f reachdists.bmp

gcc -v && g++ -v && cmake --version
 
g++ -std=c++11 -g test_main.cpp -o run -lX11 -pthread
./run
