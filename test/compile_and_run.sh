rm run
rm reachdists.bmp
g++ -std=c++11 -g test_main.cpp -o run -lX11 -pthread
./run
