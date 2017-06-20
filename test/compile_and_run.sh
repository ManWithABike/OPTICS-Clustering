rm -f run
rm -f reachdists.bmp

gcc -v && g++ -v && cmake --version
 
./compile.sh
./run
