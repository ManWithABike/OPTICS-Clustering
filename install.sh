cd include
cd optics
mkdir CImg
wget https://raw.githubusercontent.com/dtschump/CImg/master/CImg.h
cd ..
wget https://raw.githubusercontent.com/CrikeeIP/Geometry/master/include/geometry/geometry.h
cd ..
cd .
cd ..

git clone https://github.com/Dobiasd/FunctionalPlus
cd FunctionalPlus
mkdir build
cd build
cmake ..
sudo make install
cd ..
cd ..

cd OPTICS-Clustering
cd test

gcc -I../include main.cpp
./a.out
