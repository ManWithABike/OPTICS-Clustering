cd OPTICS-Clustering
cd include
cd optics
mkdir CImg
cd Cimg
wget https://raw.githubusercontent.com/dtschump/CImg/master/CImg.h
cd ..
mkdir Geometry
cd Geometry
wget https://raw.githubusercontent.com/CrikeeIP/Geometry/master/include/geometry/geometry.h
cd ..
cd ..
cd ..
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

g++ --std=c++11 -lX11 -lpthread -I../include main.cpp
./a.out