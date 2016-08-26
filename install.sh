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

mkdir dependencies
cd dependencies
git clone https://github.com/Dobiasd/FunctionalPlus
cd FunctionalPlus
mkdir build
cd build
cmake ..
sudo make install
cd ..
cd ..
cd .. 
rm dependencies

cd OPTICS-Clustering
cd test

g++ --std=c++11 -I../include main.cpp -lX11 -lpthread
./a.out