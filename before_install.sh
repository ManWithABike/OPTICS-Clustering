mkdir dependencies
cd dependencies

export GTEST_VERSION=master
export GTEST=googletest-${GTEST_VERSION}
wget https://github.com/google/googletest/archive/${GTEST_VERSION}.tar.gz
tar -xzvf *.tar.gz
cd ${GTEST}
mkdir build && cd build
cmake ..
make -j4 && sudo make install && cd ../..

mkdir include
cd include
sudo mkdir /usr/local/include/CImg
wget https://raw.githubusercontent.com/dtschump/CImg/master/CImg.h -O /usr/local/include/CImg/CImg.h
sudo mkdir /usr/local/include/Geometry
wget https://raw.githubusercontent.com/CrikeeIP/Geometry/master/include/geometry/geometry.h -O /usr/local/include/Geometry/geometry.h
cd ..

git clone https://github.com/Dobiasd/FunctionalPlus
cd FunctionalPlus
mkdir build
cd build
cmake ..
sudo make install
cd ..
cd ..

cd ..
