
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
wget https://raw.githubusercontent.com/dtschump/CImg/master/CImg.h
wget https://raw.githubusercontent.com/CrikeeIP/Geometry/master/include/geometry/geometry.h
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