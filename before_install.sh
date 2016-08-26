echo ----- google test start ----
export GTEST_VERSION=master
export GTEST=googletest-${GTEST_VERSION}
wget https://github.com/google/googletest/archive/${GTEST_VERSION}.tar.gz
tar -xzvf *.tar.gz
cd ${GTEST}
mkdir build && cd build
cmake ..
make -j4 && sudo make install && cd ../..
echo ----- google test end ----

echo ----- 1 ----
mkdir dependencies
echo ----- 2 ----
cd dependencies
echo ----- 3 ----
mkdir include
echo ----- 4 ----
cd include
echo ----- 5 ----
wget https://raw.githubusercontent.com/dtschump/CImg/master/CImg.h
wget https://raw.githubusercontent.com/CrikeeIP/Geometry/master/include/geometry/geometry.h
cd ..
echo ----- 6 ----

git clone https://github.com/Dobiasd/FunctionalPlus
cd FunctionalPlus
mkdir build
cd build
cmake ..
sudo make install
cd ..
cd ..

echo ----- 7 ----