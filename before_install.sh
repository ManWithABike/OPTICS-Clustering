mkdir dependencies
cd dependencies
	mkdir include
	cd include
	sudo mkdir /usr/local/include/CImg
	sudo wget https://raw.githubusercontent.com/dtschump/CImg/master/CImg.h -O /usr/local/include/CImg/CImg.h
	printf "--------- CImg downloaded ---------\n"
	sudo mkdir /usr/local/include/geometry
	sudo wget https://raw.githubusercontent.com/CrikeeIP/Geometry/master/include/geometry/geometry.hpp -O /usr/local/include/geometry/geometry.hpp
	printf "--------- Geometry downloaded ---------\n"
	cd ..

	git clone https://github.com/Dobiasd/FunctionalPlus
	cd FunctionalPlus
	mkdir build
	cd build
	cmake ..
	sudo make install
	cd ..
	cd ..
	printf "--------- Fplus installed ---------\n"
cd ..

sudo cp -r ./include/optics /usr/local/include/

cd test
	chmod +x ./compile_and_run.sh
	./compile_and_run.sh
	echo "--------- Tests successfully completed ---------"
cd ..