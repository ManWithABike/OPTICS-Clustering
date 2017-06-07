mkdir dependencies
cd dependencies

	mkdir include
	cd include
	sudo mkdir /usr/local/include/CImg
	sudo wget https://raw.githubusercontent.com/dtschump/CImg/master/CImg.h -O /usr/local/include/CImg/CImg.h
	echo "CImg downloaded"
	sudo mkdir /usr/local/include/geometry
	sudo wget https://raw.githubusercontent.com/CrikeeIP/Geometry/master/include/geometry/geometry.h -O /usr/local/include/geometry/geometry.h
	echo "Geometry downloaded"
	cd ..

	git clone https://github.com/Dobiasd/FunctionalPlus
	cd FunctionalPlus
	mkdir build
	cd build
	cmake ..
	sudo make install
	cd ..
	cd ..
	echo "Fplus installed"

cd ..

	cd test
	chmod +x ./compile_and_run.sh
	./compile_and_run.sh
	echo "Tests successfully completed"

cd..