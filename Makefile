	CC = g++

CFLAGS = $(shell root-config --cflags) -std=c++17 -D_GLIBCXX_USE_CXX11_ABI=0
LDFLAGS = $(shell root-config --libs)

.PHONY: build
build:
	cd libPython && CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" python3 setup.py build_ext --inplace

.PHONY: cython-build
cython-build:
	cd libPython && CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" python3 setup.py build_ext --inplace --use-cython

.PHONY: clean
clean:
	rm -f -r libPython/build/
	rm -f libPython/histUtils.so

.PHONY: clean-cpp
clean-cpp:
	rm -f libPython/histUtils.cpp

.PHONY: clean-all
clean-all: clean clean-cpp
