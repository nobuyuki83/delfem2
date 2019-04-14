#Building Instruction


## Linux and Mac OSX

In order to build the DelFEM2 python module, you need to install several packages with the command

	sudo apt install git cmake

download DelFEM2 from github 

	git clone https://github.com/nobuyuki83/delfem2.git

then, download the submodules on which the DelFEM2 depends.

	git submodule update --init --recursive


Next, build the DelFEM2 python module with the following command

    cd module_py
    mkdir buildMake
    cd buildMake
    cmake ..
    make

If everything is is successful, the binary will be placed under  ``delfem2/module_py/buildMake/``. The bindary has name ``dfm2.**.so`` where the ``**`` is the compilation environment. The compiled binary will be automatically replaced and renamed as``delfem2/module_py/dfm2/dfm2.so``.


This is not a recommended way of installation, but you may also use ``pip`` to install the module_py using the ``setup.py``

    pip install . -e


Finally, install several python modules which the DelFEM2 depends on as

	python3 -m pip install numpy PyOpenGL glfw

Try running the python examples in ``delfem2/example_py/``. If the examples work congraturations!! I hope you enjoy the library.
