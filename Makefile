.PHONY: doc
all:
	cd examples/testSuite/exec; make
	cd examples/Euler/exec; make
	cd examples/Godunov/exec; make
	cd examples/Multigrid/exec; make
	cd examples/MultigridLevelData/exec; make
	cd examples/Navier/exec; make
	cd examples/applyKernel/exec; make
	cd examples/forallKernel/exec; make
	cd examples/laplacian/exec; make
	cd examples/leveldata/exec; make
	cd examples/FourthOrder1DWENO/exec; make
clean:
	cd examples/testSuite/exec; make clean
	cd examples/Euler/exec; make clean
	cd examples/Godunov/exec; make clean
	cd examples/Multigrid/exec; make clean
	cd examples/MultigridLevelData/exec; make clean
	cd examples/Navier/exec; make clean
	cd examples/applyKernel/exec; make clean
	cd examples/forallKernel/exec; make clean
	cd examples/laplacian/exec; make clean
	cd examples/leveldata/exec; make clean
	cd examples/FourthOrder1DWENO/exec; make clean

doc:
	 doxygen doxygen.config 
	 @echo see document in doc/html/index.html