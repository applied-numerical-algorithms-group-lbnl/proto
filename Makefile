all:
	cd examples/Euler/exec; make
	cd examples/Godunov/exec; make
	cd examples/UnitTests/exec; make
	cd examples/Multigrid/exec; make
	cd examples/Navier/exec; make
	cd examples/Snippets/exec; make
clean:
	cd examples/Euler/exec; make clean
	cd examples/Godunov/exec; make clean
	cd examples/UnitTests/exec; make clean
	cd examples/Multigrid/exec; make clean
	cd examples/Navier/exec; make clean
	cd examples/Snippets/exec; make clean
