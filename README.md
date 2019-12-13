# Proto + PDFG-IR

This project contains an implementation of the Proto structured grid eDSL in C++
coupled with a Polyhedral+Dataflow graph intermediate representation (IR).

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing
purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

The project requires Git, CMake, and GCC or ICC.

```
$ git --version
git version 2.19.1

$ cmake --version
cmake version 3.12.1

$ g++ --version
g++ (Ubuntu 8.3.0-6ubuntu1~18.10) 8.3.0
```

### Installing

Clone the Proto repository and check out the 'dataflow' branch.

```
$ git clone https://bitbucket.org/berkeleylab/proto.git
$ cd proto
$ git checkout dataflow
```

Enter the Proto directory and clone the PolyhedralDataflowIR repository into 'pdfl':

```
$ git clone https://github.com/BoiseState-AdaptLab/PolyhedralDataflowIR.git ./pdfl
```

Build the polyhedral+dataflow library:

```
$ cd pdfl
$ cmake .
$ make
```

Build the Proto+Dataflow code:

```
$ cd ..
$ cmake .
$ make
```

## Running the tests

Run the test executable:

```
$ ./test/euler_par
```

## Included Packages

* [PolyhedralDataflowIR](https://github.com/BoiseState-AdaptLab/PolyhedralDataflowIR) - Polyhedral+Dataflow IR for dataflow optimizations.
* [CHiLL](https://github.com/CtopCsUtahEdu/chill-dev) - Omega+ calculator performs polyhedral code generation.
* [IEGenLib](https://github.com/CompOpt4Apps/IEGenLib) - Set and relation library with uninterpreted function symbols.
* [ISL](https://github.com/Meinersbur/isl) - Integer Set Library (required by IEGenLib).
* [GMP](https://gmplib.org) - GNU Multiple Precision Arithmetic Library (required by ISL).
* [GoogleTest](https://github.com/google/googletest) - Google Testing and Mocking Framework (for unit tests)

## Authors

* **Brian Van Straalen - *Proto DSL* - [Proto](https://bitbucket.org/berkeleylab/proto.git)
* **Eddie Davis** - *PolyhedralDataflowIR* - [PolyhedralDataflowIR](https://github.com/BoiseState-AdaptLab/PolyhedralDataflowIR)

## License

This project is licensed ... - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Insert NSF grant number...
