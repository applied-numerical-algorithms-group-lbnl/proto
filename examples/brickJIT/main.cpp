#include <iostream>
#include <fstream>
#include <memory>
#include <cxxabi.h>
#include "Proto_Point.H"
#include "Proto_Brick.H"
#include <Python.h>

typedef void (*Fun)();

int main() {
  std::cout << "Hello, World!" << std::endl;
  std::ofstream fout("test.cpp");
  fout << "#include<iostream>\n#include<string>" << std::endl;
  fout << "const char* msg_lol = \"This is a message\";" << std::endl;
  fout << "extern \"C\" void run() {" << std::endl;
  fout << "std::cout << \"io \" << msg_lol << std::endl;" << std::endl;
  fout << "}" << std::endl;
  int status;
  std::cout << abi::__cxa_demangle(typeid(double).name(), nullptr, nullptr, &status) << std::endl;
  fout.close();
  int ret = system("c++ test.cpp -o test.so -shared -fPIC");
  if (WEXITSTATUS(ret) != EXIT_SUCCESS) {
    std::cout << "compilation failed" << std::endl;
    exit(EXIT_FAILURE);
  }

  sleep(1);

  void *dynlib = dlopen("./test.so", RTLD_LAZY);
  if (!dynlib) {
    std::cerr << "error loading library:\n" << dlerror() << std::endl;
    exit(EXIT_FAILURE);
  }

  void *fun_ptr = dlsym(dynlib, "run");
  const char *dlsym_error = dlerror();
  if (dlsym_error != NULL) {
    std::cerr << "error loading symbol:\n" << dlsym_error << std::endl;
    exit(EXIT_FAILURE);
  }
  {
    auto message_ptr = (char **) dlsym(dynlib, "msg_lol");
    const char *dlsym_error = dlerror();
    if (dlsym_error != NULL) {
      std::cerr << "error loading symbol:\n" << dlsym_error << std::endl;
      exit(EXIT_FAILURE);
    }
    std::cout << *message_ptr << std::endl;
  }

  ((Fun) (fun_ptr))();

  dlclose(dynlib);

  CompilationBase cbase("./");

  std::vector<Proto::Point> pts;

  pts.emplace_back(0, 0, 0);
  pts.emplace_back(1, 0, 0);
  pts.emplace_back(-1, 0, 0);
  pts.emplace_back(0, 1, 0);
  pts.emplace_back(0, -1, 0);
  pts.emplace_back(0, 0, 1);
  pts.emplace_back(0, 0, -1);

  cbase.template getStencil<double>(static_cast<std::vector<Proto::Point> &&>(pts));
  return 0;
}
