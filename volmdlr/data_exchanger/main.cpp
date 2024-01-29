#include <pybind11/pybind11.h>
#include <string>

#include <pybind11/stl.h> // Needed for conversion between C++ containers and Python objects
#include <map>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)



int add(int i, int j) {
    return i + j;
}

std::string importFile() {
    return "The file was imported";
}

std::map<std::string, std::string> getData() {
    std::map<std::string, std::string> data;
    data["key1"] = "value1";
    data["key2"] = "value2";
    data["key3"] = "value3";
    return data;
}

namespace py = pybind11;

PYBIND11_MODULE(data_exchanger, m) { // (python_example = module name, m == py::module_)
    m.doc() = "API to acess the data exchanger function"; // optional module docstring

    m.def("add", &add, "A function that adds two numbers");
    m.def("import_file", &importFile, "A function that import a file using the data exchanger");
    m.def("get_data", &getData, "A function that returns a dictionary");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}