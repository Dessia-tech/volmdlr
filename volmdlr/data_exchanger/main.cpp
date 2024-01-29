#include <pybind11/pybind11.h>
#include <string>
#include <iostream>
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

std::vector<float> linspace(float start, float end, int num) {
    std::vector<float> result;
    if (num <= 1) {
        result.push_back(start);
        return result;
    }
    float step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result.push_back(start + i * step);
    }
    return result;
}

int factorial(int n) {
    int result = 1;
    while (n > 0) {
        result *= n;
        n--;
    }
    return result;
}

void checkNumber(int number) {
    if (number > 0) {
        std::cout << "The number " << number << " is positive." << std::endl;
    } else if (number < 0) {
        std::cout << "The number " << number << " is negative." << std::endl;
    } else {
        std::cout << "The number is zero." << std::endl;
    }
}

namespace py = pybind11;

PYBIND11_MODULE(data_exchanger, m) { // (data_exchanger = module name, m == py::module_)
    m.doc() = "API to acess the data exchanger function"; // optional module docstring

    m.def("add", &add, "A function that adds two numbers");
    m.def("import_file", &importFile, "A function that import a file using the data exchanger");
    m.def("get_data", &getData, "A function that returns a dictionary");
    m.def("linspace", &linspace, "Generate a linearly spaced array",
          py::arg("start"),
          py::arg("end"),
          py::arg("num"));
    m.def("factorial", &factorial, "Returns the factorial of a given number.");
    m.def("check_number", &checkNumber, "Checks if the number is positive, negative or zero");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}