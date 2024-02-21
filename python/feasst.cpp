#include <pybind11/pybind11.h>
#include <string>
//#include "math/include/position.h"
//#include "system/include/system.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/monte_carlo.h"
//#include "feasst/include/feasst.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    return i + j;
}

void parse(feasst::MonteCarlo * mc, const std::string& line) {
  auto parsed = feasst::parse_line(line, NULL, NULL);
  feasst::arglist list;
  list.push_back(parsed);
  mc->begin(list);
}

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    py::class_<feasst::MonteCarlo>(m, "MonteCarlo")
        .def(py::init<>());
        //.def("begin", &feasst::MonteCarlo::begin);

    m.def("parse", &parse, R"pbdoc(
        Hi
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
