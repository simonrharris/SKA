/*
 * ska_bindings.cpp
 * Python bindings for ska
 *
 */

// pybind11 headers
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void align(const std::string& files) {
  printf("align not yet implemented");
  std::cerr << files << std::endl;
}

PYBIND11_MODULE(ska, m) {
  m.doc() = "SKA functions";

  // Exported functions
  m.def("align", &align, "Align ska files",
        py::arg("files"));
  // Example args
  /*
        py::arg("db_name"), py::arg("samples"), py::arg("files"),
        py::arg("klist"), py::arg("sketch_size"),
        py::arg("codon_phased") = false, py::arg("calc_random") = true,
        py::arg("use_rc") = true, py::arg("min_count") = 0,
        py::arg("exact") = false, py::arg("num_threads") = 1,
        py::arg("use_gpu") = false, py::arg("device_id") = 0);
  */

  m.attr("version") = VERSION_INFO;
}
