/* This is an interface file for swig.
   This file is created by dev/tools/depend.py . Modifications to this
   file will likely be overwritten in the future. Thus, edit depend.py
   instead.

   usage: /path/to/feasst/py/run.sh /path/to/feasst/dev/tools/depend.py -s /path/to/feasst
 */

%module(directors="1") feasst

%feature("director:except") {
  if( $error != NULL ) {
    PyObject *ptype, *pvalue, *ptraceback;
    PyErr_Fetch( &ptype, &pvalue, &ptraceback );
    PyErr_Restore( ptype, pvalue, ptraceback );
    PyErr_Print();
    Py_Exit(1);
  }
}

%warnfilter(509);

%{
#include "utils/include/timer.h"
#include "utils/include/file.h"
#include "utils/include/argument_parse.h"
#include "utils/include/io.h"
#include "utils/include/utils.h"
#include "utils/include/custom_exception.h"
#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/checkpoint.h"
#include "utils/include/serialize.h"
#include "utils/include/progress_report.h"
#include "utils/include/cache.h"
using namespace feasst;
%}
%include "std_string.i"
%include "std_vector.i"
%include "std_shared_ptr.i"
%include "std_iostream.i"
%include "stdint.i"
%template(IntVector) std::vector<int>;
%template(Int2DVector) std::vector<std::vector<int> >;
%template(DoubleVector) std::vector<double>;
%template(Double2DVector) std::vector<std::vector<double> >;
%template(Double3DVector) std::vector<std::vector<std::vector<double> > >;
using namespace std;
%pythonnondynamic;
%include "std_map.i"
%template(args) std::map<std::string, std::string>;
%template(ArgsVector) std::vector<std::map<std::string, std::string> >;
%template(arglist) std::vector<std::pair<std::string, std::map<std::string, std::string> > >;
%shared_ptr(feasst::Timer);
%shared_ptr(feasst::ArgumentParse);
%shared_ptr(feasst::CustomException);
%shared_ptr(feasst::Checkpoint);
%shared_ptr(feasst::ProgressReport);
%shared_ptr(feasst::Cache);
%include utils/include/timer.h
%include utils/include/file.h
%include utils/include/argument_parse.h
%include utils/include/io.h
%include utils/include/utils.h
%include utils/include/custom_exception.h
%include utils/include/debug.h
%include utils/include/arguments.h
%include utils/include/checkpoint.h
%include utils/include/serialize.h
%include utils/include/progress_report.h
%include utils/include/cache.h
