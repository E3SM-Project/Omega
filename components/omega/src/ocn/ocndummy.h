// OCN dummy header file

#include <iostream>
#include "YAKL.h"

using yakl::Array;
using yakl::styleC;
using yakl::memHost;
using yakl::memDevice;
using yakl::c::parallel_for;
using yakl::c::Bounds;

typedef float real;
typedef Array<real,1,memDevice,styleC> real1d;
typedef Array<real,2,memDevice,styleC> real2d;

void die(std::string msg) {
  yakl::yakl_throw(msg.c_str());
}
