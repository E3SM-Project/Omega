(omega-user-float-except)=

## Floating-point exceptions

Omega has support for enabling floating-point exceptions. In
certain conditions floating-point exceptions may be enabled by
default in Omega unit tests. Currently, this happens only when
the GNU compiler is used. A CMake option `OMEGA_FPEXCEPTS_IN_TESTS`
can be set to enable floating-point exceptions in unit tests with
any compiler. This enables the following exceptions
- divide by zero
- invalid operation
- overflow

Note that enabling floating-point exceptions in other compilers
may require adjustment to compiler flags. For example, `clang`
has the command-line option  `-ffp-exception-behavior` which
defaults to `ignore`. Enabling floating-point exceptions without
changing this flag to either `maytrap` or `strict` can cause
false positives.
