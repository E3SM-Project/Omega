(omega-dev-logging)=

# Developing Omega Logging

Omega's logging system is built upon the
[spdlog](https://github.com/gabime/spdlog) logging tool.

Logging macros and custom formatters for Kokkos data types have been
integrated into the Omega logging system via the Logging.h header file
located in the src/infra directory. Users who wish to utilize Omega's
logging capabilities must include this header file in their code.

The `src/infra/Logging.cpp` file encompasses various functions related
to logging. As of the current version, this file contains the
`OMEGA::initLogging` function, which initialize the logging process.
It is recommended to invoke the `OMEGA::initLogging` function at the beginning
of an Omega application to properly initialize the logging system.

## Writing log messages

To write log messages, include the Logging.h header file in your code.

```c
#include "Logging.h"
```

You can then use the following macros to write logs of different severity
levels:

```c
LOG_TRACE   (logmsg);
LOG_DEBUG   (logmsg);
LOG_INFO    (logmsg);
LOG_WARN    (logmsg);
LOG_ERROR   (logmsg);
LOG_CRITICAL(logmsg);
```

The "logmsg" argument in these macros can be any string that follows the syntax
of C++'s std::format function. For example, the following code writes a
"Hello, World!" message as a TRACE level log:

```c
LOG_TRACE("Hello, {}!", "World");
```
For more information about the syntax of logmsg, please see the spdlog documentation.

Here is a table that summarizes the different log severities and their meanings:

1. `LOG_CRITICAL`: Used to indicate a critical error or issue that requires
   immediate attention and could potentially cause data loss or application
   failure.
2. `LOG_ERROR`: Used to indicate an error or issue that could cause
   application instability or incorrect behavior.
3. `LOG_WARN`: Used to indicate a warning or potential issue that should be
   addressed but does not necessarily require immediate attention.
4. `LOG_INFO`: Used to provide general informational messages about the
   application's operation or progress.
5. `LOG_DEBUG`: Used to provide detailed debugging information for
   developers to track down issues or analyze behavior.
6. `LOG_TRACE`: Used to provide the most detailed level of information,
   including function calls and variable values, and is typically used for deep
   analysis and debugging.


You can replace the string "World" with any variable whose type is
a C++ basic data type, such as an int, a float, or a string. You can
also use a limited set of Kokkos variables. For example, the following code
writes the value of the variable `MyInt` as a TRACE level log:

```c
int MyInt = 10;
LOG_TRACE("The value of MyInt is: {}", MyInt);
```

## Initializing Omega Logger

The `OMEGA::initLogging` function, located within the `src/infra/Logging.cpp`
file, serves as a pivotal component for initializing the Omega logging system.

The function establishes the default logger configuration, ensuring that
logging messages are effectively saved to a designated file. The path to
this file can be passed as the second argument to `OMEGA::initLogging`. When
it is omitted, the path comes from the `OMEGA_LOG_FILE` environment variable
if that is set and non-empty, and from `OMEGA::OmegaDefaultLogfile`
(`omega.log`) otherwise. The log file is truncated when it is opened, so it
holds only the messages from the current run.

By default, `initLogging` also redirects `std::cout` and `std::cerr` to the log
file on the ranks that log. Pass `false` as the third argument to leave them
alone. The coupled driver does this, since these streams are shared with any
other component on the same rank. In coupled runs they stay in `e3sm.log.*`,
while Omega's own messages go to `ocn.log.*`.

## Log files in the unit tests

Every unit test driver calls `OMEGA::initLogging(DefEnv)` without a file name,
so without any further arrangement all of the tests would write to a single
`omega.log` in the build's `test` directory, with nothing to say which test
wrote which line. To avoid this, `test/CMakeLists.txt` sets `OMEGA_LOG_FILE`
per test in `add_omega_test`, so that each test writes to
`test/logs/<TEST_NAME>.log`. Note that the name comes from the CTest test
name rather than the source file, because several tests are built from the
same source (for example `TEND_PLANE_TEST` and
`TEND_PLANE_SINGLE_PRECISION_TEST` both come from `TendencyTermsTest.cpp`).

A new test added with `add_omega_test` gets its own log file automatically. A
test driver that hard-codes a log file name would opt out of this, so pass no
file name to `initLogging` unless the test manages its own log files, as
`LoggingTest` does.

The environment variable relies on the MPI launcher passing its environment
to the ranks it starts, which is true of `srun`. If a launcher does not, the
tests fall back to writing `omega.log`.

## Creating Logging Macros

The Omega logging macros, denoted by the prefix `LOG_`, are defined within
the `src/infra/Logging.h` file. These macros constitute a main part of
the logging infrastructure, enabling users to seamlessly incorporate
logging functionality into their Omega applications.

Furthermore, the logging framework includes a distinct set of macros that
commence with the prefix `LOGGER_`. These macros offer enhanced
versatility by accommodating the utilization of the spdlog logger as
their first argument. This approach facilitats the integration of various
logger types.

## Customer formatter for Kokkos

Within the same header file, you will encounter specialized spdlog formatter
structs designed to accommodate Kokkos data types.

For further information on customizing the spdlog formatter, refer to
[Custom formatting](https://github.com/gabime/spdlog/wiki/3.-Custom-formatting).
