(omega-user-run-logging)=

# Logging

Omega writes its messages, at severities from trace to critical, to a log
file. How to write log messages from code is described in
{ref}`omega-dev-logging`.

## Log file name

By default Omega writes its log messages to `omega.log` in the directory it is
run from. A different name or path can be requested at run time with the
`OMEGA_LOG_FILE` environment variable:

```bash
export OMEGA_LOG_FILE=mylogs/run1.log
```

When more than one MPI rank writes a log (see `OMEGA_LOG_TASKS`), the rank
number is inserted into the name, so `omega.log` becomes `omega_0.log`,
`omega_1.log` and so on. The log file is truncated when Omega starts, so it
always contains only the messages from the current run.

## Standard output and error

In Standalone mode, Omega redirects `std::cout` and `std::cerr` to the log file
on the ranks that log.

(omega-user-run-logging-e3sm)=

## In E3SM

The E3SM driver sets the log file name, so Omega's messages go to `ocn.log.*`
in the run directory and `OMEGA_LOG_FILE` is ignored. `std::cout` and
`std::cerr` are not redirected: other output, including crash backtraces, stays
in `e3sm.log.*`. Check both files when a run fails.
