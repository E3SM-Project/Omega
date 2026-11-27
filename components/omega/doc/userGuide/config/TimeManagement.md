(omega-user-config-time-management)=

# Time Management

The calendar, and when a simulation starts and stops, are set in the
``TimeIntegration`` section of the Omega configuration file, along with the
time stepping options described in {ref}`omega-user-config-time-stepping`:
```yaml
  TimeIntegration:
    CalendarType: No Leap
    StartType: StartUp
    StartTime: 0001-01-01_00:00:00
    StopType: AtTime
    StopCriterion: 0001-01-01_02:00:00
```
In a coupled E3SM case, most of these come from the case and the coupler
instead; see {ref}`omega-user-run-e3sm`.

## Calendar

``CalendarType`` selects the calendar used for tracking simulation time:

| Config option name | Calendar type |
| ------------------ | ------------- |
| Gregorian | typical western calendar with leap years/days |
| No Leap | a Gregorian calendar without leap days |
| Julian | earlier calendar still used by Eastern Orthodox traditions |
| Julian Day | simply tracks the number of days since a reference day |
| Modified Julian Day | a modified version of Julian Day that changes the reference time |
| 360 Day | an artificial calendar with equal 30-day months |
| Custom | a calendar in which a developer can specify the number of days per month and days per year; this cannot yet be configured entirely via the configuration and would require additional code to specify those variables |
| No Calendar | for test cases in which the calendar is meaningless and time is simply tracked from time zero |

## Starting a simulation

The ``StartType`` can be one of five choices. The ``StartUp`` option is for
starting a solution from scratch from an initial state file. The ``Continue``
option is for continuing a simulation from a restart file. The ``Branch``
option will branch from an existing simulation by reading from the restart
file and is identical to ``Continue`` except that some other simulation
metadata (eg case name or simulation name) can be changed. The ``Hybrid``
option will continue a branch from a restart file but will reset the start
time and should only be used for the first step of a hybrid run, similar to the
``StartUp`` option. It is equivalent to using the ``StartUp`` option and
specifying the restart file as an init file. It is only provided in the event
a future mult-level time stepping scheme is used in which starting from a
restart with multiple time levels is an advantage over a fresh start. Finally,
a ``Coupled`` option specifies that the ``StartType`` is provided by the
coupler in a coupled configuration. The actual ``StartType`` used internally
is created by the E3SM CIME tools from a combination of the ``RUN_TYPE`` and
``CONTINUE`` options.

The ``StartTime`` refers to the starting time for the full simulation (not the
current leg of an ongoing simulation). If this is a continuation of an
existing simulation, the current time for this leg of the full simulation will
be set by the restart file.

## Stopping a simulation

A ``StopType`` determines (with the ``StopCriterion`` below) how the simulation
will be stopped. There are three options. The ``AtTime`` option will stop the
simulation at a specific time and the ``StopCriterion`` holds that specific time
as described below. The ``AfterDuration`` option runs the simulation for a fixed
time interval and the ``StopCriterion`` is used to define that interval. A final
option called ``Coupled`` is primarily for coupled simulations where the
simulation will stop after receiving a signal from the coupler.

For the ``AtTime`` stop type, the ``StopCriterion`` must be a time instant. If
the ``StopType`` is ``AfterDuration``, the ``StopCriterion`` is a time
interval (typically ``dddd_hh:mm:ss``). For the ``Coupled`` option the
``StopCriterion`` is ignored.

(omega-user-config-time-formats)=

## Time instants and intervals

A time instant, such as ``StartTime``, is a string of the form
``YYYY-MM-DD_HH:MM:SS(.sss...)`` for year, month, day, hour, minute and
second.

A time interval, such as ``TimeStep``, is a string in one of these forms:

- ``DDDD_HH:MM:SS(.sss...)``
- ``HH:MM:SS(.sss...)``
- ``MM:SS(.sss...)``
- ``SS(.sss...)``

Days, hours and minutes are optional but must be in order if included.
Fractional seconds are optional. The actual width of each unit can be
arbitrary and the separators can be any single non-numeric character.
