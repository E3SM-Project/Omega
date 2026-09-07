(omega-user-time-stepping)=

# Time stepping
Time stepper refers to the means by which a simulation is advanced in time.
The configuration of a simulation is set by the ``TimeIntegration`` section
of the Omega configuration file:
```yaml
  TimeIntegration:
    CalendarType: No Leap
    TimeStepper: Forward-Backward
    TimeStep: 0000_00:10:00
    StartType: StartUp
    StartTime: 0001-01-01_00:00:00
    StopType: AtTime
    StopCriterion: 0001-01-01_02:00:00
```
This configuration refers to the default time stepping used for the model
dynamics (momentum and continuity equations). Additional time steppers can
be used for other portions of the model (eg the barotropic mode or tracer
transport) that may have different time steps and use different algorithms.

The Calendar choice is describe in the
[Time Management](#omega-user-time-manager) section.

The TimeStepper option refers to the numerical scheme used to advance the
model in time. Omega implements a number of time-stepping schemes. The user
can select the scheme they want in the configuration file.
The following time steppers are currently available:
| Config option name | Scheme |
| ------------------- | ------- |
| Forward-Backward | forward-backward |
| RungeKutta2 | second-order two-stage midpoint Runge Kutta method |
| RungeKutta4 | classic fourth-order four-stage Runge Kutta method |

The time step refers to the main model time step used to advance the solution
forward. The time step is specified as a formatted string and can be provided
in any of the following forms:

- ``DDDD_HH:MM:SS(.sss...)``
- ``HH:MM:SS(.sss...)``
- ``MM:SS(.sss...)``
- ``SS(.sss...)``

Days, hours and minutes are optional but must be in order if included.
Fractional seconds are optional.

The StartOption can be one of three choices. The StartUp option is for starting
a solution from scratch from an initial state file. The Continue option is for
continuing a simulation from a restart file. The Branch option will branch
from an existing simulation by reading from the restart file, but it will
reset the clock to the StartTime.

The StartTime refers to the starting time for the full simulation (not the
current leg of an ongoing simulation). It is in the
format ``yyyy-mm-day_hh:mm:ss`` for year, month, day, hour, minute, second.
This refers to the initial start time; for a longer simulation, the current
time will be modified by the restart file to update to the present time for
the current segment of the simulation.

A StopType determines (with the StopCriterion below) how the simulation will
be stopped. There are three options. The AtTime option will stop the simulation
at a specific time and the StopCriterion holds that specific time as described
below. The AfterDuration option runs the simulation for a fixed time interval
and the StopCriterion is used to define that interval. A final option called
OnSignal is primarily for coupled simulation where the simulation will stop
and a signal from the coupler.

For the AtTime stop type, the StopCriterion must be a time instant in the
format ``yyyy-mm-dd_hh:mm:ss``. If the StopType is AfterDuration, the
StopCriterion is a time interval in the format described above for the time
step (but typically ``dddd_hh:mm:ss``). For the OnSignal option the
StopCriterion is ignored.
