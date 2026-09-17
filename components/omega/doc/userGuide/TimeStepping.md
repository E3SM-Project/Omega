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
    ModeSplitShare:
      BtrTimeStepper: Predictor-Corrector
      BtrTimeStep: 0000_00:00:20
      NTimeStepIteration: 2
      NBclCoriolisIteration: 2
      ReinitSplitVelocity: false
    StartTime: 0001-01-01_00:00:00
    StopTime: 0001-01-01_02:00:00
    RunDuration: none
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
| SplitExplicitRK2 | mode-split RK2 with a subcycled barotropic mode |
| UnsplitRK2 | the same RK2 scheme with mode splitting disabled |

`SplitExplicitRK2` and `UnsplitRK2` use the `ModeSplitShare` configuration
subgroup shown above. See [Split time stepping](#omega-user-split-time-stepping)
for descriptions of these schemes and their configuration options.

The time step refers to the main model time step used to advance the solution
forward. The time step is specified as a formatted string and can be provided
in any of the following forms:

- ``DDDD_HH:MM:SS(.sss...)``
- ``HH:MM:SS(.sss...)``
- ``MM:SS(.sss...)``
- ``SS(.sss...)``

Days, hours and minutes are optional but must be in order if included.
Fractional seconds are optional.

The StartTime refers to the starting time for the simulation. It is in the
format ``yyyy-mm-day_hh:mm:ss`` for year, month, day, hour, minute, second.
This refers to the initial start time; for a longer simulation, the current
time will be modified by the restart file to update to the present time for
the current segment of the simulation.

The simulation will be stopped either at a fixed StopTime if provided or
by the RunDuration. For shorter simulations, the StopTime can be used to
specify a specific time to stop. For longer simulations with multiple
segments that are restarted, the RunDuration should used and should be set
to fit within the queue time. The format for StopTime is the same as StartTime.
The format for RunDuration is the same as the TimeStep.

Only one of the StopTime or RunDuration should be specified with the other
set to either an empty string or "none". If both are specified, the
RunDuration is used instead of the StopTime.

```{toctree}
:hidden:

SplitTimeStepping
```
