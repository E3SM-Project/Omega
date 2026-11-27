(omega-user-config-time-stepping)=

# Time Stepping
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
    ModeSplitShare:
      BtrTimeStepper: Predictor-Corrector
      BtrTimeStep: 0000_00:00:20
      NTimeStepIteration: 2
      NBclCoriolisIteration: 2
      ReinitSplitVelocity: false
```
This configuration refers to the default time stepping used for the model
dynamics (momentum and continuity equations). Additional time steppers can
be used for other portions of the model (eg the barotropic mode or tracer
transport) that may have different time steps and use different algorithms.

The ``CalendarType``, ``StartType``, ``StartTime``, ``StopType`` and
``StopCriterion`` options are described in
{ref}`omega-user-config-time-management`.

The ``TimeStepper`` option refers to the numerical scheme used to advance the
model in time. Omega implements a number of time-stepping schemes. The user
can select the scheme they want in the configuration file.
The following time steppers are currently available (the schemes are
described in {ref}`omega-tech-time-stepping`):
| Config option name | Scheme |
| ------------------- | ------- |
| Forward-Backward | forward-backward |
| RungeKutta2 | second-order two-stage midpoint Runge Kutta method |
| RungeKutta4 | classic fourth-order four-stage Runge Kutta method |
| SplitExplicitRK2 | mode-split RK2 with a subcycled barotropic mode |
| UnsplitRK2 | the same RK2 scheme with mode splitting disabled |

`SplitExplicitRK2` and `UnsplitRK2` use the `ModeSplitShare` configuration
subgroup shown above. See [Split time stepping](#omega-user-config-split-time-stepping)
for descriptions of these schemes and their configuration options.

The ``TimeStep`` refers to the main model time step used to advance the solution
forward. It is a time interval, in one of the forms given in
{ref}`omega-user-config-time-formats`.
