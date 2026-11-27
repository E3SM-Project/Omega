(omega-user-run-standalone)=

# Running Standalone

:::{admonition} Partly written
:class: warning
This page does not yet walk through setting up and running a standalone case.
Progress is tracked in
[E3SM-Project/Omega#724](https://github.com/E3SM-Project/Omega/issues/724).
:::

Standalone runs are usually set up by
[Polaris](https://docs.e3sm.org/polaris/main/), whose documentation covers
creating, building and running tasks and suites. A standalone run directory
needs:

- the Omega executable, built as described in {ref}`omega-user-build-standalone`;
- the configuration file, `omega.yml`, in the same directory as the
  executable ({ref}`omega-user-config`);
- a mesh file and an initial condition, read by the streams in the
  configuration ({ref}`omega-user-io-mesh`,
  {ref}`omega-user-io-initial-conditions`).

What a run writes is described in {ref}`omega-user-run-logging`,
{ref}`omega-user-run-timing` and {ref}`omega-user-run-errors`.
