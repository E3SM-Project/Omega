(omega-dev-howto-add-analysis-group)=

# Add an Analysis Group

Analysis groups bundle operator chains under one configuration section,
such as `GlobalStats` or `MOC`; see {ref}`omega-dev-analysis` for how they
fit into the module. The code is in `src/analysis/`.

To add a new bundled group:

1. **Create the group class** in `analysisGroups/`, deriving from
   `AnalysisGroup`. In the constructor, read the group's configuration
   options, build chain stem strings (field names with spatial operators but
   without temporal operators), and call `buildTemporalChains()` to handle
   temporal operator construction, operator instantiation via
   `parseChainAndBuildOps()`, and automatic `OpChainInfo` population.
   Alternatively, for groups with custom chain building logic, manually call
   `AnalysisManager->parseChainAndBuildOps()` for each chain and record an
   `OpChainInfo` for each chain describing its output frequency and whether it
   performs temporal reduction or instantaneous output.
2. **Create the output streams** by calling `createAnalysisGroupStreams()`
   at the end of the constructor.
3. **Include the header** in `analysisGroups/Groups.h`.
4. **Add a dispatch branch** so the `Analysis` constructor constructs the new
   group when it is enabled in the configuration.
