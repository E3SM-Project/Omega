(omega-dev-howto-add-analysis-operator)=

# Add an Analysis Operator

Analysis operators are the building blocks that analysis groups chain
together; see {ref}`omega-dev-analysis` for how they fit into the module.
The code is in `src/analysis/`.

To add a new operator:

1. **Create the operator class** in `operators/` as a template on the input
   array type, deriving from `AnalysisOperator`. In the constructor, pass the
   operator type name to the base constructor, set `InputNames`, build the
   output name as `<input>_<Suffix>`, set `OutputNames` and `InstanceName`,
   allocate the output data array, and register the output Field with
   `createOutputField()` plus `attachData()`, deriving its units, standard
   name and cell methods from the input as described under
   [Output metadata](#omega-dev-analysis-output-metadata). If the operator
   depends on another operator's output (as `SpatialStdDev` depends on
   `SpatialMean`), append that dependency to `InputNames`.
2. **Implement `compute()`** to retrieve the input array(s) from the Field
   registry, perform the transformation, write the output array(s), and update
   `LastComputed`/`FieldComputed`.
3. **Implement `setPeriodAlarm()`** only if the operator is a temporal
   reduction that finalizes on a period boundary.
4. **Include the header** in `operators/Ops.h`.
5. **Register the operator** by adding a
   `AnalysisOpFactory::registerAllArrayVariants<NewOp>("NewName")` call to
   `Analysis::registerAllBaseAnalysisOperators()` in `operators/Ops.cpp`.
6. **Update the chain parser** if your operator introduces a new naming
   convention. The current parser in `Analysis::parseChainAndBuildOps()`
   identifies operators by prefix (`Spatial` for spatial reductions, `Time`
   for temporal reductions). If your operator uses a different naming
   pattern, add logic to recognize and dispatch the new token. The parser is
   preliminary and expected to evolve as new operator categories are added.
7. **Add unit tests** in `test/analysis/AnalysisOperatorTest.cpp` to verify the
   new operator implementation. Tests should validate correct computation on
   representative input data.

If your operator follows an existing naming convention, no parser changes are
needed. No changes to the dependency resolution, alarm propagation, or compute
logic are required.
