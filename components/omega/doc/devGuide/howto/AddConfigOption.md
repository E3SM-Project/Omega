(omega-dev-howto-add-config-option)=

# Add a Configuration Option

1. **Add the option and its default** to `configs/Default.yml`, in the section
   of the module that reads it. For example, the horizontal mesh's mesh
   scaling switch:

   ```yaml
     HorzMesh:
       MeshScalingEnable: true
   ```

2. **Read it** where the module initializes, with the `Config` class
   ({ref}`omega-dev-config`). Get the section from the Omega configuration,
   then the option from the section, and abort with a clear message if either
   is missing:

   ```c++
   Config *OmegaConfig = Config::getOmegaConfig();
   Config HorzMeshConfig("HorzMesh");
   Error Err = OmegaConfig->get(HorzMeshConfig);
   CHECK_ERROR_ABORT(Err, "HorzMesh: HorzMesh group not found in configuration");

   bool MeshScalingEnable;
   Err = HorzMeshConfig.get("MeshScalingEnable", MeshScalingEnable);
   CHECK_ERROR_ABORT(Err, "HorzMesh: MeshScalingEnable not found in HorzMesh "
                          "configuration");
   ```

3. **Decide what coupled cases get.** A coupled case uses the `Default.yml`
   value unless the `coupled` section of
   `cime_config/omega_buildnml/data/config_overrides.yaml` overrides it. If the
   case itself must control the option, add it to `BLOCKED_OPTIONS` in
   `cime_config/omega_buildnml/validate.py` so that `user_nl_omega` cannot set
   it. See {ref}`omega-dev-buildnml`, and run the validation it describes.

4. **Document it** on the User Guide page for the feature: add it to the YAML
   example and say what it does and what values it takes. If the meaning
   involves equations, put them on the Technical Guide page and link to it
   ({ref}`omega-dev-docs`).

5. **Tag the pull request** `[NML]`, since it adds a configuration option.
