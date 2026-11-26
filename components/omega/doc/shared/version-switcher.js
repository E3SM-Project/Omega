// Add a drop-down menu of published documentation versions to the sidebar.
//
// The versions are listed in versions.json, next to this script, by
// generate_versions_json.py. Their URLs are relative to versions.json, so the
// same files work on GitHub Pages and in a local preview.
(async function () {
  const container = document.getElementById("version-switcher");
  if (!container) return;

  // document.currentScript is only set while the script first runs
  const versionsUrl = new URL("versions.json", document.currentScript.src);

  const meta = document.querySelector('meta[name="doc-version"]');
  const currentVersion = meta ? meta.content : "";

  let versions;
  try {
    const response = await fetch(versionsUrl);
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    versions = await response.json();
  } catch (err) {
    console.warn(`Could not load ${versionsUrl}:`, err);
    return;
  }
  if (!versions.length) return;

  const select = document.createElement("select");
  select.style.marginLeft = "0.5em";
  for (const { version, url } of versions) {
    const option = document.createElement("option");
    option.value = new URL(url, versionsUrl).href;
    option.textContent = version;
    option.selected = version === currentVersion;
    select.appendChild(option);
  }
  select.addEventListener("change", () => {
    window.location.href = select.value;
  });

  const label = document.createElement("label");
  label.textContent = "Version:";
  label.style.color = "white";
  label.appendChild(select);
  container.style.marginTop = "0.5em";
  container.appendChild(label);
})();
