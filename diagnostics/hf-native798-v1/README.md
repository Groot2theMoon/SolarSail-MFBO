# Whole-sail optics diagnostic commands

Use ordinary Python and NumPy, or the installed Abaqus Python interpreter if that is the available Python/NumPy environment. `sail_optics.py` does not import `odbAccess`. It hashes only the ODBs explicitly named in the study and reads the exported geometry bundles. The exporter is the separate read-only Abaqus command; evaluation never exports or opens an ODB through Abaqus.

```text
python diagnostics/hf-native798-v1/sail_optics.py inventory --manifest STUDY --out-dir NEW_INVENTORY_DIR
python diagnostics/hf-native798-v1/sail_optics.py evaluate --manifest STUDY --out-dir NEW_EVALUATION_DIR
```

Freeze `study.initial.json` with the run, candidate, geometry basis and reviewed facts before linking exports. Save a separate `study.linked.json`: add the SHA-256 reference to the initial manifest and the case geometry links only. Both commands independently recheck the parent, current byte identities and evidence. For a real case, cite the canonical `export_verification.json` as a geometry-basis record; its record ID is arbitrary. The checked record must cover the adapter, vendor runtime, precision, mapping, element types and face-sign paths used. It does not prove historical producer linkage; cite the actual producing attempt and any byte-preserving copy evidence separately. Neither process parses solver logs for scientific meaning.

A qualified real export can be compared with the fixed flat reference even when the control or a required matching fact is unavailable. A missing control stays missing; the paired comparison stays blocked. A contrast completes only when it is eligible and both loads and their finite normalized delta pass at all three attitudes. The screen additionally requires matching clamp hardware, connections and attachment restraints apart from the two applied preloads. Case qualification, matched evidence, and numerical publication are distinct.

The evaluator uses the existing `nominal_square_v1` four rotational copies, global Sun rays and version-1 piecewise-planar faceting. Moments are about the sail origin, not a spacecraft centre of mass. If pressure is omitted, normalized outputs are in m² and m³; physical loads and angles are not inferred. Numerical agreement between integration routes is not a faceting error estimate or physical validation. Optional diagnostic gaps do not prove a settled branch or disqualify an otherwise traceable terminal shape.

Keep all real manifests, geometry bundles, ODBs and solver evidence in private directories. Use fresh output directories for both commands. The committed tests fabricate opaque ODB bytes and in-memory API fixtures; they neither inspect real ODB fields nor establish real provenance.
