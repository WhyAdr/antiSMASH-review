# ClusterBlast test fixtures

The fixture taxonomy separates compact parser contracts from reconstructed upstream
serializer shapes. None of these files is claimed to be an unmodified antiSMASH run.

## Layout

- `text/`: biologically synthetic text-sidecar input.
- `minimal/`: hand-authored parser-contract JSON. These files may omit fields that
  `antiSMASH-review` intentionally ignores.
- `golden/`: serializer-faithful reconstructed ClusterBlast module JSON. Each file
  contains every field emitted for its constructed ClusterBlast objects, embedded in
  the smallest antiSMASH result wrapper accepted by the parser.

## Compatibility matrix

| Fixture | Claimed generation | Module schema | Result schema | Historical distinction |
|---|---|---:|---:|---|
| `text/contig_1_c1.txt` | synthetic 6.x-8.x-style text | N/A | N/A | Six-column pairing table |
| `minimal/schema1_minimal.json` | 5.x/6.x parser contract | 2 | 1 | Similarity and `data_version` absent |
| `minimal/schema2_minimal.json` | 7.0.x parser contract | 2 | 2 | Optional `data_version`; similarity absent |
| `minimal/schema3_minimal.json` | 7.1.x parser contract | 2 | 3 | Similarity serialized |
| `minimal/schema5_minimal.json` | 8.x parser contract | 2 | 5 | Similarity and newer subject representation supported |
| `minimal/clusterblast_compat_module_schema1.json` | synthetic compatibility only | 1 | 1 | No released-version provenance established |

The module-schema-1 fixture only exercises an intentionally tolerated compatibility
shape. It is not evidence that any specific antiSMASH release emitted module schema 1.

## Reconstructed golden provenance

All four golden fixtures were reconstructed from the tagged upstream serializers in
`antismash/modules/clusterblast/results.py`, specifically
`ClusterBlastResults.to_json`, `GeneralResults.to_json`, and
`RegionResult.jsonify`, together with the `ReferenceCluster`, `Protein`, `Score`, and
`Subject` structures in `data_structures.py`.

Common transformations are explicit: identifiers and values are synthetic; one record,
region, ranking, pairing, and relevant reference protein are retained; unrelated
antiSMASH modules and non-ClusterBlast record fields are omitted; optional
`mibig_entries` are absent because the reconstructed objects contain none. The files
are reconstructed fixtures, not authentic antiSMASH output.

| Fixture | Upstream tag | Upstream commit | Source type | Version-specific retained shape | SHA-256 |
|---|---|---|---|---|---|
| `golden/antismash_6_1_1_clusterblast.json` | `6-1-1` | `0933904e9493eede567d65db7e20999a1225ac61` | serializer-faithful reconstruction | Complete schema-1 reference, score, subject, protein, prefix, and search-type fields | `f750d4f7049db154c147aaf67d69bc3426c70a555ff0ecd2d585b3c16ad63b7d` |
| `golden/antismash_7_0_1_clusterblast.json` | `7-0-1` | `323bf70798780f3fbaa2ac2d0d80cd9e293e4271` | serializer-faithful reconstruction | Schema 2 plus `data_version`; no similarity field | `6172e2e82587f3d2a3fd7b5b0291410f819f8e1e0abf79854505782387cd7671` |
| `golden/antismash_7_1_0_clusterblast.json` | `7-1-0` | `b291ad2d4f6b61bb4ed66abd4fd5505e696e02ca` | serializer-faithful reconstruction | Schema 3 serializer-computed similarity | `27752d74ce4b84bda8e204d12d659455c546b734b76c6990f2ed491257426ecd` |
| `golden/antismash_8_0_4_clusterblast.json` | `8-0-4` | `8ce7163dd3c2a64b654a3eca294db07289a93a4a` | serializer-faithful reconstruction | Reference coordinates, `Subject.full_name`, full protein JSON, and similarity | `9971b8bcd4b5b6eb987e2ff5065051eeabc69defe910433368a71d1c9e3aca92` |

`test_serializer_reconstructed_golden_fixtures` verifies both stable normalized parser
fields and each version-dependent raw shape. Update the provenance table and hashes
whenever a golden fixture changes intentionally. `.gitattributes` pins these JSON files
to LF line endings so the recorded byte hashes remain stable across platforms.
