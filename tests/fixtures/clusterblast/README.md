# ClusterBlast Test Fixtures

This directory contains version-faithful minimized fixtures representing raw ClusterBlast and KnownClusterBlast sidecars and JSON documents across antiSMASH generations.

## Version Matrix

| Fixture File | antiSMASH Version | Module Container Schema | Result Schema (`GeneralResults`) | Serializer Notes |
|---|---|:---:|:---:|---|
| `contig_1_c1.txt` | 6.x–8.x text | N/A | N/A | Standard 6-column tabular text format |
| `clusterblast_v6_schema1.json` | 6.x | 2 | 1 | `antismash.modules.clusterblast.results.GeneralResults` (schema 1: no `similarity`) |
| `clusterblast_v7_schema3.json` | 7.x | 2 | 3 | `antismash.modules.clusterblast.results.GeneralResults` (schema 3: added `similarity`) |
| `clusterblast_v8_schema5.json` | 8.x | 2 | 5 | `antismash.modules.clusterblast.results.GeneralResults` (schema 5: added Subject `full_name` metadata) |
| `clusterblast_legacy_module_schema1.json` | Legacy/custom | 1 | 1 | Legacy module container schema 1 compatibility check |
