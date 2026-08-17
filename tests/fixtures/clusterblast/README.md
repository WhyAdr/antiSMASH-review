# ClusterBlast Test Fixtures

This directory contains version-faithful fixtures representing raw ClusterBlast and KnownClusterBlast sidecars and JSON documents across antiSMASH generations.

## Directory Structure

- `text/`: Text sidecar format fixtures (`contig_1_c1.txt`)
- `minimal/`: Minimal valid JSON structures testing schema versions in isolation
- `golden/`: Upstream serializer-faithful JSON documents testing integration across versions

## Version Matrix

| Fixture File | antiSMASH Version | Module Container Schema | Result Schema (`GeneralResults`) | Serializer Notes |
|---|---|:---:|:---:|---|
| `text/contig_1_c1.txt` | 6.x–8.x text | N/A | N/A | Standard 6-column tabular text format |
| `minimal/schema1_minimal.json` | 5.x/6.x | 2 | 1 | Schema 1: no `similarity` |
| `minimal/schema2_minimal.json` | 7.0.x | 2 | 2 | Schema 2: data_version present, no `similarity` |
| `minimal/schema3_minimal.json` | 7.1.x | 2 | 3 | Schema 3: added `similarity` |
| `minimal/schema5_minimal.json` | 8.x | 2 | 5 | Schema 5: added Subject `full_name` metadata |
| `minimal/clusterblast_compat_module_schema1.json` | Legacy | 1 | 1 | Module container schema 1 compatibility check |
| `golden/antismash_6_1_1_clusterblast.json` | 6.1.1 | 2 | 1 | Serializer-faithful 6.1.1 output |
| `golden/antismash_7_0_1_clusterblast.json` | 7.0.1 | 2 | 2 | Serializer-faithful 7.0.1 output |
| `golden/antismash_7_1_0_clusterblast.json` | 7.1.0 | 2 | 3 | Serializer-faithful 7.1.0 output |
| `golden/antismash_8_0_4_clusterblast.json` | 8.0.4 | 2 | 5 | Serializer-faithful 8.0.4 output |
