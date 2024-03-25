# PlantandFoodResearch/pangene: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 0.2.0+dev - [25-Mar-2024]

### `Added`

1. Added changelog and semantic versioning
2. Changed license to MIT
3. Updated `.editorconfig`
4. Moved .literature to test/ branch
5. Renamed `pangene_local` to `local_pangene`
6. Renamed `pangene_pfr` to `pfr_pangene`
7. Added versioning checking
8. Updated github workflow to use pre-commit instead of prettier and editorconfig check
9. Added central singularity cache dir for pfr config
10. Added `SORTMERNA_INDEX` before `SORTMERNA`
11. Fixed sample contamination bug introduced by `file.simpleName`
12. Now using empty files for stub testing in CI
13. Now BRAKER can be skipped by including BRAKER outputs from previous runs in the `target_assemblies` param
14. Added `gffcompare` to merge liftoff annotations

### `Fixed`

### `Dependencies`

1. NextFlow!>=23.04.4

### `Deprecated`
