# Development Record - 2026-03-14

This file records development-oriented script changes and decisions for the Margie app.

## Objective
- Keep clone-and-setup workflow portable for users.
- Avoid hard-coded repository root paths.
- Preserve local setup writes to `<root-folder>/db` only.
- Support runtime read preference for shared DB path when available.

## Key Principles Applied
1. `<root-folder>` is dynamic:
- Derived from script location (`ROOT_DIR` / `SCRIPT_DIR`).
- No machine-specific absolute repo paths in active runtime scripts.

2. Shared DB precedence for runtime only:
- Prefer `/depot/lindems/data/margie/db` if present and non-empty.
- Fallback to `<root-folder>/db`.

3. Setup write policy:
- `setup_databases.sh` writes only under `<root-folder>/db`.
- Setup never writes into shared depot path.

4. Operon annotation availability:
- Operon calling is treated as a required annotation stage.
- Official runtime image source is `ghcr.io/sajalbhattarai/for_oprn:latest`.

## Scripts Touched During This Development Cycle

### Pipeline/spec docs
- `readme_under_development`
  - Expanded to detailed implementation contract.
  - Added explicit `<root-folder>` path convention note.
  - Added runtime DB precedence and setup write policy section.

### Database setup logic
- `setup_databases.sh`
  - Added dbCAN platform override (`linux/amd64`) via per-db resolver.
  - Fixed dbCAN official image invocation to use `--entrypoint dbcan_build`.
  - Added idempotent `hmmpress` handling by removing stale `.h3*` indices before press.
  - Added dbCAN fallback path using direct file downloads.
  - Added explicit output filenames for `download_file.php` URLs.
  - Added recreation of dbCAN directory after `dbcan_build --clean` removes it.

### Container setup logic
- `setup_containers.sh`
  - Added official operon image pull support for Docker and Apptainer/Singularity using `ghcr.io/sajalbhattarai/for_oprn:latest`.
  - Added pull-first strategy for prebuilt GHCR images for core callers:
    - `prodigal-annotation:1.0`
    - `rasttk-annotation:1.0`
  - Added automatic fallback to local build when prebuilt pull is unavailable.
  - Added environment overrides for forks/custom orgs:
    - `MARGIE_GHCR_OWNER`
    - `MARGIE_GHCR_REPO`
    - `MARGIE_GHCR_PRODIGAL_IMAGE`
    - `MARGIE_GHCR_RASTTK_IMAGE`

### Container publishing CI
- `.github/workflows/publish-core-images.yml`
  - Added GitHub Actions workflow to build and publish multi-arch (`linux/amd64`, `linux/arm64`) core images to GHCR.
  - Publishes:
    - `ghcr.io/<owner>/margie/prodigal-annotation:1.0` and `:latest`
    - `ghcr.io/<owner>/margie/rasttk-annotation:1.0` and `:latest`
  - Triggers on manual dispatch and on pushes that modify prodigal/rasttk container contexts.

### Runtime DB source resolution (execution path)
- `processing/scripts/lib/runtime_db_resolver.sh`
  - Added shared resolver utility used by runtime entry scripts.
  - Resolution policy:
    1) `${MARGIE_SHARED_DB_ROOT:-/depot/lindems/data/margie/db}` when non-empty.
    2) `<root-folder>/db` fallback.

- `processing/containers/rasttk/run_container.sh`
  - Added shared DB preference logic.
  - Added fallback to local `<root-folder>/db`.
  - Added `MARGIE_SHARED_DB_ROOT` env override for deployment flexibility.

- `processing/containers/rasttk/scripts/run_batch_annotations.sh`
  - Replaced machine-specific `MARGIE_DIR` with repo-derived root.
  - Added container-side DB source resolver:
    - `/depot/lindems/data/margie/db` (preferred if populated)
    - `/margie/db` fallback
    - `/margie/databases` compatibility fallback
  - Integrated host-side runtime DB resolver for consistent reporting and policy alignment.

- `processing/containers/rasttk/scripts/run_batch_formatted_annotations.sh`
  - Replaced machine-specific `MARGIE_DIR` with repo-derived root.
  - Integrated host-side runtime DB resolver for consistent policy visibility.

- `app/gui/streamlit_app.py`
  - Runtime DB inventory now follows the same precedence policy (shared depot then local db).

### Gene-caller wrappers aligned to reference implementation
- `processing/scripts/prodigal/preflight_check_prodigal.sh`
  - Added Prodigal execution preflight for host binary, Docker image, or Apptainer/Singularity SIF paths.

- `processing/scripts/prodigal/run_prodigal.sh`
  - Added host-side Prodigal wrapper matching expected output contract:
    - `output/prodigal/<genome>/<genome>.faa`
    - `output/prodigal/<genome>/<genome>.fna`
    - `output/prodigal/<genome>/<genome>.gff`
  - Supports host prodigal, Docker image, and Apptainer/Singularity fallback.

- `processing/scripts/rasttk/run_annotation.sh`
  - Added host-side RASTtk wrapper for gene-calling + annotation flow.
  - Uses runtime DB resolver (shared depot preferred, local fallback).
  - Runs incremental RASTtk container pipeline then generates consolidated `rast.tsv`.

## Paths Guaranteed by Contract
- Databases: `<root-folder>/db` (setup writes here)
- Outputs: `<root-folder>/output`
- Processing scripts/containers: `<root-folder>/processing`

## Notes for Future Merge
- This development record and `readme_under_development` are intended to be merged with top-level `README.md` later into a comprehensive one-stop guide under `docs/`.
