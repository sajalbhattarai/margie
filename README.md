# Margie: Mostly Automated Rapid Genome Inference Environment

The development of this pipeline has been possible due to combined efforts of following scientists:

**Sajal Bhattarai**, Department of Food Science, Purdue University (developer)


**Dane Deemer**, Rosen Advanced Computing CLuster, Purdue University (developer)


**Stephen R Lindemann**, Department of Food Science, Purdue University (Supervisor)

## Overview

Margie is a containerized, high-throughput functional annotation pipeline for bacterial genomes. It integrates 13 complementary annotation tools to provide context-rich gene and protein characterization from a single input genome in FASTA format. The pipeline is designed for reproducibility across workstations and HPC clusters through Docker and Apptainer (Singularity) containerization.

The pipeline performs gene calling via BV-BRC RASTtk and functional annotation against COG, Pfam, TIGRfam, KEGG, eggNOG, UniProt, TCDB, MEROPS, dbCAN, and InterPro databases. Results from all tools are consolidated into a unified, gene-centric output table with functional annotations and confidence scores.

---

## Annotation Tools and Databases

| Tool | Database | Method | Primary Function | Reference |
|------|----------|--------|-----------------|-----------|
| RASTtk | BV-BRC | Rule-based | Gene calling, initial annotation | Brettin, T., Davis, J. J., Disz, T., Edwards, R. A., Gerdes, S., Olsen, G. J., Olson, R., Overbeek, R., Parrello, B., Pusch, G. D., Shukla, M., Thomason, J. A., III, Stevens, R., Vonstein, V., Wattam, A. R., & Xia, F. (2015). RASTtk: A modular and extensible implementation of the RAST algorithm for building custom annotation pipelines and annotating batches of genomes. *Scientific Reports*, *5*, 8365. https://doi.org/10.1038/srep08365 |
| Prodigal | None (ab initio) | Ab initio ORF prediction | Alternative prokaryotic gene calling | Hyatt, D., Chen, G.-L., LoCascio, P. F., Land, M. L., Larimer, F. W., & Hauser, L. J. (2010). Prodigal: Prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics*, *11*, 119. https://doi.org/10.1186/1471-2105-11-119 |
| COGclassifier (RPS-BLAST) | COG 2024 + CDD | Official COGclassifier workflow | Orthologous group assignment | Galperin, M. Y., Wolf, Y. I., Makarova, K. S., Vera Alvarez, R., Landsman, D., & Koonin, E. V. (2021). COG database update: Focus on microbial diversity, model organisms, and widespread pathogens. *Nucleic Acids Research*, *49*(D1), D274–D281. https://doi.org/10.1093/nar/gkaa1018 |
| HMMER | Pfam 36.0 | Profile HMM | Protein family and domain annotation | Mistry, J., Chuguransky, S., Williams, L., Qureshi, M., Salazar, G. A., Sonnhammer, E. L. L., Tosatto, S. C. E., Paladin, L., Raj, S., Richardson, L. J., Finn, R. D., & Bateman, A. (2021). Pfam: The protein families database in 2021. *Nucleic Acids Research*, *49*(D1), D412–D419. https://doi.org/10.1093/nar/gkaa913 |
| HMMER | TIGRfam 15.0 | Profile HMM | Curated protein family annotation | Haft, D. H., Selengut, J. D., & White, O. (2003). The TIGRFAMs database of protein families. *Nucleic Acids Research*, *31*(1), 371–373. https://doi.org/10.1093/nar/gkg128 |
| KofamScan | KEGG KOfam | Profile HMM | KEGG Orthology assignment | Aramaki, T., Blanc-Mathieu, R., Endo, H., Ohkubo, K., Kanehisa, M., Goto, S., & Ogata, H. (2020). KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold. *Bioinformatics*, *36*(7), 2251–2252. https://doi.org/10.1093/bioinformatics/btz859 |
| DIAMOND | eggNOG 5.0 | Sequence similarity | Orthology and pathway annotation | Huerta-Cepas, J., Szklarczyk, D., Heller, D., Hernández-Plaza, A., Forslund, S. K., Cook, H., Mende, D. R., Letunic, I., Rattei, T., Jensen, L. J., von Mering, C., & Bork, P. (2019). eggNOG 5.0: A hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. *Nucleic Acids Research*, *47*(D1), D309–D314. https://doi.org/10.1093/nar/gky1085 |
| UPIMAPI + Swiss-Prot | UniProt Swiss-Prot | Official UniProt ID mapping plus Swiss-Prot sequence search | Functional protein annotation | The UniProt Consortium. (2023). UniProt: The Universal Protein Database in 2023. *Nucleic Acids Research*, *51*(D1), D523–D531. https://doi.org/10.1093/nar/gkac1052 |
| BLAST+ | TCDB | Sequence similarity | Transporter classification | Saier, M. H., Jr., Reddy, V. S., Moreno-Hagelsieb, G., Hendargo, K. J., Zhang, Y., Iddamsetty, V., Lam, K. J. K., Tian, N., Russum, S., Wang, J., & Medrano-Soto, A. (2021). The Transporter Classification Database (TCDB): 2021 update. *Nucleic Acids Research*, *49*(D1), D461–D467. https://doi.org/10.1093/nar/gkaa1004; Saier, M. H., Jr., Tran, C. V., & Barabote, R. D. (2006). TCDB: the Transporter Classification Database for membrane transport protein analyses and information. *Nucleic Acids Research*, *34*(Database issue), D181-D186. https://doi.org/10.1093/nar/gkj001 |
| DIAMOND | MEROPS | Sequence similarity | Peptidase family classification | Rawlings, N. D., Barrett, A. J., Thomas, P. D., Huang, X., Bateman, A., & Finn, R. D. (2018). The MEROPS database of proteolytic enzymes, their substrates and inhibitors in 2017 and a comparison with peptidases in the PANTHER database. *Nucleic Acids Research*, *46*(D1), D624–D632. https://doi.org/10.1093/nar/gkx1134; Rawlings, N. D., Waller, M., Barrett, A. J., & Bateman, A. (2014). MEROPS: the database of proteolytic enzymes, their substrates and inhibitors. *Nucleic Acids Research*, *42*(D1), D503-D509. https://doi.org/10.1093/nar/gkt953 |
| HMMER | dbCAN | Profile HMM | Carbohydrate-active enzyme annotation | Zhang, H., Yohe, T., Huang, L., Entwistle, S., Wu, P., Yang, Z., Busk, P. K., Xu, Y., & Yin, Y. (2018). dbCAN2: A meta server for automated carbohydrate-active enzyme annotation. *Nucleic Acids Research*, *46*(W1), W95–W101. https://doi.org/10.1093/nar/gky418 |
| Operon (custom; adapted from UniOP) | RAST GFF + protein order | Rule-based gene-neighborhood inference | Operon boundary and membership prediction | Salgado, H., Moreno-Hagelsieb, G., Smith, T. F., & Collado-Vides, J. (2000). Operons in *Escherichia coli*: Genomic analyses and predictions. *Proceedings of the National Academy of Sciences*, *97*(12), 6652–6657. https://doi.org/10.1073/pnas.110147297; Westover, B. P., Buhler, J. D., Sonnenburg, J. L., & Gordon, J. I. (2005). Operon prediction without a training set. *Bioinformatics*, *21*(7), 880–888. https://doi.org/10.1093/bioinformatics/bti123; Su, H., Zhang, R., & Söding, J. (2024). UniOP: A universal operon prediction for high-throughput prokaryotic (meta-)genomic data using intergenic distance. *bioRxiv*. https://doi.org/10.1101/2024.11.11.623000; hongsua. (n.d.). *UniOP* [Computer software]. GitHub. Retrieved March 13, 2026, from https://github.com/hongsua/UniOP |
| Mapping | InterPro | Cross-reference | Integrated domain annotation | Jones, P., Binns, D., Chang, H.-Y., Fraser, M., Li, W., McAnulla, C., McWilliam, H., Maslen, J., Mitchell, A., Nuka, G., Pesseat, S., Quinn, A. F., Sangrador-Vegas, A., Scheremetjew, M., Yong, S.-Y., Lopez, R., & Hunter, S. (2014). InterProScan 5: Genome-scale protein function classification. *Bioinformatics*, *30*(9), 1236–1240. https://doi.org/10.1093/bioinformatics/btu031 |
| Consolidation | All above | Integration | Unified gene-centric output table | — |

---

## Requirements

**Local workstation (macOS or Linux):**
- Docker Desktop 4.25 or later
- At least 16 GB RAM (32 GB recommended)
- At least 120 GB free disk space for databases
- Python 3.9+ (used for helper scripts and fallback database downloader)

**HPC cluster:**
- Apptainer 1.1+ or Singularity 3.8+
- SLURM or equivalent workload manager
- At least 120 GB available storage for databases
- Python 3.9+

**Both environments require:**
- Internet access for initial database downloads
- Host installation of `diamond`, `makeblastdb`, and `hmmpress` is optional (indexing can run inside tool containers)

---

## Installation

### 1. Clone the repository

```
git clone https://github.com/sajalbhattarai/margie.git
cd margie
```

### 2. Build containers and download databases

On a local workstation with Docker, or on an HPC node with Apptainer loaded:

```
./setup_containers.sh
```

Optional: install direct commands to `~/bin` (default), so you can run without `./`:

```
./install_commands.sh
```

Then you can use:

```
setup_containers --dry-run --cog --kegg --pfam
setup_databases --cog --pfam
annotate-genome --interactive
```

The script auto-detects the available runtime, builds all containers, and downloads databases into `db/`.
If host `diamond`/`makeblastdb`/`hmmpress` are missing, indexing falls back to the corresponding tool containers.
If host `wget` is missing, downloads fall back to a built-in Python downloader with progress output.
For `dbcan` and `eggnog`, Margie uses the official upstream container images during setup and GUI operations rather than local custom image builds.
For `prodigal` and `rasttk`, setup now attempts to pull prebuilt GHCR images first, then falls back to local Docker build automatically if pull is unavailable.

Default prebuilt image locations:

```
ghcr.io/sajalbhattarai/margie/prodigal-annotation:1.0
ghcr.io/sajalbhattarai/margie/rasttk-annotation:1.0
```

Override prebuilt source for forks or private org mirrors:

```
export MARGIE_GHCR_OWNER=<your-owner>
export MARGIE_GHCR_REPO=<your-repo>
# optional full-image overrides:
export MARGIE_GHCR_PRODIGAL_IMAGE=ghcr.io/<owner>/<repo>/prodigal-annotation:1.0
export MARGIE_GHCR_RASTTK_IMAGE=ghcr.io/<owner>/<repo>/rasttk-annotation:1.0
```

Maintainers can publish/update these prebuilt images using the workflow:

```
.github/workflows/publish-core-images.yml
```

To specify explicitly:

```
./setup_containers.sh --docker      # force Docker
./setup_containers.sh --apptainer   # force Apptainer
./setup_containers.sh --skip-db     # build containers only
./setup_containers.sh --operon      # build only the operon container path
./setup_containers.sh --prodigal --rasttk  # build only selected tool containers
./setup_containers.sh --db cog,pfam # build containers + selected databases
./setup_containers.sh --cog --kegg --pfam  # build only those tool containers and set up those DBs
./setup_containers.sh --all         # all supported databases
./setup_containers.sh --dry-run     # preview all build/download commands
./setup_containers.sh --preflight   # environment check only
```

### 2.1 Guided interactive setup (recommended for first-time users)

```
./annotate-genome --interactive
```

The interactive guide provides:
- Step-by-step beginner onboarding from a fresh clone
- Preflight check with install guidance when requirements are missing
- Runtime selection (auto, Docker, or Apptainer)
- Database profile selection (minimal, recommended, all, or custom)
- Automatic dry-run preview followed by optional full setup
- Input folder readiness check and optional GUI launch

On HPC, load the Apptainer module before running:

```
module load apptainer
./setup_containers.sh --apptainer
```

Apptainer `.sif` files are written to each container's directory under `processing/containers/`. They are excluded from git by `.gitignore`.

Each setup run now also writes a dated provenance bundle under `logs/container_setup/` using the prefix `container-setup-<YYYYMMDD-HHMMSS>`. That bundle includes the full console transcript (`.log`), a machine-readable status ledger with phase, artifact, status, detail, and command columns (`.status.tsv`), and a short human-readable conclusion summary (`.summary.txt`) listing which containers/images succeeded or failed.

### 3. Database-only mode (optional)

If you want to download databases without rebuilding containers:

```
./setup_databases.sh
./setup_databases.sh --db cog,pfam,tigrfam
./setup_databases.sh --cog --pfam --tigrfam
./setup_databases.sh --all
./setup_databases.sh --dry-run
```

Total database footprint is approximately 60-90 GB depending on which UniProt release is used. MEROPS files are pulled from the public EBI FTP release, and the pipeline requires both DIAMOND (`pepunit.lib`) and HMMER (`meropsscan.lib`) indexes for MEROPS. COG setup now follows official COGclassifier resources (`cddid.tbl.gz`, `Cog_LE.tar.gz`, `cog-24.def.tab`, `cog-24.fun.tab`) under `db/cog`, dbCAN setup requires the official `haidyi/run_dbcan:latest` build path, eggNOG setup uses the official `quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0` runtime with bacteria (taxid 2) and archaea (taxid 2157) members/annotations, KOfam uses the official Genome.jp INSTALL tarball path, TCDB uses the official BLAST+ database (`tcdb_blast`) paper-backed local CLI workflow, and UniProt uses Swiss-Prot-only resources under `db/uniprot` with official UPIMAPI (`upimapi=1.13.3`) for annotation and UniProt ID mapping.

### 4. Optional graphical interface

Margie includes an optional Streamlit GUI for users who prefer a graphical workflow.

Install GUI dependencies:

```
pip install -r requirements-gui.txt
```

Launch GUI:

```
./run_gui.sh
```

The GUI provides:
- Container and database setup with dry-run support
- Individual tool run mode (build one container and set up matching database)
- Beginner mode (simplified navigation)
- Advanced mode (database inventory and command runner)
- Database inventory viewer
- Command runner for advanced usage
- Input/output overview panel

---

## Usage

Recommended first-time workflow:

```
./annotate-genome --interactive
```

Preflight check only:

```
./annotate-genome --preflight
```

Open advanced menu:

```
./annotate-genome --menu
```

Non-interactive setup:

```
./annotate-genome --setup --docker
./annotate-genome --setup --apptainer --db cog,pfam
./annotate-genome --setup --dry-run
./annotate-genome --docker --cog --kegg --pfam
./annotate-genome --all
```

Output is written to `output/` with one subdirectory per annotation tool and a final consolidated table in `output/consolidated/`.

---

## Repository Structure

```
margie/
├── setup_containers.sh          # Builds all containers (Docker or Apptainer)
├── setup_databases.sh           # Downloads and indexes all databases
├── processing/
│   └── containers/              # One subdirectory per annotation tool
│       ├── rasttk/
│       │   ├── Dockerfile       # Docker build recipe
│       │   ├── rasttk.def       # Apptainer build recipe
│       │   └── scripts/         # Tool-specific annotation scripts
│       ├── cog/
│       ├── pfam/
│       ├── tigrfam/
│       ├── kegg/
│       ├── eggnog/
│       ├── uniprot/
│       ├── tcdb/
│       ├── merops/
│       ├── dbcan/
│       ├── interpro/
│       ├── prodigal/
│       ├── operon/
│       └── consolidation/
├── db/                          # Databases (downloaded by setup_databases.sh)
├── input/                       # Input genome FASTA files
├── output/                      # Annotation results
└── logs/                        # Execution logs
```

Databases, built container images, input genomes, and output files are excluded from version control via `.gitignore`.

---

## HPC Deployment

If building containers locally and transferring to an HPC cluster:

1. Build Docker images locally with `./setup_containers.sh --docker`
2. Export images to `.tar` archives and transfer them to the cluster
3. On the cluster, convert each archive to Apptainer format:

```
module load apptainer
apptainer build <tool>.sif docker-archive://<tool>.tar
```

Alternatively, build Apptainer images directly on the HPC from the `.def` files in each container directory, which requires no Docker:

```
module load apptainer
apptainer build processing/containers/pfam/pfam.sif processing/containers/pfam/pfam.def
```

---

## Citation

If you use this pipeline in your research, please cite the underlying tools and databases:

- Brettin T et al. (2015) RASTtk: A modular and extensible implementation of the RAST algorithm for building custom annotation pipelines and annotating batches of genomes. *Scientific Reports* 5:8365.
- Buchfink B, Xie C, Huson DH (2015) Fast and sensitive protein alignment using DIAMOND. *Nature Methods* 12:59-60.
- Eddy SR (2011) Accelerated Profile HMM Searches. *PLoS Computational Biology* 7(10):e1002195.
- Galperin MY et al. (2021) COG database update: focus on microbial diversity, model organisms, and widespread pathogens. *Nucleic Acids Research* 49:D274-D281.
- Mistry J et al. (2021) Pfam: The protein families database in 2021. *Nucleic Acids Research* 49:D412-D419.
- Kanehisa M & Goto S (2000) KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research* 28:27-30.
- Huerta-Cepas J et al. (2019) eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource. *Nucleic Acids Research* 47:D309-D314.
- Saier MH et al. (2021) The Transporter Classification Database (TCDB): 2021 update. *Nucleic Acids Research* 49:D461-D467.
- Saier MH Jr, Tran CV, Barabote RD (2006) TCDB: the Transporter Classification Database for membrane transport protein analyses and information. *Nucleic Acids Research* 34:D181-D186.
- Rawlings ND et al. (2018) The MEROPS database of proteolytic enzymes, their substrates and inhibitors in 2017. *Nucleic Acids Research* 46:D624-D632.
- Rawlings ND, Waller M, Barrett AJ, Bateman A (2014) MEROPS: the database of proteolytic enzymes, their substrates and inhibitors. *Nucleic Acids Research* 42:D503-D509.
- Zhang H et al. (2018) dbCAN2: a meta server for automated carbohydrate-active enzyme annotation. *Nucleic Acids Research* 46:W95-W101.
- Su, H., Zhang, R., & Söding, J. (2024). UniOP: A universal operon prediction for high-throughput prokaryotic (meta-)genomic data using intergenic distance. *bioRxiv*. https://doi.org/10.1101/2024.11.11.623000.
- hongsua. (n.d.) *UniOP* [Computer software]. GitHub. Retrieved March 13, 2026, from https://github.com/hongsua/UniOP.
- Jones P et al. (2014) InterProScan 5: genome-scale protein function classification. *Bioinformatics* 30:1236-1240.
- The UniProt Consortium (2023) UniProt: the Universal Protein Database in 2023. *Nucleic Acids Research* 51:D523-D531.

---

## License

MIT License. See LICENSE for details.
