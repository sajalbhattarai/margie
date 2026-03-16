import os
import shlex
import subprocess
from pathlib import Path

import streamlit as st

ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = ROOT / "input"
OUTPUT_DIR = ROOT / "output"
DB_DIR = ROOT / "db"

TOOL_DB_MAP = {
    "cog": "cog",
    "eggnog": "eggnog",
    "kegg": "kegg",
    "prodigal": "rasttk",
    "operon": "rasttk",
    "tigrfam": "tigrfam",
    "consolidation": "rasttk",
    "merops": "merops",
    "rasttk": "rasttk",
    "uniprot": "uniprot",
    "dbcan": "dbcan",
    "interpro": "interpro",
    "pfam": "pfam",
    "tcdb": "tcdb",
}

OFFICIAL_TOOL_IMAGES = {
    "dbcan": {
        "docker": ["docker", "pull", "--platform", "linux/amd64", "haidyi/run_dbcan:latest"],
        "apptainer": [
            "apptainer",
            "pull",
            str(ROOT / "processing" / "containers" / "dbcan" / "dbcan.sif"),
            "docker://haidyi/run_dbcan:latest",
        ],
    },
    "eggnog": {
        "docker": [
            "docker",
            "pull",
            "--platform",
            "linux/amd64",
            "quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0",
        ],
        "apptainer": [
            "apptainer",
            "pull",
            str(ROOT / "processing" / "containers" / "eggnog" / "eggnog.sif"),
            "docker://quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0",
        ],
    },
}


def build_command_for_tool(tool: str, runtime: str) -> tuple[list[str], str]:
    official = OFFICIAL_TOOL_IMAGES.get(tool)
    if official:
        return official[runtime], "official upstream image"

    if runtime == "docker":
        return [
            "docker",
            "build",
            "--platform",
            "linux/amd64",
            "-t",
            f"{tool}-annotation:1.0",
            str(ROOT / "processing" / "containers" / tool),
        ], "local repository recipe"

    return [
        "apptainer",
        "build",
        str(ROOT / "processing" / "containers" / tool / f"{tool}.sif"),
        str(ROOT / "processing" / "containers" / tool / f"{tool}.def"),
    ], "local repository recipe"

st.set_page_config(
    page_title="Margie Pipeline GUI",
    layout="wide",
)

st.title("Margie Pipeline GUI")
st.caption("Graphical interface for container and database setup")


def command_exists(cmd: str) -> bool:
    return subprocess.run(["bash", "-lc", f"command -v {cmd}"], capture_output=True).returncode == 0


def run_command(command: list[str], cwd: Path) -> tuple[int, str]:
    proc = subprocess.run(command, cwd=str(cwd), text=True, capture_output=True)
    output = ""
    if proc.stdout:
        output += proc.stdout
    if proc.stderr:
        output += "\n" + proc.stderr
    return proc.returncode, output.strip()


def list_visible_files(folder: Path) -> list[str]:
    if not folder.exists():
        return []
    return sorted([p.name for p in folder.iterdir() if p.is_file() and not p.name.startswith(".")])


def list_visible_dirs(folder: Path) -> list[str]:
    if not folder.exists():
        return []
    return sorted([p.name for p in folder.iterdir() if p.is_dir() and not p.name.startswith(".")])


with st.sidebar:
    st.header("Environment")
    st.write(f"Repository root: {ROOT}")

    docker_ok = command_exists("docker")
    apptainer_ok = command_exists("apptainer")
    singularity_ok = command_exists("singularity")

    st.write("Docker:", "available" if docker_ok else "not found")
    st.write("Apptainer:", "available" if apptainer_ok else "not found")
    st.write("Singularity:", "available" if singularity_ok else "not found")

    st.divider()
    st.write("Setup scripts")
    st.write("setup_containers.sh:", "found" if (ROOT / "setup_containers.sh").exists() else "missing")
    st.write("setup_databases.sh:", "found" if (ROOT / "setup_databases.sh").exists() else "missing")

    st.divider()
    ui_mode = st.radio("Interface mode", ["Beginner", "Advanced"], horizontal=True)


if ui_mode == "Beginner":
    tab_setup, tab_tool, tab_results = st.tabs([
        "Setup",
        "Tool Run",
        "Results",
    ])
else:
    tab_setup, tab_tool, tab_data, tab_run, tab_results = st.tabs([
        "Setup",
        "Tool Run",
        "Database Inventory",
        "Command Runner",
        "Results",
    ])


with tab_setup:
    st.subheader("Container and Database Setup")

    runtime = st.radio("Runtime", ["auto", "docker", "apptainer"], horizontal=True)
    dry_run = st.checkbox("Dry-run", value=True)
    skip_db = st.checkbox("Skip database setup", value=False)

    db_choices = [
        "cog",
        "pfam",
        "tigrfam",
        "dbcan",
        "kegg",
        "eggnog",
        "merops",
        "tcdb",
        "uniprot",
        "interpro",
        "rasttk",
    ]
    selected_db = st.multiselect("Database subset (optional)", options=db_choices, default=[])

    if st.button("Run setup_containers.sh", type="primary"):
        cmd = ["bash", "setup_containers.sh"]

        if runtime == "docker":
            cmd.append("--docker")
        elif runtime == "apptainer":
            cmd.append("--apptainer")

        if skip_db:
            cmd.append("--skip-db")
        elif selected_db:
            cmd.extend(["--db", ",".join(selected_db)])

        if dry_run:
            cmd.append("--dry-run")

        st.write("Command:", " ".join(shlex.quote(x) for x in cmd))
        with st.spinner("Running setup script..."):
            rc, out = run_command(cmd, ROOT)
        st.code(out or "(no output)", language="bash")
        if rc == 0:
            st.success("Completed successfully.")
        else:
            st.error(f"Command failed with exit code {rc}.")

    st.divider()
    st.subheader("Database-only setup")

    db_only = st.multiselect("Databases", options=db_choices, default=["cog", "pfam"])
    db_dry_run = st.checkbox("Dry-run (database-only)", value=True)

    if st.button("Run setup_databases.sh"):
        cmd = ["bash", "setup_databases.sh"]
        if db_only and len(db_only) < len(db_choices):
            cmd.extend(["--db", ",".join(db_only)])
        if db_dry_run:
            cmd.append("--dry-run")

        st.write("Command:", " ".join(shlex.quote(x) for x in cmd))
        with st.spinner("Running database setup..."):
            rc, out = run_command(cmd, ROOT)
        st.code(out or "(no output)", language="bash")
        if rc == 0:
            st.success("Completed successfully.")
        else:
            st.error(f"Command failed with exit code {rc}.")


with tab_tool:
    st.subheader("Individual Tool and Database Run")
    st.caption("Build one container and/or set up its matching database.")

    tool_names = list(TOOL_DB_MAP.keys())
    selected_tool = st.selectbox("Tool", options=tool_names, index=0)
    matching_db = TOOL_DB_MAP[selected_tool]
    st.write(f"Matched database: {matching_db}")
    if selected_tool in OFFICIAL_TOOL_IMAGES:
        st.info(f"{selected_tool} uses the official upstream image in setup and GUI operations.")

    tool_runtime = st.radio("Container runtime", ["docker", "apptainer"], horizontal=True, key="tool_runtime")
    tool_dry_run = st.checkbox("Dry-run (tool run)", value=True, key="tool_dry")

    col_build, col_db, col_both = st.columns(3)

    with col_build:
        if st.button("Build selected container", use_container_width=True):
            cmd, source_label = build_command_for_tool(selected_tool, tool_runtime)

            st.write("Source:", source_label)
            st.write("Command:", " ".join(shlex.quote(x) for x in cmd))
            if tool_dry_run:
                st.code("[DRY-RUN] " + " ".join(shlex.quote(x) for x in cmd), language="bash")
            else:
                with st.spinner("Building container..."):
                    rc, out = run_command(cmd, ROOT)
                st.code(out or "(no output)", language="bash")
                if rc == 0:
                    st.success("Container build completed.")
                else:
                    st.error(f"Container build failed with exit code {rc}.")

    with col_db:
        if st.button("Setup matching database", use_container_width=True):
            cmd = ["bash", "setup_databases.sh", "--db", matching_db]
            if tool_dry_run:
                cmd.append("--dry-run")

            st.write("Command:", " ".join(shlex.quote(x) for x in cmd))
            with st.spinner("Setting up database..."):
                rc, out = run_command(cmd, ROOT)
            st.code(out or "(no output)", language="bash")
            if rc == 0:
                st.success("Database setup completed.")
            else:
                st.error(f"Database setup failed with exit code {rc}.")

    with col_both:
        if st.button("Run both", type="primary", use_container_width=True):
            build_cmd, source_label = build_command_for_tool(selected_tool, tool_runtime)

            db_cmd = ["bash", "setup_databases.sh", "--db", matching_db]
            if tool_dry_run:
                db_cmd.append("--dry-run")

            st.write("Source:", source_label)
            st.write("Build command:", " ".join(shlex.quote(x) for x in build_cmd))
            st.write("DB command:", " ".join(shlex.quote(x) for x in db_cmd))

            if tool_dry_run:
                st.code("[DRY-RUN] " + " ".join(shlex.quote(x) for x in build_cmd), language="bash")
                st.code("[DRY-RUN] " + " ".join(shlex.quote(x) for x in db_cmd), language="bash")
            else:
                with st.spinner("Building container..."):
                    rc1, out1 = run_command(build_cmd, ROOT)
                st.code(out1 or "(no output)", language="bash")
                if rc1 != 0:
                    st.error(f"Container build failed with exit code {rc1}.")
                else:
                    with st.spinner("Setting up database..."):
                        rc2, out2 = run_command(db_cmd, ROOT)
                    st.code(out2 or "(no output)", language="bash")
                    if rc2 == 0:
                        st.success("Container and database setup completed.")
                    else:
                        st.error(f"Database setup failed with exit code {rc2}.")


if ui_mode == "Advanced":
    with tab_data:
        st.subheader("Database Directory Inventory")

        if not DB_DIR.exists():
            st.warning("Database directory not found.")
        else:
            db_subdirs = list_visible_dirs(DB_DIR)
            if not db_subdirs:
                st.info("No database subdirectories found.")
            else:
                for sub in db_subdirs:
                    files = list_visible_files(DB_DIR / sub)
                    with st.expander(f"{sub} ({len(files)} files)"):
                        if files:
                            st.write("\n".join(files))
                        else:
                            st.write("No files present.")


    with tab_run:
        st.subheader("Command Runner")
        st.caption("Use this for advanced commands while staying in the repository root.")

        default_cmd = "bash setup_containers.sh --dry-run"
        user_cmd = st.text_input("Command", value=default_cmd)

        if st.button("Run command"):
            try:
                cmd = shlex.split(user_cmd)
            except ValueError as exc:
                st.error(f"Invalid command syntax: {exc}")
                st.stop()

            with st.spinner("Running command..."):
                rc, out = run_command(cmd, ROOT)
            st.code(out or "(no output)", language="bash")
            if rc == 0:
                st.success("Completed successfully.")
            else:
                st.error(f"Command failed with exit code {rc}.")


with tab_results:
    st.subheader("Input and Output Overview")

    input_files = list_visible_files(INPUT_DIR)
    output_dirs = list_visible_dirs(OUTPUT_DIR)

    col1, col2 = st.columns(2)
    with col1:
        st.markdown("Input files")
        if input_files:
            st.write("\n".join(input_files))
        else:
            st.write("No input files detected.")

    with col2:
        st.markdown("Output tool directories")
        if output_dirs:
            st.write("\n".join(output_dirs))
        else:
            st.write("No output directories with data detected.")

st.divider()
st.caption("Margie GUI focuses on setup, inspection, and command execution with minimal overhead.")
