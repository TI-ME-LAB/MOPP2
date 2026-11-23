"""Run the bundled step scripts (mopp2/scripts/*.sh) with a config file."""

import logging
import os
import shutil
import subprocess
from pathlib import Path

logger = logging.getLogger("mopp2")

SCRIPTS_DIR = Path(__file__).resolve().parent / "scripts"
CONFIG_TEMPLATE = SCRIPTS_DIR / "config.template.sh"

STEP_SCRIPTS = {
    "host-deplete": "3_host_depletion.sh",
    "drop-trimmed": "drop_trimmed_metaG.sh",
    "subdb:curated": "3a_build_subdb.sh",
    "subdb:reuse": "3b_build_subdb_round2.sh",
    "subdb:globdb": "3c_build_subdb_from_globdb.sh",
    "map": "4_run_mapping.sh",
    "annotate": "5a_build_annotation.sh",
    "blacklist": "5b_make_structural_rna_blacklist.sh",
    "saf": "5c_gen_saf.sh",
    "features": "6_build_features.sh",
}


def script_path(step):
    p = SCRIPTS_DIR / STEP_SCRIPTS[step]
    if not p.is_file():
        raise FileNotFoundError(f"bundled script missing: {p}")
    return p


def run_step(step, config, extra_args=(), conda_switch=False):
    """Run one step script with ``bash``. Returns its exit code."""
    if shutil.which("bash") is None:
        raise RuntimeError("bash not found on PATH")
    config = str(Path(config).resolve())
    cmd = ["bash", str(script_path(step)), config, *extra_args]
    env = os.environ.copy()
    if not conda_switch:
        # The scripts would otherwise `conda activate` MAP_CONDA_ENV / CONDA_ENV
        # (base, features, ...). Under the mopp2 command every tool comes from
        # the environment mopp2 itself is installed in.
        env["MOPP2_NO_CONDA_SWITCH"] = "1"
    logger.info(f"[{step}] $ {' '.join(cmd)}")
    rc = subprocess.call(cmd, env=env)
    if rc != 0:
        logger.error(f"[{step}] exited with code {rc}")
    else:
        logger.info(f"[{step}] finished")
    return rc


def config_value(config, var):
    """Value of ``var`` after sourcing ``config`` in bash ('' when unset)."""
    out = subprocess.run(
        ["bash", "-c", 'set +u; source "$1" >/dev/null 2>&1; printf "%s" "${!2}"', "_",
         str(Path(config).resolve()), var],
        capture_output=True, text=True, check=False,
    )
    return out.stdout


def write_config(outpath, project=None, force=False):
    outpath = Path(outpath)
    if outpath.exists() and not force:
        raise FileExistsError(f"{outpath} exists (use --force to overwrite)")
    text = CONFIG_TEMPLATE.read_text()
    if project:
        text = text.replace("PROJECT=/path/to/project", f"PROJECT={Path(project).resolve()}", 1)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    outpath.write_text(text)
    return outpath
