import logging
import sys
import time
from pathlib import Path

import click

from mopp2 import __version__
from mopp2 import runner
from mopp2.metadata import MetadataError, generate_metadata, load_metadata
from mopp2.trim import METARS_MAX_LENGTH, METARS_MIN_LENGTH, PAIRED_MIN_LENGTH, trim_files

logger = logging.getLogger("mopp2")
_FMT = logging.Formatter("[%(asctime)s] %(levelname)s %(message)s", "%Y-%m-%d %H:%M:%S")

MSG_WELCOME = """
Welcome to MOPP2 - multi-omics processing pipeline for metagenomic (metaG),
metatranscriptomic (metaT) and metaRibo-Seq (metaRS) data.

Reads are mapped to a per-sample reference database built from the genomes
detected in that sample, and turned into gene- and genome-level count, TPM and
RPKM tables.

\b
  (1) To run every step after trimming:  mopp2 workflow --help
  (2) To run a single step:              mopp2 trim --help, mopp2 map --help, ...

\b
Steps, in order:
  metadata      create / validate the sample metadata
  trim          Trim Galore, per-omic parameters
  host-deplete  remove host reads (host-associated samples only)
  subdb         per-sample reference databases (curated | reuse | globdb)
  map           bowtie2 -> per-sample BAMs
  annotate      prodigal / eggNOG-mapper
  blacklist     ORFs inside rRNA / tRNA genes
  saf           per-sample SAF without those ORFs
  features      featureCounts -> count / TPM / RPKM tables
"""


def _setup_logging(logdir=None, name="mopp2"):
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    sh = logging.StreamHandler()
    sh.setFormatter(_FMT)
    logger.addHandler(sh)
    if logdir is not None:
        Path(logdir).mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(Path(logdir) / f"{name}_{time.strftime('%Y%m%d_%H%M%S')}.log")
        fh.setFormatter(_FMT)
        logger.addHandler(fh)


config_arg = click.argument("config", type=click.Path(exists=True, dir_okay=False))


@click.group(help=MSG_WELCOME, context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, "-v", "--version")
@click.option(
    "--conda-switch", is_flag=True, default=False,
    help="Let the step scripts `conda activate` the envs named in the config "
         "(MAP_CONDA_ENV, FEAT_CONDA_ENV, CONDA_ENV). Off by default: every tool "
         "comes from the active environment.",
)
@click.pass_context
def main(ctx, conda_switch):
    ctx.ensure_object(dict)
    ctx.obj["conda_switch"] = conda_switch


# --------------------------------------------------------------------------- metadata
@main.command(help="""Create and/or validate the sample metadata.

\b
If -i is a directory, every file ending in .fastq.gz / .fq.gz with a
metaG / metaT / metaRS token in its name is listed, and the result is
written to -o (default: <dir>/metadata_auto.tsv) and validated.
If -i is a metadata file, it is validated only.

The metadata has four tab-separated columns: sample_name, identifier, omic
(metaG/metaT/metaRS) and strand (r1/r2). See tests/test_data/metadata_template.tsv.""")
@click.option("-i", "--input", "inp", required=True, type=click.Path(exists=True),
              help="Directory of raw fastq.gz files, or a metadata tsv.")
@click.option("-o", "--outpath", type=click.Path(), help="Output metadata tsv (directory input only).")
@click.option("-p", "--paired-end", is_flag=True, help="Require r1 and r2 for every metaG / metaT library.")
@click.option("-m", "--multiomics", is_flag=True, help="Require metaG, metaT and metaRS for every sample.")
def metadata(inp, outpath, paired_end, multiomics):
    inp = Path(inp)
    if inp.is_dir():
        out = Path(outpath) if outpath else inp / "metadata_auto.tsv"
        _setup_logging(out.parent, "mopp2_metadata")
        try:
            generate_metadata(inp, out)
            load_metadata(out, paired_end=paired_end, multiomics=multiomics, indir=inp)
        except MetadataError as e:
            logger.error(str(e))
            sys.exit(1)
    else:
        if outpath:
            click.echo("note: -o is ignored when validating an existing metadata file", err=True)
        _setup_logging(inp.parent, "mopp2_metadata")
        try:
            load_metadata(inp, paired_end=paired_end, multiomics=multiomics)
        except MetadataError as e:
            logger.error(str(e))
            sys.exit(1)


# --------------------------------------------------------------------------- trim
@main.command(help=f"""Trim adapters and low-quality ends with Trim Galore.

The metadata says which omic each file is, and the parameters follow:

\b
  metaG / metaT  paired-end   --paired --length {PAIRED_MIN_LENGTH}
  metaRS         single-end   --length {METARS_MIN_LENGTH} --max_length {METARS_MAX_LENGTH}

\b
Outputs (named to match the mopp2 config defaults):
  <id>_metaG_R1.fq.gz  <id>_metaG_R2.fq.gz
  <id>_metaT_R1.fq.gz  <id>_metaT_R2.fq.gz
  <id>_metaRS.fq.gz
Reports, FastQC and per-library logs go to <outdir>/trimming_reports/.
Libraries already trimmed are skipped unless --force.""")
@click.option("-i", "--input-dir", required=True, type=click.Path(exists=True, file_okay=False),
              help="Directory holding the raw fastq.gz files.")
@click.option("-o", "--output-dir", required=True, type=click.Path(file_okay=False),
              help="Output directory, e.g. <project>/2.trim.")
@click.option("-m", "--metadata", "md_path", required=True, type=click.Path(exists=True, dir_okay=False),
              help="Sample metadata tsv (see `mopp2 metadata`).")
@click.option("-t", "--threads", default=4, show_default=True, help="Libraries trimmed at the same time.")
@click.option("--cores-per-job", default=1, show_default=True,
              help="Trim Galore --cores for each library (each core uses ~3 processes).")
@click.option("--min-length", default=PAIRED_MIN_LENGTH, show_default=True,
              help="metaG / metaT minimum read length after trimming.")
@click.option("--metars-min-length", default=METARS_MIN_LENGTH, show_default=True,
              help="metaRS minimum read length after trimming.")
@click.option("--metars-max-length", default=METARS_MAX_LENGTH, show_default=True,
              help="metaRS maximum read length after trimming (0 disables).")
@click.option("--no-fastqc", is_flag=True, help="Do not run FastQC on the trimmed reads.")
@click.option("--force", is_flag=True, help="Re-trim libraries whose outputs already exist.")
@click.option("--dry-run", is_flag=True, help="Print the trim_galore commands and exit.")
def trim(input_dir, output_dir, md_path, threads, cores_per_job, min_length,
         metars_min_length, metars_max_length, no_fastqc, force, dry_run):
    _setup_logging(None if dry_run else output_dir, "mopp2_trim")
    try:
        failed = trim_files(
            input_dir, output_dir, md_path,
            threads=threads, cores_per_job=cores_per_job,
            paired_min_length=min_length,
            metars_min_length=metars_min_length,
            metars_max_length=metars_max_length,
            fastqc=not no_fastqc, force=force, dry_run=dry_run,
        )
    except (MetadataError, RuntimeError) as e:
        logger.error(str(e))
        sys.exit(1)
    sys.exit(1 if failed else 0)


# --------------------------------------------------------------------------- config
@main.command("init-config", help="Write a copy of the config template to fill in.")
@click.option("-o", "--output", required=True, type=click.Path(dir_okay=False), help="Config file to write.")
@click.option("--project", type=click.Path(file_okay=False), help="Fill in PROJECT with this directory.")
@click.option("--force", is_flag=True, help="Overwrite an existing file.")
def init_config(output, project, force):
    try:
        p = runner.write_config(output, project=project, force=force)
    except FileExistsError as e:
        raise click.ClickException(str(e))
    click.echo(f"wrote {p}; fill in the database paths and library sections, then run e.g. `mopp2 subdb {p} --mode globdb`")


@main.command("scripts-dir", help="Print the directory holding the bundled step scripts.")
def scripts_dir():
    click.echo(str(runner.SCRIPTS_DIR))


# --------------------------------------------------------------------------- script steps
def _run(ctx, step, config, extra=()):
    _setup_logging()
    sys.exit(runner.run_step(step, config, extra, conda_switch=ctx.obj["conda_switch"]))


@main.command("host-deplete", help="""Remove host reads with bowtie2 against HUMAN_IDX (3_host_depletion.sh).

Reads TRIM_DIR, writes HF_OUT/filtered. Only for host-associated samples
(stool, skin, plant tissue); leave it out for a SynCom or an environmental
cohort. metaRS is left alone while SKIP_METARS_FILTER=1 (the default): ~28 bp
reads hit a 3 Gb host genome by chance.""")
@config_arg
@click.pass_context
def host_deplete(ctx, config):
    _run(ctx, "host-deplete", config)


@main.command("drop-trimmed", help="""Delete trimmed metaG once its host-depleted copy is verified (drop_trimmed_metaG.sh).

A disk-space utility, not a pipeline step. --watch keeps polling while
host-deplete is still running.""")
@config_arg
@click.option("--watch", is_flag=True, help="Poll until every sample is done.")
@click.pass_context
def drop_trimmed(ctx, config, watch):
    _run(ctx, "drop-trimmed", config, ["--watch"] if watch else [])


SUBDB_HELP = """Build the per-sample reference databases.

\b
  --mode curated  3a  a SynCom: two hand-curated genome pools
                      (MAGS_POOL + GLOBDB_POOL)
  --mode reuse    3b  a dropout set: reuse another project's
                      unified_db.syldb (SUBDB_SRC, OUT_DIR)
  --mode globdb   3c  a real community: discover the pool from the reads
                      with sylph against GlobDB (SYLPH_DB_R1, GLOBDB_ROOT),
                      optionally with the project's MAGs, isolates and a fungus

Writes SUBDB_ROOT/unified_db.tsv and SUBDB_ROOT/samples/<sample>/."""


@main.command(help=SUBDB_HELP)
@config_arg
@click.option("--mode", required=True, type=click.Choice(["curated", "reuse", "globdb"]),
              help="Which sub-DB builder to run.")
@click.pass_context
def subdb(ctx, config, mode):
    _run(ctx, f"subdb:{mode}", config)


@main.command("map", help="""Map metaG / metaT / metaRS to each sample's sub-DB with bowtie2 (4_run_mapping.sh).

metaG and metaT are paired-end on default settings. metaRS is single-end and
mapped twice: metaRS_strict (--score-min L,0,-0.01) and metaRS_mopp2
(-L 18 -N 1 --score-min L,0,-0.3, about one mismatch per footprint).""")
@config_arg
@click.pass_context
def map_(ctx, config):
    _run(ctx, "map", config)


@main.command(help="""Predict genes with prodigal and annotate them with eggNOG-mapper (5a_build_annotation.sh).

SKIP_EGGNOG=1 stops after prodigal; everything later steps need is still
written, and a later run with SKIP_EGGNOG=0 adds the functional layer.""")
@config_arg
@click.pass_context
def annotate(ctx, config):
    _run(ctx, "annotate", config)


@main.command(help="""List ORFs that lie inside rRNA (barrnap) or tRNA genes (5b_make_structural_rna_blacklist.sh).

Prodigal calls spurious ORFs inside structural RNA genes, and they draw a large
share of metaRS reads. Needs TRNA_IDX.""")
@config_arg
@click.pass_context
def blacklist(ctx, config):
    _run(ctx, "blacklist", config)


@main.command(context_settings={"ignore_unknown_options": True}, help="""Write per-sample ref.saf without blacklisted ORFs (5c_gen_saf.sh).

Extra options are passed to the script: --cds FILE, --mode auto|filter|rebuild,
--no-backup.""")
@config_arg
@click.argument("extra", nargs=-1, type=click.UNPROCESSED)
@click.pass_context
def saf(ctx, config, extra):
    _run(ctx, "saf", config, list(extra))


@main.command(help="""Count reads per gene with featureCounts and build count / TPM / RPKM tables (6_build_features.sh).

Refuses to run on metaRS without the structural-RNA blacklist.""")
@config_arg
@click.pass_context
def features(ctx, config):
    _run(ctx, "features", config)


# --------------------------------------------------------------------------- workflow
WORKFLOW_STEPS = ["trim", "host-deplete", "subdb", "map", "annotate", "blacklist", "saf", "features"]


@main.command(help="""Run the whole chain from one config.

\b
  trim -> [host-deplete] -> subdb -> map -> annotate -> blacklist -> saf -> features

trim runs only when -i and -m are given; it writes to --trim-dir (default:
TRIM_DIR from the config, else <PROJECT>/2.trim). host-deplete runs only with
--host-depletion. Use --from / --to to run part of the chain; every step skips
work that is already done, so a failed run can be restarted.""")
@config_arg
@click.option("--subdb-mode", required=True, type=click.Choice(["curated", "reuse", "globdb"]),
              help="Sub-DB builder (see `mopp2 subdb --help`).")
@click.option("-i", "--input-dir", type=click.Path(exists=True, file_okay=False), help="Raw fastq directory (enables trim).")
@click.option("-m", "--metadata", "md_path", type=click.Path(exists=True, dir_okay=False), help="Sample metadata tsv (enables trim).")
@click.option("--trim-dir", type=click.Path(file_okay=False), help="Where trimmed reads go.")
@click.option("-t", "--threads", default=4, show_default=True, help="Libraries trimmed at the same time.")
@click.option("--host-depletion", is_flag=True, help="Run host-deplete after trim.")
@click.option("--from", "start", type=click.Choice(WORKFLOW_STEPS), default="trim", show_default=True)
@click.option("--to", "stop", type=click.Choice(WORKFLOW_STEPS), default="features", show_default=True)
@click.pass_context
def workflow(ctx, config, subdb_mode, input_dir, md_path, trim_dir, threads, host_depletion, start, stop):
    _setup_logging()
    steps = WORKFLOW_STEPS[WORKFLOW_STEPS.index(start): WORKFLOW_STEPS.index(stop) + 1]
    if not steps:
        raise click.UsageError("--from comes after --to")
    for step in steps:
        if step == "trim":
            if not (input_dir and md_path):
                if input_dir or md_path:
                    raise click.UsageError("trim needs both -i and -m")
                logger.info("[trim] no -i/-m given, skipped")
                continue
            out = trim_dir or runner.config_value(config, "TRIM_DIR") or (
                runner.config_value(config, "PROJECT") + "/2.trim")
            _setup_logging(out, "mopp2_trim")
            try:
                failed = trim_files(input_dir, out, md_path, threads=threads)
            except (MetadataError, RuntimeError) as e:
                logger.error(str(e))
                sys.exit(1)
            if failed:
                sys.exit(1)
            continue
        if step == "host-deplete" and not host_depletion:
            logger.info("[host-deplete] --host-depletion not set, skipped")
            continue
        name = f"subdb:{subdb_mode}" if step == "subdb" else step
        rc = runner.run_step(name, config, conda_switch=ctx.obj["conda_switch"])
        if rc != 0:
            logger.error(f"workflow stopped at {step}; fix it and resume with --from {step}")
            sys.exit(rc)
    logger.info("workflow finished")


if __name__ == "__main__":
    main()
