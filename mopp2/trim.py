"""Adapter / quality trimming with Trim Galore, driven by the sample metadata.

The omic decides the parameters:

    metaG, metaT  paired-end   trim_galore --paired --length 20
    metaRS        single-end   trim_galore --length 18 --max_length 75

R1 and R2 stay separate files, because later steps map metaG and metaT as
pairs (bowtie2 -1/-2) and sketch them as pairs (sylph -1/-2). Output names
follow the mopp2 config defaults, so ``2.trim/`` can be used as METAG_DIR /
METARS_DIR / TRIM_DIR without renaming:

    <identifier>_metaG_R1.fq.gz   <identifier>_metaG_R2.fq.gz
    <identifier>_metaT_R1.fq.gz   <identifier>_metaT_R2.fq.gz
    <identifier>_metaRS.fq.gz
    <identifier>_<omic>.fq.gz     (a single-end metaG / metaT library)

Trim Galore's reports, FastQC output and one log per library go to
``<outdir>/trimming_reports/``. A library whose final files already exist is
skipped, so an interrupted run can simply be started again.
"""

import csv
import logging
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import List

from mopp2.metadata import load_metadata

logger = logging.getLogger("mopp2")

PAIRED_MIN_LENGTH = 20
METARS_MIN_LENGTH = 18
METARS_MAX_LENGTH = 75


@dataclass
class TrimJob:
    identifier: str
    omic: str
    inputs: List[Path]
    outputs: List[Path]
    command: List[str]
    produced: List[Path] = field(default_factory=list)  # trim_galore's own names
    log: Path = None

    @property
    def label(self):
        return f"{self.identifier} {self.omic}"


def build_jobs(
    indir,
    outdir,
    md,
    paired_min_length=PAIRED_MIN_LENGTH,
    metars_min_length=METARS_MIN_LENGTH,
    metars_max_length=METARS_MAX_LENGTH,
    cores_per_job=1,
    fastqc=True,
    extra_args=(),
):
    indir, outdir = Path(indir), Path(outdir)
    reports = outdir / "trimming_reports"
    jobs = []
    for ident, omics in md.items():
        for omic, files in omics.items():
            inputs = [indir / f for f in files]
            base = f"{ident}_{omic}"
            cmd = ["trim_galore", "--gzip", "--output_dir", str(reports), "--basename", base]
            if cores_per_job > 1:
                cmd += ["--cores", str(cores_per_job)]
            if fastqc:
                cmd.append("--fastqc")
            if omic == "metaRS":
                cmd += ["--length", str(metars_min_length)]
                if metars_max_length:
                    cmd += ["--max_length", str(metars_max_length)]
            else:
                cmd += ["--length", str(paired_min_length)]
            cmd += list(extra_args)
            if len(inputs) == 2:
                cmd += ["--paired", str(inputs[0]), str(inputs[1])]
                produced = [reports / f"{base}_val_1.fq.gz", reports / f"{base}_val_2.fq.gz"]
                outputs = [outdir / f"{base}_R1.fq.gz", outdir / f"{base}_R2.fq.gz"]
            else:
                cmd += [str(inputs[0])]
                produced = [reports / f"{base}_trimmed.fq.gz"]
                outputs = [outdir / f"{base}.fq.gz"]
            jobs.append(
                TrimJob(ident, omic, inputs, outputs, cmd, produced, reports / f"{base}.trim_galore.log")
            )
    return jobs


def _run_job(job):
    for f in job.produced:  # leftovers of an interrupted run
        f.unlink(missing_ok=True)
    logger.info(f"[{job.label}] trimming started")
    with open(job.log, "w") as log:
        log.write("$ " + " ".join(job.command) + "\n\n")
        log.flush()
        rc = subprocess.call(job.command, stdout=log, stderr=subprocess.STDOUT)
    if rc != 0:
        return job, f"trim_galore exited with code {rc}, see {job.log}"
    missing = [str(f) for f in job.produced if not f.is_file()]
    if missing:
        return job, f"trim_galore finished but did not write {', '.join(missing)}, see {job.log}"
    # R1 is moved last: a library counts as done only when its R1 exists.
    for src, dst in reversed(list(zip(job.produced, job.outputs))):
        shutil.move(str(src), str(dst))
    logger.info(f"[{job.label}] trimming finished -> {', '.join(o.name for o in job.outputs)}")
    return job, None


def trim_files(
    indir,
    outdir,
    md_path,
    threads=4,
    cores_per_job=1,
    paired_min_length=PAIRED_MIN_LENGTH,
    metars_min_length=METARS_MIN_LENGTH,
    metars_max_length=METARS_MAX_LENGTH,
    fastqc=True,
    force=False,
    dry_run=False,
    extra_args=(),
):
    """Trim every library in the metadata. Returns the number of failed libraries."""
    if shutil.which("trim_galore") is None and not dry_run:
        raise RuntimeError("trim_galore is not installed or not on PATH")

    md = load_metadata(md_path, indir=indir)
    outdir = Path(outdir)
    (outdir / "trimming_reports").mkdir(parents=True, exist_ok=True)

    jobs = build_jobs(
        indir, outdir, md,
        paired_min_length=paired_min_length,
        metars_min_length=metars_min_length,
        metars_max_length=metars_max_length,
        cores_per_job=cores_per_job,
        fastqc=fastqc,
        extra_args=extra_args,
    )

    todo, results = [], {}
    for job in jobs:
        if not force and all(o.is_file() for o in job.outputs):
            logger.info(f"[{job.label}] already trimmed, skip")
            results[job.label] = "skipped"
        else:
            todo.append(job)

    if dry_run:
        for job in todo:
            print(" ".join(job.command))
        return 0

    failed = 0
    if todo:
        logger.info(f"trimming {len(todo)} libraries, {threads} at a time")
        with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
            for fut in as_completed([pool.submit(_run_job, j) for j in todo]):
                job, err = fut.result()
                if err:
                    failed += 1
                    logger.error(f"[{job.label}] {err}")
                    results[job.label] = "failed"
                else:
                    results[job.label] = "trimmed"

    summary = outdir / "trim_summary.tsv"
    with open(summary, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["identifier", "omic", "status", "inputs", "outputs"])
        for job in jobs:
            w.writerow([
                job.identifier, job.omic, results.get(job.label, "failed"),
                ",".join(p.name for p in job.inputs), ",".join(p.name for p in job.outputs),
            ])
    logger.info(f"summary: {summary}")
    if failed:
        logger.error(f"{failed} of {len(jobs)} libraries failed")
    return failed
