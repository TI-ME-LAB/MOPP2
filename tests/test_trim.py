import gzip
import random
import shutil
import subprocess
import sys

import pytest

from mopp2.metadata import load_metadata
from mopp2.trim import build_jobs, trim_files

ADAPTER = "AGATCGGAAGAGC"  # Illumina universal, what Trim Galore auto-detects


def _fastq(path, n, insert_len, read_len, seed):
    rnd = random.Random(seed)
    with gzip.open(path, "wt") as fh:
        for i in range(n):
            insert = "".join(rnd.choice("ACGT") for _ in range(insert_len))
            seq = (insert + ADAPTER + "A" * read_len)[:read_len]
            fh.write(f"@r{i}\n{seq}\n+\n{'I' * read_len}\n")


@pytest.fixture
def dataset(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    _fastq(raw / "S1_metaG_R1_001.fastq.gz", 2000, 120, 150, 1)
    _fastq(raw / "S1_metaG_R2_001.fastq.gz", 2000, 120, 150, 2)
    _fastq(raw / "S1_metaRS_R1_001.fastq.gz", 2000, 28, 75, 3)
    md = tmp_path / "md.tsv"
    md.write_text(
        "sample_name\tidentifier\tomic\tstrand\n"
        "S1_metaG_R1_001.fastq.gz\tS1\tmetaG\tr1\n"
        "S1_metaG_R2_001.fastq.gz\tS1\tmetaG\tr2\n"
        "S1_metaRS_R1_001.fastq.gz\tS1\tmetaRS\tr1\n"
    )
    return raw, md, tmp_path / "2.trim"


def test_build_jobs(dataset):
    raw, md, out = dataset
    jobs = {j.omic: j for j in build_jobs(raw, out, load_metadata(md))}
    g, rs = jobs["metaG"], jobs["metaRS"]
    assert [o.name for o in g.outputs] == ["S1_metaG_R1.fq.gz", "S1_metaG_R2.fq.gz"]
    assert "--paired" in g.command and g.command[g.command.index("--length") + 1] == "20"
    assert [o.name for o in rs.outputs] == ["S1_metaRS.fq.gz"]
    assert "--paired" not in rs.command
    assert rs.command[rs.command.index("--length") + 1] == "18"
    assert rs.command[rs.command.index("--max_length") + 1] == "75"


def _count(path):
    with gzip.open(path, "rt") as fh:
        return sum(1 for _ in fh) // 4


def _lengths(path):
    with gzip.open(path, "rt") as fh:
        return [len(l.strip()) for i, l in enumerate(fh) if i % 4 == 1]


@pytest.mark.skipif(shutil.which("trim_galore") is None, reason="trim_galore not installed")
def test_trim_end_to_end(dataset):
    raw, md, out = dataset
    assert trim_files(raw, out, md, threads=2, fastqc=False) == 0
    r1, r2, rs = out / "S1_metaG_R1.fq.gz", out / "S1_metaG_R2.fq.gz", out / "S1_metaRS.fq.gz"
    assert _count(r1) == _count(r2) > 0
    lens = _lengths(rs)
    assert lens and all(18 <= l <= 75 for l in lens)
    # adapters were removed: the 28 bp inserts come back at ~28 bp
    assert sorted(lens)[len(lens) // 2] <= 30
    assert (out / "trim_summary.tsv").read_text().count("trimmed") == 2
    # a second run skips everything
    assert trim_files(raw, out, md, threads=2, fastqc=False) == 0
    assert (out / "trim_summary.tsv").read_text().count("skipped") == 2


def test_cli_help():
    r = subprocess.run([sys.executable, "-m", "mopp2", "--help"], capture_output=True, text=True)
    assert r.returncode == 0 and "trim" in r.stdout and "workflow" in r.stdout
