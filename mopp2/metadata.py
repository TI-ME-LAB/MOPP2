"""Sample metadata: create from a directory of fastq files, or load and validate.

The metadata is a tab-separated file with four columns:

    sample_name                          identifier  omic    strand
    S1_metaG_S121_L004_R1_001.fastq.gz   S1          metaG   r1
    S1_metaG_S121_L004_R2_001.fastq.gz   S1          metaG   r2
    S1_metaRS_S13_L004_R1_001.fastq.gz   S1          metaRS  r1

``sample_name`` is the raw fastq file name, ``identifier`` the sample name used
in every downstream file, ``omic`` one of metaG / metaT / metaRS, and ``strand``
r1 or r2 (single-end libraries only have r1).
"""

import csv
import logging
import re
from collections import OrderedDict
from pathlib import Path

logger = logging.getLogger("mopp2")

COLUMNS = ["sample_name", "identifier", "omic", "strand"]
OMICS = ("metaG", "metaT", "metaRS")
STRANDS = ("r1", "r2")
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz")

_OMIC_RE = re.compile(r"(?:^|[_.-])(metaG|metaT|metaRS)(?=$|[_.-])", re.IGNORECASE)
_STRAND_RE = re.compile(r"[_.-](R[12])(?=$|[_.-])")
_IDENT_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")


class MetadataError(ValueError):
    pass


def _strip_fastq_suffix(name):
    for suf in FASTQ_SUFFIXES:
        if name.endswith(suf):
            return name[: -len(suf)]
    return None


def infer_row(filename):
    """Guess (identifier, omic, strand) from a raw fastq file name.

    The omic is the first metaG / metaT / metaRS token, the identifier is
    everything before it, and the strand is the last _R1 / _R2 token after it
    (r1 when there is none). Returns None when the name has no omic token.
    """
    stem = _strip_fastq_suffix(filename)
    if stem is None:
        return None
    m = _OMIC_RE.search(stem)
    if not m:
        return None
    omic = next(o for o in OMICS if o.lower() == m.group(1).lower())
    identifier = stem[: m.start()].rstrip("_.-")
    rest = stem[m.end():]
    strands = _STRAND_RE.findall(rest)
    strand = strands[-1].lower() if strands else "r1"
    return identifier, omic, strand


def generate_metadata(indir, outpath):
    """Scan ``indir`` for fastq files and write a metadata tsv to ``outpath``."""
    indir = Path(indir)
    rows, skipped = [], []
    for f in sorted(indir.iterdir()):
        if not f.is_file() or _strip_fastq_suffix(f.name) is None:
            continue
        row = infer_row(f.name)
        if row is None or not row[0]:
            skipped.append(f.name)
            continue
        rows.append([f.name, *row])
    for name in skipped:
        logger.warning(f"could not infer identifier/omic from '{name}', left out")
    if not rows:
        raise MetadataError(f"no fastq file in {indir} has a metaG/metaT/metaRS token")
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        w.writerows(rows)
    logger.info(f"wrote {len(rows)} rows to {outpath}")
    return outpath


def load_metadata(path, paired_end=False, multiomics=False, indir=None):
    """Load and validate a metadata tsv.

    Returns an ordered dict ``{identifier: {omic: [r1_name] or [r1_name, r2_name]}}``.
    Raises MetadataError listing every problem found.
    """
    path = Path(path)
    errors = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header is None or [h.strip() for h in header] != COLUMNS:
            raise MetadataError(
                f"{path}: header must be {' / '.join(COLUMNS)} (tab-separated), got {header}"
            )
        rows = [r for r in reader if any(c.strip() for c in r)]

    md = OrderedDict()
    seen_files = set()
    for n, row in enumerate(rows, start=2):
        if len(row) != 4:
            errors.append(f"line {n}: expected 4 columns, got {len(row)}")
            continue
        name, ident, omic, strand = (c.strip() for c in row)
        if _strip_fastq_suffix(name) is None:
            errors.append(f"line {n}: '{name}' does not end in .fastq.gz or .fq.gz")
        if name in seen_files:
            errors.append(f"line {n}: '{name}' listed twice")
        seen_files.add(name)
        if not _IDENT_RE.match(ident):
            errors.append(
                f"line {n}: identifier '{ident}' may only contain letters, digits, '_' and '-'"
            )
        if omic not in OMICS:
            errors.append(f"line {n}: omic '{omic}' is not one of {', '.join(OMICS)}")
        if strand not in STRANDS:
            errors.append(f"line {n}: strand '{strand}' must be r1 or r2 (lower case)")
        if indir is not None and not (Path(indir) / name).is_file():
            errors.append(f"line {n}: '{name}' not found in {indir}")
        slot = md.setdefault(ident, OrderedDict()).setdefault(omic, {})
        if strand in slot:
            errors.append(f"line {n}: {ident} {omic} already has a {strand} file ({slot[strand]})")
        slot[strand] = name

    out = OrderedDict()
    for ident, omics in md.items():
        out[ident] = OrderedDict()
        for omic, slot in omics.items():
            if "r1" not in slot:
                errors.append(f"{ident} {omic}: r2 without r1")
                continue
            files = [slot["r1"]] + ([slot["r2"]] if "r2" in slot else [])
            if omic == "metaRS" and len(files) == 2:
                errors.append(f"{ident} metaRS: metaRS is single-end, list only r1")
            if paired_end and omic in ("metaG", "metaT") and len(files) != 2:
                errors.append(f"{ident} {omic}: --paired-end set but there is no r2")
            out[ident][omic] = files
        if multiomics:
            missing = [o for o in OMICS if o not in out[ident]]
            if missing:
                errors.append(f"{ident}: --multiomics set but {', '.join(missing)} missing")

    if errors:
        raise MetadataError(f"{path}: {len(errors)} problem(s)\n  " + "\n  ".join(errors))
    logger.info(
        f"metadata OK: {len(out)} samples, {sum(len(v) for v in out.values())} libraries"
    )
    return out
