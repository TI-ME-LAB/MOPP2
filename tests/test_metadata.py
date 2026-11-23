import pytest

from mopp2.metadata import MetadataError, generate_metadata, infer_row, load_metadata


@pytest.mark.parametrize("name,expected", [
    ("S1_metaG_S121_L004_R1_001.fastq.gz", ("S1", "metaG", "r1")),
    ("S1_metaG_S121_L004_R2_001.fastq.gz", ("S1", "metaG", "r2")),
    ("119_D0_metaRS_S13_L004_R1_001.fq.gz", ("119_D0", "metaRS", "r1")),
    ("Group_4_1_metaT_sub10_R2.fq.gz", ("Group_4_1", "metaT", "r2")),
    ("S2_metaRS.fastq.gz", ("S2", "metaRS", "r1")),
    ("undetermined_R1.fastq.gz", None),
    ("notes.txt", None),
])
def test_infer_row(name, expected):
    assert infer_row(name) == expected


def _write(tmp_path, rows):
    p = tmp_path / "md.tsv"
    p.write_text("sample_name\tidentifier\tomic\tstrand\n" + "".join("\t".join(r) + "\n" for r in rows))
    return p


def test_load_ok(tmp_path):
    p = _write(tmp_path, [
        ("a_R1.fq.gz", "S1", "metaG", "r1"),
        ("a_R2.fq.gz", "S1", "metaG", "r2"),
        ("b.fq.gz", "S1", "metaRS", "r1"),
    ])
    md = load_metadata(p, paired_end=True)
    assert md == {"S1": {"metaG": ["a_R1.fq.gz", "a_R2.fq.gz"], "metaRS": ["b.fq.gz"]}}


@pytest.mark.parametrize("rows,flags,msg", [
    ([("a.fq.gz", "S1", "metag", "r1")], {}, "omic"),
    ([("a.fq.gz", "S1", "metaG", "R1")], {}, "strand"),
    ([("a.fq.gz", "S.1", "metaG", "r1")], {}, "identifier"),
    ([("a.fq", "S1", "metaG", "r1")], {}, ".fastq.gz"),
    ([("a.fq.gz", "S1", "metaG", "r2")], {}, "r2 without r1"),
    ([("a.fq.gz", "S1", "metaRS", "r1"), ("b.fq.gz", "S1", "metaRS", "r2")], {}, "single-end"),
    ([("a.fq.gz", "S1", "metaG", "r1")], {"paired_end": True}, "no r2"),
    ([("a.fq.gz", "S1", "metaG", "r1")], {"multiomics": True}, "missing"),
])
def test_load_errors(tmp_path, rows, flags, msg):
    with pytest.raises(MetadataError, match=msg):
        load_metadata(_write(tmp_path, rows), **flags)


def test_generate(tmp_path):
    for n in ["S1_metaG_R1_001.fastq.gz", "S1_metaG_R2_001.fastq.gz", "S1_metaRS_R1_001.fastq.gz", "x.txt"]:
        (tmp_path / n).write_bytes(b"")
    out = generate_metadata(tmp_path, tmp_path / "md" / "metadata.tsv")
    md = load_metadata(out, indir=tmp_path)
    assert md["S1"]["metaG"] == ["S1_metaG_R1_001.fastq.gz", "S1_metaG_R2_001.fastq.gz"]
    assert md["S1"]["metaRS"] == ["S1_metaRS_R1_001.fastq.gz"]
