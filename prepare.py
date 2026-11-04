#!/usr/bin/env python3
"""SwiftSites: one-time setup for evaluate.py.

    python3 prepare.py            # build the reference, unpack/generate datasets, record expected outputs + timings

Steps:
  1. Build the reference: snp-sites (github.com/sanger-pathogens/snp-sites) pinned to tag v2.5.1 (commit
     657e74d, the current release on bioconda), compiled directly with gcc -O2 (autotools are often missing; the
     only generated file, config.h, is written by hand). zlib is built from GitHub if its headers are missing. The
     binary is checked against upstream's own expected outputs in tests/data/.
  2. Unpack the real alignments (data/) and generate the simulated ones (gen_data.py, fixed seeds).
  3. For each dataset run the reference `snp-sites <args> -o out IN` in a fresh directory and record the sha256 and
     size of every output file it writes (out, out.vcf, out.phylip, out.snp_sites.aln), plus stdout, and the wall
     time (best of 2).
  4. Write bench/manifest.json (expected outputs, timings, environment, script hashes).
Re-running is safe: finished steps are skipped.
"""
import gzip
import hashlib
import json
import lzma
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
EXT = ROOT / "ext"
BENCH = ROOT / "bench"
DATA = BENCH / "data"
TAG, COMMIT = "v2.5.1", "657e74d"
REF = EXT / "snp-sites" / "snp-sites-ref"
SCRIPT_FILES = ["prepare.py", "evaluate.py", "gen_data.py"]
CFLAGS = "-O2"


def sh(cmd, **kw):
    print("+", cmd if isinstance(cmd, str) else " ".join(map(str, cmd)), flush=True)
    return subprocess.run(cmd, shell=isinstance(cmd, str), check=True, **kw)


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def ensure_numpy():
    try:
        import numpy  # noqa: F401
    except ImportError:
        r = subprocess.run([sys.executable, "-m", "pip", "install", "--quiet", "numpy"])
        if r.returncode != 0:  # PEP 668 "externally managed" images
            sh([sys.executable, "-m", "pip", "install", "--quiet", "--break-system-packages", "numpy"])


def zlib_flags():
    EXT.mkdir(exist_ok=True)
    test = EXT / "zt.c"
    test.write_text("#include <zlib.h>\nint main(){return zlibVersion()[0]==0;}\n")
    if subprocess.run(["gcc", str(test), "-lz", "-o", str(EXT / "zt")], capture_output=True).returncode == 0:
        return ""
    print("system zlib headers missing; building zlib from GitHub source")
    if not (EXT / "zlib-install").exists():
        if not (EXT / "zlib").exists():
            sh(["git", "clone", "-q", "--depth", "1", "https://github.com/madler/zlib.git", str(EXT / "zlib")])
        sh("./configure --static --prefix=" + str(EXT / "zlib-install") + " >/dev/null && make -s install",
           cwd=EXT / "zlib")
    return f"-I{EXT/'zlib-install'/'include'} -L{EXT/'zlib-install'/'lib'} "


def collect(outdir, prefix="out"):
    """sha256/size of every file whose name starts with `prefix`, keyed by the suffix after the prefix."""
    res = {}
    for p in sorted(Path(outdir).iterdir()):
        if p.name.startswith(prefix):
            res[p.name[len(prefix):]] = dict(sha256=sha256_file(p), bytes=p.stat().st_size)
    return res


def build_reference(zf):
    if REF.exists():
        return
    d = EXT / "snp-sites"
    if not d.exists():
        sh(["git", "clone", "-q", "https://github.com/sanger-pathogens/snp-sites.git", str(d)])
    sh(["git", "checkout", "-q", TAG], cwd=d)
    head = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=d, capture_output=True, text=True).stdout
    assert head.strip().startswith(COMMIT), head
    (d / "src" / "config.h").write_text('#define PACKAGE_VERSION "2.5.1"\n')
    sh(f"gcc {CFLAGS} {zf}-o snp-sites-ref src/*.c -lz", cwd=d)
    out = subprocess.run([str(REF), "-V"], capture_output=True, text=True).stdout.strip()
    assert out == "snp-sites 2.5.1", out
    # upstream's own expected outputs (tests/data), compared byte for byte
    t = d / "tests" / "data"
    cases = [([], "alignment_file_one_line_per_sequence.aln", {".snp_sites.aln": "alignment_file_one_line_per_sequence.aln.snp_sites.aln",
                                                               ".vcf": "alignment_file_one_line_per_sequence.aln.vcf",
                                                               ".phylip": "alignment_file_one_line_per_sequence.aln.phylip"}),
             ([], "alignment_file_with_n.aln", {".snp_sites.aln": "alignment_file_with_n.aln.snp_sites.aln",
                                                ".vcf": "alignment_file_with_n.aln.vcf",
                                                ".phylip": "alignment_file_with_n.aln.phylip"}),
             (["-c"], "pure_mode_alignment.aln", {".snp_sites.aln": "pure_mode_alignment.aln.snp_sites.aln",
                                                  ".vcf": "pure_mode_alignment.aln.vcf",
                                                  ".phylip": "pure_mode_alignment.aln.phylip"}),
             (["-c", "-b"], "pure_mode_monomorphic_alignment.aln",
              {".snp_sites.aln": "pure_mode_monomorphic_alignment.aln.snp_sites.aln",
               ".vcf": "pure_mode_monomorphic_alignment.aln.vcf",
               ".phylip": "pure_mode_monomorphic_alignment.aln.phylip"})]
    passed, failed = [], []
    with tempfile.TemporaryDirectory() as td:
        for flags, inp, expect in cases:
            for p in Path(td).iterdir():
                p.unlink()
            subprocess.run([str(REF), "-m", "-v", "-p", *flags, "-o", "out", str(t / inp)], cwd=td, check=True,
                           capture_output=True)
            for suf, exp in expect.items():
                ok = (Path(td) / ("out" + suf)).read_bytes() == (t / exp).read_bytes()
                (passed if ok else failed).append(f"{inp} {' '.join(flags)} {suf}")
        r = subprocess.run([str(REF), "-C", "-o", "c.txt", str(t / "small_alignment.aln")], cwd=td, check=True)
        ok = (Path(td) / "c.txt").read_bytes() == (t / "small_alignment.constant_site_counts.txt").read_bytes()
        (passed if ok else failed).append("small_alignment -C")
    (EXT / "upstream_tests.json").write_text(json.dumps(dict(passed=passed, failed=failed), indent=1))
    print(f"reference snp-sites 2.5.1 built; upstream expected outputs: {len(passed)} match, {len(failed)} differ")
    assert len(passed) >= 10, failed


def materialise(name, spec):
    DATA.mkdir(parents=True, exist_ok=True)
    kind = spec["kind"]
    if kind == "real":
        path = DATA / (name + ".fa")
        if not path.exists():
            with lzma.open(ROOT / spec["src"], "rb") as fi, open(path, "wb") as fo:
                shutil.copyfileobj(fi, fo, 1 << 24)
    elif kind == "real_gz":
        path = DATA / (name + ".fa.gz")
        if not path.exists():
            with lzma.open(ROOT / spec["src"], "rb") as fi, \
                    gzip.GzipFile(path, "wb", compresslevel=6, mtime=0) as fo:
                shutil.copyfileobj(fi, fo, 1 << 24)
    elif kind == "file":
        path = DATA / (name + "".join(Path(spec["src"]).suffixes[-2:]))
        if not path.exists():
            shutil.copy(ROOT / spec["src"], path)
    else:
        import gen_data
        path = DATA / (name + (".fa.gz" if spec.get("gz") else ".fa"))
        if not path.exists():
            t0 = time.perf_counter()
            shape = gen_data.write_dataset(spec, path)
            print(f"generated {name} {shape} in {time.perf_counter()-t0:.1f}s", flush=True)
    return path


def run_ref(args, inp):
    with tempfile.TemporaryDirectory(dir=BENCH) as td:
        os.symlink(inp, Path(td) / ("in" + "".join(inp.suffixes[-2:])))
        os.sync()
        t0 = time.perf_counter()
        r = subprocess.run([str(REF), *args, "-o", "out", "in" + "".join(inp.suffixes[-2:])], cwd=td,
                           capture_output=True, check=True)
        dt = time.perf_counter() - t0
        return dt, collect(td), r.stdout


def prepare_sets(specs, kind):
    rows = []
    for name, spec in specs.items():
        path = materialise(name, spec)
        t1, outs, so = run_ref(spec["args"], path)
        t2, outs2, _ = run_ref(spec["args"], path)
        assert outs == outs2, f"reference not deterministic on {name}"
        row = dict(name=name, kind=kind, real=spec["kind"] in ("real", "real_gz", "file"),
                   path=str(path.relative_to(ROOT)), args=spec["args"], input_bytes=path.stat().st_size,
                   expected=outs, expected_stdout_sha256=sha256_bytes(so), ref_s=round(min(t1, t2), 4))
        print(json.dumps(row), flush=True)
        rows.append(row)
    return rows


def main():
    os.chdir(ROOT)
    ensure_numpy()
    sys.path.insert(0, str(ROOT))
    zf = zlib_flags()
    build_reference(zf)
    threads = os.cpu_count() or 1
    manifest_path = BENCH / "manifest.json"
    if manifest_path.exists():
        print("bench/manifest.json exists; delete bench/ to regenerate")
    else:
        import gen_data
        BENCH.mkdir(exist_ok=True)
        dev = prepare_sets(gen_data.DEV_SPECS, "dev")
        hold = prepare_sets(gen_data.HOLDOUT_SPECS, "holdout")
        cpu, flags = "", ""
        try:
            cpu = next(l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name"))
            flags = next(l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("flags"))
        except Exception:
            pass
        manifest = dict(
            created=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            reference=dict(tool="snp-sites", version="2.5.1", tag=TAG, commit=COMMIT, cflags=CFLAGS,
                           cmd="snp-sites <args> -o out IN"),
            upstream_tests=json.loads((EXT / "upstream_tests.json").read_text()),
            threads=threads, cpu=cpu,
            simd=[f for f in ["sse4_2", "popcnt", "avx2", "avx512f", "avx512bw", "bmi2"] if f in flags.split()],
            gcc=subprocess.run(["gcc", "--version"], capture_output=True, text=True).stdout.splitlines()[0],
            python=platform.python_version(),
            dev=dev, holdout=hold,
            scripts_sha256={f: sha256_file(ROOT / f) for f in SCRIPT_FILES},
        )
        manifest_path.write_text(json.dumps(manifest, indent=1))
        print(f"wrote {manifest_path}")
    print("prepare: done")


if __name__ == "__main__":
    main()
