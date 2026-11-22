"""Download the raw data and check that it is exactly what we expect.

Reads   metadata/sources.tsv   where each file comes from and its SHA-256
Writes  data/raw/...           the files listed there, made read-only (raw data is never edited)

Run it again at any time: files that are already present and correct are not downloaded
again. Use --verify-only to check without downloading. Works on Linux, macOS and Windows
(only the Python standard library).

    python scripts/download_data.py [--verify-only]
"""

import argparse
import csv
import hashlib
import http.client
import os
import stat
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCES = ROOT / "metadata" / "sources.tsv"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url, dest, attempts=3):
    """Download url to dest (through a temporary file, so a failed download leaves no stub)."""
    print(f"  downloading {url}", flush=True)
    request = urllib.request.Request(url, headers={"User-Agent": "curl/8 (cbio7501 nmtf-project)"})
    for attempt in range(1, attempts + 1):
        tmp = None
        try:
            with urllib.request.urlopen(request, timeout=60) as response, \
                    tempfile.NamedTemporaryFile(dir=dest.parent, delete=False) as out:
                tmp = Path(out.name)
                while chunk := response.read(1 << 20):
                    out.write(chunk)
            os.replace(tmp, dest)
            return
        except (OSError, http.client.HTTPException) as err:  # network hiccup: wait and retry
            if tmp is not None and tmp.exists():
                tmp.unlink()
            if attempt == attempts:
                raise
            print(f"  attempt {attempt} failed ({err}); retrying in {5 * attempt} s", flush=True)
            time.sleep(5 * attempt)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--verify-only", action="store_true", help="only check, do not download")
    args = parser.parse_args()

    with open(SOURCES, newline="") as fh:
        sources = list(csv.DictReader(fh, delimiter="\t"))

    failed = []
    for source in sources:
        dest = ROOT / source["dest"]
        if dest.exists() and sha256(dest) == source["sha256"]:
            print(f"ok       {source['dest']}")
            continue
        if args.verify_only:
            print(f"MISSING or different: {source['dest']}")
            failed.append(source["dest"])
            continue
        if dest.exists():  # wrong content: raw files are read-only, so unlock before replacing
            dest.chmod(stat.S_IWUSR | stat.S_IRUSR)
        dest.parent.mkdir(parents=True, exist_ok=True)
        download(source["url"], dest)
        got = sha256(dest)
        if got != source["sha256"]:
            failed.append(source["dest"])
            print(f"CHECKSUM MISMATCH {source['dest']}\n  expected {source['sha256']}\n  got      {got}")
            continue
        dest.chmod(stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)  # read-only
        print(f"ok       {source['dest']}")

    if failed:
        sys.exit(f"{len(failed)} file(s) failed: {', '.join(failed)}")
    print("all raw files present and verified")


if __name__ == "__main__":
    main()
