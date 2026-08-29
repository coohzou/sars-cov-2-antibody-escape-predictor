"""
Download SARS-CoV-2 variant reference genomes from NCBI.
"""

import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from utils.paths import MANIFEST_PATH, PREDICTION_DATA_DIR

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def load_manifest():
    with open(MANIFEST_PATH, "r", encoding="utf-8") as f:
        return json.load(f)


def fetch_fasta(accession):
    params = urllib.parse.urlencode({
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
    })
    url = f"{EFETCH_URL}?{params}"
    with urllib.request.urlopen(url, timeout=60) as resp:
        return resp.read().decode("utf-8")


def normalize_fasta(text, accession, lineage):
    lines = [line.strip() for line in text.strip().splitlines() if line.strip()]
    if not lines:
        raise ValueError("Empty FASTA response")

    if lines[0].startswith(">"):
        header = lines[0]
    else:
        header = f">{accession} SARS-CoV-2 {lineage} complete genome"

    sequence = "".join(line for line in lines[1:] if not line.startswith(">")).upper()
    if len(sequence) < 29000:
        raise ValueError(f"Sequence too short ({len(sequence)} bp)")

    wrapped = [header]
    for i in range(0, len(sequence), 70):
        wrapped.append(sequence[i:i + 70])
    return "\n".join(wrapped) + "\n"


def download_all(force=False):
    manifest = load_manifest()
    results = []

    for item in manifest["variants"]:
        out_path = os.path.join(PREDICTION_DATA_DIR, item["file"])

        if os.path.exists(out_path) and not force:
            results.append({"name": item["name"], "status": "skipped", "path": out_path})
            continue

        print(f"Downloading {item['name']} ({item['accession']})...")
        try:
            raw = fetch_fasta(item["accession"])
            fasta = normalize_fasta(raw, item["accession"], item["lineage"])
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            with open(out_path, "w", encoding="utf-8") as f:
                f.write(fasta)
            results.append({
                "name": item["name"],
                "status": "ok",
                "path": out_path,
                "length": len("".join(fasta.splitlines()[1:])),
            })
            print(f"  Saved {out_path}")
        except (urllib.error.URLError, ValueError) as exc:
            results.append({"name": item["name"], "status": "failed", "error": str(exc)})
            print(f"  Failed: {exc}")

        time.sleep(0.8)

    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Download variant genomes from NCBI")
    parser.add_argument("--force", action="store_true", help="Re-download even if file exists")
    args = parser.parse_args()

    summary = download_all(force=args.force)
    ok = sum(1 for r in summary if r["status"] == "ok")
    skipped = sum(1 for r in summary if r["status"] == "skipped")
    failed = [r for r in summary if r["status"] == "failed"]
    print(f"\nDone: {ok} downloaded, {skipped} skipped, {len(failed)} failed")
    if failed:
        for item in failed:
            print(f"  - {item['name']}: {item['error']}")
