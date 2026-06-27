#!/usr/bin/env python3
"""Validate the GC32 DSI reproducibility package file inventory."""
from pathlib import Path
import hashlib, csv, sys

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / 'metadata' / 'sha256_manifest.csv'


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    if not MANIFEST.is_file():
        print(f'Missing manifest: {MANIFEST}')
        return 2
    failures = []
    with MANIFEST.open(newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rel = row['path']
            p = ROOT / rel
            if not p.is_file():
                failures.append((rel, 'missing'))
                continue
            digest = sha256(p)
            if digest != row['sha256']:
                failures.append((rel, 'checksum_mismatch'))
    if failures:
        for rel, reason in failures:
            print(f'{reason}: {rel}')
        return 1
    print('All package files passed checksum validation.')
    return 0

if __name__ == '__main__':
    sys.exit(main())
