import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

RAIL = 32000
BAD_FRACTION = 0.10


def read_vhdr(path):
    txt = path.read_text(encoding="utf-8")
    assert "BinaryFormat=INT_16" in txt and "DataOrientation=MULTIPLEXED" in txt, path
    return re.findall(r"^Ch\d+=([^,]*),", txt, flags=re.M)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path, help="nm000433 dataset root")
    root = parser.parse_args().root
    report = json.loads((root / "code" / "conversion_report.json").read_text())
    saturated = {r["stem"]: r["saturated_samples"] for r in report["recordings"]}
    print("| recording | samples | channels | channels > 10%: % of samples at abs count >= 32000 "
          "| max % of the other channels | counts at 32767/-32768 = conversion_report saturated_samples |")
    print("|---|---:|---:|---|---|---|")
    summary = defaultdict(lambda: [0, 0, 0.0, 0.0])
    n_files = n_bad = 0
    for vhdr in sorted(root.glob("sub-*/ses-*/ieeg/*_ieeg.vhdr")):
        names = read_vhdr(vhdr)
        x = np.fromfile(vhdr.with_suffix(".eeg"), dtype="<i2").reshape(-1, len(names)).astype(np.int32)
        frac = (np.abs(x) >= RAIL).mean(axis=0)
        bad = frac > BAD_FRACTION
        listed = ", ".join(f"{n} {100 * f:.1f}" for n, f in zip(names, frac) if f > BAD_FRACTION)
        other = int(np.argmax(np.where(bad, -1.0, frac)))
        stem = vhdr.name.removesuffix("_ieeg.vhdr")
        if stem in saturated:
            sat = ((x == 32767) | (x == -32768)).sum(axis=0)
            same = "yes" if all(int(s) == saturated[stem][n] for n, s in zip(names, sat)) else "NO"
        else:
            same = "n/a (analog inputs)"
        worst = f"{100 * frac[other]:.2f} ({names[other]})" if frac[other] > 0 else "0"
        print(f"| {stem} | {len(x)} | {len(names)} | {bad.sum()}: {listed or '-'} | {worst} | {same} |")
        acq = re.search(r"_acq-([a-z]+)_", stem).group(1)
        s = summary[acq]
        s[0] += 1
        s[1] += int(bad.sum())
        lfp = [not b and n.startswith("LFP") for n, b in zip(names, bad)]
        rest = [not b and not n.startswith("LFP") for n, b in zip(names, bad)]
        s[2] = max(s[2], 100 * frac[lfp].max(initial=0.0))
        s[3] = max(s[3], 100 * frac[rest].max(initial=0.0))
        n_files += 1
        n_bad += int(bad.sum())
    print("\n| acq | recordings | channels > 10% | max % of the other LFP channels | max % of the other channels |")
    print("|---|---:|---:|---:|---:|")
    for acq, (n, b, lfp, rest) in sorted(summary.items()):
        print(f"| {acq} | {n} | {b} | {lfp:.2f} | {rest:.2f} |")
    print(f"\n{n_files} recordings, {n_bad} channels with more than 10% of samples at abs count >= {RAIL}.")


if __name__ == "__main__":
    main()
