"""d3-radboud-dareplane (Radboud Data Repository doi:10.34973/d214-m342, CC BY 4.0) -> BIDS-iEEG (raw).

Patient data only: sub-P001, sessions day2/day3/day4 after DBS-lead implantation (Dold et al., Dareplane proof of concept,
arXiv:2408.01242). Benchtop tests (AO/Arduino/CT/Picoscope XDF) and the c-VEP scalp-EEG speller data are not included.

Measured (probe Jobs ieeg-b3w11-rdprobe-1008004324, ieeg-b3w11-rdprobe2-1008025023; SHA-256 of all 685 downloaded files
matched the release MANIFEST.txt):
- 41 Neuro Omega .mat exports (MAT v5) in data/sub-P001_ses-day*/ao/. Each holds the 24 channels of the README CHANNEL_MAP
  (CECOG_HF_1_01..16 = LFP_1..16 = two STN DBS leads; CECOG_HF_2_01..04 = ECOG_1..4; CECOG_HF_2_05..08 = EOG_1..4) as int16
  at 22 kHz with BitResolution 38.14697265625 and Gain 55, optional analog inputs (CANALOG_IN_1..8, CADD_ANALOG_IN_1..8,
  int16, 2.75 kHz), hardware port markers CPORT__1 (44 kHz timestamps + value), stimulation markers CStimMarker_* (44 kHz),
  spike-segment channels (CSEG_*, not continuous) and SF_* system-setting tables.
- Conversion: the 24 mapped channels -> BrainVision INT_16, the int16 counts unchanged, resolution BitResolution/Gain uV per
  count (Neuro Omega convention value_uV = count * BitResolution / Gain). Analog inputs at their native 2.75 kHz go to a
  separate acq-<cond>analog file when they are not constant. Round-trip: the written .eeg is read back and compared to
  the source int16 arrays (exact equality required).
- Events: CPORT__1 hardware port values (non-zero; zero resets the port, see analysis_scripts/load_xdf.py get_AO_data) and
  CStimMarker_* pulse markers, onsets = timestamp/44 kHz - TimeBegin of the 22 kHz channels (the alignment used by the
  authors' load_xdf.py). The port values are reported as recorded: they do not equal the CopyDraw paradigm codes of the
  README (10-21, 110-121, 50) on every day; the complete CopyDraw marker stream is in the XDF files (LSL clock) in sourcedata.
- Privacy: the .mat text header holds the export date ('Created on: ... Nov 27 ... 2024'): day set to 01. XDF files: every
  date string (YYYY-MM-DD / YYYYMMDD) in header XML and in string sample chunks is set to day 01, same length; other
  files are scanned and patched the same way. Every patch is logged (file, offset, chunk type; not the old value).
Usage: python b3w11_convert_radboud.py <work_dir> <bids_root>
"""
import glob
import hashlib
import json
import os
import re
import shutil
import struct
import sys
from collections import Counter, defaultdict

import numpy as np
import scipy.io

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b3w3_common import scrub_mat_header, sha256_file, wjson, write_vhdr, wtsv  # noqa: E402

W, OUT = sys.argv[1], sys.argv[2]
DL = os.path.join(W, "sourcedata_dl")
SRC = os.path.join(DL, "data")
SUB = "sub-P001"
SD = os.path.join(OUT, "sourcedata", "radboud-d214-m342")
rep = {"recordings": [], "privacy_patches": [], "notes": []}

CHANNEL_MAP = [(f"CECOG_HF_1___{i:02d}___Array_{1 if i <= 8 else 2}___{i:02d}", f"LFP_{i}") for i in range(1, 17)] + \
              [(f"CECOG_HF_2___{i:02d}___Array_3___{i:02d}", f"ECOG_{i}") for i in range(1, 5)] + \
              [(f"CECOG_HF_2___{i:02d}___Array_3___{i:02d}", f"EOG_{i - 4}") for i in range(5, 9)]
ANALOG = [f"CANALOG_IN_{i}" for i in range(1, 9)] + [f"CADD_ANALOG_IN_{i}" for i in range(1, 9)]


def ctype(n):
    return "DBS" if n.startswith("LFP") else ("ECOG" if n.startswith("ECOG") else "EOG")


def cgroup(n):
    k = int(n.split("_")[1])
    if n.startswith("LFP"):
        return "DBS_left_STN" if k <= 8 else "DBS_right_STN"
    return "ECOG_strip_left_M1" if n.startswith("ECOG") else "EOG"


def cdesc(n):
    k = int(n.split("_")[1])
    if n.startswith("LFP"):
        side, j = ("left", k) if k <= 8 else ("right", k - 8)
        pos = " (tip, most ventral)" if j == 1 else (" (most dorsal)" if j == 8 else "")
        return f"Boston Scientific Vercise directional DBS lead in the {side} STN, contact {j} of 8{pos} (README)"
    if n.startswith("ECOG"):
        pos = {1: " (most caudal)", 4: " (most rostral, closest to the left borehole)"}.get(k, "")
        return f"4-contact Ad-Tech strip, epidural over left primary motor cortex, contact {k}{pos} (README; paper)"
    return {1: "EOG above the right eyebrow", 2: "EOG below the right eye", 3: "EOG horizontal, left of the left eye",
            4: "EOG horizontal, right of the right eye"}[k] + " (README)"


def ref_of(n):
    if n.startswith("LFP"):
        return "left mastoid"
    if n.startswith("ECOG"):
        return "left earlobe"
    return "n/a"


def scal(x):
    x = np.ravel(x)
    return x[0].item() if x.size else None


# ---------------------------------------------------------------- signals
COND = {"copydraw_off": ("stimoff", False), "copydraw_on": ("stimon", True), "clcopydraw_on": ("closedloop", True)}
TASKDESC = ("CopyDraw: the patient copies traces shown on a screen with a stylus on a Wacom tablet (Castano-Candamil et al. 2019); "
            "self-paced readiness confirmation (marker 50), then copy-drawing of a trace (markers 10-21 start, 110-121 end, ~9 s). "
            "Blocks were recorded with DBS OFF, with continuous open-loop DBS ON, or (day 4) with closed-loop adaptive DBS.")
STIMDESC = ("Paper (Dold et al.): 130 Hz, 60 us symmetric initial-negative pulses, monolateral (left or right STN) bipolar stimulation "
            "between contacts of the DBS lead chosen by a clinician; day 2: left-hemisphere amplitude reduced to 4 mA in the first "
            "stimulation block (feeling of heat); day 4: 6 mA. Closed-loop blocks (day 4): adaptive DBS switched ON/OFF by a "
            "Bollinger-band controller on an ECoG-decoded CopyDraw score.")
scans = defaultdict(list)
allmats = []
for ses_dir in sorted(glob.glob(os.path.join(SRC, f"{SUB}_ses-day*"))):
    ses = "ses-" + ses_dir.rsplit("_ses-", 1)[1]
    items = []
    for f in sorted(glob.glob(os.path.join(ses_dir, "ao", "*.mat"))):
        m = re.match(r"block(\d+)_(copydraw_off|copydraw_on|clcopydraw_on)(\d{4})?\.mat$", os.path.basename(f))
        assert m, f
        hdr = scipy.io.loadmat(f, variable_names=[CHANNEL_MAP[0][0] + "_TimeBegin"])
        items.append((float(scal(hdr[CHANNEL_MAP[0][0] + "_TimeBegin"])), f, int(m.group(1)), m.group(2), m.group(3) or "-"))
    items.sort()
    if os.environ.get("SMOKE"):
        items = items[:2]
    for run, (tb, f, block, cond, part) in enumerate(items, 1):
        acq, stim = COND[cond]
        m = scipy.io.loadmat(f, squeeze_me=False)
        src_rel = os.path.relpath(f, DL)
        D = os.path.join(OUT, SUB, ses, "ieeg")
        os.makedirs(D, exist_ok=True)
        stem = f"{SUB}_{ses}_task-copydraw_acq-{acq}_run-{run:02d}"
        k0 = CHANNEL_MAP[0][0]
        fs = float(scal(m[k0 + "_KHz"])) * 1000.0
        br, gain = float(scal(m[k0 + "_BitResolution"])), float(scal(m[k0 + "_Gain"]))
        tbegin = float(scal(m[k0 + "_TimeBegin"]))
        arrs = []
        for src, name in CHANNEL_MAP:
            a = np.ravel(m[src])
            assert a.dtype == np.int16, (f, src, a.dtype)
            assert float(scal(m[src + "_KHz"])) * 1000.0 == fs and float(scal(m[src + "_BitResolution"])) == br and float(scal(m[src + "_Gain"])) == gain, (f, src)
            assert abs(float(scal(m[src + "_TimeBegin"])) - tbegin) < 1e-9, (f, src)
            arrs.append(a)
        n = min(len(a) for a in arrs)
        assert all(len(a) == n for a in arrs), (f, [len(a) for a in arrs])
        X = np.stack(arrs)  # (24, n) int16
        X.T.astype("<i2").tofile(os.path.join(D, stem + "_ieeg.eeg"))
        res = br / gain
        names = [nm for _, nm in CHANNEL_MAP]
        write_vhdr(os.path.join(D, stem + "_ieeg"), len(names), fs, names, ["µV"] * len(names), fmt="INT_16", resolutions=[res] * len(names),
                   comment=f"b3w11_convert_radboud.py: Neuro Omega int16 counts unchanged; resolution = BitResolution {br!r} / Gain {gain!r} uV")
        back = np.fromfile(os.path.join(D, stem + "_ieeg.eeg"), dtype="<i2").reshape(-1, len(names)).T
        rt = bool(back.shape == X.shape and np.array_equal(back, X))
        del back
        # events
        ev = []
        if "CPORT__1" in m and m["CPORT__1"].size:
            P = np.asarray(m["CPORT__1"])
            pk = float(scal(m["CPORT__1_KHz"])) * 1000.0
            for t, v in zip(P[0], P[1]):
                if v != 0:
                    ev.append((t / pk - tbegin, "port_marker", int(v)))
        nstim = Counter()
        for sk in sorted(k for k in m if re.fullmatch(r"CStimMarker_\d+", k)):
            S = np.asarray(m[sk])
            if not S.size:
                continue
            sk_k = float(scal(m[sk + "_KHz"])) * 1000.0
            for t, v in zip(S[0], S[1]):
                ev.append((t / sk_k - tbegin, "stim_marker_" + sk.split("_")[1], int(v)))
                nstim[sk] += 1
        dur = n / fs
        inside = [e for e in ev if 0 <= e[0] < dur]
        outside = len(ev) - len(inside)
        inside.sort()
        if inside:
            wtsv(os.path.join(D, stem + "_events.tsv"), ["onset", "duration", "trial_type", "value", "sample"],
                 [[f"{o:.6f}", 0, tt, v, int(round(o * fs))] for o, tt, v in inside])
        # channels / sidecar
        wtsv(os.path.join(D, stem + "_channels.tsv"),
             ["name", "type", "units", "low_cutoff", "high_cutoff", "sampling_frequency", "reference", "group", "status", "status_description", "description"],
             [[nm, ctype(nm), "µV", "n/a", "n/a", fs, ref_of(nm), cgroup(nm), "good", "n/a",
               f"{cdesc(nm)}; Neuro Omega channel {src}; int16 counts x {br!r}/{gain!r} uV"] for src, nm in CHANNEL_MAP])
        sat = {nm: int(np.count_nonzero((X[i] == 32767) | (X[i] == -32768))) for i, nm in enumerate(names)}
        wjson(os.path.join(D, stem + "_ieeg.json"), {
            "TaskName": "copydraw", "TaskDescription": TASKDESC,
            "Instructions": "n/a",
            "SamplingFrequency": fs, "PowerLineFrequency": 50,
            "SoftwareFilters": "n/a", "HardwareFilters": "n/a",
            "Manufacturer": "Alpha Omega", "ManufacturersModelName": "Neuro Omega",
            "InstitutionName": "Medical Center - University of Freiburg, Department of Stereotactic and Functional Neurosurgery",
            "iEEGReference": "Paper: ECoG referenced to the left earlobe, DBS LFP referenced to the left mastoid; EOG reference not stated.",
            "iEEGGround": "right mastoid (paper)",
            "iEEGPlacementScheme": "Bilateral STN Boston Scientific Vercise directional DBS leads (8 contacts each, externalised) and one 4-contact Ad-Tech ECoG strip implanted epidurally over the left primary motor cortex through the left borehole; 4 EOG electrodes.",
            "iEEGElectrodeGroups": "DBS_left_STN: LFP_1-8; DBS_right_STN: LFP_9-16; ECOG_strip_left_M1: ECOG_1-4",
            "ElectrodeManufacturer": "Boston Scientific (DBS leads); Ad-Tech Medical Instrument Corporation (ECoG strip)",
            "DBSChannelCount": 16, "ECOGChannelCount": 4, "EOGChannelCount": 4, "SEEGChannelCount": 0,
            "RecordingType": "continuous", "RecordingDuration": round(dur, 6),
            "ElectricalStimulation": stim,
            "ElectricalStimulationParameters": (STIMDESC if stim else "Stimulation OFF block (file name 'copydraw_off').") +
                                               " Stimulation pulse markers of the Neuro Omega are in events.tsv (stim_marker_*).",
            "SubjectArtefactDescription": "Stimulation artefacts in ON/closed-loop blocks; LFP channels reach the int16 limits (saturated samples per channel are counted in code/conversion_report.json).",
        })
        scans[(SUB, ses)].append([f"ieeg/{stem}_ieeg.vhdr", "n/a", "sourcedata/radboud-d214-m342/" + src_rel,
                                  block, part, acq, f"{tbegin:.6f}"])
        rec = dict(src=src_rel, stem=stem, sfreq=fs, n_ch=len(names), n_samples=int(n), dur_s=dur, bitres=br, gain=gain, res_uV=res,
                   n_port_events=sum(1 for e in inside if e[1] == "port_marker"), n_stim_events=sum(1 for e in inside if e[1] != "port_marker"),
                   events_outside=outside, saturated_samples=sat, roundtrip_ok=rt)
        # analog inputs
        an = [a for a in ANALOG if a in m]
        if an:
            A = np.stack([np.ravel(m[a]) for a in an])
            nonconst = [a for i, a in enumerate(an) if A[i].min() != A[i].max()]
            rec["analog_nonconstant"] = nonconst
            if nonconst:
                afs = float(scal(m[an[0] + "_KHz"])) * 1000.0
                abr, ag = float(scal(m[an[0] + "_BitResolution"])), float(scal(m[an[0] + "_Gain"]))
                assert A.dtype == np.int16 and all(float(scal(m[a + "_KHz"])) * 1000.0 == afs for a in an)
                astem = f"{SUB}_{ses}_task-copydraw_acq-{acq}analog_run-{run:02d}"
                A.T.astype("<i2").tofile(os.path.join(D, astem + "_ieeg.eeg"))
                write_vhdr(os.path.join(D, astem + "_ieeg"), len(an), afs, an, [""] * len(an), fmt="INT_16", resolutions=[1.0] * len(an),
                           comment=f"b3w11_convert_radboud.py: Neuro Omega analog-input int16 counts unchanged (BitResolution {abr!r}, Gain {ag!r}; unit not documented)")
                backa = np.fromfile(os.path.join(D, astem + "_ieeg.eeg"), dtype="<i2").reshape(-1, len(an)).T
                rec["analog_roundtrip_ok"] = bool(np.array_equal(backa, A))
                rec["roundtrip_ok"] = rec["roundtrip_ok"] and rec["analog_roundtrip_ok"]
                wtsv(os.path.join(D, astem + "_channels.tsv"), ["name", "type", "units", "low_cutoff", "high_cutoff", "sampling_frequency", "status", "status_description", "description"],
                     [[a, "MISC", "n/a", "n/a", "n/a", afs, "good" if a in nonconst else "bad", "n/a" if a in nonconst else "constant over the recording",
                       f"Neuro Omega analog input {a} (content not documented in the record); raw int16 counts; BitResolution {abr!r}, Gain {ag!r}"] for a in an])
                wjson(os.path.join(D, astem + "_ieeg.json"), {
                    "TaskName": "copydraw", "TaskDescription": TASKDESC, "SamplingFrequency": afs, "PowerLineFrequency": 50,
                    "SoftwareFilters": "n/a", "HardwareFilters": "n/a", "Manufacturer": "Alpha Omega", "ManufacturersModelName": "Neuro Omega",
                    "iEEGReference": "n/a (auxiliary analog inputs)", "MiscChannelCount": len(an), "ECOGChannelCount": 0, "SEEGChannelCount": 0,
                    "RecordingType": "continuous", "RecordingDuration": round(A.shape[1] / afs, 6), "ElectricalStimulation": stim,
                    "AnalogInputTimeBeginSeconds": float(scal(m[an[0] + "_TimeBegin"])),
                    "Description": f"Auxiliary analog inputs recorded by the Neuro Omega at their native {afs} Hz during the same recording as {stem}_ieeg (Neuro Omega clock start {float(scal(m[an[0] + '_TimeBegin'])):.6f} s vs {tbegin:.6f} s for the 22 kHz channels)."})
                scans[(SUB, ses)].append([f"ieeg/{astem}_ieeg.vhdr", "n/a", "sourcedata/radboud-d214-m342/" + src_rel, block, part, acq + "analog",
                                          f"{float(scal(m[an[0] + '_TimeBegin'])):.6f}"])
                del A, backa
        rep["recordings"].append(rec)
        print(stem, X.shape, fs, "ev", len(inside), "out", outside, "rt", rec["roundtrip_ok"], "analog", rec.get("analog_nonconstant"), flush=True)
        del X, m, arrs

for (sub, ses), r in scans.items():
    wtsv(os.path.join(OUT, sub, ses, f"{sub}_{ses}_scans.tsv"), ["filename", "acq_time", "source_file", "block", "file_part", "condition", "neuro_omega_time_begin_s"], r)
wjson(os.path.join(OUT, "scans.json"), {
    "source_file": {"Description": "Original Neuro Omega .mat export in sourcedata/ (export date in the MAT header set to day 01)"},
    "block": {"Description": "CopyDraw block number from the source file name"},
    "file_part": {"Description": "Neuro Omega file index (0001, 0002, ...) from the source file name; '-' when the name has none"},
    "condition": {"Description": "stimoff / stimon (open-loop DBS) / closedloop (adaptive DBS), from the file name; suffix 'analog' = auxiliary analog-input file"},
    "neuro_omega_time_begin_s": {"Description": "TimeBegin of the channels in this file on the Neuro Omega clock (s); runs are numbered in this order within a session", "Units": "s"},
})

# ---------------------------------------------------------------- electrodes
for ses in sorted({s for _, s in scans}):
    D = os.path.join(OUT, SUB, ses, "ieeg")
    pref = f"{SUB}_{ses}"
    rows = []
    for _, nm in CHANNEL_MAP:
        if nm.startswith("EOG"):
            continue
        k = int(nm.split("_")[1])
        hemi = "L" if (nm.startswith("ECOG") or k <= 8) else "R"
        region = "subthalamic nucleus" if nm.startswith("LFP") else "primary motor cortex (epidural)"
        manu = "Boston Scientific (Vercise directional lead)" if nm.startswith("LFP") else "Ad-Tech Medical Instrument Corporation"
        rows.append([nm, "n/a", "n/a", "n/a", "n/a", cgroup(nm), hemi, region, manu, cdesc(nm)])
    wtsv(os.path.join(D, f"{pref}_space-Other_electrodes.tsv"), ["name", "x", "y", "z", "size", "group", "hemisphere", "region_stated", "manufacturer", "description"], rows)
    wjson(os.path.join(D, f"{pref}_space-Other_electrodes.json"), {
        "hemisphere": {"Description": "Hemisphere stated in the README/paper", "Levels": {"L": "left", "R": "right"}},
        "region_stated": {"Description": "Implantation target stated in the README and paper (no localisation released)"},
        "manufacturer": {"Description": "Electrode manufacturer stated in the paper"},
        "description": {"Description": "Contact position as stated in the release README"}})
    wjson(os.path.join(D, f"{pref}_space-Other_coordsystem.json"), {
        "iEEGCoordinateSystem": "Other", "iEEGCoordinateUnits": "n/a",
        "iEEGCoordinateSystemDescription": "No electrode coordinates or imaging are released; contacts are listed with x, y, z = n/a and the stated target."})

wtsv(os.path.join(OUT, "participants.tsv"), ["participant_id", "age", "sex", "handedness", "diagnosis", "paper_id"], [[SUB, 67, "F", "R", "Parkinson's disease", "PD4"]])
wjson(os.path.join(OUT, "participants.json"), {
    "age": {"Description": "Age at the recordings, as published (Dold et al., section 2.6)", "Units": "year"},
    "sex": {"Description": "As published", "Levels": {"F": "female"}},
    "handedness": {"Description": "As published", "Levels": {"R": "right"}},
    "diagnosis": {"Description": "As published"},
    "paper_id": {"Description": "Patient label in Table 5 of the paper"}})
wtsv(os.path.join(OUT, SUB, f"{SUB}_sessions.tsv"), ["session_id", "days_after_implantation", "brainamp_physiology_in_sourcedata", "closed_loop_blocks"],
     [["ses-day2", 2, "yes", "no"], ["ses-day3", 3, "yes", "no"], ["ses-day4", 4, "no", "yes"]])
wjson(os.path.join(OUT, "sessions.json"), {
    "days_after_implantation": {"Description": "Days after the first surgery (DBS-lead and ECoG-strip implantation with externalisation), from the session folder name and README", "Units": "day"},
    "brainamp_physiology_in_sourcedata": {"Description": "BrainAmp ExG heart rate, respiration and galvanic skin response (BrainVision, 5 kHz) present in sourcedata"},
    "closed_loop_blocks": {"Description": "Session contains closed-loop adaptive DBS blocks"}})
wjson(os.path.join(OUT, "task-copydraw_events.json"), {
    "trial_type": {"Description": "port_marker: hardware port CPORT__1 of the Neuro Omega (non-zero values; zero resets the port). stim_marker_<k>: Neuro Omega stimulation marker channel CStimMarker_<k> (README: arrival times of the stimulation pulses)."},
    "value": {"Description": "port_marker: the port value as recorded. The CopyDraw paradigm codes (README: 10-21 trace start, 110-121 trace end, 50 readiness) are only partly reproduced by the port values on some days; the complete marker stream (CopyDrawParadigmMarkerStream) is in the XDF files in sourcedata on the LSL clock. stim_marker: the identifier stored with the marker by the Neuro Omega."},
    "sample": {"Description": "onset x 22000 (0-based sample of the 22 kHz file)"},
    "onset": {"Description": "marker timestamp (44 kHz Neuro Omega clock) minus TimeBegin of the 22 kHz channels, as in the authors' load_xdf.py", "Units": "s"}})

# ---------------------------------------------------------------- sourcedata (privacy-scrubbed copies)
DATE_SEP = re.compile(rb"(?<![0-9])(20[0-3][0-9]|19[0-9][0-9])([-/.])(0[1-9]|1[0-2])\2(0[1-9]|[12][0-9]|3[01])(?![0-9])")
DATE_CMP = re.compile(rb"(?<![0-9.eE+-])(20[12][0-9])(0[1-9]|1[0-2])(0[1-9]|[12][0-9]|3[01])(?=[0-9]{0,9}(?![0-9]))")


def patch_dates(path, rel, regions=None):
    """In-place same-length day->01 for date strings. regions: list of (start, end, kind) byte ranges to scan (None = whole file)."""
    n = 0
    with open(path, "r+b") as fh:
        data = fh.read() if regions is None else None
        regs = [(0, len(data), "file")] if regions is None else regions
        for s, e, kind in regs:
            if data is None:
                fh.seek(s); buf = fh.read(e - s)
            else:
                buf = data[s:e]
            for rx, dpos in ((DATE_SEP, 4), (DATE_CMP, 3)):
                for mm in rx.finditer(buf):
                    day = mm.group(dpos)
                    if day == b"01":
                        continue
                    off = s + mm.start(dpos)
                    fh.seek(off); fh.write(b"01"); n += 1
                    rep["privacy_patches"].append(dict(file=rel, offset=off, kind=kind, pattern=rx.pattern[:12].decode("latin1", "replace")))
    return n


def xdf_regions(path):
    """XML header/footer chunks + string-format sample chunks."""
    regs, fmts = [], {}
    with open(path, "rb") as fh:
        assert fh.read(4) == b"XDF:"
        while True:
            pos = fh.tell()
            nb = fh.read(1)
            if not nb:
                break
            L = int.from_bytes(fh.read(nb[0]), "little")
            tag = struct.unpack("<H", fh.read(2))[0]
            start = fh.tell()
            if tag in (1, 2, 6):
                body = fh.read(L - 2)
                regs.append((start, start + L - 2, {1: "xdf_fileheader", 2: "xdf_streamheader", 6: "xdf_streamfooter"}[tag]))
                if tag == 2:
                    sid = struct.unpack("<I", body[:4])[0]
                    fm = re.search(rb"<channel_format>(.*?)</channel_format>", body)
                    fmts[sid] = fm.group(1) if fm else b""
            elif tag == 3:
                sid = struct.unpack("<I", fh.read(4))[0]
                if fmts.get(sid) == b"string":
                    regs.append((start, start + L - 2, "xdf_string_samples"))
                fh.seek(start + L - 2)
            else:
                fh.seek(start + L - 2)
    return regs


os.makedirs(SD, exist_ok=True)
for fn in ("README.md", "ABOUT.txt", "LICENSE.txt", "MANIFEST.txt"):
    shutil.copy(os.path.join(DL, fn), os.path.join(SD, fn))
os.makedirs(os.path.join(SD, "analysis_scripts"), exist_ok=True)
shutil.copy(os.path.join(DL, "analysis_scripts", "load_xdf.py"), os.path.join(SD, "analysis_scripts", "load_xdf.py"))
changed = {}
for f in sorted(glob.glob(os.path.join(SRC, f"{SUB}_ses-day*", "**", "*"), recursive=True)):
    if not os.path.isfile(f):
        continue
    if os.environ.get("SMOKE") and not (f.endswith((".vmrk", ".vhdr")) or "block1_" in f or "block_01" in f):
        continue
    rel = os.path.relpath(f, DL)
    dst = os.path.join(SD, rel)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(f, dst)
    k = 0
    if f.endswith(".mat"):
        b = open(dst, "rb").read(128)
        nb_, ch = scrub_mat_header(b)
        if ch:
            with open(dst, "r+b") as fh:
                fh.write(nb_[:116])
            k += 1
            rep["privacy_patches"].append(dict(file=rel, offset=0, kind="mat_text_header", pattern="Created on"))
    elif f.endswith(".xdf"):
        k += patch_dates(dst, rel, xdf_regions(dst))
    elif f.endswith((".vhdr", ".vmrk", ".txt", ".json", ".csv")):
        k += patch_dates(dst, rel)
    # .yaml (behaviour): only floats and LSL-clock timestamps (probe2: no date strings); checked by the residual scan below
    if k:
        changed[rel] = k
rep["sourcedata_changed_files"] = changed
# residual scan of every text-like/xml region after patching (must be empty except day 01)
resid = []
for f in sorted(glob.glob(os.path.join(SD, "**", "*"), recursive=True)):
    if os.path.isfile(f) and f.endswith((".yaml", ".vhdr", ".vmrk", ".json", ".csv")):
        t = open(f, "rb").read()
        for rx, dpos in ((DATE_SEP, 4),) + (() if f.endswith(".yaml") else ((DATE_CMP, 3),)):
            resid += [os.path.relpath(f, SD) for mm in rx.finditer(t) if mm.group(dpos) != b"01"]
rep["residual_dates_text"] = resid[:50]
wtsv(os.path.join(SD, "PRIVACY_CHANGES.tsv"), ["file", "n_patches", "sha256_after_patch"],
     [[r, n, sha256_file(os.path.join(SD, r))] for r, n in sorted(changed.items())])

os.makedirs(os.path.join(OUT, "code"), exist_ok=True)
for fn in (os.path.abspath(__file__), os.path.join(os.path.dirname(os.path.abspath(__file__)), "b3w3_common.py")):
    shutil.copy(fn, os.path.join(OUT, "code", os.path.basename(fn)))
json.dump(rep, open(os.path.join(OUT, "code", "conversion_report.json"), "w"), indent=1)

# ---------------------------------------------------------------- top-level text
AUTH = ["Matthias Dold", "Joana Pereira", "Bastian Sajonz", "Volker A. Coenen", "Jordy Thielen", "Marcus L.F. Janssen", "Michael Tangermann"]
nrun = sum(1 for r in rep["recordings"])
hrs = sum(r["dur_s"] for r in rep["recordings"]) / 3600
nan = sum(1 for r in rep["recordings"] if r.get("analog_nonconstant"))
wjson(os.path.join(OUT, "dataset_description.json"), {
    "Name": "Dareplane proof of concept: STN LFP, ECoG and EOG during the CopyDraw task with open- and closed-loop DBS (1 Parkinson's disease patient)",
    "BIDSVersion": "1.10.0", "DatasetType": "raw", "License": "CC-BY-4.0", "Authors": AUTH,
    "HowToAcknowledge": "Cite the dataset (Dold, M., Pereira, J., Sajonz, B., Coenen, V.A., Thielen, J., Janssen, M.L.F., Tangermann, M. (2024). LFP and ECoG data during CopyDraw task with deep brain stimulation - Dareplane data for proof of concept paper. Version 1. Radboud University. https://doi.org/10.34973/d214-m342) and the Dareplane paper (arXiv:2408.01242). The release asks to acknowledge: 'Data were provided (in part) by the Radboud University, Nijmegen, The Netherlands'.",
    "Funding": ["Dareplane collaboration project, co-funded by the PPP Allowance of Health Holland, Top Sector Life Sciences & Health, and the Dutch Brain Foundation (MD)",
                "Deutsche Forschungsgemeinschaft, Walter Benjamin Program, project 510112977 (JP)",
                "PD-Interaktiv I study: Bundesministerium fuer Bildung und Forschung, grant 16SV8011"],
    "EthicsApprovals": ["Ethics committee of the Medical Center - University of Freiburg, EK application 21-1545; study PD-Interaktiv I, DRKS000287039 (as stated in the paper)"],
    "ReferencesAndLinks": ["https://doi.org/10.34973/d214-m342", "https://arxiv.org/abs/2408.01242", "https://github.com/bsdlab/Dareplane"],
    "DatasetDOI": "n/a",
    "SourceDatasets": [{"DOI": "doi:10.34973/d214-m342", "URL": "https://doi.org/10.34973/d214-m342", "Version": "1"}],
    "Keywords": ["deep brain stimulation", "adaptive DBS", "closed-loop", "subthalamic nucleus", "local field potential", "ECoG", "Parkinson's disease", "CopyDraw", "Dareplane"],
    "GeneratedBy": [{"Name": "b3w11_convert_radboud.py", "Description": "Neuro Omega .mat exports to BIDS-iEEG BrainVision (int16 counts unchanged); events from the Neuro Omega port and stimulation markers; privacy-scrubbed originals in sourcedata. Code in code/."}],
})
open(os.path.join(OUT, "README"), "w").write(f"""# Dareplane proof of concept: STN LFP, ECoG and EOG during CopyDraw with open- and closed-loop DBS

Intracranial recordings from one patient with Parkinson's disease (67 years, female, right-handed; "PD4" in the paper) with
externalised bilateral subthalamic DBS leads and an epidural ECoG strip over the left primary motor cortex, recorded on
days 2, 3 and 4 after implantation at the Medical Center - University of Freiburg while she performed the CopyDraw task
with DBS OFF, with continuous open-loop DBS ON and (day 4) with closed-loop adaptive DBS run on the Dareplane platform.

Source: Dold, M., Pereira, J., Sajonz, B., Coenen, V.A., Thielen, J., Janssen, M.L.F., Tangermann, M. (2024). "LFP and ECoG
data during CopyDraw task with deep brain stimulation - Dareplane data for proof of concept paper". Version 1. Radboud
University. https://doi.org/10.34973/d214-m342 (CC BY 4.0). Paper: Dold et al., "Dareplane: A modular open-source software
platform for BCI research with application in closed-loop deep brain stimulation", arXiv:2408.01242 (J. Neural Eng.).
Please acknowledge: "Data were provided (in part) by the Radboud University, Nijmegen, The Netherlands".

## Contents
- `sub-P001` (release label), sessions `ses-day2`, `ses-day3`, `ses-day4` (days after implantation).
- {nrun} recordings, {hrs:.2f} h at 22 kHz. One recording = one Neuro Omega .mat file of the release; runs are numbered in
  recording order (Neuro Omega clock) within a session. `acq-stimoff` / `acq-stimon` / `acq-closedloop` follow the file
  names (`copydraw_off`, `copydraw_on`, `clcopydraw_on`). `scans.tsv` gives the source file, block number and file part
  for every run. Several blocks were saved in more than one file (e.g. block1 parts 0001-0004 on day 2), some only a few
  seconds long; they are kept as separate runs, unchanged.
- 24 channels per recording, names from the release README CHANNEL_MAP:
  - `LFP_1`-`LFP_8`: Boston Scientific Vercise directional lead in the left STN (`LFP_1` tip, most ventral; `LFP_8` most dorsal);
  - `LFP_9`-`LFP_16`: same for the right STN (`LFP_9` tip);
  - `ECOG_1`-`ECOG_4`: 4-contact Ad-Tech strip, epidural over left M1 (`ECOG_1` most caudal, `ECOG_4` most rostral);
  - `EOG_1` above the right eyebrow, `EOG_2` below the right eye, `EOG_3` left of the left eye, `EOG_4` right of the right eye.
  Type DBS / ECOG / EOG. The README does not say which of contacts 2-7 are the directional segments.
- Values: the Neuro Omega int16 counts are stored unchanged (BrainVision INT_16); the header resolution is
  BitResolution/Gain = 38.14697265625/55 = {38.14697265625 / 55:.10f} uV per count (Neuro Omega convention). Note that the
  authors' load_xdf.py treats the raw counts directly as microvolts (no gain); multiply by 55/38.14697265625 to reproduce that.
  LFP channels saturate (int16 limits) during stimulation; counts per channel are in code/conversion_report.json.
- Reference and ground (paper): ECoG referenced to the left earlobe, DBS recordings to the left mastoid, ground on the
  right mastoid. Hardware filter settings are not stated (n/a).
- Auxiliary Neuro Omega analog inputs (CANALOG_IN_1-8, CADD_ANALOG_IN_1-8, 2.75 kHz) are written as separate
  `acq-<condition>analog` files (MISC, raw counts, unit not documented) only where they are not constant ({nan} recordings).
- No coordinates are released. electrodes.tsv lists the 16 DBS and 4 ECoG contacts with x, y, z = n/a and the stated
  target (STN / left M1, hemisphere).

## Events
`events.tsv` holds two kinds of markers stored in the .mat files, aligned as in the authors' load_xdf.py
(onset = marker time on the 44 kHz Neuro Omega clock - TimeBegin of the 22 kHz channels):
- `port_marker`: non-zero values of the Neuro Omega hardware port CPORT__1 (zeros reset the port). The values are given as
  recorded. They do not always equal the CopyDraw codes (README: 10-21 start of copy-drawing trace k, 110-121 end, 50
  self-paced readiness): on day 3 the port only shows value 1 and on days 2 and 4 only part of the code range appears.
- `stim_marker_<k>`: the Neuro Omega stimulation marker channel CStimMarker_<k> (README: arrival times of the stimulation
  pulses); `value` is the identifier stored with the marker.
The full CopyDraw marker stream (CopyDrawParadigmMarkerStream), the LSL copy of the signals (AODataStream), the closed-loop
control streams and the mouse/keyboard streams are in the XDF files in sourcedata, on the LSL clock. The release does not
document an alignment between the LSL clock and the .mat sample clock other than the port markers, so the XDF markers are
not copied into events.tsv.

## sourcedata/radboud-d214-m342
The release files of the patient (README.md, ABOUT.txt, LICENSE.txt, MANIFEST.txt, analysis_scripts/load_xdf.py and
data/sub-P001_ses-day*/{{ao,bv,lsl,behavioral}}): Neuro Omega .mat exports, BrainAmp ExG BrainVision files (heart rate,
respiration, GSR; days 2-3), LSL XDF files and the CopyDraw behavioural YAML files (stylus traces per trial, block settings,
projection results). All SHA-256 sums matched MANIFEST.txt on download. Files changed for privacy are listed in
PRIVACY_CHANGES.tsv with their new SHA-256 (all other files are byte-identical to the release).

Not included: the benchtop tests (AO_test, Arduino_test, CT_* XDF/pickles), the c-VEP scalp-EEG speller data of three
healthy participants (data/data_p00*), and intermediate analysis results (data/*.hdf, *.npy, *.joblib). Get them from the
Radboud record.

## Privacy
No names or hospital identifiers were found. Dates were reduced to year and month (day set to 01, same byte length):
the export date in the .mat text headers and date strings in the XDF headers/marker strings and text files
(see PRIVACY_CHANGES.tsv). The BrainVision marker files already carry the placeholder date 1900-01-01. Age, sex and
handedness are those published in the paper.
""")
open(os.path.join(OUT, "CHANGES"), "w").write("1.0.0 2026-10-08\n  - BIDS-iEEG conversion of the patient recordings of Radboud Data Repository doi:10.34973/d214-m342 for NEMAR.\n")
open(os.path.join(OUT, ".bidsignore"), "w").write("code/\n")
if resid:
    print("RESIDUAL DATES", resid[:20])
    sys.exit(7)
print(json.dumps({"recordings": nrun, "hours": hrs, "privacy_changed_files": len(changed), "residual_dates_text": len(resid),
                  "rt_all": all(r["roundtrip_ok"] for r in rep["recordings"])}))
print("CONVERT_DONE")
