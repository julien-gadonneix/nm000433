# Mark channels at the int16 rail as bad

**Claim**: 26 channel entries in 10 day-4 runs have more than 10% of samples at |count| >= 32000 and are marked
`bad`: ECOG_1-4 and EOG_1-4 in stimoff run-01 (98.7%), LFP_5 and LFP_8 in the 9 closed-loop runs (17.1-92.7%).
Every other LFP channel reaches at most 3.9% with DBS off and 4.3% with open-loop DBS on.

**Source**
- PR nemarDatasets/nm000433#3, computed at head `8a066915`. Data: NEMAR copy of `origin/main` `8211a939` (only a CI
  workflow change since the PR base `344f172e`). The sha256 of each of the 164 files read (82 `.vhdr`, 82 `.eeg`)
  equals its SHA256E git-annex key at `8211a939`.
- Upstream: Radboud Data Repository, doi:10.34973/d214-m342, version 1 (Neuro Omega `.mat` exports; source file per
  run in `scans.tsv`). The converter wrote the `.mat` int16 arrays unchanged and checked exact equality on read-back
  (`roundtrip_ok` in `code/conversion_report.json`). The per-channel counts at 32767/-32768 recomputed here equal its
  `saturated_samples` in all 41 x 24 channels (last column).

**Method**: numpy reads the BrainVision INT_16 counts directly (no MNE, no scaling; the 38.14697265625/55 uV
resolution is not needed for counts). Per channel: fraction of all samples with |count| >= 32000; `bad` if above 10%.
Both cut-offs are the PR's rule, not taken from a source: 32000 is 97.7% of the int16 range, and 10% lies in the gap
between the highest unflagged channel (4.31%) and the lowest flagged one (17.1%). The README maxima are the
stimoff (DBS off) and stimon (open-loop DBS on) rows of the second table.

**Result** (script output, unedited)

| recording | samples | channels | channels > 10%: % of samples at abs count >= 32000 | max % of the other channels | counts at 32767/-32768 = conversion_report saturated_samples |
|---|---:|---:|---|---|---|
| sub-P001_ses-day2_task-copydraw_acq-stimoff_run-01 | 1208305 | 24 | 0: - | 0 | yes |
| sub-P001_ses-day2_task-copydraw_acq-stimoff_run-02 | 74497 | 24 | 0: - | 0 | yes |
| sub-P001_ses-day2_task-copydraw_acq-stimoff_run-03 | 838897 | 24 | 0: - | 0.11 (LFP_3) | yes |
| sub-P001_ses-day2_task-copydraw_acq-stimoff_run-04 | 777049 | 24 | 0: - | 0.12 (LFP_3) | yes |
| sub-P001_ses-day2_task-copydraw_acq-stimoff_run-05 | 849313 | 24 | 0: - | 0 | yes |
| sub-P001_ses-day2_task-copydraw_acq-stimoff_run-06 | 679081 | 24 | 0: - | 0 | yes |
| sub-P001_ses-day2_task-copydraw_acq-stimoff_run-07 | 1218961 | 24 | 0: - | 0 | yes |
| sub-P001_ses-day2_task-copydraw_acq-stimoff_run-08 | 4722873 | 24 | 0: - | 0 | yes |
| sub-P001_ses-day2_task-copydraw_acq-stimoff_run-10 | 7543297 | 24 | 0: - | 0 | yes |
| sub-P001_ses-day2_task-copydraw_acq-stimoff_run-12 | 8017921 | 24 | 0: - | 0 | yes |
| sub-P001_ses-day2_task-copydraw_acq-stimoff_run-14 | 6611801 | 24 | 0: - | 0 | yes |
| sub-P001_ses-day2_task-copydraw_acq-stimoffanalog_run-01 | 151039 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day2_task-copydraw_acq-stimoffanalog_run-02 | 9313 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day2_task-copydraw_acq-stimoffanalog_run-03 | 104863 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day2_task-copydraw_acq-stimoffanalog_run-04 | 97132 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day2_task-copydraw_acq-stimoffanalog_run-05 | 106165 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day2_task-copydraw_acq-stimoffanalog_run-06 | 84886 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day2_task-copydraw_acq-stimoffanalog_run-07 | 152371 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day2_task-copydraw_acq-stimoffanalog_run-08 | 590360 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day2_task-copydraw_acq-stimoffanalog_run-10 | 942913 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day2_task-copydraw_acq-stimoffanalog_run-12 | 1002241 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day2_task-copydraw_acq-stimoffanalog_run-14 | 826476 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day2_task-copydraw_acq-stimon_run-09 | 12553761 | 24 | 0: - | 2.38 (LFP_2) | yes |
| sub-P001_ses-day2_task-copydraw_acq-stimon_run-11 | 14421769 | 24 | 0: - | 1.63 (LFP_8) | yes |
| sub-P001_ses-day2_task-copydraw_acq-stimon_run-13 | 9552649 | 24 | 0: - | 2.84 (LFP_8) | yes |
| sub-P001_ses-day2_task-copydraw_acq-stimonanalog_run-09 | 1569221 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day2_task-copydraw_acq-stimonanalog_run-11 | 1802722 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day2_task-copydraw_acq-stimonanalog_run-13 | 1194082 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day3_task-copydraw_acq-stimoff_run-01 | 5653857 | 24 | 0: - | 0 | yes |
| sub-P001_ses-day3_task-copydraw_acq-stimoff_run-04 | 4933625 | 24 | 0: - | 0 | yes |
| sub-P001_ses-day3_task-copydraw_acq-stimoff_run-06 | 4814993 | 24 | 0: - | 0 | yes |
| sub-P001_ses-day3_task-copydraw_acq-stimoff_run-08 | 5024497 | 24 | 0: - | 0 | yes |
| sub-P001_ses-day3_task-copydraw_acq-stimoff_run-10 | 5310553 | 24 | 0: - | 0 | yes |
| sub-P001_ses-day3_task-copydraw_acq-stimoffanalog_run-01 | 706733 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day3_task-copydraw_acq-stimoffanalog_run-04 | 616704 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day3_task-copydraw_acq-stimoffanalog_run-06 | 601875 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day3_task-copydraw_acq-stimoffanalog_run-08 | 628063 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day3_task-copydraw_acq-stimoffanalog_run-10 | 663820 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day3_task-copydraw_acq-stimon_run-02 | 4811065 | 24 | 0: - | 3.14 (LFP_2) | yes |
| sub-P001_ses-day3_task-copydraw_acq-stimon_run-03 | 5203817 | 24 | 0: - | 2.87 (LFP_2) | yes |
| sub-P001_ses-day3_task-copydraw_acq-stimon_run-05 | 5142305 | 24 | 0: - | 2.86 (LFP_2) | yes |
| sub-P001_ses-day3_task-copydraw_acq-stimon_run-07 | 5078441 | 24 | 0: - | 2.86 (LFP_2) | yes |
| sub-P001_ses-day3_task-copydraw_acq-stimon_run-09 | 9876497 | 24 | 0: - | 1.46 (LFP_2) | yes |
| sub-P001_ses-day3_task-copydraw_acq-stimon_run-11 | 8668473 | 24 | 0: - | 1.64 (LFP_2) | yes |
| sub-P001_ses-day3_task-copydraw_acq-stimonanalog_run-02 | 601384 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day3_task-copydraw_acq-stimonanalog_run-03 | 650478 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day3_task-copydraw_acq-stimonanalog_run-05 | 642789 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day3_task-copydraw_acq-stimonanalog_run-07 | 634806 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day3_task-copydraw_acq-stimonanalog_run-09 | 1234563 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day3_task-copydraw_acq-stimonanalog_run-11 | 1083560 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day4_task-copydraw_acq-closedloop_run-02 | 5409121 | 24 | 2: LFP_5 41.4, LFP_8 75.1 | 2.06 (LFP_4) | yes |
| sub-P001_ses-day4_task-copydraw_acq-closedloop_run-03 | 4607705 | 24 | 2: LFP_5 17.1, LFP_8 31.0 | 2.30 (LFP_4) | yes |
| sub-P001_ses-day4_task-copydraw_acq-closedloop_run-04 | 4766185 | 24 | 2: LFP_5 36.7, LFP_8 68.6 | 2.29 (LFP_4) | yes |
| sub-P001_ses-day4_task-copydraw_acq-closedloop_run-05 | 4379817 | 24 | 2: LFP_5 39.2, LFP_8 66.9 | 1.72 (LFP_4) | yes |
| sub-P001_ses-day4_task-copydraw_acq-closedloop_run-10 | 5765137 | 24 | 2: LFP_5 44.9, LFP_8 78.9 | 1.99 (LFP_4) | yes |
| sub-P001_ses-day4_task-copydraw_acq-closedloop_run-11 | 3757057 | 24 | 2: LFP_5 43.4, LFP_8 82.9 | 2.43 (LFP_4) | yes |
| sub-P001_ses-day4_task-copydraw_acq-closedloop_run-12 | 5350873 | 24 | 2: LFP_5 49.5, LFP_8 90.8 | 2.14 (LFP_4) | yes |
| sub-P001_ses-day4_task-copydraw_acq-closedloop_run-13 | 4566849 | 24 | 2: LFP_5 50.7, LFP_8 91.7 | 2.06 (LFP_4) | yes |
| sub-P001_ses-day4_task-copydraw_acq-closedloop_run-14 | 4397593 | 24 | 2: LFP_5 51.0, LFP_8 92.7 | 2.14 (LFP_4) | yes |
| sub-P001_ses-day4_task-copydraw_acq-closedloopanalog_run-02 | 676141 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day4_task-copydraw_acq-closedloopanalog_run-03 | 575964 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day4_task-copydraw_acq-closedloopanalog_run-04 | 595774 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day4_task-copydraw_acq-closedloopanalog_run-05 | 547478 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day4_task-copydraw_acq-closedloopanalog_run-10 | 720643 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day4_task-copydraw_acq-closedloopanalog_run-11 | 469633 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day4_task-copydraw_acq-closedloopanalog_run-12 | 668860 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day4_task-copydraw_acq-closedloopanalog_run-13 | 570857 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day4_task-copydraw_acq-closedloopanalog_run-14 | 549700 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day4_task-copydraw_acq-stimoff_run-01 | 4247641 | 24 | 8: ECOG_1 98.7, ECOG_2 98.7, ECOG_3 98.7, ECOG_4 98.7, EOG_1 98.7, EOG_2 98.7, EOG_3 98.7, EOG_4 98.7 | 3.90 (LFP_1) | yes |
| sub-P001_ses-day4_task-copydraw_acq-stimoff_run-06 | 5059913 | 24 | 0: - | 0 | yes |
| sub-P001_ses-day4_task-copydraw_acq-stimoff_run-08 | 5112361 | 24 | 0: - | 0 | yes |
| sub-P001_ses-day4_task-copydraw_acq-stimoff_run-15 | 6998658 | 24 | 0: - | 0 | yes |
| sub-P001_ses-day4_task-copydraw_acq-stimoffanalog_run-01 | 530956 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day4_task-copydraw_acq-stimoffanalog_run-06 | 632490 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day4_task-copydraw_acq-stimoffanalog_run-08 | 639046 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day4_task-copydraw_acq-stimoffanalog_run-15 | 874833 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day4_task-copydraw_acq-stimon_run-07 | 5731257 | 24 | 0: - | 3.02 (LFP_8) | yes |
| sub-P001_ses-day4_task-copydraw_acq-stimon_run-09 | 6778105 | 24 | 0: - | 2.64 (LFP_4) | yes |
| sub-P001_ses-day4_task-copydraw_acq-stimon_run-16 | 4708041 | 24 | 0: - | 4.31 (LFP_2) | yes |
| sub-P001_ses-day4_task-copydraw_acq-stimonanalog_run-07 | 716408 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day4_task-copydraw_acq-stimonanalog_run-09 | 847264 | 16 | 0: - | 0 | n/a (analog inputs) |
| sub-P001_ses-day4_task-copydraw_acq-stimonanalog_run-16 | 588506 | 16 | 0: - | 0 | n/a (analog inputs) |

| acq | recordings | channels > 10% | max % of the other LFP channels | max % of the other channels |
|---|---:|---:|---:|---:|
| closedloop | 9 | 18 | 2.43 | 0.00 |
| closedloopanalog | 9 | 0 | 0.00 | 0.00 |
| stimoff | 20 | 8 | 3.90 | 0.00 |
| stimoffanalog | 20 | 0 | 0.00 | 0.00 |
| stimon | 12 | 0 | 4.31 | 0.00 |
| stimonanalog | 12 | 0 | 0.00 | 0.00 |

82 recordings, 26 channels with more than 10% of samples at abs count >= 32000.

**Coverage**: 82 of 82 recordings (41 Neuro Omega 22 kHz recordings and 41 analog-input recordings), all channels,
all samples. The 26 flagged entries and percentages equal the 26 `bad` rows and `status_description` values in the
10 `channels.tsv` files changed by the PR; no other row is `bad`.

Reproduce: `python code/receipts/2026-10-22_int16-rail-channels.py <nm000433 root>` (Python 3.12.15, numpy 2.5.4;
MNE not used).
