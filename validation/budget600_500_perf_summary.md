# ZINC22 RDKit/AMSOL Performance Summary

Benchmark date: 2026-05-29

Input set: 500 molecules from public ZINC22 `H24P100-N-oaa`, selected from IDs with official DB2 headers and matching tranche SMILES.

Run directory:

```text
/tmp/zinc22-official-compare/perf/run_H24P100_N_oaa_500_budget600_parallel8_20260529_01
```

Configuration:

```text
workers=8
RDKIT_CONF_BUDGET_BASE=600
RDKIT_CONF_TIMEOUT=120
pH=7.4
```

This is the full final-output path used by `submit/build-3d.bash`: tautomer/protomer generation at pH 7.4, protomer coalescing, stereo expansion, RDKit seed embedding, AMSOL solvation, RDKit conformer generation, torsion-strain annotation, `mol2db2` with terminal-H rotation enabled, and archive conversion to `solv`, `mol2`, `sdf`, `pdbqt`, and `db2.gz`.

## Topline

| Metric | Value |
|---|---:|
| Input molecules | 500 |
| Complete final outputs | 496 |
| Parallel wall time | 1184 s / 19.7 min |
| Throughput, attempted inputs | 1520 mol/hour |
| Throughput, complete outputs | 1508 mol/hour |
| Effective wall time | 2.37 s/input molecule |
| Max RSS per worker | 246 MB |

Four molecules failed during AMSOL solvation, before any final output could be written:

```text
ZINCog0000081ii0
ZINCog0000081k0G
ZINCog00000811mv
ZINCog0000080oRP
```

The AMSOL logs for those cases report an internal-coordinate failure involving atoms in an almost-straight geometry. The worker commands still exited successfully because the pipeline skips molecules without `output.solv`.

## Feature Time Breakdown

The table below uses summed worker CPU-stage time, not wall time. With 8 workers, the summed stage time is much larger than elapsed wall time. This is the right view for finding bottlenecks.

| Feature/stage | Summed worker time | Share | Approx. time per complete output | What it includes |
|---|---:|---:|---:|---|
| RDKit conformers | 8380.4 s | 93.3% | 16.90 s | ETKDG conformer generation, MMFF94s/UFF minimization, energy-window filtering, RMSD pruning, and rigid-fragment alignment/pinning before DB2. |
| Strain | 320.3 s | 3.6% | 0.65 s | Torsion-library strain lookup over the generated conformer ensemble. 23 molecules hit the nonfatal zero-strain fallback. |
| DB2 | 172.7 s | 1.9% | 0.35 s | `mol2db2` hierarchy construction, terminal-H rotation/reset expansion, clash filtering, and DB2 text generation. |
| AMSOL solvation | 89.5 s | 1.0% | 0.18 s | AMSOL water and hexane single-point solvation runs plus parsing into `output.solv` and AMSOL-charged `output.mol2`. |
| RDKit seed embed | 14.4 s | 0.16% | 0.03 s | One RDKit 3D seed conformer per protomer, written as Mol2 for AMSOL input. |
| Archive conversions | 0.6 s | 0.01% | 0.001 s | OpenBabel conversion of final Mol2 to SDF and PDBQT plus tar archive writes. |

Shell-level stages before 3D building were measured at one-second resolution and were `0s` per worker for this run: pH 7.4 tautomer/protomer generation, coalescing, and stereo expansion. They are included in the end-to-end wall time but are not material at this sample size with the current conservative RDKit protomer settings.

## Output Depth

RDKit retained conformers before DB2:

| Statistic | Count |
|---|---:|
| Mean | 499.9 |
| p50 | 549 |
| p75 | 587 |
| p90 | 596 |
| Max | 600 |

Final DB2 `sets` after terminal-H expansion:

| Statistic | Count |
|---|---:|
| Mean | 771.1 |
| p50 | 579 |
| p75 | 600 |
| p90 | 1689 |
| Max | 1797 |

`sets` is the final DOCK coordinate-state count in the DB2 file. It is not the same thing as unique heavy-atom conformers. For example, a molecule with 600 heavy conformers and one terminal hydrogen with three allowed states can produce about 1800 DB2 sets.

## Official ZINC Comparison

For the same official-covered IDs:

| Metric | Value |
|---|---:|
| Official mean DB2 sets, output IDs | 664.0 |
| Our mean DB2 sets, output IDs | 771.1 |
| Aggregate ratio, ours/official | 1.16x |
| Median per-molecule ratio, ours/official mean sets | 1.10x |
| p25/p75 per-molecule ratio | 0.93x / 2.43x |

Interpretation: `RDKIT_CONF_BUDGET_BASE=600` is close for tranche-average final DB2 depth. It slightly overproduces on aggregate, but the median molecule is near the official final set count. The wide per-molecule ratio is expected because official ZINC appears to use molecule-dependent effective caps, while this run uses one global RDKit budget.

## Practical Conclusions

1. The performance bottleneck is RDKit conformer generation. It accounts for about 93% of summed worker time at budget 600.
2. AMSOL is not the bottleneck in this configuration. It is about 1% of summed worker time, though it caused the four missing outputs.
3. DB2 construction is modest but not free. It is about 2% of summed worker time and includes terminal-H expansion, which is why final DB2 set counts can exceed retained RDKit conformer counts.
4. Budget 600 is the best current budget for matching tranche-average final ZINC-like output depth. Budget 900 overfilled the 100-molecule benchmark by about 1.5x on average.
5. The next quality issue to fix is AMSOL robustness for near-linear internal coordinates, because that is responsible for the observed missing outputs.

Parsed data files:

```text
/tmp/zinc22-official-compare/perf/run_H24P100_N_oaa_500_budget600_parallel8_20260529_01/parsed_summary.json
/tmp/zinc22-official-compare/perf/run_H24P100_N_oaa_500_budget600_parallel8_20260529_01/per_job_summary.tsv
/tmp/zinc22-official-compare/perf/run_H24P100_N_oaa_500_budget600_parallel8_20260529_01/official_set_comparison_summary.json
/tmp/zinc22-official-compare/perf/run_H24P100_N_oaa_500_budget600_parallel8_20260529_01/official_set_comparison.tsv
```
