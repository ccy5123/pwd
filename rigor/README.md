# rigor/ — partition research through the sci-adk methodology

This directory uses the [**sci-adk**](https://github.com/ccy5123/sci-adk) rigor engine
(*agents propose; the engine judges by frozen criteria; record ≠ belief; null results
are first-class*) to take a step on the partition-toolkit's open problem and record it
as a re-verifiable, tamper-evident run. Two runs:

1. **`partition-model-cv`** — the research step (a cross-validated multi-surrogate
   model), built to follow the sci-adk methodology properly (see the correction note).
2. **`partition-xtb-ufz`** — the baseline: a faithful, verifiable restatement of the
   toolkit's existing single-surrogate README numbers.

Both re-derive without an LLM (`sci-adk verify`, exit 0) and are tamper-evident.

> ### Methodology correction (v1 → v2 — recorded honestly, not hidden)
> The **first** model run set its thresholds *after* inspecting the CV numbers
> (0.60 hugged the observed albumin/muscle 0.644/0.646; a "min ≥ 0" bar was picked to
> single out muscle-chicken) and labelled the result `confirmatory`. That is **HARKing**
> — exactly what sci-adk's frozen-Spec / anti-HARKing principle forbids — and it skipped
> the prior-work checkpoint. This v2 corrects all three:
> 1. **Principled, pre-committed bars** (not tuned to the data): `useful` = mean CV R² ≥
>    **0.50** (majority of out-of-sample variance); `reliable` = the 95% across-seed CV-R²
>    interval lies entirely **above the no-skill null 0** (a textbook significance bar).
> 2. **Honest mode = `exploratory`.** The dataset has been examined, so genuine
>    confirmatory pre-registration on it is no longer possible; true confirmation needs
>    an untouched validation set (recorded as a limitation in the Spec).
> 3. **Prior work recorded** (`LITERATURE` + `PRIOR_WORK_DECISION` Evidence), reframing
>    the novelty (below).

## Prior work — what's actually novel here

A web search (recorded in the run's Evidence) shows multi-descriptor modelling of these
phases is **established PP-LFER work**, including **Endo, Brown & Goss 2013** ([ES&T,
10.1021/es401772m](https://pubs.acs.org/doi/10.1021/es401772m)) — *this dataset's own
source paper* — which models membrane/storage lipid, albumin and protein vs water at
**RMSE 0.32–0.53 log units** using measured Abraham descriptors. So "a multi-descriptor
model predicts these phases" is **not novel**, and PP-LFER is the real benchmark.

**Honest contribution of this run:** a *fully in-silico, xTB-only* surrogate model (no
measured Abraham descriptors — only structure → GFN1/2-xTB) recovers the protein-phase
predictability that the toolkit's **single uncalibrated surrogate** route reported as
unpredictable. It is **comparable-to-slightly-worse than PP-LFER** (CV RMSE ≈ 0.51–0.58
here; not a like-for-like comparison), so it is a *cheap-descriptor* alternative, not an
accuracy advance.

## Result (autonomous verdicts; `sci-adk verify` exit 0, all reproduced)

| Phase | n | Baseline R² (toolkit single surrogate) | **CV R² mean** | 95% across-seed interval | CV RMSE | useful (≥0.50) / reliable (CI>0) |
|-------|--:|--:|--:|--:|--:|---|
| membrane | 207 | 0.808 | **0.906** | [0.875, 0.921] | 0.58 | SUPPORTED / SUPPORTED |
| storage | 247 | 0.837 | **0.898** | [0.720, 0.951] | 0.57 | SUPPORTED / SUPPORTED |
| albumin | 83 | 0.307 | **0.644** | [0.585, 0.683] | 0.51 | SUPPORTED / SUPPORTED |
| muscle (fish) | 45 | −0.008 | **0.648** | [0.358, 0.767] | 0.56 | SUPPORTED / SUPPORTED |
| muscle (chicken) | 46 | −0.027 | 0.646 | **[−0.045, 0.774]** | 0.52 | SUPPORTED / **PROPOSED** |

**What the engine judged:**
- vs the toolkit's single-surrogate baseline, the in-silico multi-surrogate model makes
  **albumin** and **muscle (fish)** useful *and* reliably-above-no-skill, and lifts the
  lipids further (membrane/storage CV R² ~0.90). The README's near-zero protein R² was
  an artefact of the single uncalibrated surrogate, not a real ceiling.
- **muscle (chicken)** is `useful` on average (0.65) but its `reliable` claim is
  **PROPOSED, not SUPPORTED**: the 95% interval includes 0 (min −0.05 at n=46), so skill
  is *not reliably distinguishable from no-skill* — an honest "inconclusive", recorded as
  a first-class outcome rather than averaged away.

## How this run follows the sci-adk methodology

sci-adk is a **referee, not a player**: it builds the rigor kernel and *borrows* the
experiment as a capability from the in-session agent. This run is driven through the
kernel's **designed entry point**, honoring its core disciplines:

| Discipline | How this run honors it |
|---|---|
| **Driven via the real seam** | `run_model.py` registers the capability in sci-adk's adapter registry and invokes **`sci-adk run --capability partition-model-cv`** (`sci_adk.cli.main`) — not a hand-rolled compiler call. The engine renders every verdict. |
| **Engine = source of truth** | SUPPORTED/PROPOSED are the `DecisionEngine`'s, from frozen rules (e.g. muscle-chicken `reliable` = PROPOSED is the interval rule mapping a CI-containing-null → NEUTRAL, not a choice of mine). |
| **Anti-HARKing (S1/S3)** | Bars (0.50, null 0) are principled and pre-committed; hypotheses honestly `exploratory` (the dataset was examined). |
| **Anti-method-shopping (F3)** | The **full single-surrogate field — all 72 surrogates' R² per phase — is recorded** as an append-only OBSERVATION (`*-surrogate-scan`), so the toolkit's "best surrogate" winner cannot be cherry-picked: winners *and* losers are on the record. |
| **Prior work recorded** | `LITERATURE` + `PRIOR_WORK_DECISION` Evidence; novelty reframed against PP-LFER. |
| **Generalization gate** | This is a **second domain** (after T-1) plugged into the kernel **with no kernel edit** — the roadmap's generalization gate (`design/adoption-roadmap.md` §6). |
| **In-session only** | No autonomous LLM call; the numeric/interval rules judge for free. |
| **Docker execution seam** | The CV runs through an executor seam mirroring sci-adk's `T1Executor`: Docker-backed when a daemon + the `sci-adk-partition` image are present (provenance stamps the **image id**), else an in-process fallback. The pure science (`partition_cv.py`, no sci-adk dep) is identical either way — only the container is a seam. Each Evidence records which executor ran. |

This sandbox has no Docker daemon, so the committed run used the **in-process fallback**
— recorded honestly in provenance (`executor=in-process, fallback=docker daemon not
running`). On a Docker host (`docker build -t sci-adk-partition rigor/` first), the same
`python rigor/run_model.py` runs the CV inside the container and stamps the image id.

## Files

- `proposal.md` — four-pane pre-registration of the baseline study (EN + KR).
- `partition_capability.py` — baseline capability (→ `partition-xtb-ufz`).
- `partition_model_capability.py` — the corrected model capability (→ `partition-model-cv`):
  frozen Spec (per phase: `useful` threshold + `reliable` interval, exploratory) + an
  `ExperimentFn` that fits the ridge model, runs repeated CV on the real data, records
  prior art, and emits `data_source="measured"` Evidence.
- `run_rigor.py` / `run_model.py` — zero-LLM compile drivers.
- `runs/.../` — generated, verifiable artifacts (`spec.json` · `evidence/*.json` ·
  `claims/*.json` · `checkpoints/` · `paper/draft.md`).

## Reproduce

```bash
pip install "git+https://github.com/ccy5123/sci-adk@6f8376d407a49240153e115a94a3dbe7bc14ed0c"
pip install numpy pandas scikit-learn

python rigor/run_model.py                       # the research step
sci-adk verify rigor/runs/partition-model-cv    # REPRODUCED per claim; exit 0
```

Run the experiment inside the Docker seam (on a Docker-enabled host) so provenance
stamps the container image id:

```bash
docker build -t sci-adk-partition rigor/        # numpy/pandas/scikit-learn, version-pinned
python rigor/run_model.py                        # auto-detects the daemon + image -> runs in-container
# provenance environment then reads: executor=docker:sci-adk-partition, image_id=<sha>
```

## Remaining deviations / limitations (recorded in the Spec + here)

- **Exploratory, not confirmatory** on this dataset (it has been examined). A genuine
  confirmatory test needs an untouched external validation set or a locked hold-out.
- **Below the PP-LFER benchmark** (RMSE 0.32–0.53); the contribution is cheap in-silico
  descriptors, not accuracy.
- **Execution seam:** the Docker-backed executor is implemented (`sci-adk-partition`
  image + `partition_cv.py` running in-container, image id stamped); this sandbox has no
  Docker daemon so the committed run used the in-process fallback (recorded in
  provenance). The Docker path is validated up to the `docker run` call (sci-adk's own
  `DockerExecutor`); only an actual container launch is unexercised here.
- **muscle sets are small** (n ≈ 45) → the muscle-chicken reliability is inconclusive;
  more measured data or feature reduction is the next step.
