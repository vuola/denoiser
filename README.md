# denoiser

Non-Local Means Denoising filter implemented for CT-scan image processing.

This README is the living project document. Each development step adds or
updates a chapter here; a LaTeX/PDF report is produced when a main phase is
completed.

| Chapter | Content | Status |
|---|---|---|
| 1 | Project plan | Reviewed; question 3 of section 1.8 waiting for approval |

---

## 1. Project plan

### 1.1 Goal

Turn the 3D Non-Local Means (NLM) filter in `ORIGINAL_CODE` into a maintained,
tested and documented tool for denoising 3D CT volumes, in three main phases:

1. **Optimise** the original implementation for 3D CT data and give it a
   flexible interface for describing the input data format and pre-processing
   the data before filtering.
2. **Refactor** the filter into small, readable, fully documented subfunctions
   that are unit tested, comparable across computing platforms, and replaceable
   one by one with GPU versions later.
3. **Document** the work continuously in this README and as a LaTeX/PDF report
   at the end of each main phase.

Writing the GPU versions themselves is outside this plan; phase 2 only prepares
for them.

### 1.2 Starting point

`ORIGINAL_CODE` holds:

| File | What it is |
|---|---|
| `denoiserVersionA.cpp` | 3D NLM filter in C++ by Samuel Väisänen (2024), a port of the 2D ImageJ plugin |
| `NLMeansDenoising_.java` | The 2D ImageJ plugin by Behnel and Wagner (GPL-3) that the C++ is based on |
| `Fast_nonlocal_filtering_applied_to_electron_cryomicroscopy.pdf` | Darbon et al. 2008: the fast integral-image algorithm |
| `Non-Local Means Denoising.pdf` | Buades, Coll, Morel 2011 (IPOL): noise threshold and parameter table |
| `CORRESPONDENCE` | Emails with the usage notes of the C++ code |

Findings from the code review that shape the plan:

- **Algorithm.** The code is the Darbon algorithm (integral image of squared
  differences per search shift, constant-time patch distance, weight
  `h²/(h² + d)`, centre voxel weighted with the maximum weight) with the IPOL
  noise threshold `2σ²` and the IPOL grey-image parameter table.
- **Input.** Only one format: headerless raw, 32-bit big-endian float,
  x fastest, dimensions parsed from the file name (`name_WxHxD.raw`). Values are
  expected in 0–255 because the parameters are chosen from the absolute value of
  sigma. NaN and Inf are not handled.
- **Deviations from the papers.** The patch distance is not divided by the patch
  size, so filtering is stricter than in IPOL; the IPOL `h` values (tuned for an
  exponential kernel and 2D patches) are used with a rational kernel and 3D
  patches.
- **Numerical risk (not yet tested).** The integral image is accumulated in
  32-bit floats over a 256³ block, so patch distances may lose precision.
- **Speed.** About 14 hours for a 1000×1000×900 volume on a multi-core lab
  machine, according to the author. Known inefficiencies: all search shifts are
  computed although half would do by symmetry; threads are created per block;
  each thread holds a full 256³ result array; file I/O is voxel by voxel.

### 1.3 Repository layout

The layout follows the `NMSC-assignment` project folder.

```
denoiser/
  README.md        living project document (this file)
  makefile         build, test, benchmark and report targets
  include/         public headers (function interface documentation lives here)
  src/             implementation and command-line program
  tests/           unit tests, regression tests, cross-platform comparison tool
  tests/data/      64³ crops of the CT volumes (tracked) and the script that cuts them
  run/             benchmark scripts and their results, figures
  latex/           phase reports (.tex) and their figures
  third_party/     vendored dependencies (test framework, header-only readers)
  data/            full CT volumes, 250³ (local only; data/README.txt is tracked)
  ORIGINAL_CODE/   untouched reference material (local only, not in git)
```

Test data: four real CT volumes (samples A–D, 250×250×250, 32-bit big-endian
float, values 0–255, 62.5 MB each) live in `data/` and are ignored by git.
`python tests/data/make_crops.py` cuts the centred 64×64×64 cube of each into
`tests/data/` (1 MB each, tracked), so tests run from a fresh clone. Details
are in `data/README.txt`.

`include/` is not in the list of required folders but the NMSC project has it;
it is proposed here because phase 2 puts the interface documentation in
headers.

### 1.4 Tools and conventions

| Topic | Choice | Reason |
|---|---|---|
| Build | GNU make, one `makefile` at the root | Preferred; works on this Windows machine (MinGW g++ 16, make 4.4) and on Linux |
| CMake | Not used for now | CPU code does not need it. CUDA on Linux builds with make and `nvcc`. CMake is added only if a CUDA build on Windows is wanted, where `nvcc` needs the Microsoft compiler instead of MinGW |
| GPU target | NVIDIA/CUDA (later, outside this plan) | Decided in review |
| Reports | LaTeX built with MiKTeX on Windows | Small basic install that fetches missing packages on demand |
| Language | C++17 | Supported by current CPU and GPU compilers |
| Threads | OpenMP | Same as the NMSC project; simpler than hand-made thread pools |
| Unit tests | doctest (single header in `third_party/`) | No installation, runs with `make test` |
| Code comments | Doxygen-style blocks on every function | Phase 2 requirement; can generate a reference manual |
| Floating point | No `-ffast-math`; fixed summation order in the reference path | Needed for comparing results across platforms |

The new code derives from Samuel Väisänen's C++ and, through it, from the
GPL-3 ImageJ plugin of Behnel and Wagner. The repository therefore stays GPL-3
and the source files credit both.

### 1.5 Workflow

- Work proceeds in the numbered steps of section 1.6.
- Each step ends with a README chapter (what was done, results, open points)
  and a stop. The next step starts only after you have read the chapter and
  commented or approved.
- One commit per approved step, on a branch per main phase.
- At the end of each main phase a LaTeX report is written from the README
  chapters and built to PDF.

### 1.6 Phases and steps

#### Phase 0 — Foundation

| Step | Content | Result |
|---|---|---|
| 0.1 | Create the folder layout, `makefile`, `.gitignore` fix | `make`, `make test`, `make clean` work |
| 0.2 | Build the original code unchanged as the reference program | Reference binary |
| 0.3 | Synthetic test volumes (fibre-like phantom with known added noise) and a small brute-force NLM for tiny volumes | Ground truth for quality (PSNR) and correctness |
| 0.4 | Baseline: reference outputs, run time and memory on the test volumes | Numbers every later step is compared with |

#### Phase 1 — Optimise for 3D CT data, input interface, pre-processing

| Step | Content | Result |
|---|---|---|
| 1.1 | **Input description.** Dimensions, sample type (8/16-bit integer, 32-bit float), byte order, header offset and axis order given on the command line or in a small text file, instead of the file name | Raw CT volumes load without conversion in ImageJ |
| 1.2 | **File formats.** Raw plus the formats your CT data comes in (see open questions) | Readers and writers |
| 1.3 | **Pre-processing.** Intensity mapping to the working range and back, region-of-interest crop, NaN/Inf handling, output in the input's sample type | Documented, optional steps |
| 1.4 | **Noise estimate.** Sigma from a homogeneous region or from a histogram fit, as an alternative to typing it in | `sigma` no longer has to be found by hand |
| 1.5 | **Filter parameters.** Patch radius, search radius and `h` settable; the sigma table stays as the default. Patch-size normalisation and kernel selectable: original behaviour or as in the IPOL paper | Deviations from the papers become explicit choices |
| 1.6 | **Numerical accuracy.** Measure the float32 integral-image error; fix it if it matters | Tested accuracy statement |
| 1.7 | **Speed and memory.** Symmetric shifts, persistent threads, smaller per-thread storage, block I/O, vectorisable prefix sums | Benchmark table against the baseline |
| 1.8 | Phase 1 report | `latex/phase1.tex`, PDF |

Acceptance for phase 1: with the original settings the output matches the
reference within a stated tolerance; every intended difference is listed and
justified.

#### Phase 2 — Refactor for readability, tests and GPU replacement

| Step | Content | Result |
|---|---|---|
| 2.1 | **Subfunction design.** Proposed chain: read block → pad borders → squared difference for one shift → 3D prefix sum → patch distance from 8 corners → weight → accumulate → normalise → write block. Each works on plain contiguous arrays plus a small dimensions struct, with no hidden state | Interface chapter for review before coding |
| 2.2 | **Implementation** of the readable CPU reference with complete interface documentation; speed is secondary | Reference backend |
| 2.3 | **Unit tests** per subfunction: analytic cases, comparison with the brute-force NLM, edge cases | `make test` |
| 2.4 | **Cross-platform comparison.** Fixed test vectors; each subfunction's output written with platform tag, checksum and statistics; a tool compares two platforms with exact match for integers and a stated tolerance for floats | `make compare` |
| 2.5 | **Backend interface.** Subfunctions selected through one table so a GPU version can replace one at a time and be checked with the step 2.4 tool | Documented replacement procedure |
| 2.6 | Phase 2 report | `latex/phase2.tex`, PDF |

#### Phase 3 — Documentation

| Step | Content | Result |
|---|---|---|
| 3.1 | Continuous: README chapter per step (runs through phases 0–2) | This file |
| 3.2 | Final report: theory, differences from the papers, user manual, interface reference, test and benchmark results | `latex/report.tex`, PDF |

### 1.7 Risks

| Risk | Handling |
|---|---|
| Real CT data has no noise-free ground truth | Synthetic phantom with known noise for quality figures; real crops for regression and benchmarks |
| `ORIGINAL_CODE` is not in git, so the reference program cannot be rebuilt from a fresh clone | Reference outputs for the tracked crops are stored in `tests/data/` in step 0.4 |
| Optimisation changes results unnoticed | Baseline outputs from step 0.4 checked in every step |
| Float results differ between compilers and platforms | Tolerances defined in step 2.4, no unsafe math flags |
| Readable reference becomes too slow for full volumes | Phase 1 optimised path kept as a second backend |

### 1.8 Decisions from the plan review

| # | Question | Decision |
|---|---|---|
| 1 | CT data format | Headerless raw, 32-bit big-endian float, as in `data/`. Step 1.1 still adds the general description (sample type, byte order, header offset); other file formats in step 1.2 wait until a need appears |
| 2 | Real test data | Samples A–D in `data/`, 64³ crops in `tests/data/` |
| 3 | Default filter behaviour | Proposed: keep the original behaviour as default, see below. Waiting for approval |
| 4 | GPU family | NVIDIA/CUDA |
| 5 | LaTeX | MiKTeX on this Windows PC |
| 6 | Original code and emails | Removed from the remote repository; `ORIGINAL_CODE` is local only |

**Recommendation for question 3.** Keep the behaviour of the original code
(which is also the behaviour of the ImageJ plugin) as the default through
phases 0–2, and add the IPOL patch-size normalisation and exponential kernel as
options in step 1.5. Reasons:

- The regression baseline of step 0.4 only works if the default reproduces the
  original output.
- Sigma values and results you already have for these samples were obtained
  with this behaviour, as were those of other users of the ImageJ plugin.
- The IPOL variant is not known to be better here: its parameters were tuned
  for 2D photographs, not for 3D CT patches.

Step 1.5 then compares both variants on the phantom (PSNR against the known
noise-free volume) and on the real crops. The default changes only if that
comparison shows a clear gain, and you decide on the basis of those numbers.

### 1.9 Open questions

1. Which noise sigma did you use for samples A–D? It is needed for the baseline
   in step 0.4.
2. Does the Ubuntu server have an NVIDIA GPU? If so it becomes the second test
   platform for the cross-platform comparison of step 2.4.
