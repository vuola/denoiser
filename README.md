# denoiser

Non-Local Means Denoising filter implemented for CT-scan image processing.

This README is the living project document. Each development step adds or
updates a chapter here; a LaTeX/PDF report is produced when a main phase is
completed.

| Chapter | Content | Status |
|---|---|---|
| 1 | Project plan | Reviewed; question 3 of section 1.8 waiting for approval |
| 2 | Existing toolkit and the pipeline for new scans | Draft, waiting for review |

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
  tests/data/      64³ test scans cut from the CT scans (tracked) and the script that cuts them
  run/             benchmark scripts and their results, figures
  latex/           phase reports (.tex) and their figures
  third_party/     vendored dependencies (test framework, header-only readers)
  data/            full CT scans and historical tools (local only; data/README.txt is tracked)
  ORIGINAL_CODE/   untouched reference material (local only, not in git)
```

Test data: the full scans live in `data/RAW/` and are ignored by git.
`python tests/data/make_crops.py` cuts a 64×64×64 cube from each into
`tests/data/<scan>_c64/` as a miniature scan folder in the lab's format (`.vol`
plus `.vgi`, 0.5 MB, tracked), so tests run from a fresh clone. Details are in
`data/README.txt`.

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
| 0.2 | Build the original code unchanged as the reference program, with a small converter that feeds it a test scan as a temporary float file and discards the file afterwards | Reference binary usable on `.vol` test scans |
| 0.3 | Test volumes: the 16-bit test scans in `tests/data/` (done), a synthetic fibre-like phantom with known added noise, and a small brute-force NLM for tiny volumes | Inputs in the lab's format; ground truth for quality (PSNR) and correctness |
| 0.4 | Baseline: reference outputs, run time and memory on the test volumes | Numbers every later step is compared with |

#### Phase 1 — Optimise for 3D CT data, input interface, pre-processing

| Step | Content | Result |
|---|---|---|
| 1.1 | **Scan reader.** Read a scan folder as the CT lab delivers it: dimensions, sample type and voxel size from the `.vgi` file, voxels from the `.vol` file. The same description (dimensions, sample type, byte order, header offset) can be given by hand for other raw files | Scans load without ImageJ |
| 1.2 | **Writer.** Output as `.vol` plus `.vgi` in the input's sample type, so the result opens wherever the input does | Output in the lab's format |
| 1.3 | **Pre-processing.** Region-of-interest crop in voxel coordinates, intensity mapping to the working range and back, NaN/Inf handling | Replaces `PickCube.ijm` and `scaleToRange.py` |
| 1.4 | **Noise estimate.** Sigma estimated from the data with the method of `estimate_noise.py`, or given by hand | Replaces `estimate_noise.py` |
| 1.5 | **Filter parameters.** Patch radius, search radius and `h` settable; the sigma table stays as the default. Patch-size normalisation and kernel selectable: original behaviour or as in the IPOL paper | Deviations from the papers become explicit choices |
| 1.6 | **Numerical accuracy.** Measure the float32 integral-image error; fix it if it matters | Tested accuracy statement |
| 1.7 | **Speed and memory.** Symmetric shifts, persistent threads, smaller per-thread storage, block I/O, vectorisable prefix sums | Benchmark table against the baseline |
| 1.8 | **One-command run and run record.** `denoiser run <scan folder>` chains steps 1.1–1.5 and the filter, and writes a record of every setting next to the output | The pipeline of chapter 2 |
| 1.9 | Phase 1 report | `latex/phase1.tex`, PDF |

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
| 1 | CT data format | Revised: scans arrive as in `data/RAW`, a 16-bit unsigned little-endian `.vol` with a `.vgi` metadata file. Input and output are `.vol` only; no float files are stored. See chapter 2 |
| 2 | Real test data | Scans B (`Sample1_18p15_5mu`) and D (`Sample2_5mu`), 5 µm voxels; 64³ test scans in `tests/data/`. Samples A and C (20 µm voxels) are dropped |
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

---

## 2. Existing toolkit and the pipeline for new scans

Goal: a new scan from the CT lab is filtered with one documented command,
without ImageJ, and the run can be repeated exactly.

### 2.1 What the CT lab delivers

Each scan is a folder such as `data/RAW/Sample2_5mu/` (phoenix nanotom,
datos|x reconstruction):

| File | Content | Used by the pipeline |
|---|---|---|
| `.vol` | Voxels: headerless, 16-bit unsigned, little-endian, x fastest. 1142 × 1142 × 1000 voxels, 2.6 GB | Yes |
| `.vgi` | Text: volume size, bits per element, data type, voxel size (0.005 mm) | Yes |
| `.pcr` | Text: reconstruction settings, and the physical values that counts 0 and 65535 stand for (`Min`, `Max`) | Recorded in the run record |
| `.pca`, `.pcj`, `.pcp`, `.dtxml` | Acquisition geometry and per-projection log | No |

The `.vgi` file does not state the byte order. Little-endian is confirmed by
the check in section 2.3. The file size equals 1142 × 1142 × 1000 × 2 bytes, so
there is no header.

The two scans have different `Min`/`Max` in their `.pcr` files (upper value
0.303 against 0.242), so the same count does not mean the same attenuation in
both scans.

### 2.2 Existing tools in `data/`

| Tool | What it does | Problem for repeatability |
|---|---|---|
| ImageJ import and "Save As > Raw data" | Reads the `.vol`, converts to 32-bit float, writes big-endian raw | Menu actions; nothing records what was done |
| `macros/PickCube.ijm` | Cuts a centred cube; size and centre typed into a dialog | Dialog input; slice numbers are 1-based |
| `scripts/scaleToRange.py` | Maps an input range linearly to 0–255 and clips | Forces big-endian input, so it cannot read the `.vol` directly; needs numpy |
| `scripts/estimate_noise.py` | Estimates sigma from 20 slices with five estimators; recommends the median of the neighbour-difference and local-residual estimates | Works on the float file; needs numpy |
| `ORIGINAL_CODE/denoiserVersionA.cpp` | The filter | Float big-endian input only; dimensions from the file name |

The other macros and scripts (median, Frangi and tubeness filters, CPU/GPU
comparison scripts) are not part of the filtering chain.

### 2.3 How the existing samples were made (verified)

Samples B and D were recomputed from the two scans in `data/RAW` and compared
with the existing float files, all 15.6 million voxels each:

- **Crop:** the 250³ cube starts at voxel (446, 446, 374), counted from zero.
  That is `PickCube.ijm` with the default centre. The z start is one voxel
  lower than the x and y starts because ImageJ counts slices from 1.
- **Scaling:** value = count × 255 / 65535. The full 16-bit range was mapped,
  not a narrower window, and nothing was clipped.
- **Agreement:** the largest difference is 0.000008, which is float rounding.

So the range of samples B and D (for example 30–226) is the range of the data
itself. Sample A reaches exactly 0 and 255, which would mean saturated counts
in that scan or a different scaling. The scans of samples A and C are not in
`data/RAW`, so this could not be checked.

### 2.4 Target pipeline

One program, `denoiser`, with the steps as subcommands:

| Command | Purpose |
|---|---|
| `denoiser info <scan>` | Print what the `.vgi` says, check the file size, print value statistics |
| `denoiser noise <scan>` | Estimate sigma |
| `denoiser run <scan>` | Read, crop if asked, scale, estimate sigma unless given, filter, write the result and the run record |

`<scan>` is the scan folder or its `.vgi` file. Filtering a new scan is then:

```
1. Copy the scan folder into data/RAW/
2. denoiser run data/RAW/<scan name>
```

Options of `run`, all with defaults:

| Option | Default | Meaning |
|---|---|---|
| `--roi x0 y0 z0 w h d` | whole volume | Region to filter, zero-based voxel coordinates |
| `--centre-cube n` | off | Shortcut for a centred n³ region |
| `--scale lo hi` | `0 65535` | Counts mapped to the working range 0–255; the default reproduces the existing samples |
| `--sigma s` | estimated | Noise sigma in working-range units |
| `--out dir` | `run/<scan name>/` | Output folder |

Output folder:

| File | Content |
|---|---|
| `<scan name>_nlm.vol` | Filtered volume, 16-bit unsigned little-endian, same scaling as the input |
| `<scan name>_nlm.vgi` | Metadata for the filtered volume (size, type, voxel size) |
| `run_record.txt` | Program version, full command line, input file size and checksum, region, scaling, sigma and how it was obtained, filter parameters, run time |

Repeating a run means running the command line stored in `run_record.txt`.

### 2.5 Design choices

- **One C++ program instead of Python scripts plus ImageJ.** Reading, cropping
  and scaling are needed inside the filter anyway, because it processes the
  volume block by block. Doing them there removes numpy and ImageJ from the
  chain and leaves one tool to document and test. The existing scripts and
  macros stay in `data/` as reference.
- **No intermediate float file.** The filter reads 16-bit blocks directly from
  the `.vol`, converts each block in memory and writes 16-bit blocks. A full
  scan as float would be 5.2 GB on disk and in the old chain existed twice.
- **Zero-based voxel coordinates everywhere,** so a region means the same in
  the tool, in the record and in the test data.
- **Sigma method kept, then examined.** Step 1.4 first reproduces what
  `estimate_noise.py` recommends, so old and new runs are comparable. Its
  estimators assume noise that is independent between neighbouring voxels,
  which reconstruction filtering violates, so step 1.4 also checks it on the
  phantom with known noise.
- **Full scans need the speed work.** The original filter took about 14 hours
  for a volume of this size. `run` works on a region from the start; full-scan
  runs become practical with step 1.7.

### 2.6 Decisions from the review

| # | Question | Decision |
|---|---|---|
| 1 | Output format | `.vol` plus `.vgi` only. No float files are written or stored |
| 2 | `data/macros`, `data/scripts` | Historical reference, not tracked. Everything in `data/` except `README.txt` is ignored by git |
| 3 | Test material | Scans B and D (5 µm voxels) are enough. A and C (20 µm) are dropped. Future scans will have 1 µm voxels |

Consequences:

- **Regression baseline without stored float files.** The original program
  reads only float input. Step 0.2 therefore converts a test scan to a
  temporary float file, runs the original on it and deletes the file. The
  16-bit test scans were checked against the earlier float crops before those
  were removed: they agree to 0.000008.
- **16-bit output loses nothing that matters.** One 16-bit count is 0.004 on
  the 0–255 working scale, far below the noise.
- **1 µm scans change the scale of structures.** A fibre that is a few voxels
  wide at 5 µm is five times wider in voxels at 1 µm, and the noise per voxel
  will differ. The patch and search radii of the sigma table are in voxels, so
  the defaults may not suit 1 µm data. Step 1.5 makes the radii settable, the
  run record stores the voxel size, and the defaults are reviewed when the
  first 1 µm scan arrives.
