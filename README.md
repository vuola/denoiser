# denoiser

Non-Local Means Denoising filter implemented for CT-scan image processing.

This README is the living project document. Each development step adds or
updates a chapter here; a LaTeX/PDF report is produced when a main phase is
completed.

| Chapter | Content | Status |
|---|---|---|
| 1 | Project plan | Reviewed and approved |
| 2 | Existing toolkit and the pipeline for new scans | Reviewed and approved |
| 3 | Step 0.1: folder layout, build and tests | Done, committed |
| 4 | Step 0.2: the original filter as reference program | Done, committed |
| 5 | Research objectives and the requirements they set | Reviewed and approved |
| 6 | Step 0.3: phantom, scoring and the brute-force filter | Done, committed |
| 7 | Step 0.4: baseline (run time, memory, accuracy, phantom scores) | Done, committed |
| 8 | Step 1.1: scan reader and the command `denoiser info` | Done, committed |
| 9 | Step 1.2: scan writer | Done, committed |
| 10 | Step 1.3: pre-processing, and the brightness drift of the present scans | Done, committed |
| 11 | Step 1.4: noise estimate, and the noise of the present scans | Done, committed |
| 12 | Step 1.5: filter settings, and which of them suit fibre rods | Done, committed |
| 13 | Step 1.6: the fast filter and its accuracy | Done, committed |

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

Phase 2 prepares for GPU versions of the subfunctions. Writing them was first
left outside this plan; the size of the coming scans (section 5.6) makes them
necessary for production runs, so they are added as phase 5.

The filter serves a measurement task: finding small air bubbles in the epoxy of
glass-fibre rods. Chapter 5 describes that task and the requirements it sets,
and phase 4 adds the tools for it after the filter work.

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
  tests/reference/ scripts that run the original filter and compare volumes
  tests/phantom/   scripts that make the synthetic rod phantom and score results against it
  tests/data/      64³ test scans cut from the CT scans (tracked), the script that cuts them,
                   and the original program's results for them (original_outputs/)
  run/             benchmark scripts and their results (run/baseline/), figures
  latex/           phase reports (.tex) and their figures
  third_party/     vendored dependencies (test framework, header-only readers)
  data/            full CT scans, larger crops (data/crops/) and historical tools
                   (local only; data/README.txt is tracked)
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
| Reports | LaTeX built on the second machine (TeX Live with `pdflatex` and `latexmk`); the PDF is copied back | Already installed there, so the Windows machine needs no LaTeX. MiKTeX is dropped |
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
- One commit per approved step, on `main` (decided in review, section 3.6).
- At the end of each main phase a LaTeX report is written from the README
  chapters and built to PDF.

### 1.6 Phases and steps

#### Phase 0 — Foundation

| Step | Content | Result |
|---|---|---|
| 0.1 | Create the folder layout, `makefile`, `.gitignore` fix | `make`, `make test`, `make clean` work |
| 0.2 | Build the original code unchanged as the reference program, with a small converter that feeds it a test scan as a temporary float file and discards the file afterwards | Reference binary usable on `.vol` test scans |
| 0.3 | Test volumes: the 16-bit test scans in `tests/data/` (done); a synthetic rod phantom with known noise and a known material label for every voxel (section 5.4); a scoring script (PSNR and share of misclassified voxels); a small brute-force NLM for tiny volumes | Inputs in the lab's format; ground truth for quality, classification and correctness |
| 0.4 | Baseline: reference outputs, run time and memory on the test volumes; phantom scores of the original filter and of a median and a mean filter along z (section 5.5) | Numbers every later step is compared with |

#### Phase 1 — Optimise for 3D CT data, input interface, pre-processing

| Step | Content | Result |
|---|---|---|
| 1.1 | **Scan reader.** Read a scan folder as the CT lab delivers it: dimensions, sample type and voxel size from the `.vgi` file, voxels from the `.vol` file. The same description (dimensions, sample type, byte order, header offset) can be given by hand for other raw files | Scans load without ImageJ |
| 1.2 | **Writer.** Output as `.vol` plus `.vgi` in the input's sample type, so the result opens wherever the input does | Output in the lab's format |
| 1.3 | **Pre-processing.** Region-of-interest crop in voxel coordinates, intensity mapping to the working range and back, NaN/Inf handling. The mapping is linear and nothing is clipped by default (requirement R2). A measurement of how brightness drifts from the rod centre to the rim, as a basis for deciding on a correction (R4) | Replaces `PickCube.ijm` and `scaleToRange.py` |
| 1.4 | **Noise estimate.** Sigma estimated from the data with the method of `estimate_noise.py`, or given by hand | Replaces `estimate_noise.py` |
| 1.5 | **Filter parameters.** Patch radius and search radius settable separately for x, y and z (R5), `h` settable; the sigma table with its cubes stays as the default. Patch-size normalisation and kernel selectable: original behaviour or as in the IPOL paper | Deviations from the papers become explicit choices |
| 1.6 | **Fast filter and numerical accuracy.** The filter computed with running sums, for all settings of step 1.5, checked against the brute-force filter in double and single precision (scope extended in the review of step 1.5, section 12.9) | Tested accuracy statement |
| 1.7 | **Speed and memory.** Symmetric shifts, persistent threads, smaller per-thread storage, block I/O, vectorisable prefix sums | Benchmark table against the baseline |
| 1.8 | **One-command run and run record.** `denoiser run <scan folder>` chains steps 1.1–1.5 and the filter, and writes a record of every setting next to the output | The pipeline of chapter 2 |
| 1.9 | Phase 1 report | `latex/phase1.tex`, PDF |

Acceptance for phase 1 (revised in the review of step 0.3, section 6.5): with
the original settings the output matches the **brute-force filter** within a
tight, stated tolerance on the small test volumes. The difference from the
original program is measured and reported, and explained by the limited
accuracy of the original's 32-bit running sums. Every intended difference is
listed and justified.

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
| 3.1 | Continuous: README chapter per step (runs through phases 0–2 and 4) | This file |
| 3.2 | Final report: theory, differences from the papers, user manual, interface reference, test and benchmark results. Written after phase 4 | `latex/report.tex`, PDF |

#### Phase 4 — Material classification and bubble analysis

Added after the review of the research objectives (chapter 5). The steps are
an outline; they are detailed when phase 2 ends and the first 1 µm scan is
available.

| Step | Content | Result |
|---|---|---|
| 4.1 | **Look at the data.** Histograms of the raw and the filtered volume, brightness against distance from the rod axis, values in the bore hole and in the air outside the rod | Numbers for contrast, noise and drift |
| 4.2 | **Classification.** Label every voxel as glass, epoxy or air from the filtered volume. Output: a label volume as `.vol` plus `.vgi` | Label volume; scored on the phantom |
| 4.3 | **Masks.** Epoxy mask shrunk away from fibre edges by a settable number of voxels; reference region inside the bore hole | Voxel sets for step 4.4 |
| 4.4 | **Bubble search in the raw data.** Raw values of the masked epoxy voxels compared with the epoxy level and with the air level from the bore hole; connected groups of low voxels listed with position, size, mean value and significance | Candidate list and a record of the settings |
| 4.5 | **Validation.** Phantom with voids of known size and position: share found against void size, false alarms per volume | Smallest bubble size that can be found reliably |
| 4.6 | Phase 4 report | `latex/phase4.tex`, PDF |

#### Phase 5 — GPU versions

Added because a 1 µm scan is about three times the present size (section 5.6,
requirement R9). Outline only: the subfunctions of phase 2 are replaced one at
a time by CUDA versions on the second machine, each checked against the CPU
reference with the comparison tool of step 2.4. The volume is processed block
by block, so it need not fit in GPU memory. Whether phase 5 comes before or
after phase 4 is decided at the end of phase 2, from the run times measured by
then.

### 1.7 Risks

| Risk | Handling |
|---|---|
| Real CT data has no noise-free ground truth | Synthetic phantom with known noise for quality figures; real crops for regression and benchmarks |
| `ORIGINAL_CODE` is not in git, so the original program cannot be rebuilt from a fresh clone | The reference is the brute-force filter, which is in git (section 6.5). Outputs of the original program for the tracked crops are stored in `tests/data/` in step 0.4, for comparison only |
| Optimisation changes results unnoticed | Brute-force results on the small test volumes checked in every step |
| Float results differ between compilers and platforms | Tolerances defined in step 2.4, no unsafe math flags |
| Readable reference becomes too slow for full volumes | Phase 1 optimised path kept as a second backend |
| The 5 µm scans resolve a fibre with only 4 voxels, so they say little about how the filter behaves on 1 µm data | The phantom imitates the 1 µm case; defaults are reviewed when the first 1 µm scan arrives |
| The filter removes or blurs the small bubbles that are the target | Bubbles are searched in the raw data; the filtered volume is used only for the material labels (section 5.2) |

### 1.8 Decisions from the plan review

| # | Question | Decision |
|---|---|---|
| 1 | CT data format | Revised: scans arrive as in `data/RAW`, a 16-bit unsigned little-endian `.vol` with a `.vgi` metadata file. Input and output are `.vol` only; no float files are stored. See chapter 2 |
| 2 | Real test data | Scans B (`Sample1_18p15_5mu`) and D (`Sample2_5mu`), 5 µm voxels; 64³ test scans in `tests/data/`. Samples A and C (20 µm voxels) are dropped |
| 3 | Default filter behaviour | Approved: the original behaviour stays the default; the IPOL variant becomes an option in step 1.5. Reasons below |
| 4 | GPU family | NVIDIA/CUDA |
| 5 | LaTeX | Revised: TeX Live on the second machine replaces MiKTeX on the Windows PC |
| 6 | Original code and emails | Removed from the remote repository; `ORIGINAL_CODE` is local only |

**Decision for question 3 (approved).** Keep the behaviour of the original code
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

1. ~~Which noise sigma did you use for samples A–D?~~ Answered: 17 on the
   0–255 scale. See section 4.6.
2. ~~Is a second test platform with an NVIDIA GPU available?~~ Answered: yes, a
   Linux machine with an NVIDIA GPU and the CUDA compiler. It is the second
   platform for the cross-platform comparison of step 2.4 and the CUDA platform
   later.
3. ~~Geometry of the coming 1 µm scan.~~ Answered: a region inside the rod,
   up to three times the present voxel count. Section 5.6 recommends
   2000 × 2000 × 975 voxels and lists three points to settle with the CT lab.
   The bore hole is 900 µm in diameter (confirmed).

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

---

## 3. Step 0.1 — Folder layout, build and tests

Status: done and committed.

### 3.1 What was added

| Path | Content |
|---|---|
| `makefile` | Build, test, clean and remote-test targets |
| `.gitattributes` | Text files use LF line endings on every platform; `.vol`, `.raw`, `.pdf`, `.png` are never converted |
| `include/version.hpp`, `src/version.cpp` | The first library function, `denoiser::version()` |
| `src/main.cpp` | The command-line program. It recognises `info`, `noise`, `run`, `--version` and `--help`; the three commands only report that they are not implemented yet |
| `tests/test_main.cpp` | Test runner (doctest supplies `main`) |
| `tests/test_setup.cpp` | Two tests of the setup: the version text is well formed, and the two tracked test scans are present with the right size |
| `third_party/doctest/` | doctest 2.4.12, one header, with its MIT licence text |
| `run/`, `latex/` | Empty for now (a `.gitkeep` file keeps them in git) |

### 3.2 Commands

Run from the repository root.

| Command | Effect |
|---|---|
| `make` | Builds the program `build/denoiser` |
| `make test` | Builds and runs the unit tests (`build/run_tests`) |
| `make clean` | Removes `build/` and `build-san/` |
| `make test SANITIZE=1` | Builds into `build-san/` with address and undefined-behaviour checks, then runs the tests. Linux only |
| `make remote-test` | Sends the last commit to the second machine and runs the tests there (section 3.4) |
| `make help` | Lists the targets |

### 3.3 Build choices

- **Output goes to `build/`.** Object files and programs stay out of the source
  folders, and `build/` is already ignored by git.
- **Flags: `-O2 -std=c++17 -Wall -Wextra`, with OpenMP.** There is no
  `-march=native` and no `-ffast-math`, so that two machines compile the same
  arithmetic (section 1.4). Whether `-march=native` is worth adding for speed
  is measured in step 1.7.
- **Header dependencies are tracked automatically.** Changing a header
  rebuilds exactly the files that include it.
- **New source files need no makefile change.** Every `src/*.cpp` except
  `main.cpp` is library code and is linked into both the program and the
  tests; every `tests/*.cpp` is part of the test program.
- **On Windows the compiler's runtime is linked into the program (`-static`).**
  Without it the test program did not start from Git Bash: Git for Windows
  ships an older `libstdc++` DLL that was found before the compiler's own.
- **On Windows `make` needs the Unix tools of Git for Windows on the PATH.** It
  was tested from Git Bash and from PowerShell.

### 3.4 The second machine

A Linux machine with an NVIDIA GPU is the second platform (section 1.9). The
repository refers to it only by the SSH host alias `gpu`; the alias is defined
in `~/.ssh/config` on the development machine, outside the repository.

One-time setup:

1. Add a `Host gpu` entry with the machine's address and user name to
   `~/.ssh/config`. Key login must work without a prompt.
2. Clone the repository on that machine into `~/denoiser`.
3. On the development machine: `git remote add gpu gpu:denoiser`.

After that, `make remote-test` pushes the last commit to the branch
`remote-test` of that clone, checks it out there and runs `make test`.
Uncommitted changes are not sent, and nothing goes to GitHub. Another alias or
folder can be given as `make remote-test REMOTE=name REMOTE_DIR=folder`.

### 3.5 Results

| Platform | Compiler | `make` | `make test` | `make test SANITIZE=1` |
|---|---|---|---|---|
| Windows 11, development machine | MinGW g++ 16.2 | builds | 2 test cases, 12 checks, all pass | not available |
| Linux, second machine | g++ 15.2 | builds | 2 test cases, 12 checks, all pass | all pass, no sanitizer report |

Both builds are free of compiler warnings. After the step was committed,
`make remote-test` ran with the same result on the second machine.

### 3.6 Open points

- ~~LaTeX is not installed on either machine yet.~~ Resolved: LaTeX is now
  installed on the second machine (`pdflatex` and `latexmk` found there).
  Reports are built on that machine; a make target that sends the report
  source there and fetches the PDF is added with the first report (step 1.9).
- ~~Branch per main phase or `main`?~~ Decided: work continues on `main`.
  There is one developer and each step is reviewed before it is committed.
  Section 1.5 is changed accordingly.

---

## 4. Step 0.2 — The original filter as reference program

Status: done and committed.

### 4.1 What was added

| Path | Content |
|---|---|
| `makefile`, target `reference` | Compiles `ORIGINAL_CODE/denoiserVersionA.cpp` without any change into two programs (section 4.3) |
| `tests/reference/run_reference.py` | Runs a reference program on a scan folder (`.vgi` + `.vol`) and writes the result as `.vol` + `.vgi` |
| `tests/reference/compare_volumes.py` | Compares two volumes voxel by voxel: number of differing voxels, largest and RMS difference |
| `tests/reference/slice_images.py` | Writes PNG images of slices for visual inspection: input, filtered, difference (section 4.7) |
| `tests/reference/one_worker.cpp` | A few lines that make the original source start one worker thread instead of many |

The original source stays outside git. To build the reference on another
machine, copy that one file into `ORIGINAL_CODE/` of the clone there.

### 4.2 Commands

```
make reference
python tests/reference/run_reference.py tests/data/Sample2_5mu_c64 --sigma 10 --one-worker
python tests/reference/compare_volumes.py a.vol b.vol
python tests/reference/slice_images.py tests/data/Sample2_5mu_c64 build/reference/Sample2_5mu_c64_ref_s10.vol
```

The Python command is `python` on Windows and `python3` on Linux. On Windows
`python3` is only a placeholder that reports that Python is missing.

`run_reference.py` does what section 2.6 announced: it converts the scan to a
temporary float file (value = count × 255 / 65535), runs the original program
on it, converts the result back to 16-bit and deletes the temporary files.
`--sigma` is on the 0–255 scale. Results go to `build/reference/` unless
`--out` names another folder; `--keep-float` also keeps the unrounded result
as a `.f32` file for exact comparisons. The script holds the volume in memory
and is meant for the test scans only.

### 4.3 Finding: the original program is not repeatable

The original source compiles unchanged and without warnings with both
compilers (as C++20, which it needs). Two runs of it on the same input do not
give the same output:

| Comparison (Sample2 test scan, sigma 10) | Voxels differing | Largest difference |
|---|---|---|
| Original, run 1 against run 2, float result | 86 % | 0.0002 on the 0–255 scale |
| Original, run 1 against run 2, 16-bit result | 0.5 % | 1 count |

The cause is in the threading. The program starts one worker per processor
(minus one). The workers take search shifts from a shared counter, so which
worker handles which shift depends on timing, and each worker keeps its own
partial sums. The partial sums are added at the end, and floating-point
addition gives slightly different results for different groupings.

A program whose output changes between runs cannot serve as a regression
reference. The target `make reference` therefore builds two programs from the
same unchanged source:

| Program | Workers | Use |
|---|---|---|
| `build/denoiser_reference` | processors − 1, as the author built it | Run time of the original |
| `build/denoiser_reference_1w` | 1 | Repeatable reference for comparing results |

The one-worker program is made by linking `tests/reference/one_worker.cpp`,
which replaces the standard library's processor count with 1. The filter
source itself is not touched.

### 4.4 Results

All runs: 64³ test scans, sigma 10 (patch radius 1, search radius 10).

| Check | Result |
|---|---|
| One-worker program, run 1 against run 2, both machines, both scans | Identical, bit for bit |
| One-worker program, Windows/MinGW g++ 16.2 against Linux/g++ 15.2, both scans | Identical, bit for bit |
| One-worker against multi-worker program, float result | Largest difference 0.0055 on the 0–255 scale |
| One-worker against multi-worker program, 16-bit result | 5.8 % of voxels differ, by at most 2 counts |
| Result values outside 0–65535 after conversion, NaN or infinite values | None |

Run time of the original program on one 64³ test scan:

| Machine | Multi-worker | One worker |
|---|---|---|
| Development machine (11 workers) | 12–13 s | 38–51 s |
| Second machine (15 workers) | 7.6 s | 27–29 s |

The times are single runs, noted only to show that the reference is usable;
step 0.4 measures properly.

### 4.5 What this means for later steps

- **The one-worker program is the repeatable form of the original.** Its
  output is the same on every run and on both platforms. It was first meant
  as the regression reference; step 0.3 found that its 32-bit running sums
  are not accurate enough for that, and the brute-force filter took over that
  role (section 6.5). The one-worker program remains the way to state how far
  new code is from the original.
- **The 0.0055 difference is a first measurement for step 1.6.** One worker
  adds 9260 terms per voxel in one long chain; several workers add shorter
  chains. The gap between the two is rounding error of 32-bit sums, and it
  exceeds one 16-bit count (0.004). It is small against the noise (sigma about
  10), but it shows that the summation order must be fixed by design in the
  new code, not left to thread timing.
- **Many workers give little speed-up.** Fifteen workers are only 3.7 times
  faster than one. Likely causes are the items already listed in section 1.2
  (threads created per block, 200 MB of result storage per worker); step 1.7
  looks into it.

### 4.6 Open points

- The sigma of 10 used here was a placeholder. The earlier filtering of these
  samples used sigma 17 (answered in review). The two values fall in
  different rows of the original parameter table: up to 15 the patch radius
  is 1, from 15 to 30 it is 2. The runs above therefore did not exercise the
  patch size that was actually used. Step 0.4 takes the baseline at sigma 17
  and at sigma 10, so both rows are covered.
- A rough neighbour-difference estimate on the two test scans gives 12 and
  14, below the table's limit of 15, while 17 is above it. The estimate thus
  decides the patch size; step 1.4 looks at this when it reproduces the
  earlier estimate.
- ~~The `.vgi` writer shared with `make_crops.py` handles cubes only.~~
  Resolved in step 0.3: it accepts any volume size. The real writer is still
  step 1.2.

### 4.7 Visual inspection

`tests/reference/slice_images.py` replaces looking at results in ImageJ. It
takes the input volume and the filtered volume (each a scan folder or a `.vol`
file with its `.vgi`) and writes one PNG per axis into `build/slices/`. Each
image shows the same slice three times, left to right:

| Panel | Content |
|---|---|
| Input | The slice before filtering |
| Filtered | The slice after filtering, on the same grey scale |
| Difference | Filtered minus input; mid-grey is "no change" |

A difference panel that looks like pure noise means the filter removed only
noise. Outlines of fibres in it mean the filter also changed the structure.

Options: `--axis x|y|z|all` (default all), `--slice N` (default the middle
slice), `--window LOW HIGH` in counts (default the 0.5 and 99.5 percentiles of
the input slice), `--gain` for the difference panel (default 1), `--out`.
The script reads only the requested slices; one slice of a full
1142×1142×1000 scan takes about 4 seconds. It needs only the Python standard
library.

First look at the reference result (Sample2 test scan, sigma 10): the filtered
panel is clearly smoother and the difference panel shows noise without
visible structure. The 64³ test scans are too small to judge more than that
by eye: a slice of the full scan shows fibres running along z, but a 64-voxel
window holds only a few of them. A larger crop for visual checks, kept local
and not tracked, is proposed for step 0.4.

---

## 5. Research objectives and the requirements they set

Status: reviewed and approved. This chapter records
what the filter is for and what follows from it for the plan of section 1.6.

### 5.1 The measurement task

The samples are pultruded glass-fibre rods: glass fibres running along the
rod, with epoxy between them. The question to answer is whether the epoxy
holds small air bubbles.

| Property | Value | Source |
|---|---|---|
| Fibre diameter | about 20 µm | Project owner |
| Glass share of the volume | about 60 % (datasheet: more than 60 %) | Project owner, supplier datasheet |
| Fibre direction | Along the scan's z axis. Fibre positions drift sideways from slice to slice and form gently twisted, rope-like bundles; no sharp bends | Project owner, from the 5 µm scans |
| Earlier experience | A median filter along z over 33 voxels (165 µm), with no window in x and y, improved the contrast of the 5 µm scans clearly | Project owner; parameters from the macros in `data/macros` (section 5.5) |
| Coming scans | 1 µm voxels; rod turned down from 18.15 mm to 5 mm diameter; a hole drilled along the axis as an air reference; scanned region inside the rod, up to three times the present voxel count | Project owner; see section 5.6 |

What these numbers mean in voxels, assuming evenly spaced fibres:

| Quantity | Size | At 5 µm voxels | At 1 µm voxels |
|---|---|---|---|
| Fibre diameter | 20 µm | 4 | 20 |
| Distance between fibre centres | about 25 µm | 5 | 25 |
| Epoxy gap between two neighbouring fibres | about 5 µm | 1 | 5 |
| Largest epoxy pocket between three fibres | about 8 µm across | 1–2 | 8 |
| Window along z of the earlier median filter | 165 µm | 33 | 165 |

Two things follow directly. At 5 µm the epoxy is at most two voxels wide, so
the present scans cannot show glass and epoxy as separate regions; the 1 µm
scans can. And a bubble that sits in the epoxy without pushing fibres aside is
at most a few micrometres across: a handful of voxels even at 1 µm.

### 5.2 The analysis chain

| Stage | Input | Output | Phase |
|---|---|---|---|
| 1. Denoise | Raw scan | Filtered volume | 1–2 |
| 2. Classify | Filtered volume | Label per voxel: glass, epoxy or air | 4 |
| 3. Mask | Labels | Epoxy voxels away from fibre edges; air voxels in the bore hole | 4 |
| 4. Search | **Raw** scan, masks | Bubble candidates, judged against the epoxy level and the bore-hole air level | 4 |

The search goes back to the raw values on purpose. Non-local means averages
voxels whose surroundings look alike. A small bubble has few look-alikes, so
the filter either leaves it noisy or blends it into the epoxy. The filtered
volume is therefore used only to decide where the epoxy is.

### 5.3 Requirements

| # | Requirement | Reason | Where |
|---|---|---|---|
| R1 | Filter quality is judged by how well glass and epoxy can be told apart afterwards: share of misclassified voxels on the phantom, and whether fibre edges stay in place. PSNR is reported too but does not decide | Stage 2 is what the filter is for | Steps 0.3, 0.4, 1.5 |
| R2 | Voxel values stay linear from input to output: no clipping, no contrast stretch, same count scale in and out | The bore-hole calibration compares raw and filtered values with a reference level | Steps 1.2, 1.3 |
| R3 | The filtered volume has exactly the same voxel grid as the input, or a region of it with its offset in the run record | Masks made from the filtered volume are applied to the raw scan | Steps 1.3, 1.8 |
| R4 | Slow brightness drift across the rod is measured, and corrected if it would move the glass/epoxy threshold | A full slice of the 5 µm scans is not equally bright at centre and rim. The bore hole gives the air level at the centre only; air outside the rod gives a second point at the rim | Step 1.3 (measure), 4.1–4.2 |
| R5 | Search radius and patch radius are settable per axis | Section 5.5 | Step 1.5 |
| R6 | The epoxy mask can be shrunk by a settable number of voxels | A voxel on a fibre edge holds part glass and part epoxy and reads like "unusual epoxy", just as a bubble edge does | Step 4.3 |
| R7 | The bubble search states the smallest bubble it can find reliably | A bubble of N voxels stands out by contrast × √N / noise. Contrast comes from the bore hole, noise from the raw epoxy voxels. Single-voxel bubbles will probably stay below the limit | Steps 4.4, 4.5 |
| R8 | Filter defaults are reviewed on the first 1 µm scan | Radii are in voxels; structures are five times larger in voxels than in the scans used so far | Steps 1.5, 1.8 |
| R9 | A scan of 3.9 × 10⁹ voxels (three times the present size) is filtered in a practical time, block by block, on a GPU with less memory than the scan | Section 5.6 | Step 1.7, phase 2, phase 5 |

### 5.4 The phantom of step 0.3

The phantom imitates the 1 µm case, because the tracked 5 µm test scans cannot
show the structures that matter (section 5.1).

- Cylinders of 20 voxels diameter along z, filling about 60 % of the volume,
  with a small tilt and sideways drift so that they are not perfectly aligned.
- Epoxy between them, at a lower value than glass.
- An air-filled hole along the axis.
- Voids of several sizes in the epoxy, at known positions.
- Added noise of known strength; optionally a slow brightness drift from
  centre to rim.
- Written as `.vol` plus `.vgi`, together with the noise-free volume and the
  true label of every voxel.

A second, coarser version with fibres 4 voxels thick imitates the present
5 µm scans. Size, grey levels and noise strength are fixed in step 0.3; the
grey levels start from the histograms of the real scans.

### 5.5 Using the fibre direction in the filter

The structure changes slowly along z and quickly across it. The filter can
use that in two ways, both cheap in this algorithm:

| Change | Effect | Cost |
|---|---|---|
| Search window longer along z, shorter across | More true look-alikes per voxel | Lower than now: 11 × 11 × 41 shifts are about half of 21 × 21 × 21 |
| Patch longer along z | Patch comparison is less noisy, without mixing across fibre edges | None: the patch sum takes 8 look-ups for any box size |

Decision proposed: step 1.5 makes the radii settable per axis (R5). The
default stays the cube of the original code, so the regression baseline holds.
Whether a stretched window becomes the default is decided from the phantom
scores.

The filter used earlier on these rods is added to the baseline of step 0.4 as
a competitor. According to the macros kept in `data/macros`
(`median_filter_current_image.ijm`, `open_and_median_filter.ijm`) it was
ImageJ's "Median 3D" with radii x = 0, y = 0, z = 16 (an earlier setting,
commented out, was z = 15). That is a median over a single column of 33 voxels
along z, 165 µm at 5 µm voxels, with no window in x and y. The baseline holds
two filters along z with that window: the median, and the mean for comparison.

Both are far cheaper than NLM. Their weakness is the twist: within the window
a fibre moves sideways, so they mix glass and epoxy at fibre edges, while NLM
still finds matching voxels in the neighbouring columns. The median handles
this better than the mean, because it keeps an edge as long as the fibre
stays in the column for more than half the window. The phantom scores show
how much NLM gains over them, and at what tilt.

### 5.6 Size of the 1 µm scan

Answered in review: the scanned region is widened in x and y so that it holds
enough rod around the bore hole, and shortened in z if the voxel count gets
too large. It stays inside the 5 mm rod sideways. The total may grow to three
times the present count: 3 × 1142 × 1142 × 1000 = 3.9 × 10⁹ voxels.

The bore hole is 900 µm in diameter (confirmed in review), which is 900 voxels.

| x × y × z (voxels = µm) | Rod around the hole, each side | Rod volume in the scan | Share of voxels that is rod | `.vol` size |
|---|---|---|---|---|
| 1600 × 1600 × 1500 | 350 µm (14 fibre spacings) | 2.9 mm³ | 75 % | 7.7 GB |
| **2000 × 2000 × 975** | **550 µm (22 fibre spacings)** | **3.3 mm³** | **84 %** | **7.8 GB** |
| 2400 × 2400 × 675 | 750 µm (30 fibre spacings) | 3.5 mm³ | 89 % | 7.8 GB |

**Recommended: 2000 × 2000 × 975**, with the bore hole centred.

- The ring of rod around the hole is 22 fibre spacings wide. That holds
  several thousand fibres per slice and leaves material well away from the
  drilled surface, where drilling damage could be mistaken for bubbles.
- 975 slices are about the present depth. The filter needs depth only for its
  search along z (some tens of voxels), so more slices mainly add more of the
  same fibres, while more width adds new fibres and new epoxy pockets.
- Going wider than 2000 gains little rod volume (3.5 against 3.3 mm³) and
  costs a third of the depth.

Points to settle with the CT lab:

1. **Detector width.** The metadata of the present scans show a detector of
   1152 × 1152 pixels, which gives the 1142 voxels across. If a single scan
   cannot be wider than that, 2000 voxels need a finer detector mode or
   several scans joined together. If neither is possible, the fallback is
   about 1150 × 1150 × 2900 with the hole off-centre, so that the field of
   view holds a segment of the hole and otherwise rod. A centred 900 µm hole
   would leave a ring of only 125 µm.
2. **Scan inside the sample.** The field of view (2 mm) is smaller than the
   sample (5 mm). Reconstructions of such scans typically show a brightness
   drift towards the edge of the field of view. That makes requirement R4
   more important, and the air outside the rod is no longer available as a
   second reference: the bore hole is the only air in the picture.
3. **Hole size.** The air reference needs some thousands of voxels, not
   600 000 per slice. A smaller hole (200–300 µm), if it can be drilled,
   would leave more of the scan for rod.

Consequences for this project:

- **Run time.** The filter's time grows in proportion to the voxel count. The
  original program would need about 40 hours for such a scan (14 hours for
  the present size, according to its author). Production runs therefore need
  the speed work of step 1.7 and a GPU version (requirement R9).
- **Memory.** The scan is 7.8 GB as 16-bit counts and twice that as floats,
  more than the memory of the GPU on the second machine. The filter must work
  block by block on the GPU too; the block design of phase 2 allows that.
- **No fixed size.** Nothing in the code may assume 1142 × 1142 × 1000 or a
  cube; readers and writers take the size from the `.vgi` file.

### 5.7 Changes to the plan

| Place | Change |
|---|---|
| Section 1.1 | Pointer to this chapter |
| Step 0.3 | Phantom as in section 5.4, with true labels; scoring script |
| Step 0.4 | Phantom scores of the original filter and of a median and a mean filter along z |
| Step 1.3 | Linear, unclipped mapping stated as a requirement; brightness drift measured |
| Step 1.5 | Radii settable per axis |
| Phase 4 (new) | Classification, masks, bubble search, validation, report |
| Phase 5 (new) | GPU versions of the subfunctions, needed for production runs on 1 µm scans |
| Sections 1.9, 4.6 | Sigma 17 recorded; the baseline of step 0.4 includes it |
| Step 3.2 | Final report moves after phase 4 |
| Sections 1.7, 1.9 | Two risks and the geometry question added |

---

## 6. Step 0.3 — Phantom, scoring and the brute-force filter

Status: done and committed.

### 6.1 What was added

| Path | Content |
|---|---|
| `tests/phantom/make_phantom.py` | Makes a synthetic rod scan with known noise, the noise-free volume and the true material of every voxel (section 6.2) |
| `tests/phantom/score.py` | Scores a volume against a phantom: error, share of misclassified voxels, void contrast (section 6.3) |
| `include/nlm_bruteforce.hpp`, `src/nlm_bruteforce.cpp` | The filter written straight from its definition, in double precision, with radii per axis; and the parameter table of the original program as a function (section 6.4) |
| `tests/test_nlm_bruteforce.cpp` | Eight unit tests of the two functions |
| `tests/reference/bruteforce_main.cpp`, `make bruteforce` | The brute-force filter behind the command line of the original program -> `build/denoiser_bruteforce` |
| `tests/reference/run_reference.py` | New option `--bruteforce`; volumes no longer need to be cubes |
| `tests/data/make_crops.py` | The shared `.vgi` writer accepts any volume size |

### 6.2 The phantom

```
python tests/phantom/make_phantom.py --preset fine
python tests/phantom/make_phantom.py --preset coarse
```

Output below `build/phantom/` (not tracked; the same seed gives the same
files on every platform):

| Folder | Content |
|---|---|
| `<name>/` | The noisy scan as `.vol` + `.vgi`, 16-bit; `phantom.txt` with every setting, voxel counts and file checksums; `voids.csv` with centre, diameter and voxel count of every void |
| `<name>_clean/` | The same volume without noise |
| `<name>_labels/` | Material of every voxel, 8-bit: 0 air in the bore hole, 1 epoxy, 2 glass, 3 air in a void; 16 is added for an edge voxel that mixes two materials |

| Setting | `fine` (imitates 1 µm scans) | `coarse` (imitates 5 µm scans) |
|---|---|---|
| Size | 160 × 160 × 96 | 64 × 64 × 64 |
| Fibre diameter | 20 voxels | 4 voxels |
| Bore-hole diameter | 40 voxels | 16 voxels |
| Void diameters, 6 voids of each | 1, 3, 5, 7 voxels | 1, 3 voxels |
| Twist about the rod axis | 0.0005 rad per slice | 0.0025 rad per slice |
| Glass share | 60 % requested; 59 % and 60 % obtained outside the bore hole | |
| Sideways tilt | 0.02 in x and 0.01 in y, voxels per slice | |
| Noise sigma | 14 on the 0–255 scale (3598 counts) | |

Every setting can be changed on the command line (`--size`, `--sigma`,
`--twist`, `--void-diameters`, `--drift` for a brightness drift from centre to
rim, and so on; see `--help`).

How it is built:

- **Fibres** are circles on a hexagonal grid whose spacing gives the glass
  share. Each is moved by a random offset small enough that fibres never
  overlap. From slice to slice all fibres rotate about the rod axis and shift
  sideways together, so they stay apart at every depth. The bore hole stays on
  the axis.
- **Edges** are one voxel wide: a voxel cut by a fibre edge or by the bore
  wall gets a mixed value.
- **Voids** are balls centred on a pure-epoxy voxel. Only epoxy voxels inside
  the ball become air, so the fibres around a large void cut into it, as they
  would around a real bubble that does not push fibres aside. A void of
  7 voxels diameter keeps 120–170 of the ball's 179 voxels.
- **Noise** is independent from voxel to voxel and very close to a normal
  distribution. It is made from integer random numbers only, which is what
  makes the files identical on both platforms.

Grey levels, in counts:

| Material | Level | Where it comes from |
|---|---|---|
| Air | 18 000 | The `.pcr` files of the two scans give the attenuation that counts 0 and 65535 stand for; zero attenuation falls at count 17 900 and 15 400 |
| Epoxy | 26 000 | Assumption, see below |
| Glass | 38 000 | Assumption, see below |

The mean count of the two test scans is 32 900 and 27 600. That mean is a mix
of 60 % glass and 40 % epoxy, which the 5 µm scans do not show separately. The
split assumes that glass attenuates 2.5 times as strongly as epoxy. This ratio
is an assumption of mine, not a measurement; it sets the contrast between
glass and epoxy (12 000 counts, 3.3 times the noise sigma) and so how hard the
phantom is. It is to be replaced by the measured ratio from the first 1 µm
scan.

What the phantom leaves out: real CT noise is correlated between neighbouring
voxels, real edges are blurred over more than one voxel, and fibres vary in
diameter. These can be added when the first 1 µm scan shows how much they
matter.

### 6.3 The scores

```
python tests/phantom/score.py build/phantom/phantom_fine <volume to score>
```

| Score | Meaning |
|---|---|
| `rmse`, `psnr` | Error against the noise-free volume, on the 0–255 scale; PSNR = 20 log₁₀(255 / rmse) |
| `wrong_all` | Share of voxels classified wrongly. Each voxel is classified from its value alone, with fixed thresholds halfway between the air, epoxy and glass levels |
| `wrong_interior` | The same among voxels that hold one material only |
| `wrong_edge` | The same among edge voxels |
| `void_contrast` | How much of the step from epoxy down to air is left in the void voxels: 1 = as dark as air, 0 = filled in to the epoxy level |

A table of true material against classified material is printed as well.

First scores, with the original program (the proper baseline is step 0.4):

| Phantom | Volume scored | PSNR | `wrong_all` | `wrong_interior` | `wrong_edge` | `void_contrast` |
|---|---|---|---|---|---|---|
| fine | Noisy scan, unfiltered | 25.2 dB | 12.1 % | 10.2 % | 25.6 % | 0.98 |
| fine | Original filter, sigma 14 | 38.1 dB | 1.5 % | 0.06 % | 11.8 % | 0.07 |
| fine | Original filter, sigma 17 | 39.9 dB | 1.2 % | 0.05 % | 9.3 % | 0.25 |
| coarse | Noisy scan, unfiltered | 25.2 dB | 18.1 % | 8.2 % | 26.2 % | 1.04 |
| coarse | Original filter, sigma 14 | 31.8 dB | 8.5 % | 0.17 % | 15.2 % | −0.57 |
| coarse | Original filter, sigma 17 | 32.0 dB | 7.5 % | 0.09 % | 13.5 % | −0.23 |

What the numbers say:

- **The filter does what the classification needs.** On the fine phantom the
  share of wrongly classified voxels falls from 12 % to about 1 %, and nearly
  all remaining errors are edge voxels, which are mixtures anyway.
- **The filter removes the voids.** In the unfiltered scan the void voxels are
  as dark as air (contrast 0.98). After filtering only 7–25 % of the step is
  left. This confirms the plan of section 5.2: bubbles must be searched in the
  raw data.
- **At 5 µm the filter cannot separate the materials.** On the coarse phantom
  8 % of the voxels stay wrong, and the voids end up brighter than the epoxy
  around them (negative contrast), because the filter fills them with values
  from nearby glass edges.
- **Sigma 17 scores better than the true sigma of 14 here.** Sigma 17 uses the
  larger patch of the parameter table (radius 2 instead of 1). Step 0.4 looks
  at this over a range of sigma values.

### 6.4 The brute-force filter

`nlmBruteForce()` computes the filter exactly as its definition reads: for
every voxel and every shift in the search window it sums the squared
differences of the two patches voxel by voxel, turns the sum into a weight
and averages. The formulas, taken from the original source, are written out in
`include/nlm_bruteforce.hpp`. Differences from the original program are in
the arithmetic only:

| | Original program | Brute force |
|---|---|---|
| Patch distance | Through a running-sum table in 32-bit floats | Summed directly in double precision |
| Weighted sums | 32-bit floats | Double precision |
| Order of operations | Depends on thread timing (section 4.3) | Fixed per voxel; independent of the number of threads |
| Patch and search radius | One value for all axes | One value per axis (requirement R5) |
| Cost per voxel and shift | 8 look-ups | One operation per patch voxel |

It is slow by design and meant for volumes up to about 64³: a 64³ volume at
sigma 10 takes 34 s on the development machine.

`originalParams(sigma)` returns the patch radius, search radius and `h` that
the original program picks for a noise level.

Unit tests (`make test`, now 10 test cases with 772 checks): the parameter
table; a constant volume stays unchanged; a three-voxel case worked out by
hand; a search window of one voxel returns the input; results stay within the
input's range; mirroring the volume mirrors the result; swapping two axes of
volume and radii swaps the result; bad arguments are rejected.

### 6.5 Finding: the original program's 32-bit running sums cost accuracy

```
make bruteforce
python tests/reference/run_reference.py build/phantom/phantom_tiny --sigma 10 --bruteforce --keep-float
python tests/reference/run_reference.py build/phantom/phantom_tiny --sigma 10 --one-worker --keep-float
python tests/reference/compare_volumes.py build/reference/phantom_tiny_ref_s10.f32 build/reference/phantom_tiny_bf_s10.f32
```

Brute force against the one-worker original program, differences on the 0–255
scale:

| Volume | Sigma | Largest difference | RMS difference |
|---|---|---|---|
| 32³ phantom | 10 | 2.37 | 0.037 |
| 32³ phantom | 17 | 0.53 | 0.012 |
| 64³ test scan (Sample2) | 10 | 2.74 | 0.049 |

A difference of 2.7 is about 700 counts. Most voxels agree closely; a small
number are off by a lot, at any distance from the volume's edge.

To find the cause, a throwaway copy of the original source was changed in one
respect only: the running-sum table and the distance taken from it were
switched from 32-bit to 64-bit floats. That copy agrees with the brute-force
filter to 0.007 at the worst (RMS 0.0006) for both sigma values. So:

- **The brute-force filter computes the same filter as the original
  program.** The definition in the header is the right one.
- **The differences come from the 32-bit running sums.** The table adds up
  squared differences over the whole block, reaching values of 10⁸ and more,
  where a 32-bit float resolves only steps of 8 or 16. A patch distance is
  the small difference of eight such large numbers, so it is off by tens,
  while the filter strength `h²` is 16 (at sigma 10). For voxels with few
  good matches that error changes the weights noticeably.
- **The error grows with the block size.** The original program works in
  blocks of 256³, where the running sums are far larger than in these small
  volumes. Step 0.4 measures this on a larger crop.

This answers the "numerical risk (not yet tested)" of section 1.2 ahead of
step 1.6: the risk is real, and double precision in the running sums removes
it. The consequence for the acceptance rule of phase 1 ("the output matches
the reference within a stated tolerance") is that the tolerance against the
original program cannot be tight.

**Decided in review: the brute-force filter is the reference.** New code is
required to match it tightly on the small test volumes. Its difference from
the original program is measured, reported and explained by this finding.
The acceptance rule in section 1.6 is changed accordingly. Because the
brute-force filter is in git, the reference can now be rebuilt from a fresh
clone, which the original program could not.

### 6.6 Results on both platforms

| Check | Development machine (Windows, MinGW g++ 16.2) | Second machine (Linux, g++ 15.2) |
|---|---|---|
| `make test` | 10 test cases, 772 checks, all pass | the same |
| `make test SANITIZE=1` | not available | all pass, no sanitizer report |
| Compiler warnings | none | none |
| Phantom files (fine, coarse and a 32³ one) | | checksums identical to the development machine |
| Brute-force result on the 32³ phantom | | identical, bit for bit |
| Brute force, 32³ phantom, sigma 10 | 8 s | 1.2 s |

### 6.7 Open points

- **Glass-to-epoxy ratio of 2.5** is assumed (section 6.2).
- **The phantoms are not tracked.** They are regenerated from the seed, and
  the checksums in `phantom.txt` show that both platforms get the same files.
  Step 0.4 stores the reference outputs it needs in `tests/data/`.
- **The Python scripts are for phantoms only.** Generating and scoring the
  fine phantom take a few seconds each, but both hold the volume in memory, so
  they do not suit full scans.

---

## 7. Step 0.4 — Baseline

Status: **done, committed.**

This step records the numbers that later steps are compared with: how fast
the original program is, how much memory it takes, how far it is from the
brute-force reference, and how well it and the simple z filters score on the
phantoms. Three results matter beyond the bookkeeping, and they are stated
first.

### 7.1 Main results

1. **The original program's error grows with the block size.** Against the
   brute-force filter it is off by up to 2.7 (on the 0–255 scale) in a 64³
   volume, 7.5 in a 128³ volume and 17.6 in a 232³ volume, at sigma 10. The
   program filters real scans in blocks of 256³, so this is its working
   condition. At sigma 17, the value used on the real samples, the error is
   much smaller: up to 1.4, with an RMS of 0.04 (section 7.4).
2. **A plain mean along z classifies the fine phantom as well as the original
   filter.** A window of 17 slices leaves 0.98 % of the voxels wrong, the
   original filter at its best sigma 1.09 %. The original filter gives the
   lower noise (40.8 dB against 38.2 dB). The phantom favours the z filters,
   so this is a lead to follow, not a verdict (section 7.5).
3. **The original filter's quality depends sharply on sigma.** It is best at
   sigma 16–17 for a true noise sigma of 14, and it collapses when sigma is
   set to twice the true value (section 7.5).

### 7.2 What was added

| File | Purpose |
|---|---|
| `run/baseline.py` | Runs the whole baseline and writes the result tables |
| `run/baseline/windows/`, `run/baseline/linux/` | The result tables of the two machines: `timing.csv`, `scores.csv`, `accuracy.csv` |
| `tests/phantom/z_filter.py` | Median or mean filter along z, the kind of filter used on the rods earlier (section 5.5) |
| `tests/data/original_outputs/` | Results of the original program (one worker) for the two tracked test scans at sigma 10 and 17: four `.vol` + `.vgi` pairs, 0.5 MB each |
| `tests/reference/run_reference.py` | Now also reports the peak memory of the program it runs (measured on Linux only) |
| `tests/data/make_crops.py` | New options `--centre`, `--scan` and `--out`, for cutting larger local crops |

```
make reference bruteforce
python run/baseline.py                          everything
python run/baseline.py --parts timing accuracy --sizes 64
```

Larger crops for the accuracy part are local only (`data/crops/`, not
tracked). They are cut from the middle of scan `Sample2_5mu` with:

```
python tests/data/make_crops.py --scan Sample2_5mu --size 128 --centre --out data/crops
python tests/data/make_crops.py --scan Sample2_5mu --size 232 --centre --out data/crops
```

The 232³ crop is the largest that the original program filters as one block
(232 voxels plus a border of 12 on each side at sigma 17 make 256). It also serves for
looking at results by eye. The full baseline took about 2.5 hours on the
second machine, nearly all of it for the 232³ crop; on the development
machine only the timing and the 64³ accuracy were run.

### 7.3 Run time and memory

64³ test scans, best of three runs. The two scans give the same times within
a few per cent; the table shows `Sample1_18p15_5mu_c64`.

| Program | Sigma | Development machine | Second machine | Peak memory (second machine) |
|---|---|---|---|---|
| Original, all workers | 10 | 13.8 s | 7.6 s | 3 908 MB |
| Original, all workers | 17 | 12.6 s | 8.6 s | 3 908 MB |
| Original, one worker | 10 | 26.4 s | 27.1 s | 324 MB |
| Original, one worker | 17 | 31.4 s | 31.1 s | 324 MB |
| Brute force | 10 | 29.7 s | 11.0 s | 26 MB |
| Brute force | 17 | 75.1 s | 32.5 s | 26 MB |

Larger volumes, second machine, one run each:

| Volume | Sigma | Original, one worker | Brute force (all cores) |
|---|---|---|---|
| 64³ | 10 | 28 s | 11 s |
| 64³ | 17 | 32 s | 32 s |
| 128³ | 10 | 232 s | 90 s |
| 128³ | 17 | 233 s | 258 s |
| 232³ | 10 | 1 288 s | 541 s |
| 232³ | 17 | 1 232 s | 2 949 s |

What the numbers say:

- **Memory is the original program's weak point.** It takes 3.9 GB to filter
  a 0.5 MB volume, about 324 MB for every worker, whatever the volume size.
  The brute-force filter needs 26 MB for the same job. Step 1.7 (smaller
  per-thread storage) starts from here.
- **The original program's run time grows in proportion to the voxel count**,
  about 105 s per million voxels with one worker and about 29 s per million
  voxels with all workers on the second machine. It does not depend on the
  patch size. For a present scan (1142 × 1142 × 1000 voxels) this projects to
  roughly 10 hours, in line with the 14 hours the author reported, and to
  roughly 30 hours for the planned 2000 × 2000 × 975 scan. These are the
  figures the optimisation of step 1.7 and the GPU version of phase 5 are
  measured against.
- **The brute-force filter is usable as a reference up to about 128³.** Its
  run time grows with the patch size, so sigma 17 costs three times as much
  as sigma 10. The 64³ test scans take 11–32 s on the second machine and
  30–75 s on the development machine, which is acceptable for a check that is
  run once per step. The 232³ run at sigma 17 took 49 minutes.

### 7.4 Accuracy of the original program

Original program (one worker) against the brute-force filter on the same
volume. Float results are on the 0–255 scale; one unit is 257 counts.

| Volume | Sigma | Largest difference | RMS difference | Largest, in counts | RMS, in counts | Voxels that differ in the 16-bit result |
|---|---|---|---|---|---|---|
| `Sample1_18p15_5mu_c64` | 10 | 5.29 | 0.088 | 1 359 | 22.5 | 52 % |
| `Sample1_18p15_5mu_c64` | 17 | 1.08 | 0.022 | 276 | 5.6 | 64 % |
| `Sample2_5mu_c64` | 10 | 2.74 | 0.049 | 705 | 12.5 | 68 % |
| `Sample2_5mu_c64` | 17 | 0.24 | 0.003 | 63 | 0.9 | 31 % |
| `Sample2_5mu_c128` | 10 | 7.53 | 0.148 | 1 936 | 38.1 | 89 % |
| `Sample2_5mu_c128` | 17 | 0.39 | 0.008 | 100 | 2.0 | 49 % |
| `Sample2_5mu_c232` | 10 | 17.64 | 0.480 | 4 533 | 123.3 | 97 % |
| `Sample2_5mu_c232` | 17 | 1.43 | 0.036 | 367 | 9.2 | 71 % |

The 64³ rows are identical on both machines, digit for digit. The brute-force
results for the 64³ scans are identical bit for bit on both machines, and so
are the original program's one-worker results.

What the numbers say:

- **The error grows with the block size.** This fits the cause found in step
  0.3: the running sums are kept in 32-bit floats, and the larger the block,
  the larger the sums and the coarser their last digit. From 64³ to 232³ the
  RMS difference grows tenfold.
- **Small sigma suffers most.** The filter compares each patch distance with
  a threshold that grows with sigma² and with the patch size (5 400 at sigma
  10, 72 250 at sigma 17). The same rounding error in the distance is a much
  larger share of the smaller threshold.
- **At the setting used on the real samples the error is small.** Sigma 17
  with 256³ blocks gives an RMS difference of about 0.04 on the 0–255 scale,
  against a noise sigma of 17. Earlier results made at sigma 17 are therefore
  not called into question. A run at sigma 15 or below would be: there the
  largest errors reach 7 % of the full scale.
- **For the tolerance of phase 1** this confirms the decision of section 6.5.
  New code is compared with the brute-force filter; the table above is the
  stated and explained distance of the original program from it.

### 7.5 Phantom scores

Scores as defined in section 6.3. The phantoms have a true noise sigma of 14.
Results of the second machine; the original program ran with all workers.

**Fine phantom** (20-voxel fibres, as in the coming 1 µm scans):

| Filter | Setting | PSNR | `wrong_all` | `wrong_interior` | `wrong_edge` | `void_contrast` |
|---|---|---|---|---|---|---|
| None | | 25.2 dB | 12.05 % | 10.19 % | 25.6 % | 0.98 |
| Original | sigma 8 | 30.0 dB | 3.60 % | 0.53 % | 26.0 % | −0.35 |
| Original | sigma 10 | 31.5 dB | 2.77 % | 0.29 % | 20.9 % | −0.30 |
| Original | sigma 12 | 34.6 dB | 1.85 % | 0.10 % | 14.7 % | −0.15 |
| Original | sigma 14 | 38.1 dB | 1.47 % | 0.06 % | 11.8 % | 0.07 |
| Original | sigma 15 | 39.3 dB | 1.41 % | 0.05 % | 11.3 % | 0.16 |
| Original | **sigma 16** | **40.8 dB** | 1.09 % | 0.05 % | 8.7 % | 0.10 |
| Original | sigma 17 | 39.9 dB | 1.16 % | 0.05 % | 9.3 % | 0.25 |
| Original | sigma 20 | 36.2 dB | 1.51 % | 0.06 % | 12.1 % | 0.22 |
| Original | sigma 25 | 34.3 dB | 2.36 % | 0.22 % | 18.0 % | 0.01 |
| Original | sigma 30 | 26.8 dB | 8.61 % | 3.65 % | 44.8 % | −0.05 |
| Original | sigma 35 | 21.4 dB | 36.06 % | 33.83 % | 52.3 % | −0.48 |
| Median along z | radius 1 | 28.7 dB | 4.50 % | 2.72 % | 17.5 % | 0.87 |
| Median along z | radius 2 | 30.6 dB | 2.58 % | 0.98 % | 14.3 % | 0.77 |
| Median along z | radius 4 | 32.8 dB | 1.55 % | 0.21 % | 11.4 % | 0.50 |
| Median along z | radius 8 | 35.2 dB | 1.18 % | 0.07 % | 9.3 % | 0.19 |
| Median along z | radius 16 | 36.8 dB | 1.19 % | 0.06 % | 9.4 % | 0.06 |
| Median along z | radius 32 | 35.5 dB | 2.05 % | 0.18 % | 15.7 % | −0.02 |
| Mean along z | radius 1 | 29.9 dB | 3.05 % | 1.39 % | 15.2 % | 0.81 |
| Mean along z | radius 2 | 32.1 dB | 1.75 % | 0.33 % | 12.1 % | 0.69 |
| Mean along z | radius 4 | 34.6 dB | 1.20 % | 0.07 % | 9.5 % | 0.48 |
| Mean along z | **radius 8** | 37.0 dB | **0.98 %** | 0.06 % | 7.7 % | 0.26 |
| Mean along z | radius 16 | 38.2 dB | 1.04 % | 0.06 % | 8.2 % | 0.11 |
| Mean along z | radius 32 | 35.6 dB | 1.99 % | 0.17 % | 15.3 % | 0.02 |

The original filter took 45 s per run on this phantom (204 s at sigma 35,
which uses the larger search window), the z filters 1–7 s as Python scripts.

**Coarse phantom** (4-voxel fibres, as in the present 5 µm scans), the best
setting of each filter and the settings next to it:

| Filter | Setting | PSNR | `wrong_all` | `wrong_interior` | `wrong_edge` | `void_contrast` |
|---|---|---|---|---|---|---|
| None | | 25.2 dB | 18.1 % | 8.22 % | 26.2 % | 1.04 |
| Original | sigma 14 | 31.8 dB | 8.5 % | 0.17 % | 15.2 % | −0.57 |
| Original | **sigma 16** | 33.1 dB | 6.9 % | 0.07 % | 12.5 % | −0.31 |
| Original | sigma 17 | 32.0 dB | 7.5 % | 0.09 % | 13.5 % | −0.23 |
| Original | sigma 25 | 24.0 dB | 36.8 % | 17.6 % | 52.4 % | −0.67 |
| Median along z | radius 4 | 32.4 dB | 6.7 % | 0.16 % | 12.0 % | −0.02 |
| Median along z | radius 8 | 33.3 dB | 6.4 % | 0.06 % | 11.5 % | −0.21 |
| Mean along z | radius 4 | 34.2 dB | 5.4 % | 0.05 % | 9.8 % | 0.13 |
| Mean along z | **radius 8** | **34.4 dB** | **5.4 %** | 0.03 % | 9.8 % | −0.10 |
| Mean along z | radius 16 | 30.4 dB | 10.5 % | 0.75 % | 18.4 % | −0.26 |

All rows are in `run/baseline/linux/scores.csv`.

What the numbers say:

- **Sigma is not a safe dial on the original filter.** The best result comes
  at sigma 16, the first value that uses the larger patch (radius 2), not at
  the true sigma of 14. Above sigma 20 the result gets worse quickly, and at
  about twice the true sigma the filter treats all patches as alike and
  becomes a plain average over the whole search window: on the fine phantom
  36 % of the voxels are then wrong, three times as many as without any
  filter. The value 17 used on the real samples is close to the best one,
  provided the noise there is near 14–15 on the 0–255 scale. Step 1.4 (noise
  estimate) therefore matters more than it looked, and step 1.5 should
  separate the choice of patch size from sigma.
- **Averaging along the fibres works.** With a window of 9–33 slices the mean
  along z classifies as well as the original filter or slightly better, on
  both phantoms, at a tiny fraction of the cost. The mean does better than
  the median throughout. This supports requirement R5 (separate radii for x,
  y and z): a Non-Local Means filter that searches further along z than
  across should combine the strengths of both.
- **The window along z has an upper limit set by the twist.** Radius 32 on
  the fine phantom and radius 16 on the coarse one are clearly worse, because
  over that length the fibres have moved sideways by a good part of their
  diameter. The median with radius 16 that was used earlier at 5 µm is, by
  this measure, already past the best window (radius 8).
- **All good settings remove the voids.** Every filter that gets
  `wrong_all` near its minimum leaves less than 30 % of the void contrast.
  Short z windows keep the voids (radius 1–2: 70–87 %) but classify poorly.
  This confirms again that bubbles are searched in the raw data (section
  5.2).
- **At 5 µm no filter separates the materials well.** The best result on the
  coarse phantom leaves 5.4 % of the voxels wrong, nearly all of them edge
  voxels, which at four voxels per fibre are a large share of the volume.

How far these conclusions carry: the phantom favours the z filters. Its
noise is independent from voxel to voxel, while CT noise is correlated
between neighbours; its fibres are perfectly regular with one common twist;
and the fine phantom is only 96 slices deep, so the widest windows are cut
short at both ends. The comparison has to be repeated on real data when the
classification of phase 4 gives a way to judge it.

### 7.6 Stored outputs of the original program

`tests/data/original_outputs/` holds the one-worker results of the original
program for the two tracked test scans, at sigma 10 and sigma 17. They are
identical on both machines (same checksums). Since `ORIGINAL_CODE/` is not in
git, these files are what a fresh clone has of the original program's
behaviour. They are for stating how far new code is from the original, with
`tests/reference/compare_volumes.py`; the pass/fail reference is the
brute-force filter (section 6.5).

### 7.7 Results on both platforms

| Check | Development machine | Second machine |
|---|---|---|
| `make test` | 10 test cases, 772 checks, all pass | the same |
| Timing, 64³ test scans | measured | measured, with peak memory |
| Accuracy, 64³ test scans | measured | identical values |
| Accuracy, 128³ and 232³ crops | not run | measured |
| Phantom scores | not run | measured |
| Original program's one-worker outputs | stored in `tests/data/original_outputs/` | identical checksums |

### 7.8 Open points

- **Phase 0 is complete with this step.** Decided in review: no separate
  report for phase 0; its results go into the phase 1 report.
- **Per-axis radii come late in the plan.** Result 2 suggests that searching
  further along z is the most promising change to the filter itself. It is
  step 1.5 in the plan. It could be tried earlier, with the brute-force
  filter, which already takes separate radii for x, y and z, on the phantoms
  only. Proposal: do this as a short experiment at the start of step 1.5,
  not before, so that the reader and writer of steps 1.1–1.2 come first as
  planned.
- **The baseline is not part of `make test`.** It takes hours and is run by
  hand when a step calls for it.
- **The z filters are Python scripts for phantoms and small crops.** If a
  mean along z is to be used on full scans, it is added to the C++ program
  in phase 1.

---

## 8. Step 1.1 — Scan reader

Status: **done, committed.**

The program can now read a scan as the CT lab delivers it, without ImageJ
and without a converted copy. `denoiser info` is the first working command.

### 8.1 What was added

| File | Purpose |
|---|---|
| `include/scan_reader.hpp`, `src/scan_reader.cpp` | The scan reader |
| `include/volume.hpp` | Volume size and region, shared by all parts (the size was in `nlm_bruteforce.hpp` before) |
| `src/main.cpp` | The command `denoiser info` |
| `tests/test_scan_reader.cpp` | 10 test cases for the reader |

The program version is now 1.1.0, after the step it implements.

### 8.2 The command

```
denoiser info data/RAW/Sample2_5mu
```

```
Voxel file:   data/RAW/Sample2_5mu/Sample2_5mu_.vol
Size:         1142 x 1142 x 1000 voxels (1304164000)
Sample type:  16-bit unsigned integer, little-endian
Header:       0 bytes
Voxel size:   0.005 x 0.005 x 0.005 mm
File size:    2608328000 bytes, as expected
Minimum:      0
Maximum:      65535
Mean:         28235.24
Std. dev.:    4213.6733
```

This takes 8 seconds for the 2.6 GB scan on the development machine. The
scan may be given as its folder, its `.vgi` file or its `.vol` file.

A raw file that has no `.vgi` file is described by hand, for example one of
the float files of the old workflow:

```
denoiser info SampleB_BE255_250x250x250.raw --raw 250 250 250 --type f32 --big-endian
```

| Option | Meaning |
|---|---|
| `--raw W H D` | Volume size in voxels; marks the file as described by hand |
| `--type T` | `u8`, `u16` (default) or `f32` |
| `--big-endian` | Samples are stored highest byte first (default: lowest first) |
| `--header N` | Bytes to skip before the first voxel (default 0) |

Exit codes: 0 success, 1 wrong command line, 3 the scan cannot be read. The
message names the file and the reason, for example:

```
denoiser info: scan.vol: the file has 524288 bytes, but 64 x 64 x 63 voxels
of type 16-bit unsigned integer with a header of 0 bytes need 516096
```

### 8.3 The reader

All functions are documented in `include/scan_reader.hpp`.

| Function | What it does |
|---|---|
| `openScan(path)` | Finds the `.vgi` file, reads the description and checks the voxel file's size |
| `parseVgi(text, name)` | Reads the description from the text of a `.vgi` file |
| `describeRaw(path, size, type, byte order, header)` | The same description given by hand |
| `checkScan(info)` | Checks that the file has exactly header + voxels × bytes per voxel |
| `readRegion(info, region)` | Reads a box of voxels; only those voxels are read from the file |
| `readVolume(info)` | Reads all voxels; for small volumes |
| `scanStatistics(info)` | Minimum, maximum, mean and standard deviation, read slice by slice |

How it behaves:

- **What is taken from the `.vgi` file:** `size`, `bitsperelement`,
  `datatype`, `resolution` and `unit`, and `skipheader` if present. Keys are
  matched without regard to case, and the first occurrence counts, as in the
  Python scripts of phase 0.
- **The voxel file is found by name, not from the `.vgi` text.** It is the
  `.vgi` file's name with the extension `.vol`. The lab's `.vgi` files hold
  a path of the lab's own computer, which is of no use here.
- **Byte order.** The `.vgi` format does not state it; little-endian is
  assumed, as verified in section 2.3. Files described by hand can be
  big-endian.
- **Sample types:** 8-bit and 16-bit unsigned integers and 32-bit floats.
  Anything else is refused with a message that lists these.
- **Values are returned unchanged,** as 32-bit floats: a count of 31139
  becomes 31139.0. Mapping to the working range 0–255 is step 1.3.
- **Regions are zero-based** and must lie inside the volume; a region that
  does not fit is an error, not silently cut.
- **Large files work.** Positions are 64-bit throughout; reading single
  voxels at the end of the 2.6 GB scan gives the same values as an
  independent Python check.
- **Errors** are reported by exception with the file name in the message;
  the command line turns them into the messages above.

### 8.4 Tests

| Test case | Checks |
|---|---|
| `.vgi` text of a lab scan | Windows line ends, indented keys, a second `Size` further down, a lab path in `Name` |
| Sample types and optional keys | 8-bit, float, header, missing voxel size |
| Unreadable `.vgi` text | 11 kinds of faulty text are refused; the message names the file |
| Opening a scan | From the folder, the `.vgi` file and the `.vol` file |
| Missing and mismatching files | No folder, no `.vgi`, two `.vgi` files, no `.vol`, a `.vol` one byte short |
| Region content | Every voxel of a region, for 16-bit in both byte orders (one with a header), 8-bit, and float in both byte orders |
| Region limits | Full-width regions, the last voxel alone, five regions that do not fit |
| Description by hand | Wrong size, wrong header, missing file |
| Statistics | Known values; NaN and infinity counted and left out |
| Tracked test scans | Sum, minimum, maximum and three voxels of each scan against values computed independently in Python |

### 8.5 Results on both platforms

| Check | Development machine | Second machine |
|---|---|---|
| `make test` | 20 test cases, 902 checks, all pass | the same |
| `make test SANITIZE=1` | not available | all pass, no sanitizer report |
| Compiler warnings | none | none |
| `denoiser info` on a tracked test scan | | same figures as on the development machine |
| `denoiser info` on the full scans (2.6 GB) | 8 s each | not run (the full scans are on the development machine only) |

Statistics of the full scans, in counts:

| Scan | Minimum | Maximum | Mean | Standard deviation |
|---|---|---|---|---|
| `Sample1_18p15_5mu` | 0 | 65535 | 33 403 | 5 041 |
| `Sample2_5mu` | 0 | 65535 | 28 235 | 4 214 |

### 8.6 Open points

- **Both full scans reach 0 and 65535.** Some voxels are therefore at the
  limits of the 16-bit range, where the true value may have been cut off in
  the reconstruction. The 64³ test scans do not contain such voxels. How
  many there are and where they lie (inside the rod, or only in the air and
  the holder around it) is counted in step 1.3, because it bears on
  requirement R2 (nothing is clipped).
- **The reader returns floats.** A block of voxels takes four bytes per
  voxel in memory although the file has two. This is what the filter works
  with anyway; step 1.7 decides whether blocks are converted on the fly.
- **Step 1.2 (writer) is next.** It writes `.vol` plus `.vgi` in the input's
  sample type and uses the description of this step.

---

## 9. Step 1.2 — Scan writer

Status: **done, committed.**

The program can now write a volume in the CT lab's format, `.vol` plus
`.vgi`, in the same sample type as the input. A full scan read with the
reader of step 1.1 and written back block by block is identical to the
original, byte for byte.

### 9.1 What was added

| File | Purpose |
|---|---|
| `include/scan_writer.hpp`, `src/scan_writer.cpp` | The scan writer |
| `tests/test_scan_writer.cpp` | 8 test cases for the writer |

The program version is now 1.2.0. There is no new command: the writer is
used by `denoiser run` (step 1.8) and, before that, by the tests and by the
pre-processing of step 1.3.

### 9.2 The writer

All functions are documented in `include/scan_writer.hpp`.

| Function | What it does |
|---|---|
| `createScan(info)` | Writes the `.vgi` file and creates the `.vol` file at its full size, all voxels zero |
| `writeRegion(info, region, values)` | Writes a box of voxels into that file; everything outside the box stays as it is |
| `writeVolume(info, values)` | Both of the above for a volume held in memory |
| `vgiText(info)` | The text of the `.vgi` file |

The description `info` is the same one the reader returns. Writing a result
in the input's format is therefore: copy the input's description, change the
file name, write.

How it behaves:

- **Block by block.** A large result is created once and then filled region
  by region, in any order. This is what the filter needs, since it works on
  blocks and never holds a full scan in memory.
- **Same sample type as the input:** 8-bit or 16-bit unsigned integers or
  32-bit floats, in either byte order.
- **Rounding.** Values arrive as floats. For integer types they are rounded
  to the nearest whole number.
- **Nothing is cut off unnoticed.** A value below 0 or above the largest
  storable number (65535 for 16 bits) is stored as that limit, and NaN as 0.
  Every write returns how many values were below, above or NaN, so the
  caller can report it. With float output all values are stored unchanged.
- **The `.vgi` file** holds the size, the sample type, the voxel file's name
  without any folder, and the voxel size. For 16-bit scans it is, character
  for character, the format of the tracked test scans. No path is written,
  so nothing about the computer that made the file gets into it.
- **Guard against writing into the wrong file.** Before writing a region
  the file's size is checked against the description; a file that does not
  match is refused.
- **Files are written without a header,** as the lab's files are.

### 9.3 Tests

| Test case | Checks |
|---|---|
| `.vgi` text | Exact text for a 16-bit scan; float; no folder in the text; the reader of step 1.1 understands it |
| Creating a scan | Both files exist, all voxels zero, the reader opens it, creating again replaces it |
| Write and read back | 16-bit in both byte orders (bytes on disk checked), 8-bit, float in both byte orders |
| Values that do not fit | Rounding at the halves, values below and above the range, NaN, infinity: stored values and counts, for 16-bit, 8-bit and float |
| A region | The box holds the values, every voxel outside it is untouched |
| Blocks against the whole | A volume written in 12 blocks, last block first, and one written slice by slice are byte-identical to the same volume written at once |
| Refusals | Regions that do not fit, wrong number of values, a file of another size, a missing file, a header, an empty name, a missing folder; the file is left as it was |
| Tracked test scans | Read and written back: `.vol` and `.vgi` identical to the originals, byte for byte |

### 9.4 Results on both platforms

| Check | Development machine | Second machine |
|---|---|---|
| `make test` | 28 test cases, 1172 checks, all pass | the same |
| `make test SANITIZE=1` | not available | all pass, no sanitizer report |
| Compiler warnings | none | none |
| Full scan `Sample2_5mu` (2.6 GB) copied in slabs of 50 slices, last slab first | identical to the original, byte for byte; 18 s for reading and writing | not run (the full scans are on the development machine only) |

The full-scan copy also shows that file positions beyond 2 GB are
handled correctly when writing.

### 9.5 Open points

- **Rounding of exact halves.** A value exactly halfway between two counts
  is rounded away from zero (2.5 becomes 3). The Python scripts of phase 0
  do the same for the values that occur there. It matters only for
  bit-for-bit comparisons, and those compare programs that use this writer.
- **What the caller does with the counts** of limited values is decided in
  step 1.3 (mapping back from the working range) and step 1.8 (run record).
  Proposal: the run record states them, and the program prints a warning
  when the count is not zero.
- **Step 1.3 (pre-processing) is next:** region of interest, mapping to the
  working range and back, NaN handling, the count of voxels at the limits of
  the 16-bit range (section 8.6), and the measurement of brightness drift
  from the rod centre to the rim.

---

## 10. Step 1.3 — Pre-processing and brightness drift

Status: **done, committed.**

This step adds the pieces between reading a scan and filtering it: cutting
a region, mapping counts to the working range and back, and handling NaN
values. It also measures how brightness changes from the centre of the scan
to its rim (requirement R4). The measurement gave three results that matter
for the later phases, stated first.

### 10.1 Main results

1. **Both scans get about 10 % brighter from the centre to the rim.** The
   rise is about 3 000 counts (12 on the 0–255 scale). That is more than a
   quarter of the whole spread between dark and bright voxels, so one fixed
   threshold between glass and epoxy cannot hold across a slice. The drift
   is too slow to disturb the filter (section 10.5).
2. **The present scans contain no air.** Neither a bore hole nor air around
   the rod shows in the profiles. They therefore cannot be used to test the
   bore-hole calibration of phase 4.
3. **The usable part of a scan is a circle, not the full square.** Beyond a
   radius of about 540 voxels the brightness falls steeply; the corners of
   the square lie outside the reconstructed field of view.

The voxels at 0 and 65535 found in step 1.1 turned out to be harmless: each
scan has exactly one of each, next to the rotation axis.

### 10.2 What was added

| File | Purpose |
|---|---|
| `include/preprocess.hpp`, `src/preprocess.cpp` | Intensity mapping, centred cube, NaN replacement, cutting a region |
| `include/radial_profile.hpp`, `src/radial_profile.cpp` | Brightness against distance from the axis |
| `include/scan_reader.hpp` | The statistics now count the voxels at 0 and at the top of the range |
| `src/main.cpp` | Commands `denoiser crop` and `denoiser profile`; `denoiser info` prints the counts at the limits |
| `tests/test_preprocess.cpp` | 10 test cases |
| `run/profile/` | The profiles of the two full scans (CSV) |

The program version is now 1.3.0.

### 10.3 Commands

**Cut a region** (replaces `PickCube.ijm` and the Python crop script for
everyday use):

```
denoiser crop data/RAW/Sample2_5mu --centre-cube 232 --out data/crops/Sample2_5mu_c232
denoiser crop data/RAW/Sample2_5mu --roi 446 446 374 250 250 250 --out data/crops/old_cube
```

| Option | Meaning |
|---|---|
| `--roi X0 Y0 Z0 W H D` | First voxel, counted from 0, and size |
| `--centre-cube N` | Cube of N voxels in the middle of the scan |
| `--out DIR` | Output folder; created if missing |

The result is a scan folder in the lab's format: `DIR/<name>.vol` and
`.vgi`, where `<name>` is the last part of `DIR`, in the source's sample
type and with its voxel size. `DIR/crop_record.txt` states the source and
the first voxel, so that a result computed on the crop can be placed back
into the grid of the full scan (requirement R3):

```
denoiser 1.3.0 crop
source = data/RAW/Sample2_5mu/Sample2_5mu_.vol
source size = 1142 1142 1000
first voxel (x y z, counted from 0) = 455 455 384
size (x y z) = 232 232 232
```

**Brightness profile:**

```
denoiser profile data/RAW/Sample2_5mu --out run/profile/Sample2_5mu_radial.csv
```

| Option | Default | Meaning |
|---|---|---|
| `--centre X Y` | middle of the x-y plane | The axis, in voxels |
| `--ring W` | 8 | Ring width in voxels |
| `--slices Z0 N` | all | Use N slices starting at slice Z0 |
| `--out FILE` | screen | Where the table goes |

The voxels are sorted into rings by their distance from the axis. The table
has one line per ring: inner and outer radius, number of voxels, mean, 10th
percentile, median, 90th percentile, and the number of voxels at the bottom
and at the top of the value range. In a rod the 10th percentile follows the
darker material and the 90th the brighter one, so drift can be told from a
change in the mix of materials. A full scan takes 8 to 20 seconds.

**`denoiser info`** now ends with two more lines:

```
At 0:         1 voxels (7.67e-08 %)
At 65535:     1 voxels (7.67e-08 %)
```

### 10.4 Functions

All are documented in the headers.

| Function | What it does |
|---|---|
| `IntensityMap` with `toWork()`, `toStored()` | Linear mapping between stored values and the working range: `low` becomes 0, `high` becomes 255 |
| `defaultMap(type)` | Full range of the sample type (0–65535 for 16 bits), the scaling of the old workflow |
| `mapToWork()`, `mapToStored()` | The mapping applied to a block of values |
| `centredCube(size, n)` | Region of a cube in the middle of a volume |
| `replaceNonFinite(values, replacement)` | Replaces NaN and infinity, returns how many |
| `cropScan(source, region, path)` | Writes a region as a scan of its own, slab by slab |
| `radialProfile(info, settings)` | The brightness profile |

How the mapping behaves:

- **Linear, and nothing is clipped** (requirement R2). A value outside the
  mapped range lands outside 0–255 and comes back unchanged. Only the writer
  limits values to what the file type can hold, and it counts them (section
  9.2).
- **Same numbers as the old workflow.** With the default mapping the working
  value is count × 255 / 65535, computed in double precision and rounded to
  float once, exactly as the reference scripts of phase 0 do.
- **Loss-free there and back.** Every one of the 65 536 possible counts
  returns to itself after mapping to the working range and back, also with
  a narrower or wider mapped range.
- **The mapping decides what sigma means.** The filter picks its settings
  from sigma on the working scale. With a narrower mapped range the same
  noise is a larger sigma. The default keeps old and new runs comparable.

NaN and infinity occur only in float files. The filter cannot work with
them. The function to replace them is in place; what `denoiser run` does by
default is proposed in section 10.7.

### 10.5 Brightness drift of the two scans

Profiles around the middle of the x-y plane, all 1000 slices, rings of 8
voxels, values in counts. The full tables are in `run/profile/`.

| | `Sample1_18p15_5mu` | `Sample2_5mu` |
|---|---|---|
| Centre (radius below 96 voxels): mean | 33 022 | 27 759 |
| Centre: 10th / 90th percentile | 27 148 / 38 855 | 22 787 / 32 692 |
| Brightest ring: radius | 536–544 voxels (2.7 mm) | 536–544 voxels |
| Brightest ring: mean | 36 329 | 30 672 |
| Brightest ring: 10th / 90th percentile | 30 248 / 42 365 | 25 590 / 35 720 |
| Rise of the mean, centre to brightest ring | 3 307 counts (10.0 %) | 2 913 counts (10.5 %) |
| The same on the 0–255 scale | 12.9 | 11.3 |
| Rise as a share of the 10th–90th percentile spread | 28 % | 29 % |
| Mean at radius 571 (edge of the inscribed circle) | 33 000 | 27 900 |
| Mean in the corners (radius 800) | 26 000 | 21 900 |

Mean brightness against radius, `Sample1_18p15_5mu`:

| Radius (voxels) | 0–100 | 144 | 240 | 336 | 432 | 480 | 536 | 552 | 568 | 600 | 700 | 800 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Mean (counts) | 33 100 | 33 200 | 33 350 | 33 930 | 34 730 | 35 350 | 36 330 | 35 110 | 33 000 | 30 780 | 27 670 | 25 970 |

The drift also depends on the height in the scan. Mean in the ring of radius
480–540 voxels:

| Slices | `Sample1_18p15_5mu` | `Sample2_5mu` |
|---|---|---|
| 0–99 | 34 869 | 29 421 |
| 450–549 | 35 904 | 30 297 |
| 900–999 | 34 827 | 29 441 |

At the centre the three heights agree within 250 counts.

What the numbers say:

- **The drift is an offset, not a change of contrast.** The 10th and the
  90th percentile rise by nearly the same amount as the mean (about 3 100
  and 3 500 counts in the first scan). Dark and bright voxels are lifted
  together.
- **It is nearly flat near the centre and grows towards the rim.** Inside a
  radius of 250 voxels the mean changes by 1.5 % in the first scan and 3 %
  in the second. Most of the rise comes in the outer third.
- **It is stronger at mid-height than at the top and bottom,** by about
  1 000 counts at the rim. A correction that depends on the radius alone
  would leave this part.
- **The filter is not disturbed by it.** The steepest part of the drift is
  about 35 counts per voxel. Over the filter's search radius of 10 voxels
  that is 350 counts, or 1.4 on the 0–255 scale, against a noise sigma of
  about 17. No correction is needed before filtering, and none is applied,
  so the values stay linear (R2).
- **The classification is disturbed by it.** A threshold that separates
  glass from epoxy at the centre is off by more than a quarter of the
  dark-to-bright spread at the rim. Phase 4 (steps 4.1–4.2) has to work with
  levels that follow the position, or correct the drift first. The profile
  command gives the data for either.
- **No air in these scans.** With the lab's scaling, air would lie near
  17 900 counts in the first scan and 15 400 in the second (section 6.2).
  The 10th percentile never comes near these values at any radius inside the
  field of view; it stays at the level of the rod material. There is no bore
  hole on the axis and no air at the rim. The first scan's field of view
  (5.7 mm) lies inside the 18.15 mm rod, which explains the rim.
- **Beyond a radius of about 540 voxels the data are not reliable.** The
  brightness turns and falls by 10 000 counts towards the corners. The
  square is 1142 voxels wide, so the inscribed circle has a radius of 571;
  the last 30 voxels inside it are already affected. Results from this zone
  should not be used.
- **The voxels at 0 and 65535** are one each per scan, both within 16
  voxels of the axis. They look like marks or a reconstruction artefact on
  the rotation axis, not like clipping. No other voxel reaches a limit,
  so the default mapping loses nothing.

### 10.6 Tests and results on both platforms

| Test case | Checks |
|---|---|
| Default mapping | Same numbers as the old workflow; float scans unchanged |
| There and back | All 65 536 counts return to themselves, for three mapped ranges |
| Linear, not clipped | Values outside the mapped range continue on the same line |
| Bad ranges | Refused; the values are left untouched |
| Centred cube | The cube of the old workflow, odd leftover space, cubes that do not fit |
| NaN replacement | NaN and both infinities replaced and counted |
| Cutting a region | Right voxels, description kept, whole scan copied byte for byte, bad region and overwriting the source refused |
| Profile on rings | A volume whose value depends on the radius only: every ring's mean exact, percentiles within one histogram class, voxel counts, part of the slices, another axis, settings out of range |
| Percentiles and limits | Known 10th, 50th and 90th percentile; voxels at the limits counted |
| Statistics at the limits | Counted for integer scans, zero for the tracked scans and for float scans |

| Check | Development machine | Second machine |
|---|---|---|
| `make test` | 38 test cases, 1489 checks, all pass | the same |
| `make test SANITIZE=1` | not available | all pass, no sanitizer report |
| Compiler warnings | none | none |
| Profile of the 232³ crop | | identical output, character for character |
| `denoiser crop --centre-cube 232` of a full scan | identical to the crop made with the Python script in step 0.4, `.vol` and `.vgi` | |
| Profiles of the full scans | 8–20 s each | not run (the full scans are on the development machine only) |

### 10.7 Open points

- **Centre of a cube and the old workflow.** `--centre-cube 250` starts at
  slice 375, the true middle. The old macro started at 374 because ImageJ
  counts slices from 1 (section 2.3). To repeat an old crop exactly, use
  `--roi 446 446 374 250 250 250`.
- **NaN in float input: proposal.** `denoiser run` refuses a float scan
  that contains NaN or infinity and says how many there are, unless a
  replacement value is given with an option. Scans from the lab are 16-bit
  and cannot contain them, so this only concerns other raw files.
- **Mapped range: proposal.** `denoiser run` keeps `--scale 0 65535` as the
  default (section 2.4), for comparability with the old runs and because
  the noise estimate of step 1.4 then gives sigma on the familiar scale.
- **Drift correction is left to phase 4,** as planned in R4. This step
  shows that it is needed for the classification and that it must follow
  both radius and height. Nothing is corrected in the filter.
- **Usable region of a scan.** Results are meaningful inside a radius of
  about 530 voxels only, about 68 % of the square slice. Proposal: `denoiser
  run` gets no mask for this (the filter does no harm outside), but phase 4
  restricts its analysis to the reliable circle, and the report states the
  radius.
- **Consequence for the planned 1 µm scan (section 5.6).** If the new scan
  is also reconstructed as a square with a circular field of view, the
  recommended 2000 × 2000 voxels give a usable circle of about 1 900 voxels
  diameter, or 1.9 mm around the 0.9 mm bore hole. The rise towards the rim
  will be present there too, because that scan also lies inside the rod.
  The bore hole on the axis gives the air level where the drift is flattest,
  which is favourable. It is worth asking the CT lab whether their
  reconstruction offers a correction for scans inside the sample.
- **The test scans and crops lie near the axis** (within 165 voxels), where
  the drift is below 2 %. The baseline of step 0.4 is not affected.
- **Step 1.4 (noise estimate) is next.**

---

## 11. Step 1.4 — Noise estimate

Status: **done, committed.**

The filter needs the noise sigma: it decides the patch size, the search size
and the smoothing strength. This step adds the command `denoiser noise`. It
reproduces the estimate of the earlier script `estimate_noise.py` and checks
that estimate on the phantom and on the full scans. The check changed the
picture, so the results come first.

### 11.1 Main results

1. **The earlier method is reproduced exactly.** The new command gives the
   same numbers as the formulas of `estimate_noise.py`, to all printed
   digits (section 11.6).
2. **On the present scans the earlier method reads about 20 % low.** It
   gives 13.2 and 11.0 on the 0–255 scale. The noise is about 16.4 in the
   first scan and 13.7 in the second.
3. **The reason: the noise of neighbouring voxels is not independent.**
   Directly neighbouring voxels share 45 % of their noise, voxels two apart
   9 %, voxels three apart 1 %. The figures are the same in both scans, so
   this comes from the reconstruction, not from the samples. Every estimator
   that compares direct neighbours therefore sees too little noise.
4. **On the phantom the earlier method reads high.** There the noise is
   independent, and the method takes fibre edges for noise: 14.9 on the fine
   phantom and 19.1 on the coarse one, against a true 14.0. On the coarse
   phantom *without any noise* it still reports 14.5.
5. **Differences along z, four voxels apart, give the right value in both
   cases.** The fibres run along z, so the material hardly changes between
   slices, and at four voxels the noise is independent again. This estimate
   is 14.1 and 14.6 on the phantoms (true 14.0) and is the 16.4 and 13.7
   quoted above for the scans.
6. **The sigma of 17 used on the earlier samples was a good value for the
   first scan** (4 % above its noise) **and high for the second** (24 %
   above).

### 11.2 What was added

| File | Purpose |
|---|---|
| `include/noise_estimate.hpp`, `src/noise_estimate.cpp` | The estimators, the summary over slices, the table by distance |
| `src/main.cpp` | Command `denoiser noise`; the region options are now shared by `crop` and `noise` |
| `tests/test_noise_estimate.cpp` | 9 test cases |

The program version is now 1.4.0.

### 11.3 Command

```
denoiser noise data/RAW/Sample1_18p15_5mu --roi 196 196 0 750 750 1000
```

```
Voxel file:   data/RAW/Sample1_18p15_5mu/Sample1_18p15_5mu_.vol
Region:       first voxel 196 196 0, size 750 x 750 x 1000
Slices:       20, from 0 to 999
Scale:        0 ... 65535 stored -> 0 ... 255

Estimator      median     mean       Q1       Q3   (working scale)
std_dev        17.976   18.045   17.922   18.045
mad            17.953   18.020   17.905   18.036
diff           13.533   13.605   13.499   13.643
resid          12.781   12.852   12.751   12.901
zdiff          12.160   12.223   12.134   12.217
zdiff_robust   12.160   12.217   12.132   12.215

Sigma from differences between voxels a distance apart (working scale)
distance        x        y        z
       1   13.523   13.504   12.144
       2   17.367   17.367   15.629
       3   17.822   17.851   16.268
       4   17.844   17.859   16.384
       5   17.812   17.834   16.410
       6   17.824   17.818   16.427

Sigma, earlier method (diff and resid):   13.157   (3381.302 stored units)
Sigma, z neighbours (zdiff_robust):       12.160   (3125.151 stored units)
Sigma, z differences 4 apart:             16.384   (4210.724 stored units)
```

| Option | Default | Meaning |
|---|---|---|
| `--roi X0 Y0 Z0 W H D` or `--centre-cube N` | whole scan | Region to examine |
| `--slices N` | 20 | Number of slices examined, spread evenly from the first to the last slice of the region |
| `--scale LO HI` | `0 65535` (`0 255` for 8-bit and float scans) | Stored values that map to 0 and 255. Sigma is printed on the 0–255 scale and in stored units |
| `--lags N` | 6 | Largest distance in the second table; 0 leaves the table out |

The example takes 5 seconds; the whole scan with its corners takes 18.
A scan with NaN or infinite values is refused.

### 11.4 The estimators

Each examined slice gives one value per estimator. The first table shows
the median, mean and quartiles of these values over the slices.

| Name | Formula (s = slice, n = next slice) | What it contains besides noise |
|---|---|---|
| `std_dev` | standard deviation of s | All material contrast. An upper limit only |
| `mad` | 1.4826 × MAD(s) | The same |
| `diff` | std of all x and y neighbour differences / √2 | Material edges within the slice |
| `resid` | 1.4826 × MAD(s − mean of the 3 × 3 neighbourhood) | Material edges within the slice. Reads 5.7 % low even on pure independent noise, because the mean contains the voxel itself |
| `zdiff` | std of (n − s) / √2 | Only what changes from slice to slice |
| `zdiff_robust` | 1.4826 × MAD(n − s) / √2 | The same, less sensitive to single outliers |

MAD(v) is the median of |v − median(v)|. The first five are the estimators
of the script. Its recommendation, printed as "earlier method", is the mean
of the medians of `diff` and `resid`. `zdiff_robust` is new.

The second table repeats the robust difference estimate for voxels 1, 2, …
N apart, along each axis. How to read it:

- **All entries equal:** independent noise, no material change over these
  distances. Every entry is sigma.
- **Entries rise over the first few distances and then stay level:** the
  noise of nearby voxels is correlated. The level is sigma. The entry at
  distance d is low by the factor √(1 − correlation at distance d).
- **Entries keep rising:** the material changes over these distances.

The last line of the output, "z differences 4 apart", is the z entry of
this table at distance 4.

### 11.5 Functions

| Function | Purpose |
|---|---|
| `sliceNoise(slice, next, width, height)` | The six estimates of one slice held in memory |
| `summarise(values)` | Median, mean and quartiles; quartiles interpolated as numpy does |
| `sampleSlices(depth, count)` | The slice numbers: i × (depth − 1) / (count − 1), rounded down |
| `estimateNoise(info, region, slices)` | Reads the slices of a scan and returns all estimates and the earlier recommendation |
| `differenceByLag(info, region, slices, max_lag)` | The table by distance |

`sliceNoise` works on plain arrays and reads no files, so it can later be
replaced by a GPU version (phase 5). All results are in stored units; the
command converts them to the working scale. Memory use is two slices.

### 11.6 Check against the earlier script

The script needs numpy, which neither machine has, so it could not be run
itself. Two independent checks stand in for it:

- **Its formulas written again in plain Python** (standard library only)
  and run on the tracked scan `Sample1_18p15_5mu_c64`, in counts:

  | Estimator | Plain Python: median / mean / Q1 / Q3 | `denoiser noise` |
  |---|---|---|
  | `diff` | 3692.052 / 3696.376 / 3581.232 / 3777.587 | identical |
  | `resid` | 3456.064 / 3471.520 / 3349.729 / 3553.442 | identical |
  | `zdiff` | 3721.410 / 3701.614 / 3450.459 / 3913.973 | identical |
  | recommended | 3574.058 | identical |

- **A slice of four voxels computed by hand** is a unit test.

The rough neighbour-difference figures of section 4.6 (12 and 14 for the
two test scans) are the `diff` values found here (12.1 and 14.4).

One difference remains possible: the script picks its slices with floating
point arithmetic, the new code with whole numbers. They can differ by one
slice when i × (depth − 1) / (count − 1) is a whole number. For 20 slices
out of 1000, 250 or 64 it never is, except for the first and last slice.

### 11.7 Check on the phantom

The phantoms of chapter 6 have independent noise of sigma 14.0 on the 0–255
scale.

| Phantom | Earlier method | z neighbours | z, 4 apart |
|---|---|---|---|
| Fine (fibres 20 voxels wide) | 14.90 | 14.05 | 14.12 |
| Coarse (fibres 4 voxels wide) | 19.14 | 14.04 | 14.65 |
| Coarse, noise-free volume | 14.47 | 0.35 | 2.02 |

- **The earlier method takes fibre edges for noise.** The error is 6 % on
  the fine phantom and 37 % on the coarse one, where every other voxel is an
  edge. A noise-free volume is reported as sigma 14.5.
- **The z estimates are right** to 0.3 % at distance 1. At distance 4 they
  are 1 % and 5 % high: the phantom fibres twist around the axis, so the
  material does change a little over four slices. The coarse phantom twists
  five times more per voxel than the fine one.

### 11.8 The noise of the two full scans

Region: the square inside the reliable circle of chapter 10, all 1000
slices (`--roi 196 196 0 750 750 1000`). Values on the 0–255 scale.

| | `Sample1_18p15_5mu` | `Sample2_5mu` |
|---|---|---|
| Spread of a slice (`std_dev`) | 17.98 | 15.04 |
| Earlier method | 13.16 | 11.02 |
| z neighbours | 12.16 | 10.17 |
| **z, 4 apart** | **16.38** | **13.70** |
| In counts | 4 211 | 3 521 |

Estimate from z differences by distance:

| Distance | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| `Sample1_18p15_5mu` | 12.14 | 15.63 | 16.27 | 16.38 | 16.41 | 16.43 |
| `Sample2_5mu` | 10.16 | 13.08 | 13.61 | 13.70 | 13.72 | 13.72 |
| Share of noise in common | 45 % | 9 % | 1 % | – | – | – |

What the numbers say:

- **The noise is correlated over about two voxels.** The estimate rises by
  35 % from distance 1 to distance 4 and is level after that. The shares in
  the last row follow from the ratio to the level and are the same for both
  scans to within one percentage point. The x and y columns show the same
  rise over the same distances, so the correlation is alike in all three
  directions.
- **The level is noise, not material.** Over 40 slices the z estimate rises
  further only slowly, from 13.59 at distance 4 to 13.72 at 10, 13.91 at 20
  and 14.23 at 40 (232³ crop of the second scan). That slow rise is the
  material changing along the fibres. Carried back to distance 4 it amounts
  to less than 1 %.
- **Single voxels carry more noise than signal.** Taking the noise out of
  the slice spread leaves a material contrast of about 7.4 (first scan) and
  6.2 (second scan) as standard deviation, less than half the noise. At 5 µm
  voxels a 20 µm fibre is four voxels wide and mostly edge.
- **The noise is nearly the same everywhere in the reliable circle.**

  | Where | `Sample1_18p15_5mu` | `Sample2_5mu` |
  |---|---|---|
  | Column 200 × 200 on the axis | 16.26 | 13.59 |
  | Column at radius 250 | 16.39 | 13.70 |
  | Column at radius 400 | 16.51 | 13.78 |
  | Slices 0–99 | 17.04 | 14.25 |
  | Slices 450–549 | 16.39 | 13.70 |
  | Slices 900–999 | 16.36 | 13.69 |
  | Region of the earlier 250³ samples | 16.26 | 13.54 |

  It grows by 1.5 % from the axis to radius 400, and the first hundred
  slices are 4 % noisier than the rest. One sigma per scan is enough.
- **The whole scan including its corners reads lower** (15.70 and 13.16),
  because it includes the corners outside the field of view. The estimate
  should be taken on the region that is filtered.

### 11.9 What this means for the filter

- **Sigma decides the patch size in the original program** (section 4.6):
  up to 15 the patch radius is 1, above 15 it is 2. The earlier method would
  have put both scans in the first row (13.2 and 11.0). The noise found here
  puts the first scan in the second row (16.4) and the second scan in the
  first (13.7). The 17 used earlier is in the second row. Step 1.5 separates
  the patch size from sigma, as chapter 7 already proposed.
- **Which number the filter should get.** The filter compares patches that
  lie several voxels apart. The difference between two such patches contains
  the full noise of both, whatever the correlation between direct
  neighbours. The value at the level (16.4 and 13.7) is therefore the right
  sigma for the filter, not the neighbour value.
- **The phantom results of chapter 7 were obtained with independent
  noise.** Real noise that is shared between neighbours forms small blobs
  two to three voxels wide. A patch of radius 1 is no larger than such a
  blob, so the filter can mistake noise for structure more easily than on
  the phantom. How much this costs is not known yet.

### 11.10 Tests and results on both platforms

`make test`: 47 test cases, 1570 checks, all pass on both machines; the
build with address and undefined-behaviour checks on the second machine
reports nothing. The new cases:

| Case | What it checks |
|---|---|
| `sampleSlices` | First and last slice, rounding, fewer slices than asked for |
| `summarise` | Median, mean and quartiles as numpy gives them; NaN left out |
| Four voxels by hand | All six estimators against values computed by hand; single voxel; no next slice |
| Independent noise | All estimators within 2–3 % of the true sigma; `resid` at √(8/9) of it |
| Material equal in both slices | In-slice estimators read high, z estimators stay right |
| `estimateNoise` | Slice numbers, summaries, regions, wrong arguments |
| NaN | Refused |
| Tracked scan | Fixed values for `Sample2_5mu_c64`, the same on every platform |
| Correlated noise | Noise with known correlation along z: the table by distance shows √(1 − correlation) at distances 1 and 2 and the true sigma from 3 on |

The output of `denoiser noise` for the 232³ crop and for a tracked 64³ scan
is the same on both machines, digit for digit.

### 11.11 Open points and decisions

- **Decided in review: `denoiser run` takes "z differences 4 apart" as its sigma**
  when none is given, estimated on the region to be filtered, and writes all
  three values and the method into the run record. `--sigma` overrides it.
  The earlier method stays in the output for comparison with old runs.
- **Limit of the new estimate.** It relies on the material being nearly
  constant over four slices. That holds for rods scanned with the fibres
  along z. For a sample without such a direction the estimate reads high,
  and the table by distance shows it by rising without a level.
- **Limit of the fixed distance.** Four voxels suit these scans. A scan
  reconstructed differently may need another distance; the table shows where
  the level begins. The distance is a constant in `src/main.cpp`.
- **Decided in review: give the phantom correlated noise** (an option of
  `make_phantom.py`, 45 % and 9 % as measured) and repeat the scores of
  chapter 7 with it at the start of step 1.5, before the experiment with
  separate radii per axis. Otherwise the choice of patch size rests on a
  noise model the scans do not follow.
- **The planned 1 µm scan gives a direct check.** The air in the bore hole
  has no material contrast, so the plain standard deviation there is the
  noise. It should agree with the estimate from z differences. Whether the
  noise is correlated over two voxels there as well depends on how the lab
  reconstructs that scan.
- **Sigma differs between scans** (16.4 and 13.7 with the same nominal
  settings), so it has to be estimated for each scan, not carried over.
- **Step 1.5 (parameter choice and radii per axis) is next.** See chapter 12.

---

## 12. Step 1.5 — Filter settings

Status: **done, committed.**

This step makes every choice that the original program has built in a
setting: patch radius and search radius per axis, `h`, how the patch
distance is counted, and the weight function. It then measures on the
phantom which settings classify glass and epoxy best. As decided in the
review of step 1.4, the phantom first got noise that is shared between
neighbouring voxels, as in the real scans.

### 12.1 Main results

1. **Noise shared between neighbours costs the original filter much of its
   benefit.** Its best result on the fine phantom goes from 1.09 % wrong
   voxels to 1.93 %, on the coarse phantom from 6.9 % to 9.0 %.
2. **A patch that is a single column along z is the best patch found.** With
   a column of 17 voxels and a search window that is narrow across and long
   along z, 1.32 % of the voxels are wrong on the fine phantom and 7.0 % on
   the coarse one. Against the original filter at its best that is a third
   fewer on the fine phantom and a fifth fewer on the coarse one.
3. **That filter no longer depends on getting sigma right.** From sigma 12
   to 16 (true value 14) its result moves between 1.26 % and 1.45 %. The
   original filter changes its patch size at sigma 15 and fails above 25.
4. **It visits a third of the positions.** The window of 7 × 7 × 61 has
   2 989 positions against 9 261 of the original 21 × 21 × 21.
5. **The plain mean along z classifies equally well** (1.26–1.31 % and
   7.0 %). The column filter is ahead in the other scores: it leaves a
   smaller error against the noise-free volume and keeps more of the small
   voids.
6. **The IPOL variant helps cube patches a little and column patches not at
   all.** The original weight function stays the default.

### 12.2 What was added

| File | Purpose |
|---|---|
| `include/nlm_bruteforce.hpp`, `src/nlm_bruteforce.cpp` | Settings for the distance (sum or mean over the patch) and the weight function (rational or exponential); `nlmWeight`, `patchVoxels` |
| `tests/reference/bruteforce_main.cpp` | Options `--patch`, `--search`, `--h-factor`, `--normalise`, `--kernel` |
| `tests/reference/run_reference.py` | Passes such options on to the brute-force program |
| `tests/phantom/make_phantom.py` | Option `--correlated-noise` |
| `run/params.py` | The experiments of this chapter |
| `run/params/linux/*.csv` | Their results: 5 files, 177 rows |
| `tests/test_nlm_bruteforce.cpp` | 3 new test cases |

The program version is now 1.5.0. The command `denoiser` itself is
unchanged in this step; the settings reach it with the fast filter (section
12.9).

### 12.3 The settings

| Setting | Original program | Now |
|---|---|---|
| Patch radius | One number for all axes, chosen by sigma from a table | Three numbers: x, y, z |
| Search radius | One number, from the same table | Three numbers |
| `h` | 0.3 to 0.4 × sigma, from the table | Any factor |
| Patch distance | Sum over the patch | Sum, or mean over the patch |
| Weight | h² / (h² + excess) | That, or exp(−excess / h²) |

"Excess" is the patch distance minus what noise alone would give
(2 × sigma² per voxel). Mean distance with the exponential weight is the
variant of the IPOL paper. With the summed distance the same `h` acts more
sharply the larger the patch is; with the mean it does not.

The table of the original program remains available as `originalParams`,
and with it the brute-force filter still reproduces the original program:
on the coarse phantom both give 8.99 % wrong voxels at sigma 16.

### 12.4 The phantom with shared noise

```
python tests/phantom/make_phantom.py --preset fine --correlated-noise
```

makes `build/phantom/phantom_fine_cn`. The independent noise is averaged
with its neighbours along x, y and z with the weights 0.02, 0.245, 1,
0.245, 0.02 and scaled back to the requested sigma.

| Distance in voxels | 1 | 2 | 3 |
|---|---|---|---|
| Share of noise in common, phantom | 44.6 % | 8.9 % | 0.9 % |
| The same, measured on the real scans (section 11.8) | 45 % | 9 % | 1 % |

`denoiser noise` on the coarse phantom made this way shows the pattern of
the real scans: 10.5 from z neighbours, 14.6 from z differences four apart
(true sigma 14.0). Without the option the phantoms are the same files as
before, byte for byte, and both kinds are identical on both machines.

One thing is not known: how much noise diagonal neighbours share in the
real scans. The phantom's value follows from its construction (20 %).

### 12.5 Experiments

All runs are on the second machine with the brute-force filter, true sigma
14, sigma given to the filter 14 and `h` = 0.4 × sigma unless stated. Radii
are written x y z. "Wrong" is the share of misclassified voxels
(`wrong_all` of section 6.3); "voids" is the share of void contrast that is
left. The files in `run/params/linux/` hold all rows and all scores.

**Fine phantom, shared noise** (fibres 20 voxels wide, as in the coming
1 µm scans):

| Filter | Patch | Search | Wrong | PSNR | Voids | Time |
|---|---|---|---|---|---|---|
| None | | | 12.03 % | 25.2 dB | 1.05 | |
| Original table, sigma 14 | 1 1 1 | 10 10 10 | 2.60 % | 33.4 dB | 0.36 | 124 s |
| Original table, sigma 16 (its best) | 2 2 2 | 10 10 10 | 1.93 % | 36.8 dB | 0.44 | 413 s |
| IPOL variant | 2 2 2 | 10 10 10 | 1.79 % | 36.4 dB | 0.77 | 430 s |
| Cube patch, long window | 2 2 2 | 3 3 20 | 1.79 % | 35.6 dB | 0.58 | 90 s |
| Stretched patch | 1 1 4 | 5 5 20 | 1.61 % | 36.6 dB | 0.47 | 192 s |
| Stretched patch | 2 2 4 | 5 5 20 | 1.57 % | 36.4 dB | 0.58 | 394 s |
| Stretched patch | 1 1 8 | 5 5 20 | 1.39 % | 37.7 dB | 0.55 | 353 s |
| Column of 9 | 0 0 4 | 5 5 20 | 1.62 % | 34.4 dB | 0.45 | 67 s |
| Column of 17 | 0 0 8 | 5 5 20 | 1.33 % | 36.3 dB | 0.40 | 120 s |
| **Column of 17** | 0 0 8 | 3 3 30 | **1.32 %** | 36.4 dB | 0.42 | 73 s |
| Column of 17 | 0 0 8 | 2 2 30 | 1.34 % | 36.2 dB | 0.46 | 37 s |
| Column of 17 | 0 0 8 | 1 1 30 | 1.42 % | 35.6 dB | 0.53 | 14 s |
| Column of 25 | 0 0 12 | 3 3 30 | 1.24 % | 37.2 dB | 0.37 | 105 s |
| Column of 33 | 0 0 16 | 3 3 30 | 1.23 % | 37.6 dB | 0.33 | 142 s |
| Mean along z | | 0 0 8 | 1.31 % | 34.1 dB | 0.27 | |
| Mean along z | | 0 0 16 | 1.26 % | 35.8 dB | 0.14 | |
| Median along z | | 0 0 16 | 1.38 % | 35.0 dB | 0.08 | |

**Coarse phantom, shared noise** (fibres 4 voxels wide, as in the present
5 µm scans):

| Filter | Patch | Search | Wrong | PSNR | Voids | Time |
|---|---|---|---|---|---|---|
| None | | | 18.3 % | 25.2 dB | 0.94 | |
| Original table, sigma 14 | 1 1 1 | 10 10 10 | 12.9 % | 29.2 dB | −0.41 | 13 s |
| Original table, sigma 16 (its best) | 2 2 2 | 10 10 10 | 9.0 % | 31.3 dB | −0.13 | 44 s |
| IPOL variant | 2 2 2 | 10 10 10 | 8.6 % | 31.8 dB | 0.20 | 46 s |
| Stretched patch | 1 1 4 | 5 5 20 | 8.6 % | 31.8 dB | −0.31 | 20 s |
| Stretched patch | 1 1 8 | 5 5 20 | 7.8 % | 32.4 dB | −0.03 | 37 s |
| Column of 9 | 0 0 4 | 5 5 20 | 8.0 % | 31.9 dB | 0.06 | 7 s |
| Column of 17 | 0 0 8 | 5 5 20 | 7.0 % | 33.1 dB | 0.02 | 13 s |
| **Column of 17** | 0 0 8 | 3 3 30 | **7.0 %** | 33.0 dB | 0.02 | 8 s |
| Column of 25 | 0 0 12 | 3 3 30 | 7.1 % | 32.9 dB | −0.03 | 11 s |
| Column of 33 | 0 0 16 | 3 3 30 | 7.4 % | 32.6 dB | −0.05 | 15 s |
| Mean along z | | 0 0 8 | 7.0 % | 32.5 dB | −0.06 | |
| Mean along z | | 0 0 16 | 11.5 % | 29.9 dB | −0.28 | |

**Independent against shared noise,** wrong voxels:

| Filter | Fine, independent | Fine, shared | Coarse, independent | Coarse, shared |
|---|---|---|---|---|
| Patch 1 1 1, sigma 14 | 1.47 % | 2.60 % | 8.5 % | 12.9 % |
| Patch 2 2 2, sigma 14 | 1.36 % | 1.85 % | 8.3 % | 9.1 % |
| Patch 2 2 2, sigma 16 | 1.09 % | 1.93 % | 6.9 % | 9.0 % |
| Column of 9, search 5 5 20 | | 1.62 % | 6.0 % | 8.0 % |
| Mean along z, radius 8 | 0.98 % | 1.31 % | 5.4 % | 7.0 % |

**How the column of 17 with search 3 3 30 reacts to sigma,** and for
comparison the original table:

| Sigma given | 12 | 14 | 16 | 18 | 20 | 25 |
|---|---|---|---|---|---|---|
| Column of 17, fine | 1.45 % | 1.32 % | 1.26 % | | | |
| Column of 17, coarse | 7.5 % | 7.0 % | 6.9 % | | | |
| Column of 9, fine | 1.77 % | 1.62 % | 1.56 % | 1.54 % | | |
| Original table, fine | 2.71 % | 2.60 % | 1.93 % | | 2.75 % | 3.77 % |
| Original table, coarse | 12.9 % | 12.9 % | 9.0 % | | 16.0 % | 36.3 % |

**How the column of 9 reacts to `h`** (search 5 5 20): factors 0.4, 0.8,
1.6 and 3.2 give 1.62, 1.59, 1.57 and 1.59 % on the fine phantom and 8.0,
7.9, 8.1 and 9.3 % on the coarse one.

**The IPOL variant with long patches:** patch 0 0 4 gives 1.63 % (fine) and
8.0 % (coarse), the same as the original weight function; patch 1 1 4 gives
1.56 % and 8.4 % against 1.61 % and 8.6 %.

### 12.6 What the numbers say

- **Small cube patches are the wrong tool for this noise.** A patch of
  radius 1 is as large as one blob of shared noise, so blobs look like
  structure and survive. The damage is largest exactly there: 12.9 %
  against 8.5 % on the coarse phantom.
- **Length along z is what counts.** Along z the material is the same over
  many voxels while the noise is fresh after three, so a long column
  averages the noise out of the comparison without mixing materials. Width
  across the fibres adds cost and, on the coarse phantom, harm: 1 1 8 is
  worse than 0 0 8 there, and 2 2 4 worse than 1 1 4.
- **The best length depends on the voxel size.** On the fine phantom the
  result still improves up to a column of 33 voxels. On the coarse phantom
  17 is best and longer columns are worse, because the fibres there move
  sideways faster per voxel. Requirement R8 (defaults reviewed on the first
  1 µm scan) applies to this number in particular.
- **A narrow, long search window is enough.** Three voxels to each side
  across and 30 along z give the same result as 5 and 20, at 60 % of the
  cost. Two voxels across cost another half and lose little.
- **Sigma and `h` become uncritical.** With the column patch the result
  changes by a few tenths of a percentage point over the whole range tried.
  An error of the noise estimate of 10 to 15 % no longer matters.
- **The mean along z remains as good at classifying.** It is the limit of
  this filter family: a window of one column with all weights equal. What
  the weights add shows in the other scores (about 2 dB better PSNR on
  the fine phantom, more void contrast left), and should show on real rods,
  whose fibres are not as regular as the phantom's. The phantom cannot
  decide this.
- **The window along z has the same upper limit as before** (section 7.5):
  mean and median with radius 16 are clearly worse on the coarse phantom.
  The earlier median with radius 16 at 5 µm was past the best window.

### 12.7 Defaults for fibre rods (decided in review)

| Setting | Default | Reason |
|---|---|---|
| Patch radius | 0 0 8 | Best or equal best on both phantoms |
| Search radius | 3 3 30 | Same result as wider windows; a third of the original's positions |
| Sigma | z differences 4 apart (decided, section 11.11) | |
| `h` | 0.4 × sigma | Uncritical; the original's value |
| Distance, weight | Sum, rational (the original's) | The IPOL variant brings no gain with this patch |

- **Decided in review: this becomes the default of `denoiser run`,** in place of the original
  behaviour. Section 1.8 (decision 3) says the default changes only on a
  clear gain and by your decision. The gain on the phantom is a fifth to a third
  fewer wrong voxels than the original at its best sigma, no dependence on that
  sigma, and a third of the positions. The original table stays available
  as an option, and the regression against the original program is
  unaffected: it compares the fast filter with the brute-force filter at
  equal settings.
- **These defaults are for rods scanned with the fibres along z.** For any
  other sample the original table is the safer choice.
- **Whether the mean along z should be offered as a second method** in the
  program is left open. It costs almost nothing to compute. Phase 4 can
  compare both on real scans once a classification exists.

### 12.8 Tests and results on both platforms

`make test`: 50 test cases, 1795 checks, all pass on both machines; the
build with address and undefined-behaviour checks on the second machine
reports nothing. The new cases:

| Case | What it checks |
|---|---|
| `nlmWeight` | The four combinations of distance and weight function against values computed by hand; patch size from per-axis radii |
| Exponential weight by hand | Three voxels, result computed by hand |
| Mean distance | With the rational weight, the mean distance with `h` equals the summed distance with `h` × √(patch size); the exponential weight gives a different result |

The experiments ran on the second machine only (about two hours with the
brute-force filter).

### 12.9 Open points

- **The fast filter does not exist in the new program yet.** The plan lists
  accuracy (step 1.6) and speed (step 1.7) but no step that writes it.
  Proposal: step 1.6 writes the fast filter with all settings of this
  chapter and checks it against the brute-force filter; the accuracy of the
  summed-area arithmetic is settled there. The command-line options follow
  with `denoiser run` in step 1.8 as planned.
- **A column patch needs no summed-area table in x and y.** The patch sum
  runs along z only, which makes the fast filter simpler and its rounding
  error smaller. Step 1.6 can use that.
- **What the phantom does not show.** Its fibres are perfectly straight
  with one common twist, and its shared noise is built from the correlation
  along the axes only. The ranking of column filter and mean along z has to
  be confirmed on real scans.
- **Cost estimate.** The original program needed about 10 hours for a
  present scan (section 7.3). With a third of the positions the same code
  would need about 3 hours, before any of the speed work of step 1.7.

## 13. Step 1.6 — The fast filter and its accuracy

Status: **done, committed.**

The review of step 1.5 extended this step (section 12.9): besides settling
the accuracy of the summed-distance arithmetic, it writes the fast filter
with all settings of chapter 12 and checks it against the brute-force
filter. The fast filter computes the same filter as the brute-force filter,
but finds every patch distance with running sums instead of summing the
patch voxel by voxel.

### 13.1 Main results

1. **In double precision the fast filter gives the brute-force result
   exactly.** On five volumes with three settings each, the largest
   difference is 0. That meets the acceptance rule of phase 1 (section 1.6)
   with no tolerance at all.
2. **In single precision it is off by at most 0.004 on the 0–255 scale**,
   at most one count after rounding to 16 bits. The original program's
   32-bit running sums were off by up to 2.7 (section 6.5). The difference
   is that the running sums here never reach beyond one block plus the
   patch, so they stay small.
3. **Single precision is the default (decided in review).** It is 25–35 %
   faster on the larger volumes and its error is far below the noise.
   Double precision remains an option and is what the checks against the
   brute-force filter use (section 13.5).
4. **The cost no longer grows with the patch.** The original table at
   sigma 17 (patch 5³) takes about as long as at sigma 10 (patch 3³); the
   brute-force filter needed three to four times as long.
5. **It is not yet faster than the original program with all its workers.**
   With the original settings on the 232³ crop it takes 600 s on the second
   machine, against an estimated 360 s for the original program (section
   7.3). With the rod defaults it takes 144 s in double and 106 s in single
   precision. Per voxel and position it is slow (13.6); that is the task of
   step 1.7.

### 13.2 What was added

| File | Purpose |
|---|---|
| `include/nlm_fast.hpp`, `src/nlm_fast.cpp` | The fast filter `nlmFast` with the options block size and precision |
| `include/nlm_bruteforce.hpp`, `src/nlm_bruteforce.cpp` | The check of the arguments is now the public function `checkNlmInput`, shared by both filters |
| `tests/reference/bruteforce_main.cpp` | Options `--fast`, `--single` (default), `--double`, `--block N` |
| `run/fast_check.py` | The check of this chapter: brute force, fast double and fast single on every volume and setting |
| `run/fast/linux/check.csv` | Its results: 45 rows |
| `tests/test_nlm_fast.cpp` | 5 new test cases |

The program version is now 1.6.0. As planned, the command `denoiser` gets
the filter with `denoiser run` in step 1.8; until then it is reached through
the brute-force program with `--fast`.

### 13.3 Commands

```
make bruteforce
build/denoiser_bruteforce <input> <output> <sigma> --fast
build/denoiser_bruteforce <input> <output> <sigma> --fast --double --block 32
python run/fast_check.py
python run/fast_check.py --volumes Sample2_5mu_c64 --settings rods
```

The program takes the same input and output files as in chapter 6
(headerless big-endian floats, size at the end of the file name) and all
options of chapter 12. `run/fast_check.py` writes one row per volume, setting and computation into
`run/fast/<machine>/check.csv` as soon as it is ready, and removes the
filtered volumes afterwards.

### 13.4 How the fast filter works

The patch distance between voxel u and the voxel shifted by d is the sum of
the squared differences (v(u + q) − v(u + q + d))² over the patch offsets q.
For one shift d, the squared differences are the same numbers for all
voxels, only summed over different boxes. So the filter takes one shift at a
time:

1. Squared differences for all voxels of the block plus the patch margin.
2. Their sum over the patch box, one axis after the other: a running sum
   along x, then along y, then along z. Each costs two additions per voxel
   whatever the patch radius. An axis with patch radius 0 is skipped, so
   the column patch of the rod defaults needs the z pass only.
3. Weight of every voxel from its distance, exactly as in `nlmWeight`, and
   the weighted sums.

The volume is cut into cubic blocks (default 64³). Each block is copied
with a margin of patch plus search radius, edges replicated as in the
brute-force filter, and filtered on its own. Blocks are shared out among the
threads. The shifts run in the order of the brute-force filter, and the
accumulators are always in double precision.

Why the error is small: a running sum adds and removes squared
differences, so its rounding error grows with the size of the sums it
carries. Here a running sum covers at most one block edge plus the patch,
not the whole block as in the original program's summed-area table. In
double precision the remaining rounding never changed a result in any of
the 15 cases; in single precision (24 bits) it shows in the last digits.

Unlike the original program, the result does not depend on the number of
threads. The block size can change it at the level of rounding, because
the running sums start anew at each block.

### 13.5 Accuracy

Fast filter against the brute-force filter, all on the second machine.
"Largest" and "RMS" are differences of the float results on the 0–255
scale; "counts" is the largest difference after rounding to 16 bits.

| Volume | Setting | Double: largest | Single: largest | Single: RMS | Single: counts |
|---|---|---|---|---|---|
| `Sample1_18p15_5mu_c64` | table 10 | 0 | 0.0040 | 3.2·10⁻⁵ | 1 |
| | table 17 | 0 | 0.0011 | 1.4·10⁻⁵ | 1 |
| | rods | 0 | 0.0003 | 5.8·10⁻⁶ | 1 |
| `Sample2_5mu_c64` | table 10 | 0 | 0.0010 | 1.6·10⁻⁵ | 1 |
| | table 17 | 0 | 0.00005 | 2.9·10⁻⁶ | 1 |
| | rods | 0 | 0.00006 | 2.7·10⁻⁶ | 1 |
| `Sample2_5mu_c232` | table 10 | 0 | 0.0017 | 9.3·10⁻⁶ | 1 |
| | table 17 | 0 | 0.00008 | 1.8·10⁻⁶ | 1 |
| | rods | 0 | 0.0001 | 2.1·10⁻⁶ | 1 |
| `phantom_coarse_cn` | table 10 | 0 | 0.0029 | 4.2·10⁻⁵ | 1 |
| | table 17 | 0 | 0.0007 | 1.5·10⁻⁵ | 1 |
| | rods | 0 | 0.0002 | 5.4·10⁻⁶ | 1 |
| `phantom_fine_cn` | table 10 | 0 | 0.0035 | 2.8·10⁻⁵ | 1 |
| | table 17 | 0 | 0.0010 | 8.2·10⁻⁶ | 1 |
| | rods | 0 | 0.0003 | 3.2·10⁻⁶ | 1 |

Settings: "table 10" and "table 17" are the original table at sigma 10
(patch 3³, window 21³) and sigma 17 (patch 5³, window 21³); "rods" are the
defaults of section 12.7 at sigma 14 (patch 1 × 1 × 17, window 7 × 7 × 61).

One count is 1/257 of a step on the 0–255 scale. The noise of the present
scans (sigma 14–16 on that scale) is 3 600–4 100 counts, and the filter
moves a voxel by hundreds to thousands of counts. The classification of
phase 4 cannot see a difference of one count.

**Decided in review: single precision is the default.** The reasons:

- Its error is far below the noise, as above.
- The error does not grow with the scan. A running sum spans one block plus
  the patch, so the bound depends on the block size, not on the volume
  size; a 2000 × 2000 × 975 scan has the same bound as the 232³ crop.
- The result is repeatable: same input and settings give the same bits,
  whatever the number of threads.
- It is 25–35 % faster on the larger volumes (section 13.6) and halves the
  buffer memory.
- The GPU version of phase 5 computes in single precision; using it on the
  CPU already proves that arithmetic on real scans.

Double precision stays as an option for one purpose: it reproduces the
brute-force filter to rounding (difference 0 in all 15 cases), so a check
against the reference can demand a difference of 0, and any difference
points to a bug. The tests and `run/fast_check.py` use it for that; single
precision is checked against the same reference with a tolerance of 0.01.

### 13.6 Run time

Second machine, all threads, one run each:

| Volume | Setting | Brute force | Fast double | Fast single |
|---|---|---|---|---|
| 64³ test scan | table 10 | 13 s | 12 s | 12 s |
| | table 17 | 44 s | 15 s | 14 s |
| | rods | 8 s | 4–5 s | 4–5 s |
| `phantom_fine_cn` (160 × 160 × 96) | table 10 | 122 s | 107 s | 69 s |
| | table 17 | 413 s | 114 s | 73 s |
| | rods | 73 s | 27 s | 18 s |
| 232³ crop | table 10 | 625 s | 598 s | 402 s |
| | table 17 | 2 728 s | 623 s | 418 s |
| | rods | 373 s | 144 s | 106 s |

On the development machine (all threads), 232³ crop: rods 132 s in double,
116 s in single, 102 s in double with blocks of 32³; table 17 398 s.

What the numbers say:

- **Small volumes use one thread.** A 64³ volume is a single block, so the
  fast filter runs on one thread there while the brute-force filter uses
  all of them. Per thread it is 17 to 48 times faster than brute force, but
  the wall time hardly shows it.
- **Per voxel and position it is slow.** On the 232³ crop with the original
  settings it spends about 80 ns per voxel and position per thread on the
  second machine (about 40 ns on the development machine). The arithmetic is
  a few operations; the time goes into memory traffic. A block of 64³ with
  its margin needs well over 10 MB of buffers per thread, all of it touched
  for every shift (input, squared differences, the sums along each axis,
  three accumulators in double precision),
  which is far beyond the processor cache. Smaller blocks were already faster
  (102 s against 132 s), and single precision, with half the buffer size,
  is faster by a quarter to a third.
- **Against the original program.** With all its workers the original
  program needs about 360 s for the 232³ crop (29 s per million voxels,
  section 7.3) and 3.9 GB of memory. The fast filter is at 600 s with
  the same settings, so the original's speed is not reached yet. With the
  rod defaults, which visit a third of the positions, it is at 106–144 s.
- **Projection.** With the rod defaults, blocks of 32³ and single
  precision, the reliable region of a present scan (750 × 750 × 1000) would
  take about an hour and a quarter on the development machine; the
  original program needs about 10 hours for a whole scan.
- **Memory.** Peak 300–400 MB on the 232³ crop, of which the volume in and
  out takes about 100 MB. The original program needs 3.9 GB for any volume.

### 13.7 Functions

| Function | Purpose |
|---|---|
| `nlmFast(volume, size, params, options)` | The filter, for all settings of `NlmParams`; `NlmFastOptions` sets block size (default 64) and precision (default single) |
| `checkNlmInput(volume, size, params)` | Rejects a size, radius or setting out of range and a vector whose length does not match the size; used by both filters |

### 13.8 Tests and results on both platforms

`make test`: 55 test cases, 2489 checks, all pass on both machines; the
build with address and undefined-behaviour checks on the second machine
reports nothing. The new cases:

| Case | What it checks |
|---|---|
| Many settings | 13 combinations against brute force, in double (to 10⁻⁴) and single precision (to 0.01), on an 11 × 9 × 13 volume cut into blocks of 4: cube, column, x-only, y-only and one-voxel patches, uneven radii, a patch larger than the block, a window larger than the volume or of one voxel, mean distance, exponential weight and the IPOL variant |
| Tracked test scan | A 24³ part of `Sample2_5mu_c64` at sigma 17 with the original table, against brute force: double precision to 10⁻⁴, the default (single) to 0.01 |
| Block size and threads | In double precision, blocks of 5 and one block agree to 10⁻⁴; the same call twice gives identical bits whatever the thread timing |
| Constant volume and range | A constant volume stays exactly constant; results stay within the input's range |
| Bad arguments | A size that does not match the vector, h = 0, a negative radius and a block size of 0 are rejected |

The check of section 13.5 ran on the second machine (about two hours, most
of it the brute-force filter on the 232³ crop).

### 13.9 Open points

- **Speed is for step 1.7**, and this step shows where to start: buffers
  that fit the cache (smaller blocks, or several shifts per pass over a
  smaller tile), threads within a block so small volumes use all of them,
  and single precision for the running sums. The goal set by the baseline
  is the original program's 29 s per million voxels with the original
  settings.
- **Speed gap between the precisions.** Single precision is faster now
  mainly because its buffers are half the size and memory traffic is the
  bottleneck. If step 1.7 makes the buffers fit the cache, the gap may
  shrink; the default stays single either way, since it costs nothing in
  accuracy that matters.
