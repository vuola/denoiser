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
| 1.6 | **Numerical accuracy.** Measure the float32 integral-image error; fix it if it matters | Tested accuracy statement |
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
