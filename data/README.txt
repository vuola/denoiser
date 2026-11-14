Local CT data and historical tools
==================================

Nothing in this folder is tracked by git except this file (see .gitignore).

RAW/  -  scans as delivered by the CT lab (the pipeline's input)
----------------------------------------------------------------
One folder per scan:

  .vol  voxels, headerless, 16-bit unsigned little-endian,
        x fastest, then y, then z
  .vgi  text metadata: size, data type, voxel size
  .pcr  reconstruction settings
  .pca .pcj .pcp .dtxml  acquisition data, not used

Present scans, both 1142 x 1142 x 1000 voxels, 5 micrometre voxels, 2.6 GB:

  Sample1_18p15_5mu   (called "SampleB" in the old workflow)
  Sample2_5mu         (called "SampleD" in the old workflow)

To add a scan, copy its folder here.

Small tracked test scans
------------------------
tests/data/<scan>_c64/ holds a 64 x 64 x 64 cube of each scan as .vol + .vgi
(0.5 MB each, tracked by git). The cube starts at voxel (539, 539, 467),
zero-based. They are regenerated from RAW/ with:

  python tests/data/make_crops.py

Historical material (kept for reference, not used by the pipeline)
-------------------------------------------------------------------
macros/   ImageJ macros of the old manual workflow (PickCube.ijm etc.)
scripts/  Python scripts of the old workflow (scaleToRange.py,
          estimate_noise.py etc.); they need numpy

Sample?_BE255_250x250x250.raw
          Float samples made with the old workflow: 250^3 voxels, 32-bit
          big-endian float, values 0...255.
            SampleA = Sample1_18p15      (20 micrometre scan, not in RAW/)
            SampleB = Sample1_18p15_5mu
            SampleC = Sample2            (20 micrometre scan, not in RAW/)
            SampleD = Sample2_5mu
          B and D were verified against RAW/, all voxels:
            crop : 250^3 cube starting at voxel (446, 446, 374), zero-based
                   (PickCube.ijm with the default centre)
            scale: value = count * 255 / 65535, no clipping
          The largest difference is 0.000008 (float rounding).
          These files are no longer needed and can be deleted (239 MB).

The old workflow and the pipeline that replaces it are described in chapter 2
of the main README.md.
