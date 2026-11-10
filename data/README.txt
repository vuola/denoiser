Full CT test volumes (not tracked by git, see .gitignore)
=========================================================

Sample names
------------
SampleA = Sample1_18p15
SampleB = Sample1_18p15_5mu
SampleC = Sample2
SampleD = Sample2_5mu

Files
-----
SampleA_BE255_250x250x250.raw
SampleB_BE255_250x250x250.raw
SampleC_BE255_250x250x250.raw
SampleD_BE255_250x250x250.raw

Format
------
Headerless raw, 250 x 250 x 250 voxels, 62 500 000 bytes each.
32-bit float, big-endian ("BE"), x fastest, then y, then z.
Values scaled to the range 0...255 ("255"). No NaN or Inf.

Measured over the whole volume:

  Sample   min       max       mean      std
  A          0.000   255.000   124.690   23.791
  B         30.000   226.043   128.568   17.787
  C         97.167   249.767   182.613   13.553
  D         27.280   184.852   108.296   15.009

How to get the files on a new machine
-------------------------------------
From the repository root:

  scp "vuola@192.168.68.134:/home/vuola/matres/MicroCT/NEWFILTER/*250x250x250.raw" data/

Small tracked crops
-------------------
tests/data/ holds the centred 64 x 64 x 64 cube of each volume (1 MB each,
tracked by git). They are regenerated from the files here with:

  python tests/data/make_crops.py
