FocalMechClassifier
===================

For reproducibility, we demo the code using the dataset provide
with the USGS HASH code, http://earthquake.usgs.gov/research/software/index.php.

Download the hash.v1.2.tar.gz, unpack it and copy the north1.phase
and scsn.reverse datafiles into the "demo" directory of this project.
This will be the inputs to the demo test[2-6].py, run these from the
demo directory. 

You'll also need to compile the hash_driver1 exe is you want to run
it. You can use this to regenerate the example1.inp, example1-norev.inp
data. Do,

cat test1.out | sed -e s/*//" > example1.out

this removes the extra column that messes up genfromtxt.

We provide the input file example1.inp, example1-norev.inp input
parameters. And also data derived from them,  3146815_realizations.out,
example1.out and example1-noreverse.out. You can regenerate these using
the inputs provided.

If you get "LookupError: unknown encoding: " errors from ObsPy, set these environment variables:

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

This seemed to effect OSX: El Capitan.
