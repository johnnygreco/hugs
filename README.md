# Hunt for Ultra-LSB Galaxies (hugs) 

This repository contains a low-surface-brightness galaxy detection pipeline for the [Hyper Suprime-Cam Subaru Strategic Program (HSC-SSP)](http://hsc.mtk.nao.ac.jp/ssp/). 

### Dependencies 
- [numpy](http://www.numpy.org)
- [scipy](https://www.scipy.org)
- [pandas](http://pandas.pydata.org)
- [astropy](http://www.astropy.org)
- [LSST Science Pipelines](https://pipelines.lsst.io/install/conda.html) (for image processing)
- [SExtractor](https://www.astromatic.net/software/sextractor) version 2.19.5 (for source detection)
- [sqlalchemy](https://www.sqlalchemy.org) (for pipeline database)
- [sfdmap](http://github.com/kbarbary/sfdmap) (for extinction corrections)

## Using hugs

### For HSC-SSP team members
If you have access to the calibrated exposures from a complete HSC rerun (e.g., at Mitaka, IAA, or Princeton), which can be called by the data [butler](https://lsst-web.ncsa.illinois.edu/doxygen/x_masterDoxyDoc/classlsst_1_1daf_1_1persistence_1_1butler_1_1_butler.html), you can use the `runner.py`script to run the *hugs* pipeline on a single patch:

```shell
python scripts/runner.py -t $tract -p $patch
```

or a list of patches stored in a csv file:

```shell
python scripts/runner.py --patches_fn $filename
```

The pipeline parameters are given in a yaml configuration file. The default config file is [here](https://github.com/johnnygreco/hugs/blob/master/pipe-configs/default_config.yml). The rerun directory is given by the `data_dir` parameter. You can pass a custom config file to `runner.py` with the `--config_fn` command line argument.

The pipeline builds a catalog as a sqlite database using [sqlalchemy](https://www.sqlalchemy.org). The database will be created within the `hugs_io` directory, which is specified in the config file. The database model file is [here](https://github.com/johnnygreco/hugs/blob/master/hugs/database/tables.py).

The pipeline can be run in parallel mode using [schwimmbad](https://github.com/adrn/schwimmbad) via the `--mpi` (for mpi) or `--ncores` (for multiprocessing) arguments.

### For everyone
If you are working with HSC-SSP images, the [LSST stack](https://pipelines.lsst.io/install/conda.html) provides an extremely useful suite of image processing tools. *hugs* uses the LSST codebase to "clean" images (i.e., remove bright sources and their associated diffuse light) and *SExtractor* to extract sources. For an example of how to use *hugs* on public HSC images, check out [this notebook](https://github.com/johnnygreco/hugs/blob/master/notebooks/hugs-demo.ipynb).
