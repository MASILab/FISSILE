# FISSILE
Fiber segmentation by signal splitting and iterative least squares estimation (FISSILE) is a method for unmixing crossing white matter fibers in diffusion MRI data.

## Usage

FISSILE operates on a fiber orientation distribution (FODs) as its input. We suggest fitting using [MRtrix3](https://www.mrtrix.org/) with order 6 or order 8 spherical harmonics. The input should be a `.mat` file containing `signal`, a `N_voxels x N_sh` matrix of spherical harmonics coefficients. We suggest normalizing the FODs so that the first value is `1/sqrt(4*pi)` so the FOD integrates to 1 over the surface of the sphere.

The main code is in `fissile.m`:

```
function fissile(input_file, output_file, fracFitErrorThresh, eta2, max_fibers, n_jobs)
    ...
end
``` 

- `fracFitErrorThresh`: The threshold for the fractional fit error (5e-2).
- `eta2`: The weight of the axial symmetry constraint (5).
- `max_fibers`: The maximum number of fibers to be segmented (10).
- `n_jobs`: The number of parallel jobs to be used (10).

You could run this with something like this in bash:

```bash
matlab -batch "fissile('test_input.mat', 'test_output.mat', 5e-2, 5, 10, 10)"
```

The output file should contain `data`, a MATLAB struct with the following fields:
- `total_fod`: the original FODs
- `separate_fod`: the separated FODs in fiber frames (i.e., aligned with the z-axis)
- `est_theta`: estimated theta value in radians for each fiber
- `est_phi`: estimated phi value in radians for each fiber
- `est_v`: estimated v value for each fiber
- `rotate_separate_fod`: the rotated separated FODs in lab frame

## Containerized version

There is also a containerized version. See `INSTALL.md` for building instructions.

Run with:

```
docker run --rm -it -v "$(pwd)/data:/data" fissile:v1.1.3 /data/test_input.mat /data/test_output.mat 5e-2 5 10 10
```

Or as a singularity, available [here](https://doi.org/10.5281/zenodo.20044172) for download:

```
apptainer run -e --contain -B "$(pwd)/data:/data" fissile_v1.1.3.sif /data/test_input.mat /data/test_output.mat 5e-2 5 10 10
```

## Citation

Adam M. Saunders, Gaurav Rudravaram, Michael E. Kim, Elyssa M. McMaster, Yifan Shuai, Warren Boschen, Sophie Lin, Suzana Herculano-Houzel, Daniel C. Colvin, Yurui Gao, Kurt G. Schilling, Adam W. Anderson, and Bennett Landman. "Separating mixtures of single-fiber orientation distributions for increased sensitivity to crossing white matter fibers.” [https://arxiv.org/abs/2610.20974](https://arxiv.org/abs/2610.20974).
