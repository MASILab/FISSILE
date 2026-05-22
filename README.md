# FISSILE
The main code is in `fissile.m`:

```
function fissile(input_file, output_file, fracFitErrorThresh, eta2, max_fibers, n_jobs)
    ...
end
``` 

You could run this with something like this in bash:

```bash
matlab -batch "fissile('test_input.mat', 'test_output.mat', 5e-2, 5, 10, 10)"
```

There is also a containerized version. See `INSTALL.md` for building instructions.

Run with:

```
docker run --rm -it -v "$(pwd)/data:/data" fissile:v1.1.2 /data/test_input.mat /data/test_output.mat 5e-2 5 10 10
```

Or as a singularity:

```
apptainer run -e --contain -B "$(pwd)/data:/data" fissile_v1.1.2.sif /data/test_input.mat /data/test_output.mat 5e-2 5 10 10
```
