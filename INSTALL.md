# Containerizing MATLAB

Set up your code as a single MATLAB function with a medical image/file as the input and output. Something like this:

```
function fissile(input_file, output_file, fracFitErrorThresh, eta2, max_fibers, n_jobs)
    ...
end
``` 

You could test this with something like this in bash:

```bash
matlab -batch "fissile('test_input.mat', 'test_output.mat', 5e-2, 10)"
```

Your code cannot have any `figure`/visualizations or pre-compiled code. You also cannot use any functions that modify the MATLAB $PATH like `addpath`. After you have made the appropriate edits, compile the code in the MATLAB GUI. Use `-a` to add any folders that you need to add to the path:

```
mcc -m -d 'fissile_compiled' -a 'splitFod4_mrtrixTest_codeBundle' splitFod4_mrtrixTest_codeBundle/fissile.m
```

Now you can test with something like this in bash. You must have installed the correct MATLAB runtime as the first argument:

```
fissile_compiled/run_fissile.sh /usr/local/MATLAB/MATLAB_Runtime/R2023a test_input.mat test_output.mat 5e-2 5 10 10
```

Now write a Dockerfile based on a container that contains the correct version of MATLAB runtime:

```
FROM containers.mathworks.com/matlab-runtime:r2023a

WORKDIR /app

ENV AGREE_TO_MATLAB_RUNTIME_LICENSE='yes'
COPY fissile_compiled/ /app/fissile_compiled/

RUN chmod +x /app/fissile_compiled/run_fissile.sh

ENTRYPOINT ["/app/fissile_compiled/run_fissile.sh", "/opt/matlabruntime/R2023a/"]
```

Build with:

```
docker build -t fissile:v1.1.3 .
```

Run with:

```
docker run --rm -it -v "$(pwd)/data:/data" fissile:v1.1.3 /data/test_input.mat /data/test_output.mat 5e-2 5 10 10
```

Convert to singularity:

```
apptainer build fissile_v1.1.3.sif docker-daemon://fissile:v1.1.3
```

And run:

```
apptainer run -e --contain -B "$(pwd)/data:/data" fissile_v1.1.3.sif /data/test_input.mat /data/test_output.mat 5e-2 5 10 10
```
