FROM containers.mathworks.com/matlab-runtime:r2023a

WORKDIR /app
ENV AGREE_TO_MATLAB_RUNTIME_LICENSE='yes'

COPY fissile_compiled/ /app/fissile_compiled/
RUN chmod +x /app/fissile_compiled/run_fissile.sh

ENTRYPOINT ["/app/fissile_compiled/run_fissile.sh", "/opt/matlabruntime/R2023a"]
