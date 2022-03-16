FROM nfcore/base:1.9
LABEL authors="Andrea Telatin" \
      description="Docker image containing all software requirements for the telatin/nextflow-example project"

# Install the conda environment
COPY env.yaml /environment.yml
RUN conda env create -n env -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/env/bin:$PATH
# (alternative: replace the name with $(head -n 1 env.yaml  | cut -f 2 -d " ") or similar

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name env > workflow.yml
