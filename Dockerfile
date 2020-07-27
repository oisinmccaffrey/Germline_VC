FROM nfcore/base:1.9
LABEL authors="Barry Digby" \
      description="Docker image containing all software requirements for Germline Variant Calling Analysis"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/Germline_VC/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name Germline_VC > Germline_VC.yml
