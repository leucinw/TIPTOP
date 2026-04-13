# ============================================================
# TIPTOP – IP_PDB2txyz Docker image
#
# Dependencies:
#   obabel   – converts biomolecule PDB fragments to TINKER XYZ
#   openmm   – used by the pdb_fixer.py helper (PDB preparation)
#   pdbfixer – adds missing atoms / hydrogens before IP_PDB2txyz
#
# The TIPTOP repository is cloned from GitHub at build time so
# the image always reflects the latest main branch.  Activate
# the conda environment before running any TIPTOP script:
#
#   docker run --rm -it tiptop
#   conda activate tiptop
#   cd /opt/TIPTOP
#   python IP_PDB2txyz.py your.pdb A
# ============================================================

FROM continuumio/miniconda3:latest

LABEL maintainer="Chengwen Liu <liuchw2010@gmail.com>"
LABEL description="TIPTOP IP_PDB2txyz environment"

# ── system packages ──────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
        git \
    && rm -rf /var/lib/apt/lists/*

# ── conda environment ────────────────────────────────────────
# openbabel is installed via conda-forge to get the CLI (obabel)
RUN conda create -n tiptop -y -c conda-forge \
        python=3.11 \
        openbabel \
        openmm \
        pdbfixer \
    && conda clean -afy

# Make the tiptop env the default for all subsequent RUN/CMD layers
ENV PATH="/opt/conda/envs/tiptop/bin:$PATH"
SHELL ["conda", "run", "-n", "tiptop", "/bin/bash", "-c"]

# ── clone TIPTOP from GitHub ─────────────────────────────────
RUN git clone https://github.com/leucinw/TIPTOP.git /opt/TIPTOP

WORKDIR /opt/TIPTOP

# ── activate tiptop env automatically on interactive login ───
RUN echo "conda activate tiptop" >> /root/.bashrc

CMD ["/bin/bash", "--login"]
