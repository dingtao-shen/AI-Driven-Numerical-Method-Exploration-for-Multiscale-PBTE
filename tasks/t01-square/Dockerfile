FROM python:3.11-slim

# Everything at build time; grading runs with networking disabled.
RUN pip install --no-cache-dir \
      "numpy>=1.24" "scipy>=1.10" "pyyaml>=6.0" "numba>=0.58" "pytest>=7.0"

# Deterministic and single-threaded, so an iteration count is a fact about the
# solver rather than about the machine it ran on.
ENV OMP_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    NUMBA_NUM_THREADS=1 \
    PYTHONHASHSEED=0 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /task
COPY environment/ /task/environment/
COPY task.md      /task/task.md

# The verifier, its pristine solver, the references and the oracle are mounted
# at grading time, never baked in: the residual gate is only worth something
# if the submission never had access to the code that computes it.

WORKDIR /task/environment
CMD ["/bin/bash"]
