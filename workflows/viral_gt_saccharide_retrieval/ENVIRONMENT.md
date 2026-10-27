# Environment

Portable runtime: Python >=3.11; numpy==1.26.4; standard library otherwise. Observed: CPython 3.11.17, NumPy 1.26.4, Ubuntu WSL2 Linux x86_64. See environment/observed_validation_environment.json. An existing environment was used; no fresh environment/OS recreation is claimed.

Historical baseline inference recorded NumPy 1.26.4; later Windows paired/rescue records used NumPy 2.3.5. Explicit tolerances address cross-runtime differences. Validation limited OpenBLAS/OMP/MKL threads to one.

No torch, fair-esm, lmdb, transformers, BGC-Prophet model, GPU, structure predictor or taxonomy classifier is required. Original validate.py contains broader historical Torch checks; it is archived as evidence and not invoked.
