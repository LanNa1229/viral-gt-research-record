# Environment

Validated on Linux x86_64 using the existing recovery environment: Python 3.12.11, PyTorch 2.5.1+cpu, NumPy 1.26.4, fair-esm 2.0.0. Historical extraction additionally recorded Biopython 1.85; the portable Step 2 path reads already extracted TSV sequences and does not perform extraction. CPU float32, four intraop threads, one interop thread and seed 1729 are set by the unchanged historical torch_setup function. Classifier autograd behavior and eval mode are preserved; ESM encoding uses the original inference_mode path.

The minimal dependency recipe is environment/requirements-cpu.txt. environment/observed_validation_environment.json records the existing validation environment and its installed distributions; this is an observed snapshot, not a retrospectively recovered historical lock, wheel bundle or cross-platform guarantee. No installation, new Conda environment or GPU operation was performed.

Users creating a separate environment may install the requirements using their normal approved package workflow. Full historical CPU replay and reproduction on another machine remain untested. Model/source versions are pinned independently by resources.json; the ESM contact-regression sidecar is required by the installed loader even though contacts are not requested.
