# Environment

Python >=3.11 standard library only. Validated with CPython 3.11.17 using -I -S -B, without importing installed third-party scientific packages or creating an environment.

Runtime modules: argparse, collections, csv, hashlib, json, pathlib, re, sys. Public acquisition additionally uses urllib.request and zipfile. Validation uses unittest/importlib and an audit hook. No model or embedding libraries, random seeds, GPUs or large databases are involved.

The isolated run uses an existing Python executable; isolation is achieved through disabled user environment/site paths, fresh archive extraction and an explicit filesystem/network audit guard. It is not represented as a newly provisioned OS/container environment.
