# Reproduce the supplementary embedding figure and global clustering

Clone the [public repository](https://github.com/LanNa1229/viral-gt-research-record) and check out the full commit SHA for the published version of this document. For `RELEASE_COMMIT` entries, the checkout itself provides the immutable release revision. For any other entry, the public input manifest specifies a full commit and repository path. A regular clone contains the applicable history. A sparse or shallow clone also needs the pinned commits and paths listed in the manifest.

Run these commands from the repository root. Python 3.11 or newer is needed for checksum verification.

```sh
cd results/WP07_supplementary_global_BGC_embedding
python provenance/reproduce.py checksums
```

## Redraw the final panel

The archived plotter renders the same source coordinates, membership annotations, class counts and palette to PNG, PDF and SVG. Use Python with NumPy 2.3.5, Matplotlib 3.11.1 and Pillow 12.3.0, and specify a locally available Arial regular font file. Use a new output directory. On Windows, use a short path to stay within the filesystem path-length limit. The script compares the redrawn PNG to the archived PNG. Font revisions, operating system, Matplotlib dependencies and export metadata can affect exact equality.

```sh
python provenance/reproduce.py figure --font /path/to/arial.ttf --output /path/to/new/figure-replay
```

## Verify or replay global Leiden

The public inputs provide 2,502 reference vectors, 556 viral representative vectors, the array of 737 viral window vectors from which the representatives were selected, and the row order of the 287 selected representatives. The checker verifies the input digests, joins by identifiers, rebuilds the weighted cosine union-kNN graph with k=20, and compares the edges, weights, partition objective, saved assignments, reference class counts and figure coordinates with the historical tables.

```sh
python provenance/reproduce.py graph
```

To replay the original Leiden optimizer on the reconstructed graph and verify its per-sample baseline labels, install the archived graph stack in an isolated environment and pass `--replay-leiden`. The recorded versions are NumPy 2.3.5, SciPy 1.18.1, scikit-learn 1.9.0, igraph 1.0.0 and leidenalg 0.12.0. The replay uses the archived floating-point edge weights after verifying the graph from the embeddings. The recorded baseline uses the RB configuration partition, resolution 1.0 and seed 42. The graph and saved assignments are verified without an optimizer replay during publication.

```sh
python provenance/reproduce.py graph --replay-leiden
```

## Restore historical inputs and inspect upstream limitations

```sh
python provenance/reproduce.py restore --output /path/to/new/workspace
```

The restorer checks published hashes and writes the documented relative paths into a new workspace. It extracts and verifies bundled numeric inputs. Consult [the input manifest](provenance/INPUT_MANIFEST.tsv) before invoking historical scripts, some of which use their original workstation configuration. The [scientific bundle manifest](provenance/BUNDLED_SCIENTIFIC_INPUTS.tsv) records checksums for individual archive members. The [software versions](provenance/SOFTWARE_VERSIONS.tsv) are a record of the historical runtimes.

A complete upstream replay also needs the externally distributed [15.1 GB BGC-Prophet protein-vector LMDB, model checkpoints and versioned GenBank records](provenance/EXTERNAL_DEPENDENCIES.tsv). The LMDB’s archived hash is historical. GenBank annotations can change even for the same accession/version, so compare any downloaded bytes with the recorded checksum. Frozen sequence outputs are supplied where used. These inputs are unnecessary for redrawing the final figure or rebuilding the graph from saved vectors.

The effective epoch count of the historical UMAP fit was unrecorded. An estimate from the library default is not an execution record. The exact fit-time script bytes could not be established. Some initial commands, logs and runtime details are incomplete. The frozen numerical outputs and surviving execution records are retained. See the [Methods record](METHODS_RECORD.md) and [validation report](provenance/VALIDATION_REPORT.md).
