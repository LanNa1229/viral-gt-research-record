#!/bin/bash

# ===== USER SETTINGS =====

INPUT_PDB_DIR="/home/ln72030/virus_pdb/virus_initial_hits"  # <<< Change this to your pdb/cif directory
OUTPUT_DIR="/home/ln72030/foldseek_virus_cluster_clean"     # <<< Change this to your desired output directory
THREADS=64

# ==========================

mkdir -p "$OUTPUT_DIR"/{database,cluster,align}

DB_DIR="$OUTPUT_DIR/database"
CLUSTER_DIR="$OUTPUT_DIR/cluster"
ALIGN_DIR="$OUTPUT_DIR/align"

DB_NAME="$DB_DIR/protein_db"
CLUSTER_TMP="$CLUSTER_DIR/tmp"

echo "📦 Step 1: Creating Foldseek database..."
foldseek createdb "$INPUT_PDB_DIR" "$DB_NAME"

echo "📊 Step 2: Clustering proteins..."
foldseek cluster "$DB_NAME" "$CLUSTER_DIR/cluster" "$CLUSTER_TMP" --min-seq-id 0.0 -c 0.05 --cov-mode 0 --threads $THREADS

echo "🔗 Step 3: Aligning representative vs members..."
foldseek align "$DB_NAME" "$DB_NAME" "$CLUSTER_DIR/cluster" "$ALIGN_DIR/aln_result" -a --threads $THREADS

echo "📄 Step 4: Converting alignments to m8 format..."
foldseek convertalis "$DB_NAME" "$DB_NAME" "$ALIGN_DIR/aln_result" "$ALIGN_DIR/aln_result.m8"

echo "✅ ALL DONE!"
echo "Cluster file: $CLUSTER_DIR/cluster_cluster.tsv"
echo "Alignment file (for Cytoscape): $ALIGN_DIR/aln_result.m8"

