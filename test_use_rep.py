#!/usr/bin/env python
"""
Test script for use_rep parameter functionality.
Tests that ESCHR can accept user-supplied dimensional reductions.
"""

import numpy as np
import pandas as pd
import anndata
import eschr as es
import os
import shutil

print("Testing use_rep parameter with pre-computed dimensional reduction...")
print("=" * 70)

# Create test data
np.random.seed(42)
n_cells = 200
n_genes = 1000

# Generate synthetic gene expression data
data = np.random.rand(n_cells, n_genes) * 10

# Create AnnData object
adata = anndata.AnnData(X=data)
print(f"Created AnnData with shape: {adata.X.shape}")

# Compute PCA manually (simulating user-supplied reduction)
from sklearn.decomposition import PCA
pca = PCA(n_components=30, random_state=42)
X_pca = pca.fit_transform(adata.X)
adata.obsm['X_pca'] = X_pca
print(f"Added pre-computed PCA to adata.obsm['X_pca'] with shape: {X_pca.shape}")

# Test 1: Run ESCHR with default (raw data)
print("\n" + "=" * 70)
print("TEST 1: Running ESCHR with raw data (use_rep='X')")
print("=" * 70)
zarr_loc_1 = "./test_zarr_raw.zarr"
if os.path.exists(zarr_loc_1):
    shutil.rmtree(zarr_loc_1)

try:
    adata_test1 = adata.copy()
    adata_test1 = es.tl.consensus_cluster(
        adata=adata_test1,
        zarr_loc=zarr_loc_1,
        use_rep='X',  # Use raw data
        ensemble_size=10,  # Small for testing
        nprocs=1  # Serial execution to avoid multiprocessing issues
    )
    print(f"✓ SUCCESS: Clustering completed with {len(adata_test1.obs['hard_clusters'].unique())} clusters")
    print(f"  Hard clusters shape: {adata_test1.obs['hard_clusters'].shape}")
    print(f"  Soft membership matrix shape: {adata_test1.obsm['soft_membership_matrix'].shape}")
except Exception as e:
    print(f"✗ FAILED: {str(e)}")
    import traceback
    traceback.print_exc()
finally:
    if os.path.exists(zarr_loc_1):
        shutil.rmtree(zarr_loc_1)

# Test 2: Run ESCHR with pre-computed PCA
print("\n" + "=" * 70)
print("TEST 2: Running ESCHR with pre-computed PCA (use_rep='X_pca')")
print("=" * 70)
zarr_loc_2 = "./test_zarr_pca.zarr"
if os.path.exists(zarr_loc_2):
    shutil.rmtree(zarr_loc_2)

try:
    adata_test2 = adata.copy()
    adata_test2 = es.tl.consensus_cluster(
        adata=adata_test2,
        zarr_loc=zarr_loc_2,
        use_rep='X_pca',  # Use pre-computed PCA
        ensemble_size=10,  # Small for testing
        nprocs=1  # Serial execution to avoid multiprocessing issues
    )
    print(f"✓ SUCCESS: Clustering completed with {len(adata_test2.obs['hard_clusters'].unique())} clusters")
    print(f"  Hard clusters shape: {adata_test2.obs['hard_clusters'].shape}")
    print(f"  Soft membership matrix shape: {adata_test2.obsm['soft_membership_matrix'].shape}")
    print(f"  Used pre-computed reduction with {X_pca.shape[1]} components instead of {n_genes} genes")
except Exception as e:
    print(f"✗ FAILED: {str(e)}")
    import traceback
    traceback.print_exc()
finally:
    if os.path.exists(zarr_loc_2):
        shutil.rmtree(zarr_loc_2)

# Test 3: Test error handling for invalid use_rep
print("\n" + "=" * 70)
print("TEST 3: Testing error handling for invalid use_rep")
print("=" * 70)
zarr_loc_3 = "./test_zarr_invalid.zarr"
if os.path.exists(zarr_loc_3):
    shutil.rmtree(zarr_loc_3)

try:
    adata_test3 = adata.copy()
    adata_test3 = es.tl.consensus_cluster(
        adata=adata_test3,
        zarr_loc=zarr_loc_3,
        use_rep='X_nonexistent',  # Invalid key
        ensemble_size=10,
        nprocs=1  # Serial execution to avoid multiprocessing issues
    )
    print(f"✗ FAILED: Should have raised ValueError for invalid use_rep")
except ValueError as e:
    print(f"✓ SUCCESS: Correctly raised ValueError: {str(e)}")
except Exception as e:
    print(f"✗ FAILED: Wrong exception type: {str(e)}")
finally:
    if os.path.exists(zarr_loc_3):
        shutil.rmtree(zarr_loc_3)

print("\n" + "=" * 70)
print("All tests completed!")
print("=" * 70)
