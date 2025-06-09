import h5py
import anndata
import scipy.sparse as sp
import numpy as np
import pandas as pd
import sys

# --- Configuration ---
h5seurat_path = 'SD-A1_processed.h5Seurat'
h5ad_path = 'SD-A1_processed.h5ad'
main_assay_name = 'Spatial.008um' # The name of your assay to use for X

# --- Helper function for printing messages ---
def log_message(message):
    print(f"[INFO] {message}")

def log_warning(message):
    print(f"[WARNING] {message}", file=sys.stderr)

def log_error(message):
    print(f"[ERROR] {message}", file=sys.stderr)
    sys.exit(1) # Exit on critical errors

log_message(f"Starting conversion of {h5seurat_path} to {h5ad_path}")

try:
    with h5py.File(h5seurat_path, 'r') as f_h5seurat:
        log_message("H5Seurat file opened successfully.")

        # --- Step 1 & 2: Verify top-level and assay layers ---
        log_message(f"Top-level keys: {list(f_h5seurat.keys())}")
        
        layers_group_path = f'assays/{main_assay_name}/layers'
        if layers_group_path not in f_h5seurat:
            log_error(f"Error: Layers group '{layers_group_path}' not found in h5Seurat file.")
        
        spatial_layers_group = f_h5seurat[layers_group_path]
        log_message(f"Layers within {main_assay_name}: {list(spatial_layers_group.keys())}")

        # --- Step 3: Identify the primary expression matrix (X) and its dimensions ---
        x_layer_subgroup_name = None
        
        # Prioritize 'data' (normalized) over 'counts' (raw)
        for layer_name in ['data', 'counts']:
            if layer_name in spatial_layers_group:
                layer_group = spatial_layers_group[layer_name]
                if all(key in layer_group for key in ['data', 'indices', 'indptr']):
                    x_layer_subgroup_name = layer_name
                    log_message(f"Using '{x_layer_subgroup_name}' layer as primary X matrix.")
                    break
        
        if x_layer_subgroup_name is None:
            log_error(f"Could not find a valid sparse matrix (data or counts) within {layers_group_path}/.")

        x_layer_full_path = f"{layers_group_path}/{x_layer_subgroup_name}"
        sparse_matrix_group = f_h5seurat[x_layer_full_path]
        
        if 'dims' not in sparse_matrix_group.attrs:
            log_error(f"Error: 'dims' attribute not found on {x_layer_full_path}.")
            
        matrix_dims = sparse_matrix_group.attrs['dims']
        n_features, n_cells = matrix_dims[0], matrix_dims[1] # Seurat (features, cells)

        log_message(f"Dimensions of X matrix (Seurat format: features x cells): {n_features} x {n_cells}")
        log_message(f"Will be transposed to (cells x features) for AnnData: {n_cells} x {n_features}")

        # --- Step 4: Load and reconstruct the X sparse matrix ---
        sp_data = sparse_matrix_group['data'][:]
        sp_indices = sparse_matrix_group['indices'][:]
        sp_indptr = sparse_matrix_group['indptr'][:]

        # Reconstruct CSC (column-major) matrix from Seurat's format
        X_csc = sp.csc_matrix((sp_data, sp_indices, sp_indptr), shape=(n_features, n_cells))
        
        # Transpose to CSR (row-major) for AnnData (cells x features)
        X = X_csc.T.tocsr()
        log_message(f"Shape of reconstructed X matrix (AnnData format: cells x features): {X.shape[0]} x {X.shape[1]}")

        # --- Step 5: Extract Cell Metadata (.obs) ---
        obs_data_dict = {}
        if 'meta.data' in f_h5seurat:
            meta_data_group = f_h5seurat['meta.data']
            for key in meta_data_group.keys():
                item = meta_data_group[key]
                if isinstance(item, h5py.Group):
                    # Handle categorical factors (nested groups with 'values' and 'levels')
                    if "values" in item and "levels" in item:
                        values = item['values'][:].astype(np.int64)
                        levels = item['levels'][:].astype(str)
                        if len(values) == n_cells:
                            obs_data_dict[key] = levels[values]
                            log_message(f"Loaded meta.data/{key} as categorical data.")
                        else:
                            log_warning(f"Metadata '{key}' values length mismatch ({len(values)} vs {n_cells}). Skipping.")
                    else:
                        log_warning(f"Skipping nested group /meta.data/{key} (no 'values' and 'levels').")
                elif isinstance(item, h5py.Dataset):
                    # Handle simple datasets
                    val = item[()]
                    if isinstance(val, np.ndarray) and val.dtype.type is np.bytes_:
                        obs_data_dict[key] = val.astype(str)
                    elif isinstance(val, bytes):
                        obs_data_dict[key] = val.decode('utf-8')
                    else:
                        obs_data_dict[key] = val
                    
                    if np.isscalar(obs_data_dict[key]) and n_cells > 1: # Scalar applied to all cells
                        obs_data_dict[key] = np.full(n_cells, obs_data_dict[key])
                    elif isinstance(obs_data_dict[key], np.ndarray) and obs_data_dict[key].shape != (n_cells,):
                        log_warning(f"Metadata '{key}' shape mismatch ({obs_data_dict[key].shape} vs ({n_cells},)). Skipping.")
                        del obs_data_dict[key] # Remove it if it doesn't fit
                    else:
                        log_message(f"Loaded meta.data/{key} as direct dataset.")
                else:
                    log_warning(f"Skipping unknown type for /meta.data/{key}.")
        
        cell_names = None
        if 'cell.names' in f_h5seurat:
            cell_names_read = f_h5seurat['cell.names'][:].astype(str)
            if len(cell_names_read) == n_cells:
                cell_names = cell_names_read
            else:
                log_warning(f"Length of cell.names ({len(cell_names_read)}) does not match n_cells ({n_cells}). Using default indices.")
        
        if cell_names is None:
            cell_names = np.arange(n_cells).astype(str) # Fallback to generic names

        obs_df = pd.DataFrame(index=pd.Index(cell_names, name='cell_id'))
        for k, v in obs_data_dict.items():
            if isinstance(v, np.ndarray) and v.shape == (n_cells,):
                obs_df[k] = v
            elif np.isscalar(v) and n_cells > 0:
                obs_df[k] = v
            else:
                log_warning(f"Skipping metadata '{k}' due to shape mismatch or unexpected type after processing: {type(v)} {getattr(v, 'shape', 'N/A')}")

        log_message(f"Shape of obs (cell metadata) DataFrame: {obs_df.shape[0]} x {obs_df.shape[1]}")

        # --- Step 6: Extract Feature Metadata (.var) ---
        var_df = pd.DataFrame(index=pd.RangeIndex(n_features, name='gene_id')) # Default
        
        if f'assays/{main_assay_name}/features' in f_h5seurat:
            feature_group = f_h5seurat[f'assays/{main_assay_name}/features']
            
            feature_names = None
            if "name" in feature_group:
                feature_names_read = feature_group['name'][:].astype(str)
                if len(feature_names_read) == n_features:
                    feature_names = feature_names_read
                    var_df = pd.DataFrame(index=pd.Index(feature_names, name='gene_id'))
                    log_message("Using extracted feature names for var.")
                else:
                    log_warning(f'Length of feature names ({len(feature_names_read)}) does not match n_features ({n_features}). Using default indices for var.')
            else:
                log_warning('No "name" dataset found in features group. Using default gene IDs for var.')

            # Add other feature attributes like 'key', 'id'
            for attr_name in ['key', 'id']: # Common Seurat feature attributes
                if attr_name in feature_group:
                    attr_values = feature_group[attr_name][:].astype(str)
                    if len(attr_values) == n_features:
                        var_df[attr_name] = attr_values
                        log_message(f"Added feature attribute '{attr_name}' to var.")
                    else:
                        log_warning(f"Feature attribute '{attr_name}' length mismatch ({len(attr_values)} vs {n_features}). Skipping.")
        else:
            log_warning('No features group found for main assay. Using default gene IDs for var.')

        log_message(f"Shape of var (feature metadata) DataFrame: {var_df.shape[0]} x {var_df.shape[1]}")

        # --- Step 7: Create the AnnData object ---
        adata = anndata.AnnData(X=X, obs=obs_df, var=var_df)
        log_message("AnnData object created.")
        log_message(f"AnnData object info:\n{adata}")

        # --- Step 8: Add other layers (e.g., raw counts, scaled data) ---
        if x_layer_subgroup_name == 'data' and 'counts' in spatial_layers_group:
            counts_group_in_layers = spatial_layers_group['counts']
            if all(key in counts_group_in_layers for key in ['data', 'indices', 'indptr']):
                sp_data_counts = counts_group_in_layers['data'][:]
                sp_indices_counts = counts_group_in_layers['indices'][:]
                sp_indptr_counts = counts_group_in_layers['indptr'][:]
                
                counts_csc = sp.csc_matrix((sp_data_counts, sp_indices_counts, sp_indptr_counts), shape=(n_features, n_cells))
                adata.layers['counts'] = counts_csc.T.tocsr() # Transpose for AnnData
                log_message("Added 'counts' layer to AnnData object.")
            else:
                log_warning("Counts layer found, but not in expected sparse matrix format. Skipping.")

        if 'scale.data' in spatial_layers_group:
            scaled_data = spatial_layers_group['scale.data'][:]
            # scale.data is dense and has dimensions (n_features, n_cells) from Seurat
            adata.layers['scaled_data'] = scaled_data.T # Transpose for AnnData
            log_message("Added 'scaled_data' layer to AnnData object.")

        # --- Step 9: Add Spatial Embeddings/Reductions (.obsm) ---
        if 'reductions' in f_h5seurat:
            reductions_group = f_h5seurat['reductions']
            for reduction_name in reductions_group.keys():
                reduction_path = f'reductions/{reduction_name}/cell.embeddings'
                if reduction_path in f_h5seurat:
                    embeddings = f_h5seurat[reduction_path][:]
                    # Embeddings are typically (n_cells, n_dimensions) already in Seurat
                    # Verify shape before adding to avoid issues
                    if embeddings.shape[0] == n_cells:
                        adata.obsm[f'X_{reduction_name.lower()}'] = embeddings
                        log_message(f"Added reduction '{reduction_name}' to .obsm as X_{reduction_name.lower()}.")
                    else:
                        log_warning(f"Reduction '{reduction_name}' embeddings shape mismatch ({embeddings.shape[0]} vs {n_cells}). Skipping.")
                else:
                    log_warning(f"Skipping reduction '{reduction_name}': no 'cell.embeddings' found.")

        # --- Step 10: Save the AnnData object ---
        adata.write(h5ad_path)
        log_message(f"Successfully converted {h5seurat_path} to {h5ad_path}")

except Exception as e:
    log_error(f"An unexpected error occurred during conversion: {e}")

# The h5py file object 'f_h5seurat' is automatically closed when exiting the 'with' block.
log_message("H5Seurat file connection closed.")