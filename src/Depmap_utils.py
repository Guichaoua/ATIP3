import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests
import gseapy as gp




def get_mutated_genes(row):
    """Extract mutated gene names from a DataFrame row.
    
    Args:
        row (pd.Series): A row from the DataFrame containing genetic feature information.

    Returns:
        list: A list of mutated gene names.
    """
    mutated_genes = []

    # Check if the feature is marked as mutated
    if row['IS Mutated'] == 1:
        
        # Check if the gene name contains "_mut" in the 'Genetic Feature' column
        if '_mut' in row['Genetic Feature']:
            gene_name = row['Genetic Feature'].replace('_mut', '')
            mutated_genes.append(gene_name)
    
    # Return a list of unique gene names (remove duplicates)
    return list(set(mutated_genes))

def load_metadata():
    """
    Load and preprocess metadata for cell lines.
    
    Args:
        None

    Returns:
        pd.DataFrame: A DataFrame containing the processed metadata for cell lines.
    """
    try: 
        Genetic_features = pd.read_csv("../data/Cell-lines/PANCANCER_Genetic_features_2025_04_07.csv")
        Genetic_features.rename(columns={'COSMIC ID': 'COSMIC_ID'}, inplace=True)
        # Apply the function to each row of the DataFrame to create the 'Mutated Genes' column
        Genetic_features['Mutated Genes'] = Genetic_features.apply(get_mutated_genes, axis=1)
        grouped_Genetic_features = Genetic_features.groupby(['Cell Line Name', 'COSMIC_ID']).agg({
            'GDSC Desc1': 'first',
            'GDSC Desc2': 'first',
            'TCGA Desc': 'first',
            'Mutated Genes': lambda x: list(set([gene for sublist in x for gene in sublist]))
        }).reset_index()
    except Exception as e:
        print(f"Error loading PANCANCER_Genetic_features: {e}")
        return pd.DataFrame()

    try:
        meta = pd.read_csv("../data/Cell-lines/Model.csv") 
        meta.rename(columns={'ModelID': 'DepMap_ID', 'COSMICID':'COSMIC_ID' }, inplace=True)
    except Exception as e:
        print(f"Error loading Model.csv: {e}")
        return pd.DataFrame()  # Return empty DataFrame on error

    return meta, grouped_Genetic_features



def load_merge_data( dict_filters, list_mutation, drop_mutation=False, keep_only_mutation=False, expr_keep=[]):
    """
    Load and merge omics and mutation data according to input parameters.

    Args:
        dict_filters (dict): Filtering conditions as {column_name: lambda}.
        list_mutation (list): Genes to use for mutation-based filtering.
        drop_mutation (bool): If True, keep only non-mutated samples.
        keep_only_mutation (bool): If True, keep only mutated samples.
        expr_keep (list): List of expression columns to retain.

    Returns:
        pd.DataFrame: Merged and filtered dataset ready for analysis.
    """

    # 1. Load expression data depending
    # RNA-seq data (log2(TPM+1)), filter for selected genes
    try : 
        df_1 = pd.read_csv("../data/Cell-lines/OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv", index_col=0) #DepMap Public 24Q4
        df_1.index.name = 'DepMap_ID'
        df_1 = df_1.reset_index()
        df_1.columns = df_1.columns.str.extract(r'^([^\(]+)', expand=False).str.strip()
        df_1.rename(columns={'COSMICID': 'COSMIC_ID'}, inplace=True)
    except Exception as e:
        print(f"Error loading RNA-seq data: {e}")
        return pd.DataFrame()  # Return empty DataFrame on error
    if expr_keep:
        df_1 = df_1[['DepMap_ID'] + expr_keep]

    # 2. Load metadata
    meta, grouped_Genetic_features = load_metadata()
    meta_grouped_Genetic_features = meta.merge(grouped_Genetic_features, on='COSMIC_ID', how='inner')#,how = 'left')
    #Add metadata (e.g., COSMIC_ID, tissue type, etc.)
    df_1 = df_1.merge(meta_grouped_Genetic_features, on='DepMap_ID', how='inner')
    #df_1.drop_duplicates(inplace=True) (I can't do this because there are lists with mutated_genes)
    df_1.drop_duplicates(subset=['DepMap_ID', 'COSMIC_ID'],inplace=True)
    
     
    # 2. Apply user-defined filters (e.g., TCGA type, lineage)
    for col, condition in dict_filters.items():
        if col in df_1.columns:
            df_1 = df_1[condition(df_1[col])]
    
    # 3. Mutation filtering (via grouped genetic features)
    mask = grouped_Genetic_features["Mutated Genes"].apply(
        lambda x: any(gene in x for gene in list_mutation)
    )
    grouped_Genetic_features_filtered = grouped_Genetic_features.copy()
    df_2 = df_1.copy()

    if drop_mutation:
        # Keep only non-mutated samples
        grouped_Genetic_features_filtered = grouped_Genetic_features[~mask].copy()
        df_2 = df_1.merge(grouped_Genetic_features_filtered, on='COSMIC_ID')

    if keep_only_mutation:
        # Keep only mutated samples
        grouped_Genetic_features_filtered = grouped_Genetic_features[mask].copy()
        df_2 = df_1.merge(grouped_Genetic_features_filtered, on='COSMIC_ID')
    
    # 4. Final cleanup
    try:
        df_2.drop_duplicates(inplace=True)
    except:
        pass
    
    # 5. Filter expression_df to keep only relevant samples
    depmap_id_list = df_2['DepMap_ID'].unique().tolist()
    expression_df = df_1[df_1['DepMap_ID'].isin(depmap_id_list)][['DepMap_ID'] + expr_keep]
    
    return df_2, expression_df

def load_prism_data():
    """
    Load PRISM drug response data.

    Args:
        None

    Returns:
        pd.DataFrame: Merged PRISM dataset.
    """
    try:
        lfc = pd.read_csv("../data/Cell-lines/Repurposing_Public_24Q2_LFC_COLLAPSED.csv")
        cell_meta = pd.read_csv("../data/Cell-lines/Repurposing_Public_24Q2_Cell_Line_Meta_Data.csv")
        drug_meta = pd.read_csv("../data/Cell-lines/Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv")
    except Exception as e:
        print(f"Error loading PRISM data: {e}")
        return pd.DataFrame()

    # Keep only relevant columns from cell and drug metadata
    cell_meta = cell_meta[['row_id', 'ccle_name', 'depmap_id']]
    drug_meta = drug_meta[['IDs', 'Drug.Name', 'repurposing_target', 'MOA']]

    # Clean the 'IDs' column by removing the 'BRD:' prefix to match with 'broad_id' in lfc
    drug_meta['IDs'] = drug_meta['IDs'].str.replace('^BRD:', '', regex=True)

    # Merge drug response data with cell line metadata to get cell line names
    df = lfc.merge(cell_meta, on='row_id', how='left')

    # Merge with drug metadata to get drug names, targets, and MOA
    df = df.merge(drug_meta, left_on='broad_id', right_on='IDs', how='left')

    # Drop entries with missing drug IDs (unmatched merges)
    df = df.dropna(subset=['IDs'])

    # Add a constant column indicating the dataset source
    df['DATASET'] = 'PRISM'

    # Rename depmap_id to match GDSC2 naming convention
    df = df.rename(columns={'depmap_id': 'DepMap_ID'})

    # Extract cell line name and cancer origin from the 'ccle_name' field
    df[['CELL_LINE_NAME', 'Origin']] = df['ccle_name'].str.split('_', n=1, expand=True)

    # Remove duplicate rows
    df.drop_duplicates(inplace=True)

    # Rename columns to harmonize with the GDSC2 format
    df_prism = df.rename(columns={
        'Drug.Name': 'DRUG_NAME',
        'repurposing_target': 'PUTATIVE_TARGET',
        'MOA': 'PATHWAY_NAME',
        #'LFC': 'AUC',  # In PRISM, LFC is inversely related to sensitivity: more negative = more sensitive
        'COSMICID': 'COSMIC_ID'  # This column may not be available unless merged externally
    })

    # Reorder and select final columns for the standardized output
    df_prism = df_prism[[
        'DATASET', 'ccle_name', 'CELL_LINE_NAME', 'Origin', 'DepMap_ID', 'IDs', 'DRUG_NAME',
        'PUTATIVE_TARGET', 'PATHWAY_NAME', 'LFC'
    ]]
    return df_prism

def load_data_sensitivity( dict_filters, list_mutation, target_genes, 
              gdsc_n='gdsc2', drop_mutation=False, keep_only_mutation=False, expr_keep=[]):
    """
    Load and merge omics, mutation, and drug sensitivity data according to input parameters.

    Args:
        dict_filters (dict): Filtering conditions as {column_name: lambda}.
        list_mutation (list): Genes to use for mutation-based filtering.
        target_genes (list): Drug target gene list.
        gdsc_n (str): Which dataset to use for drug response ('gdsc2' or 'prism').
        drop_mutation (bool): If True, keep only non-mutated samples.
        keep_only_mutation (bool): If True, keep only mutated samples.
        expr_keep (list): List of expression columns to retain.

    Returns:
        pd.DataFrame: Merged and filtered dataset ready for analysis.
    """
    # 1. Load drug response datasets
    if gdsc_n == 'gdsc2':
        try : 
            gdsc = pd.read_csv("../data/Cell-lines/GDSC2_fitted_dose_response_27Oct23.csv", sep = ';')
        except Exception as e:
            print(f"Error loading GDSC2 data: {e}")
            return pd.DataFrame()  # Return empty DataFrame on error
    elif gdsc_n == 'prism':
        gdsc = load_prism_data()

    # 2. load data
    df_2,_ = load_merge_data(dict_filters, list_mutation, drop_mutation=drop_mutation, keep_only_mutation=keep_only_mutation, expr_keep=expr_keep)

    # 3. Filter drug response data by target genes (if any)
    if target_genes:
        def target_match(x):
            if pd.isna(x):
                return False
            return any(gene in x for gene in target_genes)
        gdsc_filtered = gdsc[gdsc['PUTATIVE_TARGET'].apply(target_match)].copy()
    else:
        gdsc_filtered = gdsc.copy()

    # 4. Merge with drug data on appropriate key
    if gdsc_n == 'prism':
        merged_df = pd.merge(df_2, gdsc_filtered, on='DepMap_ID', how='inner')
    else:
        merged_df = pd.merge(df_2, gdsc_filtered, on='COSMIC_ID', how='inner')

    # 5. Cleanup and sort
    try:
        merged_df.drop_duplicates(inplace=True)
    except:
        pass
    try:
        return merged_df.sort_values(by='MTUS1')
    except:
        return merged_df
    

def correl_expression_measure_sensitivity(df, measure_col, expr_gene='MTUS1'):
    """
    Computes the Pearson correlation between gene expression (e.g., MTUS1) and a drug response measure.

    Args:
        df (pd.DataFrame): DataFrame containing expression and measure columns.
        measure_col (str): Column name for drug sensitivity (for gdsc2 : 'AUC' or 'LN_IC50', for prism : 'LFC').
        expr_gene (str): Column name for gene expression (default is 'MTUS1').

    Returns:
        tuple: (correlation coefficient, p-value)
    """
    # Clean missing or improperly formatted values
    df[measure_col] = df[measure_col].astype(str).str.replace(',', '.').astype(float)

    sub_df = df[[expr_gene, measure_col]].dropna().drop_duplicates(subset=[expr_gene, measure_col])

    # Compute Pearson correlation
    corr, pval = pearsonr(sub_df[expr_gene], sub_df[measure_col])
    return corr, pval


def correlations_expression_measure_all_drugs_sensitivity(df, measure_col, expr_gene='MTUS1'):
    """
    Computes Pearson correlation between gene/protein expression and drug response for each drug.

    Args:
        df (pd.DataFrame): DataFrame with expression, response, and drug columns.
        measure_col (str): Column name for drug sensitivity (for gdsc2 : 'AUC' or 'LN_IC50', for prism : 'LFC').
        expr_gene (str): Expression column (e.g. 'MTUS1').


    Returns:
        pd.DataFrame: With columns for drug name, correlation, p-value, FDRand sample size.
    """
    results = []

    for drug in df['DRUG_NAME'].unique():
        sub = df[df['DRUG_NAME'] == drug][[expr_gene, measure_col, 'PUTATIVE_TARGET', 'PATHWAY_NAME']].dropna()
        sub = sub.drop_duplicates(subset=[expr_gene, measure_col])
        sub[measure_col] = sub[measure_col].astype(str).str.replace(',', '.')
        sub[measure_col] = pd.to_numeric(sub[measure_col], errors='coerce')
        sub = sub.dropna()

        if len(sub) >= 3:
            x = sub[expr_gene].to_numpy()
            y = sub[measure_col].to_numpy()
            r, p = pearsonr(sub[expr_gene], sub[measure_col])

            # Optionally, perform permutation test for empirical p-value
            #n_perm = 10_000
            #perm_r = np.empty(n_perm)

            #rng = np.random.default_rng(seed=42)      
            #for i in range(n_perm):
            #    y_perm = rng.permutation(y)           # permute AUC
            #    perm_r[i], _ = pearsonr(x, y_perm)    # new correlation

            # empirical two-tailed p-value
            #p_perm = (np.sum(np.abs(perm_r) >= abs(r)) + 1) / (n_perm + 1)

            results.append({
                'DRUG_NAME': drug,
                'Putative_Target': sub['PUTATIVE_TARGET'].iloc[0] if 'PUTATIVE_TARGET' in sub.columns else None,
                'PATHWAY_NAME': sub['PATHWAY_NAME'].iloc[0] if 'PATHWAY_NAME' in sub.columns else None,
                'correlation': r,
                'p_value': p,
            #    'p_value_perm': p_perm,
                'n': len(sub)
            })

    res_df = pd.DataFrame(results)

    # Apply Benjamini-Hochberg correction
    if not res_df.empty:
        _, fdr_corrected_pvals, _, _ = multipletests(res_df['p_value'], method='fdr_bh')
        res_df['p_value_FDR'] = fdr_corrected_pvals
    try: 
        return res_df.sort_values('p_value')
    except:
        return res_df


def plot_expr_vs_measure_sensitivity(
    df,
    drug_name,
    measure_col,
    expr_gene='MTUS1',
    ):
    """
    Plots the correlation between gene expression and drug response for a specific drug.

    Args:
        df (pd.DataFrame): Must contain drug name, expression, and response columns.
        drug_name (str): Drug to focus on.
        measure_col (str): Column name for drug sensitivity (for gdsc2 : 'AUC' or 'LN_IC50', for prism : 'LFC').
        expr_gene (str): Column name for gene expression.
        type_data (str): Either 'TPMLogp1' (transcriptomic). Used for x-axis labeling.

    Returns:
        tuple: (correlation coefficient, p-value)
    """

    plt.figure(figsize=(6, 5))

    # Filter and clean data
    sub = df[df['DRUG_NAME'] == drug_name].dropna(subset=[expr_gene, measure_col]).drop_duplicates(subset=[expr_gene, measure_col])

    # Conversion ici (si jamais des virgules traînent)
    sub[measure_col] = sub[measure_col].astype(str).str.replace(',', '.')
    sub[measure_col] = pd.to_numeric(sub[measure_col], errors='coerce')

    if len(sub) < 3:
        print(f"Not enough data points for {drug_name}")
        return None, None

    # Compute Pearson correlation
    r, p = pearsonr(sub[expr_gene], sub[measure_col])

    # Plot scatter and regression line
    sns.regplot(data=sub, x=expr_gene, y=measure_col, scatter_kws={'alpha': 0.6})

    # Highlight TNBC samples if available
    if 'ModelSubtypeFeatures' in sub.columns:
        sns.scatterplot(
            data=sub[sub['ModelSubtypeFeatures'].str.contains('TNBC', na=False)],
            x=expr_gene, y=measure_col,
            label='TNBC', color='r', alpha=0.6
        )

    x_label = f"{expr_gene} (log2(TPM+1))"

    # Annotate the plot
    plt.title(f"{drug_name}\nCorrelation {expr_gene} vs {measure_col}\nr = {r:.2f}, p = {p:.1e}")
    plt.xlabel(x_label)
    plt.ylabel(measure_col)
    plt.axhline(0, linestyle='--', color='gray', linewidth=0.5)
    plt.tight_layout()
    plt.show()

    return r, p

#################
def _ensure_depmap_id_columns(df):
    """
    Ensure 'DepMap_ID' is a column of type str (not just an index).
    
    Args:
        df (pd.DataFrame): DataFrame to ensure 'DepMap_ID' column.

    Returns:
        pd.DataFrame: DataFrame with 'DepMap_ID' as a column of type str.
    """
    if df.index.name == "DepMap_ID" or "DepMap_ID" not in df.columns:
        df = df.reset_index()
    if "DepMap_ID" not in df.columns:
        raise ValueError("DepMap_ID not found as column or index.")
    df = df.copy()
    df["DepMap_ID"] = df["DepMap_ID"].astype(str)
    return df


def correlate_gene_effect_with_expression(effect_df, expression_df, effect_gene='WEE1', expr_gene='MTUS1'):
    """
    Correlates the dependency on a gene (CRISPRGeneEffect) with the expression of another gene (e.g. MTUS1)
    using Pearson correlation.

    Args:
        effect_df (pd.DataFrame): DataFrame containing gene effect data.
        expression_df (pd.DataFrame): DataFrame containing gene expression data.
        effect_gene (str): Name of the gene for which to assess knockout effect.
        expr_gene (str): Name of the gene for which to assess expression.

    Returns:
        tuple: Pearson correlation coefficient and p-value.
    """
    effect_df = _ensure_depmap_id_columns(effect_df)
    expression_df = _ensure_depmap_id_columns(expression_df)

    # Ensure DepMap_ID is of type string
    effect_df['DepMap_ID'] = effect_df['DepMap_ID'].astype(str)
    expression_df['DepMap_ID'] = expression_df['DepMap_ID'].astype(str)

    # Merge and clean data
    sub = expression_df[['DepMap_ID', expr_gene]].merge(effect_df[['DepMap_ID', effect_gene]], on='DepMap_ID')
    sub = sub.dropna()
    corr, pval = pearsonr(sub[expr_gene], sub[effect_gene])
    return corr, pval

def correlate_selected_genes_with_expression(effect_df, expression_df, expr_gene='MTUS1', genes=[]):
    """
    Correlates the expression of a selected gene with the effect of other genes. 

    Args:
        effect_df (pd.DataFrame): DataFrame containing gene effect data.
        expression_df (pd.DataFrame): DataFrame containing gene expression data.
        expr_gene (str): Name of the gene for which to assess expression.
        genes (list): List of genes to correlate with the expression of expr_gene.

    Returns:
        pd.DataFrame: DataFrame containing Pearson correlation coefficients, p-values, FDR (Benjamini-Hochberg) corrected p-values.

    """
    effect_df = _ensure_depmap_id_columns(effect_df)
    expression_df = _ensure_depmap_id_columns(expression_df)

    if expr_gene not in expression_df.columns:
        raise KeyError(f"{expr_gene} not in expression_df.")

    if genes is None or len(genes) == 0:
        # assume effect_df has DepMap_ID + gene columns
        genes = [c for c in effect_df.columns if c not in {"DepMap_ID"}]
        if expr_gene in genes:
            genes.remove(expr_gene)

    rows = []
    expr = expression_df[["DepMap_ID", expr_gene]].dropna()
    for g in genes:
        if g not in effect_df.columns:
            continue
        merged = expr.merge(effect_df[["DepMap_ID", g]], on="DepMap_ID", how="inner").dropna()
        if len(merged) >= 3:
            r, p = pearsonr(merged[expr_gene].to_numpy(), merged[g].to_numpy())
            rows.append(
                {
                    "Gene": g,
                    "Correlation": float(r),
                    "P-value": float(p),
                    "n": int(len(merged)),
                    "min_score": float(merged[g].min()),
                    "max_score": float(merged[g].max()),
                }
            )

    res = pd.DataFrame(rows)
    if res.empty:
        return res

    # Benjamini–Hochberg FDR
    res["P-value_FDR"] = multipletests(res["P-value"].to_numpy(), method="fdr_bh")[1]
    res_df = res.reset_index(drop=True)
    try: 
        return res_df.sort_values('P-value') 
    except: 
        return res_df

def plot_expr_vs_effect(effect_df, expression_df, expr_gene='MTUS1', effect_gene='MYC'):
    """
    Plots the expression of a gene against the effect of another gene.

    Args:
        effect_df (pd.DataFrame): DataFrame containing gene effect data.
        expression_df (pd.DataFrame): DataFrame containing gene expression data.
        expr_gene (str): Name of the gene for which to assess expression.
        effect_gene (str): Name of the gene for which to assess effect.

    Returns:
        None
    """

    plt.figure(figsize=(6, 5))

    sub = pd.merge(
            expression_df[['DepMap_ID', expr_gene]],
            effect_df[['DepMap_ID', effect_gene]],
            on='DepMap_ID'
        ).dropna()

    # Check if there are enough points
    if len(sub) < 3:
        print(f"Not enough points")
        return

    # Calculate correlation
    r, p = pearsonr(sub[expr_gene], sub[effect_gene])
    sns.regplot(data=sub, x=expr_gene, y=effect_gene, scatter_kws={'alpha':0.6})
    #sns.regplot(data=sub[sub['ModelSubtypeFeatures'].str.contains('TNBC', na=False)],x=expr_gene, y=effect_gene, scatter_kws={'alpha': 0.6}, color='r', label='TNBC')
    #sns.scatterplot(data=sub[sub['ModelSubtypeFeatures'].str.contains('TNBC', na=False)],x=expr_gene, y=effect_gene, label='TNBC', color='r', alpha=0.6)
    

    plt.title(f"Correlation {expr_gene} vs {effect_gene}\nr = {r:.2f}, p = {p:.1e}")
    plt.xlabel(f"{expr_gene} (log2(TPM+1))")
    plt.ylabel(f"{effect_gene}effect (CRISPRGeneEffect)")
    plt.axhline(0, linestyle='--', color='gray', linewidth=0.5)
    plt.tight_layout()
    plt.show()

def plot_expr_vs_expression( expression_df,  expr_gene='MTUS1', expression_col='MYC'):
    """
    Plots the expression of one gene against another.
    
    Args:
        expression_df (pd.DataFrame): DataFrame containing gene expression data.
        expr_gene (str): Name of the gene for which to assess expression.
        expression_col (str): Name of the gene for which to assess expression.

    Returns:
        None
    """
    plt.figure(figsize=(6, 5))

    sub = pd.merge(
            expression_df[['DepMap_ID', expr_gene]],
            expression_df[['DepMap_ID', expression_col]],
            on='DepMap_ID'
        ).dropna()
    
    # Check if there are enough points
    if len(sub) < 3:
        print(f"Not enough points")
        return

    # Calculate correlation
    r, p = pearsonr(sub[expr_gene], sub[expression_col])
    sns.regplot(data=sub, x=expr_gene, y=expression_col, scatter_kws={'alpha':0.6})
    #sns.regplot(data=sub[sub['ModelSubtypeFeatures'].str.contains('TNBC', na=False)],x=expr_gene, y=expression_col, scatter_kws={'alpha': 0.6}, color='r', label='TNBC')
    #sns.scatterplot(data=sub[sub['ModelSubtypeFeatures'].str.contains('TNBC', na=False)],x=expr_gene, y=expression_col, label='TNBC', color='r', alpha=0.6)
    
    

    # Annotation
    plt.title(f"Correlation {expr_gene} vs {expression_col}\nr = {r:.2f}, p = {p:.1e}")
    plt.xlabel(f"{expr_gene} (log2(TPM+1))")
    plt.ylabel(f"{expression_col} (log2(TPM+1))")
    plt.axhline(1, linestyle='--', color='gray', linewidth=0.5)
    plt.tight_layout()
    plt.show()