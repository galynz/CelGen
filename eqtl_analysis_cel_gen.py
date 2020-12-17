# Imports
import argparse
import numpy as np
import pandas as pd
import statsmodels.api as sm

from tree_tools import generate_tree

GENE_CHR = "chromosome_name"
GENE_START = "start_position"
GENE_END = "end_position"
GENE_NAME_INDEX = 0
SNP_NAME = "marker"
SNP_CHR = "chromosome"
SNP_POS = "bp"
SNP_NAME_INDEX = 0


def filter_cis(snp_id, snp_positions, gene_expression, gene_positions, cis_add=5000000):
    try:
        snp_chr = snp_positions.loc[snp_id, SNP_CHR]
        snp_pos = snp_positions.loc[snp_id, SNP_POS]
    except KeyError:
        return pd.DataFrame()

    genes = gene_positions[(gene_positions[GENE_CHR] == snp_chr) & (gene_positions[GENE_START] > snp_pos - cis_add) & (
                gene_positions[GENE_END] < snp_pos + cis_add)]
    return gene_expression.loc[genes.index.intersection(gene_expression.index)]


def filter_snp(snp_series, cur_gene_exp, cell_fractions, maf=0.1):
    not_null_snps = snp_series.dropna().index
    snp_series = snp_series[not_null_snps]
    if snp_series.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    snp_counts = snp_series.value_counts(normalize=True)
    if len(snp_counts) <= 1 or snp_counts.iloc[1] < maf:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    return snp_series[not_null_snps], cur_gene_exp[not_null_snps], cell_fractions.loc[not_null_snps]


def test_gene_snp(cur_gene_exp, cur_genotype, cell_fractions, nodes_dict):
    cur_df = pd.concat([cur_gene_exp, cur_genotype], axis=1, ignore_index=True)
    cur_df.columns = ["gene_exp", "genotype"]
    cur_df = cur_df.join(cell_fractions)
    genotype_cell_prop = np.multiply(cell_fractions.T, cur_genotype.values).T
    genotype_cell_prop.columns = genotype_cell_prop.columns.map(lambda x: x+"_genotype")
    cur_df = cur_df.join(genotype_cell_prop)

    res_restricted = sm.OLS(cur_df["gene_exp"], sm.add_constant(cur_df[["genotype"]], has_constant='add')).fit()

    full_res = {"gene_id":                 cur_gene_exp.name,
                "snp_id":                  cur_genotype.name,
                "generic_log_likelihood":  res_restricted.llf,
                "generic_ssr":             res_restricted.ssr,
                "generic_ess":             res_restricted.ess,
                "generic_centered_tss":    res_restricted.centered_tss,
                "generic_r2":              res_restricted.rsquared,
                "generic_r2_adj":          res_restricted.rsquared_adj,
                "generic_p_value":         -np.log10(res_restricted.f_pvalue),
                "generic_aic":             res_restricted.aic,
                "generic_bic":             res_restricted.bic,
                "generic_genotype_coeff":  res_restricted.params["genotype"],
                "generic_intercept_coeff": res_restricted.params["const"]
                }

    for node in cell_fractions.columns:
        cells = nodes_dict[node].get_subtree_names()
        cur_df[node+"_root_genotype"] = cur_df[[cell+"_genotype" for cell in cells]].sum(axis=1)
        res_full = sm.OLS(cur_df["gene_exp"],
                          sm.add_constant(cur_df[["genotype", node + "_root_genotype"]], has_constant='add')).fit()
        f_test = -np.log10(res_full.compare_f_test(res_restricted)[1])
        full_res.update({node+"_restricted_full_f_test": f_test,
                         node+"_full_log_likelihood":    res_full.llf,
                         node+"_ssr":                    res_full.ssr,
                         node+"_ess":                    res_full.ess,
                         node+"_centered_tss":           res_full.centered_tss,
                         node+"_r2":                     res_full.rsquared,
                         node+"_r2_adj":                 res_full.rsquared_adj,
                         node+"_p_value":                -np.log10(res_full.f_pvalue),
                         node+"_full_aic":               res_full.aic,
                         node+"_full_bic":               res_full.bic,
                         node+"_genotype_coeff":         res_full.params["genotype"],
                         node+"_genotype_cell_coeff":    res_full.params[node+"_root_genotype"],
                         node+"_intercept_coeff":        res_full.params["const"]})

    return pd.DataFrame(full_res, index=[cur_gene_exp.name])


def main():
    # Read inputs
    gene_expression = pd.read_csv(args.gene_expression_path, index_col=0)
    genotype = pd.read_csv(args.genotype_path, index_col=0)
    cell_fractions = pd.read_csv(args.cell_fractions_path, index_col=0)
    cell_fractions.columns = cell_fractions.columns.str.upper()
    cell_fractions.index = cell_fractions.index.astype(str)

    shared_columns = gene_expression.columns.intersection(genotype.columns).intersection(
        cell_fractions.index.astype(str))
    gene_expression = gene_expression[shared_columns]
    genotype = genotype[shared_columns]
    cell_fractions = cell_fractions.loc[shared_columns]

    # Read tree's structure
    tree_df = pd.read_csv(args.tree_path, index_col=0)
    root_name = tree_df.loc[0, 'name']
    tree, nodes_dict = generate_tree(args.tree_path, root_name)

    gene_positions = pd.read_csv(args.gene_positions, index_col=GENE_NAME_INDEX)
    snp_positions = pd.read_csv(args.snp_positions, index_col=SNP_NAME_INDEX)

    gene_expression = gene_expression.loc[gene_positions.index.intersection(gene_expression.index)]
    gene_positions = gene_positions.loc[gene_positions.index.intersection(gene_expression.index)]

    all_res = []
    for snp_id, cur_genotype in genotype.iterrows():
        cis_gene_expression = filter_cis(snp_id, snp_positions, gene_expression, gene_positions)
        for gene_id, cur_gene_expression in cis_gene_expression.iterrows():
            cur_genotype, cur_gene_expression, cur_cell_fractions = filter_snp(cur_genotype,
                                                                               cur_gene_expression,
                                                                               cell_fractions)

            if cur_gene_expression.empty or np.isclose(cur_gene_expression.var(), 0):
                continue

            all_res.append(test_gene_snp(cur_gene_expression, cur_genotype, cur_cell_fractions, nodes_dict))

    all_res_df = pd.concat(all_res)
    all_res_df.to_csv(args.output_path, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run Cel-Gen eQTL analysis')
    parser.add_argument('output_path', help='where to save the output file')
    parser.add_argument('--genotype', dest='genotype_path', help="genotype file path")
    parser.add_argument('--gene_expression', dest='gene_expression_path', help="gene_expression file path")
    parser.add_argument('--cell_fractions', dest='cell_fractions_path', help="cell fractions file path")
    parser.add_argument('--tree_structure', dest='tree_path', help="tree structure file path")
    parser.add_argument('--gene_positions', dest="gene_positions",
                        help="a csv file with chromosome and gene positions (start and end) for each gene")
    parser.add_argument('--snp_positions', dest="snp_positions",
                        help="a csv file with chromosome and snp position for each marker")

    args = parser.parse_args()
    main()
