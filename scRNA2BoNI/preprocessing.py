#!/usr/bin/env python3
# -*- coding: utf-8 -*

try:
    import csv
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy import stats
    import json
    from functools import partial
    import pkg_resources
    import clyngor


    from .utils import  read_file, load_data, get_classes, save_to_file

except ImportError as E:
    print(E)
    print('Verify your dependancies.')
    exit()

def get_extrem_expr_values(mtx: pd.DataFrame, readouts:list, cell_types_selected:list)-> float:
    """Get the maximum value of a subset of a dataframe, where subset columns are given in entry.

    Parameters
    ----------
    mtx : pd.DataFrame
        Dataframe to treat
    readouts : list
        Columns line of the subset

    Returns
    -------
    float
        Maximum value
    """
    # TODO : Change function description !!!
    reduced_mtx = mtx[mtx['clusterUmap'].isin(cell_types_selected)] 
    reduced_mtx = reduced_mtx[readouts]
    max_ = reduced_mtx.max(numeric_only=True).max()
    min_ = reduced_mtx.min(numeric_only=True).min()
    return max_, min_



def reduce_and_binarize_matrix(mtx, cell_types_selected, readouts, inputs, intermediates, annotations_len):
    def binarize(x):
        if x < 2:
            return 0
        else:
            return 1

    def normalize(x, max_expr_value, min_expr_value):
        return (x-min_expr_value)/(max_expr_value - min_expr_value)

    # Maybe useles, because file containing the genes are already tested if present in matrix
    # If yes, remove all ?
    i_o_genes = list(set(mtx.columns.tolist()).intersection(set(inputs+intermediates)))
    readouts = list(set(mtx.columns.tolist()).intersection(set(readouts)))
    max_expr_value, min_expr_value = get_extrem_expr_values(mtx, readouts, cell_types_selected)
    
    


    annotations = mtx.columns[:annotations_len].to_list()
    
    reduced_expr_mtx = mtx.loc[mtx.clusterUmap.isin(cell_types_selected)]
    reduced_expr_mtx = reduced_expr_mtx[annotations+i_o_genes+readouts]

    bin_reduced_expr_mtx = reduced_expr_mtx.copy()
    for gene in bin_reduced_expr_mtx.columns[annotations_len:]:
        if gene in readouts:
            bin_reduced_expr_mtx[gene] = bin_reduced_expr_mtx[gene].apply(partial(normalize,max_expr_value=max_expr_value, min_expr_value=min_expr_value))
        elif gene in i_o_genes:
            bin_reduced_expr_mtx[gene] = bin_reduced_expr_mtx[gene].apply(binarize)

    return bin_reduced_expr_mtx, i_o_genes, readouts

def encode(to_encode:list)-> dict:
    encoded = {
        'decode':dict(),
        'encode':dict()
    }
    for i in range(len(to_encode)):
        encoded['decode'][i+1] = to_encode[i]
        encoded['encode'][to_encode[i]] = i+1
    return encoded



def to_lp_file(to_transform:list, predicate_name:str):
    instance = str()
    for item_to_transform in to_transform:
        items = item_to_transform.split(',')
        transformed = str()
        for i in range(len(items)):
            transformed += f'"{items[i]}",'
        transformed = transformed[:-1] # Remove last useless ','
        instance += f'{predicate_name}({transformed}).\n'        
    return instance

def get_parents(encoding, out_dir):
    print('GET PARENTS')
    instance_path = pkg_resources.resource_filename(__name__, 'data/processing/get_parents.lp')
    intermediates_path = f'{out_dir}/intermediates_instance.lp'
    programs = [instance_path, intermediates_path]
    print(programs)
    answers = clyngor.solve(programs, inline=encoding)
    print(answers)
    for answer in answers.by_predicate:
        parents_answer = list(answer['parent'])
        parents = list()
        for tuple in parents_answer:
            parents.append(','.join([elem.strip('"') for elem in tuple]))
        # parents = [x.replace('"','') for x in parents]
        print(parents)
        return parents



def create_expression_instance(matrix, readouts, annotations_len):
    instance = str()
    cells_list = dict()
    for class_ in get_classes(matrix):
        cells_list[class_] = list()

    for i in list(matrix.index.values):
        cell = i
        class_ = matrix.at[i, 'clusterUmap']
        cells_list[class_].append(cell)
        for gene in matrix.columns[annotations_len:]:
            value = matrix.at[i, gene]
            if gene in readouts:
                pass
                # TODO: deal both classes for the readouts. Add the class in the predicate maybe ? 2 files (one for each class) ?
                # instance += f'r("{cell}","{gene}","{value}").\n'
            else:
                instance += f'pert("{cell}","{gene}",{value},"{class_}").'
        instance += '\n'
    return instance, cells_list


def run_preprocessing(config):

    out_dir = config['output_dir']
    cell_types_selected = config['class_types']



    # matrix = load_data(f"{out_dir}/bin_reduced_matrix.csv")
    matrix = load_data(f"{out_dir}/pkn_gene_reduced_expr_mtx.csv", index_name='Name')
    readouts = read_file(f"{out_dir}/no_successors_in_the_matrix.txt")
    inputs = read_file(f"{out_dir}/no_predecessors_in_the_matrix.txt")
    intermediates = read_file(f"{out_dir}/intermediates_in_the_matrix.txt")
    annotations_len = config['annotation_len']
    bin_reduced_matrix, io_genes_intersectPKN_matrix, readouts_genes_intersectPKN_matrix = reduce_and_binarize_matrix(matrix, cell_types_selected, readouts, inputs, intermediates, annotations_len)
    # bin_reduced_matrix = bin_reduced_matrix.set_index('Name')
    bin_reduced_matrix.to_csv(f'{out_dir}/bin_reduced_matrix.csv', index=True)
    # list_to_file(io_genes_intersectPKN_matrix, f'{out_dir}/io_genes_intersectPKN_matrix.txt')
    # list_to_file(readouts_genes_intersectPKN_matrix, f'{out_dir}/readouts_genes_intersectPKN_matrix.txt')

    # genes_list = matrix.columns[annotations_len:]
    # genes_hash_map = encode(genes_list)
    # json.dump(genes_hash_map, open(f"{out_dir}/genes_hash_map.json", "w"), sort_keys=True, indent=4)
    # cells_list = list(matrix.index.values)
    # cells_hash_map = encode(cells_list)
    # json.dump(cells_hash_map, open(f"{out_dir}/cells_hash_map.json", "w"), sort_keys=True, indent=4)
    # classes_list = config['class_types']
    # classes_hash_map = encode(classes_list)
    # json.dump(classes_hash_map, open(f"{out_dir}/classes_hash_map.json", "w"), sort_keys=True, indent=4)


    inputs_instance = to_lp_file(inputs, 'input')
    save_to_file(inputs_instance, f'{out_dir}/inputs_instance.lp')
    intermediates_instance = to_lp_file(intermediates, 'intermediate')
    save_to_file(intermediates_instance, f'{out_dir}/intermediates_instance.lp')
    
    pkn_encoding = open(f'{out_dir}/pkn.lp').read()
    parents = get_parents(pkn_encoding, out_dir)
    parents_instance = to_lp_file(parents,'parent')

    expression_instance, cells_list = create_expression_instance(bin_reduced_matrix, readouts, annotations_len)
    
    
    save_to_file(parents_instance, f'{out_dir}/parents_instance.lp')
    save_to_file(expression_instance, f'{out_dir}/expression_instance.lp')
    with open(f'{out_dir}/cells_list_by_classes.json', 'w') as f:
        json.dump(cells_list, f)

