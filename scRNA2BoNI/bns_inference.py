#!/usr/bin/env python3
# -*- coding: utf-8 -*

try:
    import csv
    import pandas as pd
    import numpy as np 
    import re
    import json
    from argparse import Namespace
    import pydotplus
    import os
    import networkx as nx
    import matplotlib.pyplot as plt

    from .utils import load_data, read_file, get_cell_class
    # from caspo.caspo.console.handlers
    from caspo.console.handlers import learn_handler, classify_handler, visualize_handler


except ImportError as E:
    print(E)
    print('Verify your dependancies.')
    exit()


def find_last_index(search_list, starting_item):
    starting_item = '^'+starting_item
    indices = [ i for i, word in enumerate(search_list) if re.search(starting_item, word) ]
    return indices[-1]

def get_value(atom):
    return re.findall('\"([A-Za-z0-9\.\-]*)\"', atom)


# def get_value(atom):
#     print('atom:', atom)
#     print(re.findall('\"([A-Z0-9\.]*)\"', atom))
#     return re.findall('([0-9]+)', atom)


def load_answer_set(filename: str) -> set:
    sel_genes = list()
    affinities = list()
    file_list = read_file(filename)
    index = find_last_index(file_list, 'Optimization:')-1 # -1 to have the answer-set of the last optimization
    answer_set_raw = file_list[index]
    answer_set_list = answer_set_raw.split(' ')
    for atom in answer_set_list:
        if re.search('^selgene', atom):
            encoded_gene = get_value(atom)[0]
            # gene = genes_hash_map['decode'][encoded_gene]
            gene = encoded_gene
            sel_genes.append(gene)
        elif re.search('^affinity', atom):
            cur_cells = get_value(atom)
            encoded_cell_i = cur_cells[0]
            # cell_i = cells_hash_map['decode'][encoded_cell_i]
            cell_i = encoded_cell_i
            encoded_cell_j = cur_cells[1]
            # cell_j = cells_hash_map['decode'][encoded_cell_j]
            cell_j = encoded_cell_j
            affinities.append([cell_i, cell_j])
            # if cell_i in affinities.keys():
            #     affinities[cell_i].append(cell_j)
            # else:
            #     affinities[cell_i] = [cell_j]
            

    return sel_genes, affinities

def get_classes_cells(matrix, classes):
    # Not necessary, see config dict
    classes_cells = dict()
    for class_ in classes:
        classes_cells[class_] = matrix[matrix['clusterUmap']==class_].index.tolist()
    return classes_cells


def get_expr_vector(cell:str, genes:list, matrix:pd.DataFrame)-> list:
    expr_vector = list()
    for gene in genes:
        if gene in matrix.columns.to_list(): # TO REMOVE
            expr_vector.append(matrix.at[cell, gene])
    return expr_vector

def get_redundancies(studied_cell:str, selgenes:list, class_:str, matrix:pd.DataFrame)-> list:
    redundancies = list()
    redundancies.append(studied_cell)
    cells = matrix.index.tolist()
    studied_cell_vector = get_expr_vector(studied_cell, selgenes, matrix)
    for cell in cells: 
        if cell != studied_cell:

            vector = get_expr_vector(cell, selgenes, matrix)
            if vector == studied_cell_vector:
                redundancies.append(cell)
    
    return redundancies

def redundancies_calculation(affinities:dict, selgenes:list, class_names:list, matrix:pd.DataFrame)-> dict:
    redundancy_cells = {
        class_names[0]: {},
        class_names[1]: {}
    }
    redundancy_vectors = {}
    c1_class_name = class_names[0]
    c1_submatrix = matrix[matrix['clusterUmap']==c1_class_name]

    c2_class_name = class_names[1]
    c2_submatrix = matrix[matrix['clusterUmap']==c2_class_name]
    for i in range(len(affinities)):
        c1 = affinities[i][0]
        c2 = affinities[i][1]
        redundancy_vectors[i] = list()
        c1_redundancies = get_redundancies(c1, selgenes, c1_class_name, c1_submatrix)
        redundancy_cells[c1_class_name][i] = c1_redundancies
        c2_redundancies = get_redundancies(c2, selgenes, c2_class_name, c2_submatrix)
        redundancy_cells[c2_class_name][i] = c2_redundancies
        
        for c1_red in c1_redundancies:
            for c2_red in c2_redundancies:
                redundancy_vectors[i].append((c1_red, c2_red))
    
    sum_ = 0
    for values in redundancy_cells[c1_class_name].values():
        sum_ += len(values)
    redundancy_cells[c1_class_name]['nb'] = sum_
    sum_ = 0
    for values in redundancy_cells[c2_class_name].values():
        sum_ += len(values)
    redundancy_cells[c2_class_name]['nb'] = sum_
    return redundancy_cells, redundancy_vectors


    # for aff_class, aff_cells in affinities.items():
    #     for aff_cell in aff_cells:
    #         aff_expr = get_expr_vector(aff_cell, selgenes, matrix)
    #         for red_cell in class_cells[aff_class]:
    #             red_expr = get_expr_vector(red_cell, selgenes, matrix)
    #             if aff_expr == red_expr:
    #                 redundancies[aff_class]['redundancies'].append(red_cell)
    #     redundancies[aff_class]['nb'] = len(redundancies[aff_class]['redundancies'])
    # return redundancies

# def duplications_calculation(affinities:dict, selgenes:list, matrix:pd.DataFrame) -> dict:
#     duplications = dict()
#     aff_i_list = list(affinities.keys())
#     for aff in aff_i_list:
#         duplications[aff] = {
#             'to_treat':True,
#             'duplications': [aff]
#         }
#     for i in aff_i_list[:-1]:
#         for j in aff_i_list[1:]:
#             if get_expr_vector(i, selgenes, matrix) == get_expr_vector(j, selgenes, matrix) and duplications[i]['to_treat']==True:
#                 duplications[i]["duplications"].append(j)
#                 duplications[j]["to_treat"] = False
#     return duplications


def readouts_diff_calculation(aff_i, aff_j, readouts, matrix) -> float:
    aff_i_vector = get_expr_vector(aff_i,readouts, matrix)
    aff_j_vector = get_expr_vector(aff_j,readouts, matrix)
    diff = 0
    for idx in range(len(aff_i_vector)):
        diff += abs(aff_i_vector[idx]-aff_j_vector[idx])
    return diff

def readouts_maximization(redundancy_vector:dict, readouts:list, matrix:pd.DataFrame) -> list:
    perturbations = list()
    for aff in range(len(redundancy_vector)):
        max_pair = 0
        max_idx = 0
        for idx in range(len(redundancy_vector[aff])):
            aff_i = redundancy_vector[aff][idx][0]
            aff_j = redundancy_vector[aff][idx][1]
            readout_diff = readouts_diff_calculation(aff_i, aff_j, readouts, matrix)
            if readout_diff > max_pair:
                max_pair = readout_diff
                max_idx = idx 
        perturbations.append(redundancy_vector[aff][max_idx])
    return perturbations
    

def readouts_mean(redundancy_vector:dict, readouts:list, matrix:pd.DataFrame) -> list:
    readouts_i = dict()
    readouts_j = dict()
    mean_i = 0
    mean_j = 0
    for r_gene in readouts:
        readouts_i[r_gene] = list()
        readouts_j[r_gene] = list()
        for aff in range(len(redundancy_vector)):
            redundancies_nb = len(redundancy_vector[aff])
            for red in redundancy_vector[aff]:
                mean_i += matrix.at[red[0], r_gene]
                mean_j += matrix.at[red[1], r_gene]
            readouts_i[r_gene].append(mean_i/redundancies_nb)
            readouts_j[r_gene].append(mean_j/redundancies_nb)
            
    return readouts_i, readouts_j



def create_midas_file(expr_data: pd.DataFrame, readout_genes: str, ii_genes: list, cells: list, out_dir: str, class_name: str, intermediate_genes:list) :

    # Clean lists # TODO Check if it's really necessary
    gene_mtx = set(expr_data.columns)
    readout_genes = list(set(readout_genes).intersection(gene_mtx))
    ii_genes = list(set(ii_genes).intersection(gene_mtx))



    inputs_dict = dict()
    DA_readouts_dict = dict()
    DV_readouts_dict = dict()

    class_len = len(cells)

    inputs_dict[f'TR:{class_name}:CellLine'] = [1 for x in range(class_len*2)]

    for ii_gene in ii_genes:
        current_column = list()
        for cell in cells:
            # idx = expr_data.index[expr_data.index == cell].tolist()[0]
            expr_value = expr_data.at[cell, ii_gene]
            # if intermediates genes (stimuli for caspo) reverse the expression value
            if ii_gene in intermediate_genes:
                expr_value = 1 if expr_value==0 else 0
            current_column.append(expr_value)
        inputs_dict[f'TR:{ii_gene}'] = current_column+current_column

    for r_gene in readout_genes:
        DA_readouts_dict[f'DA:{r_gene}'] = [0 for x in range(class_len)]+[10 for x in range(class_len)]

        current_column = list()
        for cell in cells:
            # idx = expr_data.index[expr_data['Name'] == cell].tolist()[0]
            current_column.append(expr_data.at[cell, r_gene])

        DV_readouts_dict[f'DV:{r_gene}'] = [0 for x in range(class_len)]+current_column

    inputs_df = pd.DataFrame(inputs_dict)
    DA_readouts_df = pd.DataFrame(DA_readouts_dict)
    DV_readouts_df = pd.DataFrame(DV_readouts_dict)

    # concatenate the 3 dataframe
    midas_df = pd.concat([inputs_df, DA_readouts_df, DV_readouts_df], axis=1)
    midas_df.to_csv(f'{out_dir}/{class_name}_midas.csv', index=False)
 
def MinMaxNormalization(df):
    min_val = df.min().min()
    max_val = df.max().max()
    print(min_val, max_val)
    normalized_df = (df - min_val) / (max_val - min_val)
    normalized_df.to_csv('./DFFFF_TODEL.csv')
    df.to_csv('./DF_TODEL.csv')
    return normalized_df

def create_midas_file2(expr_data: pd.DataFrame, readout_genes: str, ii_genes: list, cells: list, out_dir: str, class_name: str, intermediate_genes:list, curent_readout_values:list) :

    # Clean lists # TODO Check if it's really necessary
    gene_mtx = set(expr_data.columns)
    readout_genes = list(set(readout_genes).intersection(gene_mtx))
    ii_genes = list(set(ii_genes).intersection(gene_mtx))



    inputs_dict = dict()
    DA_readouts_dict = dict()
    DV_readouts_dict = dict()

    class_len = len(cells)

    inputs_dict[f'TR:{class_name}:CellLine'] = [1 for x in range(class_len*2)]

    for ii_gene in ii_genes:
        current_column = list()
        for cell in cells:
            # idx = expr_data.index[expr_data.index == cell].tolist()[0]
            expr_value = expr_data.at[cell, ii_gene]
            # if intermediates genes (stimuli for caspo) reverse the expression value
            if ii_gene in intermediate_genes:
                expr_value = 1 if expr_value==0 else 0
            current_column.append(expr_value)
        inputs_dict[f'TR:{ii_gene}'] = current_column+current_column

    for r_gene in readout_genes:
        DA_readouts_dict[f'DA:{r_gene}'] = [0 for x in range(class_len)]+[10 for x in range(class_len)]

        # current_column = list()
        # for readout_val in curent_readout_values:
        #     # idx = expr_data.index[expr_data['Name'] == cell].tolist()[0]
        #     current_column.append(expr_data.at[cell, r_gene])

        DV_readouts_dict[f'DV:{r_gene}'] = [0 for x in range(class_len)]+curent_readout_values[r_gene]

    inputs_df = pd.DataFrame(inputs_dict)
    DA_readouts_df = pd.DataFrame(DA_readouts_dict)
    DV_readouts_df = pd.DataFrame(DV_readouts_dict)
    print(DV_readouts_df)

    DV_readouts_df = MinMaxNormalization(DV_readouts_df)

    print("++++++++++++++++++++++++++")
    print(DV_readouts_df)
    print("===============================")
    # concatenate the 3 dataframe
    midas_df = pd.concat([inputs_df, DA_readouts_df, DV_readouts_df], axis=1)
    midas_df.to_csv(f'{out_dir}/{class_name}_midas.csv', index=False)
 

def create_setup_file(ii_genes:list, readouts_genes:list, input_genes:list, intermediates_genes:list, out_dir:str):
    setup = {
        'stimuli' : list(),
        'inhibitors' : list(),
        'readouts' : list()
    }
    setup['readouts'] = readouts_genes
    for gene in ii_genes:
        if gene in input_genes:
            setup['stimuli'].append(gene)
        elif gene in intermediates_genes:
            setup['inhibitors'].append(gene)
    json.dump(setup, open(f"{out_dir}/setup.json", 'w'))

def transform_sif_file(out_dir:str):
    with open(f'{out_dir}/reduced_pkn.sif') as f:
        content = f.read()
    to_replace = {
        'PART_OF' : '1',
        'ACTIVATION' : '1',
        'INHIBITION' : '-1',
        '\"' : '',
        ' ': '',
        '.*\sUNKNOWN\s.*\n':''
    }
    for key,value in to_replace.items():
        # content = content.replace(key, value)
        content = re.sub(key, value, content)
    with open(f'{out_dir}/reduced_transformed_pkn.sif', 'w') as f:
        f.write(content)


def learn_format_args(args:dict, class_:str)-> Namespace:
    formatted_args = dict()
    out_dir = f"{args['output_dir']}"
    formatted_args['out'] = f'{out_dir}/{class_}'
    formatted_args['pkn'] = f'{out_dir}/reduced_transformed_pkn.sif'
    formatted_args['midas'] = f'{out_dir}/{class_}_midas.csv'
    formatted_args['time'] = 10
    formatted_args['fit'] = args['fit'] if args['fit'] else 0
    formatted_args['size'] = args['size'] if args['size'] else 0
    formatted_args['length'] = args['length'] if args['length'] else 0
    formatted_args['optimum'] = args['optimum'] if args['optimum'] else 0
    # Non-used params - Caspo default value
    formatted_args['discretization'] = 'round'  
    formatted_args['factor'] = 100             
    formatted_args['threads'] = 1              
    formatted_args['conf'] = 'many'            
    return Namespace(**formatted_args)


def classify_format_args(args:dict, class_:str)-> Namespace:
    formatted_args = dict()
    out_dir = f"{args['output_dir']}"
    formatted_args['out'] = f'{out_dir}/{class_}'
    formatted_args['pkn'] = f'{out_dir}/reduced_transformed_pkn.sif'
    formatted_args['setup'] = f'{out_dir}/setup.json'
    formatted_args['networks'] = f'{out_dir}/{class_}/networks.csv'
    formatted_args['midas'] = [f'{out_dir}/{class_}_midas.csv', 10]
    # Non-used params - Caspo default value   
    formatted_args['threads'] = 1              
    formatted_args['conf'] = 'many'         
    return Namespace(**formatted_args)


def visualize_format_args(args:dict, class_:str)-> Namespace:
    formatted_args = dict()
    out_dir = f"{args['output_dir']}"
    formatted_args['out'] = f'{out_dir}/{class_}'
    formatted_args['pkn'] = f'{out_dir}/reduced_transformed_pkn.sif'
    formatted_args['setup'] = f'{out_dir}/setup.json'
    formatted_args['midas'] = [f'{out_dir}/{class_}_midas.csv', 10]
    formatted_args['networks'] = f'{out_dir}/{class_}/networks.csv'
    formatted_args['sample'] = args['sample'] if args['sample'] else -1
    # Non-used params  
    formatted_args['stats_networks'] = None              
    formatted_args['stats_strategies'] = None           
    formatted_args['strategies'] = None 
    formatted_args['designs'] = None 
    formatted_args['behaviors'] = None 
    formatted_args['predictions'] = None 
    return Namespace(**formatted_args)



def dot_to_pdf(config:dict, class_:str):
    out_dir = f'{config["output_dir"]}/{class_}'
    # TODO Handle sample = [0,N]
    # if config['sample'] != '-1' :
    #     to_convert = ['networks-union', 'pkn', 'pkn-zip']
    # else :
    #     to_convert = ['pkn', 'pkn-zip']
    to_convert = ['networks-union', 'pkn', 'pkn-zip']
    for file in to_convert:
        in_path = f'{out_dir}/{file}.dot'
        out_path = f'{out_dir}/{file}.pdf'
        graph = pydotplus.graphviz.graph_from_dot_file(in_path)
        graph.write_pdf(out_path)


def run_bns_inference(config):
    out_dir = config['output_dir']
    classes = config['class_types']
    matrix_filename = f"{out_dir}/bin_reduced_matrix.csv"
    readout_filename = f"{out_dir}/no_successors_in_the_matrix.txt"
    readout_genes = read_file(readout_filename)

    # genes_hash_map = json.load(open(f'{out_dir}/genes_hash_map.json'))
    # cells_hash_map = json.load(open(f'{out_dir}/cells_hash_map.json'))
    # classes_hash_map = json.load(open(f'{out_dir}/classes_hash_map.json'))

    # genes_hash_map = ''
    # cells_hash_map = ''
    # classes_hash_map = ''

    answer_set_filename = f"{out_dir}/pseudo_perturbation_answer_sets.txt"
    sel_genes, affinities = load_answer_set(answer_set_filename)

    input_genes = read_file(f"{out_dir}/no_predecessors_in_the_matrix.txt")
    intermediate_genes = read_file(f"{out_dir}/intermediates_in_the_matrix.txt")
    expr_data = load_data(matrix_filename, index_name='Name')
    # expr_data
    
    # Redundancies calculation 
    classes = list()
    for c in affinities[0]:
        classes.append(get_cell_class(c, expr_data))
    classes_cells = get_classes_cells(matrix=expr_data, classes=classes)
    # redundancies = redundancies_calculation(affinities=affinities, selgenes=sel_genes, class_names=classes, class_cells=classes_cells,matrix=expr_data)
    # json.dump(redundancies, open(f'{out_dir}/redundancies.json', 'w'), indent=4, sort_keys=True)

    # Maximization of readouts difference
    json.dump(affinities, open(f"{out_dir}/affinities.json",'w'))
    
    redondancy_cells, redondancy_vectors = redundancies_calculation(affinities, sel_genes, classes, expr_data)
        # duplications = duplications_calculation(affinities=affinities, selgenes=sel_genes, matrix=expr_data)
    json.dump(redondancy_cells, open(f"{out_dir}/redondancy_cells.json",'w'))
    json.dump(redondancy_vectors, open(f"{out_dir}/redondancy_vectors.json",'w'))
    readouts = readouts_mean(redondancy_vectors, readout_genes, expr_data)
    # perturbations = affinities
    # json.dump(perturbations, open(f"{out_dir}/perturbations.json",'w'))



    # Preprocessing step
    for i  in range(len(classes)):
        class_ = classes[i]
        if not os.path.exists(f'{out_dir}/{class_}'): os.makedirs(f'{out_dir}/{class_}')
        # TO TREAT !!!
        # curent_cells = [item[i] for item in perturbations] 
        curent_cells = [item[i] for item in affinities]
        curent_readout_values = readouts[i]
        create_midas_file2(ii_genes=sel_genes, expr_data=expr_data, cells=curent_cells,
                        readout_genes=readout_genes, out_dir=out_dir, class_name=class_, intermediate_genes=intermediate_genes, 
                        curent_readout_values=curent_readout_values)
    
    transform_sif_file(out_dir)
    create_setup_file(sel_genes, readout_genes, input_genes, intermediate_genes, out_dir)

    for class_ in classes:
        print(f"--> Class: {class_}")
        # Learning step
        learn_args = learn_format_args(config, class_)
        learn_handler(learn_args)

        # Classify step
        # classify_args = classify_format_args(config, class_)
        # classify_handler(classify_args)

        # Visualize step
        visualize_args = visualize_format_args(config, class_)
        visualize_handler(visualize_args)
        dot_to_pdf(config, class_)

    