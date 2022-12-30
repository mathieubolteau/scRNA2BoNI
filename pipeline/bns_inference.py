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

    from .utils import load_data, read_file
    from pipeline.caspo.caspo.console.handlers import learn_handler, classify_handler, visualize_handler


except ImportError as E:
    print(E)
    print('Verify your dependancies.')
    exit()


def find_last_index(search_list, starting_item):
    starting_item = '^'+starting_item
    indices = [ i for i, word in enumerate(search_list) if re.search(starting_item, word) ]
    return indices[-1]

def get_value(atom):
    return re.findall('\"([A-Z0-9\.]*)\"', atom)


def load_answer_set(filename: str, classes:list) -> set:
    sel_genes = list()
    cells = {
        classes[0] : list(),
        classes[1] : list()
    }
    file_list = read_file(filename)
    index = find_last_index(file_list, 'Optimization:')-1 # -1 to have the answer-set of the last optimization
    answer_set_raw = file_list[index]
    answer_set_list = answer_set_raw.split(' ')
    for atom in answer_set_list:
        if re.search('^selprot', atom):

            sel_genes.append(get_value(atom)[0])
        elif re.search('^affinite', atom):
            cur_cells = get_value(atom)
            cells[classes[0]].append(cur_cells[0])
            cells[classes[1]].append(cur_cells[1])
            

    return sel_genes, cells


def create_midas_file(matrix_filename: str, readout_genes: str, ii_genes: list, cells: list, out_dir: str, class_name: str):
    expr_data = load_data(matrix_filename)

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
            idx = expr_data.index[expr_data['Name'] == cell].tolist()[0]
            current_column.append(expr_data.loc[idx, ii_gene])
        inputs_dict[f'TR:{ii_gene}'] = current_column+current_column

    for r_gene in readout_genes:
        DA_readouts_dict[f'DA:{r_gene}'] = [0 for x in range(class_len)]+[10 for x in range(class_len)]

        current_column = list()
        for cell in cells:
            idx = expr_data.index[expr_data['Name'] == cell].tolist()[0]
            current_column.append(expr_data.loc[idx, r_gene])

        DV_readouts_dict[f'DV:{r_gene}'] = [0 for x in range(class_len)]+current_column

    inputs_df = pd.DataFrame(inputs_dict)
    DA_readouts_df = pd.DataFrame(DA_readouts_dict)
    DV_readouts_df = pd.DataFrame(DV_readouts_dict)

    # MAX_VALUE = max(DV_readouts_df.max())

    # def normalize(x):
    #     return x/MAX_VALUE

    # DV_readouts_df = DV_readouts_df.apply(np.vectorize(normalize))
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
    with open(f'{out_dir}/modified_pkn.sif') as f:
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
    with open(f'{out_dir}/tranformed_pkn.sif', 'w') as f:
        f.write(content)


def learn_format_args(args:dict, class_:str)-> Namespace:
    formatted_args = dict()
    out_dir = args['output_dir']
    formatted_args['out'] = f'{out_dir}/{class_}'
    formatted_args['pkn'] = f'{out_dir}/tranformed_pkn.sif'
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
    out_dir = args['output_dir']
    formatted_args['out'] = f'{out_dir}/{class_}'
    formatted_args['pkn'] = f'{out_dir}/tranformed_pkn.sif'
    formatted_args['setup'] = f'{out_dir}/setup.json'
    formatted_args['networks'] = f'{out_dir}/networks.csv'
    formatted_args['midas'] = [f'{out_dir}/{class_}_midas.csv', 10]
    # Non-used params - Caspo default value   
    formatted_args['threads'] = 1              
    formatted_args['conf'] = 'many'         
    return Namespace(**formatted_args)


def visualize_format_args(args:dict, class_:str)-> Namespace:
    formatted_args = dict()
    out_dir = args['output_dir']
    formatted_args['out'] = f'{out_dir}/{class_}'
    formatted_args['pkn'] = f'{out_dir}/tranformed_pkn.sif'
    formatted_args['setup'] = f'{out_dir}/setup.json'
    formatted_args['midas'] = [f'{out_dir}/{class_}_midas.csv', 10]
    formatted_args['networks'] = f'{out_dir}/networks.csv'
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
    answer_set_filename = f"{out_dir}/pseudo_perturbation_answer_sets.txt"
    sel_genes, cells = load_answer_set(answer_set_filename, classes)

    # Preprocessing step
    for class_ in classes:
        curent_cells = cells[class_]
        print(class_)
        create_midas_file(ii_genes=sel_genes, matrix_filename=matrix_filename, cells=curent_cells,
                          readout_genes=readout_genes, out_dir=out_dir, class_name=class_)
    
    input_genes = read_file(f"{out_dir}/no_predecessors_in_the_matrix.txt")
    intermediate_genes = read_file(f"{out_dir}/intermediates_in_the_matrix.txt")
    transform_sif_file(out_dir)
    create_setup_file(sel_genes, readout_genes, input_genes, intermediate_genes, out_dir)

    
    for class_ in classes:
        print(config)
        # Learning step
        learn_args = learn_format_args(config, class_)
        learn_handler(learn_args)

        # Classify step
        classify_args = classify_format_args(config, class_)
        classify_handler(classify_args)

        # Visualize step
        visualize_args = visualize_format_args(config, class_)
        visualize_handler(visualize_args)
        dot_to_pdf(config, class_)

