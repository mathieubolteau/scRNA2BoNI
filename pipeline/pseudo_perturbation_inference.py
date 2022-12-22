#!/usr/bin/env python3
# -*- coding: utf-8 -*

try:
    import csv
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy import stats
    import json
    import os
    from sys import argv
    from argparse import Namespace
    import configparser
    from pipeline.pyBRAvo.src import pyBravo
    import shutil
    from math import ceil


except ImportError as E:
    print(E)
    print('Verify your dependancies.')
    exit()




def run_pseudo_perturbation_inference(config):
    out_dir = config['output_dir']
    k = config['k']
    i = ceil(k*config['gene_inputs_pourcentage']/100)
    expression_instance_filename = f'{out_dir}/expression_instance.lp'
    input_instance_filename = f'{out_dir}/inputs_instance.lp'
    problem_instance = "./pipeline/data/problem.lp"
    timeout = config['timeout'] if config['timeout'] != '' else 0

    cmd = f'clingo --const k={k} --const i={i} --time-limit={timeout} {expression_instance_filename} {input_instance_filename} {problem_instance}'
    print(f'CMD:  {cmd}')
    clingo_output = os.popen(cmd).read()
    print(clingo_output)
