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
    from scRNA2BoNI.pyBRAvo.src import pyBravo
    import shutil
    from math import ceil
    import pkg_resources


except ImportError as E:
    print(E)
    print('Verify your dependancies.')
    exit()




def run_pseudo_perturbation_inference(config):
    out_dir = config['output_dir']
    k = config['k']
    expression_instance_filename = f'{out_dir}/expression_instance.lp'
    input_instance_filename = f'{out_dir}/inputs_instance.lp'
    problem_instance = pkg_resources.resource_filename(__name__, 'data/pseudo_perturbation_inference/problem.lp')
    timeout = config['timeout'] if config['timeout'] != '' else 0

    cmd = f'clingo --const k={k} --time-limit={timeout} {expression_instance_filename} {input_instance_filename} {problem_instance} > {out_dir}/pseudo_perturbation_answer_sets.txt'
    print(f'CMD:  {cmd}')
    clingo_output = os.popen(cmd).read()
