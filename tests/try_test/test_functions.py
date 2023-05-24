import os
import sys
import fileinput
import json
import subprocess as sp
import threading as th
import pandas as pd

def make_config_file(mk : float,dmk:float,fractions_filename:str,new_filename : str):
    with open('config_template.txt', 'r') as file :
        filedata = file.read()
    filedata = filedata.replace('*mk*', f'{mk}')
    filedata = filedata.replace('*dmk*', f'{dmk}')
    filedata = filedata.replace('*fractions_filename*', f'{fractions_filename}')
    with open(new_filename, 'w') as file:
        file.write(filedata)


m_progname = os.path.relpath( "../../build-Distribution-Desktop_Qt_6_4_2_MinGW_64_bit-Debug/distrib_vector_grid")
print(m_progname)

def run_prog(prog_name :str,config_name:str):
    sp.run([prog_name ,'-config', config_name])

def get_fractions(fractions_filename:str):
    with open(fractions_filename) as f:
        data = json.load(f)
    return data

def get_fraction_filename(mk : float,dmk:float):
    return f'fractions_{mk}_{dmk}.txt'


def get_config_filename(mk : float,dmk:float):
    return f'config_{mk}_{dmk}.txt'

def make_files(mk : float,dmk:float):
    fractions_filename = get_fraction_filename(mk,dmk)
    config_filename = get_config_filename(mk,dmk)
    make_config_file(mk,dmk,fractions_filename,config_filename)
    return (config_filename,fractions_filename)

def get_fractions_cycle(mk : float,dmk:float,prog_name = m_progname):
    fractions_filename = f'fractions_{mk}_{dmk}.txt'
    config_filename = f'config_{mk}_{dmk}.txt'
    make_config_file(mk,dmk,fractions_filename,config_filename)
    run_prog(prog_name,config_filename)
    return get_fractions(fractions_filename)

def test_parametr_array(parametrs_list:list,prog_name = m_progname):
    fractions_list = list(range(len(parametrs_list)))
    
    for (idx,params) in enumerate(parametrs_list):
        mk = params['mk']
        dmk = params['dmk']
        frac = get_fractions_cycle(mk,dmk,prog_name)
        frac.update({'mk':mk,'dmk':dmk})
        fractions_list[idx] = frac
    return fractions_list

def test_parametr_array_p(parametrs_list:list,prog_name = m_progname):
    fractions_list = list(range(len(parametrs_list)))
    threads = list(range(len(parametrs_list)))
    fraction_filename_list = list(range(len(parametrs_list)))
    for (idx,params) in enumerate(parametrs_list):
        mk = params['mk']
        dmk = params['dmk']
        (config_filename,fractions_filename) = make_files(mk,dmk)
        fraction_filename_list[idx] = fractions_filename
        threads[idx] = th.Thread(target = run_prog,args = (prog_name,config_filename))
        threads[idx].start()

    for (idx,params) in enumerate(parametrs_list):
        threads[idx].join()
        frac = get_fractions(fraction_filename_list[idx])
        mk = params['mk']
        dmk = params['dmk']
        frac.update({'mk':mk,'dmk':dmk})
        fractions_list[idx] = frac
    return fractions_list

def test_parametr_array_mks(mk_list:list,dmk:float,prog_name = m_progname,parallel = False):
    parametrs_list = [{'mk':mk,'dmk':dmk} for mk in mk_list]
    if(parallel):
        return  test_parametr_array(parametrs_list,prog_name)
    else:
        return  test_parametr_array_p(parametrs_list,prog_name)

def test_parametr_array_dmks(mk:float,dmk_list:list,prog_name = m_progname,parallel = False):
    parametrs_list = [{'mk':mk,'dmk':dmk} for dmk in dmk_list]
    if(parallel):
        return  test_parametr_array_p(parametrs_list,prog_name)
    else:
        return  test_parametr_array(parametrs_list,prog_name)

def load_fractions(parametrs_list:list):
    fractions_list = list(range(len(parametrs_list)))
    for (idx,params) in enumerate(parametrs_list):
        mk = params['mk']
        dmk = params['dmk']
        frac = get_fractions(get_fraction_filename(mk,dmk))
        frac.update({'mk':mk,'dmk':dmk})
        fractions_list[idx] = frac
    return fractions_list

def load_list_mk(mk_list:list,dmk:float):
    parametrs_list = [{'mk':mk,'dmk':dmk} for mk in mk_list]
    return load_fractions(parametrs_list)

def load_list_fmk(mk:float,dmk_list:list):
    parametrs_list = [{'mk':mk,'dmk':dmk} for dmk in dmk_list]
    return load_fractions(parametrs_list)

def plotting_command(columns,filename:str,x_column = 'mk'):
    cols = pd.Index(columns)
    x_col = cols.get_loc(x_column)
    return 'plot '+ ', '.join([f'\"{filename}\" using {x_col+1}:{cols.get_loc(item)+1} title \"{item}\" with lp' for item in cols if(item != 'mk' and item != 'dmk')])


def dict_concat(a,b):
    c = dict(a)
    c.update(b)
    return c

def map_dict(values:dict,func):
    return dict({key:func(values[key]) for key in values})