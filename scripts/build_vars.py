#!/usr/bin/env python3

import sys
import traceback as tb
import json
from shutil import copy
import subprocess as sub

LIKWID_PERF = 0
VTUNE_PERF = 0

def usage():
    print("build_vars <main_file> <header_name> <var_1> ... [<var_n>]")

def command(cmd):
    return sub.check_output(cmd, stderr=sub.STDOUT, shell=True).decode('utf-8')

def compile_flags(compiler):
    flags = ' -g '
    version = command(compiler + ' --version')    
    if 'Intel' in version:
        flags += '-dynamic -qopt-report=5 -qopenmp -restrict'
    else:
        flags += '-fopenmp -ftree-vectorize -funroll-all-loops'
        for lib in ['m', 'omp']:
            flags += ' -l' + lib
    flags += ' -O3'
    return flags

def main():
    if len(sys.argv) < 2:
        usage()
        return

    main_file = sys.argv[1]
    header_file = sys.argv[2]
    variants = sys.argv[3:]

    nruns = 7
    nretries = 3
    nthreads = [1, 2, 4, 8, 16, 24, 32, 64]
    
    compiler = 'cc'
    compile_stub = compiler + ' '
    
    consts = {'DIM': 3, 'NUMCELLS': 64, 'DATAFLOW_CODE': 1, 'MPI_ENABLE': 1}
    if LIKWID_PERF:
        consts['LIKWID_PERF'] = 1
    elif VTUNE_PERF:
        consts['VTUNE_PERF'] = 1

    for const in consts:
        compile_stub += ' -D' + const + '=' + str(consts[const])

    compile_stub += compile_flags(compiler)

    includes = ['include', 'examples/Euler/src', 'pdfl/src', 'out']
    if VTUNE_PERF:
        includes.append(os.environ['VTUNE_AMPLIFIER_XE_2018_DIR'] + '/include')
    for include in includes:
        compile_stub += ' -I' + include

    links = []
    if VTUNE_PERF:
        links.append(os.environ['VTUNE_AMPLIFIER_XE_2018_DIR'] + '/lib64')
    for link in links:
        compile_stub  += ' -L' + link    

    libs = []
    if VTUNE_PERF:
        libs.append('ittnotify')
    for lib in libs:
        compile_stub  += ' -l' + lib

    if VTUNE_PERF:
        command('module load vtune/2018.up2')

    compile_stub += ' ' + main_file
    run_data = {}

    #srun -n $n -c $t $SDE -d -iform 1 -omix sde_${suffix}.out -i -top_blocks 500 -global_region -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -- $exe
    fp = open("results.json", "a")
    for variant in variants:
        #Compile the variant
        print('%s => %s' % (variant + '.h', header_file + '.h'))
        copy(variant + '.h', header_file + '.h')
        exec_file = variant.split('/')[-1] + '.x'
        
        compile_cmd = compile_stub + ' -o ' + exec_file
        if 'qopt' in compile_cmd:
            compile_cmd += ' -qopt-report-file=' + exec_file.replace('.x', '.rpt')
        print(compile_cmd)
        out = command(compile_cmd)
        print(out)
  
        # Run the variant
        run_stub = './' + exec_file +' ' + str(nruns)
        run_data = {}
        for n in nthreads:
            run_cmd = run_stub
            if n > 1:
                srun_cmd = 'srun -n ' + str(n) + ' '
                if VTUNE_PERF and n == nthreads[-1]:
                    srun_cmd += 'sde64 -hsw -d -iform 1 -omix ' + exec_file.replace('.x', '.sde')
                    srun_cmd += ' -i -top_blocks 500 -global_region -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -- '
                run_cmd = srun_cmd + run_cmd
            print(run_cmd)

            retry = 0
            while retry < nretries:
                try:
                    out = command(run_cmd)
                    run_data[exec_file] = {}
                    result = out.split(':')[1].strip()
                    if len(run_data[exec_file]) < 1:
                        for elem in result.split(','):
                            (key, val) = elem.split('=')
                            run_data[exec_file][key] = [val]
                    else:
                        for elem in result.split(','):
                            (key, val) = elem.split('=')
                            run_data[exec_file][key].append([val])
                    break
                except Exception as e:
                    retry += 1

            print(out)
        print(json.dumps(run_data, sort_keys=True, indent=4), file=fp) 
    fp.close()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt as e: # Ctrl-C
        print("Closing gracefully on keyboard interrupt...")
    except Exception as e:
        print('ERROR: ' + str(e))
        tb.print_exc()
