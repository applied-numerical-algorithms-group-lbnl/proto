#!/usr/bin/env python3

import sys
import traceback as tb
from shutil import copy
import subprocess as sub

def usage():
    print("build_vars <main_file> <header_name> <var_1> ... [<var_n>]")

def command(cmd):
    return sub.check_output(cmd, stderr=sub.STDOUT, shell=True).decode('utf-8')

def compile_flags(compiler):
    flags = ' -g '
    version = command(compiler + ' --version')    
    if 'Intel' in version:
        flags += '-dynamic -qopt-report=5 -qopenmp'
    else:
        flags += '-fopenmp -ftree-vectorize -funroll-all-loops -I/usr/include/x86_64-linux-gnu/mpich -L/usr/lib/x86_64-linux-gnu'
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

    nruns = 9
    nthreads = [1, 2, 4, 8, 16, 24, 32, 64]
    
    compiler = 'cc'
    compile_stub = compiler + ' '
    
    consts = {'DIM': 3, 'NUMCELLS': 64, 'DATAFLOW_CODE': 1, 'MPI_ENABLE': 1}
    for const in consts:
        compile_stub += ' -D' + const + '=' + str(consts[const])

    compile_stub += compile_flags(compiler)

    includes = ['include', 'examples/Euler/src', 'pdfl/src', 'out']
    for include in includes:
        compile_stub += ' -I' + include

    links = []
    for link in links:
        compile_stub  += ' -L' + link    

    compile_stub += ' ' + main_file
    run_data = {}

    for variant in variants:
        #Compile the variant
        print('%s => %s' % (variant + '.h', header_file + '.h'))
        exec_file = variant.split('/')[-1] + '.x'
        compile_cmd = compile_stub + ' -o ' + exec_file
        if 'qopt' in compile_cmd:
            compile_cmd += ' -qopt-report-file=' + exec_file.replace('.x', '.rpt')
        print(compile_cmd)
        out = command(compile_cmd)
        print(out)
  
        # Run the variant
        run_cmd = './' + exec_file +' ' + str(nruns)
        run_data = {}
        for n in nthreads:
            if n > 1:
                run_cmd = 'srun -n ' + str(n) + ' ' + run_cmd
            print(run_cmd)
            out = command(run_cmd)
            print(out)

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

    import json
    print(json.dumps(run_data, sort_keys=True, indent=4)) 

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt as e: # Ctrl-C
        print("Closing gracefully on keyboard interrupt...")
    except Exception as e:
        print('ERROR: ' + str(e))
        tb.print_exc()
