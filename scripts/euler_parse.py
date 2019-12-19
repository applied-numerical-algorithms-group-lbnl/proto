#!/usr/bin/env python3

import os
import sys
import traceback as tb
import json
import csv
import subprocess as sub

def usage():
    print("euler_parse.py <output_file>")

def command(cmd):
    return sub.check_output(cmd, stderr=sub.STDOUT, shell=True).decode('utf-8')

def main():
    if len(sys.argv) < 2:
        usage()
        return

    out_file = open(sys.argv[1])
    csv_file = open(sys.argv[1].replace('.out', '.csv'), 'w', newline='')
    writer = None

    for line in out_file:
        (variant, results) = line.split(':')
        pos = variant.rfind('_')
        variant = variant[2:pos]
        items = results.strip().split(',')

        par_type = 'Boxes'
        if '_hyb' in variant:
            par_type = 'Hybrid'
        elif '_pt' in variant:
            par_type = 'Tiles'

        fuse_type = 'Serial'
        if '_pf' in variant:
            fuse_type = 'Partial'
        elif '_ff' in variant:
            fuse_type = 'Full'

        tiled = ('_ts' in variant)

        data = {'variant': variant, 'par_type': par_type, 'fuse_type': fuse_type, 'tiled': tiled}
        for item in items:
            (key, val) = item.split('=')
            data[key] = val

        if writer is None:
            writer = csv.DictWriter(csv_file, fieldnames=data.keys())
            writer.writeheader()
        writer.writerow(data)

    out_file.close()
    csv_file.close()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt as e: # Ctrl-C
        print("Closing gracefully on keyboard interrupt...")
    except Exception as e:
        print('ERROR: ' + str(e))
        tb.print_exc()
