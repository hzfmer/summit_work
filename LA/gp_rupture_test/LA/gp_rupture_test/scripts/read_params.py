#!/usr/bin/env python
''' Read parameters in param.sh
'''

import argparse

def _convert_arg_line_to_args(arg_line):
    for arg in arg_line.split():
        if not arg.strip():
            continue
        if arg[0] == '#':
            break
        yield arg

def read_params():
    parser = argparse.ArgumentParser(description="Read the parameters in param.sh, output the list contains the options and values", fromfile_prefix_chars="@")
    parser.convert_arg_line_to_args = _convert_arg_line_to_args

    parser.add_argument('-X', type=int)
    parser.add_argument('-Y', type=int)
    parser.add_argument('-Z', type=int)
    parser.add_argument('--NBGX', type=int)
    parser.add_argument('--NBGY', type=int)
    parser.add_argument('--NBGZ', type=int)
    parser.add_argument('--NEDX', type=int)
    parser.add_argument('--NEDY', type=int)
    parser.add_argument('--NEDZ', type=int)
    parser.add_argument('--NSKPX', type=int, default=1)
    parser.add_argument('--NSKPY', type=int, default=1)
    parser.add_argument('--NSKPZ', type=int, default=1)
    parser.add_argument('--DT', type=float)
    parser.add_argument('--TMAX', type=float)
    parser.add_argument('--NSRC', type=int)
    parser.add_argument('--NST', type=int)
    parser.add_argument('--NTISKP', type=int)
    parser.add_argument('--SXRGO', default=repr('output_sfc/SX_0_'))
    parser.add_argument('--SYRGO', default=repr('output_sfc/SY_0_'))
    parser.add_argument('--SZRGO', default=repr('output_sfc/SZ_0_'))
    parser.add_argument('--READ_STEP', type=int)
    parser.add_argument('--WRITE_STEP', type=int)
    parser.add_argument('--IVELOCITY', type=int, default=0)

    args = vars(parser.parse_known_args(['@param.sh'])[0])

    return args

if __name__ == "__main__":
    args = read_params()
    for key, value in args.items():
        print(key, value, '\n')
