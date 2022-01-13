#!/usr/bin/env python3

"""Extract thermo data from a Lammps output file.

Python function and CLI based on argeparse.
"""

import sys, argparse, io, logging
import pandas as pd

def thermo_from_stream(in_stream,
                       start_flg='Per MPI rank memory', end_flg='Loop time of',
                       header_char='',
                       debug=False):
    """Parse stream of Lammps output to extract thermo data.

    Beginning and end of thermo loops are given by start and end flags.
    Thermo is converted into Pandas Dataframe using CVS utility. Header of thermo labels Dataframe columns.
    A list of Dataframes is returned.
    """

    # Set up LOGGER
    c_log = logging.getLogger('thermo_from_stream')
    # Adopted format: level - current function name - mess. Width is fixed as visual aid
    std_format = '[%(levelname)5s - %(funcName)10s] %(message)s'
    logging.basicConfig(format=std_format)
    c_log.setLevel(logging.INFO)
    # Set debug option
    if debug: c_log.setLevel(logging.DEBUG)

    parse_flg = False
    c_thermo_output = []
    thermo_data = []
    i=0
    for line in in_stream.readlines():
        if end_flg in line:
            parse_flg = False
            c_log.debug("----End parse at line %5i", i)
            # Join all lines and parse as a CSV file, but with whitesapce as delimiter instead of comma.
            # Assume floating point data. First row will automatically be header.
            thermo_data.append(pd.read_csv(io.StringIO('\n'.join(c_thermo_output)),
                                           delim_whitespace=True, dtype=float)
            )
            c_thermo_output = []

        if parse_flg:
            if header_flg:
                line = header_char + line
                header_flg = False
            c_thermo_output.append(line)


        if start_flg in line:
            parse_flg = True
            c_log.debug("++Start parse at line %5i", i)
            header_flg = True
        i+=1

    # Deal with incomplete runs
    if len(c_thermo_output):
        c_log.warning("Incomplete run of %i timesteps", len(c_thermo_output))
        try:
            thermo_data.append(pd.read_csv(io.StringIO('\n'.join(c_thermo_output)),delim_whitespace=True, dtype=float))
        except ValueError or TypeError:
            c_log.warning("Skip invalid last line %s" %  c_thermo_output[-1])
            thermo_data.append(pd.read_csv(io.StringIO('\n'.join(c_thermo_output[:-1])),delim_whitespace=True, dtype=float))

    # Return list of dataframes with runs output
    return thermo_data

def extract_thermo(argv):
    #-------------------------------------------------------------------------------
    # Argument parser
    #-------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description=__doc__)
    # Positional arguments
    parser.add_argument('filename',
                        type=str, nargs='?',
                        help='input file. If not given stdin is used.')
    # Optional args
    parser.add_argument('--start_flg', '-s',
                        dest='start_flg', type=str, default='Per MPI rank memory',
                        help='flag string to start parsing. This line is not included in output. Default "Per MPI rank memory".')
    parser.add_argument('--end_flg', '-e',
                        dest='end_flg', type=str, default='Loop time of',
                        help='flag string to end parsing. This line is not included in output. Defaul "Loop time of"')
    parser.add_argument('--to_files', '-f',
                        dest='basename', default=None, type=str,
                        help='save each run thermo to different file: <basename>-1, <basename>-2, etc.')
    parser.add_argument('--head',
                        dest='header_char', default='', type=str,
                        help='trail every header line with this string, e.g. "#" for gnuplot.')
    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        help='show debug informations.')

    #-------------------------------------------------------------------------------
    # Initialize and check variables
    #-------------------------------------------------------------------------------
    args = parser.parse_args(argv)

    # Set up LOGGER
    c_log = logging.getLogger(__name__)
    # Adopted format: level - current function name - mess. Width is fixed as visual aid
    std_format = '[%(levelname)5s - %(funcName)10s] %(message)s'
    logging.basicConfig(format=std_format)
    c_log.setLevel(logging.INFO)
    # Set debug option
    if args.debug: c_log.setLevel(logging.DEBUG)
    c_log.debug(args)

    if args.filename:
        in_stream = open(args.filename, 'r')
    else:
        in_stream = sys.stdin

    thermo_data = thermo_from_stream(in_stream, args.start_flg, args.end_flg,
                                     header_char=args.header_char,
                                     debug=args.debug)
    if args.filename: in_stream.close()


    if args.basename != None:
        for i, data in enumerate(thermo_data):
            c_log.debug("Print run %5i", i)
            with open("%s-%i.dat" % (args.basename, i), 'w') as out_stream:
                print(data.to_string(), file=out_stream)
    else:
        for i, data in enumerate(thermo_data):
            c_log.debug("Print run %5i", i)
            print(data.to_string())


if __name__ == "__main__":
    # From https://github.com/python/mypy/issues/2893
    # "Python interpreter [...] overrides the default handling of SIGPIPE; specifically, it ignores this signal."
    # Restore it to the default handler, SIG_DFL
    import signal
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    extract_thermo(sys.argv[1:])
