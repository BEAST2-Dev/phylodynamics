#!/usr/bin/env python

from argparse import ArgumentParser, FileType
from sys import argv, exit, stdin, stdout

if __name__ == "__main__":

    parser = ArgumentParser(description="Convert file with data along rows into file with data down columns.")
    parser.add_argument("infile", type=FileType('r'), default=stdin, nargs="?",
                        help="file to convert (default stdin)")
    parser.add_argument("-d", type=str, default=" ",
                        help="delimiter")
    parser.add_argument("-o", dest="outfile", metavar="outfile", type=FileType('w'), default=stdout,
                        help="output file (default stdout)")

    args = parser.parse_args(argv[1:])

    # Read input
    data = {}
    N = None
    for line in args.infile.readlines():
        thisline = line.strip().split(args.d)
        if len(thisline)>0:
            data[thisline[0]] = thisline[1:]
            if N == None:
                N = len(thisline[1:])
            else:
                if len(thisline[1:]) != N:
                       print "Rows have different numbers of elements. Aborting."
                       exit(1)


    # Generate output
    for key in data.keys():
        args.outfile.write("{}{}".format(key, args.d))
    args.outfile.write("\n")

    for i in range(N):
        for key in data.keys():
            args.outfile.write("{}{}".format(data[key][i], args.d))
        args.outfile.write("\n")
