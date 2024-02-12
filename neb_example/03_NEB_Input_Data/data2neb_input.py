#!/usr/bin/env python3
"""Convert a data file to a file which can be read
by the LAMMPS neb command.

Synopsis
--------
The neb command of Lammps uses the following format:
<N>
<id1> <x1> <y1> <z1>
(...)
<idN> <xN> <yN> <zN>
where id is the atom ID, x,y, and z are the coordinates,
and N is the number of atoms in the file.
See also the documentation of the neb command.

Author
------
Wolfram Noehring
Fri Jul 17 19:30:56 CEST 2015
"""

import re
import argparse
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__)
    parser.add_argument('infile',
        type=str, metavar = 'DATA',
        help='LAMMPS data file (input)')
    parser.add_argument('outfile',
        type=str, metavar = 'NEB-DATA',
        help='File which can be read by the Lammps neb command')
    args=parser.parse_args()

    data = to_list(args.infile, split_lines=False)
    print("Reading input configuration")
    header_section_or_comment = re.compile('.*?[a-zA-Z](?!\+|-)')
    comment = re.compile('^\#')
    for line_number, line in enumerate(data):
        if header_section_or_comment.match(line):
            # Strip comments
            words = line.rstrip().split()
            comment_pos = [i for i, word in enumerate(words)
                if comment.match(word)]
            if comment_pos:
                words = words[0:comment_pos[0]]
            is_float = [read_as_float(w) for w in words]
            if all(is_float):
                continue
            # Divide into numbers and words
            first_string_pos = is_float.index(False)
            # Join strings to form the keyword
            keyword = ' '.join(words[first_string_pos:])

            if keyword == 'atoms':
                natoms = int(words[:first_string_pos][0])
            if keyword == 'Atoms':
                atoms_section_start = line_number + 2
                break
    atoms_section = data[atoms_section_start:atoms_section_start + natoms]
    atoms_section = [line.rstrip().split() for line in atoms_section]
    # Remove the IDs and the image flags
    for i,line in enumerate(atoms_section):
        del(line[1])
        del(line[4:])
        atoms_section[i] = ' '.join(line) + '\n'
    # Write the outfile
    with open(args.outfile, 'w') as file:
        file.write('{:d}\n'.format(natoms))
        file.writelines(atoms_section)

def read_as_float(my_obj):
    """ Check if an object can be converted to a number
    See also http://stackoverflow.com/q/354038
    """
    try:
        float(my_obj)
        return True
    except ValueError:
        return False

def to_list(infile, split_lines=False):
    """ Reads a file and returns the contents as a list or list of lists."""
    with open(infile, 'r') as file:
        list_of_lines = file.readlines()
    if split_lines:
        list_of_lines = [line.rstrip().split() for line in list_of_lines
            if not line.startswith('#')]
    return list_of_lines

if __name__ == '__main__':
    main()