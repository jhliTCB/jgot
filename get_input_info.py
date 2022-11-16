#!/usr/bin/env python3
import argparse
import re
import jgot

def parser_func():
    parser = argparse.ArgumentParser(description='Update the coordinates of a gaussian input template',
    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i','--input_com',type=str,required=True,help=
    '''    A com input file that has correct atom types and charges
    str, Required, e.g. -i ts_clac_IRC_RC.com
    ''')
    parser.add_argument('-measure',dest='measure',default=False,help=
    '''
        Atom numbers or atom PDB information for measurements, seperated by "-"
    default is False, Examples (using Amber-like atom mask)
        -measure 7192-198-1784 ;
        -measure :UPG@O3B-:UPG@C1p-:UNK@O5
        -measure :HIS:16@NE2-:UNK:465@O5
    ''')
    parser.add_argument('-mf','--mask_file',default='None',help=
    '''
        A file contains one or multiple lines atom list for measurement, '#' for comment
        e.g.:
        -mf active_coords.txt
        # Active coordinates for TS
        :UPG@C1p-:UNK@O5          # the CV1 distance
        :UNK@O5-:UPG@C1p-:UPG@O3B # the CV2 angle
        :UNK@O5-:HID:16@NE2       # comment3
    ''')
    parser.add_argument('-v','--verbose',action='store_true',help=
    '''    verbose output information for link0 and title
    ''')
    return parser

parsers = parser_func()
args = parsers.parse_args()
comInfo = jgot.get_com_info(args.input_com, args.verbose)
print('\n'.join(comInfo))

comlst = jgot.com2list(args.input_com)
blanks = jgot.get_blank_line_ind(comlst)
coords = jgot.get_coord(comlst,blanks)
connections = jgot.get_connection(comlst, blanks)
others = jgot.get_others(comlst, blanks)

if args.verbose:
    print("")
    modred = jgot.parse_com_end(others, "modredundant")
    if len(modred) > 0:
        for i in range(len(modred)):
            non_num = re.findall(r'[A-Za-z]+', modred[i])
            if len(non_num) == 1:
                print("Warning: no action for {}".format(modred[i]))
                num_lst = modred[i][modred[i].index(non_num[0][-1])+1:]
            elif len(non_num) == 2:
                num_lst = modred[i][modred[i].index(non_num[0][-1])+1:modred[i].index(non_num[1][-1])]
            else:
                print("Warning: only numbers in modredundant!{}".format(''.join(modred[i])))
                break
            print("Modredundant info: {}".format(''.join(modred[i])))
            num_lst = num_lst.split()
            xyzs = []
            for n in num_lst:
                coord = coords[int(n)-1]
                print("{} is: {}".format(n,coord.split()[0]))
                xyzs.append([float(coord.split()[x]) for x in [2,3,4]])

            if len(num_lst) == 2:
                print("initial distance of {} {} is: {:.3f}".format(*num_lst, jgot.measure(xyzs)))
            elif len(num_lst) == 3:
                print("initial angle of {} {} {} is : {:.3f}".format(*num_lst, jgot.measure(xyzs)))
            elif len(num_lst) == 4:
                print("initial dihedral of {} {} {} {} is : {:.3f}".format(*num_lst, jgot.measure(xyzs)))
            print("")

    print("")
    print("Number of Connection Records: {}".format(len(connections)))

if args.measure:
    pdb_info, xyzs = jgot.get_pdb_info_and_xyzs(args.measure, coords)
    print("Measurement for {}: {:.3f}".format(pdb_info, jgot.measure(xyzs)))
elif args.mask_file != "None":
    masks = jgot.parse_mask_file(args.mask_file)
    for mask in masks:
        pdb_info, xyzs = jgot.get_pdb_info_and_xyzs(mask, coords)
        print("Measurement for {}: {:.3f}".format(pdb_info, jgot.measure(xyzs)))
elif args.measure and args.mask_file != "None":
    print("--measure and --mask_file can not be used together")
    print("not going to measure anything")