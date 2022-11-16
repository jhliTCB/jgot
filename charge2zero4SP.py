#!/usr/bin/env python3
import argparse
import jgot
default_route = "#T oniom(b3lyp/6-311++g(2d,p):amber=SoftFirst)=embedcharge geom=connectivity"

def parser_func():
    parser = argparse.ArgumentParser(description='Change the charges to zero for a set of residues in an ONIOM input',
    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i','--input_com',type=str,required=True,help=
    '''    A ONIOM input (.com) that has correct atom types and charges
    str, Required, e.g. -gv embBS2_ts_clac_IRC_RC.com
    ''')
    parser.add_argument('-z','--zerocharges',type=str,required=True,help=
    '''    A list of residues that need to change to zero charge 1 by 1,
           format example: "HID 16"
    str, Required, e.g. -z zero_charge.lst
    ''')
    parser.add_argument('-o','--output_com',type=str,required=True,help=
    '''    output .com basename, e.g. -o embBS2_TS  Required option
    ''')
    parser.add_argument('-p','--processors',type=str,default='32',help=
    '''    number of CPUs, e.g. -p 32''')

    parser.add_argument('-m','--mem',type=str,default='64GB',help=
    '''    memory size, e.g. -m 64GB''')

    parser.add_argument('-l','--lindaworker',action='store_true',help=
    '''    if LindaWorkers are needed for parallelization, -e.g. -l''')

    parser.add_argument('-r','--routine',type=str,default=default_route,help=
    '''    Statement of routine in a pair of quota,
    e.g. "#T oniom(b3lyp/6-31g(d):amber=SoftFirst)=EmbedCharge geom=connectivity"
    e.g. -r keep # will keep the routine used in input .com
    ''')
    parser.add_argument('--title',type=str,default='no title',help=
    '''    Title section of the output com file, e.g. --title "no title"
           Or keep it from template: --title keep
    ''')
    return parser

parsers = parser_func()
args = parsers.parse_args()
outcom = args.output_com
comlst = jgot.com2list(args.input_com)
empind = jgot.get_blank_line_ind(comlst)
crglst = jgot.parse_crglst(args.zerocharges)

new_link0, new_routine = [], []
new_link0.append('%nprocshared={}\n'.format(args.processors))
new_link0.append('%mem={}\n'.format(args.mem))

if args.lindaworker:
    new_link0.append('%LindaWorkders=node1,node2\n')

new_link0.append('%chk={}.chk\n'.format("CHK_REPLACE"))
if args.routine == "keep":
    new_routine = jgot.get_routine_section(comlst, empind)
else:
    new_routine = [args.routine+'\n']

if args.title == 'keep':
    title = jgot.get_title_section(comlst, empind)
else:
    title = [args.title + '\n']

chg_mul = jgot.get_chg_multiplicity(comlst, empind)
coords1 = jgot.get_coord(comlst, empind)

for res in crglst:
    new_coords = jgot.zero_charges(coords1, res)
    if len(new_coords) == 0:
        continue
    connectivity = jgot.get_connection(comlst, empind)
    others = jgot.get_others(comlst, empind)
    outname = "{}_{}.com".format(outcom, ''.join(res))
    new_link0[-1] = '%chk={}.chk\n'.format(outname[0:-4])
    outputs = new_link0+new_routine+['\n']+title+['\n']+chg_mul+new_coords+['\n']+connectivity+others
    outcom2 = open(outname, 'w')
    outcom2.write(''.join(outputs))
    outcom2.close()