#!/usr/bin/env python3
import argparse
import jgot

def parser_func():
    parser = argparse.ArgumentParser(description='Update the coordinates of a gaussian input template',
    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i',dest='input_com',type=str,required=True,help=
    '''    A com input exported from GaussView or other place, 
    which does not have correct atom type and charge or incomplete QM region settings
    str, Required, e.g. -gv ts_clac_IRC_RC.com
    ''')
    parser.add_argument('-t','--template_com',type=str,required=True,help=
    '''    A com template input, with full MM atomtype and charge, connections, etc...
    str, Required, e.g. -t mod_rc1.com
    ''')
    parser.add_argument('-o','--output_com',type=str,required=True,help=
    '''    output com filename, e.g. -o system_ircRC.com  Required option
    ''')
    parser.add_argument('-p','--processors',type=str,default='32',help=
    '''    number of CPUs, e.g. -p 32, default is 32''')

    parser.add_argument('-m','--mem',type=str,default='64GB',help=
    '''    memory size, e.g. -m 64GB, default is 64GB''')

    parser.add_argument('-l','--lindaworker',action='store_true',help=
    '''    if LindaWorkers are needed for parallelization, -e.g. -l''')

    parser.add_argument('-r','--routine',type=str,default="tmpl",help=
    '''    Statement of routine in a pair of quota, default is from template
    e.g. -r "#T oniom(b3lyp/6-31g(d):amber=SoftOnly) geom=connectivity"
    e.g. -r keep  : keep the routine section from input;
    e.g. -r tmpl  : use the routine section from template''')
    parser.add_argument('-title',dest='title',type=str,default='no title',help=
    '''    Title section of the output com file, default is: "no title"
    e.g. -title keep  : keep the title section from input
    e.g. -title tmpl  : use the title section from template
    ''')
    return parser

parsers = parser_func()
args = parsers.parse_args()
outcom = args.output_com
comlst1 = jgot.com2list(args.template_com)
comlst2 = jgot.com2list(args.input_com)
blank_ind1 = jgot.get_blank_line_ind(comlst1)
blank_ind2 = jgot.get_blank_line_ind(comlst2)

new_link0, new_routine = [], []

new_link0.append('%nprocshared={}\n'.format(args.processors))
new_link0.append('%mem={}\n'.format(args.mem))
if args.lindaworker:
    new_link0.append('%LindaWorkders=node1,node2\n')
name = '.'.join(outcom.split('.')[:-1])
if '/' in name:
    name = name.split('/')[-1]
new_link0.append('%chk={}.chk\n'.format(name))

if args.routine == 'keep':
    new_routine = jgot.get_routine_section(comlst2, blank_ind2)
elif args.routine == 'tmpl':
    new_routine = jgot.get_routine_section(comlst1, blank_ind1)
else:
    new_routine = [args.routine+'\n']

if args.title == 'keep':
    title = jgot.get_title_section(comlst2, blank_ind2)
elif args.title == 'tmpl':
    title = jgot.get_title_section(comlst1, blank_ind2)
else:
    title = [args.title + '\n']
chg_mul = jgot.get_chg_multiplicity(comlst1, blank_ind1)
coords1 = jgot.get_coord(comlst1, blank_ind1)
coords2 = jgot.get_coord(comlst2, blank_ind2)

#if "(" not in new_coords[0].split()[0]:
#    new_coords = jgot.add_pdb_info(coords1, coords2, False)
#else:
new_coords = jgot.add_pdb_info(coords1, coords2)

connectivity = jgot.get_connection(comlst1, blank_ind1)
others = jgot.get_others(comlst1, blank_ind1)

outputs = new_link0+new_routine+['\n']+title+['\n']+chg_mul+new_coords+['\n']+connectivity+others
outcom2 = open(outcom, 'w')
outcom2.write(''.join(outputs))
outcom2.close()