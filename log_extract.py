#!/usr/bin/env python3
import argparse, jgot

def parser4LogExtract(): 
    parser = argparse.ArgumentParser(description=
    'Extract something from a given ONIOM output (.log/.out)',
    formatter_class=argparse.RawTextHelpFormatter)
    basic = parser.add_argument_group(description="Required and basic options:")
    basic.add_argument('-g',dest='gau_log',type=str,required=True,help=
    '''Required! A gaussian .log/.out files with at least one geometry and energy
       e.g. -g ts_clac_IRC_RC.log/.out
    ''')
    basic.add_argument('-o',dest='out_file_name',type=str,default="Empty",help=
    '''optional. The output filename with supported suffix (currently: .com/.gjf/.pdb), 
       e.g. -o system_ircRC.com/.gjf
       e.g. -o scanned_point7.pdb
       when -ener is not given, -o is required!
    ''')
    basic.add_argument('-ener',dest='energy', action='store_true',help=
    '''optional. The energy details will be printed,             
        e.g. -ener
    when -o is not given, -ener is required!
    ''')
    basic.add_argument('-t',dest='input_template',type=str,default="FromLog",help=
    '''optional. A .com/.gjf template input, with full MM atomtype and charge, connections
       default is taken from the .com/.gjf file that generates the .log/.out file
       When .log/.out does not contain PDBinfo, -t must be given!
       e.g. -t mod_rc1.com/.gjf
    ''')
    
    extract = parser.add_argument_group(description="Options for .log/.out extracting:")
    extract.add_argument('-s',dest='step',type=str,default='last',help=
    '''optional. Step number for geometry extraction, default is the last stationary point.
Examples: (OP: optimization; SP: stationary point; PN: IRC point number)
    -s 1      >> extract the first OP geometry
    -s 7      >> extract the 7th geometry, SP prioritize OP, if no SP, will be the 7th OP
    -s lastpt >> extract the last geometry in OP list (Maybe not a SP)
    -s 7p     >> extract the 7 OP geometry despites that there are more than 7 SPs
# The following examples extract multple geometries to multiple outputs:
    -s all    >> extract all the geometries, SP/PN prioritize OP if there is any
    -s allp   >> extract all the OP geometries, SP/PN included''')
    extract.add_argument('-measure',dest='measure',default=False,type=str,help=
    '''optional. Default is not print or print scanned coords if it is a scan job.
    e.g. ":UPG@C1p-UNK@O5"''')
    extract.add_argument('-mf',dest="mask_file",default=False,type=str,help=
    '''optional. A mask file list.                               e.g. -mf masks.lst''')
    extract.add_argument('-pkl',dest="pickle_act",action='store_true',help=
    '''optional. dump/load coordinates to/from a pickle (.pkl) file automatically
    (with the same prefix as .log/.out file)
    by given the -pkl option, will decide if it is dump or load''')
    extract.add_argument('-qm',dest='qm_extract',action='store_true',help=
    '''optional. flag to extract the high-lyaer atoms (and linked hydrogen),
    used with -o -s ; if -qm is given, -o outputs only the HL atoms''')
    extract.add_argument("-v",dest='verbose',action='store_true',help=
    '''optional. verbose output, print more information.         e.g. -v''')

    com_gjf= parser.add_argument_group(description="Options for getting the .com/.gjf output(s):")
    com_gjf.add_argument('-linda',dest='lindaworker',action='store_true',help=
    '''optional. If LindaWorkers are needed for parallelization, e.g. -l''')
    com_gjf.add_argument('-p',dest='processors',type=str,default='SameAsLog',help=
    '''optional. Number of CPUs, default is SameAsLog;           e.g. -p 32''')
    com_gjf.add_argument('-m',dest='mem',type=str,default='SameAsLog',help=
    '''optional. Memory size, default is SameAsLog;              e.g. -m 108GB''')
    com_gjf.add_argument('-r',dest='routine',type=str,default='SameAsLog',help=
    '''optional. Statement of routine in a pair of quota, default is SameAsLog
    e.g. "#T oniom(b3lyp/6-31g(d):amber=SoftOnly)\\ngeom=connectivity"
    e.g. -r template # will keep the routine from template .com/.gjf''')
    com_gjf.add_argument('-title',dest='title',type=str,default='SameAsLog',help=
    '''optional. Title section of the output com file, Default is SameAsLog
    e.g.  --title "no title
    e.g.  --title "template", get it from a template .com''')
    com_gjf.add_argument('-gencon',dest="GenConnect",action='store_true',help=
    '''optional. Generate new connectivity for the first extracted geometry, 
    if given, the "$GAUSS_EXEDIR/newzmat" will be used, please ensure it exists
    e.g. -gencon
    ''')

    pdb_out = parser.add_argument_group(description="Options for getting the .pdb output(s):")
    pdb_out.add_argument('-shift',dest='resnum_shift',type=int,default=0,help=
    '''optional. number of shift for restore the PDB residue number, default is 0
    e.g. -shift 2''')
    return parser

IRC, newzmat = False, False
parsers = parser4LogExtract()
args = parsers.parse_args()
OutName = args.out_file_name
if not args.energy and OutName == "Empty":
    print("at least use -ener or -o, or use both of them!")
    exit(1)
logInit = jgot.read_log_inits(args.gau_log,args.verbose)

TemplateInput = jgot.get_template_name(args.input_template,args.gau_log)    
comlst1 = jgot.com2list(TemplateInput)
blank_ind1 = jgot.get_blank_line_ind(comlst1)
chg_mul = jgot.get_chg_multiplicity(comlst1, blank_ind1)
coords1 = jgot.get_coord(comlst1, blank_ind1)
others = jgot.get_others(comlst1, blank_ind1)
if "PDBName" not in coords1[0] or "(" not in coords1[0]:
    print("{} does not contain PDB info.!".format(TemplateInput))
    print("Error, we need a .com/.gjf file with PDB info!")
    exit(1)
if args.verbose:
    print("Number of atoms in tmeplate com: {}".format(len(coords1)))

### Get more log info, either from log or previous dumpped pickle
# [[indicies of Stationary points],[lists of coordinates[+energy],
# irc_pt_list, links_scale]]    <- logInfo
# new coord structure: coord[-1] = energy_geom
# lne(coord) = len(coords1) or len(coords1)+1
logInfo = jgot.pickle_io(args.gau_log,args.pickle_act,len(coords1))

print("Total number of geometries (OP): {}".format(len(logInfo[1])))
if len(logInfo[0]) == 0:
    print("Warning: No Stationary/Path geometry found!")
else:
    print("Number of Stationary/Path points: {}".format(len(logInfo[0])))

if args.verbose:
    if len(logInfo[0]) != 0:
        print("Stationary indices in all geometries: "+str(logInfo[0]))
    if len(logInfo[1]) != logInfo[2]:
        print("Not all the geometries have energy record")

if args.measure:
    measures = jgot.get_index_by_mask(args.measure,coords1)
    measures = [measures]
elif args.mask_file:
    masks = jgot.parse_mask_file(args.mask_file)
    measures = [jgot.get_index_by_mask(x,coords1) for x in masks]
elif logInit[-1][0] != "NOREDUNDANT":
    measures = jgot.get_index_from_modredundant(logInit[-1])
else:
    measures = []

# if it is an IRC log
if len(logInfo[2]) > 0:
    path = {"0":"TS","1":"FORWARD","2":"REVERSE"}
    IRC = [[x[0],path[x[1]]] for x in logInfo[2]]
    IRC[0] = ['0', 'TS']
else:
    IRC = False

if args.energy:
    jgot.print_log_ener(logInfo,measures,coords1,IRC)
if OutName == "Empty":
    exit(0)

### Finally extracted the desired geometries
OutPrefix = OutName.split('.')[0:-1]
out_steps = jgot.step_parser(args.step,logInfo[0],len(logInfo[1]))
if len(out_steps) > 1:
    print("going to extract all (stationary) geometries")
for ii in range(len(out_steps)):
    out_step = out_steps[ii]
    if len(out_steps) > 1:
        if IRC:
            NN = "pt{:03d}_{}_{:03d}".format(ii,IRC[ii][1],int(IRC[ii][0]))
        else:
            NN = "{:03d}".format(ii+1)
        tmp =  OutPrefix+[NN,OutName.split('.')[-1]]
        OutName = '.'.join(tmp)
    else:
    # print which IRC point will be extracted if it is an IRC log
        if IRC:
            if out_step == 0:
                print("The starting geometry (TS) will be extracted")
            else:
                point_path = IRC[int(args.step)-1]
                print("The {} th point in the {} path will be extracted".format(*point_path))

    if len(logInfo[1][out_step])-len(coords1) == 1:
        coord2 = logInfo[1][out_step][:-1]
    else:
        coord2 = logInfo[1][out_step] # this is a geometry without energy
        print("This geometry has 0 energy records, if it's not the last, please check!")

    new_coords = jgot.add_pdb_info(coords1, coord2)

    if args.qm_extract:
        new_coords = jgot.get_qm_coords(new_coords,logInfo[-1])

    if args.GenConnect:
        connectivity = jgot.newzmat_conn(new_coords,chg_mul)    
    else:
        connectivity = jgot.get_connection(comlst1, blank_ind1)

    if OutName.split('.')[-1] in ['com', 'gjf']:
        headings = jgot.get_newcom_headings(args,logInit,OutName,comlst1,blank_ind1)
        if args.qm_extract:
            chg_mul2 = chg_mul[0].split()[-2:]
            chg_mul2 = ' '.join(chg_mul2) + '\n'
            outputs = headings+['\n',chg_mul2]+new_coords+['\n','\n']
        else:
            outputs = headings+['\n']+chg_mul+new_coords+['\n']+connectivity+others
    elif OutName.split('.')[-1] == 'pdb':
        PDB = jgot.coords2pdb(new_coords,args.resnum_shift)
        outputs = ''.join(PDB)

    outFile = open(OutName, 'w')
    outFile.write(''.join(outputs))
    outFile.close()