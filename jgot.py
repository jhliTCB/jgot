# jgot: JY Gaussian ONIOM Toolkit, rc0
# By Junhao Li, PhD, jun828hao@gmail.com
# TODO merge charge2zero.py get_input_info.py gvcom2molUPcom.py into oniom_input.py
# TODO add -modredundat option, add modredundant section according to given mask
# TODO read frequencies information
# TODO print/modifie freezed residues, atom ids, QM residues, atom ids in com file
# TODO to use class??? logger?
import math
import pickle
import re, os

def Num2Element(num):
    tupl_table = ('H', 'He', 'Li', 'Be', 'B' , 'C' , 'N' , 'O', 'F', 'Ne', 'Na',
            'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar','K', 'Ca', 'Sc', 'Ti', 'V',
            'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se',
            'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Te', 'Ru', 'Rh', 'Pd',
            'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce',
            'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm','Yb',
            'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb',
            'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu',
            'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm','Md', 'No', 'Lr', 'Rf', 'Db','Sg',
            'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og', 'Uue')
    elem_number = tuple([x for x in range(1,len(tupl_table)+1)])
    Element_dic = dict(zip(elem_number,tupl_table))
    return Element_dic[num]

def get_template_name(templ,log):
    TemplateInput = '.'.join(log.split('.')[0:-1])
    if templ == "FromLog":
        Input_isfile = False
        for in_suffix in ['.com','.gjf']:
            if os.path.isfile(TemplateInput+in_suffix):
                TemplateInput += in_suffix
                Input_isfile = True
                break
        if not Input_isfile:
            err="Error, the input file with the same prefix of {} does not exist"
            print(err.format(log))
            print("Please use -t name.com to provide a template oniom input file")
            print("(which should have the same MM settings with the log file)")
            exit(1)
    else:
        if os.path.isfile(templ):
            TemplateInput = templ
        else:
            print("{} does not exist!".format(templ))
            exit(1)
    return TemplateInput

def step_parser(args_step,stat_list,op_len):
    out_step = ['last']
    if "all" in args_step:
        if len(stat_list) > 0:
            out_step = stat_list
        else:
            out_step = [i for i in range(op_len)]
    ## ":" and "," are not practical in most situations !
    ##elif ":" in args_step:
    ##    out_step = args_step.split(':')
    ##elif "," in args_step:
    ##    out_step = args_step.split(',')
    elif "allp" in args_step:
        out_step = [i for i in range(op_len)]
    else:
        if args_step == 'last' and len(stat_list) != 0:
            out_step[0] = stat_list[-1]
            print("The last Stationary geometry will be extracted")
        elif args_step == 'last' and len(stat_list) == 0:
            out_step[0] = -1
            print("The geometry will be extracted, but it's not a stationary point")
        elif len(stat_list) != 0 and 'p' not in args_step:
            out_step[0] = stat_list[int(args_step)-1]
            info = "Step_ID {} given, extracting the {} th stationary/IRC point"
            print(info.format(args_step, int(args_step)))
        elif 'p' in args_step and args_step != 'lastpt':
            out_step[0] = int(args_step.split('p')[0])-1
            info = "Step_ID {} given, extracting the {} th optimization point"
            print(info.format(args_step, out_step[0]+1))
        elif args_step == 'lastpt':
            out_step[0] = -1
            print("The last optimization, maybe stationary, will be extracted")
        else:
            out_step[0] = int(args_step)-1
    return out_step

def com2list(com):
    return [line for line in open(com, 'r')]

def get_blank_line_ind(comlst):
    i=0
    blank_line_ind = []
    for line in comlst:
        # Some blank lines are not completely blank!
        if line.strip(" ") == '\n':
            blank_line_ind.append(i)
        i += 1
    return blank_line_ind

def get_link0_section(comlst, blank_line_ind):
    return [x for x in comlst[0:blank_line_ind[0]] if x.startswith('%')]

def get_routine_section(comlst, blank_line_ind):
    return [x for x in comlst[0:blank_line_ind[0]] if not x.startswith('%')]

def get_title_section(comlst, blank_line_ind):
    return comlst[blank_line_ind[0]+1:blank_line_ind[1]]

def get_chg_multiplicity(comlst, blank_line_ind):
    return [comlst[blank_line_ind[1]+1]]

def get_coord(comlst, blank_line_ind):
    return comlst[blank_line_ind[1]+2:blank_line_ind[2]]

def get_connection(comlst, blank_line_ind):
    return comlst[blank_line_ind[2]+1:blank_line_ind[3]]

def get_others(comlst, blank_line_ind):
    return comlst[blank_line_ind[3]:]

def add_pdb_info(coords1, coords2):
    # coords1: xyzs of all atoms with PDBinfo
    # coords2: xyzs without PDBinfo
    new_xyz = []
    pdbForm = " {:60s} {:4d} {:15.7f} {:15.7f} {:15.7f} "
    if len(coords1) != len(coords2):
        print("n_atoms differed:{} {}\n".format(len(coords1), len(coords2)))
        exit(1)
    for i in range(0, len(coords1)):
        lle1 = coords1[i].split()
        lle2 = coords2[i].split()
        if "(" in lle2[0]:
            element2 = lle2[0].split('-')[0]
        elif len(re.findall(r'[A-Z][a-z]', lle2[0])) == 0:
            element2 = Num2Element(int(lle2[0].strip()))
        else:
            element2 = lle2[0]
        if lle1[0].strip().split('-')[0] != element2:
            print("Error, element not match! something wrong??")
            exit(1)
        for j in range(2, 5):
            lle1[j] = lle2[j]
        flx = [float(x) for x in lle1[2:5]]
        line2 = pdbForm.format(lle1[0],int(lle1[1]),flx[0],flx[1],flx[2])
        line2 += ' '.join(lle1[5:])+'\n'
        new_xyz.append(line2)
    return(new_xyz)

def remove_pdb_info(coords,chg_mul):
    # return a gaussian input file without PDB info for newzmat
    outName = ".tmp_new_coords.gjf"
    outFile = open(outName,'w')
    nzLink0 = "%mem=16GB\n%nprocshared=8\n# hf/6-31g(d)\n\n No Title\n\n{}\n"
    outFile.write(nzLink0.format(chg_mul))
    for line in coords:
        elemt = line.strip()[0:2].strip('-') # no 3-leter element
        coord = line.split()[2:5]
        outFile.write(" {} {} {} {}\n".format(elemt,*coord))
    outFile.write('\n\n\n')
    outFile.close()
    return outName

def scale_xyz(H_xyz, L_xyz,s):
    # s is the last scaling factor read from log file (scaling for high layer)
    # For a two-layer calculation, the scale factors specify the link atom bond distance 
    # in the model system when calculated at the low and high levels (respectively)
    #return scaled L_xyz' [x3, y3, z3] for making it available in new high calc
    return [s*(L_xyz[i]-H_xyz[i])+H_xyz[i] for i in range(len(L_xyz))]

def get_qm_coords(new_coords,link_lines):
    #new coords come from extraction
    #link_indx is a list output from read log
    # out is a coords with each line seperated
    #exmaple link statement:
    #" ONIOM: Cut between C /H  1368 and C  1370 factor= 0.723886 0.723886"
    if len(link_lines) > 0:
        linkAtoms = [x.split()[4].strip('/') for x in link_lines]
        link_indx = [re.findall(r'[0-9].[0-9]+',x) for x in link_lines]
    out = []
    for line in new_coords:
        llist = line.split()
        if llist[5] == 'H':
            #elem = five_repl(llist[0])[0]
            #out.append(" ".join([elem]+llist[2:5]+['\n']))
            out.append(line)
        elif len(llist) > 6 and llist[6][0] == 'H':
            linIndx = llist[7]
            linLine = new_coords[int(linIndx)-1] #HL line
            for i in range(len(link_lines)):
                if linIndx in link_indx[i]:
                    scaling = link_indx[i][-1]
                    linkAtm = linkAtoms[i]
                    pdbInfo = new_coords[int(link_indx[i][1])-1].split()[0]
                    indx = get_chg_indx(pdbInfo)
                    dash = indx[2]
                    pdbInfo = five_repl(pdbInfo)
                    break
            # using "dash" is much more better than using repl_five!!!!!
            newInfo = [linkAtm,linkAtm,dash+pdbInfo[2+indx[1]],pdbInfo[3+indx[1]],
                        linkAtm+pdbInfo[4+indx[1]]] + pdbInfo[5+indx[1]:]
            newInfo = ' {}-{}-{}({}={},{}={},{}={}) 0 '.format(*newInfo)
            newList = scale_xyz([float(linLine.split()[x]) for x in [2,3,4]],
                        [float(llist[x]) for x in [2,3,4]], float(scaling))
            #outputs = [linkAtm+llist[1]]+[str(x) for x in newList]+[llist[5],'\n']
            #outputs = [pdbInfo,llist[1]]+[str(x) for x in newList]+['H', '\n']
            #outputs = [linkAtm]+[str(x) for x in newList]+['\n']
            outputs = [newInfo]+[str(x) for x in newList]+['H', '\n']
            out.append(" ".join(outputs))
    return out

def five_repl(a):
    # retrun a list of the first column of the input
    # e.g. a = 'N-N--0.415700(PDBName=N,ResName=GLY,ResNum=450)'
    #      b = 'C-CA-0.807600(PDBName=CZ,ResName=ARG,ResNum=449)'
    #            0    1    2   3          4          5      6         7      8         9     10
    # a return ['N', 'N', '', '0.415700', 'PDBName', 'N', 'ResName', 'GLY', 'ResNum', '450', '']
    #            0    1        2          3          4      5         6      7         8      9 
    # b return ['C', 'CA',    '0.807600', 'PDBName', 'CZ', 'ResName', 'ARG', 'ResNum', '449', '']
    _r = a.replace('-','_').replace('(','_').replace(')','_').replace(',','_').replace('=','_')
    return _r.split('_')

def repl_five(a, b):
        ind = b[1]
        tmp = tuple(a[0:2]+a[2+ind:9+ind])
        return ' {}-{}-{}({}={},{}={},{}={})'.format(*tmp)

def zero_charges(coords, res):
    # residue: list, e.g. [HID, 16]
    new_coords = []
    orig_chg = []
    Charged, HL = False, False
    if res[0] in ['ASP', 'GLU', 'ARG', 'LYS', 'HIP']:
        Charged = True
        Total = []
    for i in range(len(coords)):
        linelist = coords[i].split()
        atominfo = linelist[0]
        indx = get_chg_indx(atominfo)
        if Charged:
            tmp_lst2 = five_repl(atominfo)
            charges2 = float(indx[0])*float(tmp_lst2[2+indx[1]])
            Total.append(charges2)

        pdb_info = [x.split('=') for x in atominfo.split('(')[1].strip(')').split(',')]
        if pdb_info[1][1] == res[0] and pdb_info[2][1] == res[1]:
            if linelist[5] == 'H':
                HL = True
                print("residue {}{} contains high-layer atoms, will be skipped".format(*res))
                break #not going to change High-layer residues
            otherInf = coords[i][coords[i].index(')')+1:]
            tmp_lst = five_repl(atominfo)
            charge = float(indx[0])*float(tmp_lst[2+indx[1]])
            orig_chg.append(charge)
            tmp_lst[2+indx[1]] = '0.000000'
            new_coords.append(' {:<60s} {}'.format(repl_five(tmp_lst, indx), otherInf))
        else:
            new_coords.append(coords[i])
    if HL:
        return []
    if len(orig_chg) == 0:
        print("Warning: residue {}{} not Found, do nothing!".format(*res))
    else:
        print("Total charge of residue {}{} is {}".format(res[0],res[1],sum(orig_chg)))
        if round(sum(orig_chg)) != 0:
            if Charged:
                nchg = sum(Total)-sum(orig_chg)
                print("Charged, SUM charge of the new system is {} ({})".format(round(nchg), nchg))
            else:
                print("Warning: sum charge of a neutralized residue is not Zero!, please Check!!")
        print("Charges of all atoms in residue {}{} has changed to zero".format(*res))
    print('')
    return new_coords

def read_oniom_log(log, natoms): 
    # Extract geometries and energies; Convergence follows the ONIOM energy flag
    OP,joa, SP, jsa, conv = 0, 0, 0, 0, 0
    OP_coords, SP_coord_ind = [], []
    flag, flag2 = 'NULL', 'NULL'
    links_scale = []
    irc_pt_list = [] # same length with SP_coord_ind
    energy_flag = False
    conver_flag = False # " Converged?" in line
    for line in open(log, 'r'):
        if "ONIOM: Cut between" in line and line not in links_scale:
        #will be risky if Gaussian can addapt the scaling factor during opt
            links_scale.append(line)
        elif "Path Number:" in line:
            irc_pt_list.append(re.findall(r"[0-9]+",line))

        if len(line) > 50:
            if line[43:66] == 'Coordinates (Angstroms)':
                OP += 1
                op_coord, joa = [], 0
                #abc.append('%5s%9s\n' % tuple(['MODEL', OP]))
                if not flag == 'STATIONARY_POINT':
                    flag = 'COORD_START'
                else:
                    jsa = 0
                    SP += 1
                    flag = 'STAT_COORD_START'
            if flag2 == 'op_collecting' and "--------" in line:
                if joa > 0 and joa != natoms:
                    print("Error reading: Atom number ({}) not match the given ({})!".format(joa, natoms))
                    exit(1)
                else:
                    flag = 'COORD_END'
                    flag2 = 'NULL'
            elif flag2 == 'sp_collecting' and "--------" in line:
                if jsa >0 and jsa != natoms:
                    print("Error reading: Atom number ({}) not match the given ({})!".format(jsa, natoms))
                    exit(1)
                else:
                    flag = 'STAT_COORD_END'
                    flag2 = 'NULL'
        # IRC log won't have "Stationary"; opt log won't have "Path Number", using the same flag
        elif line[4:17] == '-- Stationary' or 'Path Number:' in line:
            flag = 'STATIONARY_POINT'
        # make flags for the energies
        elif "ONIOM: calculating energy" in line:
            energy_flag = True
            energy_geom = []

        # Collect the coordinates 
        # OP or SP statment won't show up at the same time
        if flag == 'COORD_START':
            if line.strip()[0] in [str(x) for x in range(10)]:
                flag2 = 'op_collecting'
                op_coord.append(line[10:].strip()) # atom ID is not needed!
                joa += 1
        elif flag == 'STAT_COORD_START':
            if line.strip()[0] in [str(x) for x in range(10)]:
                op_coord.append(line[10:].strip()) # atom ID is not needed!
                flag2 = 'sp_collecting'
                jsa += 1
        elif flag == 'COORD_END':
            OP_coords.append(op_coord)
            flag = "NULL"
        elif flag == 'STAT_COORD_END':
        # OP-2: OP has firstly += 1 and we need coordinates just before "Stationary"
            SP_coord_ind.append(OP-2) 
            OP_coords.append(op_coord)
            flag = "NULL"
        if energy_flag: 
        # The followed 4 lines are energy records for:
        # low_model, high_model, low_read, and total
        # The last geometry usually has no energies, or has the same with previous
            if "calculating energy" not in line:
                energy_geom.append(line.split()[-1])
            if "ONIOM: extrapolated" in line:
                OP_coords[-1].append(energy_geom) 
                energy_flag = False

        if " Converged?" in line:
            conver_list = []
            conver_flag = True
        elif conver_flag:
            conver_list.append([x.strip() for x in line[22:].split()]) 
            conv += 1
            if conv == 5:
                OP_coords[-1].append(conver_list) # array of 6 x 3
                conver_flag = False
                conv = 0

    return [SP_coord_ind, OP_coords, irc_pt_list, links_scale]

def read_log_inits(log,verbose):
    # get the initial data and then break!
    # return [[Gau_version,mem,nproc,chk,routine,title,charge,multiplicity],
    #        [coords],[modredundant],
    outInfo, coord, mod_red = [], [], []
    routine, title, chg_mul = "", "", ""
    star, flag, mord = False, False, False
    dash, chrg, gauv = 0, 0, 0
    link0 = []
    for line in open(log, 'r'):
        if "**********" in line:
            star = True
        if star and line[0:9] == " Gaussian":
            gauv = line.strip()
            star = False
        if line[0:2] == " %" and 'linda' not in line.lower():
            # linda parallelization is not needed in most cases
            link0.append(line.strip())
        if '---' in line:
            dash += 1
        if dash == 3 and '---' not in line:
            routine += line.strip()
        elif dash == 5 and '---' not in line:
            title += line.strip()
        elif dash == 6 and '---' not in line:
            if "Charge =" in line:
                chrg += 1
                chg_mul += " "+" ".join(re.findall(r'(-?\d+\.?\d*)', line))
            elif chrg == 3 and len(coord) == 0:
                flag = True
            if flag and len(line.strip()) == 0:
                flag = False
            elif flag:
                coord.append(line)
            if "following ModRedundant input" in line:
                mord = True
            if mord and "ModRedundant" not in line and len(line.strip())>0:
                mod_red.append(line.strip())
            if mord and len(line.strip()) == 0:
                mord = False
            if "Read MM parameter" in line:
                break
    if len(mod_red) == 0:
        mod_red = ["NOREDUNDANT"]
    if "PDBName" not in coord[0] or "(" not in coord[0]:
        print("log is ran from com without PDB info!")
        print("-t is needed to specify a .com/.gjf file with PDB info!")
    outInfo = [gauv,link0,routine,title,chg_mul.strip(),coord,mod_red]
    if verbose:
        print("Basic .log/.out information:")
        print("GauVersion: {}".format(gauv))
        print("Link0: {}".format(sorted(link0)))
        print("Routine: {}".format(routine))
        print("Charge and multiplicity: {}".format(chg_mul.strip()))
        print("modredundant info: {}".format("\n".join(mod_red)))
        print("Number of atoms from log file:   {}".format(len(coord)))
    return outInfo        

def parse_crglst(crglst):
    out = []
    for line in open(crglst, 'r'):
        out.append([x.strip() for x in line.split()])
    return out

def get_chg_indx(AtmInf):
    out = []
    if '--' in AtmInf:
        out = [-1.0, 1, "-"]
    else:
        out = [1.0,  0, ""]
    return out

def get_com_info(com, verbose):
    #parse a com input file!
    #return multiplicity, charges, atom numbers, routine, scan info, title, etc
    # with -v: link0, routine, title
    # without -v: chg_mul settings, charges in each layer, 
    #             numbers of total atoms/freezed atoms/non-freezed atoms/H/L atoms
    # layer/freezes: #atoms #charge_sum
    outinf = []
    comlst = com2list(com)
    blanks = get_blank_line_ind(comlst)
    mulchg = get_chg_multiplicity(comlst, blanks) # system_low model_low model_high
    coords = get_coord(comlst,blanks)
    HL, LL, FF, NF, LA, SC = [], [], [], [], [], []
    for i in range(len(coords)):
        linelist = coords[i].split()
        atominfo = linelist[0]
        indx     = get_chg_indx(atominfo)
        atomlist = five_repl(atominfo)
        atcharge = float(indx[0])*float(atomlist[2+indx[1]])
        SC.append(atcharge)
        if linelist[1] == "0":
            NF.append(atcharge)
        elif linelist[1] == "-1":
            FF.append(atcharge)
        if linelist[5] == "H":
            HL.append(atcharge)
        if linelist[5] == "L":
            LL.append(atcharge)
        if len(linelist) > 6:
            #PDBName, ResName, ResID, 4,6,8
            atom_tmp = [atomlist[x+indx[1]] for x in [4,6,8]]
            link_at1 = "{}@{}{}".format(*atom_tmp)
            linkingH = linelist[6]
            atom_id2 = int(linelist[7])-1
            atmList2 = five_repl(coords[atom_id2].split()[0])
            indx2    = get_chg_indx(coords[atom_id2].split()[0])
            atom_tmp = [atmList2[x+indx2[1]] for x in [4,6,8]] # this is a way to get pdbInfo
            link_at2 = "{}@{}{}".format(*atom_tmp)
            LLinks = [link_at1, linelist[5], str(i), link_at2, coords[atom_id2].split()[5], str(atom_id2)]
            formIn = "{:12s} ({:2s}, index: {:6s}) Links to {:12s} ({:2s}, index: {:6s})"
            LinkInfo = formIn.format(*LLinks)
            LinkInfo += " by Atom-AtomType: {}".format(linkingH)
            LA.append(LinkInfo)
    if verbose:
        outinf.append(' '.join([x.strip() for x in get_link0_section(comlst, blanks)]))
        outinf.append(''.join(get_routine_section(comlst, blanks)))
        outinf.append(''.join(get_title_section(comlst, blanks)))
        outinf.append("")
    outinf.append("Charges and Multiplicities")
    outinf.append("real_low, model_high, model_low): {}".format(" ".join([x.strip() for x in mulchg])))
    outinf.append("Number of atoms:       {:6d}".format(len(coords)))
    outinf.append("Num_free_atoms =       {:6d}, Num_Freezed_atoms = {:6d}".format(len(NF),len(FF)))
    outinf.append("Num_High_layer =       {:6d}, Num_Low_layer     = {:6d}".format(len(HL),len(LL)))
    if verbose and len(LA) > 0:
        outinf.append('\n'.join(LA))
    outinf.append("Sum of atomic Charges in real         = {:.3f}".format(sum(SC)))
    outinf.append("Sum of atomic Charges in high (model) = {:.3f}".format(sum(HL)))
    outinf.append("Sum of atomic Charges in Low-Layer    = {:.3f}".format(sum(LL)))
    return outinf

def parse_com_end(com_endings, what_to_return):
    # return a list of modredundant atom number
    modredun = ['B','R','A','D','L','X','O']
    abbremod = ['Bond', 'Stretch', 'Angle', 'Bend', 'Dihedral', 'Torsion', 'Linear', 'LinearBend',]
    if what_to_return == "modredundant":
        out = []
        flag = 0
        for line in com_endings:
            #print(line.split())
            if line != '\n' and line.split()[0] in modredun+abbremod:
                out.append(" ".join([x.strip() for x in line.split()]))
                flag = 1
            elif line == '\n' and flag == 1:
                break
        return out

def get_index_from_modredundant(modredundant):
    out = [] #2D list
    for line in modredundant:
        #for T in ['S','F','A','B','K','R','D','H']:
        #    if T in line and 'S' in line:
        #        out.append([int(x)-1 for x in line.split('S')[0].split()[1:]])
        #        break
        if 'S' in line: T = 'S'
        if 'F' in line: T = 'F'
        if T:
            out.append([int(x)-1 for x in line.split(T)[0].split()[1:]])
    return(out)

def vector(a, b):
    return tuple([b[0]-a[0], b[1]-a[1], b[2]-a[2]])

def dot(a, b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def cross(a, b):
    return tuple([a[1]*b[2]-a[2]*b[1],  a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]])

def measure(xyzs):
    if len(xyzs) == 2: # xyzs[0]-xyzs[1] distance
        return math.sqrt(sum([(xyzs[0][i]-xyzs[1][i])**2 for i in range(3)]))
    elif len(xyzs) == 3: # xyzs[0]-xyzs[1]-xyzs[2] angle
        BA, BC = vector(xyzs[1],xyzs[0]), vector(xyzs[1],xyzs[2])
        cosB = dot(BA,BC) / math.sqrt(dot(BA, BA)) / math.sqrt(dot(BC,BC))
        rads = math.acos(cosB)
        return math.degrees(rads)
    elif len(xyzs) == 4: # xyzs[0]-xyzs[1]-xyzs[2]-xyzs[3] dihedral
        BA, BC = vector(xyzs[1],xyzs[0]), vector(xyzs[1],xyzs[2])
        CB, CD = vector(xyzs[2],xyzs[1]), vector(xyzs[2],xyzs[3])
        n1, n2 = cross(BA,BC), cross(CB,CD)
        cosD = dot(n1, n2) / math.sqrt(dot(n1,n1)) / math.sqrt(dot(n2,n2))
        rads = math.acos(cosD)
        if dot(n1, CD) < 0:
            return math.degrees(rads)
        else:
            return -1.0*math.degrees(rads)

def get_mask_by_index(atom, coords):
    # atom index pairs as input, e.g. 0-7128-7124
    # or atom index list directly: [0, 7128, 7124]
    # return amber-like mask, e.g. :PRO:2@N-:UPG:452@C1p-:UNK:453@O5
    j = 0
    if type(atom) is list:
        atoms = atom
    elif type(atom) is str:
        atoms = atom.split('-')
    for i in atoms:
        line = coords[int(i)]
        ind1 = get_chg_indx(line)
        info = five_repl(line.split()[0])
        ind2 = [info[x+ind1[1]] for x in [6,8,4]]
        text = ":{}:{}@{}".format(*ind2)
        if j == 0:
            out = text
        else:
            out = out + "-" + text
        j += 1
    return out

def parse_mask(mask):
    # example1 :UPG@O3B-:UPG@C1p-:UNK@O5
    # example2 :HIS:16@NE2-:UNK:465@O5
    # return: 
    atomlist = [] # PDBName-ResName/ResNum
    for x in mask.split('-'):
        if len(x.split(':')) == 2: # Res_name/Res_num not recomanded
            tmp = x.split(':')[1].split('@')
            if len(re.findall(r'[A-Za-z]', tmp[0])) == 0:
                ResName = ""
                ResNum  = tmp[0]
            else:
                ResName = tmp[0]
                ResNum  = ""
        elif len(x.split(':')) == 3:
            ResNum  = re.findall(r'[0-9]+', x)[0]
            ResName = re.findall(r'[A-Za-z]+', x)[0]
        PDBName = x.split(':')[-1].split('@')[-1]
        atomlist.append([PDBName, ResName, ResNum])
    return atomlist

def get_index_by_mask(mask, coords):
    # atom masks as input, e.g. :PRO:2@N-:UPG:452@C1p
    # return indicies of these atoms, e.g. 0-7128-7124
    inds = []
    risk = {1:2, 2:1} # risky list for risky mask
    atomlist = parse_mask(mask) #[PDBName, ResName, ResNum]
    for atom in atomlist:
        tmps = []
        if "" in atom:
            shift = atom.index("")
        else:
            shift = 0
        for i in range(len(coords)):
            line = coords[i]
            ind1 = get_chg_indx(line)
            info = five_repl(line.split()[0])
            ind2 = [info[x+ind1[1]] for x in [4,6,8]]
            if ind2[0] == atom[0]:
                if shift > 0:
                    if ind2[risk[shift]] == atom[risk[shift]]:
                        tmps.append(str(i))
                elif ind2[1] == atom[1] and ind2[2] == atom[2]:
                    tmps.append(str(i))
        if len(tmps) == 0:
            print("Error, {} is not found in the atom section, is it a correct mask?".format(atom))
            print("Out-of-range error will be triggered")
        elif len(tmps) > 1:
            print("Warning! more than one atoms were found for {}, will use the first".format(atom))
        inds.append(tmps[0])
    return "-".join(inds)

def get_xyz_by_index(atom, coords):
    # atom index pairs as input, e.g. 0-7128-7124
    # or a direct list of atom indices, e.g. [0, 7128, 7124]
    # return xyzs list of these atoms
    out = []
    if type(atom) is list:
        atoms = atom
    elif type(atom) is str:
        atoms = atom.split('-')
    for i in atoms:
        out.append([float(coords[int(i)].split()[x]) for x in [2,3,4]])
    return out

def get_pdb_info_and_xyzs(mask,coords):
    # the 'mask' is eigther amber-like mask, or a given atom indices
    if len(re.findall(r'[A-Za-z]', mask)) == 0:
        pdb_info = get_mask_by_index(mask,coords)
        xyz_indices = mask
    else:
        pdb_info = mask
        xyz_indices = get_index_by_mask(mask,coords)
    xyzs = get_xyz_by_index(xyz_indices,coords)
    return pdb_info, xyzs

def parse_mask_file(mf):
    # parse the masks in a mask_file list
    out = [] # 'mask1', 'mask2', ...
    for line in open(mf, 'r'):
        mask = line.split('#')[0].strip()
        if len(mask) != 0:
            out.append(mask)
    return out

def coords2pdb(coords,shift):
    # coords must have PDB infor record! TODO add remark info to the head
    #"{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}\
    #{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(...)
    PDBForm = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}" + \
            "{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n"
    outList = []
    f_xyzs  = [[float(xx) for xx in RR.split()[2:5]] for RR in coords ]
    for i in range(len(coords)):
        record = coords[i]
        indx = get_chg_indx(record.split()[0])
        reclst = five_repl(record)
        elemt  = reclst[0].strip()
        pdbNam = reclst[4+indx[1]]
        resNam = reclst[6+indx[1]]
        resNum = int(reclst[8+indx[1]])+shift
        if resNum > 9999:
            resNum -= 9999
        xyz = f_xyzs[i]
        #charge = float(indx[0])*float(reclst[2+indx[1]])
        atmlst = ['ATOM',i+1,pdbNam," ",resNam,"X",resNum," ",
                xyz[0],xyz[1],xyz[2],1.0,1.0,elemt,"  "]
        outList.append(PDBForm.format(*atmlst))
    return outList

def newzmat_conn(new_coords,chg_mul):
    if "GAUSS_EXEDIR" in os.environ:
        if os.name == 'posix':
            newzmat = [x for x in os.environ["GAUSS_EXEDIR"].split(":") if "bsd" not in x]
            newzmat = newzmat[0]+"/newzmat"
        elif os.name == 'nt':
            newzmat = os.environ["GAUSS_EXEDIR"]+os.sep+"newzmat.exe"
        tmp_com = ".tmp_new_connection.gjf"
        tmp_inp = remove_pdb_info(new_coords,' '.join(chg_mul[0].split()[0:2]))
        formats = [newzmat,tmp_inp,tmp_com]
        if "16" in os.environ["GAUSS_EXEDIR"]:
            cmd_run = "{} -icart -gencon {} {}".format(*formats)
        elif "09" in os.environ["GAUSS_EXEDIR"]:
            cmd_run = "{} {} {}".format(*formats)
        print(cmd_run)        
        os.system(cmd_run)
        with open(tmp_com, 'a') as tmporary_com:
            tmporary_com.write("\n\n\n")
        return get_connection(com2list(tmp_com),get_blank_line_ind(com2list(tmp_com)))
    else:
        print("The GAUSS_EXEDIR environment is not set, please load the gaussian module!")
        exit(1)

def pickle_io(logName,action,N_atoms):
    pkl_name = "{}.pkl".format(".".join(logName.split('.')[0:-1]))
    if not os.path.isfile(pkl_name):
        info = read_oniom_log(logName,N_atoms)
        print("The list logInfo has been read from {}".format(logName))
        if action:
            with open(pkl_name, 'wb') as f_pkl:
                pickle.dump(info, f_pkl)
            print("The list logInfo has dumpped to {}".format(pkl_name))
    else:
        with open(pkl_name, 'rb') as f_pkl:
            info = pickle.load(f_pkl)
        print("The list logInfo has been loaded from {}".format(pkl_name))
    return info

def conv_minus(conv):
    # conv: 6 x 3
    out = []
    for c in conv:
        yesno = c[-1]
        if yesno == 'NO':
            if '*' in c[0]:
                yesno += '_********'
            else:
                distc = abs(float(c[0])-float(c[1]))
                yesno += str("_{:.6f}".format(distc))
        out.append("{:8s}".format(yesno))
    return out

def print_log_ener(logInfo,measure_list,coords1,IRC,conv_flag):
    if len(logInfo[0]) != 0:
        ener_tile = ["StationaryID","Low_model","High_model","Low_real","Extrapolated"]
        coods_out = [logInfo[1][i] for i in logInfo[0]]
    else:
        ener_tile = ["OptimizationID","Low_model","High_model","Low_real","Extrapolated"]
        coods_out = logInfo[1]

    if len(measure_list) > 0:
        for indes in measure_list:
            mask = get_mask_by_index(indes,coords1)
            ener_tile.append(mask)

    if IRC: # make sure IRC and coods_out have same length
        print("This is an IRC log file")
        print("The IRC information will be printed in the end")
        ener_tile += ["#IRC_Step","point_number","path_direction"]
        IDnumP = [[str(i+1)]+IRC[i] for i in range(len(IRC))]
    
    if conv_flag:
        ener_tile += ["d_Max_Force","d_RMS_Force","d_Max_Displ",
                "d_RMS_Displ","d_Max_MM_Force","d_RMS_MM_Force"]

    print(','.join(ener_tile))

    for i in range(len(coods_out)):
        if len(coods_out[i][-2]) < 5:
            eners = [str(i+1)]+coods_out[i][-2]
        else:
            eners = [str(i+1)]+['NaN']*4
        if len(measure_list) > 0:
            for indes in measure_list:
                c = measure(get_xyz_by_index(indes,coods_out[i][0:-1]))
                eners.append("{:.4f}".format(c))
        if IRC:
            eners += IDnumP[i]
        if conv_flag:
            eners += conv_minus(coods_out[i][-1])

        print(','.join(eners))

def get_newcom_headings(args,logInit,OutName,comlst1,blank_ind1):
    new_link0, new_routine = [], []
    if args.processors == "SameAsLog":
        new_link0.append(sorted(logInit[1])[2]+'\n')
    else:
        new_link0.append('%nprocshared={}\n'.format(args.processors))
    if args.mem == "SameAsLog":
        new_link0.append(sorted(logInit[1])[1]+'\n')
    else:
        new_link0.append('%mem={}\n'.format(args.mem))
    if args.lindaworker:
        new_link0.append('%LindaWorkders=node1,node2\n')
    name = '.'.join(OutName.split('.')[:-1])
    if '/' in name:
        name = name.split('/')[-1]
    new_link0.append('%chk={}.chk\n'.format(name))

    if args.routine == "template":
        new_routine = get_routine_section(comlst1, blank_ind1)
    elif args.routine == "SameAsLog":
        new_routine = [logInit[2]+'\n']
    else:
        new_routine = [args.routine+'\n']

    if args.title == 'template':
        title = get_title_section(comlst1, blank_ind1)
    elif args.title == "SameAsLog":
        title = [logInit[3]+'\n']
    else:
        title = [args.title + '\n']
    return new_link0+new_routine+['\n']+title