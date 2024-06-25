# jgot
A toolkit of a set of functions used for the manipulations of Gaussian 16 ONIOM input (.com/.gjf) and output (.log/.out).  
  
Requirement: python3, pickle, argparse, re, os, math  
  
Initial version: 2022-11-16


压缩文件里有如下文件（所有的脚本不依赖numpy/scipy，可用系统默认的python3运行）：

**jgot.py** ----------------------- 包含了大部分函数，供其他脚本调用

**log_extract.py** -------------- 主要脚本，用以提取ONIOM log的一些信息

**get_input_info.py** ---------- 快速查看一个oniom 输入文件的基本信息

**charge2zero4SP.py** -------- 把QM区域之外的残基的电荷手动设置为0再算SP的能量



### log_extract.py的用法：


    usage: log_extract.py [-h] -g GAU_LOG [-o OUT_FILE_NAME] [-ener] [-conv]
                      [-t INPUT_TEMPLATE] [-s STEP] [-measure MEASURE]
                      [-mf MASK_FILE] [-pkl] [-qm] [-v] [-linda]
                      [-p PROCESSORS] [-m MEM] [-r ROUTINE] [-title TITLE]
                      [-gencon] [-shift RESNUM_SHIFT]

	Extract something from a given ONIOM output (.log/.out)

	optional arguments:
	  -h, --help           show this help message and exit

	  Required and basic options:

	  -g GAU_LOG           Required! A gaussian .log/.out files with at least one geometry and energy
                               e.g. -g ts_clac_IRC_RC.log/.out
                           
	  -o OUT_FILE_NAME     optional. The output filename with supported suffix (currently: .com/.gjf/.pdb), 
                               e.g. -o system_ircRC.com/.gjf
                               e.g. -o scanned_point7.pdb
                               when -ener is not given, -o is required!
                           
	  -ener                optional. The energy details will be printed,             
                               e.g. -ener
                           when -o is not given, -ener is required!
                           
	  -conv                optional, but use together with -ener
                               e.g. -ener -conv
	  -t INPUT_TEMPLATE    optional. A .com/.gjf template input, with full MM atomtype and charge, connections
                               default is taken from the .com/.gjf file that generates the .log/.out file
                               When .log/.out does not contain PDBinfo, -t must be given!
                               e.g. -t mod_rc1.com/.gjf
                           

	  Options for .log/.out extracting:

	  -s STEP              optional. Step number for geometry extraction, default is the last stationary point.
                       Examples: (OP: optimization; SP: stationary point; PN: IRC point number)
                           -s 1      >> extract the first OP geometry
                           -s 7      >> extract the 7th geometry, SP prioritize OP, if no SP, will be the 7th OP
                           -s lastpt >> extract the last geometry in OP list (Maybe not a SP)
                           -s 7p     >> extract the 7 OP geometry despites that there are more than 7 SPs
                       # The following examples extract multple geometries to multiple outputs:
                           -s all    >> extract all the geometries, SP/PN prioritize OP if there is any
                           -s allp   >> extract all the OP geometries, SP/PN included
	  -measure MEASURE     optional. Default is not print or print scanned coords if it is a scan job.
                           e.g. ":UPG@C1p-UNK@O5"
	  -mf MASK_FILE        optional. A mask file list.                               e.g. -mf masks.lst
	  -pkl                 optional. dump/load coordinates to/from a pickle (.pkl) file automatically
                           (with the same prefix as .log/.out file)
                           by given the -pkl option, will decide if it is dump or load
	  -qm                  optional. flag to extract the high-lyaer atoms (and linked hydrogen),
                           used with -o -s ; if -qm is given, -o outputs only the HL atoms
	  -v                   optional. verbose output, print more information.         e.g. -v

	  Options for getting the .com/.gjf output(s):

	  -linda               optional. If LindaWorkers are needed for parallelization, e.g. -l
	  -p PROCESSORS        optional. Number of CPUs, default is SameAsLog;           e.g. -p 32
	  -m MEM               optional. Memory size, default is SameAsLog;              e.g. -m 108GB
	  -r ROUTINE           optional. Statement of routine in a pair of quota, default is SameAsLog
                           e.g. "#T oniom(b3lyp/6-31g(d):amber=SoftOnly)\ngeom=connectivity"
                           e.g. -r template # will keep the routine from template .com/.gjf
	  -title TITLE         optional. Title section of the output com file, Default is SameAsLog
                           e.g.  --title "no title
                           e.g.  --title "template", get it from a template .com
	  -gencon              optional. Generate new connectivity for the first extracted geometry, 
                           if given, the "$GAUSS_EXEDIR/newzmat" will be used, please ensure it exists
                           e.g. -gencon
                           

	  Options for getting the .pdb output(s):

	  -shift RESNUM_SHIFT  optional. number of shift for restore the PDB residue number, default is 0
                           e.g. -shift 2

结果示例：


- **a.查看任务失败的opt优化信息:**
  
><jhli@beta ~/bin/jgot> *log_extract.py -g gmGlc_tsall_ircRC.log -ener -mf active_coords.mf -conv*
>
>The list logInfo has been read from gmGlc_tsall_ircRC.log
>Total number of geometries (OP): 448
>Warning: No Stationary/Path geometry found!
>OptimizationID,Low_model,High_model,Low_real,Extrapolated,:UNK:452@O5-:UPG:451@C1p,:UPG:451@C1p-:UPG:451@O3B,:UNK:452@O5-:UPG:451@C1p-:UPG:451@O3B,:UNK:452@H40-:UNK:452@O5,d_Max_Force,d_RMS_Force,d_Max_Displ,d_RMS_Displ,d_Max_MM_Force,d_RMS_MM_Force
>1,-0.141130657285,-3043.092223175412,-39.755405587340,-3082.706498105467,2.6813,1.6962,152.1920,1.0201,NO_0.027215,NO_0.002335,NO_0.176915,NO_0.028893,NO_0.005085
>2,-0.171740852221,-3043.100435509996,-39.787949926670,-3082.716644584444,2.7522,1.5698,151.7727,1.0178,NO_0.016415,NO_0.001834,NO_0.456080,NO_0.084267,NO_0.000045
>3,-0.180967966441,-3043.094359502330,-39.800000719274,-3082.713392255163,2.9217,1.3263,150.5819,1.0066,NO_0.096185,NO_0.004898,NO_0.193759,NO_0.044483,YES   
>4,-0.183036662582,-3043.104711843111,-39.802053959104,-3082.723729139633,2.8767,1.5190,149.4648,1.0095,NO_0.019133,NO_0.002324,NO_0.111319,NO_0.023062,NO_0.000475
>5,-0.186603726790,-3043.106147782981,-39.806251051752,-3082.725795107943,2.9146,1.4890,148.2524,1.0035,NO_0.019291,NO_0.001955,NO_0.196738,NO_0.036378,NO_0.000132
>6,-0.189458892785,-3043.106992912125,-39.807919486604,-3082.725453505944,2.9340,1.4790,147.8165,1.0018,NO_0.061589,NO_0.005258,NO_0.116321,NO_0.049401,YES   
>7,-0.189478298133,-3043.106359065329,-39.809528815007,-3082.726409582204,2.9614,1.4763,146.4037,1.0028,NO_0.029412,NO_0.002663,NO_0.146629,NO_0.034121,YES
>......
>442,-0.190281977052,-3043.107785716332,-39.815572992174,-3082.733076731455,2.9929,1.4686,140.0788,1.0001,NO_0.042955,NO_0.003668,YES    ,YES     ,YES     
>443,-0.190284357885,-3043.107787086557,-39.815577039091,-3082.733079767764,2.9929,1.4686,140.0794,1.0000,NO_0.042884,NO_0.003661,YES    ,YES     ,YES     
>444,-0.190286697014,-3043.107788432298,-39.815581013341,-3082.733082748625,2.9929,1.4686,140.0800,1.0000,NO_0.042797,NO_0.003653,YES    ,YES     ,YES     
>445,-0.190289010804,-3043.107789763321,-39.815584947013,-3082.733085699530,2.9929,1.4686,140.0806,1.0000,NO_0.042649,NO_0.003642,YES    ,YES     ,YES     
>446,-0.190291314781,-3043.107791089126,-39.815588872296,-3082.733088646642,2.9929,1.4686,140.0812,1.0000,NO_0.042663,NO_0.003642,YES    ,YES     ,YES     
>447,-0.190293611782,-3043.107792409389,-39.815592764480,-3082.733091562087,2.9929,1.4686,140.0818,1.0001,NO_0.042577,NO_0.003633,YES    ,YES     ,YES     
>448,NaN,NaN,NaN,NaN,2.9929,1.4686,140.0818,1.0001,1      ,        ,       ,        ,        ,2      , ....
>
active_coords.mf的内容如下:

>\# UNK : ligand residue name
>
>\# UPG : cofactor residue name
>
>:UNK@O5-:UPG@C1p          # CV1
>
>:UPG@C1p-:UPG@O3B         # leaving bond length 1
>
>:UNK@O5-:UPG@C1p-:UPG@O3B # CV2
>
>:UNK@H40-:UNK@O5          # leaving bond length 2
>


上述信息中，如第445个点，备注CV1的距离是1.4686 Å,备注CV2的角度是140.0806 °。从后面剩下的两个NO的值可以看到这个log文件只要提取最后一个坐标继续优化就应该能收敛了。



- **b. 提取opt的log文件的某一个点到新的输入文件（.com/.gjf)**
  
><jhli@beta ~/bin/jgot> log_extract.py -g gmGlc_tsall_ircRC.log -o gmGlc_tsall_ircRC_conti.com -s 445p -v
>
>Basic .log/.out information:
>
>GauVersion: Gaussian 16: ES64L-G16RevC.01  3-Jul-2019
>
>Link0: ['%chk=gmGlc_tsall_ircRC.chk', '%mem=108GB','%nprocshared=28']
>
>Routine: #T oniom(b3lyp/6-31g(d):amber=SoftFirst)opt=calcfc freq geom=connectivity
>
>Charge and multiplicity: 1 1 -1 1 -1 1
>
>modredundant info: NOREDUNDANT
>
>Number of atoms from log file:   12801
>
>Number of atoms in tmeplate com: 12801
>
>The list logInfo has been read fromgmGlc_tsall_ircRC.log
>
>Total number of geometries (OP): 448
>
>Warning: No Stationary/Path geometry found!
>
>Not all the geometries have energy record
>
>Step_ID 445p given, extracting the 445 thoptimization point
>
><jhli@beta ~/bin/jgot> grep -n b3lyp gmGlc_tsall_ircRC_conti.com
>
>4:#T oniom(b3lyp/6-31g(d):amber=SoftFirst)opt=calcfc freq geom=connectivity
>
><jhli@beta ~/bin/jgot> log_extract.py -g gmGlc_tsall_ircRC.log -o gmGlc_tsall_ircRC_conti.com -s 445p -r "#Toniom(b3lyp/6-31g(d):amber=SoftFirst) opt=(calcfc,maxstep=2,notrust) freqgeom=connectivity"
>
>The list logInfo has been read fromgmGlc_tsall_ircRC.log
>
>Total number of geometries (OP): 448
>
>Warning: No Stationary/Path geometry found!
>
>Step_ID 445p given, extracting the 445 thoptimization point
>
><jhli@beta ~/bin/jgot> grep -n b3lyp gmGlc_tsall_ircRC_conti.com
>
>4:#T oniom(b3lyp/6-31g(d):amber=SoftFirst)opt=(calcfc,maxstep=2,notrust) freq geom=connectivity
>
><jhli@beta ~/bin/jgot> log_extract.py -g gmGlc_tsall_ircRC.log -o gmGlc_tsall_ircRC_conti.pdb -s 445p
>
>The list logInfo has been read fromgmGlc_tsall_ircRC.log
>
>Total number of geometries (OP): 448
>
>Warning: No Stationary/Path geometry found!
>
>Step_ID 445p given,extracting the 445 th optimization point
>


- **c. 从opt优化好的log文件提取pdb**
  
><jhli@beta ~/bin/jgot> log_extract.py -g xmXyl_tsall_ircPC.log -o xmXyl_PC.pdb
>
>The list logInfo has been read from xmXyl_tsall_ircPC.log
>Total number of geometries (OP): 126
>Number of Stationary/Path points: 1
>The last Stationary geometry will be extracted
>
><jhli@beta ~/bin/jgot> head -4 xmXyl_PC.pdb
>ATOM      1  N  LPROX   1     -22.543  -5.247  22.387 -1.00-0.20           N
>ATOM      2  H2 LPRO X  1     -22.174  -4.532  21.775 -1.00 0.31           H
>ATOM      3  H3 LPRO X  1     -21.853  -5.557  23.056 -1.00 0.31           H
>ATOM      4  CD LPRO X  1     -23.733  -4.873  23.131 -1.00 -0.01          C
>
><jhli@beta ~/bin/jgot> log_extract.py -g xmXyl_tsall_ircPC.log -o xmXyl_PC.pdb -shift 3
>
>......
>
><jhli@beta ~/bin/jgot> head -4 xmXyl_PC.pdb
>
>ATOM      1  N  LPROX   4     -22.543  -5.247  22.387 -1.00-0.20           N
>ATOM      2  H2 LPRO X  4     -22.174  -4.532  21.775 -1.00 0.31           H
>ATOM      3  H3 LPRO X  4     -21.853  -5.557  23.056 -1.00 0.31           H
>ATOM     4  CD LPRO X   4     -23.733  -4.873 23.131 -1.00 -0.01           C
>


- **d. 快速查看柔性扫描任务的信息和提取结构**
  
><jhli@beta~/bin/jgot> log_extract.py -g gmApi_sc0.log -ener -conv -mf active_coords.mf
>
>The list logInfo has been read from gmApi_sc0.log
>Total number of geometries (OP): 192
>Number of Stationary/Path points: 8
>
>StationaryID,Low_model,High_model,Low_real,Extrapolated,:UNK:452@O5-:UPG:451@C1p,:UPG:451@C1p-:UPG:451@O3B,:UNK:452@O5-:UPG:451@C1p-:UPG:451@O3B,:UNK:452@H40-:UNK:452@O5,d_Max_Force,d_RMS_Force,d_Max_Displ,d_RMS_Displ,d_Max_MM_Force,d_RMS_MM_Force

>1,-0.050459278095,-2928.541423078239,-39.123555739087,-2967.614519539231,3.1039,1.4681,162.0988,1.0009,YES,YES ,YES ,YES ,YES
>2,-0.049817740458,-2928.541274987220,-39.122354168198,-2967.613811414960,3.0039,1.4714,161.7257,1.0012,YES,YES ,YES ,YES ,YES
>3,-0.048129948031,-2928.540657237018,-39.120114556347,-2967.612641845333,2.9039,1.4751,161.4714,1.0024,YES,YES ,YES ,YES ,YES
>4,-0.045041807195,-2928.539514840195,-39.116434829366,-2967.610907862365,2.8039,1.4795,161.4069,1.0035,YES,YES ,YES ,YES ,YES
>5,-0.040496414320,-2928.538506991298,-39.110526207568,-2967.608536784546,2.7039,1.4849,160.9737,1.0054,YES,YES ,YES ,YES ,YES
>6,-0.032587808584,-2928.535398987689,-39.102431173197,-2967.605242352303,2.6039,1.4920,161.4012,1.0060,YES,YES ,YES ,YES ,YES
>7,-0.020098158073,-2928.531168614695,-39.089627740066,-2967.600698196688,2.5039,1.5032,161.9769,1.0064,YES,YES ,YES ,YES ,YES
>8,-0.000329283282,-2928.525397024345,-39.069532066256,-2967.594599807319,2.4039,1.5186,162.7589,1.0068,YES,YES ,YES ,YES ,YES
>


- **e. 提取PC的QM区域的原子，并更新坐标的键连信息**
  
><jhli@beta~/bin/jgot> log_extract.py -g xmXyl_tsall_ircPC.log -o xmXyl_PC.qm.com -qm-gencon
>
>The list logInfo has been read from xmXyl_tsall_ircPC.log
>Total number of geometries (OP): 126
>Number of Stationary/Path points: 1
>The last Stationary geometry will be extracted


- **f. 查看IRC任务信息**
  
><jhli@beta ~/bin/jgot> log_extract.py -g xmXyl_tsall_irc.log  -ener -mf active_coords.mf
>
>The list logInfo has been read fromxmXyl_tsall_irc.log
>Total number of geometries (OP): 83
>Number of Stationary/Path points: 57
>This is an IRC log file
>The IRC information will be printed in the end
>
>StationaryID,Low_model,High_model,Low_real,Extrapolated,:UNK:452@O5-:UPG:451@C1p,:UPG:451@C1p-:UPG:451@O3B,:UNK:452@O5-:UPG:451@C1p-:UPG:451@O3B,:UNK:452@H40-:UNK:452@O5,#IRC_Step,point_number,path_direction
>
>1,0.777636138275,-2928.528467994788,-38.829465016202,-2968.135569149265,2.2072,2.4800,161.8979,1.4432,1,0,TS
>2,0.853861006521,-2928.528573163492,-38.753177365702,-2968.135611535716,2.1798,2.5045,161.8962,1.4684,2,1,FORWARD
>3,0.924542199921,-2928.528724784458,-38.682466861936,-2968.135733846315,2.1525,2.5293,161.8900,1.4868,3,2,FORWARD
>4,0.997123436066,-2928.528948275011,-38.609857077093,-2968.135928788170,2.1249,2.5536,161.8894,1.5038,4,3,FORWARD
>5,1.072981232576,-2928.529237016522,-38.533973526996,-2968.136191776094,2.0974,2.5778,161.8905,1.5199,5,4,FORWARD
>6,1.152452340458,-2928.529589663570,-38.454477928533,-2968.136519932561,2.0697,2.6016,161.8949,1.5349,6,5,FORWARD
>7,1.236999782726,-2928.530004719950,-38.369906306137,-2968.136910808813,2.0422,2.6251,161.8980,1.5495,7,6,FORWARD
>8,1.326516450095,-2928.530479176480,-38.280367414376,-2968.137363040952,2.0146,2.6483,161.9037,1.5630,8,7,FORWARD
>9,1.423937453447,-2928.531015904567,-38.182921697055,-2968.137875055069,1.9870,2.6712,161.9118,1.5764,9,8,FORWARD
>10,1.530135619925,-2928.531612828581,-38.076697779454,-2968.138446227960,1.9593,2.6938,161.9228,1.5892,10,9,FORWARD
>11,1.647882145177,-2928.532270528608,-37.958922511056,-2968.139075184841,1.9316,2.7161,161.9351,1.6020,11,10,FORWARD
>......
>28,8.185097856532,-2928.552350635909,-31.420015819673,-2968.157464312113,1.5614,3.0055,162.6574,1.7765,28,27,FORWARD
>29,8.841689008827,-2928.553777825276,-30.763567011569,-2968.159033845673,1.5544,3.0118,162.7209,1.7800,29,28,FORWARD
>30,0.707094790506,-2928.528454504224,-38.900067670351,-2968.135616965081,2.2347,2.4555,161.8967,1.4179,30,1,REVERSE
>31,0.626690367818,-2928.528503684380,-38.980592937792,-2968.135786989989,2.2615,2.4319,161.8949,1.3795,31,2,REVERSE
>32,0.536462994184,-2928.528676900208,-39.071064278291,-2968.136204172682,2.2851,2.4120,161.8802,1.3174,32,3,REVERSE
>33,0.457016247082,-2928.529428793143,-39.150883488689,-2968.137328528913,2.3019,2.3985,161.8608,1.2321,33,4,REVERSE
>......
>54,-0.083329803418,-2928.555947062072,-39.691704943813,-2968.164322202468,2.6876,1.8244,163.1090,1.0139,54,25,REVERSE
>55,-0.096794171130,-2928.557488260972,-39.705130723363,-2968.165824813205,2.7057,1.7937,163.1213,1.0131,55,26,REVERSE
>56,-0.109105485375,-2928.559048317587,-39.717401006499,-2968.167343838711,2.7237,1.7632,163.1308,1.0124,56,27,REVERSE
>57,-0.120247697075,-2928.560614241086,-39.728499888268,-2968.168866432279,2.7414,1.7330,163.1365,1.0117,57,28,REVERSE
>


- **g. 根据IRC的方向和反应坐标的值提取RC和PC的初始结构进行下一步的优化**
  
><jhli@beta ~/bin/jgot> log_extract.py -g xmXyl_tsall_irc.log -o xmXyl_irc_RC_opt.com -s 57 -r "#Toniom(b3lyp/6-31g(d):amber=SoftFirst) opt=(calcfc,maxstep=2,notrust) freqgeom=connectivity"
>
>The list logInfo has been read from xmXyl_tsall_irc.log
>Total number of geometries (OP): 83
>Number of Stationary/Path points: 57
>Step_ID 57 given, extracting the 57 th stationary/IRC point
>The 28 th point in the REVERSE path will be extracted
>This geometry has 0 energy records, if it's not the last, please check!
>
><jhli@beta ~/bin/jgot> log_extract.py -g xmXyl_tsall_irc.log  -o xmXyl_irc_PC_opt.com -s 29 -r "#Toniom(b3lyp/6-31g(d):amber=SoftFirst) opt=(calcfc,maxstep=2,notrust) freqgeom=connectivity"
>The list logInfo has been read fromxmXyl_tsall_irc.log
>Total number of geometries (OP): 83
>Number of Stationary/Path points: 57
>Step_ID 29 given, extracting the 29 thstationary/IRC point
>The 28 th point in the FORWARD path will beextracted
>This geometry has 0 energyrecords, if it's not the last, please check!
>



## get_input_info.py

这个脚本就是快速查看ONIOM输入文件的信息，可以结合-mf选项读取mask（amber格式）文件，主要检查一下反应坐标、扫描坐标等是否正确。



## charge2zero4SP.py

这个脚本的功能是根据给定的非QM区域残基的列表批量生成某个残基电荷为0的EMB单点能计算输入文件，目的比较鸡肋。用法如下:


	<jhli@beta ~/bin/jgot> charge2zero4SP.py -h
	usage: charge2zero4SP.py [-h] -i INPUT_COM -z ZEROCHARGES -o OUTPUT_COM
                         [-p PROCESSORS] [-m MEM] [-l] [-r ROUTINE]
                         [--title TITLE]

	Change the charges to zero for a set of residues in an ONIOM input

	optional arguments:
	  -h, --help            show this help message and exit
	  -i INPUT_COM, --input_com INPUT_COM
                            A ONIOM input (.com) that has correct atom types and charges
                            str, Required, e.g. -gv embBS2_ts_clac_IRC_RC.com
                            
	  -z ZEROCHARGES, --zerocharges ZEROCHARGES
                            A list of residues that need to change to zero charge 1 by 1,
                                   format example: "HID 16"
                            str, Required, e.g. -z zero_charge.lst
                            
	  -o OUTPUT_COM, --output_com OUTPUT_COM
                            output .com basename, e.g. -o embBS2_TS  Required option
                            
	  -p PROCESSORS, --processors PROCESSORS
                            number of CPUs, e.g. -p 32
	  -m MEM, --mem MEM         memory size, e.g. -m 64GB
	  -l, --lindaworker         if LindaWorkers are needed for parallelization, -e.g. -l
	  -r ROUTINE, --routine ROUTINE
                            Statement of routine in a pair of quota,
                            e.g. "#T oniom(b3lyp/6-31g(d):amber=SoftFirst)=EmbedCharge geom=connectivity"
                            e.g. -r keep # will keep the routine used in input .com
                            
	  --title TITLE             Title section of the output com file, e.g. --title "no title"
                                   Or keep it from template: --title keep

目前还没有频率分析的功能，单点SP任务的geometry提取有bug，以后有需要再加入和更改。

