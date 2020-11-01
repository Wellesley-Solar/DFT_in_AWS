import numpy as np
import re
np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:10f}'.format})
def extract_ABX(chem_formula):
    """
    
    Args:
        chem_formula: perovskite composition in (A)(B)(X)3 format. 
          Examples include (Cs)(Pb0.5Sn0.5)(I)3
        
    Returns:
        A nested list whose first list contains tuples of the A elements 
        and their percentage, second list tuples of B elements and their 
        percentages, and third list tuples of C elements and their 
        percentages.
        
        [[(Cs, 1)], [(Pb, 0.5), (Sn, 0.5)], [(I, 1)]]
    
    Raises:
        None
    
    <<< extract_ABX('(Cs.75Fa.1Ma.14Rb.01)(Pb.25Sn.75)(I.33Br.66Cl.01)3')
    
    >>> [[('Cs', '.75'), ('Fa', '.1'), ('Ma', '.14'), ('Rb', '.01')],
         [('Pb', '.25'), ('Sn', '.75')],
         [('I', '.33'), ('Br', '.66'), ('Cl', '.01')]]
    
    """
    temp_formula = chem_formula
    ABX3=[]
    for i in range(3):
        Ei = temp_formula.find('(')
        Ef = temp_formula.find(')')
        ABX3.append(temp_formula[Ei+1:Ef])
        temp_formula=temp_formula[Ef+1:]
        
    for i,v in enumerate(ABX3):
        peri0 = v.find('.')
        if peri0 == -1:
            ABX3[i]=[(v,'1')]
            continue   # if there is only one element for a ion position, continue to the next ion/end loop
            
        chem_index=None
        for p,c in enumerate(v[peri0:]):
            if str.isupper(c):
                chem_index=p+peri0
                continue
        

        ABX3[i] = separate_elements(v)#ABX3_element#[(temp_element,temp_perc),(v[chem_index:chem_index+peri1],perc)]
        
    
    return ABX3


def separate_elements(subchem_formula, debug=False):
    """
    
    Args:
        subchem_formula: a string specifying the percentage of elements
        in a given A, B, or X, such as 'Pb0.5Sn0.5'
        
    Returns:
        A nested list whose first list contains tuples of the A elements 
        and their percentage, second list tuples of B elements and their 
        percentages, and third list tuples of C elements and their 
        percentages.

    
    Raises:
        None
    
    <<< separate_elements('Cs.75Fa.1Ma.14Rb.01')
    
    >>> [('Cs', '.75'), ('Fa', '.1'), ('Ma', '.14'), ('Rb', '.01')]
    """
    accum_index=0
    elements_list=[]
    perc_list=[]
    temp_formula=subchem_formula
    
    while True:
        temp_zero = temp_formula.find('.',1)
        
        if debug: print('temp zero',temp_zero)
        # when there are no more zero,stop searching for additional elements
        if temp_zero == -1:
            perc_list.append(temp_formula)
            break
        for p,c in enumerate(temp_formula):
            if str.isupper(c):
                accum_index=p
                break
        if debug: print('accum_index',accum_index)
        percentage = temp_formula[:accum_index]
        if percentage!='':
            perc_list.append(percentage)

        elements_list.append(temp_formula[accum_index:temp_zero])
        element = temp_formula[temp_zero:]

        temp_formula = temp_formula[temp_zero:]
        if debug: print('end loop temp_formula',temp_formula)
        if debug: print("percentage list",perc_list,'element list',elements_list)

    if debug: print(elements_list, perc_list)

    ABX3_element=[(e, p) for e,p in zip(elements_list, perc_list)]
    
    perc_sum = sum([float(percent) for percent in perc_list])
    if perc_sum != 1:
        print('Warning: the composition', subchem_formula, 'do not sum up to one. (', perc_sum,'instead)')
    if debug: print(ABX3_element)
    return ABX3_element


def sqs_percent(ABX3_list):
    """
    
    Args:
        subchem_formula: a string specifying the percentage of elements in a given A,
        B, or X, such as 'Pb0.5Sn0.5'
        
    Returns:
        A nested list whose first list contains tuples of the A elements 
        and their percentage, second list tuples of B elements and their 
        percentages, and third list tuples of C elements and their 
        percentages.
    
    Raises:
        None
    
    <<< separate_elements('Cs.75Fa.1Ma.14Rb.01')
    
    >>> [('Cs', '.75'), ('Fa', '.1'), ('Ma', '.14'), ('Rb', '.01')]
    """
    sqs_list=[]
    for row in ABX3_list:
        first=True
        temp=''
        for pair in row:
            if not first: 
                temp +=', '
            if pair[1]=='1':
                temp += pair[0]
            else:
                temp += pair[0]+'='+pair[1]
            first=False
        sqs_list.append(temp)
    return sqs_list


def removeChar(myString, remove="()[]"):
    """
    Remove parentatheses in a string
    
    Arguement: 
        myString: a string of type str
        
    Return:
        a string with parentatheses removed
    """
    tempStr = myString
    for each in remove:
        # print(type(each))
        tempStr=tempStr.replace(str(each),"")
    return tempStr


class Molecules:
    def __init__(self, chem):
        self.db={
            "FA":[
                ("C", np.array([0.000001125,0.06614,0])),
                ("N", np.array([0.181743125,-0.02871,0])),
                ("N", np.array([-0.181744875,-0.02871,0])),
                ("H", np.array([0.000000125,0.241308,0])),
                ("H", np.array([0.311782125,0.066036,0])),
                ("H", np.array([0.201975125,-0.19105,0])),
                ("H", np.array([-0.201974875,-0.19105,0])),
                ("H", np.array([-0.311781875,0.066036,0]))                
            ]
        }
        self.chem = chem
        self.coor = self.db[chem]


class Control:
    def __init__(self, calculation="'scf'", outdir="'./out'",prefix="'none'",\
                 pseudo_dir="'../../pseudo'", verbosity="'high'", tprnfor=".true.",\
                 tstress=".true.",etot_conv_thr="5.0000000000d-05",\
                 forc_conv_thr="1.0000000000d-04", **kwargs):
        """
        
        """
        self.block = "&CONTROL"
        self.params={'calculation':calculation, 'outdir':outdir, 'prefix':prefix,\
                       'pseudo_dir':pseudo_dir, 'verbosity':verbosity, 'tprnfor':tprnfor,\
                       'tstress':tstress,'etot_conv_thr':etot_conv_thr,\
                       'forc_conv_thr':forc_conv_thr, **kwargs
                      }
        
        
        
class System:
    def __init__(self, ibrav="0",nat="5",ntyp="3",ecutwfc="0",ecutrho="0",\
                input_dft="'pbe'",occupations="'smearing'",smearing="'cold'",\
                degauss="0.005d0", noncolin = ".true.",**kwargs):
        self.block = "&SYSTEM"
        self.params={'ibrav':ibrav, 'nat':nat, 'ntyp':ntyp,\
                       'ecutwfc':ecutwfc, 'ecutrho':ecutrho, 'input_dft':input_dft,\
                       'occupations':occupations,'smearing':smearing,\
                       'degauss':degauss, "noncolin":noncolin, **kwargs
                    }
        
class Electrons:
    def __init__(self, conv_thr="1d-08", mixing_beta="0.7d0", electron_maxstep="500", **kwargs):
        """
        
        """
        self.block = "&ELECTRONS"
        self.params={
            'conv_thr':conv_thr, 'mixing_beta':mixing_beta,\
            'electron_maxstep':electron_maxstep, **kwargs
        }

        
        
class Ions:
    def __init__(self,ion_dynamics="'bfgs'", **kwargs):
        self.block = "&IONS"
        self.params={
            'ion_dynamics':ion_dynamics, **kwargs
        }
        
        
class Cell:
    def __init__(self,cell_dynamics="'bfgs'", press="0.1d0", press_conv_thr="0.5d0", **kwargs):
        self.block = "&CELL"
        self.params={
            'cell_dynamics':cell_dynamics, 'press':press,\
            'press_conv_thr':press_conv_thr, **kwargs
        }
        
        
class sqs2QE():
    def __init__(self, fullrel=False, celldim=1):
        """
        Args:
            fullrel: boolean argument selecting 'fully relativistic PP' when true and
                 'scalar relativistic PP' when false
            *args: custom objects specifying the contorl parameters of QE calculation
                   Available class options: {Control, System, Electrons, Ions, Cells}

        Returns:
            None
        """
        self.params = None
        self.dirpath = None
        self.sqsLines = None
        self.atoms = None
        self.elements = []
        self.mass = {
            "Cs": 132.90545, "Pb": 207.20000, "I": 126.90447,\
            "Br": 79.90400, "Cl": 35.45, "Sn": 10086,\
            "C": 12.011, "N": 14.007, "H": 1.008
        }

        # self.pseudo = {
        #     "scalar":{
        #         "Cs": "Cs_pbe_v1.uspp.F.UPF", "Pb": "Pb.pbe-dn-kjpaw_psl.0.2.2.UPF",\
        #         "I": "I.pbe-n-kjpaw_psl.0.2.UPF", "Br": "br_pbe_v1.4.uspp.F.UPF",\
        #         "Cl": "Cl.pbe-n-rrkjus_psl.1.0.0.UPF","Sn": "Sn_not_relativistic",\
        #         "C": "C.pbe-n-kjpaw_psl.1.0.0.UPF", "H": "H_ONCV_PBE-1.0.oncvpsp.upf",\
        #         "N": "N.oncvpsp.upf"
        #     },\
        #     "rel":{
        #         "Cs": "Cs.rel-pbe-spnl-rrkjus_psl.1.0.0.UPF", "Pb": "Pb.rel-pbe-dn-rrkjus_psl.1.0.0.UPF",\
        #         "I": "I.rel-pbe-n-rrkjus_psl.1.0.0.UPF", "Br": "Br.rel-pbe-n-rrkjus_psl.1.0.0.UPF",\
        #         "Cl": "Cl.rel-pbe-n-rrkjus_psl.1.0.0.UPF", "Sn": "Sn_very_relativistic",\
        #         "C": "C.rel-pbe-n-rrkjus_psl.1.0.0.UPF", "H": "H.rel-pbe-rrkjus_psl.1.0.0.UPF",\
        #         "N": "N.rel-pbe-n-rrkjus_psl.1.0.0.UPF"
        #     }
        #               }

        self.precision={
            "Cs": {
                "cutoff": 30.0,
                "dual": 8.0,
                "filename": "Cs_pbe_v1.uspp.F.UPF",
                "md5": "3476d69cb178dfad3ffaa59df4e07ca4",
                "pseudopotential": "GBRV-1.2",
                "rho_cutoff": 240.0
            },
            "Pb": {
                "cutoff": 45.0,
                "dual": 8.0,
                "filename": "Pb.pbe-dn-kjpaw_psl.0.2.2.UPF",
                "md5": "9d431e6316058b74ade52399a6cf67da",
                "pseudopotential": "031PAW",
                "rho_cutoff": 360.0
            },
            "I": {
                "cutoff": 45.0,
                "dual": 8.0,
                "filename": "I.pbe-n-kjpaw_psl.0.2.UPF",
                "md5": "d4ef18d9c8f18dc85e5843bca1e50dc0",
                "pseudopotential": "031PAW",
                "rho_cutoff": 360.0
            },
            "Br": {
                "cutoff": 90.0,
                "dual": 8.0,
                "filename": "br_pbe_v1.4.uspp.F.UPF",
                "md5": "d3ffb7b29f6225aa16fe06858fb2a80b",
                "pseudopotential": "GBRV-1.4",
                "rho_cutoff": 720.0
            },
            "Cl": {
                "cutoff": 100.0,
                "dual": 8.0,
                "filename": "Cl.pbe-n-rrkjus_psl.1.0.0.UPF",
                "md5": "18cc83b4be324290a879bee3176034ba",
                "pseudopotential": "100US",
                "rho_cutoff": 800.0
            },
            "C": {
                "cutoff": 45.0,
                "dual": 8.0,
                "filename": "C.pbe-n-kjpaw_psl.1.0.0.UPF",
                "md5": "5d2aebdfa2cae82b50a7e79e9516da0f",
                "pseudopotential": "100PAW",
                "rho_cutoff": 360.0
            },
            "H": {
                "cutoff": 80.0,
                "dual": 4.0,
                "filename": "H_ONCV_PBE-1.0.oncvpsp.upf",
                "md5": "1790becc920ee074925cf490c71280fe",
                "pseudopotential": "SG15",
                "rho_cutoff": 320.0
            },
            "N": {
                "cutoff": 80.0,
                "dual": 4.0,
                "filename": "N.oncvpsp.upf",
                "md5": "563d65bfb082928f0c9eb97172f6c357",
                "pseudopotential": "Dojo",
                "rho_cutoff": 320.0
            },
        "Sn": {
                "cutoff": 70.0,
                "dual": 8.0,
                "filename": "Sn_pbe_v1.uspp.F.UPF",
                "md5": "4cf58ce39ec5d5d420df3dd08604eb00",
                "pseudopotential": "GBRV-1.2",
                "rho_cutoff": 560.0
            },
        }

        self.efficiency={
            "Cs": {
                "cutoff": 30.0,
                "dual": 8.0,
                "filename": "Cs_pbe_v1.uspp.F.UPF",
                "md5": "3476d69cb178dfad3ffaa59df4e07ca4",
                "pseudopotential": "GBRV-1.2",
                "rho_cutoff": 240.0
            },
            "C": {
                "cutoff": 45.0,
                "dual": 8.0,
                "filename": "C.pbe-n-kjpaw_psl.1.0.0.UPF",
                "md5": "5d2aebdfa2cae82b50a7e79e9516da0f",
                "pseudopotential": "100PAW",
                "rho_cutoff": 360.0
            },
            "N": {
                "cutoff": 60.0,
                "dual": 8.0,
                "filename": "N.pbe-n-radius_5.UPF",
                "md5": "16739722b17309cd8fe442a2ace49922",
                "pseudopotential": "THEOS",
                "rho_cutoff": 480.0
            },
            "H": {
                "cutoff": 60.0,
                "dual": 8.0,
                "filename": "H.pbe-rrkjus_psl.1.0.0.UPF",
                "md5": "f52b6d4d1c606e5624b1dc7b2218f220",
                "pseudopotential": "100US",
                "rho_cutoff": 480.0
            },
            "I": {
                "cutoff": 35.0,
                "dual": 8.0,
                "filename": "I.pbe-n-kjpaw_psl.0.2.UPF",
                "md5": "d4ef18d9c8f18dc85e5843bca1e50dc0",
                "pseudopotential": "031PAW",
                "rho_cutoff": 280.0
            },
            "Br": {
                "cutoff": 30.0,
                "dual": 8.0,
                "filename": "br_pbe_v1.4.uspp.F.UPF",
                "md5": "d3ffb7b29f6225aa16fe06858fb2a80b",
                "pseudopotential": "GBRV-1.4",
                "rho_cutoff": 240.0
            },
            "Cl": {
                "cutoff": 40.0,
                "dual": 8.0,
                "filename": "cl_pbe_v1.4.uspp.F.UPF",
                "md5": "fc6f6913ecf08c9257cb748ef0700058",
                "pseudopotential": "GBRV-1.4",
                "rho_cutoff": 320.0
            },
            "Pb": {
                "cutoff": 40.0,
                "dual": 8.0,
                "filename": "Pb.pbe-dn-kjpaw_psl.0.2.2.UPF",
                "md5": "9d431e6316058b74ade52399a6cf67da",
                "pseudopotential": "031PAW",
                "rho_cutoff": 320.0
            },
        }

        if fullrel:
            self.pseudo="Not Valid"
        else:
            self.pseudo=self.efficiency
        
    def load_sqs(self, dirpath=".", isqs='bestsqs0.out'):
        """
        Reads the lines in an sqs output file
        
        Args:
            path: file path to the sqs output file from ATAT package

        Returns:
            None
        """
        self.dirpath=dirpath
        self.sqsf = open(self.dirpath+"/"+isqs)
        self.sqsLines = self.sqsf.readlines()
        self.sqsf.close()
        
        
    def format_chem2QE(self, line):
        """
        Converts the atomic position in ATAT format to QE format
        If molecules such as FA were given, it will convert it to its constituents
        elements.
        
        Args:
            line: a line specifying the element names and its coordinate position in the ATAT
                  package format

        Returns:
            A (list of) string containing the element names and its coordinate position.
        """
        temp_line = line.replace('\n','')
        element_index=None
        for i,c in enumerate(temp_line):
            if c.isupper():
                element_index=i
                break
        QEline = []
        
        if temp_line[element_index:].find("FA") != -1: 
            # Split "FA" into constituent components
            A_FA = Molecules("FA")
            coord_FA = np.array([2.0,2.0,2.0])-np.array([float(coor) for coor in temp_line[:element_index-1].split(' ')])
            for element in A_FA.coor:
                if element[0] not in self.elements: self.elements.append(element[0])
                # print(' '.join(map(str, element[1]+coord_FA)))
                # QEline.append(element[0]+" "+' '.join(map(str, element[1]+coord_FA)))
                actual_coord=(element[1]+coord_FA)
                actual_coord=actual_coord/np.array([2.0,2.0,2.0])
                QEline.append(element[0]+" "+removeChar(np.array2string(actual_coord),"[]"))
                # QEline.append(element[0]+" "+removeChar(np.array2string(element[1]+coord_FA),"[]"))
        else:
            tempElement= temp_line[element_index:]
            coord = np.array([2.0,2.0,2.0])-np.array([float(coor) for coor in temp_line[:element_index-1].split(' ')])
            coord = coord/np.array([2.0,2.0,2.0])
            QEline.append(tempElement+" "+removeChar(np.array2string(coord),"[]"))
            if tempElement not in self.elements: self.elements.append(tempElement)
            # QEline = tempElement+' '+temp_line[:element_index-1]

        return QEline

    def format_chem2GPAW(self, line):
        """
        Converts the atomic position in ATAT format to QE format
        If molecules such as FA were given, it will convert it to its constituents
        elements.
        
        Args:
            line: a line specifying the element names and its coordinate position in the ATAT
                  package format

        Returns:
            A (list of) string containing the element names and its coordinate position.
        """
        temp_line = line.replace('\n','')
        element_index=None
        for i,c in enumerate(temp_line):
            if c.isupper():
                element_index=i
                break
        GPAWline = []
        
        if temp_line[element_index:].find("FA") != -1: 
            # Split "FA" into constituent components
            A_FA = Molecules("FA")
            coord_FA = np.array([2.0,2.0,2.0])-np.array([float(coor) for coor in temp_line[:element_index-1].split(' ')])
            for element in A_FA.coor:
                if element[0] not in self.elements: self.elements.append(element[0])
                # print(' '.join(map(str, element[1]+coord_FA)))
                # QEline.append(element[0]+" "+' '.join(map(str, element[1]+coord_FA)))
                actual_coord=(element[1]+coord_FA)
                actual_coord=actual_coord/np.array([2.0,2.0,2.0])
                GPAWline.append((element[0], actual_coord))
                # QEline.append(element[0]+" "+removeChar(np.array2string(element[1]+coord_FA),"[]"))
        else:
            tempElement= temp_line[element_index:]
            coord = np.array([2.0,2.0,2.0])-np.array([float(coor) for coor in temp_line[:element_index-1].split(' ')])
            coord = coord/np.array([2.0,2.0,2.0])
            GPAWline.append((tempElement,coord))
            if tempElement not in self.elements: self.elements.append(tempElement)
            # QEline = tempElement+' '+temp_line[:element_index-1]

        return GPAWline 
    
        
    def toMatrix(self, lines, start_index):
        """
        A (list of) string containing the element names and its coordinate position
        
        Args:
            path: file path to the sqs output file from ATAT package

        Returns:
            None
        """
        return np.matrix(lines[start_index].replace("\n",";")+\
                         lines[start_index+1].replace("\n",";")+
                         lines[start_index+2].replace("\n",""))

    def sqschem2QE(self):
        """
        Convert all lines of atomic positions from the SQS output in ATAT format
        to QE format and store at self.atom
        
        Args:
            None

        Returns:
            None
        """
        self.atoms = [self.format_chem2QE(i) for i in self.sqsLines[6:]]

    def sqschem2GPAW(self):
        """
        Convert all lines of atomic positions from the SQS output in ATAT format
        to QE format and store at self.atom
        
        Args:
            None

        Returns:
            None
        """
        self.atoms = [self.format_chem2GPAW(i) for i in self.sqsLines[6:]]
        
    
    def flattenList(self, aList, outList=None):
        """
        Recursively flatten a list of QE atomic position in case nested lists
        are created by converting molecules
        
        Args:
            None

        Returns:
            None
        """
        if outList==None:
            outList=[]
        for each in aList:
            if type(each)==list:
                self.flattenList(each, outList)
            else:
                outList.append(each)
        return outList
    
    def load_scf(self):
        f = open(self.dirpath+"/"+self.dirpath+"scf.qeout",'r')
        fermi_val = np.zeros(0)
        for i in f:
            if "Kohn-Sham states" in i:
                fermi_line=re.findall(r"[-+]?\d*\.\d+|\d+", i)
                self.KSnum=float(fermi_line[0])
                break
        f.close()
        return None

    def load_vcrelax(self):
        f = open(self.dirpath+"/"+self.dirpath+"vc-relax.qeout",'r')
        coor=[]
        addc=False
        for i,v in enumerate(f):
            if "Begin final coordinates" in v:
                addc=True
            elif "End final coordinates" in v:
                addc=False
                break
            if addc:
                coor.append(v)
        f.close()
        celli=0
        for i,v in enumerate(coor):
            if "CELL_PARAMETERS" in v:
                print(v)
                celli=i
                break
        self.vccellshape=coor[celli:celli+4]
        # print(self.vccellshape)

        atomi=0
        for i,v in enumerate(coor):
            if "ATOMIC_POSITIONS" in v:
                atomi=i
                break
        self.vcatoms=coor[atomi:]
        # print(self.atoms)
        return None

    def writeQE(self, calc="scf", fromsqs=True):
        """
        Write all the control information and atomic information into
        a Quantum-ESPRESSO input file
        
        Args:
            chem: perovskite composition in (A)(B)(X)3 format. 
                  Examples include (Cs)(Pb0.5Sn0.5)(I)3

        Returns:
            None
        """
        self.calc=calc
        self.dirpath = self.dirpath
        print(self.dirpath)
        self.sqschem2QE()
        # flatten the list in case pseudo ions like "FA" are replaced wit a list of their constituents
        self.atoms = self.flattenList(self.atoms)
        
        # Configure DFT Controls appropriately
        self.wfc=0
        self.rho=0
        # if self.params[1]
        for each in self.elements:
            print(type(self.pseudo[each]["cutoff"]))
            if self.pseudo[each]["cutoff"]>self.wfc:
                self.wfc=self.pseudo[each]["cutoff"]
            if self.pseudo[each]["rho_cutoff"]>self.rho:
                self.rho=self.pseudo[each]["rho_cutoff"]
        self.Con = Control(calculation="'"+calc+"'",prefix="'"+self.dirpath+"'")
        self.Sys = System(ecutwfc=str(int(self.wfc)),ecutrho=str(int(self.rho)), \
                          ntyp=str(len(self.elements)),nat=str(len(self.atoms)))
        if calc == "bands":
            self.Sys = System(ecutwfc=str(int(self.wfc)),ecutrho=str(int(self.rho)), \
                            ntyp=str(len(self.elements)),nat=str(len(self.atoms)),\
                              nbnd=str(self.KSnum))
        self.Eon = Electrons()
        self.Ions = Ions()
        self.Cell = Cell()

        if calc == "scf" or calc == "bands":
            self.params = [self.Con, self.Sys, self.Eon]
        elif calc == "relax":
            self.params = [self.Con, self.Sys, self.Eon, self.Ions]
        elif calc == "vc-relax":
            self.params = [self.Con, self.Sys, self.Eon, self.Ions, self.Cell]
        
        # write all available Control blocks
        self.qefpath=self.dirpath+calc+'.qein'
        self.qef = open(self.dirpath+"/"+self.qefpath,"w+")
        for param in self.params:
            self.qef.write(param.block+"\n")
            for v in param.params:
                self.qef.writelines("  "+v+"="+param.params[v]+",\n")
            self.qef.write("/  \n\n")           
        
        # Write Cell parameter lines
        self.coor = self.toMatrix(self.sqsLines, 0)
        self.lvec = self.toMatrix(self.sqsLines, 3)
        self.cellshape = self.coor * self.lvec
        cells = None
        nline = "\n"
        if self.calc == "vc-relax":
            self.qef.write("CELL_PARAMETERS angstrom\n")
            cells=self.cellshape
        elif self.calc == "scf" or self.calc=="bands":
            cells=self.vccellshape
            nline=""
        
        for line in cells:
            self.qef.write(str(line).replace("[[","").replace("]]","")+nline)
        self.qef.write("\n")
        

        # Convert atomic position from SQS to QE
        self.qef.write("ATOMIC_SPECIES\n")
        for each in self.elements:
            self.qef.write("  "+each+"   "+str(self.mass[each])+"  "+self.pseudo[each]["filename"]+'\n')
        self.qef.write("\n")
        
        atom=None
        nline="\n"
        if self.calc == "vc-relax":
            self.qef.write("ATOMIC_POSITIONS crystal\n")
            atoms=self.atoms
        elif self.calc == "scf" or self.calc=="bands":
            # do nothing
            atoms=self.vcatoms
            nline=""
        # self.atoms.sort()
        for line in atoms:
            self.qef.write(line+nline)
        self.qef.write("\n")
        
        ### Writing the K Point Mesh
        self.qef.write("K_POINTS automatic\n")
        self.qef.write("4 4 4 0 0 0\n\n")
        
        # close qe.in file
        self.qef.close()


    def writeGPAW(self, calc="scf", fromsqs=True):
        """
        Write all the control information and atomic information into
        a Quantum-ESPRESSO input file
        
        Args:
            chem: perovskite composition in (A)(B)(X)3 format. 
                  Examples include (Cs)(Pb0.5Sn0.5)(I)3

        Returns:
            None
        """
        # Write Cell parameter lines
        self.coor = self.toMatrix(self.sqsLines, 0)
        self.lvec = self.toMatrix(self.sqsLines, 3)
        self.gcell = self.coor*self.lvec

        self.calc=calc
        self.dirpath = self.dirpath
        print(self.dirpath)
        self.sqschem2GPAW()
        # flatten the list in case pseudo ions like "FA" are replaced wit a list of their constituents
        self.atoms = self.flattenList(self.atoms)
        self.atoms.sort(key = lambda x: x[0])
        
        # make symbol for 
        self.elements.sort()
        self.symlist=[[element, 0] for element in self.elements]
        for atom in self.atoms:
        #     print(atom)
            for each in self.symlist:
                if atom[0]==each[0]:
                    each[1]=each[1]+1
        self.symstr=""
        for each in self.symlist:
            self.symstr=self.symstr+each[0]+str(each[1])

        self.pos=[]
        for atom in self.atoms:
            self.pos.append(np.asarray(atom[1]*self.gcell)[0].tolist())

        # self.pos = self.coor
        # # Configure DFT Controls appropriately
        # self.wfc=0
        # self.rho=0
        # # if self.params[1]
        # for each in self.elements:
        #     print(type(self.pseudo[each]["cutoff"]))
        #     if self.pseudo[each]["cutoff"]>self.wfc:
        #         self.wfc=self.pseudo[each]["cutoff"]
        #     if self.pseudo[each]["rho_cutoff"]>self.rho:
        #         self.rho=self.pseudo[each]["rho_cutoff"]
        # self.Con = Control(calculation="'"+calc+"'",prefix="'"+self.dirpath+"'")
        # self.Sys = System(ecutwfc=str(int(self.wfc)),ecutrho=str(int(self.rho)), \
        #                   ntyp=str(len(self.elements)),nat=str(len(self.atoms)))
        # if calc == "bands":
        #     self.Sys = System(ecutwfc=str(int(self.wfc)),ecutrho=str(int(self.rho)), \
        #                     ntyp=str(len(self.elements)),nat=str(len(self.atoms)),\
        #                       nbnd=str(self.KSnum))
        # self.Eon = Electrons()
        # self.Ions = Ions()
        # self.Cell = Cell()

        # if calc == "scf" or calc == "bands":
        #     self.params = [self.Con, self.Sys, self.Eon]
        # elif calc == "relax":
        #     self.params = [self.Con, self.Sys, self.Eon, self.Ions]
        # elif calc == "vc-relax":
        #     self.params = [self.Con, self.Sys, self.Eon, self.Ions, self.Cell]
        
        # # write all available Control blocks
        # self.qefpath=self.dirpath+calc+'.qein'
        # self.qef = open(self.dirpath+"/"+self.qefpath,"w+")
        # for param in self.params:
        #     self.qef.write(param.block+"\n")
        #     for v in param.params:
        #         self.qef.writelines("  "+v+"="+param.params[v]+",\n")
        #     self.qef.write("/  \n\n")           
        
        # # Write Cell parameter lines
        # self.coor = self.toMatrix(self.sqsLines, 0)
        # self.lvec = self.toMatrix(self.sqsLines, 3)
        # self.cellshape = self.coor * self.lvec
        # cells = None
        # nline = "\n"
        # if self.calc == "vc-relax":
        #     self.qef.write("CELL_PARAMETERS angstrom\n")
        #     cells=self.cellshape
        # elif self.calc == "scf" or self.calc=="bands":
        #     cells=self.vccellshape
        #     nline=""
        
        # for line in cells:
        #     self.qef.write(str(line).replace("[[","").replace("]]","")+nline)
        # self.qef.write("\n")
        

        # # Convert atomic position from SQS to QE
        # self.qef.write("ATOMIC_SPECIES\n")
        # for each in self.elements:
        #     self.qef.write("  "+each+"   "+str(self.mass[each])+"  "+self.pseudo[each]["filename"]+'\n')
        # self.qef.write("\n")
        
        # atom=None
        # nline="\n"
        # if self.calc == "vc-relax":
        #     self.qef.write("ATOMIC_POSITIONS crystal\n")
        #     atoms=self.atoms
        # elif self.calc == "scf" or self.calc=="bands":
        #     # do nothing
        #     atoms=self.vcatoms
        #     nline=""
        # # self.atoms.sort()
        # for line in atoms:
        #     self.qef.write(line+nline)
        # self.qef.write("\n")
        
        # ### Writing the K Point Mesh
        # self.qef.write("K_POINTS automatic\n")
        # self.qef.write("4 4 4 0 0 0\n\n")
        
        # # close qe.in file
        # self.qef.close()
