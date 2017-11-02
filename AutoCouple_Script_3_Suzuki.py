## Python Scripts for AutoCouple
##
## Author:  Laurent Batiste
##
## Affiliation:  A. Caflisch' group at the Department of Biochemistry of the University of Zurich
##
## Date:  October 31, 2017
##--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##
##
import sys, os, rdkit, rdkit.Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdmolops
from rdkit.Chem import EditableMol
 
#Class1 : Phosporus_ylids / Class2 : Phosporus_ylids_poly / Class3 : Boronic_acid_ester / Class4 : Boronic_acid_ester_poly / Class5 : I_II_amines / Class6 : I_II_amines_poly / Class7 : alcohols / Class8 : alcohols_poly / Class9 : C+III_carbonyls / Class10 : C+III_carbonyls_poly / Class11 : sulfonyl_chlorides / Class12 : sulfonyl_chlorides_poly / Class13 : isocyanates_thio / Class14 : isocyanates_thio_poly / Class15 : alkyl_halides / Class16 : aryl_vinyl_halides / Class17 : Grignard / Class18 : ketons_aldehydes / Class19 : ketons_aldehydes_poly / Class20 : epoxydes / Class21 : terminal_alkenes / Class22 : terminal_alkenes_poly

try:
  Bromids_file = open(sys.argv[1],'r')
  Nu_file = open(sys.argv[2],'r')
  Title = sys.argv[3]
except:
  print 'Usage \'python script.py bromid_sdf nucleophile_sdf output_name\''
  exit()
NuBr_Mol_list = []
Bromids_string = Bromids_file.readlines()
Nu_string = Nu_file.readlines()
Bromids_file.close()
Nu_file.close()
count=0
linecount=0
sectioncount=0
outputstring = ''
Nu_mol_list = {}
Br_mol_list = {}
competitivecount=0
dictcompetitive = {'halide':[' 1 ',' 2 ',' 3 ',' 4 ',' 9 ',' 10 ',' 11 ',' 12 ',' 13 ',' 14 ',' 20 '],'Boronic':[' 15 ',' 16 ',' 1 ',' 2 ']}

"""
1 = 'Phosporus_ylids';class
2 = 'Phosporus_ylids_poly';class
3 = 'Boronic_acid_ester';class
4 = 'Boronic_acid_ester_poly';class
5 = 'I_II_amines';class
6 = 'I_II_amines_poly';class
7 = 'alcohols';class
8 = 'alcohols_poly';class
9 = 'C+III_carbonyls';class
10 = 'C+III_carbonyls_poly';class
11 = 'sulfonyl_chlorides';class
12 = 'sulfonyl_chlorides_poly';class
13 = 'isocyanates_thio';class
14 = 'isocyanates_thio_poly';class
15 = 'alkyl_halides';class
16 = 'aryl_vinyl_halides';class
17 = 'Grignard';class
18 = 'ketons_aldehydes';class
19 = 'ketons_aldehydes_poly';class
20 = 'epoxydes';class
21 = 'terminal_alkenes';class
22 = 'terminal_alkenes_poly';
"""
#######################################DICTIONARY_BROMID_REACTANTS############################################################
output=[]
atomoutput=[]
bondoutput=[]
atomidx=unknowCAS=0
atomsection=False
bondsection=False
CAS_nb=None
reactivecenter=''
competitive_Nu =False
outputlist=[]
nitril_flag=CIII_flag=False
br_count=0
br_allcount=0
flag_Mg=False
competitive_nucleophile=[]
flag_Mg = False
atomcount=0
activehalogen=activecarbon=atomidx=0
CIII_flag=False
nitril_flag=False
keton_aldh_flag = False

for n in range(0,len(Bromids_string)):
        try:
            outputlist.append(Bromids_string[n].split('\n')[0])
        except:
            print 'Line refuses command \"outputlist.append(Bromids_string[n].split(\'\\n\')[0])\"' + Bromids_string[n]
        if 'REACTIVE_CENTERS' in Bromids_string[n]:
            try:
                reactivecenter=Bromids_string[n].split('<REACTIVE_CENTERS>')[1].split('[')[1].split(']')[0].split(',')
            except:
                print 'Line refuses command Bromids_string[n].split(<REACTIVE_CENTERS>)' + Bromids_string[n]
        elif 'CLASS' in Bromids_string[n]:
            classes=Bromids_string[n].split('<CLASS>')[1]
            if ('3' in classes and '15' in classes) or ('4' in classes and '15' in classes) or ('3' in classes and '16' in classes) or ('4' in classes and '16' in classes) or ('15' in classes and '16' in classes):
                  competitive_Nu = True
        #    for function in ['halide','sulfonyl_chlorides','C+III','keton']:
        #        if function in sys.argv[1]:
        #            competitive_nucleophile = dictcompetitive[function]
        #    for classe in competitive_nucleophile:
        #        if classe in classes:
        #             competitive_Nu = True
        elif 'CAS' in Bromids_string[n]:
            CAS_nb = Bromids_string[n].split('<CAS>')[1]
        elif '$$$$' in Bromids_string[n]:
            activecarbon=int(reactivecenter[0])+1
            activehalogen=int(reactivecenter[1])+1
            newatomidx = {}
            atomcount=bond=atomidx=None
            atomidx=0
            for m in range(0,len(outputlist)):
                if '99 V2000' in outputlist[m]:
                    atomcount = int(outputlist[m].split()[0])
                    bond= int(outputlist[m].split()[1])-1
                    atomsection=True
                    continue
                elif atomsection==True:
                    if 'Mg' in outputlist[m] and 'halide' in sys.argv[1]:
                        flag_Mg=True
                    atomidx+=1
                    if len(outputlist[m].split()) != len(outputlist[m+1].split()):
                        bondsection=True
                        atomsection=False
                    atomoutput.append(outputlist[m])
                elif bondsection==True:
                    if len(outputlist[m].split()) != len(outputlist[m+1].split()) or 'M' in outputlist[m+1] or 'CHG' in outputlist[m+1]:
                         bondsection=False
                    if (activehalogen == int(outputlist[m].split()[0]) or activehalogen == int(outputlist[m].split()[1])) and (activecarbon == int(outputlist[m].split()[0]) or activecarbon == int(outputlist[m].split()[1])):
                        if nitril_flag == True:
                             bondoutput.append(u"{0:>3s}{1:>3s}{2:>3s}  0\n".format(outputlist[m].split()[0],outputlist[m].split()[1],str(2)))
                             bond=bond+1
                        pass
                        if keton_aldh_flag == True:
                             bondoutput.append(u"{0:>3s}{1:>3s}{2:>3s}  0\n".format(outputlist[m].split()[0],outputlist[m].split()[1],str(1)))
                             bond=bond+1
                        pass
                    else:
                        bondoutput.append(u"{0:>3s}{1:>3s}{2:>3s}  0\n".format(outputlist[m].split()[0],outputlist[m].split()[1],outputlist[m].split()[2]))
            if CAS_nb == None or len(list(CAS_nb)) < 4:
                unknowCAS +=1
                CAS_nb = 'no_CAS_'+str(unknowCAS)
            if competitive_Nu == False and flag_Mg==False:# and nitril_flag==False:
                Br_mol_list[CAS_nb] = [atomoutput,bondoutput,atomcount,bond,activecarbon]
                br_count+=1
                if br_count == 100:
                  br_allcount +=br_count 
                  br_count=0
            else:
                competitivecount+=1
            reactivecenter=None
            atomoutput=[]
            bondoutput=[]
            competitive_nucleophile=[]
            CAS_nb=None
            competitive_Nu = False
            flag_Mg = False
            atomsection=bonsection=False
            nitril_flag=False
            outputlist=[]
            atomcount=0
            activehalogen=activecarbon=atomidx=0
            CIII_flag=False
            nitril_flag=False
            keton_aldh_flag = False

br_allcount +=br_count
print '-------------REACTANT 1--------------'
print 'Non-tolerated Building Blocks     :  '+str(competitivecount)
print 'Building Blocks kept for coupling :  '+str(br_allcount)

#######################################DICTIONARY_NUCLEOPHILE_REACTANTS############################################################
mol_count=0
outputstring = []
CAS_nb=None
classes=CAS_nb=smarts=reactive_center=''
competitive_Nu=flag_halogen=False
nu_allcount=0
nu_count=0
flagmultiplectr=competitive_Nu=flag_halogen=False
mol_count=0
headoutput=''
outputstring = []
atomoutput=[]
bondoutput=[]
flag_halogen=False
CAS_nb=None
classes=CAS_nb=smarts=reactive_center=''
newatomidx={}
goodcount=0
competitivecount=0
if 'Boronic' in sys.argv[2] or 'Grignard' in sys.argv[2]:
    for n in range(0,len(Nu_string)):
        try:
            outputlist.append(Nu_string[n].split('\n')[0])
        except:
            print 'Line : \"%\" refuses command \"outputlist.append(Nu_string[n].split(\'\\n\')[0])\"' %Nu_string[n]
        if 'REACTIVE_CENTERS' in Nu_string[n]:
            try:
                reactive_center=Nu_string[n].split('<REACTIVE_CENTERS>')[1].split('[')[1].split(']')[0].split(',')
            except:
                print Nu_string[n]
        elif 'FRAG>' in Nu_string[n]:
            smarts = Nu_string[n].split('<FRAG>')[1]
        elif 'CLASS' in Nu_string[n]:
            classes=Nu_string[n].split('<CLASS>')[1]
            competitive_nucleophile = []
            for function in ['Boronic','Grignard']:
                if function in sys.argv[1]:
                    competitive_nucleophile = dictcompetitive[function]
            for classe in competitive_nucleophile:
                if classe in classes:
                     competitive_Nu = True
            if ('3' in classes and '15' in classes) or ('4' in classes and '15' in classes) or ('3' in classes and '16' in classes) or ('4' in classes and '16' in classes) or ('15' in classes and '16' in classes):
                  competitive_Nu = True
        elif 'CAS' in Nu_string[n]:
            CAS_nb = Nu_string[n].split('<CAS>')[1]
        elif '$$$$' in Nu_string[n]:
            activecarbon=int(reactive_center[0])+1
            activehalogen=int(reactive_center[1])+1
            newatomidx = {}
            atomcount=bond=atomidx=None
            atomidx=0
            for m in range(0,len(outputlist)):
                if '99 V2000' in outputlist[m]:
                    atomcount = int(outputlist[m].split()[0])
                    bond= int(outputlist[m].split()[1])-1
                    headoutput=' '+str(atomcount)+' '+str(bond)+'  0  0  0  0  0  0  0  0999 V2000'
                    atomsection=True
                    continue
                elif atomsection==True:
                    atomidx+=1
                    if int(len(outputlist[m].split())) != int(len(outputlist[m+1].split())):
                        bondsection=True
                        atomsection=False
                    atomoutput.append(outputlist[m])
                elif bondsection==True:
                    if len(outputlist[m].split()) != len(outputlist[m+1].split()) or 'M' in outputlist[m+1]:
                         bondsection=False
                    if (activehalogen == int(outputlist[m].split()[0]) or activehalogen == int(outputlist[m].split()[1])) and (activecarbon == int(outputlist[m].split()[0]) or activecarbon == int(outputlist[m].split()[1])):
                        pass
                    else:
                        bondoutput.append(u"{0:>3s}{1:>3s}{2:>3s}  0".format(outputlist[m].split()[0],outputlist[m].split()[1],outputlist[m].split()[2]))
            if CAS_nb == None or len(list(CAS_nb)) < 4:
                unknowCAS +=1
                CAS_nb = 'no_CAS_'+str(unknowCAS)
            if competitive_Nu == False:
                goodcount+=1
                Nu_mol_list[CAS_nb] = [headoutput,atomoutput,bondoutput,activecarbon,smarts]
            else:
                competitivecount+=1
            reactive_center=None
            CAS_nb=None
            classes=CAS_nb=smarts=reactive_center=''
            flagmultiplectr=competitive_Nu=False
            competitive_Nu = False
            atomsection=bonsection=False
            headoutput=''
            outputlist=[]
            atomoutput=[]
            bondoutput=[]
            newatomidx = {}
            atomcount=0
            activehalogen=activecarbon=atomidx=0
print '-------------REACTANT 2--------------'
print 'Non-tolerated Building Blocks     :  '+str(competitivecount)
print 'Building Blocks kept for coupling :  '+str(goodcount)
########################################REACTION_SECTION#################################################################
print '\n\nNo more than 10\'000 ligands will be generated - if you wish to change the number of generated ligands go to line 245 in the Script \"if threshold\"\n\n'
threshold=0
for CAS in (Nu_mol_list.keys()):
  for Br_CAS in (Br_mol_list.keys()):
    threshold+=1
    if threshold==10000:
        exit()
    output=''
    atomsection=False
    bondsection=False
    endbondsection=False
    entry=False
    output+=CAS+' '+Br_CAS+'\nNucleophile Coupling\n\n'
    atomsNu=int(Nu_mol_list[CAS][0].split()[0])
    atomcount = int(Nu_mol_list[CAS][0].split()[0])+Br_mol_list[Br_CAS][2]
    bond= int(Nu_mol_list[CAS][0].split()[1])+Br_mol_list[Br_CAS][3]+1
    newline=' '+str(atomcount)+' '+str(bond)+'  0  0  0  0  0  0  0  0999 V2000\n'
    output+=newline
    for n in range(0,len(Nu_mol_list[CAS][1])):
            output+=Nu_mol_list[CAS][1][n]+'\n' 
    for t in range(0,len(Br_mol_list[Br_CAS][0])):
            output+=Br_mol_list[Br_CAS][0][t]+'\n'
    for f in range(0,len(Nu_mol_list[CAS][2])):
            output+=Nu_mol_list[CAS][2][f]+'\n'
    bondstoadd=Br_mol_list[Br_CAS][1]
    for bondline in bondstoadd:
                     atom1=str(int(bondline.split()[0])+atomsNu)
                     atom2=str(int(bondline.split()[1])+atomsNu)
                     order=bondline.split()[2]
                     try:
                         output+= u"{0:>3s}{1:>3s}{2:>3s}  0\n".format(atom1,atom2,order) 
                     except:
                         output+=' '+atom1+' '+atom2+' '+order+'\n'
                #try:
    output+= u"{0:>3s}{1:>3s}{2:>3s}  0\n".format(str(Nu_mol_list[CAS][3]),str(int(Br_mol_list[Br_CAS][4])+atomsNu),'1')
                #except:
                #    output+=' '+str(int(Nu_mol_list[CAS][1].split()[0])+1)+' '+str(int(Br_mol_list[Br_CAS][4])+atomsNu+1)+' '+'1'
    atomsection=bondsection=entry=False
        #elif entry==False:
            #output+=filelines[n]
    output+='M  END\n\n$$$$'
    os.system(" echo \'"+output+"\' >> "+Title+".sdf")
    
