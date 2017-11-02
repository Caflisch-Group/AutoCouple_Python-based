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
import sys
import sys, os, rdkit, rdkit.Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdmolops
from rdkit import Chem

try:
    sys.argv[1]
except:
    print "    Usage \'python script.py building-blocks_library.sdf\'"
    exit()


# SEE REACTANTS DESCRIPTION AT THE END OF PAGE
filein = open(sys.argv[1],'r')
filelines = filein.readlines()
filein.close()
class_dic = {}
class_name={}
class_name={}
class_name[1] = 'Phosporus_ylids';class_name[2] = 'Phosporus_ylids_poly';class_name[3] = 'Boronic_acid_ester';class_name[4] = 'Boronic_acid_ester_poly';class_name[5] = 'I_II_amines';class_name[6] = 'I_II_amines_poly';class_name[7] = 'alcohols';class_name[8] = 'alcohols_poly';class_name[9] = 'C+III_carbonyls';class_name[10] = 'C+III_carbonyls_poly';class_name[11] = 'sulfonyl_chlorides';class_name[12] = 'sulfonyl_chlorides_poly';class_name[13] = 'isocyanates_thio';class_name[14] = 'isocyanates_thio_poly';class_name[15] = 'alkyl_halides';class_name[16] = 'aryl_vinyl_halides';class_name[17] = 'Grignard';class_name[18] = 'ketons_aldehydes';class_name[19] = 'ketons_aldehydes_poly';class_name[20] = 'epoxydes';class_name[21] = 'terminal_alkenes';class_name[22] = 'terminal_alkenes_poly';
count=0
linecount=0
sectioncount=0
dumped=0
CAS_list=[]
reactive_wittig=[]
reactive_boron=[]
reactive_CIII=[]
reactive_SO2Cl=[]
reactive_isocyanate=[]
reactive_epox=[]
redundant=0
mol_list=[]
unknownCAS=0
mol_count=0
outputstring = ''
CAS_nb = None
for i in range(0,len(filelines)):
    outputstring += filelines[i]
    if '<CAS>' in filelines[i]:
        try:
            CAS_nb = filelines[i].split('<CAS>')[1].split('\n')[0] 
        except:
            try:
                 CAS_nb = filelines[i].split()[1]
            except:
                 print filelines[i]
    if 'M  END' in filelines[i]:
        limitedoutput = outputstring
    if "$$$$" in filelines[i]:
        if CAS_nb == None:
            unknownCAS+=1
            CAS_nb ='no_CAS_X'+str(unknownCAS)
            print 'warning'
        sdmol = rdkit.Chem.MolFromMolBlock(outputstring,sanitize=False)
        sdmolsaninititzed = rdkit.Chem.MolFromMolBlock(outputstring,sanitize=True)
        if sdmol != None:
            mol_count+=1
            mol_list.append([sdmol,CAS_nb,outputstring,limitedoutput,sdmolsaninititzed])
        outputstring = ''
        limitedoutput = None
        sdmol = None
        CAS_nb = None
        sdmolsaninititzed = None
print mol_count 

for item in mol_list:
  sdmol = item[0]
  CAS_nb = item[1]
  line = item[1]+'\n'
  sdf =  item[3]
  sdmolsan = item[4]
  linecount+=1
  sectioncount+=1
  if sectioncount==1000:
      print linecount
      print dumped
      print redundant
      sectioncount=0
  if CAS_nb not in CAS_list: 
    CAS_list.append(CAS_nb)
  else:
    redundant+=1
    continue
# Try to split molecule input into fragments (based on salts)
  if  sdmol == None:
      print 'nothing in %s' %CAS_nb
  if sdmol != None:
    fragcount=0
    fragatoms=0
    dictfrag = {}
    try:
        sdmolfrags = rdkit.Chem.rdmolops.GetMolFrags(sdmol, asMols=True)
    except:
        sdmolfrags = [ sdmol, ]
    for sdmol in sdmolfrags:
     fragcount+=1
     fragatoms=len(sdmol.GetAtoms())
     dictfrag[fragcount] = fragatoms
     count+=1
# Set all variables and lists to zero
     class1=class2=class3=class4=class5=class6=class7=class8=class9=class10=class11=class12=class13=class14=class15=class16=class17=class18=class19=class20=class21=class22=False
     halogenated_carbonyl = False
     list_nitrogen=list_oxygen=O_bond=N_bond=None
     all_carbo_halides=[]
     all_halides=[]
     all_nitrogens=[]
     conjug_nitrogens=[]
     all_oxygens=[]
     conjug_oxygens=[]
     all_CII_carbonyl=[]
     all_alkenes = []
     reactive_wittig=[]
     reactive_boron=[]
     reactive_CIII=[]
     reactive_SO2Cl=[]
     reactive_isocyanate=[]
     reactive_epox=[]
     reactivedictio = {}
     for n in range(1,23):
        class_dic[n] = False
# Exclude molecule/fragment with heavy atoms more than 26 or less than 3
     if sdmol.GetNumHeavyAtoms() <= 3 or sdmol.GetNumHeavyAtoms() > 35:
         dumped+=1
         continue
# Exclude molecule/fragment with heavy metals
     HeavyMetal = False
     for atom in sdmol.GetAtoms():
         if atom.GetAtomicNum() in [13,20,31,32,33,49,50,51,81,82,83,84,21,22,23,24,25,26,27,28,29,30,39,40,41,42,43,44,45,46,47,48,72,73,74,75,76,77,78,79,80,104,105,106,107,108,109,110,111,112,57,58,59,60,61,62,63,64,65,66,67,68,69,70,7189,90,91,92,93,94,95,96,97,98,99,100,101,102,103]:
             HeavyMetal = True
     if HeavyMetal == True:
         dumped+=1
         continue 
# Try to Sanitize otherwise run the search without 

     # Wittig precursors 
     # [#15]([#6])([#6])([#6])[#6], [#15X4](~[#8])([#8])([#8])[#6]
     for function in ['P(C)(c1ccccc1)(c1ccccc1)(c1ccccc1)','P(C)(~O)(O)(O)','P(C)(C=C(-C))(C=C(-C))(C=C(-C))']:
        for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function, sanitize=False)):
            if substruct[1] not in [x[0] for x in reactive_wittig]:
                reactive_wittig.append([substruct[1],substruct[0]])
     if len(reactive_wittig)==1:
            class1 = True
            class_dic[1] = True
            reactivedictio[1]=reactive_wittig 
     elif len(reactive_wittig) > 1:
            class2 = True
            class_dic[2] = True
            reactivedictio[2]=reactive_wittig
     # boronic acid, ester                                   |
     # [#5]([#8])([#8])[#6]
     for function in ['B(C)(O)(O)','B(c1ccccc1)(O)(O)','B(c1****1)(O)(O)']:
        for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function, sanitize=False)):
            if substruct[1] not in [x[0] for x in reactive_boron]:
                reactive_boron.append([substruct[1],substruct[0]])
     if len(reactive_boron)==1:
            class3 = True
            class_dic[3] = True
            reactivedictio[3]=reactive_boron
     elif len(reactive_boron) > 1:
            class4 = True
            class_dic[4] = True
            reactivedictio[4]=reactive_boron
     # acids, acyl chlorides, anhydrides, nitriles           | 
     # [#8X1][#6X3]([#6])[#8!H0]','[#8X1][#6X3]([#6])[#17]','[#8X1][#6X3]([#6])[#8][#6X3]([#8X1])[#6]','[#7]#[#6X2]']
     for function in ['C(=O)(OC)C','C(=O)(Cl)','C(=O)(Br)','C(=O)(OC(=O))','C#N']:
        for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function, sanitize=False)):
            if substruct[0] not in [x[0] for x in reactive_CIII]:
                try:
                     reactive_CIII.append([substruct[0],substruct[2],substruct[1]])
                except:
                     reactive_CIII.append([substruct[0],substruct[1]])
     if len(reactive_CIII)==1:
            class9 = True
            class_dic[9] = True
            reactivedictio[9]=reactive_CIII
     elif len(reactive_CIII) > 1:
            class10 = True
            class_dic[10] = True
            reactivedictio[10]=reactive_CIII

     # sulfonyl chlorides                                    |
     # [#16X4](~[#8])(~[#8])(~[#6])[#17]
     for function in ['S(~O)(~O)(Cl)','S(=O)(=O)(Cl)']:
        for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function, sanitize=False)):
            if substruct[0] not in [x[0] for x in reactive_SO2Cl]:
                reactive_SO2Cl.append([substruct[0],substruct[3]])
     if len(reactive_SO2Cl)==1:
            class11 = True
            class_dic[11] = True
            reactivedictio[11]=reactive_SO2Cl
     elif len(reactive_SO2Cl) > 1:
            class12 = True
            class_dic[12] = True
            reactivedictio[12]=reactive_SO2Cl 
     # isocyanates, isothiocyanates                          |
     #[#7X2]~[#6X2]~[#8X1,#16X1]
     for function in ['N=C=O','N=C=S']:
        for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function, sanitize=False)):
            if substruct[1] not in [x[0] for x in reactive_isocyanate]:
                reactive_isocyanate.append([substruct[1],substruct[0]])
     for function in ['C(=N)(=O)','C(=N)(=S)']:
        for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function, sanitize=False)):
            if substruct[0] not in [x[0] for x in reactive_isocyanate]:
                reactive_isocyanate.append([substruct[0],substruct[0]])
     if len(reactive_isocyanate)==1:
            class13 = True
            class_dic[13] = True
            reactivedictio[13]=reactive_isocyanate
     elif len(reactive_isocyanate) > 1:
            class14 = True
            class_dic[14] = True
            reactivedictio[14]=reactive_isocyanate
     # alkyl-halides except Fluorine	                   | no boronic acid, ester
     # [#6X4]~[#9,#17,#35,#53] 
     for function in ['CBr','CCl','CI','c1(Br)*****1','c1(Cl)*****1','c1(I)*****1','C(Br)=C','C(Cl)=C','C(I)=C','C(Br):C','C(Cl):C','C(I):C']:
         for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function, sanitize=False)):
              if substruct[1] not in all_halides:
                  all_halides.append(substruct[1])
              pass
              if substruct[0] not in [x[0] for x in all_carbo_halides]:
                      all_carbo_halides.append([substruct[0],substruct[1]])
     for function in ['[#6]~[#17,#35,#53]','[#6]:1(~[#17,#35,#53]):[#6]:[#6]:[#6]:[#6]:[#6]1']:
         for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmarts(function)):
              if substruct[1] not in all_halides:
                  all_halides.append(substruct[1])
              pass
              if substruct[0] not in [x[0] for x in all_carbo_halides]:
                      all_carbo_halides.append([substruct[0],substruct[1]])
     if len(all_carbo_halides) == 1:
        for atom in sdmol.GetAtomWithIdx(all_carbo_halides[0][0]).GetNeighbors():
           if 8==atom.GetAtomicNum() or 16==atom.GetAtomicNum():
               if sdmol.GetAtomWithIdx(all_carbo_halides[0][0]).IsInRing() ==False:
                   halogenated_carbonyl = True
        if halogenated_carbonyl == False:
           class15 = True
           class_dic[15] = True
           for bond in sdmol.GetAtomWithIdx(all_carbo_halides[0][0]).GetBonds():
             bontype = str(bond.GetBondType())
             for typ in ['AROM','TRIPL','DOUBL']:
               if typ in bontype:
                 class15 = False
                 class_dic[15] = False
     # aryl-halides, vinyl-halides                           | no boronic acid, ester
     # [*]1~[*]~[*]~[*]~[#6X3]~1~[#9,#17,#35,#53], [*]1~[*]~[*]~[*]~[*]~[#6X3]~1~[#9,#17,#35,#53], [#6X3]~[#6X3]~[#9,#17,#35,#53] 
                 classe16 = True
                 class_dic[16] = True

     # I or II aryl-, vinyl-amines                           | no amid, thioamid, sulfonamid
     # [#6X4]~[N!H0] but not [#7,#8,#16]=[#6,#16]~[N!H0]
     list_nitrogen = sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmarts('[#7]'))
     list_nitrogen += sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmiles('N'))
     for atom in list_nitrogen:
       if atom[0] not in all_nitrogens:
           all_nitrogens.append(atom[0])
       else:
           continue
     for atom in all_nitrogens:
       flag_multiplebond=False
       for N_bond in sdmol.GetAtomWithIdx(atom).GetBonds():
            if N_bond.GetBondTypeAsDouble() > 1.0:
                flag_multiplebond= True
       if flag_multiplebond==True:
           if atom not in conjug_nitrogens:
               conjug_nitrogens.append(atom)
     for function in ['N(N)(C)(C)','C(N(C)(C))','ON','BN','C(N)=S','C(N)=O','P(N)=O','S(N)(=O)(=O)','C(N)(=N)','C-N-N','N-N-C','C(N)(:N)']:
        for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function, sanitize=False)):
            if substruct[1] not in conjug_nitrogens:        
                 conjug_nitrogens.append(substruct[1])
     for function in ['[#6][#7]([#6])([#6])','[#7][#7]([#6])([#6])','[#8][#7]([#6])([#6])']:
        for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmarts(function)):
            if substruct[1] not in conjug_nitrogens:
                 conjug_nitrogens.append(substruct[1])
     for atom in all_nitrogens:
        condition1=False
        condition2=False
        if sdmol.GetAtomWithIdx(atom).IsInRingSize(5) and sdmol.GetAtomWithIdx(atom).GetIsAromatic():
            condition1=True
            for function in ['C(=O)N','S(=O)N','N:C:O','N:O','N:C:S','NC=N','N-N','NC:N']:
              for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function, sanitize=False)):
                 if atom in substruct:
                     condition2=True
            for N_bond in sdmol.GetAtomWithIdx(atom).GetBonds():
                if N_bond.GetBondTypeAsDouble() > 1.5:
                     condition2=True
        if condition1 and (condition2==False):
            try:
                conjug_nitrogens.remove(atom)
            except:
                print 'failed removing atoms in list conjug nitrogen XXXXXXXXXXXXXXXXXXXXXXXXXXX' 
     if len(all_nitrogens) - len(conjug_nitrogens) == 1:
          class5 = True
          class_dic[5] = True
     if len(all_nitrogens) - len(conjug_nitrogens) > 1:
          class6 = True
          class_dic[6] = True

     # alcools                                               | no ester, thioester, sulfonic acid
     # [#6]~[O!H0] but not [#7,#8,#16]=[#6,#16]~[O!H0]
     list_oxygen = sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmarts('[#6][#8]'))
     list_oxygen += sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmiles('CO', sanitize=False))
     for substruct in list_oxygen:
          if substruct[1] not in all_oxygens:
              all_oxygens.append(substruct[1])
          else:
              continue
     for atom in all_oxygens:
       flag_multiplebond=False
       for O_bond in sdmol.GetAtomWithIdx(atom).GetBonds():
           if O_bond.GetBondTypeAsDouble() > 1.0:
               flag_multiplebond=True
       if flag_multiplebond==True:
           if atom not in conjug_oxygens:
               conjug_oxygens.append(atom)
           else:
               continue
     for function in ['NO','COC','BO','C(O)=S','C(O)=O','P(O)=O','S(O)(=O)(=O)']:
        for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function, sanitize=False)):
            if substruct[1] not in conjug_oxygens:
                conjug_oxygens.append(substruct[1])
            else:
                continue
     for function in ['[#14][#8]','[*]1-[#8]-[*]=[*]-[*]=1','[#6][#8][#6]','[#6](-[#8]-[#1])(=[#8])','[#6](-[#8]-[#6])(=[#8])']:
        for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmarts(function)):
            if substruct[1] not in conjug_oxygens:
                conjug_oxygens.append(substruct[1])
            else:
                continue
     if len(all_oxygens) - len(conjug_oxygens) ==1:
        class7 = True
        class_dic[7] = True
     elif len(all_oxygens) - len(conjug_oxygens) > 1:
           class8 = True
           class_dic[8] = True

     # ketones, aldehydes            for reduc elimination   | no halide, ester, acid, anhydride, acyl chloride, isocyanate, amid 
     # [#6][#6X3](~[#8X1])[#1,#6]
     for function in ['CC(=O)(C)']:#,'CC(=O)(c)','cC(=O)(C)','cC(=O)(c)']:
         for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function, sanitize=False)):
              if substruct[1] not in [x[0] for x in all_CII_carbonyl]:
                  all_CII_carbonyl.append([substruct[1],substruct[2]])
     for function in ['[#6]-[#6](=[#8])',]:
         for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmarts(function)):
              if len(sdmol.GetAtomWithIdx(substruct[1]).GetNeighbors()) == 2:
                  if substruct[1] not in [x[0] for x in all_CII_carbonyl]:
                      all_CII_carbonyl.append([substruct[1],substruct[2]])
     if len(all_CII_carbonyl) == 1:
        class18 = True
        class_dic[18] = True
        reactivedictio[18]=all_CII_carbonyl
     elif len(all_CII_carbonyl) > 1:
        class19 = True
        class_dic[19] = True
        reactivedictio[19]=all_CII_carbonyl

     # epoxides                      for Grignards AN        | no halide, ester, acid, anhydride, acyl chloride, isocyanate, amid 
     # [#1,#6:1][#6X4:2]1[#8X2,#7X3:3][#6X4:4]1[#1,#6:5]
     for function in ['C1OC1']:
        for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function, sanitize=False)):
                  if len(sdmol.GetAtomWithIdx(substruct[0]).GetNeighbors()) < len(sdmol.GetAtomWithIdx(substruct[2]).GetNeighbors()):
                       if substruct[0] not in [x[0] for x in reactive_epox]:
                            reactive_epox.append([substruct[0],substruct[1]])
                  elif len(sdmol.GetAtomWithIdx(substruct[2]).GetNeighbors()) < len(sdmol.GetAtomWithIdx(substruct[0]).GetNeighbors()):
                       if substruct[2] not in [x[0] for x in reactive_epox]:
                            reactive_epox.append([substruct[2],substruct[1]])
     for function in ['[#6]-1-[#8]-[#6]1']:
        for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmarts(function)):
                  if len(sdmol.GetAtomWithIdx(substruct[0]).GetNeighbors()) < len(sdmol.GetAtomWithIdx(substruct[2]).GetNeighbors()):
                       if substruct[0] not in [x[0] for x in reactive_epox]:
                            reactive_epox.append([substruct[0],substruct[1]])
                  elif len(sdmol.GetAtomWithIdx(substruct[2]).GetNeighbors()) < len(sdmol.GetAtomWithIdx(substruct[0]).GetNeighbors()):
                       if substruct[2] not in [x[0] for x in reactive_epox]:
                            reactive_epox.append([substruct[2],substruct[1]])

     if len(reactive_epox)==1:
            class20 = True
            class_dic[20] = True
            reactivedictio[20]=reactive_epox
     elif len(reactive_epox) > 1:
            class21 = True
            class_dic[20] = True
            reactivedictio[20]=reactive_epox
     # alkyl-, aryl-, vinyl-halides  for Grignards           | no protic H, ester, multihalides, anhydride, isocyanates, amid
     # [#6]~[#9,#17,#35,#53] but not [N!H0,O!H0,S!H0], [#7,#8,#16]=[#6,#16]~[#7,#8]with-bond1-2-inring, [#7X2]~[#6X2]~[#8X1]  
     if (class15 == True or class16 == True) and (( class1 or class2 or class3 or class4 or class5 or class6 or class7 or class8 or class9 or class10 or class11 or class12 or class13 or class14 or class20 or sdmol.HasSubstructMatch(rdkit.Chem.MolFromSmiles('C(=O)N')) or sdmol.HasSubstructMatch(rdkit.Chem.MolFromSmarts('[#6](-[#8]-[#1])(=[#8])')) or sdmol.HasSubstructMatch(rdkit.Chem.MolFromSmiles('[H]OC(=O)')) or sdmol.HasSubstructMatch(rdkit.Chem.MolFromSmarts('[#6](-[#7])(=[#8])')) or sdmol.HasSubstructMatch(rdkit.Chem.MolFromSmiles('[C](=O):[N]'))) == False):
        class17 = True
        class_dic[17] = True
        reactivedictio[17]=all_carbo_halides
     pass   
     if len(sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmarts('[#6]~[#12]'))) >= 1:
        class17 = True
        class_dic[17] = True
        reactivedictio[17]=[[sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmarts('[#6]~[#12]'))[0][0],sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmarts('[#6]~[#12]'))[0][1]],]
     # terminal alkenes              for Metathesis, Heck    |  
     # [#6][#6]=[#6H2] 
     for function in ['C-C=C','C-C:C','c1(-C=C)*****1','c1(-C=C)****1','c1(-C:C)*****1','c1(-C:C)****1']:
         for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function, sanitize=False)):
              if substruct[2] not in [x[0] for x in all_alkenes]:
                  if len(sdmol.GetAtomWithIdx(substruct[2]).GetNeighbors()) == 1:
                      all_alkenes.append([substruct[2],substruct[1],substruct[0]])
     for function in ['[#6][#6]=[#6]','[#6][#6]:[#6]']:
         for substruct in sdmol.GetSubstructMatches(rdkit.Chem.MolFromSmarts(function)):
             if substruct[2] not in [x[0] for x in all_alkenes]:
                  if len(sdmol.GetAtomWithIdx(substruct[2]).GetNeighbors()) == 1:
                      all_alkenes.append([substruct[2],substruct[1],substruct[0]])
     if len(all_alkenes)==1:
        class21 = True
        class_dic[21] = True
        reactivedictio[21]=all_alkenes
     elif len(all_alkenes) >1:
        class21 = True
        class_dic[21] = True
        reactivedictio[21]=all_alkenes

  if sdmolsan != None: 
     if class18 == False and class19 == False:
        all_CII_carbonyl=[]
        try:
           for function in ['CC(=O)(C)']:#,'CC(=O)(c)','cC(=O)(C)','cC(=O)(c)']:
               for substruct in sdmolsan.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function)):
                    if substruct[1] not in [x[0] for x in all_CII_carbonyl]:
                        all_CII_carbonyl.append([substruct[1],substruct[2]])
           for function in ['[#6]-[#6](=[#8])',]:
               for substruct in sdmolsan.GetSubstructMatches(rdkit.Chem.MolFromSmarts(function)):
                    if len(sdmolsan.GetAtomWithIdx(substruct[1]).GetNeighbors()) == 2:
                        if substruct[1] not in [x[0] for x in all_CII_carbonyl]:
                            all_CII_carbonyl.append([substruct[1],substruct[2]])
           if len(all_CII_carbonyl) == 1:
              class18 = True
              class_dic[18] = True
              reactivedictio[18]=all_CII_carbonyl
           elif len(all_CII_carbonyl) > 1:
              class19 = True
              class_dic[19] = True
              reactivedictio[19]=all_CII_carbonyl
        except:
           print CAS_nb
     if class7 == False and class8 == False:
        try:
         all_oxygens = []
         conjug_oxygens = []
         list_oxygen = sdmolsan.GetSubstructMatches(rdkit.Chem.MolFromSmarts('[#6][#8]'))
         list_oxygen += sdmolsan.GetSubstructMatches(rdkit.Chem.MolFromSmiles('CO'))
         for substruct in list_oxygen:
              if substruct[1] not in all_oxygens:
                  all_oxygens.append(substruct[1])
              else:
                  continue
         for atom in all_oxygens:
           flag_multiplebond=False
           for O_bond in sdmolsan.GetAtomWithIdx(atom).GetBonds():
               if O_bond.GetBondTypeAsDouble() > 1.0:
                   flag_multiplebond=True
           if flag_multiplebond==True:
               if atom not in conjug_oxygens:
                   conjug_oxygens.append(atom)
               else:
                   continue
         for function in ['NO','COC','BO','C(O)=S','C(O)=O','P(O)=O','S(O)(=O)(=O)']:
            for substruct in sdmolsan.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function)):
                if substruct[1] not in conjug_oxygens:
                    conjug_oxygens.append(substruct[1])
                else:
                    continue
         for function in ['[#14][#8]','[*]1-[#8]-[*]=[*]-[*]=1','[#6][#8][#6]','[#6](-[#8]-[#1])(=[#8])','[#6](-[#8]-[#6])(=[#8])']:
            for substruct in sdmolsan.GetSubstructMatches(rdkit.Chem.MolFromSmarts(function)):
                if substruct[1] not in conjug_oxygens:
                    conjug_oxygens.append(substruct[1])
                else:
                    continue
         if len(all_oxygens) - len(conjug_oxygens) ==1:
            class7 = True
            class_dic[7] = True
         elif len(all_oxygens) - len(conjug_oxygens) > 1:
               class8 = True
               class_dic[8] = True
        except:
          print CAS_nb
     if class5 ==False and class6 == False:
       try:
         all_nitrogens = []
         conjug_nitrogens = []
         list_nitrogen = sdmolsan.GetSubstructMatches(rdkit.Chem.MolFromSmarts('[#7]'))
         list_nitrogen += sdmolsan.GetSubstructMatches(rdkit.Chem.MolFromSmiles('N'))
         for atom in list_nitrogen:
           if atom[0] not in all_nitrogens:
               all_nitrogens.append(atom[0])
           else:
               continue
         for atom in all_nitrogens:
           flag_multiplebond=False
           for N_bond in sdmolsan.GetAtomWithIdx(atom).GetBonds():
                if N_bond.GetBondTypeAsDouble() > 1.0:
                    flag_multiplebond= True
           if flag_multiplebond==True:
               if atom not in conjug_nitrogens:
                   conjug_nitrogens.append(atom)
         for function in ['N(N)(C)(C)','C(N(C)(C))','ON','BN','C(N)=S','C(N)=O','P(N)=O','S(N)(=O)(=O)','C(N)(=N)','C-N-N','N-N-C','C(N)(:N)']:
            for substruct in sdmolsan.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function)):
                if substruct[1] not in conjug_nitrogens:
                     conjug_nitrogens.append(substruct[1])
         for function in ['[#6][#7]([#6])([#6])','[#7][#7]([#6])([#6])','[#8][#7]([#6])([#6])']:
            for substruct in sdmolsan.GetSubstructMatches(rdkit.Chem.MolFromSmarts(function)):
                if substruct[1] not in conjug_nitrogens:
                     conjug_nitrogens.append(substruct[1])
         for atom in all_nitrogens:
            condition1=False
            condition2=False
            if sdmolsan.GetAtomWithIdx(atom).IsInRingSize(5) and sdmolsan.GetAtomWithIdx(atom).GetIsAromatic():
                condition1=True
                for function in ['C(=O)N','S(=O)N','N:C:O','N:O','N:C:S','NC=N','N-N','NC:N']:
                  for substruct in sdmolsan.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function)):
                     if atom in substruct:
                         condition2=True
                for N_bond in sdmolsan.GetAtomWithIdx(atom).GetBonds():
                    if N_bond.GetBondTypeAsDouble() > 1.5:
                         condition2=True
            if condition1 and (condition2==False):
                try:
                    conjug_nitrogens.remove(atom)
                except:
                    print 'failed removing atoms in list conjug nitrogen XXXXXXXXXXXXXXXXXXXXXXXXXXX'
         if len(all_nitrogens) - len(conjug_nitrogens) == 1:
              class5 = True
              class_dic[5] = True
         if len(all_nitrogens) - len(conjug_nitrogens) > 1:
              class6 = True
              class_dic[6] = True
       except:
          print CAS_nb
     if class9 ==False and class10 == False:
        try:
         reactive_CIII=[]
         for function in ['C(=O)(OC)C','C(=O)(Cl)','C(=O)(Br)','C(=O)(OC(=O))','C#N']:
            for substruct in sdmolsan.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function)):
                if substruct[0] not in [x[0] for x in reactive_CIII]:
                    try:
                         reactive_CIII.append([substruct[0],substruct[2],substruct[1]])
                    except:
                         reactive_CIII.append([substruct[0],substruct[1]])
         if len(reactive_CIII)==1:
                class9 = True
                class_dic[9] = True
                reactivedictio[9]=reactive_CIII
         elif len(reactive_CIII) > 1:
                class10 = True
                class_dic[10] = True
                reactivedictio[10]=reactive_CIII
        except:
         print CAS_nb
     if class20 ==False:
       try:
         reactive_epox = []
         for function in ['C1OC1']:
            for substruct in sdmolsan.GetSubstructMatches(rdkit.Chem.MolFromSmiles(function)):
                      if len(sdmolsan.GetAtomWithIdx(substruct[0]).GetNeighbors()) < len(sdmolsan.GetAtomWithIdx(substruct[2]).GetNeighbors()):
                           if substruct[0] not in [x[0] for x in reactive_epox]:
                                reactive_epox.append([substruct[0],substruct[1]])
                      elif len(sdmolsan.GetAtomWithIdx(substruct[2]).GetNeighbors()) < len(sdmolsan.GetAtomWithIdx(substruct[0]).GetNeighbors()):
                           if substruct[2] not in [x[0] for x in reactive_epox]:
                                reactive_epox.append([substruct[2],substruct[1]])
         for function in ['[#6]-1-[#8]-[#6]1']:
            for substruct in sdmolsan.GetSubstructMatches(rdkit.Chem.MolFromSmarts(function)):
                      if len(sdmolsan.GetAtomWithIdx(substruct[0]).GetNeighbors()) < len(sdmolsan.GetAtomWithIdx(substruct[2]).GetNeighbors()):
                           if substruct[0] not in [x[0] for x in reactive_epox]:
                                reactive_epox.append([substruct[0],substruct[1]])
                      elif len(sdmolsan.GetAtomWithIdx(substruct[2]).GetNeighbors()) < len(sdmolsan.GetAtomWithIdx(substruct[0]).GetNeighbors()):
                           if substruct[2] not in [x[0] for x in reactive_epox]:
                                reactive_epox.append([substruct[2],substruct[1]])

         if len(reactive_epox)==1:
                class20 = True
                class_dic[20] = True
                reactivedictio[20]=reactive_epox
         elif len(reactive_epox) > 1:
                class21 = True
                class_dic[20] = True
                reactivedictio[20]=reactive_epox
       except:
         print CAS_nb    

  if sdmol != None or sdmolsan!=None:
     for n in [1,2,3,4,9,10,11,12,13,14,17,18,19,20,21,22]:
        if class_dic[n]:
            outputline = sdf
            outputline+=' \n><CAS>'+CAS_nb+'<CAS>'
            outputline+=' \n><REACTIVE_CENTERS>'
            for j in reactivedictio[n]:
                outputline+=' '+str(j)+' '
            outputline+='<REACTIVE_CENTERS> \n><CLASS>'
            for m in range(1,23):
                if class_dic[m]:
                    outputline+=' '+str(m)+' '
            outputline+='<CLASS>\n\n$$$$\n'
            smilesout = open("Rx_"+class_name[n]+".sdf",'a')
            smilesout.write(outputline)
            smilesout.close() 
     for n in [5,6]:    
        if class_dic[n]:
            outputline = sdf
            outputline+=' \n><CAS>'+CAS_nb+'<CAS>' 
            for i in conjug_nitrogens:
                try:
                    all_nitrogens.remove(i)
                except:
                    os.system("echo \'problem with "+CAS_nb+"\' >> Rx_misstreatment.dat")
            outputline+=' \n><REACTIVE_CENTERS>'
            for j in all_nitrogens:
                outputline+=' '+str(j)+' '
            outputline+='<REACTIVE_CENTERS> \n><CLASS>'
            for m in range(1,23):
                if class_dic[m]:
                    outputline+=' '+str(m)+' '
            outputline+='<CLASS>\n\n$$$$\n'
            smilesout = open("Rx_"+class_name[n]+".sdf",'a')
            smilesout.write(outputline)
            smilesout.close() 
     for n in [7,8]:
        if class_dic[n]:
            outputline = sdf
            outputline+=' \n><CAS>'+CAS_nb+'<CAS>'
            for i in conjug_oxygens:
                try:
                    all_oxygens.remove(i)
                except:
                    os.system("echo \'problem with "+CAS_nb+"\' >> Rx_misstreatment.dat")
            outputline+=' \n><REACTIVE_CENTERS>'
            for j in all_oxygens:
                outputline+=' '+str(j)+' '
            outputline+='<REACTIVE_CENTERS> \n><CLASS>'
            for m in range(1,23):
                if class_dic[m]:
                    outputline+=' '+str(m)+' '
            outputline+='<CLASS>\n\n$$$$\n'
            smilesout = open("Rx_"+class_name[n]+".sdf",'a')
            smilesout.write(outputline)
            smilesout.close()
     for n in [15,16]:
       if class_dic[n]:
            outputline = sdf
            outputline+=' \n><CAS>'+CAS_nb+'<CAS>'
            outputline+=' \n><REACTIVE_CENTERS>'
            outputline+=' '+str(all_carbo_halides[0])
            outputline+='<REACTIVE_CENTERS> \n><CLASS>'
            for m in range(1,23):
                if class_dic[m]:
                    outputline+=' '+str(m)+' '
            outputline+='<CLASS>\n\n$$$$\n'
            smilesout = open("Rx_"+class_name[n]+".sdf",'a')
            smilesout.write(outputline)
            smilesout.close()
     if (class1 or class2 or class3 or class4 or class5 or class6 or class7 or class8 or class9 or class10 or class11 or class12 or class13 or class14 or class15 or class16 or class17 or class18 or class19 or class20 or class21 or class22) ==False:
         smilesout=open('Rx_failed.sdf','a')
         smilesout.write(sdf+'\n$$$$\n')
         smilesout.close() 
         dumped+=1
print count
print dumped
exit()
"""
#REACTANT GROUPS
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
