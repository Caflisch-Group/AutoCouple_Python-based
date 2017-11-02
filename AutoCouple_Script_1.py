## Python Scripts for AutoCouple
##
## Author:  Laurent Batiste
##
## Affiliation:  A. Caflisch' group at the Department of Biochemistry of the University of Zurich
##
## Date:  October 31, 2017
##--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#!bin/python
import sys, os, rdkit, rdkit.Chem
from rdkit.Chem import rdmolfiles
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit.Chem import Lipinski
from rdkit.Chem import rdmolops
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from datetime import datetime
import re

try:
    sys.argv[1]
except:
    print "    Usage \'python script.py building-blocks_library.sdf\'"
    exit()
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Standart command line requires the providers' libraries as an SDF files:
# > python script.py library_provider1.sdf library_provider2.sdf library_provider3.sdf 

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#FIRST SECTION : EACH PROVIDER HAS ITS OWN LAYOUT FOR SDF-FILE, THE FOLLOWING SCRIPT REQUIRES SOME MANUAL INPUT BEFORE TO PROCEED

# dicofile is a dictionary storing the library itself
# dicoprovider is a dictionary storing the provider's name (given by user)
# the dictionaries : dicounit;dicoprice;dicoqty;;dicocaskey;dicocatalogkey;dicocurrency store, for each different layout, the indicators preceding information about respectively the unit (gramm,kg,liter), the price, the quantity, the CAS key, the catalogue number, the currenry 
dicounit={};dicoprice={};dicoqty={};dicofile={};dicoprovider={};dicocaskey={};dicocatalogkey={};dicocurrency={}
for i in range(1, len(sys.argv)):
    dicoqty[i]={}
    dicoprice[i]={}
    print "-"*200+"\nFIRST SECTION : EACH PROVIDER HAS ITS OWN LAYOUT FOR MOL-FILE, THE FOLLOWING SCRIPTS REQUIRE SOME MANUAL INPUT BEFORE TO PROCEED\n\n"+"-"*200+"\n\t\tChemicals provider no."+str(i)+" in file :"+sys.argv[i]
    provider= raw_input('\n\n\t\tPlease write below the name of the chemicals provider in lowercase without \"_\" :\n\t\t')
    filein = open(sys.argv[i],'r')
    dicofile[i] = filein.readlines()
    filein.close()
    output=''
    items= dicofile[i].count('$$$$\n')
    count=0
    for line in dicofile[i]:
        output+=line
        if '$$$$' in line:
          count+=1
          if count== int(items/9): 
            output=''
          if count== int(int(items/9)+1):
            count=0
            print output
            CAS_key=raw_input('\n\n\t\tPlease read the Molfile shown above and Copy/Paste below the line preceding the CAS number (should be <CAS>, <rn> or <ref>) :\n\t\t(enter \'none\' if the information cannot be found)\n\t\t(press ENTER in order to skip this section, NO information about the chemicals\' CAS and catalog numbers or the prices will be stored)\n\t\t')
            if 'none' in CAS_key:
                 output=''
                 continue
            Catalog_key=raw_input('\n\n\t\tThank you, now please read the Molfile show above and Copy/Paste below the line preceding the CATALOGUE number (looks like <ID>, <CAT> or any other):\n\t\t(enter \'none\' if the information cannot be found)\n\t\t')
            if 'none' in Catalog_key:
                 output=''
                 continue
            qtyinfo = raw_input('\n\n\t\tThank you, now please copy/paste below the line preceding the quantity or size no.1 (looks like <Q1> <quantity1> or <size_1>):\n\t\t(enter \'none\' if the information cannot be found)\n\t\t(press ENTER in order to skip this section, NO information about the chemicals\' prices will be stored)\n\t\t')
            dicoqty[i]=qtyinfo.replace('1','IDX')
            if 'none' in dicoqty[i]:
                        output=''
                        continue
            pricesinfo = raw_input('\n\n\t\tThank you, now please copy/paste below the line preceding the price no.1 that is related to this quantity (looks like <P1> or <price_1>):\n\t\t')
            dicoprice[i]=pricesinfo.replace('1','IDX')
            unitinfo = raw_input('\n\n\t\tThank you, now please copy/paste below the line preceding the unit that is related to this quantity (looks like <unit>):\n\t\t')
            dicounit[i]=unitinfo.replace('1','IDX')
            currencyinfo = raw_input('\n\n\t\tThank you, now please specify the currency used for prices in this catalog as its international symbol (EUR, CHF, USD, GBP):\n\t\t')
            dicocurrency[i]=currencyinfo;currencyinfo =''
            dicoprovider[i]=provider;provider=''
            dicocaskey[i]=CAS_key;CAS_key=''
            dicocatalogkey[i]=Catalog_key;Catalog_key=''
            break

# print out a summary of all keyword entered by the user:
for i in range(1, len(sys.argv)):
      print dicoprovider[i]
      print dicocaskey[i]
      print dicocatalogkey[i]
      print dicoqty[i]
      print dicoprice[i]
      print dicounit[i]
      print dicocurrency[i]
      
print "-"*200+"\nSECOND SECTION : THE LIBRARIES ARE SCREENED AND SORTED OUT ACCORDING TO CAS NUMBER - 450'000 mol files requires approx. 5 hours. \n\n"+"-"*200+"\n" 

dic_data = {} ; dico_provider = {} ; dico_commentary={} 
# the dictionary nddictioprices will store the quantity, and related price and unit (g,kg,l,ml) for each molecule
nddictioprices = {'QT':{},'PR':{},'UN':{}} ; rddictioprices = {}
unknowCAS = 0
unknowcatalog = 0
list_CAS_nb = []
advance = 0
safetycount = 0
redundantCAS=good=0
# Through filters specified later molecules with resin, metal, too many rotatable bond, multiple chiral centers , unwanted features or molecules that have already been stored will be excluded.
# dumpedresin is the count of how many molecules were discarded because they contained resin ---- same applies for dumpedmetal ; dumpedHAC ; dumpedfastlooping ; dumpedrdkit ; dumpedrotatable ; dumpedchiral ; dumpedtrimethylene
dumpedresin=dumpedmetal=dumpedHAC=dumpedfastlooping=dumpedrdkit=dumpedrotatable=dumpedchiral=dumpedtrimethylene=0

failist=[]
# Each provider library is loaded and parsed; each molecule of the library is screened one after the other.
for l in range(1, len(dicofile.keys())+1):
    for p in range(1,6):
                        nddictioprices['PR'][p]=''
                        nddictioprices['QT'][p]=''
                        nddictioprices['UN'][p]=''
    outputstring = ''
    CAS_nb = ''
    catalog_nb = ''
    commentary_flag=False
    commentary=''
    Amids_string=dicofile[l]
    chem_provider=dicoprovider[l]
    currency = dicocurrency[l]
    if len(list(currency)) == 0:
        currency = 'notspe'
    for i in range(0,len(Amids_string)):
        if i < len(Amids_string)-1 and ('POA' in Amids_string[i+1] or 'enq' in Amids_string[i+1].lower() or 'demand' in Amids_string[i+1].lower()):
            continue
        outputstring += Amids_string[i]
# The commercial information regarding the molecule is stored in "commentary"
        if commentary_flag == True and '$$$$' not in Amids_string[i]:
            commentary+= Amids_string[i]
# The crutial information about atoms and bonds is stored in "limitedoutput"
        if 'M  END' in Amids_string[i]:
            commentary_flag = True
            flagunit=False
            limitedoutput = outputstring
# Detect the line containing the catalog number:
        if len(list(dicocatalogkey[l])) !=0 and dicocatalogkey[l] in Amids_string[i]:
              try:
                  catalog_nb = Amids_string[i+1].split('[')[1].split(']')[0]
              except:
                  try:
                      catalog_nb = Amids_string[i+1].split('\n')[0]
                  except:
                      unknowcatalog+=1
                      catalog_nb = 'noCatNb-'+str(unknowcatalog)
# Detect the line containing the CAS number:
        if len(list(dicocaskey[l])) !=0 and dicocaskey[l] in Amids_string[i]:
            try:
                CAS_nb = Amids_string[i+1].split('[')[1].split(']')[0]
            except:
                try:
                    CAS_nb = Amids_string[i+1].split('\n')[0]
                except:
                    unknowCAS +=1
                    CAS_nb = 'noCAS-'+chem_provider+catalog_nb
# Each lines of the prices' section are parse and the price, the quantity and the related unit are stored independently 
        if len(list(dicoqty[l]))!=0:
          for h in range(1,6):
            if dicoqty[l].replace('IDX',str(h)) in Amids_string[i]:
                for unit_ref in ['g','mg','ug','kg','gr','l','ml','lt','ul']:
                    if unit_ref in Amids_string[i+1].lower():
                        nddictioprices['UN'][h] = unit_ref
                        flagunit=True
                qty = (Amids_string[i+1].lower()).replace(nddictioprices['UN'][h],'')
                if 'x' in Amids_string[i+1].lower():
                    try:
                         qty = float(qty.split('x')[0]) * float(qty.split('x')[1])
                    except:
                         failist.append('multiply'+ Amids_string[i+1].split('\n')[0])
                else:
                    try:
                         qty= float(qty.split('\n')[0]) 
                    except:
                         failist.append('floatcv'+ Amids_string[i+1].split('\n')[0])
                nddictioprices['QT'][h]=qty
        if len(list(dicoprice[l]))!=0:
          for h in range(1,6):
            if dicoprice[l].replace('IDX',str(h)) in Amids_string[i]:
                nddictioprices['PR'][h] = Amids_string[i+1].split('\n')[0]
                try:
                    nddictioprices['PR'][h] = float(Amids_string[i+1].split('\n')[0])
                except:
                    try:
                        nddictioprices['PR'][h] = float(re.findall(r"[-+]?\d*\.\d+|\d+",Amids_string[i+1])[0])
                    except:
                        failist.append('price' +Amids_string[i+1])
        if len(list(dicounit[l]))!=0:
          for h in range(1,6):
            if dicounit[l].replace('IDX',str(h)) in Amids_string[i]:
                nddictioprices['UN'][h] = Amids_string[i+1].split('\n')[0].lower()
# The line $$$$ indicate the end of a molecule and the start of a new one, the former molecule is analysed using chemical toolkit RDKIT:
        if "$$$$" in Amids_string[i]:
          commentary_flag = False
          sanitizeflag=False
          resineflag=False
#################################################################################################
# Prints a report every 5000 molecules
          safetycount+=1
          if safetycount==5000:
             advance += safetycount
             print 'Molfile items screened %s' %advance
             print 'Dumped resin:%s fastloop:%s rdkit:%s' %(dumpedresin,dumpedfastlooping,dumpedrdkit)
             print 'Dumped metal:%s HAC:%s rotablebonds:%s CH2-X-CH2:%s chiral:%s'  %(dumpedmetal,dumpedHAC,dumpedrotatable,dumpedtrimethylene,dumpedchiral)
             print 'Redundant via CAS %s, Stored %s' %(redundantCAS,good)
             safetycount = 0
             fileout = open('Global_Library_Reactants.sdf','w')
             for CAS in list_CAS_nb:
                 try:
                     fileout.write(dic_data[CAS]+dico_provider[CAS]+rddictioprices[CAS]+'\n\n$$$$\n')#+dico_commentary[CAS]+'\n\n$$$$\n')
                 except:
                     fileout.write('\n$$$$\n\n'+CAS+'\n\n$$$$\n')
             fileout.close()
###################################################################################################
# Assigns a Catalog number and CAS number if the information was not available
          if len(list(catalog_nb))==0:
                  unknowcatalog+=1
                  catalog_nb = 'noCatNb-'+str(unknowcatalog)
          if len(list(CAS_nb))==0 or len(CAS_nb.split('-')) < 2:
                  unknowCAS +=1
                  CAS_nb = 'noCAS-'+chem_provider+catalog_nb
          pass 
# If a molecule with the exact same CAS number was already stored, the present molecule information is appended to it
          if CAS_nb in list_CAS_nb:
                dico_provider[CAS_nb] =dico_provider[CAS_nb]+' ; '+chem_provider+':'+catalog_nb
                dico_commentary[CAS_nb]=dico_commentary[CAS_nb]+'\n\n> <SECTION FROM PROVIDER :>'+chem_provider+'\n\n'+commentary+'\n'
                rddictioprices[CAS_nb]=rddictioprices[CAS_nb]+'\n> <PRICE SECTION FROM CATALOG COMPOUND '+chem_provider+':'+catalog_nb+'>'
                for p in range(1,6):
                    if len(list(str(nddictioprices['PR'][p]))) > 0:
                        rddictioprices[CAS_nb]=rddictioprices[CAS_nb]+'\nPr.'+str(p)+'\t'+str(nddictioprices['PR'][p])+'\tCurr. '+currency+'\tQt.'+str(p)+'\t'+str(nddictioprices['QT'][p])+'\tunits. '+str(nddictioprices['UN'][p])
                dumpedfastlooping+=1
# Otherwise the molecule is analysed using chemical toolkit RDKIT 
          else:
           try:
                sdmol = None
                sdmol = rdkit.Chem.MolFromMolBlock(outputstring)
                sanitizeflag=True
##############################################################################################
           except:
            try:
             sdmol = None
             sdmol = rdkit.Chem.MolFromMolBlock(outputstring,sanitize=False)
             print 'NO SANITIZATION'
            except:
               dic_data[CAS_nb] = limitedoutput+'\n\n> <PROVIDER>\n '+chem_provider+':'+catalog_nb+'\n\n> <RDKIT READING FAILED>'+'\n\n> <SECTION FROM PROVIDER :>'+chem_provider+'\n\n'+commentary+'\n'
               dumpedrdkit+=1
               CAS_nb = None
               catalog_nb = None
               outputstring = ''
               continue
##############################################################################################
           if ' resin' in outputstring:
                  dumpedresin+=1
                  resineflag = True
           if sdmol != None and resineflag==False:
               outputstring = ''
               for sdmol in [sdmol,]:
                    rdkit.Chem.rdmolops.AssignAtomChiralTagsFromStructure(sdmol)
# Using RDKit, excludes all metal-containing molecules:
                    HeavyMetal = False
                    trimethyleneflag=False
                    for atom in sdmol.GetAtoms():
                        if atom.GetAtomicNum() in [13,20,31,32,33,49,50,51,81,82,83,84,21,22,23,24,25,26,27,28,29,30,39,40,41,42,43,44,45,46,47,48,72,73,74,75,76,77,78,79,80,104,105,106,107,108,109,110,111,112,57,58,59,60,61,62,63,64,65,66,67,68,69,70,7189,90,91,92,93,94,95,96,97,98,99,100,101,102,103]:
                            HeavyMetal = True
                    if HeavyMetal == True:
                        dumpedmetal+=1
                        continue
# Using RDKit, excludes all molecules with less than 3 or more than 35 non-hydrogen atoms:
                    if sdmol.GetNumHeavyAtoms() <= 3 or sdmol.GetNumHeavyAtoms() > 35:
                        dumpedHAC+=1
                        continue
# Using RDKit, excludes all molecules with more than 10 Rotatable Bonds:
                    if rdkit.Chem.Lipinski.NumRotatableBonds(sdmol) > 5:
                        dumpedrotatable+=1
                        continue
                    sdmolH=None
# Add Hydrogens to the molecular framwork:
                    try:
                        sdmolH=rdkit.Chem.rdmolops.AddHs(sdmol)
                    except:
                        sdmolH=None
                    if sdmolH!=None:          
# Using RDKit, excludes all molecules which contain the following substructures: CH2-Xsp3-CH2 , CH2-CH=CH-CH2 , S-S
                        for substruc in sdmolH.GetSubstructMatches(rdkit.Chem.MolFromSmarts('[CH2]-[CH2,O,NH]-[CH2]')):
                            if sdmolH.GetAtomWithIdx(substruc[1]).IsInRing():
                                pass
                            else:
                                trimethyleneflag=True
                        for substruc in sdmolH.GetSubstructMatches(rdkit.Chem.MolFromSmarts('[CH2]-[CH]=[CH]-[CH2]')):
                            if sdmolH.GetAtomWithIdx(substruc[1]).IsInRing() or sdmolH.GetAtomWithIdx(substruc[2]).IsInRing():
                                pass
                            else:
                                trimethyleneflag=True
                        if len(sdmolH.GetSubstructMatches(rdkit.Chem.MolFromSmarts('[#16]~[#16]'))) > 0:
                                trimethyleneflag=True
                    else:
                        dumpedrdkit+=1
                    if trimethyleneflag==True:
                        dumpedtrimethylene+=1
                        smile=rdkit.Chem.MolToSmiles(sdmol)
                        continue
# Using RDKit, excludes all molecules with more than 1 chiral center:
                    if len(rdkit.Chem.FindMolChiralCenters(sdmol, force=True, includeUnassigned=True)) > 1:
                        dumpedchiral+=1
                        continue
                    dic_data[CAS_nb] = limitedoutput+'\n\n> <CAS>'+CAS_nb
                    dico_provider[CAS_nb] ='\n\n> <PROVIDER>\n'+chem_provider+':'+catalog_nb
                    dico_commentary[CAS_nb]='\n\n> <SECTION FROM PROVIDER :>'+chem_provider+'\n\n'+commentary+'\n'
                    rddictioprices[CAS_nb]='\n\n\n> <PRICE SECTION FROM CATALOG COMPOUND '+chem_provider+':'+catalog_nb+'>'
                    for p in range(1,6):
                        if len(list(str(nddictioprices['PR'][p]))) > 0:
                            rddictioprices[CAS_nb]=rddictioprices[CAS_nb]+'\nPr.'+str(p)+'\t'+str(nddictioprices['PR'][p])+'\tCurr. '+currency+'\tQt.'+str(p)+'\t'+str(nddictioprices['QT'][p])+'\tunits. '+str(nddictioprices['UN'][p])
                        nddictioprices['PR'][p]=''
                        nddictioprices['QT'][p]=''
                        nddictioprices['UN'][p]=''    
                    list_CAS_nb.append(CAS_nb)
                    good+=1
                    CAS_nb = ''
                    catalog_nb = ''
           else:
                dumpedrdkit+=1
          limitedoutput=outputstring=commentary=''                          
          
failout=open('Failing_outputs','w')
for line in failist: 
    failout.write(line+'\n')
failout.close()
 
fileout= open('Global_Library_Reactants.sdf','w')
    
for CAS in list_CAS_nb:
      try:
        fileout.write(dic_data[CAS]+dico_provider[CAS]+rddictioprices[CAS]+'\n\n$$$$\n')
      except:
        fileout.write('\n$$$$\n\n'+CAS+'\n\n$$$$\n')

fileout.close()
print 'Molfile items screened %s' %advance
print 'Dumped resin:%s fastloop:%s rdkit:%s' %(dumpedresin,dumpedfastlooping,dumpedrdkit)
print 'Dumped metal:%s HAC:%s rotablebonds:%s CH2-X-CH2:%s chiral:%s'  %(dumpedmetal,dumpedHAC,dumpedrotatable,dumpedtrimethylene,dumpedchiral)
print 'Redundant via CAS %s, Stored %s' %(redundantCAS,good)
    
