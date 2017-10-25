# AutoCouple_Python-based
Python Scripts used within the AutoCouple methodology developped by the groups of Prof. Caflisch and Prof. Nevado

A report on AutoCouple's application to the development of Bromodomain binders is in press at ACS Central Science.

The methodology requires the use of three Python-encoded scripts:

- the first type (AutoCouple_Script_1.py) retreaves building-blocks libraries from various chemicals providers (provided by user as inputs) and generates a uniquified library of building-blocks based on their CAS number (heavy metals-, resin-containing compounds are discarded). 
The output is Global_Library_Reactants.sdf.


- the second type (AutoCouple_Script_2.py) filters and sorts out each building-block based on the chemical functionalities it contains ( for instance amines, boronic acids, halides...). 
The output files are libraries of reactants called Rx_"name_of_the_functionality".sdf.


- the third type of script (AutoCouple_Script_3_Buchwald-Hartwig.py or AutoCouple_Script_3_Suzuki.py) perform a virtual library of new ligand by coupling two libraries of reactants via either Buchwald-Hartwig or Suzuki. 
The output file name is specified as last argument (+ ".sd").
