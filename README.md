# AutoCouple_Python-based

Python Scripts for AutoCouple

Author:  Laurent Batiste

Affiliation:  A. Caflisch' group at the Department of Biochemistry of the University of Zurich

Date:  October 31, 2017

- AutoCouple_Script_1.py retrieves building-blocks libraries from various chemicals providers (provided by user as inputs) and generates a uniquified library of building-blocks based on their CAS number (heavy metals-, resin-containing compounds are discarded).      The output is Global_Library_Reactants.sdf.

- AutoCouple_Script_2.py filters and sorts out each building-block based on the chemical functionalities it contains (for instance amines, boronic acids, halides...).      The output files are libraries of reactants called Rx_"name_of_the_functionality".sdf.

- AutoCouple_Script_3_Buchwald-Hartwig.py or AutoCouple_Script_3_Suzuki.py perform a virtual library of new ligand by coupling two libraries of reactants via either Buchwald-Hartwig or Suzuki, respectively.      The output file name is specified as last argument (+ ".sd").
