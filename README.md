# DIST
Generate atomic structures for common defects in materials

One of the very powerful aspects is the colorful examples. Before you start to use the toolkit, please first check if the examples include the structures that you want. For examples, we already have:
(1) the {111}<11-2>/6 edge dislocation for fcc;
(2) the <111> screw dislocation for bcc;
(3) the magic angle bilayer graphene (1.1);
....

If you cannot find the structures that meets your requirements, you can follow the howto.dat file in the example folders to creat them by yourself.

This is the first part of the DIST toolkit. The second part is about multi-scale modeling of dislocations and strains. This part will be online soon.

Please write me if you have problems with the toolkit. The email address can be found from the DIST paper below.

If you find this toolkit is helpful to you, please kindly cite the following paper:
Zongrui Pei, DIST: A dislocation-simulation toolkit, Computer Physics Communications 233(2018)44-50.

---------------bibtex----------------------------

@article{PEI201844,
title = "DIST: A dislocation-simulation toolkit",
journal = "Computer Physics Communications",
volume = "233",
pages = "44 - 50",
year = "2018",
issn = "0010-4655",
doi = "https://doi.org/10.1016/j.cpc.2018.06.021",
url = "http://www.sciencedirect.com/science/article/pii/S0010465518302297",
author = "Zongrui Pei",
keywords = "Dislocation, Stacking fault, Simulation, Toolkit",
abstract = "Dislocations are important defects determining the mechanical properties in metallic materials for structural applications. In order to simulate dislocations, either generalized stacking faults (GSFs) or atomistic dislocation structures need to be constructed firstly. So far there is a lack of light toolkit to help easily generate supercells with GSFs or dislocations. Easy-to-use independent tools to generate any screw or edge dislocations in any crystal structures, will reduce the barrier for dislocation simulations and probably attract more beginners to embark on such mechanical-properties related simulations. Based on an effective algorithm we develop a toolkit for dislocation simulations using Python, a computer language that is readily available in almost all operation systems. Besides, this toolkit also includes tools for multi-scale Peierls–Nabarro modeling of dislocations with gamma surfaces as the key input information. This part is written in C++ language."
}


#-----------------generation of a screw dislocation--------------------

The code is not smart enough. You have to (i) use cartesian coordinate and (ii) make sure the supercell size is the same as one with dislocation. It cannot magnify and generate dislocation at the same time. Assume you prepare your bcc_super_cart, dislocation will be inserted without changing its size. You can do the following to get a dislocation:
python ../../../gen_screwDislocation.py bcc_super_cart > bcc.screw 
A screw dislocation can be found in bcc_screw
You can find these files in examples/BCC/bcc_screw_dislocation/
Please do not foget to cite my paper if you find the toolkit is useful. Thanks!

