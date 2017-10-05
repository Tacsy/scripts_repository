#How to use pathway analysis scripts

I create this subfolder for futhure use of pathway analysis using PATHWAYS plugin in VMD.
The following scripts are described in sequential order.

##pw_analysis.inp
In this slurm submission script, the working directory should be decleared first and then
a findpathways.tcl vmd script should also be defined. PDB files are in frame#.pdb format, and
the psf file is named as complex.psf

##parsepathway.inp
In this bash shell script, for each specific raw-output PATHWAYS file, the atom-formed pathways
are extracted and saved.

##findpathways.tcl
findpathways.tcl vmd script. One can change the number of pathways to determine in -p option,
and donor and acceptor index are atom number -1.

##findpathways.py
Transform atom-formed pathways to residue formed pathways for future analysis.

##analyzepathways.py
Statistically measure the frequency of every pathway in a trajectory.

