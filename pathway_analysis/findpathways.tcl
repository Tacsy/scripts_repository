mol load psf complex.psf pdb frame0.pdb
package require pathways
pathways -d "index 6377" -a "index 6304" -p "10" -withh "1" -exph "1.146" -expts "1.146" -cda "1" -o pw_frame0

exit
