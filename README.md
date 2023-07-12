# MDProcessing
MDProcessing is written in Python for pre and post-processing of MD trajectories.
Packages required: MDAnalysis, Numpy, Pandas, Matplotlib

Analysis that can be done:

import MDProcessing as MDP

1. To find the upper and lower leaflet of a membrane:
   struct=MDP.leaflet('structure.gro'); #gro file is the structure file
   print(struct.finder());
   
2. To find the change in secondary structure content with time:
   First, generate "ssdump.dat" file using Gromacs do_dssp command,
   ss=MDP.secondary_structure('ssdump.dat');
   print(ss.time_ss(number-of-residues, number of frames));
   
3. To find the average change in secondary structure content for each residue:
   First, generate "ssdump.dat" file using Gromacs do_dssp command,
   ss=MDP.secondary_structure('ssdump.dat');
   print(ss.residue_ss(number-of-residues));

4. To calculate the area per lipid of the membrane:
   res1=MDP.post_processing('structure.gro', 'trajectory.xtc');
   print(res1.apl(number-of-lipids));

5. To calculate the residue-wise depth of protein in the membrane:
   depth=MDP.post_processing('structure.gro', 'trajectory.xtc');
   print(depth.depth(initial-residue,final-residue)); # initial and final residue to define the range

6. To compute the end-to-end distance of the protein:
   e2e=MDP.post_processing('structure.gro', 'trajectory.xtc');
   print(e2e.e2e(initial-residue,final-residue));

7. To compute the lateral mobility of protein in the membrane:
   lm=MDP.post_processing('structure.gro', 'trajectory.xtc');
   print(lm.lateral_mobility());

8. To compute the tilt angle of protein in the membrane:
   titl=MDP.post_processing('structure.gro', 'trajectory.xtc');
   print(tilt.tilt_angle(initial-residue,final-residue)); #initial and final residue to make a vector

9. To compute the hydration shell around the protein:
   hydrate=MDP.post_processing('structure.gro', 'trajectory.xtc');
   print(hydrate.hydration('molecule-name',cut-off));

10. To calculate the contact map:
    cmap=MDP.post_processing('structure.gro', 'trajectory.xtc');
    print(cmap.contact_map(number-of-residues));

11. To plot the upper and lower membrane surfaces: This requires MembraneCurvature python package
    mem=MDP.post_processing('structure.gro', 'trajectory.xtc');
    print(mem.membrane_surface(box-x,box-y)); #box-x and box-y is the x and y length of box

12. To plot the 2D membrane thickness: This requires MembraneCurvature python package
    mem=MDP.post_processing('structure.gro', 'trajectory.xtc');
    print(mem.membrane_thickness(box-x,box-y)); #box-x and box-y is the x and y length of box

13. To compute the mean curvature: This requires MembraneCurvature python package
    mem=MDP.post_processing('structure.gro', 'trajectory.xtc');
    print(mem.membrane_mean_curvature(box-x,box-y)); #box-x and box-y is the x and y length of box

14. To compute the Gaussian curvature: This requires MembraneCurvature python package
    mem=MDP.post_processing('structure.gro', 'trajectory.xtc');
    print(mem.membrane_gaussian_curvature(box-x,box-y)); #box-x and box-y is the x and y length of box

15. To compute the autocorrelation of data:
    AC=MDP.autocorrelation();
    print(AC.autocorr(y)); #y is an array

16. To compute the occupancy change with time of a molecule on a protein residue:
    res1=MDP.post_processing('structure.gro', 'trajectory.xtc');
    print(res1.time_occupancy(residue-id, 'molecule-name', cutoff)); #This will calculate the 'molecule-name' occupancy on the protein residue-id within a given cutoff.

17. To compute the average occupancy of a molecule on a protein residue:
    res1=MDP.post_processing('structure.gro', 'trajectory.xtc');
    print(res1.avg_occupancy(residue-id, 'molecule-name', cutoff));


Let me know if you find any bugs or want me to add additional analysis. You can write to me at avijeetkulshrestha@gmail.com

If you are using the written code, please cite the following papers:
1. A. Kulshrestha, S. Maurya, T. Gupta, R. Roy, S. N. Punnathanam, K. G. Ayappa, "Conformational flexibility is a key determinant for the lytic activity of the pore-forming protein, Cytolysin A” J. Phys. Chem. B (2022), 127 (1), 69-84
2. A. Kulshrestha, S. N. Punnathanam, K. G. Ayappa, "Finite temperature string method with umbrella sampling using path collective variables: Application to secondary structure change in a protein” Soft Matter (2022), 18 (39), 7593-7603
3. A. Kulshrestha, S. N. Punnathanam, R. Roy, K. G. Ayappa, "Cholesterol catalyzes unfolding in membrane inserted motifs of the pore forming protein cytolysin A" BioRxiv

Thanks,
Avijeet Kulshrestha

    












