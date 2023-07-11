#*******#
#***Code written by Avijeet Kulshrestha*******
#structure-file is *.gro or *.tpr
#trajectory-file is *.xtc
#lipid name is the resname (CHL1, POPC, DOPE, etc.)- should match the force field
#cutoff in angstrom (float)
#Residue number in integer
#******

import MDAnalysis as MDA
import numpy as np
from MDAnalysis.analysis.distances import distance_array
import pandas as pd
import matplotlib.pyplot as plt

#structure=argv[1];traj=argv[2]; lipid=argv[3];cutoff=argv[4];residue=argv[5];

class leaflet:
    def __init__(self, structure):
        self.u=MDA.Universe(str(structure))
        
    def finder(self):
        #*****select all phosphorous atoms
        self.p=self.u.select_atoms('name P*')
        
        #**Center of mass of lipid's phosphorous atoms in the system
        center=self.p.center_of_mass() 
        c_z=center[2]
        upper=self.p.select_atoms('prop z>' +str(c_z)) 
        lower=self.p.select_atoms('prop z<'+str(c_z))
        
        upper.residues.atoms.write('upper_leaflet.ndx')
        lower.residues.atoms.write('lower_leaflet.ndx')
        
        return("upper_leaflet.ndx contains upper leaflet atoms and lower_leaflet.ndx conpatins lower leaflet atoms")
        
        
class post_processing:
    def __init__(self, structure, traj):
        self.structure=structure
        self.traj=traj
        self.u=MDA.Universe(str(structure), str(traj))
        
    def time_occupancy(self, res, lip, cut):
        self.residue = res
        self.lipid = lip
        self.cutoff=float(cut)
        print ("Occupany is for residue:", self.residue ,", and lipid:", self.lipid,", with cutoff value:", self.cutoff,"angstrom")
        self.p=self.u.select_atoms("resid "+ str(self.residue))
        sele_lipid=self.u.select_atoms("resname "+self.lipid)
        
        tot_avg=[];
        for tr in self.u.trajectory:        
            s2=sele_lipid.residues[0];
            lipid_id="resid "+str(s2.resid); #+" and not name H*"
            s3=self.u.select_atoms(lipid_id)
            #Distance 2D-array
            di=distance_array(self.p.positions, s3.positions, box=self.u.dimensions)
            tot=len(np.where( di <= self.cutoff )[0])
            
            for k in range(1,len(sele_lipid.residues)):
                s2=sele_lipid.residues[k]#.select_atoms(sel_atoms)
                lipid_id="resid "+str(s2.resid) #+" and not name H*"
                s3=self.u.select_atoms(lipid_id);
                #Distance 2D-array
                di2=distance_array(self.p.positions, s3.positions, box=self.u.dimensions)
                tot=tot+len(np.where( di2 <= self.cutoff )[0])
            
            tot_avg.append(np.mean(tot))
        
        self.time=tot_avg
        df=pd.DataFrame()
        df['occupancy']=tot_avg;
        df.to_csv('occupancy_time.dat',sep=' ', index=True, header=True)
        return ('occupancy_time.dat file contains the occupancy data with frames')
    
    def avg_occupancy(self, res, lip, cut):
        #time_occupancy(res,lip,cut)
        return ('Average and standard deviation:',np.mean(self.time), np.std(self.time))
        
    def apl(self, lipids):
        apl=[]
        for ts in self.u.trajectory:
            apl.append((self.u.dimensions[0]*self.u.dimensions[1]))
        
        apl=np.array(apl); apl=apl/(lipids) #apl is the area per lipid
        df=pd.DataFrame()
        df['APL']=apl;
        df.to_csv('APL.dat',sep=' ', index=True, header=True)
        return ('APL.dat file contains area per lipid change with frames')
        
    def depth(self, initial_residue, final_residue):
        res=np.arange(initial_residue,final_residue+1,1) #Residue ids
        d=[] #This is to add mean depth values
        std=[] #This is to add standard deviation
        for i in range(initial_residue,final_residue+1): #length of residues 
            m=self.u.select_atoms("not protein and not name SOL")
            p=self.u.select_atoms("resid "+str(i))
            dr=[]
            for ts in self.u.trajectory:
                m_com=m.center_of_mass(); 
                p_com=p.center_of_mass(); 
                dr.append(abs(m_com[2]-p_com[2]))
                #dr.append(p.positions[0][2])   
            d.append(np.average(dr)/10)
            std.append(np.std(dr)/10)
            del dr[:]
            
        df=pd.DataFrame()
        df['residue']=res; df['depth_mean']=d; df['depth_std']=std;
        df.to_csv('residue_depth.dat',sep=' ', index=True, header=True)
        return ('residue_depth.dat file contains the residue-wise depth in the membrane')
        
    def e2e_dist(self,initial, final):
        from numpy.linalg import norm
        dist=[] #This is to add end to end distance
        for ts in self.u.trajectory:
            first=self.u.select_atoms('resid '+ str(initial) +' and name CA')
            last=self.u.select_atoms('resid '+str(final)+' and name CA')
            
            dist.append(norm(first.positions-last.positions))
            
        df=pd.DataFrame()
        df['e2e_dist']=dist; 
        df.to_csv('e2e_dist.dat',sep=' ', index=True, header=True)
        return ('e2e_dist.dat file contains the end to end distance')
        
    def lateral_mobility(self):
        p=self.u.select_atoms("protein")
        p_com=p.center_of_geometry()
        xi=p_com[0]; yi=p_com[1]
        m=[] #variable to add lateral mobility data
        for ts in self.u.trajectory:
            p=self.u.select_atoms("protein")
            p_com=p.center_of_geometry()
            xn=p_com[0]; yn=p_com[1]
            m.append(np.sqrt(pow((xn-xi),2)+pow((yn-yi),2)))
            xi=xn; yi=yn 
            
        df=pd.DataFrame()
        df['lateral_mobility']=m; 
        df.to_csv('lateral_mobility.dat',sep=' ', index=True, header=True)
        return ('lateral_mobility.dat file contains the lateral mobility of protein')
        
    def tilt_angle(self,initial, final):
        t1=[]
        for ts in self.u.trajectory:
            f=self.u.select_atoms('resid '+ str(initial) +' and name CA')
            f_x=f.positions[0][0]
            f_y=f.positions[0][1]
            f_z=f.positions[0][2]
            
            l=self.u.select_atoms('resid '+str(final)+' and name CA')
            l_x=l.positions[0][0]
            l_y=l.positions[0][1]
            l_z=l.positions[0][2]
            
            a1=(np.arccos((f_z-l_z)/np.sqrt(pow((f_x-l_x),2)+pow((f_y-l_y),2)+pow((f_z-l_z),2))))
            t1.append(a1)
        
        t1 = [i * 57.2958 for i in t1] #Tilt angle array
        df=pd.DataFrame()
        df['tilt_angle']=t1; 
        df.to_csv('tilt_angle.dat',sep=' ', index=True, header=True)
        return ('tilt_angle.dat file contains the tilt angle (degree) between +Z and a vector from inital to final residue')
        
    def hydration(self,molecule, cutoff):
        shell=[] #Array to add numbers in hydration shell
        for ts in self.u.trajectory:
            #Below select is for water, lipids molecules can a;so be selected
            w=self.u.select_atoms('resname '+str(molecule)+ " and around "+str(cutoff)+" protein")
            shell.append(len(w.residues))
            
        df=pd.DataFrame()
        df['tilt_angle']=shell; 
        df.to_csv('hydration.dat',sep=' ', index=True, header=True)
        return ('hydration.dat file contains the hydration of given molecule around the specified cutoff')
        
    def contact_map(self,residues):
        import MDAnalysis.analysis.distances as dist
        d = np.zeros(shape=(residues,residues)) #38 is the number of residues
        for ts in self.u.trajectory:
            #C-alpha contact map of 1 nm of cutoff
            res=self.u.select_atoms('protein and name CA')
            d1=dist.contact_matrix(res.positions, cutoff=10.0, returntype='numpy', box=self.u.dimensions)
            d1=np.where(d1=='True', 1, d1)
            d=d+d1
            
        d=np.divide(d,len(self.u.trajectory)); 
        
        #**Adding a condition for upper triangle to avoid redundancy of the symmetry
        triu=np.triu_indices(38); #print(triu)
        triu_a=triu[0]; #triu_a=triu_a[1:-1]
        triu_b=triu[1]; #triu_b=triu_b[1:-1]
        for i in range(0,len(triu_a)):
            if(triu_a[i]==triu_b[i]):
                d[triu_a[i]][triu_b[i]] = d[triu_a[i]][triu_b[i]]
            else:
                d[triu_a[i]][triu_b[i]] = 0
                
        df=pd.DataFrame(d)
        df.to_csv('contact_map.dat',sep=' ', index=True, header=True)
        return ('contact_map.dat file contains the data')
        
    def membrane_surface(self, box_x, box_y):
        from membrane_curvature.base import MembraneCurvature
        from scipy import ndimage
        from matplotlib import ticker
        #***Select phosphorous atoms from all types of lipids
        p=self.u.select_atoms('name P*');
        center=p.center_of_mass()        
        #**Identfying the upper and lower leaflets
        c_z=center[2]        
        mc_upper = MembraneCurvature(p, select='prop z>' +str(c_z), n_x_bins=6, n_y_bins=6).run()
        mc_lower = MembraneCurvature(p, select='prop z<' +str(c_z), n_x_bins=6, n_y_bins=6).run()
        
        # surface
        surf_upper = mc_upper.results.average_z_surface
        surf_lower = mc_lower.results.average_z_surface
        
        #*****Interpolating
        # Resample your data grid by a factor of 3 using cubic spline interpolation.
        surf_upper = ndimage.zoom(surf_upper, 3); surf_upper=surf_upper/10
        surf_lower = ndimage.zoom(surf_lower, 3); surf_lower=surf_lower/10
        
        phi=np.linspace(0,box_x,surf_upper.shape[0])
        psi=np.linspace(0,box_y,surf_upper.shape[1])
        phi, psi= np.meshgrid(phi,psi)
        
        #****Calculating center and subtracting
        c_z=[]
        for tr in self.u.trajectory:
            p=self.u.select_atoms('resname POPC and name P')
            center=(p.center_of_mass())
            c_z.append(center[2])
        
        center=np.mean(c_z)/10
        surf_upper=surf_upper-center
        surf_lower=surf_lower-center
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        surf_u=ax.plot_surface(phi, psi, surf_upper, cmap="bwr", vmin=-3.0, vmax=3.0)
        surf_l=ax.plot_surface(phi, psi, surf_lower, cmap="bwr", vmin=-3, vmax=3)
        cbar1=fig.colorbar(surf_l, shrink=0.5, aspect=10, location='right')
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar1.locator = tick_locator
        cbar1.update_ticks()
        cbar1.ax.tick_params(labelsize=10)
        cbar1.ax.set_ylabel(r'$Z_{center}-Z_{leaflet}$ ($nm$)', fontsize=10, fontname='Arial')
        
        ax.w_zaxis.line.set_lw(0.)
        ax.set_zticks([])
        ax.set_xlabel(r"X(nm)", fontsize=12, fontname='Arial')
        ax.set_ylabel(r"Y(nm)", fontsize=12, fontname='Arial')
        plt.savefig('mean-surface.png', bbox_inches = 'tight', pad_inches = 0.1,dpi=300)
        
    def membrane_thickness(self,box_x, box_y):
        from membrane_curvature.base import MembraneCurvature
        from scipy import ndimage
        from matplotlib import ticker
        #***Select phosphorous atoms from all types of lipids
        p=self.u.select_atoms('name P*'); 
        center=p.center_of_mass()        
        #**Identfying the upper and lower leaflets
        c_z=center[2]        
        mc_upper = MembraneCurvature(p, select='prop z>' +str(c_z), n_x_bins=6, n_y_bins=6).run()
        mc_lower = MembraneCurvature(p, select='prop z<' +str(c_z), n_x_bins=6, n_y_bins=6).run()
        
        # surface
        surf_upper = mc_upper.results.average_z_surface
        surf_lower = mc_lower.results.average_z_surface
        
        #*****Interpolating
        # Resample your data grid by a factor of 3 using cubic spline interpolation.
        surf_upper = ndimage.zoom(surf_upper, 3); surf_upper=surf_upper/10
        surf_lower = ndimage.zoom(surf_lower, 3); surf_lower=surf_lower/10
        
        thickness= surf_upper- surf_lower

        #***Define the X and Y range
        X=np.linspace(0,box_x,thickness.shape[0])
        Y=np.linspace(0,box_y,thickness.shape[1])
        X, Y= np.meshgrid(X,Y)
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.patch.set_edgecolor('black')  ##Border color
        ax.patch.set_linewidth(2.5)  #Border linewidth
        levels=np.linspace(1,5,20)
        plt.contourf(X,Y,thickness, cmap='PRGn', levels=levels)
        cbar=plt.colorbar()
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()
        
        cbar.ax.tick_params(labelsize=20)
        cbar.ax.set_ylabel(r'Thickness (nm)', fontsize=20, fontname='Arial')
        
        plt.locator_params(axis='y', nbins=5)
        plt.locator_params(axis='x', nbins=5)
        
        ax.set_xlabel(r"X(nm)", fontsize=20, fontname='Arial')
        ax.set_ylabel(r"Y(nm)", fontsize=20, fontname='Arial')
        plt.xticks(fontsize=20, fontname='Arial') #Xtick size
        plt.yticks(fontsize=20, fontname='Arial') #Ytick size
        plt.savefig('thickness.png', bbox_inches = 'tight', pad_inches = 0.1,dpi=300)
        
    def membrane_mean_curvature(self, box_x, box_y):
        from membrane_curvature.base import MembraneCurvature
        from scipy import ndimage
        from matplotlib import ticker
        #***Select phosphorous atoms from all types of lipids
        p=self.u.select_atoms('name P*');
        center=p.center_of_mass()        
        #**Identfying the upper and lower leaflets
        c_z=center[2]        
        mc_upper = MembraneCurvature(p, select='prop z>' +str(c_z), n_x_bins=6, n_y_bins=6).run()
        mc_lower = MembraneCurvature(p, select='prop z<' +str(c_z), n_x_bins=6, n_y_bins=6).run()
        
        # mean curvature
        mean_upper = mc_upper.results.average_mean
        mean_lower = mc_lower.results.average_mean
        
        #******UPPER leaflet********
        # Resample your data grid by a factor of 3 using cubic spline interpolation.
        mean_upper = ndimage.zoom(mean_upper, 3)
        phi=np.linspace(0,box_x,mean_upper.shape[0])
        psi=np.linspace(0,box_y,mean_upper.shape[1])
        phi, psi= np.meshgrid(phi,psi)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.patch.set_edgecolor('black')  ##Border color
        ax.patch.set_linewidth(2.5)  ##Border linewidth
        levels=np.linspace(-1,1,20)
        plt.contourf(phi,psi,mean_upper, 20, cmap='bwr', levels=levels)
        cbar=plt.colorbar()
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()
        cbar.ax.tick_params(labelsize=20)
        cbar.ax.set_ylabel(r'Mean curvature ($nm^{-1}$)', fontsize=20, fontname='Arial') 
        plt.locator_params(axis='y', nbins=5)
        plt.locator_params(axis='x', nbins=5)
        ax.set_xlabel(r"X(nm)", fontsize=20, fontname='Arial')
        ax.set_ylabel(r"Y(nm)", fontsize=20, fontname='Arial')
        plt.xticks(fontsize=20, fontname='Arial') #Xtick size
        plt.yticks(fontsize=20, fontname='Arial') #Ytick size
        plt.savefig('mean-upper-curvature.png', bbox_inches = 'tight', pad_inches = 0.1,dpi=300)
        plt.clf()
        
        #******Lower leaflet********
        # Resample your data grid by a factor of 3 using cubic spline interpolation.
        mean_lower = ndimage.zoom(mean_lower, 3)
        phi=np.linspace(0,box_x,mean_lower.shape[0])
        psi=np.linspace(0,box_y,mean_lower.shape[1])
        phi, psi= np.meshgrid(phi,psi)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.patch.set_edgecolor('black')  ##Border color
        ax.patch.set_linewidth(2.5)  ##Border linewidth
        levels=np.linspace(-1,1,20)
        plt.contourf(phi,psi,mean_lower, 20, cmap='bwr', levels=levels)
        cbar=plt.colorbar()
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()
        cbar.ax.tick_params(labelsize=20)
        cbar.ax.set_ylabel(r'Mean curvature ($nm^{-1}$)', fontsize=20, fontname='Arial') 
        plt.locator_params(axis='y', nbins=5)
        plt.locator_params(axis='x', nbins=5)
        ax.set_xlabel(r"X(nm)", fontsize=20, fontname='Arial')
        ax.set_ylabel(r"Y(nm)", fontsize=20, fontname='Arial')
        plt.xticks(fontsize=20, fontname='Arial') #Xtick size
        plt.yticks(fontsize=20, fontname='Arial') #Ytick size
        plt.savefig('mean-lower-curvature.png', bbox_inches = 'tight', pad_inches = 0.1,dpi=300)
        plt.clf()
        
    def membrane_gaussian_curvature(self, box_x, box_y):
        from membrane_curvature.base import MembraneCurvature
        from scipy import ndimage
        from matplotlib import ticker
        #***Select phosphorous atoms from all types of lipids
        p=self.u.select_atoms('name P*');
        center=p.center_of_mass()        
        #**Identfying the upper and lower leaflets
        c_z=center[2]        
        mc_upper = MembraneCurvature(p, select='prop z>' +str(c_z), n_x_bins=6, n_y_bins=6).run()
        mc_lower = MembraneCurvature(p, select='prop z<' +str(c_z), n_x_bins=6, n_y_bins=6).run()
        
        # Gaussian curvature
        gauss_upper = mc_upper.results.average_gaussian
        gauss_lower = mc_lower.results.average_gaussian
        
        #******UPPER leaflet********
        # Resample your data grid by a factor of 3 using cubic spline interpolation.
        gauss_upper = ndimage.zoom(gauss_upper, 3)
        phi=np.linspace(0,box_x,gauss_upper.shape[0])
        psi=np.linspace(0,box_y,gauss_upper.shape[1])
        phi, psi= np.meshgrid(phi,psi)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.patch.set_edgecolor('black')  ##Border color
        ax.patch.set_linewidth(2.5)  ##Border linewidth
        levels=np.linspace(-1,1,20)
        plt.contourf(phi,psi,gauss_upper, 20, cmap='bwr', levels=levels)
        cbar=plt.colorbar()
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()
        cbar.ax.tick_params(labelsize=20)
        cbar.ax.set_ylabel(r'Gaussian curvature ($nm^{-2}$)', fontsize=20, fontname='Arial') 
        plt.locator_params(axis='y', nbins=5)
        plt.locator_params(axis='x', nbins=5)
        ax.set_xlabel(r"X(nm)", fontsize=20, fontname='Arial')
        ax.set_ylabel(r"Y(nm)", fontsize=20, fontname='Arial')
        plt.xticks(fontsize=20, fontname='Arial') #Xtick size
        plt.yticks(fontsize=20, fontname='Arial') #Ytick size
        plt.savefig('gaussian-upper-curvature.png', bbox_inches = 'tight', pad_inches = 0.1,dpi=300)
        plt.clf()
           
        #******Lower leaflet********
        # Resample your data grid by a factor of 3 using cubic spline interpolation.
        gauss_lower = ndimage.zoom(gauss_lower, 3)
        phi=np.linspace(0,box_x,gauss_lower.shape[0])
        psi=np.linspace(0,box_y,gauss_lower.shape[1])
        phi, psi= np.meshgrid(phi,psi)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.patch.set_edgecolor('black')  ##Border color
        ax.patch.set_linewidth(2.5)  ##Border linewidth
        levels=np.linspace(-1,1,20)
        plt.contourf(phi,psi,gauss_lower, 20, cmap='bwr', levels=levels)
        cbar=plt.colorbar()
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()
        cbar.ax.tick_params(labelsize=20)
        cbar.ax.set_ylabel(r'Mean curvature ($nm^{-2}$)', fontsize=20, fontname='Arial') 
        plt.locator_params(axis='y', nbins=5)
        plt.locator_params(axis='x', nbins=5)
        ax.set_xlabel(r"X(nm)", fontsize=20, fontname='Arial')
        ax.set_ylabel(r"Y(nm)", fontsize=20, fontname='Arial')
        plt.xticks(fontsize=20, fontname='Arial') #Xtick size
        plt.yticks(fontsize=20, fontname='Arial') #Ytick size
        plt.savefig('gaussian-lower-curvature.png', bbox_inches = 'tight', pad_inches = 0.1,dpi=300)
        plt.clf()
                
                
        

#**********input****************** 
#ssdump.day number-of-residues total-number-of-frames
class secondary_structure:
    def __init__(self, filename):
        self.file=str(filename)
        
    def time_ss(self, residue, frames):
        df=pd.read_fwf(self.file, widths=[1] * residue, header=None)
        df=df.T
        rows=['~','H','B','E','G','I','T','S']
        prob=pd.DataFrame(index=rows)
        p1=df[0].value_counts(normalize=True)
        result = pd.concat([prob, p1], axis=1, sort=False)
        
        for i in range(1,frames): #frames is the number of row in ssdump file
            p=df[i].value_counts(normalize=True)
            result = pd.concat([result, p], axis=1, sort=False)
        
        result=result.fillna(0); result=result.T
        result.to_csv('secondary_structure_time.dat',sep=' ', index=True, header=True)
        return('secondary_structure_time.dat file contains the secondary structure change with time for full protein')
        
    def residue_ss(self, residue):
        df=pd.read_fwf(self.file, widths=[1] * residue, header=None)
        rows=['~','H','B','E','G','I','T','S']
        prob=pd.DataFrame(index=rows)
        
        p1=df[0].value_counts(normalize=True)
        result = pd.concat([prob, p1], axis=1, sort=False)
        
        for i in range(1,residue):
            p=df[i].value_counts(normalize=True)
            result = pd.concat([result, p], axis=1, sort=False)
        
        result=result.fillna(0)
        result=result.T
        #print(result['H'])
        #plt.plot(result['H'])
        
        result.to_csv('secondary_structure_residue.dat',sep=' ', index=True, header=True)
        return('secondary_structure_residue.dat file contains the average secondary structure for each residue')



