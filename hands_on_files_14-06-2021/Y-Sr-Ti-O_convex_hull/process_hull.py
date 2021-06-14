import sys
import os
import pandas
import re
import glob
import shutil
from ase.io import *
from ase.visualize import *
from decimal import *
import pickle
from fuse_stable import run_gulp

### script for processing results from FUSE when the lowest n energy structures have
### been extracted and recomputed using VASP, with the presence of structures from ICSD
### script will then produce a folder "pymatgen" which should contain all of the input 
### files for pymat gen for each composition from FUSE plus and extras from ICSD
### will then try to add function to then execute & process the results from pymatgen
#
# Expeected data format
#		current working dir (inc. location of this file)
#				    |	
#                                   |
#	 ----------------------------------------------------
#        |                                                   |
#        |                                                   |
#    ncomp folders                                 (optional) icsd folder
#        |                                                   |
#        |                                                   |
# one folder per FUSE calc                             one folder per icsd 
#                                                      VASP calc
#
# where ncomps indicates the number of compositions for probe structures
# 
######## inputs ################################################################

comps=[] # list of the directory names for the compositions (we'll fill this out later)
icsd = 'ICSD' # is there a folder containing relaxed ICSD / reference structures?, leave as '' if not
norm_to=['Sr','Ti','Y'] # atomic species to use with pymatgen
all_elements=['Sr','Ti','Y','O'] # All of the elements present in the system
system_python_command='python' # path to python exe to use for pymatgen
exclude=[icsd,"pymatgen"]

################################################################################
# use GULP to get energies from cifs
kwds=['opti conj conp c6','opti conp c6'] # keywords for the gulp input, 
#multiple strings indicate running GULP in multiple stages
gulp_opts=[
['\nlibrary lib2.lib\ndump temp.res\ntime 15 minuets\nmaxcyc 50\nstepmax 0.1\ntime 5 minutes'],
['\nlibrary lib2.lib\ndump temp.res\ntime 15 minutes\nmaxcyc 1500\nlbfgs_order 5000\ntime 5 minutes'],
]	# options for gulp, must include one line per set of inputs in "kwds"
lib='lib2.lib' # library file for interatomic potentials
shel=['']	# species using shells in gulp

################################################################################

################################################################################
#work out the number of compositions
ncomp=0
files=glob.glob("*")
for i in range(len(files)):
	if os.path.isdir(files[i]):
		if not files[i] in exclude:
			ncomp+=1
			comps.append(files[i])

################################################################################
compositions=[]
entries_for_pymatgen={}
data_table={}

#create empty dictionary which we will later convert into a pandas data frame
# which we will then use to create an csv file containing the hull and / or 
# plot a colour map from

convex_hull={}
for i in range(len(all_elements)):
	convex_hull[all_elements[i]]=[]

convex_hull['energy from hull (meV/atom)']=[]


# check if pymatgen directory exisits & create it

if not os.path.isdir("pymatgen"):
	os.mkdir("pymatgen")

for i in range(ncomp):
	
	os.chdir(comps[i])
	shutil.copy(str("../"+lib),".")
	print ("\ncomposition: " + comps[i])
	# extract the raw energy from the FUSE outputs
	eng=0
	runs=[]
	
	files=glob.glob("*.cif")
	engs=[] 
	for i in files:
		atoms = read(i)
		atoms,energy,converged=run_gulp(atoms=atoms,shel=shel,kwds=kwds,opts=gulp_opts,lib=lib)
		engs.append(energy/len(atoms))
	#print(min(engs))	
	energy=min(engs)
	min_engs=engs.index(min(engs))
	#print(min_engs)
	## now create the entry in the entries for pymatgen dictionary
	#print(os.getcwd())
	#print(runs[min_engs])
	atoms=read(files[min_engs])
	syms=atoms.get_chemical_symbols()
	#print(syms)
	#print(atoms.get_chemical_formula())
	counts={}
	for j in range(len(all_elements)):
		if syms.count(all_elements[j]) != 0:
			counts[all_elements[j]]=syms.count(all_elements[j])
	
	keys=list(counts.keys())
	values=list(counts.values())
	
	values_norm=[]
	for j in range(len(norm_to)):
		if syms.count(norm_to[j]) != 0:
			values_norm.append(syms.count(norm_to[j]))
	
	
	norm=min(values_norm)
	for j in range(len(values)):
		if values[j] != 0:
			values[j]=values[j]/norm
	
	for j in range(len(values_norm)):
		if values_norm[j] != 0:
			values_norm[j]=values_norm[j]/norm
		
	# check that the normalisation produced formlulas containing all integer values
	# if not, correct it.
	
	all_integers=True
	indicies=[]
	for j in range(len(values)):
		value=values[j]
		test=float(value).is_integer()
		if test == False:
			all_integers=False
			indicies.append(j)
	
	if all_integers == False:
		#print(indicies)
		import math
		norms=[]
		for j in range(len(indicies)):
			num=values[indicies[j]]
			frac,whole=math.modf(num)
			#print(frac,whole)
			#print(1/frac)
			for k in range(len(values)):
				if frac != 0:
					values[k]=Decimal(values[k]*(1/frac))
					values[k]=float(values[k].quantize(Decimal('1e-6')))
			all_integers = True
			
			for k in range(len(values)):
				#print("hello")
				test=values[k].is_integer()
				if test == False:
					all_integers = False

			while all_integers == False:
				
				for j in range(len(values)):
					value=values[j]
					test=float(value).is_integer()
					if test == False:
						all_integers=False
						indicies.append(j)
				
				num=values[indicies[j]]
				frac,whole=math.modf(num)
				if frac != 0:
					#print("frac whole: ", frac,whole)
					#print("1/frac: ", 1/frac)
					for k in range(len(values)):
						values[k]=Decimal(values[k]*(1/frac))
						values[k]=float(values[k].quantize(Decimal('1e-6')))
					for k in range(len(values)):
						test=values[k].is_integer()
						if test == True:
							all_integers = True

	
	for j in range(len(keys)):
		counts[keys[j]]=values[j]

	## now create the entry in the entries for pymatgen dictionary
	
	#print("counts ",counts)
	#energy=energy*(sum(values))
	#print(energy)
	comp_label=''
	for j in range(len(keys)):
		comp_label += str(keys[j])
		comp_label += str(int(values[j]))
	#print (comp_label)
	
	norm_label=''
	norm_counts={}
	for j in range(len(norm_to)):
		if norm_to[j] in counts.keys():
			norm_label+=str(norm_to[j])
			norm_label+=str(int(counts[norm_to[j]]))
			norm_counts[norm_to[j]]=counts[norm_to[j]]
			
	#print(norm_label)
	#print(norm_counts)
	#print(energy)
	#print(counts)
	#print(values)
	norm_values=list(norm_counts.values())
	energy = energy * sum(values)
	#print(energy)
	#sys.exit()
	
	test = comp_label in list(entries_for_pymatgen.keys())
	if test == False:
		entries_for_pymatgen[str(comp_label)]={'energy_per_norm_FU':[energy],'pymat_comp':{norm_label:norm_counts},'comp_dict':counts}
	if test == True:
		entries_for_pymatgen[str(comp_label)]['energy_per_norm_FU'].append(energy)
	
	#print(os.getcwd())
	
	os.chdir("../")

#temp_keys=list(entries_for_pymatgen.keys())
#print(len(temp_keys))
#for i in range(len(temp_keys)):
#	print(temp_keys[i],entries_for_pymatgen[temp_keys[i]]['energy_per_norm_FU'])
#
#sys.exit()
#print ("\nNow moving on to reference calculations...\n")
#sys.exit()
# now need to go get ICSD structures if they exist - should be similar to the above

if icsd != '':
	known_materials={}
	for i in range(len(all_elements)):
		known_materials[all_elements[i]]=[]
	os.chdir(icsd)
	shutil.copy(str("../"+lib),".")
	structs=0
	folders=[]
	files=glob.glob("*")
	for i in range(len(files)):
		if os.path.isdir(files[i]):
			structs+=1
			folders.append(files[i])
	
	#print("folders: ",folders)
	for i in range(len(folders)):
		pymatgen_keys=list(entries_for_pymatgen.keys())
		#print("\n")
		os.chdir(str(folders[i]))
		shutil.copy(str("../"+lib),".")
		cif=glob.glob("*.cif")[0]
		atoms1 = read("POSCAR")
		atoms.rattle(0.1)
		atoms1,energy,converged=run_gulp(atoms=atoms1,shel=shel,kwds=kwds,opts=gulp_opts,lib=lib)
		energy=energy/len(atoms1)
		#label=glob.glob("*.txt")
		#label=label[0][:-4]+".cif"
		#write(label,atoms1)
		#print("atoms1 ",atoms1)
		syms=atoms1.get_chemical_symbols()
		#print("syms: ",syms)
		counts={}
		for j in range(len(all_elements)):
			if syms.count(all_elements[j]) != 0:
				counts[all_elements[j]]=syms.count(all_elements[j])
		#print(counts)
		keys=list(counts.keys())
		values=list(counts.values())
		#print("keys: ",keys)
		#print("values: ",values)
		values_norm=[]
		
		for j in range(len(norm_to)):
			if syms.count(norm_to[j]) !=0:
				values_norm.append(syms.count(norm_to[j]))
		
		norm=min(values_norm)
		
		for j in range(len(values)):
			if values[j] != 0:
				values[j]=values[j]/norm

		for j in range(len(values_norm)):
			if values_norm[j] !=0:
				values_norm[j] = values_norm[j]/norm
				
		# check that the normalisation produced formlulas containing all integer values
		# if not, correct it.

		all_integers=True
		indicies=[]
		
		#print (values)
		
		for j in range(len(values)):
			value=values[j]
			test=float(value).is_integer()
			if test == False:
				all_integers=False
				indicies.append(j)
				
		
		if all_integers == False:
			#print (counts)
			#print(indicies)
			import math
			norms=[]
			for j in range(len(indicies)):
				num=values[indicies[j]]
				frac,whole=math.modf(num)
				if frac != 0:
					#print("frac whole: ", frac,whole)
					#print("1/frac: ", 1/frac)
					for k in range(len(values)):
						values[k]=Decimal(values[k]*(1/frac))
						values[k]=float(values[k].quantize(Decimal('1e-6')))
					all_integers = True
				
				for k in range(len(values)):
					#print("hello")
					test=values[k].is_integer()
					if test == False:
						all_integers = False
				
				while all_integers == False:
					
					for j in range(len(values)):
						value=values[j]
						test=float(value).is_integer()
						if test == False:
							all_integers=False
							indicies.append(j)
					
					num=values[indicies[j]]
					frac,whole=math.modf(num)
					if frac != 0:
						#print("frac whole: ", frac,whole)
						#print("1/frac: ", 1/frac)
						for k in range(len(values)):
							values[k]=Decimal(values[k]*(1/frac))
							values[k]=float(values[k].quantize(Decimal('1e-6')))
						for k in range(len(values)):
							test=values[k].is_integer()
							if test == True:
								all_integers = True
					
					
		#print(all_integers)
		
		for j in range(len(keys)):
			counts[keys[j]]=values[j]
		
		comp_label=''
		for j in range(len(keys)):
			comp_label += str(keys[j])
			comp_label += str(int(values[j]))
		
		norm_label=''
		norm_counts={}
		for j in range(len(norm_to)):
			if norm_to[j] in counts.keys():
				norm_label+=str(norm_to[j])
				norm_label+=str(int(counts[norm_to[j]]))
				norm_counts[norm_to[j]]=counts[norm_to[j]]
		
		norm_values=list(norm_counts.values())
		energy = energy * sum(values)

		for j in range(len(all_elements)):
			if all_elements[j] in list(counts.keys()):
				known_materials[all_elements[j]].append(counts[all_elements[j]])
			else:
				known_materials[all_elements[j]].append(0.0)
				
		
		#print(comp_label)
		
		if comp_label in pymatgen_keys:
			entries_for_pymatgen[comp_label]['energy_per_norm_FU'].append(energy)
			
		if not comp_label in pymatgen_keys:
			entries_for_pymatgen[str(comp_label)]={'energy_per_norm_FU':[energy],'pymat_comp':{norm_label:norm_counts},'comp_dict':counts}
			pymatgen_keys.append(comp_label)	
					
		os.chdir("../")

	os.chdir("../")
#print(entries_for_pymatgen)
#print(len(entries_for_pymatgen.keys()))
#print(known_materials)
## now go on to write out pymatgen files

pymatgen_keys=list(entries_for_pymatgen.keys())
print("Unique compositions for pymatgen: \n",pymatgen_keys)
#sys.exit()
os.chdir("pymatgen")

for i in range(len(pymatgen_keys)):
	#print(pymatgen_keys[i])
	#print(entries_for_pymatgen[pymatgen_keys[i]])
	pymat_label=list(entries_for_pymatgen[pymatgen_keys[i]]['pymat_comp'].keys())[0]
	pymat_energy=min(entries_for_pymatgen[pymatgen_keys[i]]['energy_per_norm_FU'])
	#print(pymat_label)
	#print(pymat_energy)
	f=open(str(i+1)+".py",'w')
	f.write(
	"import pylab\nimport pymatgen\nfrom pymatgen.entries.computed_entries import ComputedEntry\nfrom pymatgen.ext.matproj import MPRester\nfrom pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter\nimport sys\n"	
	)
	
	f.write("entries = []\n")
	
	f.write(str("entries.append(ComputedEntry('"+str(pymat_label)+"',"+str(pymat_energy)+"))\n"))
	
	
	for j in range(len(pymatgen_keys)):
		pymat_label=list(entries_for_pymatgen[pymatgen_keys[j]]['pymat_comp'].keys())[0]
		pymat_energy=min(entries_for_pymatgen[pymatgen_keys[j]]['energy_per_norm_FU'])
		f.write(str("entries.append(ComputedEntry('"+str(pymat_label)+"',"+str(pymat_energy)+"))\n"))
			
	f.write("\n\n")
	
	f.write("pd=PhaseDiagram(entries)\n")
	f.write(str("fi=open('pymat_"+str(i+1)+"_out.txt','w')"))
	f.write('\nimport collections\ndecomp,ehull=pd.get_decomp_and_e_above_hull(entries[0])\nprint ("Decomposition:")\nprint(decomp)\nprint ("Energy above hull:")\nprint(ehull)\nfi.write("Composition"+str(entries[0])+"\\n")\nfi.write(str(pd.get_decomp_and_e_above_hull(entries[0])))\nfi.write(str("\\nEnergy above hull:\\n"))\nfi.write(str(ehull))\n')
	f.close()


# now need to cycle through running the pymatgen files
hull=open("hull.txt",'w')
files_to_run=glob.glob("*.py")
for i in range(len(files_to_run)):
	try:
		os.system(str(system_python_command +" "+ str(i+1) + ".py > out.txt"))
		fi=open(str("pymat_"+str(i+1)+"_out.txt"),'r')
		lines=fi.readlines()
		hull.write(str(str(i+1).rjust(3)+"      "+str(lines[-1]).ljust(10)+"\n"))
	except:
		pass
hull.close()

## now need to assemble the data table required to plot a convex hull 

#print(convex_hull)

hull=open("hull.txt",'r').readlines()
for i in range(len(pymatgen_keys)):
	fi=open(str("pymat_"+str(i+1)+"_out.txt"),'r')
	lines=fi.readlines()
	pymat_energy_raw=float(lines[-1])
	entries_for_pymatgen[pymatgen_keys[i]]['pymat_energy_raw']=pymat_energy_raw
	comp_dict=entries_for_pymatgen[pymatgen_keys[i]]['comp_dict']
	key=list(entries_for_pymatgen[pymatgen_keys[i]]['pymat_comp'])
	norm_dict=entries_for_pymatgen[pymatgen_keys[i]]['pymat_comp'][key[0]]
	
	#print(comp_dict)
	#print(norm_dict)
	
	comp_value=sum(list(comp_dict.values()))
	norm_value=sum(list(norm_dict.values()))
		
	pymat_energy=pymat_energy_raw*norm_value
	pymat_energy=pymat_energy/comp_value
	pymat_energy=pymat_energy*1000
	
	entries_for_pymatgen[pymatgen_keys[i]]['pymat_energy']=pymat_energy
	
	convex_hull['energy from hull (meV/atom)'].append(pymat_energy)
	
	keys=list(comp_dict.keys())
	
	for j in range(len(all_elements)):
		if all_elements[j] in keys:
			convex_hull[all_elements[j]].append(comp_dict[all_elements[j]])
		else:
			convex_hull[all_elements[j]].append(0.0)
	
data=pandas.DataFrame.from_dict(convex_hull)
print(data)
data.to_csv("../hull.csv",index=False)
	
## create a data table containing the compositions found in the icsd folder if used

if icsd != '':
	data2=pandas.DataFrame.from_dict(known_materials)
	data2.to_csv("../known_materials.csv",index=False)
	
	
	
	
	
	
	
	
	
	
	
		
