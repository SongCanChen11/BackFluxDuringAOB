import re
import math
import numpy as np

def parse_one_side(text):
	components={}
	a=re.split(r' \+ ',text.strip())
	for i in a:
		match=re.search(r'^(\d)\s*(.*?)$',i.strip())
		if match:
			stiochiometry=match.group(1)
			component=match.group(2)
			components[component]=float(stiochiometry)
		else:
			components[i.strip()]=1
	return components


def parse_equation(text):
	left,right=re.split(r' == ',text.strip())
	reactants=parse_one_side(left)
	products=parse_one_side(right)
	reaction={'reactants':reactants,'products':products}
	return reaction


def get_all_reactions(infile,dG0_column):
	reactions=[]
	num=0
	for l in open(infile).readlines()[1:]:
		a=re.split(r'\t',l.strip())
		if num<=16:
			reaction=parse_equation(a[1])
			reaction['dG0']=float(a[dG0_column])
			reaction['dGdis']=0#float(a[3])
			reaction['sto']=float(a[3])
			reactions.append(reaction)
		elif num==18:
			dG0_SR=float(a[dG0_column])
		elif num==17:
			dG0_overall=float(a[dG0_column])
		num+=1



	return [reactions,dG0_overall,dG0_SR]


def get_all_substrates(reactions):
	all_substrates=[]
	for r in reactions:
		for item in reactions:
			reactants=item['reactants']
			products=item['products']
			for component in reactants:
				if not component in all_substrates:
					all_substrates.append(component)
				
			for component in products:
				if not component in all_substrates:
					all_substrates.append(component)
				#all_substrates.add(component)
	return all_substrates

def linear_solver(reactions,all_substrates,concentration_constraints):
	A=[]
	y=[]
	index=0
	for r in reactions:
		if index!=reaction_for_global_dG_calculation:
			stiochiometry=[]
			for s in all_substrates:
				if s in r['products']:
					stiochiometry.append(float(r['products'][s]))
				elif s in r['reactants']:
					stiochiometry.append(float((-1)*r['reactants'][s]))
				else:
					stiochiometry.append(0)
			A.append(stiochiometry)
			y.append((r['dGdis']-r['dG0'])/R/T)
		index+=1

	for i in concentration_constraints:
		stiochiometry=[]
		for s in all_substrates:
			if s==i:
				stiochiometry.append(1)
			else:
				stiochiometry.append(0)
		A.append(stiochiometry)
		y.append(float(concentration_constraints[i]))

	
	#print(np.array(A,dtype=float).shape)
	#print(np.array(y,dtype=float).shape)
	#print(A)
	#print(y)
	x = np.linalg.solve(A, y)
	return x
	#print(x)


def update(all_substrates,current_logC):
	total_C={}
	for i in maxC_constraints:
		total_C[i]=0
		for s in all_substrates:
			if re.search(r'%s'%(re.sub(r'HS-',r'',i)),s):
				#print(i,s)
				indice=all_substrates.index(s)
				total_C[i]+=math.exp(current_logC[indice])

	new_logC_coenzymes={}
	current_logC_coenzymes={}
	for s in all_substrates:
		if s in logC_coenzymes:
			indice=all_substrates.index(s)
			current_logC_coenzymes[s]=current_logC[indice]
			try:
				new_logC_coenzymes[s]=math.log(math.exp(current_logC[indice])-math.exp(current_logC[indice])*(total_C[s]-maxC_constraints[s])/total_C[s])
			except ValueError:
				new_logC_coenzymes[s]=-100

	#print(new_logC_coenzymes)
	#print(current_logC_coenzymes)
	#print('\n')
	return new_logC_coenzymes

def calc_global_dG(current_logC):
	r=reactions[reaction_for_global_dG_calculation]
	dG=r['dG0']
	for i in r['reactants']:
		indice=all_substrates.index(i)
		dG-=R*T*current_logC[indice]*r['reactants'][i]
	
	for i in r['products']:
		indice=all_substrates.index(i)
		dG+=R*T*current_logC[indice]*r['products'][i]

	return dG

def check_dG(reactions,current_logC):
	global_dG=0
	print('\n\nCheck dG of the reactions:\n\n')
	for r in reactions:
		dG=r['dG0']
		for i in r['reactants']:
			indice=all_substrates.index(i)
			dG-=R*T*current_logC[indice]*r['reactants'][i]
		
		for i in r['products']:
			indice=all_substrates.index(i)
			dG+=R*T*current_logC[indice]*r['products'][i]
		global_dG+=dG*r['sto']

		print(' + '.join(r['reactants']) + ' == ' + ' + '.join(r['products'])+'\ndG=%.2f\n'%(dG))
	print('global_dG: %.1f'%(global_dG))
	print('catalytic energy: %.3f'%(calc_dGcat_Ex(logC_substrates)))

#def calc_dGcat(logC_substrates):
	# C4H10 (g) + 8 H2O → 4 CO2 (g) + 13 H2 (g) dG0=337.2
#	dG=337.2
#	dG+=R*T*(13*logC_substrates['H2 (g)']+4*logC_substrates['CO2 (g)']-logC_substrates['C4H10 (g)'])
#	return dG


def calc_dGcat_Ex(logC_substrates):
	# C4H10 (g) + 8 H2O → 4 CO2 (g) + 13 H2 (g) dG0=337.2
	#dG=-199.9
	dG=dG0_overall#-249.88
	#dG=-250.08
	dG+=R*T*(13*logC_substrates['XH2']+4*logC_substrates['CO2 (g)']-logC_substrates['C4H10 (g)'])
	return dG

def calc_dG_SR(logC_substrates):
	#4 XH2 (g) + SO42- + H+ == HS- + 4H2O + 4 X
	#dG=28.198
	#dG=-152.45
	dG=dG0_SR
	sulfide=0.005
	sulfate=0.028
	dG+=R*T*(math.log(sulfide)-math.log(sulfate)-4*logC_substrates['XH2'])
	return dG



def calc_limited_conc(concentration_constraints):
	for cycle in range(30):
		current_logC=linear_solver(reactions,all_substrates,concentration_constraints)
		logC_coenzymes=update(all_substrates,current_logC)
		concentration_constraints={**logC_substrates,**logC_constant,**logC_coenzymes}

	return current_logC
		
def summarize_energy_coupling_sites(reactions):
	count=0
	for r in reactions:
		if r['dGdis']>0:
			print('Energy re-invested at site %s\t%.3f'%(count+1,r['dGdis']))
		if r['dGdis']<0:
			print('Energy dissipation at site %s\t%.3f'%(count+1,r['dGdis']))
		count+=1


def calc_dG0_overall(reactions):
	dG0=0
	for r in reactions:
		dG0+=r['dG0']*r['sto']
	return dG0

def check_reactions(reactions):
	def add(x,y,sto):
		z={}
		for i in list(x.keys()) + list(y.keys()):
			z[i] = x.get(i,0) + y.get(i,0)*sto
		return z

	def substract(x,y):
		z={}
		for i in list(x.keys()) + list(y.keys()):
			z[i] = x.get(i,0) - y.get(i,0)
		return z

	def print_reaction(d):
		reactants = []
		products = []
		for i in d:
			if d[i]>0:
				reactants.append(str(d[i])+i)
			elif d[i] <0:
				products.append(str(int(abs(d[i])))+i)
		print(' + '.join(reactants) + ' == ' + ' + '.join(products)+' dG0 = %.2f'%(calc_dG0_overall(reactions)))


	reactants = {}
	products = {}
	#sto=[1]*7+[3,2]+[2]*6+[4,2,2]
	for i in range(len(reactions)):
		r = reactions[i]['reactants']
		p = reactions[i]['products']
		reactants = add(reactants, r, reactions[i]['sto'])
		products = add(products, p, reactions[i]['sto'])
	#print(reactants)
	#print(products)

	print_reaction(substract(reactants,products))



def is_feasible(current_logC):
	report=['C4H9-S-CoM','Butyryl-CoA','Crotonyl-CoA','(S)-3-Hydroxybutyryl-CoA','Acetoacetyl-CoA','Acetyl-CoA','CH3-H4MPT','CH2=H4MPT','CH≡H4MPT+','HCO-H4MPT','HCO-MF']#,'Acetyl-CoA','CO2 (g)']
	for s in report:
		indice=all_substrates.index(s)
		value=current_logC[indice]
		log10C=math.log10(math.exp(value))
		global_dG=calc_global_dG(current_logC)
		if log10C<=-6.5 or global_dG > 0 or log10C>-2:
			return False
	return True




def find_optimal_setting(feasible_settings):
	if len(feasible_settings)==0:
		print('No feasible setting was found!!\n')
		return False
		#exit()
	else:
		minimum=1000000
		#print(feasible_settings)
		for s in feasible_settings:
			dG_set=s['dG_set']
			total_energy_dissipated=0
			for step in dG_set:
				name=step['name']
				index=step['reaction index']
				dG=step['dG']
				sto=reactions[index]['sto']
				if re.search('dissipation',name):
					total_energy_dissipated-=1.0*float(sto)*float(dG)


			if minimum>total_energy_dissipated:
				minimum=total_energy_dissipated
				optimal_setting=s

		return optimal_setting

		#current_logC=optimal_setting['current_logC']
		#dG_set=optimal_setting['dG_set']
		#print_log10C(current_logC,all_substrates)
		#print_dG_set(dG_set)
		#print('energy global energy\t%.3f'%(calc_global_dG(current_logC)))
		#print('total energy dissipated: %.3f'%(minimum))
		#print('energy\tremained for SR\t%.3f\n'%(calc_dG_SR(logC_substrates)))


def find_optimal_none_feasible_setting(none_feasible_settings):
	if len(none_feasible_settings)==0:
		print('No non-feasible setting was found!!\n')
		return False
		#exit()
	else:
		minimum=1000000
		#print(feasible_settings)
		for s in none_feasible_settings:
			current_logC=s['current_logC']
			dist=0
			report=['C4H9-S-CoM','Butyryl-CoA','Crotonyl-CoA','(S)-3-Hydroxybutyryl-CoA','Acetoacetyl-CoA','Acetyl-CoA','CH3-H4MPT','CH2=H4MPT','CH≡H4MPT+','HCO-H4MPT','HCO-MF']
			for metabolite in report:
				indice=all_substrates.index(metabolite)
				value=current_logC[indice]
				if value >=-6 and value <= -2:
					dist+=0
				elif value<-6:
					dist+=abs(value+6)
				elif value>-2:
					dist+=abs(value+2)

			if minimum>dist:
				minimum=dist
				optimal_setting=s

		return optimal_setting




def search_for_feasible_dG_set(dG,step,dG_set):
	if step>=len(dG):
		for j in dG_set:
			index=j['reaction index']
			reactions[index]['dGdis']=j['dG']
		concentration_constraints={**logC_substrates,**logC_constant,**logC_coenzymes}
		current_logC=calc_limited_conc(concentration_constraints)
		if is_feasible(current_logC):

			feasible_settings.append({'current_logC':current_logC.copy(),'dG_set':dG_set.copy()})
			#print_log10C(current_logC,all_substrates)
			#print_dG_set(dG_set)
			#print('energy\tglobal energy\t%.3f'%(calc_global_dG(current_logC)))
			#print('total energy dissipated: %.3f'%(dG_set[]))
			#print('energy\tremained for SR\t%.3f\n\n\n'%(calc_dG_SR(logC_substrates)))
			#exit()
		else:
			global_dG=calc_global_dG(current_logC)
			if global_dG<0:
				none_feasible_settings.append({'current_logC':current_logC.copy(),'dG_set':dG_set.copy()})
	else:
		for i in dG[step]['dG_range']:
			dG_set.append({'reaction index':dG[step]['reaction index'],'dG':i,'name':dG[step]['name']})
			search_for_feasible_dG_set(dG,step+1,dG_set)
			dG_set.pop()


def calc_net_H_quantum(dG_set,global_dG,global_dG_sto):
	n = int(global_dG * global_dG_sto / (-20))
	for i in dG_set:
		if re.search('proton gradient',i['name']):
			index=i['reaction index']
			n+=int(i['dG']*reactions[index]['sto']/(-20))
	return n

def calc_total_energy_dissipated(dG_set):
	dG_dis = 0
	for i in dG_set:
		if re.search('dissipation',i['name']):
			index=i['reaction index']
			dG_dis+=i['dG']*reactions[index]['sto']
	return dG_dis

def write_output(conc_XH2,optimal_setting,fout1,fout2):
	current_logC=optimal_setting['current_logC']
	dG_set=optimal_setting['dG_set']
	dG_overall=calc_dGcat_Ex(logC_substrates)
	global_dG = calc_global_dG(current_logC)
	global_dG_sto = reactions[reaction_for_global_dG_calculation]['sto']
	net_H_quantum=calc_net_H_quantum(dG_set,global_dG,global_dG_sto)
	dG_dis=calc_total_energy_dissipated(dG_set)
		
	for i in range(len(substrates_to_report)):
		s=substrates_to_report[i]
		indice=all_substrates.index(s)
		value=current_logC[indice]
		
		fout1.write('%s\t%.1e\t%s\t%s\t%.3f\n'%(EX,conc_XH2,i,s,math.log10(math.exp(value))))

	for i in dG_set:
		index=i['reaction index']
		fout2.write('%s\t%.1e\t%.2f\t%.2f\t%i\t%s\t%.2f\t%s\n'%(EX,conc_XH2,dG_overall,dG_dis,net_H_quantum,i['name'],i['dG'],reactions[index]['sto']))

	fout2.write('%s\t%.1e\t%.2f\t%.2f\t%i\t%s\t%.2f\t%s\n'%(EX,conc_XH2,dG_overall,dG_dis,net_H_quantum,name_of_reaction_for_global_dG_calc,global_dG,global_dG_sto))



## Main function
#logC_constant={'H+':0,'H2O':0,'F420':0,'Fdox':0,'FAD':0,'NAD+':0,'X':0}
R=0.008314
T=298
# CoM to CoA conversion: 2 (index); 4 (line number in the table)
dG_CoM2CoA=[40,60,80,100,120,140]
# NADH oxidation: 7 (index); 9 (line number in the table)
dG_NADH_ox=[0,-5,-10,-15]
# F420H2 oxidation: 14 (index); 16 (line number in the table)
dG_F420H2_ox=[0,-10,-20,-30]
# MQH2 oxidation: 16 (index): 18 (line number in the table)
dG_MQH2_ox=[0,10,20,30]
# Fdred oxidation: 15 (index); 17 (line number in the table)
dG_Fd_ox=[0,-20,-40,-60,-80,-100]
# methylene-H4MPT oxidation: 10 (index); 12 (line number in the table)
dG_CH2_ox=[0,-5,-10,-15]


dG=[{'name':'CoM2CoA (consumption of proton gradient)','dG_range':dG_CoM2CoA,'reaction index':2},
	{'name':'NADH_ox (energy dissipation)','dG_range':dG_NADH_ox,'reaction index':7},
	#{'name':'F420H2_ox','dG_range':dG_F420H2_ox,'reaction index':14},
	{'name':'MQH2_ox (consumption of proton gradient)','dG_range':dG_MQH2_ox,'reaction index':16},
	{'name':'Fd_ox (build-up of proton gradient)','dG_range':dG_Fd_ox,'reaction index':15},
	{'name':'CH2_ox (energy dissipation)','dG_range':dG_CH2_ox,'reaction index':10}]
logC_constant={'H+':0,'H2O':0,'F420':0,'Fdox':0,'NAD+':0,'X':0,'MQ':0}
logC_coenzymes={'HS-CoB':0,'HS-CoM':0,'CoA':0,'H4MPT':0,'MF':0}
maxC_constraints={'CoA':0.01,'HS-CoB':0.01,'H4MPT':0.01,'HS-CoM':0.01,'MF':0.01}
substrates_to_report=['C4H10 (g)','C4H9-S-CoM','Butyryl-CoA','Crotonyl-CoA','(S)-3-Hydroxybutyryl-CoA','Acetoacetyl-CoA','Acetyl-CoA','CH3-H4MPT','CH2=H4MPT','CH≡H4MPT+','HCO-H4MPT','HCO-MF','CO2 (g)','F420H2','NADH','Fdred']



model = 'AOB_M1'
infile='17_AOB_M1.txt'
EX = 'Ex = -0.220V'
dG0_column = 6
conc_XH2_list=[2,4,8]
name_of_reaction_for_global_dG_calc='F420H2_ox (build-up of proton gradient)'
reaction_for_global_dG_calculation=14
output_prefix='17_AOB_M1_SRB__res'





output_log10C = output_prefix+'_log10C.txt'
output_dG_set = output_prefix+'_dG.txt'
non_feasible_log10C = output_prefix+'_log10C_nonfeasible.txt'
non_feasible_dG_set = output_prefix+'_dG_nonfeasible.txt'
fout1=open(output_log10C,'w')
fout2=open(output_dG_set,'w')
fout3=open(non_feasible_log10C,'w')
fout4=open(non_feasible_dG_set,'w')
reactions,dG0_overall,dG0_SR=get_all_reactions(infile,dG0_column)
all_substrates=get_all_substrates(reactions)



for conc_XH2 in conc_XH2_list:#[50,100,500,1000]:
	logC_substrates={'C4H10 (g)':math.log(0.3),'CO2 (g)':math.log(0.1),'XH2':math.log(conc_XH2)}
	check_reactions(reactions)
	dG0_overall=calc_dG0_overall(reactions)
	dG_overall=calc_dGcat_Ex(logC_substrates)
	print('Actual gibbs free energy of AOB: %.2f'%(dG_overall))
	print('Gibbs free energy of SR: %.2f'%(calc_dG_SR(logC_substrates)))
	print('XH2 concentration: %.1e\n'%(conc_XH2))
	feasible_settings=[]
	none_feasible_settings=[]
	search_for_feasible_dG_set(dG,0,[])
	optimal_setting=find_optimal_setting(feasible_settings)
	optimal_none_feasible_setting=find_optimal_none_feasible_setting(none_feasible_settings)

	if optimal_setting:
		write_output(conc_XH2,optimal_setting,fout1,fout2)
	else:
		write_output(conc_XH2,optimal_none_feasible_setting,fout3,fout4)






