#!/usr/bin/python3

from loaddata import * 
from filternplot import *

import sys
args = sys.argv
alldata, cndf = LoadTables()


if args[1]=="jku":
	kth = 21.5
	jku(alldata, kth)
elif args[1]=="stern2005_1x1":
	stern2005_1x1_alluse(alldata, kth=21.0)
elif args[1]=="stern2005_1x1":
	stern2005_1x1(alldata)
elif args[1]=="stern2005_3x2":
	stern2005_3x2(alldata)
elif args[1]=="lacy2007_1x1":
	s = lacy2007_1x1(alldata, kth=21.0)
	#uvj_3x2(kth=21.0,selection=s)
	#plot_fracs(kth=21.0,selection=s)
elif args[1]=="color_color_1x1":
	s = color_color_1x1(alldata, kth=21.0)
elif args[1]=="color_color_3x2":
	s = color_color_3x2(alldata, kth=21.0)
elif args[1]=="uvj":
	kth = 23.4
	uvj_3x2(alldata, kth)
elif args[1]=="plot_fracs":
	plot_fracs(alldata)
elif args[1]=="smf":
	kth = 23.4
	s = mass_selection(alldata)
	#s2 = lacy2007_1x1(kth=kth,selection=s)
	#print(len(s2))
	uvj_3x2(alldata, selection=s)
	plot_fracs(alldata, selection=s)
elif args[1]=="xray":
	#s_cm = cm('all')
	#exit(1)
	s_cm = crossmatch(alldata, cat='all',cclass22=[])
	#stern2005_3x2(alldata, selection=s_cm)
	#labelx=r'$log(J/H)$'
	#labely=r'$log(K_S/f_{3.6\mu m})$'
	#labelx=None
	#labely=None
	labelx=r'$log(f_{5.8\mu m}/f_{3.6\mu m})$'
	labely=r'$log(f_{8.0\mu m}/f_{4.5\mu m})$'
	s2 = lacy2007_1x1(alldata, selection=s_cm)
	s2 = color_color_3x2(alldata, selection=s_cm, labelx=labelx, labely=labely)
	exit(1)
	s_cm = crossmatch(alldata, cat='xmm')
	s2 = lacy2007_1x1(alldata, selection=[s_cm[0],])
	s_cm = crossmatch(alldata, cat='chandra',cclass18=[1,])
	s2 = lacy2007_1x1(alldata, selection=[s_cm[0],])
	s_cm = crossmatch(alldata, cat='chandra',cclass18=[2,])
	s2 = lacy2007_1x1(alldata, selection=[s_cm[0],])
elif args[1]=="split":
	# First we obtain the mass cut and returns selection list 
	# differenciated in (color_)IRAGN and (color_)not_IRAGN
	s_smf = mass_selection(alldata, cndf)
	# From here we have to separate the not_IRAGN to feed it
	# into the crossmatch:
	print('len(s_smf)=%d' % len(s_smf[0]))
	s_smf_niragn = [
			[a for (a,b) in zip(s_smf[0],s_smf[1]) if b == '[not IR-AGN]'],
			[b for (a,b) in zip(s_smf[0],s_smf[1]) if b == '[not IR-AGN]']
			]
	print('len(s_smf_niragn)=%d' % len(s_smf_niragn[0]))
	s_smf_iragn = [
			[a for (a,b) in zip(s_smf[0],s_smf[1]) if b == '[IR-AGN]'],
			[b for (a,b) in zip(s_smf[0],s_smf[1]) if b == '[IR-AGN]']
			]
	print('len(s_smf_iragn)=%d' % len(s_smf_iragn[0]))
	# The result is fed into the crossmatch
	# with X-ray sources which returns those matching (color_)XAGN
	# and (color_)not_AGN.
	# Several catalog to choice from here:
	#s_cm = crossmatch(cat='all',selection=s_smf_niragn)
	print(Counter(s_smf[1]))
	s_cm = crossmatch(alldata, cat='all',cclass22=[],selection=s_smf)
	s_cm_xagn = [
			[a for (a,b) in zip(s_cm[0],s_cm[1]) if b in mycolors_xagn],
			[b for (a,b) in zip(s_cm[0],s_cm[1]) if b in mycolors_xagn]
			]
	print('len(s_cm_xagn)=%d' % len(s_cm_xagn[0]))
	s_cm_iragn = [
			[a for (a,b) in zip(s_cm[0],s_cm[1]) if b == '[IR-AGN]'],
			[a for (a,b) in zip(s_cm[0],s_cm[1]) if b == '[IR-AGN]']
			]
	print('len(s_cm_iragn)=%d' % len(s_cm_iragn[0]))
	s_cm_not_agn = [
			[a for (a,b) in zip(s_cm[0],s_cm[1]) if b == '[not IR-AGN]'],
			[a for (a,b) in zip(s_cm[0],s_cm[1]) if b == '[not IR-AGN]']
			]
	print('len(s_cm_not_agn)=%d' % len(s_cm_not_agn[0]))
	#s_tot = [list(s_smf_iragn[0]) + list(s_cm[0]), list(s_smf_iragn[1]) + list(s_cm[1])]
	s_tot = [list(s_cm[0]), list(s_cm[1])]
	labelx=r'$log(f_{5.8\mu m}/f_{3.6\mu m})$'
	labely=r'$log(f_{8.0\mu m}/f_{4.5\mu m})$'
	s2 = color_color_3x2(alldata, selection=s_tot, labelx=labelx, labely=labely)
	plot_fracs(alldata, selection=s_tot)
	def xray(c):
		if c in mycolors_xagn: return '[Xray AGN]'
		else: return c
	s_cmn = [xray(k) for k in s_cm[1]]
	s_tot = [list(s_cm[0]), list(s_cmn)]
	uvj_3x2(alldata, selection=s_tot)

