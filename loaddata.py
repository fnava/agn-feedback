#!/usr/bin/python3

kth = 21.0
photo_keys=["id","ra","dec","xpix","ypix","Ks_tot","eKs_tot","Ks","eKs","H","eH","J","eJ","Y","eY","ch4","ech4","ch3","ech3","ch2","ech2","ch1","ech1","zp","ezp","ip","eip","rp","erp","V","eV","gp","egp","B","eB","u","eu","IA484","eIA484","IA527","eIA527","IA624","eIA624","IA679","eIA679","IA738","eIA738","IA767","eIA767","IB427","eIB427","IB464","eIB464","IB505","eIB505","IB574","eIB574","IB709","eIB709","IB827","eIB827","fuv","efuv","nuv","enuv","mips24","emips24","K_flag","K_star","K_Kron","apcor","z_spec","z_spec_cc","z_spec_id","star","contamination","nan_contam","orig_cat_id","orig_cat_field","USE"]
redshift_keys=["id","z_spec","z_a","z_m1","chi_a","z_p","chi_p","z_m2","odds","l68","u68","l95","u95","l99","u99","nfilt","q_z","z_peak","peak_prob","z_mc"]
restframe_uv=["id","z","DM","nfilt_fit","chi2_fit","L153","L155"]
restframe_vj=["id","z","DM","nfilt_fit","chi2_fit","L155","L161"]
masses=["id","z","ltau","metal","lage","Av","lmass","lsfr","lssfr","la2t","chi2"]
#import libraries
import numpy as np
import pandas as pd
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.interpolate import interp1d
from astropy.table import Table, join
from astropy.io import ascii
from collections import Counter

mycolors = {
		'[AGN]'				: ["black",		"AGN",				's',"black"],
		'[Xray AGN]'		: ["crimson",	"X-Ray-AGN",		'o',"crimson"],
		'[XMM]'				: ["magenta",	"XMM",				'^',"magenta"],
		#'[XMM]'				: ["green",	"XMM",				'^',"green"],
		'[Chandra Star]'	: ["grey",		"Chandra (Star)",	'o',"grey"],
		'[Chandra BLAGN]'	: ["blue",		"Chandra (BLAGN)",	's',"blue"],
		#'[Chandra not-BLAGN]': ["orangered","Chandra (not-BLAGN)",'s',"orangered"],
		'[Chandra not-BLAGN]': ["orange","Chandra (not-BLAGN)",'s',"orange"],
		'[Chandra unobs]'	: ["blue",		"Chandra (Unobscured)",	's',"blue"],
		'[Chandra obsc]'	: ["orangered",	"Chandra (Obscured)",'s',"orangered"],
		'[Chandra galaxy]'	: ["cyan",		"Chandra (Galaxy)"	,'o',"cyan"],
		'[IR-AGN]'			: ["green",		"IR-AGN",			'x',"green"],
		#'[IR-AGN]'			: ["gold",		"IR-AGN", 'x',"gold"],
		'[not IR-AGN]'		: ["maroon",	"Not AGN",			'.',"none"],
		'[exc IR-AGN]'		: ["turquoise",	"Excluded IR-AGN",	'+',"none"],
		'[not AGN]'			: ["mediumblue","Not AGN",			'.',"none"]
		}
mycolors_xagn = ['[XMM]', '[Chandra BLAGN]', '[Chandra not-BLAGN]', '[Chandra obsc]', '[Chandra unobs]']
mycolors_agn = mycolors_xagn + ['[IR-AGN]',]

def WriteTable():
	# This function is supossed to be used only once to generate the 
	# single FITS compressed file to be loaded. It should weight less than 100MB
	# so it can be uploaded to github.
	print("WriteTable")
	import zipfile
	with zipfile.ZipFile('data/UVISTA_final_v4.1.cat2.zip') as myzip:
		with myzip.open('UVISTA_final_v4.1.cat2') as f1:
				d1 = f1.read().decode() 
	#f1 = open('./data/UVISTA_final_v4.1.cat2', 'r')
	#d1 = f1.read()
	f1.close()
	photo = ascii.read(d1, guess=False, format='basic')
	f2 = open('./data/UVISTA_final_v4.1b.zout', 'r')
	d2 = f2.read()
	redshift = ascii.read(d2, guess=False, format='basic')
	redshift.rename_column('z_spec', 'z_spec2')	
	f2.close()
	f3 = open('./data/UVISTA_final_v4.1.153-155b.rf', 'r')
	d3 = f3.read()
	restframe_uv = ascii.read(d3, guess=False, format='basic')
	f3.close()
	f4 = open('./data/UVISTA_final_v4.1.155-161b.rf', 'r')
	d4 = f4.read()
	restframe_vj = ascii.read(d4, guess=False, format='basic')
	f4.close()
	f5 = open('./data/UVISTA_final_BC03_v4.1b.fout', 'r')
	#f5 = open('../K-selected FAST Stellar Masses and SPS parameters/UVISTA_final_M05_v4.1b.fout', 'r')
	d5 = f5.read()
	masses = ascii.read(d5, guess=False, format='basic')
	masses.rename_column('z', 'z_mass')	
	f5.close()

	join1 = join(photo, redshift, keys="id")
	join2 = join(join1, restframe_uv, keys="id")
	join3 = join(join2, restframe_vj, keys="id")
	alldata = join(join3, masses, keys="id")
	alldata.write('./data/alldata.fits')

def LoadTables():
	print("LoadTables")
	f6 = open('./data/cndf_in.txt', 'r')
	d6 = f6.read()
	cndf = ascii.read(d6, guess=False, format='basic')
	alldata = Table.read('./data/alldata.fits.gz')

	return alldata, cndf

def FilterByFields (alldata, *fields):
	result = []
	for key in fields:
		result.append( alldata[key] )
	return(list(fields), list(zip(*result)))

def FilterByConds (fields, data, funcs, auxcols = None):
	table = []
	for datarow in data:
		condrow = []
		for k,f in funcs.items():
			if auxcols is None:
				condrow.append( f(datarow[fields.index(k)]) )
			else:
				inputargs = [datarow[fields.index(k)], ]
				for a in auxcols:
					inputargs.append(datarow[fields.index(a)])
				condrow.append( f(*inputargs) )
		table.append([datarow,condrow])
	result = [ row[0] for row in table if all(i for i in row[1]) ]
	return result

# Returns a table with a row for each sample and a column for each filter
# fields: list of strings identifying columns in the same order
# data: input values table
# funs: dictonary of functions to be used to transform the respective column (filter)
# auxcols: list of strings identifying columns (filters) to be provided as additional
#       input to the transform function
def FilterByFuncs (fields, data, funcs, auxcols = None):
	table = []
	for datarow in data:
		newrow = []
		for k in fields:
			idx = fields.index(k)
			oldval = datarow[idx]
			if k in funcs.keys():
				if auxcols is None:
					newval = funcs[k](oldval)
				else:
					inputargs = [oldval, ]
					for a in auxcols:
						inputargs.append(datarow[fields.index(a)])	
					newval = funcs[k](*inputargs)
			else:
				newval = oldval
			newrow.append(newval) 
		table.append(newrow)
	return table

# conds
def pos(v): return v>0
def is1(v): return v==1
def is0(v): return v==0
def zcut(z): 
	return z<3.0
def kscut(ks, th=None): 
	if th is None:
		return True
	elif not ks>0: 
		return False
	else:
		return 25-2.5*np.log10(ks)<th
def sncut(sn): return sn>5
import random
def kscut2(flux, ks, ks_tot, th=21.5): 
	#if not ks>0: return False
	return random.random()<0.3 and 25.0-2.5*np.log10(flux*ks_tot/ks)<th
#	else: return random.random()<0.2
# funcs
def colorjku(i):
	if i>0.90: j=1
	else: j=0
	return {0:'red', 1:'blue'}[j]
def color2(i):
	if i==1: return 'red'
	else: return 'blue'
def magAB(flux, ks, ks_tot):
	return 25.0-2.5*np.log10(flux*ks_tot/ks)
	#return 25.0-2.5*np.log10(flux)
def log(flux, ks, ks_tot):
	#return np.log10(flux*ks_tot/ks)
	return np.log10(flux)
def sn1(flux, ks, ks_tot, c1, e1, c2, e2, c3, e3, c4, e4):
	return abs(c1/e1)
def sn2(flux, ks, ks_tot, c1, e1, c2, e2, c3, e3, c4, e4):
	return abs(c2/e2)
def sn3(flux, ks, ks_tot, c1, e1, c2, e2, c3, e3, c4, e4):
	return abs(c3/e3)
def sn4(flux, ks, ks_tot, c1, e1, c2, e2, c3, e3, c4, e4):
	return abs(c4/e4)

def color(i):
	return {0:'[not IR-AGN]', 1:'[IR-AGN]'}[i]

def lacy2007_func(flux, c1, c2, c3, c4, z):
	# z is not used, but included to allow the bariadj variation below
	# and avoid touching too much the call to FilterByFunc
	x = np.log10(c3/c1) #log_58 - log_36
	y = np.log10(c4/c2) #log_80 - log_45
	if (x > -0.1 and 
			y > -0.2 and
			y <= 0.8*x + 0.5): 
		return color(1) 
	else:
		return color(0)

def lacy2007_cond(flux, c1, c2, c3, c4):
	x = c3 - c1 #log_58 - log_36
	y = c4 - c2 #log_80 - log_45
	if (x > -0.1 and 
			y > -0.2 and
			y <= 0.8*x + 0.5 and
			c2 > c1 and c3 > c2 and c4 > c3):
		return True 
	else:
		return False

def donley2012_func(flux, c1, c2, c3, c4, z):
	# z is not used, but included to allow the bariadj variation below
	# and avoid touching too much the call to FilterByFunc
	#x = np.log10(c3/c1) #log_58 - log_36
	#y = np.log10(c4/c2) #log_80 - log_45
	x = c3 - c1
	y = c4 - c2
	if (x >= 0.08 and 
			y >= 0.15 and
			y <= 1.21*x + 0.27 and
			y >= 1.21*x - 0.27):
		return color(1) 
	else:
		return color(0)

def donley2012_cond(flux, c1, c2, c3, c4):

	x = c3 - c1 #log_58 - log_36
	y = c4 - c2 #log_80 - log_45
	if (x >= 0.08 and 
			y >= 0.15 and
			y <= 1.21*x + 0.27 and
			y >= 1.21*x - 0.27 and
			c2 > c1 and c3 > c2 and c4 > c3):
		return True 
	else:
		return False
	
bariadj_obsagn=[
		(0.0, 0.2, 0.0, 0.0),
		(0.2, 0.5, -0.12334819085155703, 0.09894967234429784),
		(0.5, 1.0, -0.13185353031352148, -0.0671463347445458), 
		(1.0, 1.5, -0.07989849010137863, -0.07561034249145582), 
		(1.5, 2.0, 0.10125731448524965, -0.04382036450512951), 
		(2.0, 2.5, 0.23099877403033042, 0.08357985231782124), 
		(2.5, 3.0, 0.26954172677464505, 0.16298817160981743),
		(3.0, 10.0, 0.0, 0.0)	
		]
bariadj_obsagn_ref = (-0.12334819085155703, 0.09894967234429784)
bariadj_obsagn_ref1 = (-0.07989849010137863, -0.07561034249145582) 


bariadj_blagn=[
		(0.0, 0.2, 0.0, 0.0),
		(0.2, 0.5, -0.12334819085155703, 0.09894967234429784), 
		(0.5, 1.0, -0.13185353031352148, -0.0671463347445458), 
		(1.0, 1.5, -0.07989849010137863, -0.07561034249145582), 
		(1.5, 2.0, 0.10125731448524965, -0.04382036450512951), 
		(2.0, 2.5, 0.23099877403033042, 0.08357985231782124), 
		(2.5, 3.0, 0.26954172677464505, 0.16298817160981743),
		(3.0, 10.0, 0.0, 0.0)
		]
bariadj_blagn_ref = (-0.13185353031352148, -0.0671463347445458)

bariadj_all=[
		(0.0, 0.2, 0.0, 0.0),
		(0.2, 0.5, -0.12334819085155703, 0.09894967234429784), 
		(0.5, 1.0, -0.13185353031352148, -0.0671463347445458), 
		(1.0, 1.5, -0.07989849010137863, -0.07561034249145582), 
		(1.5, 2.0, 0.10125731448524965, -0.04382036450512951), 
		(2.0, 2.5, 0.23099877403033042, 0.08357985231782124), 
		(2.5, 3.0, 0.26954172677464505, 0.16298817160981743),
		(3.0, 10.0, 0.0, 0.0)
		]
bariadj_all_ref1 = (-0.07989849010137863, -0.07561034249145582) 

bariadj = bariadj_all
bariadj_ref = bariadj_all_ref1
#bariadj = bariadj_blagn
#bariadj_ref = bariadj_blagn_ref

def donley2012_baricut_func(flux, c1, c2, c3, c4, z):
	adj = [(bariadj_ref[0]-x, bariadj_ref[1]-y) 
		for (z0, z1, x, y) 
		in bariadj if z > z0 and z <= z1]
	#x = np.log10(c3/c1) #log_58 - log_36
	#y = np.log10(c4/c2) #log_80 - log_45
	x = c3 - c1
	y = c4 - c2
	if z >= 1.0:
		oy = -adj[0][1]
	else:
		oy = 0.0
	if (x >= 0.08 and 
			y >= 0.15+oy and
			y <= 1.21*x + 0.27 and
			y >= 1.21*x - 0.27 and
			c2 > c1 and c3 > c2 and c4 > c3):
		return '[IR-AGN]'
	elif (x >= 0.08 and 
			y < 0.15+oy and
			y >= 0.15 and
			y <= 1.21*x + 0.27 and
			y >= 1.21*x - 0.27 and
			c2 > c1 and c3 > c2 and c4 > c3):
		return '[exc IR-AGN]'
	else:
		return '[not IR-AGN]'

def donley2012_bariadj_func(flux, c1, c2, c3, c4, z):
	adj = [(bariadj_ref[0]-x, bariadj_ref[1]-y) 
		for (z0, z1, x, y) 
		in bariadj if z > z0 and z <= z1]
	x = np.log10(c3/c1)+adj[0][0] #log_58 - log_36
	y = np.log10(c4/c2)+adj[0][1] #log_80 - log_45
	if (x >= 0.08 and 
			y >= 0.15 and
			y <= 1.21*x + 0.27 and
			y >= 1.21*x - 0.27 and
			c2 > c1 and c3 > c2 and c4 > c3):
		return '[IR-AGN]'
	else:
		return '[not IR-AGN]'

def uplogb(flux, ks, ks_tot):
	#return 25-2.5*np.log10(flux*ks_tot/ks)
	return 25-2.5*np.log10(flux)
def uplog1b(flux, ks, ks_tot, c1, e1, c2, e2, c3, e3, c4, e4):
	return np.log10(c1*ks_tot/ks) if c1>2*e1 else np.log10(2*e1*ks_tot/ks)
def uplog2b(flux, ks, ks_tot, c1, e1, c2, e2, c3, e3, c4, e4):
	return np.log10(c2*ks_tot/ks) if c2>2*e2 else np.log10(2*e2*ks_tot/ks)
def uplog3b(flux, ks, ks_tot, c1, e1, c2, e2, c3, e3, c4, e4):
	return np.log10(c3*ks_tot/ks) if c3>2*e3 else np.log10(2*e3*ks_tot/ks)
def uplog4b(flux, ks, ks_tot, c1, e1, c2, e2, c3, e3, c4, e4):
	return np.log10(c4*ks_tot/ks) if c4>2*e4 else np.log10(2*e4*ks_tot/ks)

def uplog1(flux, ks, ks_tot, c1, e1, c2, e2, c3, e3, c4, e4):
	return np.log10(c1) if c1>2*e1 else np.log10(2*e1)
def uplog2(flux, ks, ks_tot, c1, e1, c2, e2, c3, e3, c4, e4):
	return np.log10(c2) if c2>2*e2 else np.log10(2*e2)
def uplog3(flux, ks, ks_tot, c1, e1, c2, e2, c3, e3, c4, e4):
	return np.log10(c3) if c3>2*e3 else np.log10(2*e3)
def uplog4(flux, ks, ks_tot, c1, e1, c2, e2, c3, e3, c4, e4):
	return np.log10(c4) if c4>2*e4 else np.log10(2*e4)

def upsn1(flux, ks, ks_tot, c1, e1, c2, e2, c3, e3, c4, e4):
	return c1>2*e1
def upsn2(flux, ks, ks_tot, c1, e1, c2, e2, c3, e3, c4, e4):
	return c2>2*e2
def upsn3(flux, ks, ks_tot, c1, e1, c2, e2, c3, e3, c4, e4):
	return c3>2*e3
def upsn4(flux, ks, ks_tot, c1, e1, c2, e2, c3, e3, c4, e4):
	return c4>2*e4

def upsncut(a, e1, e2, e3, e4):
	if a == color(0): return False
	elif (not e1 and not e3): return False
	elif (not e2 and not e4): return False
	else: return True
	
def upsarrow(a, e1, e2, e3, e4):
	# right, up, left, down
	reply = {
			(False, True,  True,  True ) : 'r',
			(False, False, True,  True ) : 'ur',
			(False, True,  True,  False) : 'dr',
			(True,  True,  False, True ) : 'l',
			(True,  False, False, True ) : 'ul',
			(True,  True,  False, False) : 'dl',
			(True,  False, True,  True ) : 'u',
			(False, False, True,  True ) : 'ur',
			(True,  False, True,  False) : 'ul',
			(True,  True,  True,  False) : 'd',
			(False, True,  True,  False) : 'dr',
			(True,  True,  False, False) : 'dl'
			}
	args = (e1,e2,e3,e4)
	if args in reply:
		return reply[args]
	else:
		return '.'

# https://stackoverflow.com/questions/16861339/down-arrow-symbol-in-matplotlib
arrows = {
		'l':u'$\u2190$',
		'u':u'$\u2191$',
		'r':u'$\u2192$',
		'd':u'$\u2193$',
		'ul':u'$\u2196$',
		'ur':u'$\u2197$',
		'dr':u'$\u2198$',
		'dl':u'$\u2199$',
		'.':'.'
		}

arrows2 = {
		'dr':[(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, -1.0)],
		'ul':[(-1.0, -1.0), (1.0, 1.0), (-1.0, 1.0), (-1.0, -1.0)],
		'dl':[(-1.0, -1.0), (-1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)],
		'ur':[(-1.0, 1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)],
		'u':'^',
		'd':'v',
		'l':'<',
		'r':'>',
		'.':'.'
		}

def rscolor(z):
	if z>2: return "red"
	elif z>1.5: return "orange"
	elif z>1.0: return "yellow"
	elif z>0.5: return "green"
	elif z>0.0: return "blue"
	else: return "black"

def quiescent(z,flux_u1,flux_v1,flux_v2,flux_j2):
	y = -2.5*np.log10(flux_u1/flux_v1)
	x = -2.5*np.log10(flux_v2/flux_j2)
	if x > 1.6 or y < 1.3: return False
	if z > 1.0:
		y0 = 0.49
	elif z > 0.5:
		y0 = 0.59
	else:
		y0 = 0.69
	return y > 0.88 * x + y0

def errbars(r, n1s, nts):
	import math
	# Following formula (26) of Gehrels1986
	def h1(n2): 
		return (1/((2*n2)-1))	
	def h2(n1): 
		return (1/((2*n1)+1))
	def h(n1,n2): 
		return 2/(h1(n2)+h2(n1))
	def w(n1,n2):
		ht = h(n1,n2)
		t1 = math.sqrt(ht-(1/3)) / ht
		t2 = (h1(n2)-h2(n1))
		t3 = (1/2)-(2/(3*ht))
		return t1+(t2*t3)
	CL = 0.8413
	def p1u_b(n1,n2):
		if n2 == 0:
			return 1
		elif n2 == 1:
			return math.pow(CL, 1/(n1+n2))
		elif n1 == 0:
			return 1-(math.pow(1-CL,(1/n2)))
		num = (n1+1)*math.exp(2*w(n1,n2))
		den = n2+num
		return num/den
	def p1l_b(n1,n2):
		return 1-p1u_b(n2,n1)
	y = [n1/nt for (y,n1,nt) in zip(r,n1s,nts) ]
	p2u = [p1u_b(n1,nt-n1)-rt for (rt, n1,nt) in zip(y, n1s,nts)]
	p2l = [rt-p1l_b(n1,nt-n1) for (rt, n1,nt) in zip(y, n1s,nts)]
	p2u_lim = [n1 == nt for (n1,nt) in zip(n1s,nts)]
	p2l_lim = [n1 == 0 for (n1,nt) in zip(n1s,nts)]
	print(*n1s)
	print(*nts)
	#print(r)
	print(*y)
	print(*p2u)
	print(*p2l)
	return (y, p2u, p2l, p2u_lim, p2l_lim)

def plotsteps(plotlist = [], labelx = "", labely = ""):

	fig, ax = plt.subplots()
	for (x1, x2, r, n1s, nts, c, offset, e) in plotlist:
		rmax = 0.0
		for i in range(0, len(n1s)-1):
			rmax = r[i] if r[i]>rmax else rmax
		ax.plot([*x1,x2[-1]],[*r,r[-1]], c=c[0], drawstyle='steps-post', label=e[0])
		x3 = [o+(i+j)/2 for (i,j,o) in zip(x1, x2, offset)]
		print(e[0])
		print(*x1)
		print(*x2)
		(y, p2u, p2l, p2u_lim, p2l_lim) = errbars(r, n1s, nts)
		#ax.errorbar(x3, r, yerr=[p2l, p2u], uplims=p2u_lim, lolims=p2l_lim, ecolor=c[0], fmt='none')
		ax.errorbar(x3, r, yerr=[p2l, p2u], ecolor=c[0], fmt='none')
		x4  = [(i+j)/2+o for (i,j,o, k,l) in zip(x1, x2, offset, p2l_lim, p2u) if k]
		p2u = [l for (i,j,k,l) in zip(x1, x2, p2l_lim, p2u) if k]
		ax.errorbar(x4, p2u, yerr=p2u, uplims=True, ecolor=c[0], fmt='none')
		x4  = [(i+j)/2+o for (i,j,o, k,l) in zip(x1, x2, offset, p2u_lim, p2l) if k]
		p2l = [l for (i,j,k,l) in zip(x1, x2, p2u_lim, p2l) if k]
		ax.errorbar(x4, p2l, yerr=p2l, lolims=True, ecolor=c[0], fmt='none')
	legend = ax.legend(loc='upper right', shadow=True)
	ax.set_xlim(0.2,3.0)
	print("rmax=%f" % rmax)
	#ax.set_ylim(0.0,rmax*1.2)
	ax.set_ylim(0.0,1.19)
	#ax.set_ylim(0.0,rmax+0.1)
	#ax.axis([0.2,3.0,0.0,rmax*1.1])
	ax.set_xlabel(r'$%s$' % labelx)
	ax.set_ylabel(r'$%s$' % labely)
	plt.show()

def plotwbars(x1, x2, r, n1s, nts, xyrange, labelx="", labely=""):
	fig, ax = plt.subplots()
	fit = np.polyfit(x1, r, deg=1)
	line = []
	#for i in x:
	#	line.append(fit[0]*float(i)+fit[1])
	#ax.plot(x, line, color='black')

	(y, p2u, p2l, p2u_lim, p2l_lim) = errbars(r, n1s, nts)
	#p1u = [(np.sqrt(n)/nt) for (rt,n,nt) in zip(y,n1s,nts) ]
	#p1l = [(np.sqrt(n)/nt) for (rt,n,nt) in zip(y,n1s,nts) ]
	#ax.errorbar(x, y, yerr=[p1u, p1l], marker='s')
	ax.plot([*x1,x2[-1]],[*r,r[-1]], drawstyle='steps-post')
	x3 = [(i+j)/2 for (i,j) in zip(x1, x2)]
	ax.errorbar(x3, r, yerr=[p2l, p2u], uplims=p2u_lim, lolims=p2l_lim, fmt='none')
	#ax.scatter(x, r, marker='s')
	ax.annotate(r'$N_{%s}/N_{Tot} = %.2f %.2f\%s,$' % (labely, fit[1], fit[0], labelx), xy=xyrange )
	ax.set_xlabel(r'$%s$' % labelx)
	ax.set_ylabel(r'$%s$' % labely)
	plt.show()

from astroML.crossmatch import crossmatch_angular
from astroML.plotting import hist

def crossmatch(alldata, kth=None, cat='all', cclass18=None, cclass22=None, selection=None):
	print("begin crossmatch")
	
	selset1 = None
	if selection is not None:
		selset0 = set(selection[0])
		if len(selection) > 1:
			print(Counter(selection[1]))
			selset1 = True

	def inselection(x): return selection is None or x in selset0
	def selectifcolor(x): return x == 'r' or x == 'y'
	def selectioncolor(x,i): 
		if selset1 is not None:
			return selection[1][selection[0].index(i)]
		else: 
			return x

	fields, data = FilterByFields(alldata, "id","lmass","z_peak","ra","dec","Ks_tot","USE","star")
	conds = { 
			"id":inselection,
			#"Ks_tot": lambda k: kscut(k, th=kth),
			#"z_peak":zcut,
			"USE":is1
			}
	r1 = FilterByConds(fields, data, conds)
	r2 = FilterByFuncs(fields, r1, {"star":selectioncolor}, ["id"])
	[idx, lmass, zpeak, ra, dec, kstot, use, star] = zip(*r2)
	imX = np.empty((len(idx), 2), dtype=np.float64)
	imX[:, 0] = ra
	imX[:, 1] = dec

	#f2 = open('../xray/ccosmos_bright_source_catalog_v2.1b.tbl', 'r')
	f2 = open('data/C_COSMOS_identification_catalog_irsab.dat', 'r') 
	d2 = f2.read()
	f2.close()
	chandra = ascii.read(d2, guess=False, format='basic')
	if cclass18 is not None:
		if len(cclass18) == 0:
			cclass18 = [0,1,2]
			#cclass18 = [1,2]
		mycolors_chandra = ['[Chandra Star]', '[Chandra BLAGN]', '[Chandra not-BLAGN]']
		chandra_ra = [i for (i,s) in zip(chandra['RA'], chandra['z_spec_class']) if s in cclass18 ]
		chandra_dec = [i for (i,s) in zip(chandra['Dec'], chandra['z_spec_class']) if s in cclass18 ]
		chandra_col = [mycolors_chandra[s] for (i,s) in zip(chandra['Dec'], chandra['z_spec_class']) if s in cclass18 ]
	elif cclass22 is not None:
		if len(cclass22) == 0:
			#cclass22 = [1,2,3]
			cclass22 = [1,2]
		mycolors_chandra = ['', '[Chandra unobs]', '[Chandra obsc]', '[Chandra galaxy]']
		chandra_ra = [i for (i,s) in zip(chandra['RA'], chandra['z_phot_class']) if s in cclass22 ]
		chandra_dec = [i for (i,s) in zip(chandra['Dec'], chandra['z_phot_class']) if s in cclass22 ]
		chandra_col = [mycolors_chandra[s] for (i,s) in zip(chandra['Dec'], chandra['z_phot_class']) if s in cclass22 ]
	#chandra_ra = chandra['RA']
	#chandra_dec = chandra['Dec']
	#chandra_ra = chandra['ra']
	#chandra_dec = chandra['dec']
	stX_c = np.empty((len(chandra_ra), 3), dtype=np.float64)
	stX_c[:, 0] = chandra_ra
	stX_c[:, 1] = chandra_dec 
	stX_cc = chandra_col

	f3	 = open('data/cosmos_xmm_200811_v2b.tbl', 'r')
	d3 = f3.read()
	f3.close()
	xmm = ascii.read(d3, guess=False, format='basic')
	#xmm_flux = [s for (i,s) in zip(xmm['ra'], xmm['s210']) if s >0.0]
	xmm_ra = xmm['ra']
	xmm_dec = xmm['dec']
	stX_x = np.empty((len(xmm_ra), 2), dtype=np.float64)
	stX_x[:, 0] = xmm_ra
	stX_x[:, 1] = xmm_dec 
	stX_xc = len(xmm_ra) * [mycolors_xagn[0]]
	#stX_xc = len(xmm_ra) * ['[Xray AGN]']

	if cat == 'all':
		print('All chandra and XMM ')
		stX = np.empty((len(xmm_ra)+len(chandra_ra), 2), dtype=np.float64)
		all_ra = list(xmm_ra) + list(chandra_ra)
		stX[:, 0] = all_ra 
		all_dec = list(xmm_dec) + list(chandra_dec)
		stX[:, 1]= all_dec
		stXc = stX_xc + stX_cc
	elif cat == 'xmm':
		print('XMM only')
		stX=stX_x
		stXc = stX_xc
	elif cat == 'chandra':
		print("cclass18=%s, cclass22=%s" % (cclass18, cclass22))
		stX=stX_c
		stXc = stX_cc
	elif cat == 'none':
		stX = np.empty((0, 2), dtype=np.float64)
	
	print("Chandra, XMM, Total (in catalogs)")
	print(len(stX_c),len(stX_x),len(stX))
	print(Counter(stXc))

	# crossmatch catalogs
	max_radius = 1.0 / 3600  # 1 arcsec
	dist, ind = crossmatch_angular(stX, imX, max_radius)
	match = ~np.isinf(dist)
	print(Counter(match))	

	seen = {}
	deg_d = []
	if selection is not None:
		selection = idx
		if selset1:
			colors = list(star)
		else:
			colors = len(idx) * ['[not AGN]']
		#print(len(selection),len(colors))
		dists = len(idx) * [0.0,]
		for i in range(0,len(ind)):
			if not np.isinf(dist[i]):
				i2 = idx[ind[i]]
				if colors[ind[i]] not in mycolors_xagn:
					seen[i2] = [colors[ind[i]],]
				else:
					seen[i2].append(colors[ind[i]])
					deg_d.append(dist[i])
				#print(colors[ind[i]], stXc[i])
				colors[ind[i]] = stXc[i]
				#print(i, stX[i], ind[i], ra[ind[i]], dec[ind[i]])
				#print(i2, alldata["id"][i2-1], alldata["ra"][i2-1], alldata["dec"][i2-1])
				#colors[ind[i]] = mycolors_xagn
				dists[ind[i]] = dist[i]
		xagn_d = [ 3600*d for d in deg_d]
	else:
		selection = []
		colors = []
		dists = []
		for i in range(len(ind)-1,-1,-1):
			if not np.isinf(dist[i]):
				i2 = idx[ind[i]]
				if i2 in selection:
					deg_d.append(dist[i])
				else:	
					colors.append(stXc[i])
					selection.append(i2)
					dists.append(dist[i])
		dupes = []
		for x in range(0,len(selection)):
			if selection[x] not in seen:
				seen[selection[x]] = [colors[x],]
			else:
				seen[selection[x]].append(colors[x])
				dupes.append(dists[x])
		xagn_d = [ 3600*d for d in deg_d]
	
	#print(len([seen[s] for s in seen.keys()]))
	#print(len([seen[s] for s in seen.keys() if len(seen[s])>1]))
	#print(len([seen[s] for s in seen.keys() if len(seen[s])==1]))
	
	#dist_match = dist[match]
	#dist_match *= 3600
	
	if not selset1:
		ax = plt.axes()
		nbins = 20
		for k in sorted(mycolors.keys()):
		#for k in mycolors_xagn:
			x = [ 3600*d for (d,s,c) in zip(dists,selection,colors) if c == k]
			if len(x)>0:
				ax.hist(x, bins=nbins, color=mycolors[k][0], histtype='step', label="%s, N=%d" % (mycolors[k][1], len(x))) 
		#ax.hist(xmm_d, bins=nbins, color=mycolors['[XMM]'][0], histtype='step', label="%s, N=%d" % ('[XMM]', len(xmm_d))) 
		#ax.hist(cha0_d, bins=nbins, color=mycolors['[Chandra Star]'][0], histtype='step', label="%s, N=%d" % ('[Chandra Star]', len(cha0_d)))
		#ax.hist(cha1_d, bins=nbins, color=mycolors['[Chandra BLAGN]'][0], histtype='step', label="%s, N=%d" % ('[Chandra BLAGN]', len(cha1_d)))
		#ax.hist(cha2_d, bins=nbins, color=mycolors['[Chandra not-BLAGN]'][0], histtype='step', label="%s, N=%d" % ('[Chandra not-BLAGN]', len(cha2_d)))
		#ax.hist(xagn_d, bins=nbins, color='[Xray AGN]', histtype='step', label="%s, N=%d" % ("Multiple match", len(xagn_d)))
		#hist(dist_match, bins='knuth', ax=ax,
		#     histtype='stepfilled', ec='k', fc='#AAAAAA')
		ax.set_xlabel('radius of match (arcsec)')
		ax.set_ylabel('N(r, r+dr)')
		#ax.text(0.95, 0.95, "Total objects: %i\nNumber with match: %i" % (imX.shape[0], np.sum(match)),
	 #    ha='right', va='top', transform=ax.transAxes)
		ax.legend(loc=1, ncol=1) 
		ax.set_xlim(0, 1.1)
		plt.show()
	
	print(Counter(colors))
	print("Duplicates: %d" % len(xagn_d))
	print("end crossmatch")
	return [selection, colors]
