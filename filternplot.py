#!/usr/bin/python3

from loaddata import * 

def stern2005_1x1_alluse(alldata, kth=None, selection=None):
	print("Plotting Stern+2005 space for all redshifts, all USE")
	selset1 = None
	if selection is not None:
		selset0 = set(selection[0])
		if len(selection) > 1:
			print(Counter(selection[1]))
			selset1 = True

	def inselection(x): return selection is None or x in selset0
	def selectioncolor(x,i): 
		if selset1 is not None:
			return selection[1][selection[0].index(i)]
		else: 
			return x
	
	fields, data = FilterByFields(alldata, "id","ch1","ech1","ch2","ech2","ch3","ech3","ch4","ech4","contamination","star","Ks","Ks_tot","USE")
	conds = {
		"id":inselection,
		"ch1":pos,
		"ch2":pos,
		"ch3":pos,
		"ch4":pos,
		#"USE":is1,
		"contamination":is0,
		"Ks_tot": lambda k: kscut(k, th=kth)
		#"Ks_tot": lambda k: kscut2(k, th=kth)
		}
	r1 = FilterByConds(fields, data, conds)
	r2 = FilterByFuncs(fields, r1, 
			{
				"ch1":lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: magAB(f,ks,kst)-2.78,
				"ch2":lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: magAB(f,ks,kst)-3.26,
				"ch3":lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: magAB(f,ks,kst)-3.75,
				"ch4":lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: magAB(f,ks,kst)-4.38,
				"ech1":sn1,
				"ech2":sn2,
				"ech3":sn3,
				"ech4":sn4,
				"star":lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: color2(f)
				}, 
			["Ks","Ks_tot","ch1","ech1","ch2","ech2","ch3","ech3","ch4","ech4"])
	#r3 = FilterByConds(fields, r2, {"ech1":sncut, "ech2":sncut, "ech3":sncut, "ech4":sncut})
	r3 = FilterByFuncs(fields, r2, {"star":selectioncolor}, ["id"])
	[idx,mag_36,e1,mag_45,e2,mag_58,e3,mag_80,e4,contamination,star,ks,ks_tot,use] = zip(*r3)
	
	m58_80 = [ i-j for (i,j) in zip(mag_58,mag_80)]
	m36_45 = [ i-j for (i,j) in zip(mag_36,mag_45)]
	x = m58_80
	y = m36_45
	if selset1 is not None:
		u = [ mycolors[z][0] for (u,z) in zip(use,star)] # if z > b["z0"] and z <= b["z1"] ]
		plt.scatter(x, y, marker=".", c=u, edgecolors='none')
	else:
		plt.scatter(x, y, marker=".", c=star, edgecolors='none')
	plt.xlim(xmin=-1.0,xmax=3.5)
	plt.ylim(ymin=-0.4,ymax=1.5)
	plt.xlabel("[5.8]-[8.0] (Vega)")
	plt.ylabel("[3.6]-[4.5] (Vega)")
	plt.plot([0.6, 0.6], [0.3,1.5], color='gray', linestyle='-', lw=2)
	plt.plot([0.6, 1.6], [0.3,0.5], color='gray', linestyle='-', lw=2)
	plt.plot([1.6, 2.0], [0.5,1.5], color='gray', linestyle='-', lw=2)
	plt.show()

def stern2005_1x1(alldata):
	print("Plotting Stern+2005 space for all redshifts, USE=1")
	fields, data = FilterByFields(alldata,"ch1","ech1","ch2","ech2","ch3","ech3","ch4","ech4","USE","z_spec","Ks","Ks_tot")
	conds = {
		"ch1":pos,
		"ch2":pos,
		"ch3":pos,
		"ch4":pos,
		"USE":is1,
		"Ks_tot":kscut
		}
	r1 = FilterByConds(fields, data, conds)
	r2 = FilterByFuncs(fields, r1, 
			{
				"ch1":lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: magAB(f,ks,kst)-2.78,
				"ch2":lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: magAB(f,ks,kst)-3.26,
				"ch3":lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: magAB(f,ks,kst)-3.75,
				"ch4":lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: magAB(f,ks,kst)-4.38,
				"ech1":sn1,
				"ech2":sn2,
				"ech3":sn3,
				"ech4":sn4,
				#"z_spec":lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: color(f)
				}, 
			["Ks","Ks_tot","ch1","ech1","ch2","ech2","ch3","ech3","ch4","ech4"])
	r3 = FilterByConds(fields, r2, {"ech1":sncut, "ech2":sncut, "ech3":sncut, "ech4":sncut})
	[mag_36,e1,mag_45,e2,mag_58,e3,mag_80,e4,use,z_spec,ks,ks_tot] = zip(*r3)
	
	m58_80 = [ i-j for (i,j) in zip(mag_58,mag_80)]
	m36_45 = [ i-j for (i,j) in zip(mag_36,mag_45)]
	x = m58_80
	y = m36_45
	print(len(x))
	#plt.scatter(x, y, marker=".", c=z_spec)
	plt.scatter(x, y, marker=".", c='b')
	plt.xlim(xmin=-1.0,xmax=3.5)
	plt.ylim(ymin=-0.4,ymax=1.5)
	plt.xlabel("[5.8]-[8.0] (Vega)")
	plt.ylabel("[3.6]-[4.5] (Vega)")
	plt.plot([0.6, 0.6], [0.3,1.5], color='darkgray', linestyle='-', lw=2)
	plt.plot([0.6, 1.6], [0.3,0.5], color='darkgray', linestyle='-', lw=2)
	plt.plot([1.6, 2.0], [0.5,1.5], color='darkgray', linestyle='-', lw=2)
	plt.show()

def stern2005_3x2(alldata, kth=None, selection=None, labelx=None, labely=None):
	print("Plotting Stern+2005 space for each redshift, USE=1")
	selset1 = None
	if selection is not None:
		selset0 = set(selection[0])
		if len(selection) > 1:
			print(Counter(selection[1]))
			selset1 = True

	def inselection(x): return selection is None or x in selset0
	#def selectifcolor(x): return x == 'r' or x == 'y'
	def selectioncolor(x,i): 
		if selset1 is not None:
			return selection[1][selection[0].index(i)]
		else: 
			return x


	#(f1,ef1,f2,ef2,f3,ef3,f4,ef4) = ("ch3","ech3","ch1","ech1","ch4","ech4","ch2","ech2")
	(f1,ef1,f2,ef2,f3,ef3,f4,ef4) = ("ch1","ech1","ch2","ech2","ch3","ech3","ch4","ech4")
	if labelx is None: labelx = r'$%s/%s$' % (f1,f2)
	if labely is None: labely = r'$%s/%s$' % (f3,f4)
	fields0, data0 = FilterByFields(alldata, "id",f1,ef1,f2,ef2,f3,ef3,f4,ef4,"z_peak","USE","star","Ks","Ks_tot")
	r10 = FilterByConds(fields0, data0, {
		f1:pos,
		f2:pos,
		f3:pos,
		f4:pos,
		"Ks_tot": lambda k: kscut(k, th=kth),
		"USE":is1,
		})
	r20 = FilterByFuncs(fields0, r10, 
			{
				f1:lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: magAB(f,ks,kst)-2.78,
				f2:lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: magAB(f,ks,kst)-3.26,
				f3:lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: magAB(f,ks,kst)-3.75,
				f4:lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: magAB(f,ks,kst)-4.38,
				ef1:sn1,
				ef2:sn2,
				ef3:sn3,
				ef4:sn4,
				}, 
			["Ks","Ks_tot",f1,ef1,f2,ef2,f3,ef3,f4,ef4])
	#r30 = FilterByConds(fields0, r20, { "USE":upsncut }, [ef1,ef2,ef3,ef4])
	r30 = FilterByConds(fields0, r20, {ef1:sncut, ef2:sncut, ef3:sncut, ef4:sncut})
	[idx0,log_10,e10,log_20,e20,log_30,e30,log_40,e40,zpeak0,use0,star0,ks0,ks_tot0] = zip(*r30)
	fields, data = FilterByFields(alldata, "id",f1,ef1,f2,ef2,f3,ef3,f4,ef4,"z_peak","USE","star","Ks","Ks_tot")
	r1 = FilterByConds(fields, data, {
		f1:pos,
		f2:pos,
		f3:pos,
		f4:pos, 
		"id":inselection, 
		#"Ks_tot": lambda k: kscut(k, th=kth),
		#"USE":is1, 
	})
	r2 = FilterByFuncs(fields, r1, 
			{
				f1:lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: magAB(f,ks,kst)-2.78,
				f2:lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: magAB(f,ks,kst)-3.26,
				f3:lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: magAB(f,ks,kst)-3.75,
				f4:lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: magAB(f,ks,kst)-4.38,
				ef1:sn1,
				ef2:sn2,
				ef3:sn3,
				ef4:sn4,
				}, 
			["Ks","Ks_tot",f1,ef1,f2,ef2,f3,ef3,f4,ef4])
	#r31 = FilterByConds(fields, r2, { "USE":upsncut }, [ef1,ef2,ef3,ef4])
	#r5 = FilterByFuncs(fields, r31, { "USE":upsarrow, }, [ef1,ef2,ef3,ef4])
	r5 = FilterByConds(fields, r2, {ef1:sncut, ef2:sncut, ef3:sncut, ef4:sncut})
	#r6 = FilterByFuncs(fields, r5, {"star":color2})
	r6 = FilterByFuncs(fields, r5, {"star":selectioncolor}, ["id"])
	[idx,log_1,e1,log_2,e2,log_3,e3,log_4,e4,zpeak,use,star,ks,ks_tot] = zip(*r6)
	ntot = len(idx)
	print("Ntot=%d" % ntot)

	l12 = [ i-j for (i,j) in zip(log_1,log_2)]
	l34 = [ i-j for (i,j) in zip(log_3,log_4)]

	fig, axis = plt.subplots(2,3)
	binres = (100, 100)
	zbins = [ 
			{"z0":0.2, "z1":0.5, "ax":[0,0], "y0":0.69},
			{"z0":0.5, "z1":1.0, "ax":[0,1], "y0":0.59},
			{"z0":1.0, "z1":1.5, "ax":[0,2], "y0":0.49},
			{"z0":1.5, "z1":2.0, "ax":[1,0], "y0":0.49},
			{"z0":2.0, "z1":2.5, "ax":[1,1], "y0":0.49},
			{"z0":2.5, "z1":3.0, "ax":[1,2], "y0":0.49}
			]	
	for b in zbins:
		ax = axis[b["ax"][0],b["ax"][1]]
		y0 = [ i-j for (i,j,k,l,z) in zip(log_10,log_20,log_30,log_40,zpeak0) if 
				not np.isnan(i+j+k+l) and z >= b["z0"] and z < b["z1"]]
		x0 = [ k-l for (i,j,k,l,z) in zip(log_10,log_20,log_30,log_40,zpeak0) if
				not np.isnan(i+j+k+l) and z >= b["z0"] and z < b["z1"]]
		#ax.hist2d(x0, y0, bins=binres, cmap=plt.cm.Greys,zorder= 0)
		if selset1 is None:
			y = [ i for (i,j,z) in zip(l12,use,zpeak) if 
					z >= b["z0"] and z < b["z1"]]
			x = [ i for (i,j,z) in zip(l34,use,zpeak) if 
					z >= b["z0"] and z < b["z1"]]
			c = [ j for (i,j,z) in zip(l34,star,zpeak) if
					z >= b["z0"] and z < b["z1"]]        
			ax.scatter(x, y, marker='.', c=c, edgecolors="none", alpha=0.5, label="")
				#ax.scatter(x, y, marker='.', c='b', edgecolors='none')
		else:
			for k in mycolors.keys():
				y = [ i for (i,j,z) in zip(l12,star,zpeak) if 
						j == k  and z >= b["z0"] and z < b["z1"]]
				x = [ i for (i,j,z) in zip(l34,star,zpeak) if 
						j == k  and z >= b["z0"] and z < b["z1"]]
				if len(x) > 0:
					ax.scatter(x, y, marker=mycolors[k][2], c=mycolors[k][0], edgecolors='none', label=mycolors[k][1], alpha=0.8)
					#ax.scatter(x, y, marker='.', c=mycolors[k][0], edgecolors='none', label=mycolors[k][1])
		ax.set_title(r'$%.1f\leq z<%.1f$' % (b["z0"],b["z1"]))
		#ax.axis([xmin,xmax,ymin,ymax])
		#ax.set_aspect("equal")
		legend = ax.legend(loc='lower right', fontsize='small', shadow=True)
		ax.axis([-1.0,3.5,-0.4,1.5])
		ax.set_xlabel("[5.8]-[8.0] (Vega)")
		ax.set_ylabel("[3.6]-[4.5] (Vega)")
		ax.plot([0.6, 0.6], [0.3,1.5], color='darkgray', linestyle='-', lw=2)
		ax.plot([0.6, 1.6], [0.3,0.5], color='darkgray', linestyle='-', lw=2)
		ax.plot([1.6, 2.0], [0.5,1.5], color='darkgray', linestyle='-', lw=2)
	plt.show()
	selection = [ i for (i,j) in zip(idx,use) if j == '.' ]
	print("selection=%d" % len(selection))
	return [selection, ] 
#end stern2005_3x2

def lacy2007_1x1(alldata, kth=None, selection=None):    
	print("Plotting Lacy+2007 space for all redshifts, USE=1")    
	selset1 = None    
	if selection is not None:
		selset0 = set(selection[0])
		if len(selection) > 1:
			print(Counter(selection[1]))
			selset1 = True

	def inselection(x): return selection is None or x in selset0
	#def selectifcolor(x): return x == 'r' or x == 'y'
	def selectioncolor(x,i): 
		if selset1 is not None:
			return selection[1][selection[0].index(i)]
		else: 
			return x
	
	fields, data = FilterByFields(alldata, "id","ch1","ech1","ch2","ech2","ch3","ech3","ch4","ech4","USE","star","Ks","Ks_tot")
	conds = {
			"id":inselection,
			#"ch1":pos,
			#"ch2":pos,
			#"ch3":pos,
			#"ch4":pos,
			"USE":is1,
			#"Ks_tot": lambda k: kscut(k, th=kth)
			}
	r1 = FilterByConds(fields, data, conds)
	r2 = FilterByConds(fields, r1,
			{"Ks": lambda f, ks, kstot: kscut2(f, ks, kstot, th=kth)},
			["Ks","Ks_tot"]
			)
	r3 = FilterByFuncs(fields, r2, 
			{
				"ch1":uplog1b,
				"ch2":uplog2b,
				"ch3":uplog3b,
				"ch4":uplog4b,
				"ech1":upsn1,
				"ech2":upsn2,
				"ech3":upsn3,
				"ech4":upsn4,
				#"z_spec":lambda f,ks,kst,c1,e1,c2,e2,c3,e3,c4,e4: color(f)
				}, 
			["Ks","Ks_tot","ch1","ech1","ch2","ech2","ch3","ech3","ch4","ech4"])
	r31 = FilterByConds(fields, r3, 
			{
				"USE":upsncut
				},
			["ech1","ech2","ech3","ech4"])
	r32 = FilterByConds(fields, r31, 
			{
				"USE":donley2012_cond
				#"USE":lacy2007_cond
				},
			["ch1","ch2","ch3","ch4"])
	[idxw,log_36,e1,log_45,e2,log_58,e3,log_80,e4,use,star,ks,ks_tot] = zip(*r32)
	nsel = len(idxw)
	print("Nsel=%d" % nsel)
	r5 = FilterByFuncs(fields, r31, 
			{
				"USE":upsarrow,
				},
			["ech1","ech2","ech3","ech4"])
	r6 = FilterByFuncs(fields, r5, {"star":selectioncolor}, ["id"])
	#r4 = FilterByConds(fields, r3, 
	#		{
	#			"ch1":donley2012_cond,
	#			"ch2":donley2012_cond,
	#			"ch3":donley2012_cond,
	#			"ch4":donley2012_cond
	#			},
	#		["ch1","ch2","ch3","ch4"]
	#		)
	[idx,log_36,e1,log_45,e2,log_58,e3,log_80,e4,use,star,ks,ks_tot] = zip(*r6)
	ntot = len(idx)
	print("Ntot=%d" % ntot)

	l58_36 = [ i-j for (i,j) in zip(log_58,log_36)]
	l80_45 = [ i-j for (i,j) in zip(log_80,log_45)]
	#u = [ color(s) for (u,s) in zip(use,star)] # if z > b["z0"] and z <= b["z1"] ]
	u = [ s for (u,s) in zip(use,star)] # if z > b["z0"] and z <= b["z1"] ]
	#plt.scatter(x, y, marker=".", c=use)
	#x = np.array(l58_36)
	#y = np.array(l80_45)
	if not selset1:
		for k in arrows.keys():
			x = [ i for (i,j) in zip(l58_36,use) if j == k ]
			y = [ i for (i,j) in zip(l80_45,use) if j == k ]
			#plt.scatter(x, y, marker=arrows[k], c='b', edgecolors="none")
			plt.scatter(x, y, marker='.', c='b', edgecolors='none')
	else:
		for k in sorted(mycolors.keys()):
			x = [ i for (i,j) in zip(l58_36,star) if j == k ]
			y = [ i for (i,j) in zip(l80_45,star) if j == k ]
			if len(x) > 0:
				plt.scatter(x, y, marker=mycolors[k][2], c=mycolors[k][0], edgecolors='none', label=mycolors[k][1])
				#plt.scatter(x, y, marker=mycolors[k][2], c=mycolors[k][0], edgecolors='none', label=mycolors[k][1])
		legend = plt.legend(loc='upper right', shadow=True)
		#plt.title(r'$K_{S_{tot}}\,<\,%.1f\,\,N_{tot}\,=\,%d\,\,N_{sel}\,=\,%d$' % (kth, ntot, nsel))
		#plt.title(r'$N_{tot}\,=\,%d\,\,N_{sel}\,=\,%d$' % (ntot, nsel))
	xmax = 1.0
	plt.xlim(xmin=-0.6,xmax=xmax)
	plt.ylim(ymin=-0.8,ymax=1.4)
	plt.xlabel(r'$log(f_{5.8\mu m}/f_{3.6\mu m})$', fontsize="x-large")
	plt.ylabel(r'$log(f_{8.0\mu m}/f_{4.5\mu m})$', fontsize="x-large")
	# Original wedge
	plt.plot([-0.1, -0.1], [-0.2,0.4], color='gray', linestyle='--', lw=2)
	plt.plot([-0.1, 1.0], [-0.2,-0.2], color='gray', linestyle='--', lw=2)
	plt.plot([-0.1, 1.0], [0.4,1.3], color='gray', linestyle='--', lw=2)
	# Power law
	plt.plot([0.1, 0.6], [0.125,0.75], color='black', linestyle='-', lw=2)
	# Revised wedge
	plt.plot([0.08, 0.3471], [0.15,0.15], color='gray', linestyle='-', lw=2)
	plt.plot([0.08, 0.08], [0.15,0.3668], color='gray', linestyle='-', lw=2)
	plt.plot([0.3471, xmax], [0.15,1.21*xmax-0.27], color='gray', linestyle='-', lw=2)
	plt.plot([0.08, xmax], [0.3668,1.21*xmax+0.27], color='gray', linestyle='-', lw=2)
	plt.plot([xmax, xmax], [1.21*xmax-0.27,1.21*xmax+0.27], color='gray', linestyle='-', lw=2)
	plt.show()
	selection = [ i for (i,j) in zip(idxw,use) if j == '.' ]
	print("selection=%d" % len(selection))
	return [selection, ] 
#end lacy2007_1x1

def color_color_1x1(alldata, kth=None, selection=None):
	print("Plotting color_color space for all redshifts, USE=1")
	(f1,ef1,f2,ef2,f3,ef3,f4,ef4) = ("J","eJ","H","eH","Ks","eKs","ch1","ech1")
	print("Plotting: %s-%s vs %s-%s" % (f3,f4,f1,f2))

	selset1 = None
	if selection is not None:
		selset0 = set(selection[0])
		if len(selection) > 1:
			print(Counter(selection[1]))
			selset1 = True

	def inselection(x): return selection is None or x in selset0
	#def selectifcolor(x): return x == 'r' or x == 'y'
	def selectioncolor(x,i): 
		if selset1 is not None:
			return selection[1][selection[0].index(i)]
		else: 
			return x

	fields0, data0 = FilterByFields(alldata, "id",f1,ef1,f2,ef2,f3,ef3,f4,ef4,"USE","star","Ks","Ks_tot")
	r10 = FilterByConds(fields0, data0, {"USE":is1})
	r20 = FilterByFuncs(fields0, r10, 
			{
				f1:uplog1b,
				f2:uplog2b,
				f3:uplog3b,
				f4:uplog4b,
				ef1:upsn1,
				ef2:upsn2,
				ef3:upsn3,
				ef4:upsn4,
				}, 
			["Ks","Ks_tot",f1,ef1,f2,ef2,f3,ef3,f4,ef4])
	r30 = FilterByConds(fields0, r20, { "USE":upsncut }, [ef1,ef2,ef3,ef4])
	[idx0,log_10,e10,log_20,e20,log_30,e30,log_40,e40,use0,star0,ks0,ks_tot0] = zip(*r30)
	
	fields, data = FilterByFields(alldata, "id",f1,ef1,f2,ef2,f3,ef3,f4,ef4,"USE","star","Ks","Ks_tot")
	r1 = FilterByConds(fields, data, { "id":inselection, "USE":is1, })
	r2 = FilterByFuncs(fields, r1, 
			{
				f1:uplog1b,
				f2:uplog2b,
				f3:uplog3b,
				f4:uplog4b,
				ef1:upsn1,
				ef2:upsn2,
				ef3:upsn3,
				ef4:upsn4,
				}, 
			["Ks","Ks_tot",f1,ef1,f2,ef2,f3,ef3,f4,ef4])
	r31 = FilterByConds(fields, r2, { "USE":upsncut }, [ef1,ef2,ef3,ef4])
	
	r5 = FilterByFuncs(fields, r31, { "USE":upsarrow, }, [ef1,ef2,ef3,ef4])
	r6 = FilterByFuncs(fields, r5, {"star":selectioncolor}, ["id"])
	[idx,log_1,e1,log_2,e2,log_3,e3,log_4,e4,use,star,ks,ks_tot] = zip(*r6)
	ntot = len(idx)
	print("Ntot=%d" % ntot)

	l12 = [ i-j for (i,j) in zip(log_1,log_2)]
	l34 = [ i-j for (i,j) in zip(log_3,log_4)]
	
	binres = (100, 100)
	x0 = [ i-j for (i,j,k,l) in zip(log_10,log_20,log_30,log_40) if not np.isnan(i+j+k+l)]
	y0 = [ k-l for (i,j,k,l) in zip(log_10,log_20,log_30,log_40) if not np.isnan(i+j+k+l)]
	#(xmin, xmax, ymin, ymax) = (min(l12+x0), max(l12+x0), min(l34+y0), max(l34+y0))
	(xmin, xmax, ymin, ymax) = (min(l12), max(l12), min(l34), max(l34))
	plt.hist2d(x0, y0, bins=binres, cmap=plt.cm.Greys,zorder= 0)
	
	if not selset1:
		for k in arrows.keys():
			x = [ i for (i,j) in zip(l12,use) if j == k ]
			y = [ i for (i,j) in zip(l34,use) if j == k ]
			plt.scatter(x, y, marker=arrows[k], c='b', edgecolors="none", alpha=0.5)
			#plt.scatter(x, y, marker='.', c='b', edgecolors='none')
	else:
		for k in mycolors.keys():
			x = [ i for (i,j) in zip(l12,star) if j == k ]
			y = [ i for (i,j) in zip(l34,star) if j == k ]
			if len(x) > 0:
				#plt.scatter(x, y, marker=mycolors[k][2], c=mycolors[k][0], edgecolors='none', label=mycolors[k][1])
				plt.scatter(x, y, marker='.', c=mycolors[k][0], edgecolors='none', label=mycolors[k][1])
	plt.xlim(xmin=xmin,xmax=xmax)
	plt.ylim(ymin=ymin,ymax=ymax)
	plt.xlabel(r'$%s/%s$' % (f1,f2))
	plt.ylabel(r'$%s/%s$' % (f3,f4))
	plt.show()
	selection = [ i for (i,j) in zip(idx,use) if j == '.' ]
	print("selection=%d" % len(selection))
	return [selection, ] 
#end color_color_1x1

def color_color_3x2(alldata, kth=None, selection=None, labelx=None, labely=None):
	print("Plotting color_color space at each redshift range, USE=1")
	(f1,ef1,f2,ef2,f3,ef3,f4,ef4) = ("ch3","ech3","ch1","ech1","ch4","ech4","ch2","ech2")
	print("Plotting: %s-%s vs %s-%s" % (f3,f4,f1,f2))
	selset1 = None
	if selection is not None:
		selset0 = set(selection[0])
		if len(selection) > 1:
			print(Counter(selection[1]))
			selset1 = True

	def inselection(x): return selection is None or x in selset0
	#def selectifcolor(x): return x == 'r' or x == 'y'
	def selectioncolor(x,i): 
		if selset1 is not None:
			return selection[1][selection[0].index(i)]
		else: 
			return x

	#(f1,ef1,f2,ef2,f3,ef3,f4,ef4) = ("lmass","chi2","H","eH","L153_uv","eu","L153_vj","eB")
	#(f1,ef1,f2,ef2,f3,ef3,f4,ef4) = ("J","eJ","H","eH","Ks","eKs","ch1","ech1")
	(f1,ef1,f2,ef2,f3,ef3,f4,ef4) = ("ch3","ech3","ch1","ech1","ch4","ech4","ch2","ech2")
	if labelx is None: labelx = r'$%s/%s$' % (f1,f2)
	if labely is None: labely = r'$%s/%s$' % (f3,f4)
	fields0, data0 = FilterByFields(alldata, "id",f1,ef1,f2,ef2,f3,ef3,f4,ef4,"z_peak","USE","star","Ks","Ks_tot")
	r10 = FilterByConds(fields0, data0, {"USE":is1})
	r20 = FilterByFuncs(fields0, r10, 
			{
				f1:uplog1b,
				f2:uplog2b,
				f3:uplog3b,
				f4:uplog4b,
				ef1:upsn1,
				ef2:upsn2,
				ef3:upsn3,
				ef4:upsn4,
				}, 
			["Ks","Ks_tot",f1,ef1,f2,ef2,f3,ef3,f4,ef4])
	r30 = FilterByConds(fields0, r20, { "USE":upsncut }, [ef1,ef2,ef3,ef4])
	[idx0,log_10,e10,log_20,e20,log_30,e30,log_40,e40,zpeak0,use0,star0,ks0,ks_tot0] = zip(*r30)
	
	fields1, data1 = FilterByFields(alldata, "id",f1,ef1,f2,ef2,f3,ef3,f4,ef4,"z_peak","USE","star","Ks","Ks_tot")
	r11 = FilterByConds(fields1, data1, {"star":is1})
	r21 = FilterByFuncs(fields1, r11, 
			{
				f1:uplog1b,
				f2:uplog2b,
				f3:uplog3b,
				f4:uplog4b,
				ef1:upsn1,
				ef2:upsn2,
				ef3:upsn3,
				ef4:upsn4,
				}, 
			["Ks","Ks_tot",f1,ef1,f2,ef2,f3,ef3,f4,ef4])
	#r31 = FilterByConds(fields0, r20, { "USE":upsncut }, [ef1,ef2,ef3,ef4])
	[idx1,log_11,e11,log_21,e21,log_31,e31,log_41,e41,zpeak1,use1,star1,ks1,ks_tot1] = zip(*r21)
	
	fields, data = FilterByFields(alldata, "id",f1,ef1,f2,ef2,f3,ef3,f4,ef4,"z_peak","USE","star","Ks","Ks_tot")
	r1 = FilterByConds(fields, data, { "id":inselection, "USE":is1, })
	r2 = FilterByFuncs(fields, r1, 
			{
				f1:uplog1b,
				f2:uplog2b,
				f3:uplog3b,
				f4:uplog4b,
				ef1:upsn1,
				ef2:upsn2,
				ef3:upsn3,
				ef4:upsn4,
				}, 
			["Ks","Ks_tot",f1,ef1,f2,ef2,f3,ef3,f4,ef4])
	r31 = FilterByConds(fields, r2, { "USE":upsncut }, [ef1,ef2,ef3,ef4])
	r5 = FilterByFuncs(fields, r31, { "USE":upsarrow, }, [ef1,ef2,ef3,ef4])
	r6 = FilterByFuncs(fields, r5, {"star":selectioncolor}, ["id"])
	[idx,log_1,e1,log_2,e2,log_3,e3,log_4,e4,zpeak,use,star,ks,ks_tot] = zip(*r6)
	ntot = len(idx)
	print("Ntot=%d" % ntot)

	l12 = [ i-j for (i,j) in zip(log_1,log_2)]
	l34 = [ i-j for (i,j) in zip(log_3,log_4)]
	x0 = [ i-j for (i,j,k,l,z) in zip(log_10,log_20,log_30,log_40,zpeak0) if not np.isnan(i+j+k+l)]
	y0 = [ k-l for (i,j,k,l,z) in zip(log_10,log_20,log_30,log_40,zpeak0) if not np.isnan(i+j+k+l)]
	#(xmin, xmax, ymin, ymax) = (min(l12+x0), max(l12+x0), min(l34+y0), max(l34+y0))
	(xmin, xmax, ymin, ymax) = (min(l12), max(l12), min(l34), max(l34))
	(xmin, xmax, ymin, ymax) = (-0.6, 1.0, -0.8, 1.4)
	(xmin, xmax, ymin, ymax) = (-0.6, 1.0, -1.0, 1.0)

	fig, axis = plt.subplots(2,3)
	binres = (100, 100)
	zbins = [ 
			{"z0":0.2, "z1":0.5, "ax":[0,0], "y0":0.69},
			{"z0":0.5, "z1":1.0, "ax":[0,1], "y0":0.59},
			{"z0":1.0, "z1":1.5, "ax":[0,2], "y0":0.49},
			{"z0":1.5, "z1":2.0, "ax":[1,0], "y0":0.49},
			{"z0":2.0, "z1":2.5, "ax":[1,1], "y0":0.49},
			{"z0":2.5, "z1":3.0, "ax":[1,2], "y0":0.49}
			]	
	centroids = []
	do_centroids = True
	for b in zbins:
		ax = axis[b["ax"][0],b["ax"][1]]
		x0 = [ i-j for (i,j,k,l,z) in zip(log_10,log_20,log_30,log_40,zpeak0) if 
				not np.isnan(i+j+k+l) and z >= b["z0"] and z < b["z1"]]
		y0 = [ k-l for (i,j,k,l,z) in zip(log_10,log_20,log_30,log_40,zpeak0) if
				not np.isnan(i+j+k+l) and z >= b["z0"] and z < b["z1"]]
		ax.hist2d(x0, y0, bins=binres, cmap=plt.cm.Greys,zorder= 0)
		centroids.append((b["z0"], b["z1"], np.mean(x0), np.mean(y0)))
		#x1 = [ i-j for (i,j,k,l,z) in zip(log_11,log_21,log_31,log_41,zpeak1) if 
		#		not np.isnan(i+j+k+l) and z >= b["z0"] and z < b["z1"]]
		#y1 = [ k-l for (i,j,k,l,z) in zip(log_11,log_21,log_31,log_41,zpeak1) if
		#		not np.isnan(i+j+k+l) and z >= b["z0"] and z < b["z1"]]
		#ax.scatter(x1, y1, marker='.', c='red', edgecolors="none", alpha=0.3, label="stars")
		x = [ i-j for (i,j,k,l,z) in zip(log_10,log_20,log_30,log_40,zpeak0) if 
				not np.isnan(i+j+k+l) and z >= b["z0"] and z < b["z1"]]
		y = [ k-l for (i,j,k,l,z) in zip(log_10,log_20,log_30,log_40,zpeak0) if
				not np.isnan(i+j+k+l) and z >= b["z0"] and z < b["z1"]]
		if selset1 is None:
			for k in arrows.keys():
				x = [ i for (i,j,z) in zip(l12,use,zpeak) if 
						j == k and z >= b["z0"] and z < b["z1"]]
				y = [ i for (i,j,z) in zip(l34,use,zpeak) if 
						j == k and z >= b["z0"] and z < b["z1"]]
				ax.scatter(x, y, marker=arrows[k], c='b', edgecolors="none", alpha=0.5, label="")
				#ax.scatter(x, y, marker='.', c='b', edgecolors='none')
		else:
			for k in sorted(mycolors.keys()):
				x = [ i for (i,j,z) in zip(l12,star,zpeak) if 
						j == k  and z >= b["z0"] and z < b["z1"]]
				y = [ i for (i,j,z) in zip(l34,star,zpeak) if 
						j == k  and z >= b["z0"] and z < b["z1"]]
				if len(x) > 0:
					ax.scatter(x, y, marker=mycolors[k][2], c=mycolors[k][0], edgecolors='none', label=mycolors[k][1])
					#ax.scatter(x, y, marker='.', c=mycolors[k][0], edgecolors='none', label=mycolors[k][1]) 
		ax.set_title(r'$%.1f\leq z<%.1f$' % (b["z0"],b["z1"]))
		ax.axis([xmin,xmax,ymin,ymax])
		#ax.set_aspect("equal")
		legend = ax.legend(loc='lower right', fontsize='small', shadow=True)
		ax.set_xlabel(labelx, fontsize="large")
		ax.set_ylabel(labely, fontsize="large")
		# Original wedge
		ax.plot([-0.1, -0.1], [-0.2,0.4], color='gray', linestyle='--', lw=2)
		ax.plot([-0.1, 1.0], [-0.2,-0.2], color='gray', linestyle='--', lw=2)
		ax.plot([-0.1, 1.0], [0.4,1.3], color='gray', linestyle='--', lw=2)
		# Power law
		ax.plot([0.1, 0.6], [0.125,0.75], color='black', linestyle='-', lw=2)
		# Revised wedge
		ax.plot([0.08, 0.3471], [0.15,0.15], color='gray', linestyle='-', lw=2)
		ax.plot([0.08, 0.08], [0.15,0.3668], color='gray', linestyle='-', lw=2)
		ax.plot([0.3471, xmax], [0.15,1.21*xmax-0.27], color='gray', linestyle='-', lw=2)
		ax.plot([0.08, xmax], [0.3668,1.21*xmax+0.27], color='gray', linestyle='-', lw=2)
	if do_centroids:
		print(centroids)
		for b in zbins:
			ax = axis[b["ax"][0],b["ax"][1]]
			bx = [bx for (z0, z1,bx, by) in centroids ]
			by = [by for (z0, z1,bx, by) in centroids ]
			bxx = [bx for (z0, z1, bx, by) in centroids if z0 == b["z0"]]
			byy = [by for (z0, z1, bx, by) in centroids if z0 == b["z0"]]
			#ax.plot(bx, by, marker='o', c='red', linewidth=2, markersize=3)
			#ax.plot(bxx[0], byy[0], marker='*', c='red', linewidth=5, markersize=20)
			# Revised wedge
			if b["z0"] > 1.0:
				ox = bxx[0] - bariadj_obsagn_ref1[0]
				oy = byy[0] - bariadj_obsagn_ref1[1]
				print(oy, 0.15+oy)
				left_bar_size = (0.3668-0.15)
				# lower horiz bar
				if oy < left_bar_size:
					ax.plot([0.08, 0.3471+oy/1.21], [0.15+oy,0.15+oy], color='red', linestyle='-', lw=2)
				else:
					ax.plot([0.08+(oy-left_bar_size)/1.21, 0.3471+oy/1.21], [0.15+oy,0.15+oy], color='red', linestyle='-', lw=2)
				# left vert bar
				#ax.plot([0.08, 0.08], [0.15,0.3668], color='red', linestyle='-', lw=2)
				# lower diag bar
				#ax.plot([0.3471+oy/1.21, xmax], [0.15+oy,1.21*xmax-0.27], color='red', linestyle='-', lw=2)
				# upper diag bar
				#ax.plot([0.08, xmax], [0.3668,1.21*xmax+0.27], color='red', linestyle='-', lw=2)
	plt.show()
	selection = [ i for (i,j) in zip(idx,use) if j == '.' ]
	print("selection=%d" % len(selection))
	return [selection, ] 
#end color_color_3x2


def mass_selection(alldata, cndf, kth=None):
	print("mass selection")
	plotit = True
	fields, data = FilterByFields(alldata, "id","lmass","z_mass","z_peak","ch1","ech1","ch2","ech2","ch3","ech3","ch4","ech4","Ks","Ks_tot","USE")
	#fields, data = FilterByFields(alldata, "id","lmass","z_mass","z_peak","ch1","ch2","ch3","ch4","Ks","Ks_tot","USE")
	conds = { 
			"Ks_tot": lambda k: kscut(k, th=kth),
			"lmass":pos,
			"z_peak":pos,
			"USE":is1
			}
	r1 = FilterByConds(fields, data, conds)
	r21 = FilterByFuncs(fields, r1, 
			{
				"ch1":uplog1b,
				"ch2":uplog2b,
				"ch3":uplog3b,
				"ch4":uplog4b,
				"ech1":upsn1,
				"ech2":upsn2,
				"ech3":upsn3,
				"ech4":upsn4,
				}, 
			["Ks","Ks_tot","ch1","ech1","ch2","ech2","ch3","ech3","ch4","ech4"])
	r2 = FilterByFuncs(fields, r21, 
			{
				#"USE":lacy2007_func
				#"USE":donley2012_func
				#"USE":donley2012_bariadj_func
				"USE":donley2012_baricut_func
				},
			["ch1","ch2","ch3","ch4","z_peak"])
	[idx, lmass, zmass, zpeak, c1, e1, c2, e2, c3, e3, c4, e4, ks, kstot, use] = zip(*r2)

	(xmin, ymin) = (0.0, 0.0)
	(xmax, ymax) = (2.5, 2.5)
	#(binmin, binmax, binstep) = (4.5, 12.5, 0.001)
	(binmin, binmax, binstep) = (9.5, 11.5, 0.01)
	binres = []
	binidx = binmin
	while binidx < binmax:
		binres.append(binidx)
		binidx += binstep
	zbins = [ 
			{"z0":0.2, "z1":0.5, "ax":[0,0], "y0":0.69},
			{"z0":0.5, "z1":1.0, "ax":[0,1], "y0":0.59},
			{"z0":1.0, "z1":1.5, "ax":[0,2], "y0":0.49},
			{"z0":1.5, "z1":2.0, "ax":[1,0], "y0":0.49},
			{"z0":2.0, "z1":2.5, "ax":[1,1], "y0":0.49},
			{"z0":2.5, "z1":3.0, "ax":[1,2], "y0":0.49},
			{"z0":3.0, "z1":4.0, "ax":[2,0], "y0":0.49}
			]	
	csmf_u = {}
	csmf_n = {}
	for b in range(0,len(zbins)):
		V = 1.0 + ( (zbins[b]["z0"]+zbins[b]["z1"])/2.0 )**3.0 
		m = [l for (i, l,z) in zip(idx, lmass,zpeak) if z >= zbins[b]["z0"] and z < zbins[b]["z1"]]
		# using trick for truly empirical CDF in
		# https://stackoverflow.com/questions/24575869/read-file-and-plot-cdf-in-python
		#sm = np.sort(m[::-1])
		#y1 = np.arange(len(sm))/float(len(sm))
		#csmf_u[b] = y1[::-1]
		#csmf_n[b] = y1[::-1]/float(V)
		
		# using good'ol histogram + cumsum
		counts, bin_edges = np.histogram(m, bins=binres)
		cdf = np.cumsum(counts[::-1])
		csmf_u[b] = cdf[::-1]
		csmf_n[b] = cdf[::-1]/float(V)

	cndf_m = [8,9,10,11,11.5]
	cndf_n = {}
	cndf_r = {}
	for row in cndf:
		zrow = []
		for col in ["M8","M9","M10","M11","M11.5"]:
			zrow.append(row[col])
		cndf_n[row["zbin"]]=interp1d(cndf_m,zrow, kind='cubic')(binres[0:-1])
		#cndf_n[row["zbin"]]=interp1d(cndf_m,zrow, kind='quadratic')(binres[0:-1])
		cndf_r[row["zbin"]]=[cndf_m, zrow]

	# Choose which version to use of the SMFs
	# normalized cmsf_n 
	# non normalized cmsf_u
	#csmf = csmf_u
	#csmf = csmf_n
	csmf = cndf_n

	def findMz(cz, z):
		for i in range(0,len(binres)-1):
			czz = csmf[z][i]
			if czz < cz:
				#Mz = binmin+(i-1)*binstep
				Mz = binres[i]
				return Mz, czz
		return False # did not converge

	#cdf3 = binres
	#for k in range(0,len(binres)-100,100):
	#	count = len([i for (i,l,z) in zip(idx, lmass, zpeak) if 
	#			z > zbins[b]["z0"] and z <= zbins[b]["z1"] and l >= binres[k] ])
	#	for n in range(0,100):
	#		cdf3[k+n] = count

	Mz3 = 10.86
	cz3 = csmf[5][int((Mz3-binmin)/binstep)]
		
	fig, ax = plt.subplots()
	for c in csmf.keys():
		ax.plot(binres[0:-1],csmf[c],label=("%.1f <= z < %.1f" % (zbins[c]["z0"], zbins[c]["z1"])))
		ax.scatter(cndf_r[c][0],cndf_r[c][1],marker='o',c='white',facecolors='black')
	#ax.set_yscale('log')
	#ax.axis([9.75,11.75,0.1*int(cz3),10*int(cz3)])
	ax.plot([9.75,11.75], [cz3,cz3], color='darkgray', linestyle='--', lw=2)
	#ax.plot([Mz3,Mz3], [0.1*int(cz3),10*int(cz3)], color='darkgray', linestyle='--', lw=2)
	ax.set_xlabel(r'$\log(M_\star\,/\,M_\odot)$',fontsize="large")
	ax.set_ylabel(r'$log\,n(> M_\star)\,[\mathrm{Mpc}^{-3}]$',fontsize="large")
	ax.set_xlim(9.7, 11.7)
	ax.set_ylim(-5.5, -2.9)

	selection = []
	colors = []
	selm = {}
	mvsz = []
	frac = []
	for b in range(0,len(zbins)):
		(Mz, cz) = findMz(cz3, b)
		mvsz.append(((zbins[b]["z0"]+zbins[b]["z1"])/2,Mz))
		cz = csmf[b][int((Mz-binmin)/binstep)]
		czu = csmf_u[b][int((Mz-binmin)/binstep)]
		ax.plot(Mz, cz, 'k^')
		ax.annotate(r'$n=%d$' % czu, xy=(Mz,cz), xytext=(Mz*(1+0.01*((-1)**b)), cz*(1-0.15*((-1)**b))), arrowprops=dict(facecolor='black', arrowstyle='->'),fontsize="large")
		print("[%.1f:%.1f): Mz=%.2f, N=%d" % (zbins[b]["z0"],zbins[b]["z1"], Mz, czu), end="")
		sel = [i for (i,m,z) in zip(idx, lmass, zpeak) if 
				z >= zbins[b]["z0"] and z < zbins[b]["z1"] and m >= Mz ]
		selm[b] = [m for (i,m,z) in zip(idx, lmass, zpeak) if 
				z >= zbins[b]["z0"] and z < zbins[b]["z1"] and m >= Mz ]
		selection += sel
		col = [i for (i,m,z) in zip(use, lmass, zpeak) if 
				z >= zbins[b]["z0"] and z < zbins[b]["z1"] and m >= Mz ]
		colors += col
		agn = len([c for c in col if c == '[IR-AGN]'])
		nagn = len([c for c in col if c == '[not IR-AGN]'])
		print("; Ntot=%d, %s" % (len(sel),Counter(col)))
		frac.append( (zbins[b]["z0"], zbins[b]["z1"], (agn/(nagn+agn)), agn, agn+nagn ) )
		#print("             Ntot=%d, IR-AGN=%d, not IR-AGN=%d" % (len(sel),agn, nagn))
		#print(stats.describe(selm[b]))
	ax.legend(loc=3, ncol=2)
	if plotit: plt.show()


	binmin = 10
	binmax = 12
	binstep = 0.05
	binres = []
	binidx = binmin
	while binidx < binmax:
		binres.append(binidx)
		binidx += binstep

	if plotit:
		fig, ax = plt.subplots()
		bintot,edg = np.histogram(selm[0], bins=binres)
		bintot = 0.0
		for b in range(0,len(zbins)-1):
			Phi,edg = np.histogram(selm[b], bins=binres)
			selm_key = "[%.1f,%.1f]" % (zbins[b]["z0"],zbins[b]["z1"]) 
			ax.plot(binres[0:-1],Phi,label=selm_key) 
			bintot = bintot + Phi
			#V = float(Vphys) * ( (1.0 + (b["z0"]+b["z1"])/2.0)**3.0 )
			#edg[0] = edg[0:-1] + dM/2
			#Phi_n = Phi_u / float(V)
		ax.plot(binres[0:-1],bintot,label="Total") 
		ax.set_yscale('log') 
		ax.set_xlabel(r'$\log(M_\star\,/\,M_\odot)$') 
		#ax.set_ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$') 
		ax.legend(loc=2, ncol=2) 
		plt.show()

	#if plotit:
	#	print("IR AGN Fractions")
	#	(z1,z2,r,agn,tot) = zip(*frac[:-1])
	#	plotwbars(z1,z2, r, agn, tot, (0.5,10.8), "z", "frac_{AGN}")
	#	
	#	fig, ax = plt.subplots()
	#	(z,m) = zip(*mvsz[:-1])
	#	fit = np.polyfit(z, m, deg=1)
	#	line = []
	#	for i in z:
	#		line.append(fit[0]*float(i)+fit[1])
	#	ax.plot(z, line, color='black')
	#	ax.scatter(z, m, marker='o',c='black')
	#	tt = ax.text(0.5, 10.9, r'$M_z = %.2f %.2f\,z$' % (fit[1], fit[0]), size=16)
	#	#ax.annotate(tt, xy=(0.5,10.8))
	#	ax.set_ylabel(r'$\log(M_n)\,(M_\odot)$') 
	#	ax.set_xlabel(r'$z$') 
	#	plt.show()

	#for b in range(0,len(zbins)):
	#	for m in np.arange(10,12,0.5):
	#		s = [i for (i,l,z) in zip(idx, lmass, zpeak) 
	#				if z > zbins[b]["z0"] and z <= zbins[b]["z1"] and l >= m ]
	#		print("[%d] Mz = %.2f, len = %d, cum = %d" % (b, m, len(s),
	#			csmf[b][int((m-binmin)/binstep)] ) )
	return [selection,colors]
# end mass_selection

def plot_fracs(alldata, kth=None, selection=None):
	print("plot_fracs")

	selset1 = None
	if selection is not None:
		selset0 = set(selection[0])
		if len(selection) > 1:
			print(Counter(selection[1]))
			selset1 = True

	def inselection(x): return selection is None or x in selset0
	def selectioncolor(x,i): 
		if selset1 is not None:
			return selection[1][selection[0].index(i)]
		else: 
			return x

	fields, data = FilterByFields(alldata, "id","L153_uv","L155_uv","L155_vj","L161_vj", "z_peak","Ks_tot","star","USE")
	conds = { "id":inselection,
			"Ks_tot": lambda k: kscut(k, th=kth),
			"USE":is1,
			}
	r1 = FilterByConds(fields, data, conds)
	r2 = FilterByFuncs(fields, r1, {"USE":selectioncolor}, ["id"])
	# In input selection both quiescent and not quiescent
	[idx1, flux_u1,flux_v1,flux_v2,flux_j2,zspec1,kstot,star1,use1] = zip(*r2)

	r3 = FilterByConds(fields, r1, 
			{"star": quiescent},
			["L153_uv","L155_uv","L155_vj","L161_vj"])
	r4 = FilterByFuncs(fields, r3, {"USE":selectioncolor}, ["id"])
	# In input selection but only quiescent
	[idx2, flux_u1,flux_v1,flux_v2,flux_j2,zspec2,kstot,star2,use2] = zip(*r4)
	
	print("len(idx1) = %d, len(idx2) = %d"  % (len(idx1), len(idx2)) )
	#zlist = [(0.2,0.5),(0.5,1.0),(1.0,1.5),(1.5,2.0),(2.0,2.5),(2.5,3.0)]
	zlist = [(0.2,0.5),(0.5,1.0),(1.0,1.5),(1.5,2.0),(2.0,2.5),(2.5,3.0),(3.0,4.0)]

	fracs_agn = []
	frac_agn_tot = []
	frac_xagn_tot = []
	frac_iragn_tot = []
	for z1,z2 in zlist:
		if selset1 is not None:
			z1_tot = len([i for (i,z,u) in zip(idx1,zspec1,use1) if z >= z1 and z < z2])
			z1_agn = len([i for (i,z,u) in zip(idx1,zspec1,use1) if z >= z1 and z < z2 and u in mycolors_agn])
			z1_xagn = len([i for (i,z,u) in zip(idx1,zspec1,use1) if z >= z1 and z < z2 and u in mycolors_xagn])
			z1_iragn = len([i for (i,z,u) in zip(idx1,zspec1,use1) if z >= z1 and z < z2 and u in ['[IR-AGN]',]])
			r_agn = ( 0 if z1_agn == 0 else z1_agn/z1_tot )
			r_xagn = ( 0 if z1_xagn == 0 else z1_xagn/z1_tot )
			r_iragn = ( 0 if z1_iragn == 0 else z1_iragn/z1_tot )
			frac_agn_tot.append( (z1, z2, r_agn, z1_agn, z1_tot, mycolors['[AGN]'][0], -0.03, "AGN fraction of total") )
			frac_xagn_tot.append( (z1, z2, r_xagn, z1_xagn, z1_tot, mycolors['[Xray AGN]'][0], 0.03, "X-ray AGN fraction of total") )
			frac_iragn_tot.append( (z1, z2, r_iragn, z1_iragn, z1_tot, mycolors['[IR-AGN]'][0], 0.05, "IR AGN fraction of total") )
	(z1,z2,r,agn,tot,c,offset,e) = zip(*frac_agn_tot)
	fracs_agn.append((z1,z2,r,agn,tot,c,offset,e))
	(z1,z2,r,agn,tot,c,offset,e) = zip(*frac_xagn_tot)
	fracs_agn.append((z1,z2,r,agn,tot,c,offset,e))
	(z1,z2,r,agn,tot,c,offset,e) = zip(*frac_iragn_tot)
	fracs_agn.append((z1,z2,r,agn,tot,c,offset,e))
	print("AGN fractions of total")
	plotsteps(fracs_agn, "z", "frac_{AGN}")

	fracs_agn = []
	frac_agn_tot = []
	frac_agn_qui = []
	for z1,z2 in zlist:
		if selset1 is not None:
			z1tot = len([i for (i,z,u) in zip(idx1,zspec1,use1) if z >= z1 and z < z2])
			z1agnt = len([i for (i,z,u) in zip(idx1,zspec1,use1) if z >= z1 and z < z2 and u in mycolors_agn])
			z1qui = len([i for (i,z,u) in zip(idx2,zspec2,use2) if z >= z1 and z < z2])
			z1agnq = len([i for (i,z,u) in zip(idx2,zspec2,use2) if z >= z1 and z < z2 and u in mycolors_agn])
			r_agn_tot = ( 0 if z1agnt == 0 else z1agnt/z1tot )
			r_agn_qui = ( 0 if z1agnq == 0 else z1agnq/z1qui )
			frac_agn_tot.append( (z1, z2, r_agn_tot, z1agnt, z1tot, mycolors['[Xray AGN]'][0], -0.03, "AGN fraction of total") )
			frac_agn_qui.append( (z1, z2, r_agn_qui, z1agnq, z1qui, mycolors['[not AGN]'][0], 0.03, "AGN fraction of quiescents") )
	(z1,z2,r,agn,tot,c,offset,e) = zip(*frac_agn_tot[:-1])
	fracs_agn.append((z1,z2,r,agn,tot,c,offset,e))
	(z1,z2,r,agn,tot,c,offset,e) = zip(*frac_agn_qui[:-1])
	fracs_agn.append((z1,z2,r,agn,tot,c,offset,e))
	print("AGN fractions")
	plotsteps(fracs_agn, "z", "frac_{AGN}")

	edges_qui = {
			'Quiescent fraction of Not AGNs':
				[mycolors['[not AGN]'][0],	['[not IR-AGN]',],	-0.03],
			'Quiescent fraction of AGNs (IR+X Ray)':
				[mycolors['[Xray AGN]'][0],		mycolors_agn,				0.03]
			}
	fracs_qui = []
	print("zbin		q	nq	tot	q/tot	comment")
	for e,l in edges_qui.items():
		frac_qui = []
		for z1,z2 in zlist:
			if selset1 is not None:
				z1tot = len([i for (i,z,u) in zip(idx1,zspec1,use1) if z >= z1 and z < z2 and u in l[1]])
				z1qui = len([i for (i,z,u) in zip(idx2,zspec2,use2) if z >= z1 and z < z2 and u in l[1]])
				#print("z = [%.1f,%.1f) q=%d, nq=%d q/tot=%f, agn=%d"  % (z1, z2, z1qui, z1tot-z1qui, z1qui/z1tot, z1agn) )
				r_qui = ( 0 if z1qui == 0 else z1qui/z1tot )
				#print("z = [%.1f,%.1f) q=%d, nq=%d, tot=%d, q/tot=%.2f, %s"  % (z1, z2, z1qui, z1tot-z1qui, z1tot, r, e) )
				print("[%.1f,%.1f)	%d	%d	%d	%.2f	%s"  % (z1, z2, z1qui, z1tot-z1qui, z1tot, r_qui, e) )
				frac_qui.append( (z1, z2, r_qui, z1qui, z1tot, l[0], l[2], e) )
		(z1,z2,r,qui,tot,c,offset,e) = zip(*frac_qui[:-1])
		fracs_qui.append((z1,z2,r,qui,tot,c,offset,e))
	print("QUI fractions")
	plotsteps(fracs_qui, "z", "frac_{QUI}")
# end plot_fracs

def uvj_3x2(alldata, kth=None, selection=None):
	print("uvj_3x2")

	selset1 = None
	if selection is not None:
		selset0 = set(selection[0])
		if len(selection) > 1:
			print(Counter(selection[1]))
			selset1 = True

	def inselection(x): return selection is None or x in selset0
	#def selectifcolor(x): return x == 'r' or x == 'y'
	def selectioncolor(x,i): 
		if selset1 is not None:
			return selection[1][selection[0].index(i)]
		else: 
			return x

	fields0, data0 = FilterByFields(alldata, "id","L153_uv","L155_uv","L155_vj","L161_vj", "z_peak","Ks_tot","USE")
	conds0 = { 
			"Ks_tot": lambda k: kscut(k, th=kth),
			"USE":is1
			}
	r0 = FilterByConds(fields0, data0, conds0)
	[idx0, flux_u10,flux_v10,flux_v20,flux_j20,zspec0,kstot0,use0] = zip(*r0)

	fields, data = FilterByFields(alldata, "id","L153_uv","L155_uv","L155_vj","L161_vj", "z_peak","Ks_tot","USE")
	conds = { "id":lambda s: inselection(s),
			"Ks_tot": lambda k: kscut(k, th=kth),
			"USE":is1
			}
	r1 = FilterByConds(fields, data, conds)
	r2 = FilterByFuncs(fields, r1, {"USE":selectioncolor}, ["id"])
	#r3 = FilterByConds(fields, r2, {"USE":selectifcolor})
	[idx, flux_u1,flux_v1,flux_v2,flux_j2,zspec,kstot,use] = zip(*r2)
	print("z = [:] %d" % len(idx))

	fig, axis = plt.subplots(2,3)
	(xmin, ymin) = (0.0, 0.0)
	(xmax, ymax) = (2.5, 2.5)
	binres = (50, 50)
	zbins = [ 
			{"z0":0.2, "z1":0.5, "ax":[0,0], "y0":0.69},
			{"z0":0.5, "z1":1.0, "ax":[0,1], "y0":0.59},
			{"z0":1.0, "z1":1.5, "ax":[0,2], "y0":0.49},
			{"z0":1.5, "z1":2.0, "ax":[1,0], "y0":0.49},
			{"z0":2.0, "z1":2.5, "ax":[1,1], "y0":0.49},
			{"z0":2.5, "z1":3.0, "ax":[1,2], "y0":0.49}
			]	
	for b in zbins:
		x = [ -2.5*np.log10(i/j) for (i,j,z) in zip(flux_v2,flux_j2,zspec) if z > b["z0"] and z <= b["z1"] ]
		y = [ -2.5*np.log10(i/j) for (i,j,z) in zip(flux_u1,flux_v1,zspec) if z > b["z0"] and z <= b["z1"] ]
		u = [ u for (u,z) in zip(use,zspec) if z > b["z0"] and z <= b["z1"] ]
		x0 = [ -2.5*np.log10(i/j) for (i,j,z) in zip(flux_v20,flux_j20,zspec0) if z > b["z0"] and z <= b["z1"] ]
		y0 = [ -2.5*np.log10(i/j) for (i,j,z) in zip(flux_u10,flux_v10,zspec0) if z > b["z0"] and z <= b["z1"] ]
		print("z = [%.1f,%.1f] %d, %d" % (b["z0"],b["z1"],len(x),len(y)))
		ax = axis[b["ax"][0],b["ax"][1]]
		x1 = 1.3/0.88-b["y0"]
		x2 = 0.88*1.6+b["y0"]
		ax.plot([xmin, x1], [1.3,1.3], color='darkgray', linestyle='--', lw=2)
		ax.plot([x1, 1.6], [1.3,x2], color='darkgray', linestyle='--', lw=2)
		ax.plot([1.6, 1.6], [x2,xmax], color='darkgray', linestyle='--', lw=2)
		
		ax.hist2d(x0, y0, bins=binres, cmap=plt.cm.Greys,zorder= 0)
		#ax.scatter(x,y,c=u,marker='.',edgecolors='none', alpha=0.5, zorder = 1)
		#edges = {
		#		'not_agn':	[ ['[not IR-AGN]',], mycolors['[not AGN]'][0], "none", '.'],
		#		'xagn':		[ mycolors_xagn,			mycolors['[Xray AGN]'][0], mycolors['[Xray AGN]'][0], 'o'],
		#		'iragn':	[ ['[IR-AGN]',],		mycolors['[IR-AGN]'][0], mycolors['[IR-AGN]'][0], '^']
		#		}
		#for e,l in edges.items():
		#	xe = [i for (i,j) in zip(x,u) if j in l[0] ]
		#	ye = [i for (i,j) in zip(y,u) if j in l[0] ]
		#	#print(e,l,xe,ye)
		#	ax.scatter(xe, ye, c=l[1], edgecolors=l[2], marker=l[3]) #zorder=1.0)
		#if selset1 is None:
		#	ax.hist2d(x, y, bins=binres, cmap=plt.cm.Greys)
		for k in sorted(mycolors.keys()):
			xe = [i for (i,j) in zip(x,u) if j == k ]
			ye = [i for (i,j) in zip(y,u) if j == k ]
			if len(xe)>0:
				ax.scatter(xe, ye, c=mycolors[k][0], edgecolors=mycolors[k][3], marker=mycolors[k][2], label=mycolors[k][1]) #zorder=1.0)
		if selset1 is None:
			ax.hist2d(x, y, bins=binres, cmap=plt.cm.Greys)
		ax.axis([xmin,xmax,ymin,ymax])
		legend = ax.legend(loc='lower right', fontsize='small', shadow=True)
		ax.set_aspect("equal")
		if selset1 is None:
			ax.set_title(r'$%.1f\leq z<%.1f\,\,N=%d$' % (b["z0"],b["z1"],len(x)))
		ax.set_xlabel(r'$V-J_{rest}$')
		ax.set_ylabel(r'$U-V_{rest}$')
	plt.show()
# end uvj_3x2


def jku(alldata, kth=None):

	fields, data = FilterByFields(alldata, "J","Ks","u","Ks_tot","nan_contam","K_star")
	conds = {
		#"J":pos,
		#"Ks":pos,
		#"u":pos,
		#"nan_contam":is0,
		}
	r1 = FilterByConds(fields, data, {"K_star":is1})
	r2 = FilterByConds(fields, data,
			{"Ks": lambda f, ks, kstot: kscut2(f, ks, kstot, th=kth)},
			["Ks","Ks_tot"]
			)
	r3 = FilterByFuncs(fields, r2, 
			{
				"J":uplogb,
				"Ks":uplogb,
				"u":uplogb,
				}, 
			["Ks","Ks_tot"])
	[flux_j,flux_ks,flux_u,kstot,nancon,star] = zip(*r3)
	print(star[:10])
	r4 = FilterByFuncs(fields, r3, {"K_star":colorjku})
	[flux_j,flux_ks,flux_u,kstot,nancon,star] = zip(*r4)
	print(star[:10])
	
	for co in ["red","blue"]:
		uj = [ i-j for (i,j,c) in zip(flux_u,flux_j,star) if c == co]
		jk = [ i-j for (i,j,c) in zip(flux_j,flux_ks,star) if c == co]
		x = uj
		y = jk
		plt.scatter(x, y, c=co, marker=".", edgecolors="none")
	#plt.scatter(-2.5*np.log10(x), -2.5*np.log10(y), marker=".", edgecolors="none", c=star)
	plt.xlim(xmin=-2,xmax=12)
	plt.ylim(ymin=-1,ymax=1.5)
	plt.xlabel(r'$u^*\,-\,J$')
	plt.ylabel(r'$J\,-\,K_s$')
	plt.plot([-1.38, 3], [-1,-0.21], color='darkgray', linestyle='-', lw=2)
	plt.plot([3, 12], [-0.21,0.51], color='darkgray', linestyle='-', lw=2)
	plt.show()

