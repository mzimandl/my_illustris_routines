print("V Illustrisu 1 chybí snapshoty 53 a 55! Berte to v úvahu při užívání funkcí (je tam -2)")

import img_manip as manip

import matplotlib.pyplot as plt
import matplotlib.colors as colors

import h5py
import numpy as np
import requests

import os.path

baseUrl = 'http://www.illustris-project.org/api/'
headers = {"api-key":"b90b416c1d9ec1456eb15582b4f26dd4"}

PartType = {
		'gas':'PartType0',
		'dm':'PartType1',
		'stars':'PartType4',
		'bhs':'PartType5'	
	}

def get(path, params=None, filename=None):

	r = requests.get(path, params=params, headers=headers)
	r.raise_for_status()

	if r.headers['content-type'] == 'application/json':
		return r.json() # parse json responses automatically

	if 'content-disposition' in r.headers:
		if not filename: filename = r.headers['content-disposition'].split("filename=")[1]
		with open(filename, 'wb') as f:
			f.write(r.content)
		return filename # return the filename string
	return r

r = get(baseUrl)
sim = get(r['simulations'][0]['url'])
snaps = get(sim['snapshots'])

############################################### funkce ####################################################

#-----------------------downloading illustris data
def down_sub(simulace, snimek, subhalo, typ = '', filename=''):
	if typ == "stars":
		cutout_request = {'stars':'Coordinates,Masses'}
	elif typ == "gas":
		cutout_request = {'gas':'Coordinates,Masses'}
	elif typ == "dm":
		cutout_request = {'dm':'Coordinates'}
	elif typ == "gs":
		cutout_request = {'stars':'Coordinates,Masses', 'gas':'Coordinates,Masses'}
	else:
		cutout_request = {'all'}

	num_snapshots = r['simulations'][simulace]['num_snapshots']

	sim = get(r['simulations'][simulace]['url'])
	sub = get(sim['snapshots']+str(snimek)+'/subhalos/'+str(subhalo)+'/')
	
	if len(filename)==0:	
		filename = 'cutoutSub-'+str(snimek)+'-'+str(subhalo)

	info = get(sub['meta']['info'])
	np.save(filename, info)
	soubor = filename + '.hdf5'
	if not os.path.isfile(soubor):
		soubor = get(sub['cutouts']['subhalo'], cutout_request, soubor)
	
	return filename

#-----------------------making 2D histogram, calculating surface brightness
Ill_mass_unit = 10**10 #solar masses
Ill_spatial_unit = 10**3 #parsec

#surface density per parsec squared (Msol/pc**2)
def hist_dens(filename, typ, x_res, y_res, x_range, y_range, projekce):
	if not os.path.isfile(filename+'.npy') or not os.path.isfile(filename+'.hdf5'):
		return False
	else:	
		info = np.load(filename+'.npy').item()
	
		with h5py.File(filename+'.hdf5', 'r') as f:
			x = f[PartType[typ]]['Coordinates'][:,projekce[0]] - info['Subhalo']['SubhaloPos'][projekce[0]]
			y = f[PartType[typ]]['Coordinates'][:,projekce[1]] - info['Subhalo']['SubhaloPos'][projekce[1]]
		
			if typ == 'stars':		
				dens = f[PartType[typ]]['Masses'][:]
				grid, xedges, yedges = np.histogram2d(x, y, bins=(x_res, y_res), range=[x_range, y_range], weights = dens)
				
				#calculating surface density
				np.multiply(grid, (Ill_mass_unit)/(Ill_spatial_unit**2*np.diff(x_range)*np.diff(y_range)/(x_res*y_res)), grid) # density per parsec**2

			elif typ == 'dm':
				grid, xedges, yedges = np.histogram2d(x, y, bins=(x_res, y_res), range=[x_range, y_range])
		return grid, xedges, yedges

MassToLightV = 2 #McGaugh 2016
MassToLightI = 1
MV = 4.8

#magnitude per arcsec
def hist_surf(grid):
	np.divide(grid, MassToLightV, grid) # luminosity per parsec**2
	np.log10(grid, grid)
	np.multiply(grid,2.5,grid)
	np.subtract(MV+21.572,grid,grid) # surface brightness in mag
	#grid[np.isinf(grid)] = float('Inf') # nahradí NaN za Inf

#-----------------------drawing histogram
def save_hist(grid, xedges, yedges, label_x='', label_y='', name='test', log=False, bad='black'):
	plt.clf()
	plt.xlabel(label_x)
	plt.ylabel(label_y)
	
	if log:
		img = plt.imshow(grid, interpolation='none', cmap='Greys', norm=colors.LogNorm(), origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
	else:
		img = plt.imshow(grid, interpolation='none', cmap='Greys', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
	
	img.cmap.set_bad(color=bad) # set NaN values color
	plt.colorbar(img)

	plt.savefig(name + '.png', bbox_inches='tight')
	
	return

def show_hist(grid, xedges, yedges, log=False, bad='black', colormap='Greys'):
	plt.clf()
	
	if log:
		img = plt.imshow(grid, interpolation='none', cmap=colormap, norm=colors.LogNorm(), origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
	else:
		img = plt.imshow(grid, interpolation='none', cmap=colormap, origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
	
	img.cmap.set_bad(color=bad) # set NaN values color
	plt.colorbar(img)

	plt.show()
	
	return
