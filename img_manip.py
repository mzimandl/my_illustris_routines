from scipy.ndimage.filters import gaussian_filter
from scipy.signal import medfilt # MEDIAN FILTER

import numpy as np

def MedianFilter(grid, rozsah):
	grid = medfilt(grid, rozsah)
	return grid

def MinMasking(grid, rozsah):
	t_grid = grid
	for i in range(rozsah,len(grid[0])-rozsah):
		for j in range(rozsah,len(grid)-rozsah):
			temp = t_grid[i-rozsah:i+rozsah,j-rozsah:j+rozsah]				
			lmin = np.amin(temp)
			grid[i][j] -= lmin
	return

def UnsharpMasking(grid, sigma, norm=1):
	t_grid = gaussian_filter(grid, sigma)
	
	orig = 0
	temp = 0

	for i in range(0,len(t_grid[0])):
		for j in range(0,len(t_grid)):
			orig += grid[i][j]
			temp += t_grid[i][j]

	np.multiply(t_grid, norm*orig/temp, t_grid)
	np.subtract(grid, t_grid, grid)

	for i in range(0,len(grid[0])):
		for j in range(0,len(grid)):
			if grid[i][j] < 0: grid[i][j]=0;
	return

def MaxContrast(grid):
	lmin = np.amin(grid)
	np.subtract(grid,lmin,grid)

	lmax = np.amax(grid)
	np.divide(grid,lmax,grid)
	return

def SelectiveContrast(grid, rozsah):
	t_grid = grid
	for i in range(rozsah,len(grid[0])-rozsah):
		for j in range(rozsah,len(grid)-rozsah):
			temp = t_grid[i-rozsah:i+rozsah,j-rozsah:j+rozsah]				
			lmin = np.amin(temp)
			for k in range(0, len(temp[0])):
				for l in range(0, len(temp)):
					temp[k][l]-=lmin
			lmax = np.amax(temp)
			grid[i][j] = (grid[i][j]-lmin)/lmax
	return

def Gamma(grid):
	for i in range(0,len(grid[0])):
		for j in range(0,len(grid)):
			#H[i][j] = np.log10(H[i][j])
			grid[i][j] = np.power(grid[i][j],1/6)
