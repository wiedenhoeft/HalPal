from __future__ import division 
from colormath.color_objects import LabColor, sRGBColor
from colormath.color_conversions import convert_color
from colormath.color_diff import delta_e_cie2000
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
import scipy.stats

def HaltonSequence(base, index=0):
	"""Halton sequence generator for a given base (must be prime)."""
	while True:
		i = index
		result = 0
		f = 1
		while i>0:
			f = f/base
			result += f *(i % base);
			i = np.floor(i / base)
		yield result
		index += 1
		

def PaletteGenerator(
	maxNrColors = None, 
	L_min = 50, 
	L_max = 95, 
	a_min = -60, 
	a_max = 80, 
	b_min=-60, 
	b_max=80, 
	L_base=11, 
	a_base=2, 
	b_base=3,
	L_index = 0,
	a_index = 0,
	b_index = 0,
	minDistance=0,	# color distance must be greater than this
	maxFails=1000,
	gradient=False):
	
	"""Generate a color palette using Halton sequences in CIE L*a*b* space."""
	
	
	# traverse L*a*b* color space using a low-discrepancy sequence (Halton sequence) for each dimension, and reject a color if it is outside the RGB gamut
	allhex=[]
	alllab = []
	nrFails = 0
	HL = HaltonSequence(L_base, L_index)	# the first argument controls the number of lighness steps
	Ha = HaltonSequence(a_base, a_index)
	Hb = HaltonSequence(b_base, b_index)
	if gradient:
		z = HL.next()
	while True:
		x = Ha.next()
		y = Hb.next()
		if not gradient:
			z = HL.next()
		L= z*(L_max-L_min)+L_min
		a = x*(a_max-a_min)+a_min	
		b = y*(b_max-b_min)+b_min
		labcolor = LabColor(L, a, b)
		rgbcolor = convert_color(labcolor, sRGBColor)
		rgbhex = rgbcolor.get_rgb_hex()
		valid=True
		
		# check if within gamut
		for v in rgbcolor.get_upscaled_value_tuple():
			if v <=0 or v>=255:	# colormath keeps values out of gamut; it does not use negative values though, so any 00 is potentialy out of gamut as well
				valid = False
				break
		
		if valid:
			if rgbhex in allhex:
				valid=False
		
		# check if too close to a color in the palette
		if valid and minDistance > 0:
			for al in alllab:	# TODO fast spatial query structure (kd-tree etc., or binning)
				colorDist = delta_e_cie2000(al, labcolor)
				if colorDist < minDistance:
					valid = False
					break
		if valid:
			allhex.append(rgbhex)
			alllab.append(labcolor)
			if gradient:
				z = HL.next()	# change the lightness only if we found a valid color, which leads to higher contrast in adjacent colors
			nrFails = 0
			yield rgbhex
		else:
			nrFails += 1
			if nrFails >= maxFails:
				print "[ERROR] Could not find a new color after %d iterations!" % maxFails
				break
		if len(alllab) == maxNrColors:
			break





class Palette:
	"""This class implements an \"infinite\" color palette, i.e. when querying p[i] new colors will be created until this query can be satisfied, the palette then contains len(p) = i+1 colors. Notice that this means the length potentially changes with each query! This class is usefull if you have plotting routine that cannot predict how many colors you are going to need."""
		
	def __init__(
		self, 
		nrColors = None, 
		L_min = 40, 
		L_max = 95, 
		a_min = -60, 
		a_max = 80, 
		b_min=-60, 
		b_max=80, 
		L_base=11, 
		a_base=2, 
		b_base=3,
		L_index = 0,
		a_index = 0,
		b_index = 0,
		minDistance = 0, 
		maxFails=1000,
		gradient=False):
		
		self.generator = PaletteGenerator(None, L_min, L_max, a_min, a_max, b_min, b_max, L_base, a_base, b_base, L_index, a_index, b_index, minDistance, maxFails, gradient)
		self.allhex = []
		self.extendTo(nrColors)
		
	
	def extendTo(self, nrColors):
		"""Extends the palette to <nrColors> colors."""
		while len(self.allhex) < nrColors:
			self.allhex.append(self.generator.next())
	
	
	def extendBy(self, nrColors):
		"""Extends the palette by <nrColors> colors."""
		for i in xrange(nrColors):
			self.allhex.append(self.generator.next())
	
	
	def __getitem__(self, i, j=True):
		"""self[i] returns the i-th color if it exists, and an IndexError if we are out of bounds. self[[i]] always returns a color for i>=0, as it extends the palette to the necessary number of colors (i+1)."""
		if type(i)==list:
			self.extendTo(i[0]+1)
			return self.allhex[i[0]]
		return self.allhex[i]


	def __len__(self):
		return len(self.allhex)


	def plotPalette(self, filename, firstcolor=0, lastcolor=None):
		
		allhex=self.allhex[firstcolor:lastcolor]
		if len(allhex)>0:
			# plot the resulting palette for reference
			rowstretch=2
			nrRows = int(np.ceil(np.sqrt(len(allhex))))
			nrCols = int(np.ceil(np.sqrt(len(allhex))))
			fig = plt.figure(figsize=(nrCols, rowstretch*nrRows))
			ax = fig.add_subplot(111, aspect='equal')
			ax.set_xticks([])
			ax.set_yticks([])
			i=0
			for r in reversed(xrange(nrRows)):
				for c in xrange(nrCols):
					ax.add_patch(patches.Rectangle((c, rowstretch*r+0.5), 1, 	1, facecolor=allhex[i]))
					ax.annotate(str(str(i+firstcolor)), xy=(c+0.5, rowstretch*r), xycoords='data', va="bottom", ha="center", fontsize=24)
					i += 1
					if i==len(allhex):
						break
				if i==len(allhex):
					break
			plt.xlim([0, nrCols])
			plt.ylim([r, rowstretch*nrRows])
			plt.tight_layout()
			fig.savefig(filename, dpi=90, bbox_inches="tight")
			plt.close()
		else:
			print "Palette empty, nothing to plot!"


	def savePalette(self, filename, firstcolor=0, lastcolor=None):
		# save the RGB values to a text file
		f = file(filename, "w")
		f.write("\n".join(self.allhex[firstcolor:lastcolor]))
		f.close()


	def rgblist(self, firstcolor=0, lastcolor=None):
		"""Returns a simple list of colors as RGB hex strings."""
		return self.allhex[firstcolor:lastcolor]
	
	

if __name__=="__main__":

	pal = Palette(20, minDistance=0, maxFails=2000)
	#pal = Palette(200, minDistance=0, maxFails=2000, L_min=20, L_max=100, a_min=-128, a_max=128, b_min=-128, b_max = 128)
	#pal = Palette(10, L_index=59, a_index=59, b_index=59)
	#pal = Palette(10, L_base=2, a_base=11, b_base=13)
	print len(pal)
	pal.extendTo(200)
	print len(pal)
	pal.extendBy(3)
	print len(pal)
	pal.extendTo(10)
	print len(pal)
	print pal[200]	# print the 200-th color, create new colors if necessary
	print pal[[300]]	# print the 300-th color, create new colors if necessary
	print len(pal)
	pal.plotPalette("halpal.pdf")
	pal.savePalette("palette.txt")
	
