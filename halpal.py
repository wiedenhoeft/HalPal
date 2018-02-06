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


# Get ranges for LAB values from a set of RGB hexes
def getRanges(hexlist):
	labcolors = [convert_color(sRGBColor.new_from_rgb_hex(h), LabColor) for h in hexlist]
	labtuples = np.array([t.get_value_tuple() for t in  labcolors])
	maxima =  labtuples.max(axis=0)
	minima =  labtuples.min(axis=0)
	L = len(labcolors)
	minDist=np.inf
	for i in  xrange(L-1):
		for j in xrange(i+1, L):
			minDist = min(minDist, delta_e_cie2000(labcolors[i], labcolors[j]))
	return [minima, maxima, minDist]




def PaletteGenerator(
	maxNrColors = None, 
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
	minDistance=0,	# color distance must be greater than this
	adaptToKeepInput=False,	# adapt the parameters so as to keep all input colors in RGBhexes
	maxFails=1000,
	RGBhexes=[],	# a set of colors in RGB hex format. These are the first emitted from the generator, excluding those that violate the minDistance constraint ith respect to the previously generated ones
	gradient=False
	):
	
	if adaptToKeepInput:
		minLab, maxLab, mind = getRanges(RGBhexes)
		L_min = min(minLab[0], L_min)
		a_min = min(minLab[1], a_min)
		a_min = min(minLab[1], b_min)
		L_max = max(maxLab[0], L_max)
		a_max = max(maxLab[1], a_max)
		a_max = max(maxLab[1], b_max)
		minDistance = min(minDistance, mind)

	
	print "L\t%f\t%f" % (L_min, L_max)
	print "a\t%f\t%f" % (a_min, a_max)
	print "b\t%f\t%f" % (b_min, b_max)
	print "d\t%f" % minDistance
	
	"""Generate a color palette using Halton sequences in CIE L*a*b* space."""
	assert minDistance>0, "Minimum distance must be greater than 0!"
	
	alllab = []	# all previously yielded colors in Lab format
	nrFails = 0
	
	# traverse L*a*b* color space using a low-discrepancy sequence (Halton sequence) for each dimension, and reject a color if it is outside the RGB gamut
	HL = HaltonSequence(L_base, L_index)	# the first argument controls the number of lighness steps
	Ha = HaltonSequence(a_base, a_index)
	Hb = HaltonSequence(b_base, b_index)

	i = 0
	while True:
		
		# Yield RGBhexes first, then the ones generated by Halton sequences
		processingInput = (i < len(RGBhexes))
		if processingInput:
			rgbhex = RGBhexes[i]
		else:
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
		
		# check if RGB is within gamut
		if not processingInput:	# manually input colors are always within RGB gamut
			for v in rgbcolor.get_upscaled_value_tuple():
				if v <=0 or v>255:	# colormath keeps values out of gamut; it does not use negative values though, so any 00 is potentialy out of gamut as well
					valid = False
					break

		
		# check if too close to a color in the palette
		nearestDist = np.inf
		if valid:
			
			# Round the color to RGB integer precision, otherwise the minimum distance might be violated when converting RGB hexes to Lab
			# NOTE This is not done earlier since conversions to hex fail if RGB is out of gamut
			rgbcolor = sRGBColor.new_from_rgb_hex(rgbhex)
			labcolor = convert_color(rgbcolor, LabColor)
			
			# check if the minimum distance is violated
			for al in alllab:	# TODO fast spatial query structure (kd-tree etc., or binning)
				colorDist = delta_e_cie2000(al, labcolor)
				nearestDist = min(nearestDist, colorDist)
				if colorDist < minDistance:
					valid = False
					break
		
		if valid:
			alllab.append(labcolor)
			if gradient:
				z = HL.next()	# change the lightness only if we found a valid color, which leads to higher contrast in adjacent colors
			print "%s\tFound color %d of minimum distance %f after %d iteration%s" % (rgbhex, i,  nearestDist, nrFails+1, min(nrFails, 1)*"s")
			nrFails = 0
			i += 1
			yield rgbhex
		else:
			if processingInput:
				print "       \tDropping input color %s at distance %f" % (rgbhex, nearestDist)
			nrFails += 1
			if nrFails >= maxFails:
				print "       \t[ERROR] Could not find a new color after %d iterations!" % maxFails
				yield None
				break
		if len(alllab) == maxNrColors:
			yield None
			break
		




class Palette:
	"""This class implements an \"infinite\" color palette, i.e. when querying p[i] new colors will be created until this query can be satisfied, the palette then contains len(p) = i+1 colors. Notice that this means the length potentially changes with each query! This class is useful if you have plotting routine that cannot predict how many colors you are going to need."""
		
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
		adaptToKeepInput=False,
		maxFails=1000,
		RGBhexes = [],	# a set of colors in RGB hex format. These are the first emitted from the
		gradient=False):
		
		self.generator = PaletteGenerator(None, L_min, L_max, a_min, a_max, b_min, b_max, L_base, a_base, b_base, L_index, a_index, b_index, minDistance, adaptToKeepInput,  maxFails, RGBhexes, gradient)
		self.allhex = []
		self.extendTo(nrColors)
		
	
	def extendTo(self, nrColors):
		"""Extends the palette to <nrColors> colors."""
		while len(self.allhex) < nrColors:
			color = self.generator.next()
			if color is not None:
				self.allhex.append(color)
			else:
				break
	
	
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
					ax.annotate(allhex[firstcolor+i], xy=(c+0.5, rowstretch*r+1.5), xycoords='data', va="bottom", ha="center", fontsize=12, family="monospace")
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

	paired12 = [
		#ColorBrewer's Paired12, see http://colorbrewer2.org/
		"#a6cee3",
		"#1f78b4",
		"#b2df8a",
		"#33a02c",
		"#fb9a99",
		"#e31a1c",
		"#fdbf6f",
		"#ff7f00",
		"#cab2d6",
		"#6a3d9a",
		"#ffff99",
		"#b15928"
			]
	
	
	pal = Palette(
		56, 
		L_base=7,
		minDistance=20, 
		adaptToKeepInput=True,  
		maxFails=500000, 
		RGBhexes = paired12,
		L_min=90, 
		L_max=95, 
		a_min=-128, 
		a_max=128, 
		b_min=-128, 
		b_max = 128
		)
	
	pal.plotPalette("palette.pdf")
	pal.savePalette("palette.txt")
	
