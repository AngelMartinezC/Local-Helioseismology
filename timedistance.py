#-*- coding:utf-8 -*-

"""
	Program to make the time-distance plot of a given array.

	It is necessary to have the location x0 and y0 of the suspected sunquake,
	as well as two angles (in degrees) in which it is located.
	If a .sav file exists, it is possible to make the distribution of times
	and distance automatically, otherwise it will plot only the time 
	(in every 45 seconds) and the distance in pixels
	
	Package polarTransform is needed. This package can be downloaded from 
		https://github.com/addisonElliott/polarTransform
		-----
			pip install PolarTransform
		-----
	Installing process is detailed on the page
	
	Date created: July 10 2018
	Date last modified: July 14 2018
	
"""

from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sp
import time
import polarTransform
from datetime import datetime
import matplotlib.patches as mpatches



class cubo:
	
	def getdata(self,array,savfile):
		
		if isinstance(array,(str)):
			cube = fits.getdata(array,0)
			if savfile:
				head = sp.readsav(savfile)
				header = head.aiatimes
				return cube, header
			else:
				return cube, None
		elif isinstance(array,(int,float)):
			return "Not a valid array"
		else: 
			return NotImplemented



class td:

	rsun_m = 695.5e6
	rsun_pix = 1884 # Value for HMI
	
	def __init__(self,array, x0=148, y0=140, theta0=192, theta1=260, savfile=None, path=None,\
	   radius=None, rad0=0, images=None, cartesian=None, xticks=5, yticks=6):
	
		self.array = array
		self.images = images
		self.x0 = x0
		self.y0 = y0
		self.theta0 = theta0
		self.theta1 = theta1
		self.radius = radius
		self.rad0 = rad0
		self.cartesian = cartesian
		self.xticks = xticks
		self.yticks = yticks
		self.sav = savfile
		self.path = path


	################################################################################################
	
	
	# Obtain the information of the header as a readable file
	
	def values(self):
		
		if self.sav:
			headi = cubo().getdata(self.array,self.sav)
			head = headi[1]
			frame = int(len(head)/2)
			datef = str(head[frame][0:10])
			date = datef[2:12]
			d = datetime.strptime(date, '%Y-%m-%d')
			data = d.strftime(r"%B %d$^{\mathrm{th}}$, %Y")
			value = []
			for i in range(len(head)):
				value.append(head[i][11:16].decode("utf-8")) 
			return data, value
		else: 
			pass

	################################################################################################

	
	# Function wich makes the time-distance plot
	
	def tdplot(self,colorbar=None, plot=None, save=None, noprint=None, notime=None):
		
		self.colorbar = colorbar
		self.plot = plot
		self.save = save
		self.noprint = noprint
		self.notime = notime
		
		# Call the array
		cube = cubo().getdata(self.array,self.sav)
		s = time.time()
		
		# Convert to radians
		t0 = self.theta0*np.pi/180
		t1 = self.theta1*np.pi/180
		
		# To ensure the minimum radius in the plot
		if self.radius == None:
			mx = cube[0][0].shape[0]/2
			my = cube[0][0].shape[1]/2
			if self.x0<=mx and self.y0<=my:
				rad = min(self.x0,self.y0)
			elif self.y0>my and self.x0<=mx:
				rad = min(self.x0,my*2 -self.y0)
			elif self.x0>mx and self.y0<=my:
				rad = min(self.y0,mx*2 -self.x0)
			elif self.x0>mx and self.y0>my:
				rad = min(mx*2 -self.x0,my*2 -self.y0)
		else: 
			rad = self.radius
		
		if self.noprint == None:
			print("\n Radius in pixels is {}".format(rad))
		else: 
			pass

		# If put, the final image will be up the selected one, 
		#  otherwise it uses the entire range-1
		if self.images == None:
			image = len(cube[0])-1
		else: image = self.images
		
		# td graph
		image_td = []
		for ij in range(image):
			polarImage, ptSettings = polarTransform.convertToPolarImage(cube[0][ij],\
								center=[self.x0, self.y0], initialRadius = self.rad0, \
								finalRadius = rad, initialAngle = t0, finalAngle = t1)
			slices = []
			for i in range(len(polarImage[:,0])):
				slices.append(np.mean(polarImage[i,:]))
			image_td.append(slices)
		e = time.time()
		
		"""-----------------------------------------------------------------------------
		   Variables for the graph. The final graph wil be shown if is
		   selected the keyword "graph", and will be saved if the keyword
		   "save" is selected
		"""
		
		if self.sav:
			
			# Call the function which has information of the flare
			main = self.values()

			# Labels for x range
			xdiv = self.xticks
			xrangex = np.linspace(0,len(image_td[0]),xdiv)
			x2Mm1 = rad*695.5/1884 # Rsun[Mm]=695.5 & Rsun[pixels]=1884
			x2Mm0 = self.rad0*695.5/1884
			xlabel = np.linspace(x2Mm0,x2Mm1,xdiv,dtype=int)
			xlabel = list(map(str,xlabel))

			# Labels for y range
			ydiv = self.yticks
			yrange = np.linspace(0,len(image_td)-1,ydiv)
			ysec = np.array(yrange,dtype=int)
			ylabel = [main[1][i] for i in ysec]
			
			plt.title(main[0])
			plt.xticks(xrangex,xlabel)
			plt.yticks(yrange,ylabel)
			plt.xlabel("Distance (Mm)")
			plt.ylabel("Time [UT]")
			# Image of the variables with all labels defined
	
		else:
			plt.xlabel("Distance [pixels]")
			plt.ylabel("Time [45 seconds]")
		
		plt.imshow(image_td,origin='lower',cmap='Greys_r',interpolation = 'spline36')
		
		
		
		if self.colorbar == None:
			plt.colorbar(label='Velocity (m/s)')
		else: 
			pass
		
		"""
		   Final of the graph 
		   -----------------------------------------------------------------------------
		"""
		
		# Save the graph
		if self.save == None: 
			pass
		else:
			if self.path == None:
				plt.savefig("Sunquake_"+str(self.x0)+"_"+str(self.y0)+"_"+ \
				str(self.theta0)+"_"+str(self.theta1)+"_"+str(rad)+".png",format="png",dpi=450)
				plt.savefig("Sunquake_"+str(self.x0)+"_"+str(self.y0)+"_"+ \
				str(self.theta0)+"_"+str(self.theta1)+"_"+str(rad)+".eps",format="eps",dpi=450)
			else:	
				plt.savefig(self.path+"Sunquake_"+str(self.x0)+"_"+str(self.y0)+"_"+ \
				str(self.theta0)+"_"+str(self.theta1)+"_"+str(rad)+".png",format="png",dpi=450)
				plt.savefig(self.path+"Sunquake_"+str(self.x0)+"_"+str(self.y0)+"_"+ \
				str(self.theta0)+"_"+str(self.theta1)+"_"+str(rad)+".eps",format="eps",dpi=450)
		
		# Show the time spent making the graph
		if self.notime == None:
			print(" The time spent is {:0.4f} s".format(e-s))
		else: 
			pass
		
		# Show the graph
		if self.plot == None:
			plt.show()
		else: 
			pass

		# Return a cartesian version of the set 
		if self.cartesian == None:
			return image_td
		else:
			cartesianImage = ptSettings.convertToCartesianImage(polarImage)
			return cartesianImage
			

		return image_td
	
	################################################################################################


	# Makes a test of many pixels around (5 columns and 4 rows)
	# Reference to IDL procedure (Martinez-Oliveros)
	
	def test(self,columns=5, rows=4):
		
		s = time.time()
		self.x = columns
		self.y = rows
		fig = plt.figure(figsize=(38, 38))
		n=0
		print("\n")
		for j in range(self.y,0,-1):
			for i in range(0,self.x):
				n+=1
				plt.subplot(self.y,self.x,n)
				print("   Plotting image {} of {}".format(n,self.x*self.y))
				Y0 = self.y0-int(self.y/2)+j
				X0 = self.x0-int(self.x/2)+i
				rad = self.radius
				
				image = td(self.array, x0=X0, y0=Y0, theta0=self.theta0, theta1=self.theta1, \
				radius=rad, savfile=self.sav, path=self.path,\
				rad0=self.rad0)
				
				final = image.tdplot(plot=False, colorbar=False,noprint=False,notime=False)
				plt.title(" ")
				plt.xlabel(" ")
				plt.ylabel(" ")
				plt.imshow(final,cmap='Greys_r',origin='lower',interpolation='spline36')
				red_patch = mpatches.Patch(color='red', label="x="+str(X0)+"\ny="+str(Y0))
				plt.legend(handles=[red_patch],bbox_to_anchor=(1, .33))
		e = time.time()
		print("\n Time spent in the test {:0.4f} s".format(e-s))	
		plt.show()		
		
	################################################################################################
	
	
	def slider(self):
		
		from matplotlib.widgets import Slider, Button, RadioButtons

		# Image configuration and first call
		fig, ax = plt.subplots()
		fig.subplots_adjust(bottom=0.15)
		
		def data(x0,y0,t0,t1):
			image = td(self.array, x0=int(x0), y0=int(y0), theta0=int(t0), theta1=int(t1), \
				radius=self.radius, savfile=self.sav, path=self.path, rad0=int(self.rad0), images = self.images)
			im = image.tdplot(plot=False,colorbar=False)
			return im
		
		im = data(self.x0,self.y0,self.theta0,self.theta1)
		image = ax.imshow(im, cmap='Greys_r', interpolation='spline36',origin='lower')

		# Remove actual image to not be shown
		plt.show(block=False)
		plt.pause(0.0001)
		plt.close('all')

		# Set up sliders
		axt1 = plt.axes([0.1, 0.06, 0.34, 0.04], facecolor='lightcyan')
		axt0 = plt.axes([0.1, 0.01, 0.34, 0.04], facecolor='pink')
		axy0 = plt.axes([0.55, 0.06, 0.35, 0.04], facecolor='lightcyan')
		axx0 = plt.axes([0.55, 0.01, 0.35, 0.04], facecolor='pink')
		axrd = plt.axes([0.1,0.85,0.34,0.04])
		
		slider_ra = Slider(axrd,"radio",0,20,valinit = self.rad0)
		slider_t1 = Slider(axt1,r'$\theta _1(^\circ)$',0, 359, valinit = self.theta1)
		slider_t0 = Slider(axt0,r'$\theta _0(^\circ)$',0, 359, valinit = self.theta0)
		slider_x0 = Slider(axy0,r'x',0,200, valinit = self.x0) 
		slider_y0 = Slider(axx0,r'y',0,200, valinit = self.y0)

		plt.axes([0.25, 0.25, 0.6, 0.6])

		# Update the image 
		def update(val):
			theta0 = slider_t0.val
			theta1 = slider_t1.val
			x0 = slider_x0.val
			y0 = slider_y0.val
			r0 = slider_ra.val
			image.set_data(data(int(x0),int(y0),int(theta0),int(theta1)))
			fig.canvas.draw()
			ax.clear()
			plt.draw()
		
		# Change sliders
		slider_t0.on_changed(update)
		slider_t1.on_changed(update)
		slider_x0.on_changed(update)
		slider_y0.on_changed(update)
		slider_ra.on_changed(update)

		plt.show()
	
	################################################################################################

	"""
	# Visualization of the time distance plot cahnging (not necessary)
	def cbar_slider(self):
		
		from matplotlib.widgets import Slider, Button, RadioButtons
		
		fig = plt.figure()
		ax = fig.add_subplot(111)
		fig.subplots_adjust(left=0.25, bottom=0.25)
		
		image = td(self.array, x0=self.x0, y0=self.y0, theta0=self.theta0, theta1=self.theta1,\
		radius = self.radius, savfile=self.sav, path=self.path, rad0=self.rad0)
		im = np.array(image.tdplot(plot=False, colorbar=False, noprint=False)) 
		
		im1 = ax.imshow(im,cmap='Greys_r',origin='lower')
		fig.colorbar(im1)
		axcolor = 'gainsboro'
		axmin = fig.add_axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
		axmax  = fig.add_axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
		
		min0 = im.min()
		max0 = im.max()
		smin = Slider(axmin, 'Min', im.min()*1.5, 0, valinit=min0, color='dimgray')
		smax = Slider(axmax, 'Max', 0, im.max()*1.5, valinit=max0, color='dimgray')
		
		def update(val):
			im1.set_clim([smin.val,smax.val])
			fig.canvas.draw()
		smin.on_changed(update)
		smax.on_changed(update)
		
		plt.show()
	"""
	################################################################################################

# testing branch


if __name__ == "__main__":

	# Paths where the files are located
	# Menawhile this path will be set
	path = "/home/angel/IDLWorkspace/Python/Codes/Examples/"
	savfile = path+"HMIDoppler.difference_coord.sav"
	flare = path+"HMIDoppler.difference.fits"

	# Values (parameters) which are set by default
	image = td(flare,savfile=savfile,rad0=15,path=path,radius=140) 
	
	# Plot the td image
	final = image.tdplot()
	
	###plt.imshow(final)
	###plt.show()
	
	# Test to make many plots around the given center (x0, y0)
	#testing = image.test(rows=3,columns=3)
	
	# Just for change the visualization
	# Uncomment cbar_slider function in td class
	#colorbar = image.cbar_slider()
	
	# Slider across the angles and pixels
	# slide = image.slider()

