#!/usr/bin/env python

'''
    TimeDistance
    ============

    Program to make the time-distance plot of a given array.

    It calculates angular averages from a given x and y location in pixels
    between two angles (in degrees) along a radius whose size is also defined
    as an argument.

    The package polarTransform is needed. This package converts a cartesian
    image into its polar representation. This module can be downloaded from
    https://github.com/addisonElliott/polarTransform

        >>> pip install PolarTransform

    Installing process is detailed on the page.

    Date created: July 10 2018
    Date last modified: Mar 2023
'''


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
# import scipy.io as sp
import time
import polarTransform
# from datetime import datetime
# import matplotlib.patches as mpatches


class cubo:
    def getdata(self, array):
        if isinstance(array, (str)):
            cube = fits.getdata(array, 0)
            return cube, None
        elif isinstance(array, (int, float)):
            return "Not a valid array"
        else:
            return NotImplemented


class td:

    rsun_m = 695.5e6
    rsun_pix = 1884  # Value for HMI

    def __init__(self, array, x0=148, y0=140, theta0=192, theta1=260,
                 path=None, radius=None, radius0=0, time1=False,
                 time0=False, cartesian=None, xticks=10, yticks=6):
        '''
            timedistance.td
            ===============

            Python class that calls the initial parameters for the TD plot.

            The class functions are tdplot, test, slider, and cbar_slider.

                tdplot:
                =======

                Computes the TD diagram provided the datacube in fits format.

                Example: A datacube `DOPP.fits` has a suspected sunquake in the
                coordinates x:520 y:210 between two arcsections th0:45deg and
                th1:72deg, extending out to 150pix in the image. To calculate
                the TD associated the procedure is:

                    >>> from timedistance import td
                    >>> data = td('DOPP.fits', x0=520, y0=210, \
                    >>>     theta0=45, theta1=72, radius=150)
                    >>> data.plot()

                Additional arguments can be parsed to plot such as
                norm=True      - To normalize intensity images,
                intensity=True - To show intensity images
                magneto=True   - Shows colorbar for magnetograms
                colorbar=True  - To show colorbar
        '''

        self.array = array
        self.time1 = time1
        self.time0 = time0
        self.x0 = x0
        self.y0 = y0
        self.theta0 = theta0
        self.theta1 = theta1
        self.radius = radius
        self.radius0 = radius0
        self.cartesian = cartesian
        self.xticks = xticks
        self.yticks = yticks
        self.path = path

    # Function that makes the time-distance plot

    def tdplot(self, colorbar=False, params=True, norm=False, cmap='gray',
               interpolation='sinc', **kwargs):

        self.colorbar = colorbar
        # self.plot = plot
        self.norm = norm
        # self.save = save
        # self.notime = notime
        # self.title = title
        self.params = params
        # self.magneto = magneto
        # self.intensity = intensity
        self.cmap = cmap

        # Call the array
        cube = cubo().getdata(self.array)
        s = time.time()

        # Convert the angles to radians
        t0 = self.theta0*np.pi/180
        t1 = self.theta1*np.pi/180

        # To ensure the minimum radius in the plot
        # This is to set finalRadius in converting to polar coordinates
        if self.radius is None:
            mx = cube[0][0].shape[0]/2
            my = cube[0][0].shape[1]/2
            if (self.x0 <= mx) and (self.y0 <= my):
                rad = min(self.x0, self.y0)
            elif (self.y0 > my) and (self.x0 <= mx):
                rad = min(self.x0, my*2 - self.y0)
            elif (self.x0 > mx) and (self.y0 <= my):
                rad = min(self.y0, mx*2 - self.x0)
            elif (self.x0 > mx) and (self.y0 > my):
                rad = min(mx*2 - self.x0, my*2 - self.y0)
        else:
            rad = self.radius

        # If put, the final image will be up the selected one,
        # otherwise it uses the entire range-1
        if self.time0:
            image0 = self.time0
        else:
            image0 = 0

        if self.time1:
            image1 = self.time1
        else:
            image1 = len(cube[0])-1

        # td graph
        image_td = []
        for ij in range(image0, image1):
            polarImage, ptSettings = polarTransform.convertToPolarImage(
                cube[0][ij], center=[self.x0, self.y0],
                initialRadius=self.radius0, finalRadius=rad,
                initialAngle=t0, finalAngle=t1)
            slices = []
            for i in range(len(polarImage[0, :])):
                slices.append(np.mean(polarImage[:, i]))
            image_td.append(slices)
        image_td = np.array(image_td)
        e = time.time()

        # Variables for the graph. The final graph will be shown if is
        # selected the keyword "graph", and will be saved if the keyword
        # "save" is selected
        if params:
            plt.rcParams.update({'font.size': 12})
            plt.subplots_adjust(right=0.96, bottom=0.12, top=0.93, left=0.11)
        else:
            pass

        plt.xlabel("Distance (Mm)")
        plt.ylabel("Time (min)")
        # Every 4 pixels in PT correstpond to 1 pix in Postel
        rad_Pos = rad*335/240  # Here I convert from postel pix to PT
        pix_PT = np.arange(0, rad_Pos, self.xticks*4)
        pix_Mm = np.arange(0, rad_Pos/4, self.xticks, dtype=int)
        plt.xticks(pix_PT, pix_Mm)
        plt.yticks(np.arange(0, 240, 40/3), np.arange(0, 180, 10))
        plt.title(self.title)

        if norm:
            image_td /= np.max(image_td)

        plt.imshow(image_td, origin='lower', interpolation=interpolation,
                   aspect='auto', cmap=self.cmap, **kwargs)
        plt.tick_params(direction='out', length=6, width=1.0, colors='k',
                        grid_color='yellow', grid_alpha=0.99)

        # if self.colorbar:
        #     if self.intensity:
        #         plt.colorbar(label=r'$\Delta$ Counts', pad=0.03, aspect=18)
        #     elif self.magneto:
        #         plt.colorbar(label=r'$\Delta B_{LOS}$ (Gauss)', pad=0.03, \
        #                 aspect=18) #, **kwargs)
        #     else:
        #         plt.colorbar(label=r'$\Delta V_{LOS}$ (m/s)', pad=0.03, \
        #                 aspect=18) #, **kwargs)
        # else:
        #     pass

        # ----------------------------------------------------------------
        # -- Final of the graph

        # Save the graph
        # if self.save:
        #     if self.path is None:
        #         plt.savefig("Sunquake_x"+str(self.x0)+"_y"+str(self.y0)+\
        #             "-"+ str(self.theta0)+"_"+str(self.theta1)+\
        #             "deg.png",format="png") #,dpi=350)
        #         #plt.savefig("Sunquake_x"+str(self.x0)+"_y"+str(self.y0)+\
        #         #    "-"+str(self.theta0)+"_"+str(self.theta1)+\
        #         #    "deg.pdf",format="pdf") #,dpi=350)
        #     else:
        #         plt.savefig(self.path+"Sunquake_x"+str(self.x0)+"_y"+\
        #             str(self.y0)+"-"+str(self.theta0)+"_"+\
        #             str(self.theta1)+"deg.png",format="png")
        #         plt.savefig(self.path+"Sunquake_x"+str(self.x0)+"_y"+\
        #             str(self.y0)+"-"+str(self.theta0)+"_"+\
        #             str(self.theta1)+"deg.pdf",format="pdf")
        # else:
        #     pass

        # Show the time spent making the graph
        # if self.notime == None:
        print(" Time spent: {:0.2f} s".format(e-s))

        # Show the graph
        # if self.plot:
        #     plt.show()
        # else:
        #     pass

        # Return a cartesian version of the set
        if self.cartesian is None:
            pass
        else:
            cartesianImage = ptSettings.convertToCartesianImage(polarImage)
            return cartesianImage

        # At the end the function returns the array of the time-distance plot
        return image_td

    """
        The next functions ar not required to show the td plot.

        Functions:
        ----------

        -- slider: --
            Useful function to slide between pixels and angles without need
            to compile by hand new values of pixels and angles.

        -- cbar_slider: --
            To handle color contrast (min and max of colorbar)

    """

    # Makes a test of many pixels around (5 columns and 4 rows)
    # Reference to IDL procedure (Martinez-Oliveros)

    def test(self, columns=5, rows=4):
        '''
        Given the number of columns and rows, it plot as many time distance
        plots as are in columns*rows.

        Params:
        ------
          - columns `int`: Number of columns in the plot grid
          - rows `int`: Number of rows
        '''

        s = time.time()
        self.x = columns
        self.y = rows
        fig = plt.figure(figsize=(18, 12))
        fig.subplots_adjust(bottom=0.07, left=0.10, right=0.96, top=0.95,
                            wspace=0.05, hspace=0.05)
        n = 0
        print("\n")
        for j in range(self.y, 0, -1):
            for i in range(0, self.x):
                n += 1
                ax = fig.add_subplot(self.y, self.x, n)
                print("   Plotting image {} of {}".format(n, self.x*self.y))
                Y0 = self.y0-int(self.y/2)+j
                X0 = self.x0-int(self.x/2)+i
                rad = self.radius

                image = td(self.array, x0=X0, y0=Y0, theta0=self.theta0,
                           theta1=self.theta1, time0=self.time0,
                           time1=self.time1, radius=rad, path=self.path,
                           radius0=self.radius0)

                final = image.tdplot(plot=False, colorbar=False, notime=False)
                plt.title(" ")
                plt.xlabel(" ")
                plt.ylabel(" ")
                plt.imshow(final, cmap='Greys_r', origin='lower',
                           interpolation='spline36', aspect='auto')
                ax.text(0.45, 0.1, "x="+str(X0)+"\ny="+str(Y0),
                        bbox=dict(boxstyle="round", ec='k', fc=(1, 1, 1, 0.6)),
                        transform=ax.transAxes)
        e = time.time()
        print("\n Time spent in test: {:0.2f} s".format(e-s))
        plt.show()

    def slider(self, **kwargs):

        from matplotlib.widgets import Slider  # , Button, RadioButtons

        # Image configuration and first call
        # plt.figure(figsize=(13,9))
        # fig, ax = plt.subplots()
        # fig.subplots_adjust(bottom=0.15)

        def data(x0, y0, t0, t1, radius0):
            image = td(self.array, x0=int(x0), y0=int(y0), theta0=int(t0),
                       theta1=int(t1), radius0=int(radius0),
                       radius=self.radius, path=self.path,
                       time0=self.time0, time1=self.time1)
            im = image.tdplot(plot=False, colorbar=False, params=False,
                              **kwargs)
            return im

        im = data(self.x0, self.y0, self.theta0, self.theta1, self.radius0)
        image = plt.imshow(im, cmap='Greys_r', interpolation='spline36',
                           origin='lower')

        # Remove actual image to not be shown
        # plt.show(block=False)
        # plt.pause(0.0001)
        plt.close('all')

        plt.figure(figsize=(10, 7))
        plt.rcParams.update({'font.size': 12})
        # Set up sliders
        axt1 = plt.axes([0.10, 0.06, 0.33, 0.04], facecolor='lightcyan')
        axt0 = plt.axes([0.10, 0.11, 0.33, 0.04], facecolor='pink')
        axy0 = plt.axes([0.59, 0.11, 0.33, 0.04], facecolor='lightcyan')
        axx0 = plt.axes([0.59, 0.06, 0.33, 0.04], facecolor='pink')
        axrd = plt.axes([0.10, 0.89, 0.34, 0.04], facecolor='pink')

        slider_r0 = Slider(axrd, 'radio', 0, 50, valinit=self.radius0)
        slider_t1 = Slider(axt1, r'$\theta _1(^\circ)$', self.theta1-30,
                           self.theta1+30, valinit=self.theta1)
        slider_t0 = Slider(axt0, r'$\theta _0(^\circ)$', self.theta0-30,
                           self.theta0+30, valinit=self.theta0)
        slider_x0 = Slider(axy0, 'x (pix)', self.x0-30, self.x0+30,
                           valinit=self.x0, dragging=False)
        slider_y0 = Slider(axx0, 'y (pix)', self.y0-30, self.y0+30,
                           valinit=self.y0)

        plt.axes([0.25, 0.25, 0.6, 0.6])

        # Update the image
        def update(val):
            global image
            theta0 = slider_t0.val
            theta1 = slider_t1.val
            x0 = slider_x0.val
            y0 = slider_y0.val
            radius0 = slider_r0.val
            image.set_data(data(int(x0), int(y0), int(theta0), int(theta1),
                                int(radius0)))
            # fig.canvas.draw()
            # plt.clear()
            plt.draw()

        # Change sliders
        slider_r0.on_changed(update)
        slider_t0.on_changed(update)
        slider_t1.on_changed(update)
        slider_x0.on_changed(update)
        slider_y0.on_changed(update)

        plt.show()
        del image

    # Visualization of the time distance plot cahnging (not necessary)

    def cbar_slider(self, **kwargs):

        from matplotlib.widgets import Slider  # , Button, RadioButtons

        fig = plt.figure()
        ax = fig.add_subplot(111)
        fig.subplots_adjust(left=0.15, bottom=0.25)

        image = td(self.array, x0=self.x0, y0=self.y0, theta0=self.theta0,
                   theta1=self.theta1, time0=self.time0, time1=self.time1,
                   radius=self.radius, path=self.path, radius0=self.radius0)
        im = np.array(image.tdplot(plot=False, colorbar=False,
                                   params=False))

        im1 = ax.imshow(im, cmap='Greys_r', origin='lower', aspect='auto',
                        interpolation='spline36', **kwargs)
        fig.colorbar(im1)
        axcolor = 'gainsboro'
        axmin = fig.add_axes([0.1, 0.05, 0.65, 0.03], facecolor=axcolor)
        axmax = fig.add_axes([0.1, 0.1, 0.65, 0.03], facecolor=axcolor)

        min0 = im.min()
        max0 = im.max()
        smin = Slider(axmin, 'Min', -1*abs(im.min()*1.5), abs(im.max()*1.5),
                      valinit=min0, color='dimgray')
        smax = Slider(axmax, 'Max', -1*abs(im.min()*1.5), abs(im.max()*1.5),
                      valinit=max0, color='dimgray')

        def update(val):
            im1.set_clim([smin.val, smax.val])
            fig.canvas.draw()
        smin.on_changed(update)
        smax.on_changed(update)

        plt.show()


if __name__ == "__main__":

    #  Initial parameters
    name = 'hmi_data/DOPP_20110730-M9_3-DIFF.fits'  # Name of datacube
    x0 = 495   # x-center
    y0 = 510   # y-center
    th0 = 0    # Initial angle
    th1 = 40   # End angle

    rad0 = 25   # Skip these pixels in the TD plot
    rad1 = 150  # Final radius (distance) in pixels

    time0 = 5   # Initial frame for calculations
    time1 = 90  # Final frame for calculations

    # Calculates the TD diagram (calls the td class instance with init values)
    data = td(name, x0=x0, y0=y0, theta0=th0, theta1=th1,
              radius0=rad0, radius=rad1, time0=time0, time1=time1)

    # Plots the TD diagram
    plot = data.tdplot(colorbar=True, plot=True, vmin=-510, vmax=510,
                       title='2011-07-30 M9.3 T02:09')

    # Slider to change the colorbar limits (improve contrast)
    data.cbar_slider()

    # Interactive slides to best fit the TD
    data.slider()

    # Plot many around
    # data.test(rows=3,columns=5)
