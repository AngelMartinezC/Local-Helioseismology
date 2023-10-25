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
from datetime import datetime
import polarTransform
import logging


class TimeDistance:

    r_sun = 696  # Mm

    def __init__(self, array, x0, y0, theta0, theta1, radius0=0, radius=None,
                 time0=0, time1=None, readcube=True):
        '''
        Arc section parameters for the TD plot.

        Parameters:
        ----------

        array: `str`
            Specifies the file name from relative path.

        x0, y0: `float`
            Arc section projection center in pixels.

        theta0, theta1: `float`
            Start and end angles in degrees of the arc section.

        radius: `float`, optional
            Outer radius of the arc section. If not set, then the radius is
            calculated to be the minimun distance from (`x0`, `y0`) to one of
            the sides of the datacube. Radius in pixels.

        radius0: `float`, optional
            Inner radius of the section. If not set, it defaults to 0.

        time0, time1: `float`, optional
            Start and end frames of the datacube. If not set it defaults to
            the first and last frame of the datacube.

        Example: A datacube `DOPP.fits` has a suspected sunquake in the
        coordinates (520, 210) between two arcsections (45°, 72°) extending
        out to 150 pix in the image. To calculate and plot the td, we can do:

        >>> from timedistance import td
        >>> import matplotlib.pyplot as td
        >>> data = td('DOPP.fits', x0=520, y0=210, theta0=45, theta1=72,
        >>>           radius=150)
        >>> data.plot()
        >>> plt.show()
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
        self.read_cube = readcube

        if readcube:
            global data, hdr
            data, hdr = read_cube(array)
        self.data = data
        self.hdr = hdr
        if time1 is None:
            time1 = data.shape[0]
        image_td, polarImage, ptSettings = calc_td(data, hdr, x0, y0, theta0,
                                                   theta1, radius0, radius,
                                                   time0, time1)
        self.image_td = image_td
        self.polarImage = polarImage
        self.ptSettings = ptSettings

        # Check existence of some keywords
        self.t_rec = datetime.now().strftime(r"%Y-%m-%d %H:%M:%S")
        try:
            self.t_rec = self.hdr['T_REC']
        except KeyError:
            logging.warning('T_REC keyword not found in header.'
                            f'Set to {self.t_rec}')

        self.daxis = 1
        try:
            self.daxis = self.hdr['DAXIS1']
        except KeyError:
            logging.warning('DAXIS1 keyword not found in header.'
                            f'Set to {self.daxis}')

        self.daxis3 = 1
        try:
            self.daxis3 = self.hdr['DAXIS3']
        except KeyError:
            logging.warning('DAXIS3 keyword not found in header.'
                            f'Set to {self.daxis3}')

    def plot(self, colorbar=True, cmap='gray', interpolation='sinc',
             **kwargs):
        r'''
        Calculates the TD diagram from a given location (`x0`, `y0`) in a
        circular arc enclosed by the angles :math:`\theta_1` and
        :math:`\theta_2`.
        '''

        self.colorbar = colorbar
        self.cmap = cmap

        plt.gcf()
        extent = [0, (self.radius-self.radius0)*self.daxis*TimeDistance.r_sun,
                  0, (self.time1-self.time0)*45/60]
        im = plt.imshow(self.image_td, origin='lower', aspect='auto',
                        interpolation=interpolation, cmap=self.cmap,
                        extent=extent, **kwargs)
        plt.tick_params(direction='out', length=6, width=1.0, colors='k',
                        grid_color='yellow', grid_alpha=0.99)
        plt.title(self.t_rec)
        plt.xlabel("Distance (Mm)")
        plt.ylabel("Time (min)")

        if self.colorbar:
            plt.colorbar(im, label=r'$\Delta$ V$_{LOS}$ (m/s)')

        return self.image_td

    @property
    def toCartesian(self):
        cartesianImage = self.ptSettings.convertToCartesianImage(
            self.polarImage)
        return cartesianImage

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

    def slider(self, **kwargs):

        from matplotlib.widgets import TextBox

        def data_params(x0, y0, t0, t1, radius):
            td_im = TimeDistance(self.array, x0=x0, y0=y0, theta0=t0,
                                 theta1=t1, radius0=self.radius0,
                                 radius=radius, time0=self.time0,
                                 time1=self.time1, readcube=False)
            extent = [0, (radius-self.radius0) * self.daxis *
                      TimeDistance.r_sun, 0, (self.time1-self.time0)*45/60]
            return td_im.image_td, extent  # , td_im.toCartesian

        if self.radius is None:
            self.radius = 200
        first_image = data_params(self.x0, self.y0, self.theta0, self.theta1,
                                  self.radius)

        fig = plt.figure(figsize=(7, 5))
        # plt.subplot(111)
        ax = plt.axes([0.15, 0.30, 0.6, 0.6])
        im_plot = ax.imshow(first_image[0], cmap='gray', origin='lower',
                            interpolation='spline36', **kwargs)
        plt.rcParams.update({'font.size': 12})
        plt.xlabel("Distance (Mm)")
        plt.ylabel("Time (min)")
        ax.set_aspect('auto')

        axx0 = plt.axes([0.15, 0.11, 0.10, 0.04], facecolor='pink')
        axy0 = plt.axes([0.15, 0.06, 0.10, 0.04], facecolor='lightcyan')
        axt0 = plt.axes([0.35, 0.11, 0.10, 0.04], facecolor='lightcyan')
        axt1 = plt.axes([0.35, 0.06, 0.10, 0.04], facecolor='pink')
        axr1 = plt.axes([0.58, 0.06, 0.10, 0.04], facecolor='pink')
        axmin = plt.axes([0.85, 0.11, 0.10, 0.04], facecolor='pink')
        axmax = plt.axes([0.85, 0.06, 0.10, 0.04], facecolor='pink')
        axcb = plt.axes([0.8, 0.30, 0.025, 0.6], facecolor='pink')
        plt.colorbar(im_plot, cax=axcb, label=r'$\Delta $ v$_{LOS}$ (m/s)')

        text_box_x0 = TextBox(axx0, r"$x_0$ (pix)  ", textalignment="center")
        text_box_y0 = TextBox(axy0, r"$y_0$ (pix)  ", textalignment="center")
        text_box_t0 = TextBox(axt0, r"$\theta_0$ (°)  ",
                              textalignment="center")
        text_box_t1 = TextBox(axt1, r"$\theta_1$ (°)  ",
                              textalignment="center")
        text_box_r1 = TextBox(axr1, r"$\rho_1$ (pix)  ",
                              textalignment="center")
        text_box_min = TextBox(axmin, r"cbar$_{vmin}$  ",
                               textalignment="center")
        text_box_max = TextBox(axmax, r"cbar$_{vmax}$  ",
                               textalignment="center")

        def submit(val):
            value = [tb.text for tb in [text_box_x0, text_box_y0,
                                        text_box_t0, text_box_t1,
                                        text_box_r1]]
            if value[0] == '':
                value[0] = self.x0
            if value[1] == '':
                value[1] = self.y0
            if value[2] == '':
                value[2] = self.theta0
            if value[3] == '':
                value[3] = self.theta1
            if value[4] == '':
                value[4] = self.radius
            func1, ext = data_params(float(value[0]), float(value[1]),
                                     float(value[2]), float(value[3]),
                                     float(value[4]))
            im_plot.set_data(func1)
            im_plot.set_extent(ext)
            fig.canvas.draw_idle()

        def submit_cbar(var):
            value = [tb.text for tb in [text_box_min, text_box_max]]
            if value[0] == '':
                value[0] = np.min(first_image[0])
            if value[1] == '':
                value[1] = np.min(first_image[1])
            im_plot.set_clim(int(value[0]), int(value[1]))
            fig.canvas.draw_idle()

        for tb in [text_box_x0, text_box_y0, text_box_t0, text_box_t1,
                   text_box_r1]:
            tb.on_submit(submit)

        text_box_x0.set_val(self.x0)  # Trigger submit with initial string.
        text_box_y0.set_val(self.y0)
        text_box_t0.set_val(self.theta0)
        text_box_t1.set_val(self.theta1)
        text_box_r1.set_val(self.radius)

        for tb in [text_box_min, text_box_max]:
            tb.on_submit(submit_cbar)
        text_box_min.set_val(int(0.5*np.min(first_image[0])))
        text_box_max.set_val(int(0.5*np.max(first_image[0])))

        plt.show()
        # del im_plot

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

        self.x = columns
        self.y = rows
        fig = plt.figure(figsize=(12, 10))
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

                image = TimeDistance(self.array, x0=X0, y0=Y0,
                                     theta0=self.theta0, theta1=self.theta1,
                                     time0=self.time0, time1=self.time1,
                                     radius=rad, radius0=self.radius0,
                                     readcube=False)

                final = image.plot(colorbar=False)
                plt.title(" ")
                plt.xlabel(" ")
                plt.ylabel(" ")
                plt.imshow(final, cmap='Greys_r', origin='lower',
                           interpolation='spline36', aspect='auto')
                ax.text(0.45, 0.1, f"x={X0}\ny={Y0}", transform=ax.transAxes,
                        bbox=dict(boxstyle="round", ec='k',
                                  fc=(1, 1, 1, 0.6)))
        plt.show()

    def cbar_slider(self, **kwargs):

        from matplotlib.widgets import Slider  # , Button, RadioButtons

        fig = plt.figure()
        ax = fig.add_subplot(111)
        fig.subplots_adjust(left=0.15, bottom=0.25)

        image = TimeDistance(self.array, x0=self.x0, y0=self.y0,
                             theta0=self.theta0, theta1=self.theta1,
                             time0=self.time0, time1=self.time1,
                             radius=self.radius, radius0=self.radius0,
                             readcube=False)
        im = np.array(image.plot(colorbar=False))

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


def read_cube(array):
    hdul = fits.open(array)
    hdr = hdul[0].header
    data = hdul[0].data
    print('hello, Cube read')
    return data, hdr


def calc_td(data, hdr, x0, y0, theta0, theta1, radius0, radius, time0, time1):
    # Convert the angles to radians
    t0 = theta0*np.pi/180
    t1 = theta1*np.pi/180

    # To ensure the minimum radius in the plot
    # This is to set finalRadius in converting to polar coordinates
    naxis1 = hdr['NAXIS1']
    naxis2 = hdr['NAXIS2']
    if radius is None:
        mx = naxis1/2
        my = naxis2/2
        if (x0 <= mx) and (y0 <= my):
            rad = min(x0, y0)
        elif (y0 > my) and (x0 <= mx):
            rad = min(x0, my*2 - y0)
        elif (x0 > mx) and (y0 <= my):
            rad = min(y0, mx*2 - x0)
        elif (x0 > mx) and (y0 > my):
            rad = min(mx*2 - x0, my*2 - y0)
    else:
        rad = radius
    rad = int(rad)

    image_td = []
    for ij in range(time0, time1):
        polarImage, ptSettings = polarTransform.convertToPolarImage(
            data[ij], center=[x0, y0],
            initialRadius=radius0, finalRadius=rad,
            initialAngle=t0, finalAngle=t1,
            radiusSize=rad-radius0, order=0)
        slices = []
        for i in range(len(polarImage[0, :])):
            slices.append(np.mean(polarImage[:, i]))
        image_td.append(slices)
    image_td = np.array(image_td)
    return image_td, polarImage, ptSettings


if __name__ == "__main__":

    from timedistance import TimeDistance as td
    import os

    os.system('clear')

    name = './DOPP_DIFF.fits'  # Name of datacube
    x0 = 503  # x-center
    y0 = 512  # y-center
    th0 = 194  # Initial angle
    th1 = 258  # End angle

    rad0 = 0  # Skip these pixels in the TD plot
    rad1 = 200  # Final radius (distance) in pixels

    time0 = 100  # Initial frame for calculations
    time1 = 200  # Final frame for calculations

    # Calculates the TD diagram (calls the td class instance with init values)
    data = td(name, x0=x0, y0=y0, theta0=th0, theta1=th1, radius0=rad0,
              radius=rad1, time0=time0, time1=time1)

    # plt.figure()
    # plt.subplot(111)
    # image = data.plot(colorbar=True, vmin=-300, vmax=300,
    #                   interpolation='sinc')
    # image_cart = data.toCartesian
    # plt.show()

    data.slider()

    # data.cbar_slider()
    # data.test()
