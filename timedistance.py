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


class td:

    r_sun = 696  # Mm
    rsun_pix = 1884  # Value for HMI

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

    def plot(self, colorbar=True, cmap='gray', interpolation='sinc',
             **kwargs):
        r'''
        Calculates the TD diagram from a given location (`x0`, `y0`) in a
        circular arc enclosed by the angles :math:`\theta_1` and
        :math:`\theta_2`.
        '''

        self.colorbar = colorbar
        self.cmap = cmap

        # Check existence of some keywords
        t_rec = datetime.now().strftime(r"%Y-%m-%d %H:%M:%S")
        try:
            t_rec = self.hdr['T_REC']
        except KeyError:
            print(f'T_REC keyword not found in header. Set to {t_rec}')

        daxis = 1
        try:
            daxis = self.hdr['DAXIS1']
        except KeyError:
            print(f'DAXIS1 keyword not found in header. Set to {daxis}')

        daxis3 = 1
        try:
            daxis3 = self.hdr['DAXIS3']
        except KeyError:
            print(f'DAXIS3 keyword not found in header. Set to {daxis3}')

        plt.gcf()
        extent = [0, (self.radius-self.radius0)*daxis*td.r_sun, 0,
                  (self.time1-self.time0)*45/60]
        im = plt.imshow(self.image_td, origin='lower', aspect='auto',
                        interpolation=interpolation, cmap=self.cmap,
                        extent=extent, **kwargs)
        plt.tick_params(direction='out', length=6, width=1.0, colors='k',
                        grid_color='yellow', grid_alpha=0.99)
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

        # from matplotlib.widgets import Slider
        from matplotlib.widgets import TextBox

        def data_params(x0, y0, t0, t1, radius0):
            td_im = td(self.array, x0=x0, y0=y0, theta0=t0, theta1=t1,
                       radius0=radius0, radius=self.radius,
                       time0=self.time0, time1=self.time1, readcube=False)
            return td_im.image_td

        first_image = data_params(self.x0, self.y0, self.theta0, self.theta1,
                                  self.radius0)

        print('\n\n ', self.radius0, self.radius)
        fig = plt.figure(figsize=(10, 7))
        ax = plt.subplot(121)
        # ax = plt.axes([0.25, 0.25, 0.6, 0.6])
        global im_plot
        im_plot = ax.imshow(first_image, cmap='gray', origin='lower',
                            interpolation='spline36', **kwargs)
        print(first_image.shape)
        plt.rcParams.update({'font.size': 12})
        # Set up sliders
        # axt1 = plt.axes([0.10, 0.06, 0.33, 0.04], facecolor='lightcyan')
        # axt0 = plt.axes([0.10, 0.11, 0.33, 0.04], facecolor='pink')
        # axy0 = plt.axes([0.59, 0.11, 0.33, 0.04], facecolor='lightcyan')
        # axx0 = plt.axes([0.59, 0.06, 0.33, 0.04], facecolor='pink')
        # axrd = plt.axes([0.10, 0.89, 0.34, 0.04], facecolor='pink')

        # slider_t1 = Slider(axt1, r'$\theta _1(^\circ)$', self.theta1-30,
        #                    self.theta1+30, valinit=self.theta1)
        # slider_t0 = Slider(axt0, r'$\theta _0(^\circ)$', self.theta0-30,
        #                    self.theta0+30, valinit=self.theta0)
        # slider_x0 = Slider(axy0, 'x (pix)', self.x0-30, self.x0+30,
        #                    valinit=self.x0, dragging=False)
        # slider_y0 = Slider(axx0, 'y (pix)', self.y0-30, self.y0+30,
        #                    valinit=self.y0)
        # slider_r0 = Slider(axrd, 'radio', 0, 50, valinit=self.radius0)

        # # Update the image
        # def update(val):
        #     radius0 = slider_r0.val
        #     theta0 = slider_t0.val
        #     theta1 = slider_t1.val
        #     x0 = slider_x0.val
        #     y0 = slider_y0.val
        #     im_plot.set_data(data_params(int(x0), int(y0), int(theta0),
        #                                  int(theta1), int(radius0)))
        #     fig.canvas.draw_idle()
        #     # fig.canvas.draw()
        #     # plt.clear()
        #     plt.draw()

        def submit(expression):
            theta0 = eval(expression, {'theta0': self.theta0})
            x0 = eval(expression, {'x0': str(self.y0)})
            # y0 = eval(expression, {'y0': self.y0})
            # theta1 = eval(expression, {'t1': self.theta1})
            print(theta0, x0)
            print('\n \n \n \n')
            im_plot.set_data(data_params(x0, self.y0, theta0,
                                         self.theta1, 0))
            fig.canvas.draw_idle()
            # ax.relim()
            # ax.autoscale_view()
            plt.draw()

        # axbox = fig.add_axes([0.1, 0.05, 0.8, 0.075])
        axx0 = plt.axes([0.10, 0.11, 0.33, 0.04], facecolor='pink')
        axy0 = plt.axes([0.10, 0.06, 0.3, 0.04], facecolor='lightcyan')
        axt0 = plt.axes([0.59, 0.06, 0.33, 0.04], facecolor='pink')
        axt1 = plt.axes([0.59, 0.11, 0.33, 0.04], facecolor='lightcyan')
        text_box_x0 = TextBox(axx0, r"$x_0$ ", textalignment="center")
        text_box_y0 = TextBox(axy0, r"$y_0$ ", textalignment="center")
        text_box_t0 = TextBox(axt0, r"$\theta_0$ ", textalignment="center")
        text_box_t1 = TextBox(axt1, r"$\theta_0$ ", textalignment="center")
        text_box_x0.on_submit(submit)
        # text_box_y0.on_submit(submit)
        text_box_t0.on_submit(submit)
        # text_box_t1.on_submit(submit)
        text_box_x0.set_val(self.x0)  # Trigger submit with initial string.
        # text_box_y0.set_val(self.y0)
        text_box_t0.set_val(self.theta0)
        # text_box_t1.set_val(self.theta1)

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
                           time1=self.time1, radius=rad, radius0=self.radius0)

                final = image.tdplot(plot=False, colorbar=False, notime=False)
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

        image = td(self.array, x0=self.x0, y0=self.y0, theta0=self.theta0,
                   theta1=self.theta1, time0=self.time0, time1=self.time1,
                   radius=self.radius, radius0=self.radius0)
        # im = np.array(image.tdplot(plot=False, colorbar=False,
        im = np.array(image.tdplot(colorbar=False))
        #                            params=False))

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

    from timedistance import td
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
    # plt.subplot(121)
    # image = data.plot(colorbar=True, vmin=-300, vmax=300,
    #                   interpolation='sinc')
    #
    # plt.show()
    # exit()
    image_cart = data.toCartesian
    data.slider()
    plt.show()
    exit()

    print('End')
    # #  Initial parameters
    # # name = 'hmi_data/DOPP_20110730-M9_3-DIFF.fits'  # Name of datacube
    # name = './DOPP_DIFF.fits'  # Name of datacube
    # x0 = 503  # x-center
    # y0 = 512  # y-center
    # th0 = 194  # Initial angle
    # th1 = 258  # End angle

    # rad0 = 0  # Skip these pixels in the TD plot
    # rad1 = 200  # Final radius (distance) in pixels

    # time0 = 100  # Initial frame for calculations
    # time1 = 200  # Final frame for calculations

    # # Calculates the TD diagram (calls the td class inst with init values)
    # data = td(name, x0=x0, y0=y0, theta0=th0, theta1=th1, radius0=rad0,
    #           radius=rad1, time0=time0, time1=time1)

    # fig = plt.figure(figsize=(12, 5))
    # plt.subplots_adjust(wspace=0.2, left=0.1, right=0.95)

    # plt.subplot(121)
    # image = data.tdplot(colorbar=False, vmin=-300, vmax=300,
    #                     interpolation='sinc')

    # plt.subplot(122)
    # angle1, rho1, time1, ampl1 = np.loadtxt(
    #     './../../../test_bc_filter/td/OUT', unpack=True)
    # r_sun = 696  # Mm
    # frame0 = 30
    # image = data.tdplot(colorbar=False, vmin=-300, vmax=300,
    #                     interpolation='sinc')
    # # print(im)
    # plt.ylabel(' ')
    # plt.plot(rho1*r_sun, (time1 + frame0*45)/60)
    # plt.xlim(0, image.shape[1]*0.0005*r_sun)
    # plt.ylim(0, image.shape[0]*45/60)
    # # cb_ax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
    # # cbar = fig.colorbar(im, cax=cb_ax)
    # plt.show()

    # plt.show()
    # exit()
    # Plots the TD diagram

    # Slider to change the colorbar limits (improve contrast)
    # data.cbar_slider()

    # Interactive slides to best fit the TD
    # data.slider()

    # Plot many around
    # data.test(rows=3,columns=5)
