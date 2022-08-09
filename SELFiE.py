#STANDARD LIBRARIES
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib import style
import matplotlib.image as mpimg
import scipy.optimize as sciopt
import scipy.interpolate as interp
import pandas as pd

#THIRD PARTY LIBRARIES
import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import askopenfilename
import tkinter.messagebox as tkMessageBox
import tkinter.simpledialog as tkSimpleDialog
from tkinter.simpledialog import Dialog
import astropy.io.ascii as ascii
from astropy.table import Table

#LOCAL APPLICATIONS
import linefittingfunctions as lff
import spectrumstitchingfunctions as ssf

#CONSTANTS
C_kms = 2.99792458e5  #Speed of light [km/s]

#DEFINE FONTS TO USE ON BUTTONS/LABELS
XL_FONT = ("Verdana", 22)
LARGE_FONT = ("Verdana", 12)
NORM_FONT = ("Verdana", 10)
SMALL_FONT = ("Verdana", 8)
style.use("ggplot")

#FUNCTIONS
def reset():
    """Opens a new instance of the GUI
    :return: None
    """
    global app

    app.destroy()  #Destroys (but doesn't close old window... closes it if "Reset" button is used again though)
    app = SELFiE()
    app.geometry("1280x720")
    app.mainloop()

#MAIN APP IS DEFINED HERE
class SELFiE(tk.Tk):

    def __init__(self, *args, **kwargs):

        tk.Tk.__init__(self, *args, **kwargs)

        #tk.Tk.iconbitmap(self, default="duck.ico") #Need an actual icon file here :( Try to make a duckface in photoshop!
        tk.Tk.wm_title(self, "SELFiE")

        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}

        for F in (HomePage, ContinuumPage, STISPage, COSPage):

            frame = F(container, self)

            self.frames[F] = frame

            frame.grid(row=0, column=0, sticky="nsew") #Display pages when app is initialized

        self.show_frame(HomePage) #This is the front page of the app

    def show_frame(self, cont): #Pulls up frame when button is clicked

        frame = self.frames[cont]
        frame.tkraise()

############DEFINES METHODS FOR THE FIRST PAGE###########################################################
class HomePage(tk.Frame):

    def __init__(self, parent, controller):

        tk.Frame.__init__(self, parent, bg="gray")

        label = tk.Label(self, text="Welcome to SELFiE: STIS/COS Emission Line FItting and Continuum Extraction", font=XL_FONT) #Title of app
        label.pack(pady=10, padx=10)

        button = ttk.Button(self, text="Extract Continuum",
                            command=lambda: controller.show_frame(ContinuumPage)) #Need to use lambda function with "command=" to pass function with parameters (otherwise nothing will happen); this is the button you click to go to the "Extrapolate Continuum" page

        button.pack()

        button2 = ttk.Button(self, text="Measure STIS Emission Lines",
                             command=lambda: controller.show_frame(STISPage)) #This is the button you click to get to the "Measure STIS Emission Lines" page

        button2.pack()

        button3 = ttk.Button(self, text="Measure COS Emission Lines",
                             command=lambda: controller.show_frame(COSPage)) #This is the button you click to get to the "Measure COS Emission Lines" page

        button3.pack()

        #####READS IN DUCK MASCOT PHOTO!#######################################
        g = Figure(figsize=(5,5), dpi=100)
        b = g.add_subplot(111)
        img = mpimg.imread("duck.png")
        b.imshow(img)
        #####################################################################

        canvas = FigureCanvasTkAgg(g, self) #Allows image to display in interactive plotting window
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self) #Initialize matplotlib toolbar
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

#########DEFINES METHODS FOR THE "EXTRACT CONTINUUM" PAGE - NONE ADDED YET######################
class ContinuumPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent, bg="gray")
        label = tk.Label(self, text="Extract Continuum", font=LARGE_FONT) #Title of page
        label.pack(pady=10, padx=10)

        button1 = ttk.Button(self, text="Back to Home",
                             command=lambda: controller.show_frame(HomePage)) #Click this button to return to the start page

        button1.pack()

        button2 = ttk.Button(self, text="Measure STIS Emission Lines",
                             command=lambda: controller.show_frame(STISPage)) #Click this button to measure STIS emission lines

        button2.pack()

        button3 = ttk.Button(self, text="Measure COS Emission Lines",
                             command=lambda: controller.show_frame(COSPage)) #Click this button to measure COS emission lines

        button3.pack()

        button4 = ttk.Button(self, text="Reset",
                             command=lambda: reset()) #Click this button to reset the GUI - it will open a new instance, and the old one can be closed
        button4.pack()
        button4.place(x=1175, y=0)

        #####ENTRY BOX FOR USER TO ENTER FILE NAME OF SPECTRUM WITH LINES TO FIT#######################
        label1 = ttk.Label(self, text="Spectrum Filename: ")
        label1.pack(side=tk.TOP, pady=10)
        E1 = ttk.Entry(self)
        E1.pack(side=tk.TOP, pady=10)

        submit = ttk.Button(self, text="Enter", command=lambda: getfilename("Data File with Lines to Fit:")) #Get list of lines to fit! Should have same format as France2012_Table2.csv to work correctly
        submit.pack(side=tk.TOP, pady=10)

        #####CAN USE THIS PAGE TO FIT CONTINUUM AROUND CO ABSORPTIONS OR TO FIT ENTIRE FUV CONTINUUM USING KEVIN'S ALGORITHM - THAT'S WHY THERE'S NO OPTION TO READ IN A LIST OF LINES################

        #Function that will be repeated to allow user to keep picking regions of the continuum to fit - can't get rid of a region, but there is an option to re-do it at the end
        def fitdisplay():

            global cont_label, cont_bttn, done_bttn

            cont_label.destroy()
            cont_bttn.destroy()
            done_bttn.destroy()

            global continuumfit_label, continuumfit_bttn
            continuumfit_label = ttk.Label(self, text="Navigate to a nearby region free of emission lines, then click 'Enter':") #Fit continuum to put upper limit on RMS flux
            continuumfit_label.place(bordermode=tk.OUTSIDE, x=0, y=0)

            continuumfit_bttn = ttk.Button(self, text="Enter", command = lambda: pickcontinuum())
            continuumfit_bttn.place(x=0, y=30)

        #Re-plots full spectrum, prompts user to enter central wavelength - then will plot +/- 20 Angstroms around central wavelength, but can still pan across full spectrum
        def restartfit():

            global happy_bttn_yes, happy_bttn_no, start_label, start_entry, start_bttn, ax1, new_bttn, saved_label

            happy_bttn_yes.destroy()
            happy_bttn_no.destroy()
            new_bttn.destroy()
            saved_label.destroy()

            ax1.clear()
            ax1.plot(wavelength, flux, "k") #First, plot full spectrum - user will be able to interact with it before fitting emission lines, but should add option to go back to this
            ax1.set_xlabel("Wavelength (Angstroms)")
            ax1.set_ylabel("Flux")
            canvas.draw()

            global continuum_wavelength_list, continuum_flux_list, continuum_error_list
        #Will store user-selected ranges for continuum fits in these lists (ranges don't need to be the same length)
            continuum_wavelength_list = []
            continuum_flux_list = []
            continuum_error_list = []

            start_label = ttk.Label(self, text="Enter a central wavelength (Angstroms):") #Set up widgets to guide user through fitting process
            start_label.pack(side=tk.TOP)
            start_label.place(bordermode=tk.OUTSIDE, x=0, y=0)

            start_entry = ttk.Entry(self)
            start_entry.pack(side=tk.TOP)
            start_entry.place(x=0, y=30)

            start_bttn = ttk.Button(self, text="Enter", command = lambda: choose_regions())
            start_bttn.pack(side=tk.TOP)
            start_bttn.place(x=0, y=60)

        #If the continuum fit is good, now normalize fluxes, save figure and file - will need to adjust this to save full FUV continuum
        def savefit():

            global central_wavelength, spline_fit, continuum_wavelengths_flattened, ax1, canvas, f, happy_bttn_yes, happy_bttn_no

            #Destroy buttons that asked whether or not to keep the fit
            happy_bttn_yes.destroy()
            happy_bttn_no.destroy()

            #Get range of data within +10/-1.5 Angstroms from line center (for CO absorptions)
            line_wavelengths, line_fluxes, line_errors = lff.get_range(wavelength, flux, flux_err, np.amin(continuum_wavelengths_flattened), np.amax(continuum_wavelengths_flattened))

            #Divide data by continuum (obtained from spline fit)
            norm_line_fluxes = line_fluxes / spline_fit(line_wavelengths)
            #norm_line_fluxes = line_fluxes / 1.0 #For IR CO lines that are already normalized well (i.e. not embedded in emission lines)

            #Plot normalized flux - still optimized for CO lines
            ax1.clear()
            ax1.set_title("Normalized Flux: " + str(central_wavelength) + " Angstroms")
            ax1.plot(line_wavelengths, norm_line_fluxes, "k")
            ax1.hlines(1.0, np.amin(line_wavelengths), np.amax(line_wavelengths), "r", linewidth = 4., linestyle = "--")
            ax1.set_xlim(np.amin(line_wavelengths), np.amax(line_wavelengths))
            ax1.set_xlabel(r"Wavelength $\left( \AA \right)$")
            plt.gca().xaxis.set_ticks(np.array([central_wavelength-20., central_wavelength-15., central_wavelength-10., central_wavelength-5., central_wavelength, central_wavelength+5., central_wavelength+10., central_wavelength+15., central_wavelength+20.]))
            plt.gca().set_xticklabels(np.array([str(central_wavelength-20.)[0:6], str(central_wavelength-15.)[0:6], str(central_wavelength-10.)[0:6], str(central_wavelength-5.)[0:6], str(central_wavelength)[0:6], str(central_wavelength+5.)[0:6], str(central_wavelength+10.)[0:6], str(central_wavelength+15.)[0:6], str(central_wavelength+20.)[0:6]]))
            #ax1.set_ylim(0.0, 1.5)
            #ax1.set_ylim(0.9, 1.05)
            ax1.set_ylim(0.0, np.amax(norm_line_fluxes)+0.1)
            ax1.set_ylabel("Normalized Flux")
            ax1.minorticks_on()
            canvas.draw()
            f.savefig("continuumnormalized_" + str(central_wavelength) + ".pdf", format = "pdf")

            #Write out parameters to file - can make new plots with this if necessary
            linefit_table = Table([line_wavelengths, line_fluxes, line_errors, spline_fit(line_wavelengths), norm_line_fluxes], names = ("Wavelength (Angstroms)", "Line Flux (erg/s/cm^2/Angstrom)", "Line Flux Errors", "Continuum Flux (erg/s/cm^2/Angstrom)", "Normalized Line Flux (erg/s/cm^2/Angstrom)"))
            #linefit_table = Table([line_wavelengths, line_fluxes, line_errors, np.ones(np.size(line_wavelengths)), norm_line_fluxes], names = ("Wavelength (Angstroms)", "Line Flux (erg/s/cm^2/Angstrom)", "Line Flux Errors", "Continuum Flux (erg/s/cm^2/Angstrom)", "Normalized Line Flux (erg/s/cm^2/Angstrom)")) #Just write out 1s for continuum flux if you don't want to normalize (for IR CO lines that aren't embedded in emissions)

            ascii.write(linefit_table, "continuumfit_" + str(central_wavelength) + ".txt")

            global saved_label
            saved_label = ttk.Label(self, text="Saved! Fit another region?")
            saved_label.place(bordermode=tk.OUTSIDE, x=0, y=0) #Label will be placed in upper left corner of GUI
            global new_bttn
            new_bttn = ttk.Button(self, text="Enter", command=lambda: restartfit())
            new_bttn.place(bordermode=tk.OUTSIDE, x=0, y=30)

            #Define buttons, but don't place - will need to be destroyed by restartfit() function
            happy_bttn_yes = ttk.Button(self, text="Keep Fit", command=lambda: savefit())
            happy_bttn_no = ttk.Button(self, text="Discard Fit", command=lambda: restartfit())

        #After all the clean continuum regions are selected, concatenate the ranges/apply cubic spline fit
        def done():

            global cont_label, cont_bttn, done_bttn, continuum_wavelengths_flattened

            cont_label.destroy()
            cont_bttn.destroy()
            done_bttn.destroy()

            #Flatten list of lists into 1D array
            continuum_wavelengths_flattened = np.concatenate(continuum_wavelength_list)
            continuum_fluxes_flattened = np.concatenate(continuum_flux_list)
            #Smooth fluxes by 7-pixel boxcar to remove noise spikes before fitting
            #lowres_continuumfluxes = smooth(continuum_wavelengths_flattened, continuum_fluxes_flattened, 1., "Boxcar") #Don't need to do this for CRIRES spectra
            continuum_errors_flattened = np.concatenate(continuum_error_list)

            nan_indices = np.where(np.isnan(continuum_fluxes_flattened) == True)[0]
            continuum_wavelengths_flattened = np.delete(continuum_wavelengths_flattened, nan_indices)
            lowres_continuumfluxes = np.delete(continuum_fluxes_flattened, nan_indices)
            continuum_errors_flattened = np.delete(continuum_errors_flattened, nan_indices)

            global spline_fit
            #Cubic spline fit to "clean" continuum regions
            spline_fit = interp.UnivariateSpline(continuum_wavelengths_flattened, lowres_continuumfluxes, k=2) #Use k=2 for CO absorptions
            #spline_fit = interp.interp1d(continuum_wavelengths_flattened, lowres_continuumfluxes, kind="slinear")

            #spline_fit = interp.UnivariateSpline(continuum_wavelengths_flattened, continuum_fluxes_flattened, k=2)
            #Model array to see how fit compares to data
            x_test = np.arange(np.amin(continuum_wavelengths_flattened), np.amax(continuum_wavelengths_flattened), 0.01)

            #Plot spline fit over data that fit was derived from/decide whether or not to keep it
            global ax1, canvas, f
            ax1.clear()
            ax1.set_title("Continuum Fit: " + str(central_wavelength) + " Angstroms")
            for j in range(0, np.size(continuum_wavelength_list)):
                ax1.plot(continuum_wavelength_list[j], continuum_flux_list[j], "k")
            ax1.plot(x_test, spline_fit(x_test), "r", lw=4.)
            ax1.set_xlim(central_wavelength - 20., central_wavelength + 20.)
            plt.gca().xaxis.set_ticks(np.array([central_wavelength-20., central_wavelength-15., central_wavelength-10., central_wavelength-5., central_wavelength, central_wavelength+5., central_wavelength+10., central_wavelength+15., central_wavelength+20.]))
            plt.gca().set_xticklabels(np.array([str(central_wavelength-20.)[0:6], str(central_wavelength-15.)[0:6], str(central_wavelength-10.)[0:6], str(central_wavelength-5.)[0:6], str(central_wavelength)[0:6], str(central_wavelength+5.)[0:6], str(central_wavelength+10.)[0:6], str(central_wavelength+15.)[0:6], str(central_wavelength+20.)[0:6]]))
            ax1.set_ylim(np.amin(continuum_fluxes_flattened), np.amax(continuum_fluxes_flattened))
            canvas.draw()
            f.savefig("continuumfit_" + str(central_wavelength) + ".pdf", format = "pdf")

            global happy_bttn_yes, happy_bttn_no

            #Clicking this will normalize the fluxes and save the output to a plot/text file
            happy_bttn_yes = ttk.Button(self, text="Keep Fit", command=lambda: savefit())
            happy_bttn_yes.place(bordermode=tk.OUTSIDE, x=500, y=625)

            #Clicking this will discard the fit and re-plot the original spectrum - can start over from there
            happy_bttn_no = ttk.Button(self, text="Discard Fit", command=lambda: restartfit())
            happy_bttn_no.place(bordermode=tk.OUTSIDE, x=700, y=625)

        #Get x-coordinate at blue side of clean continuum region
        def getleftside(event1):

            global continuum_left, guide_label_blue, canvas, cid_blue

            continuum_left = event1.xdata #Save x-coordinate of button click event (event3)

            guide_label_blue.destroy() #Destroy label prompting user to click left side of line
            canvas.mpl_disconnect(cid_blue) #Disconnect mouse from canvas

        #Get x-coordinate at red side of clean continuum region/ask to add more regions
        def getrightside(event2):

            global continuum_right, guide_label_red, canvas, cid_red, continuum_left

            continuum_right = event2.xdata #Save x-coordinate of button click event (event3)

            guide_label_red.destroy() #Destroy label prompting user to click left side of line
            canvas.mpl_disconnect(cid_red) #Disconnect mouse from canvas

            #Extract clean continuum data from full spectrum/append to lists
            continuum_wavelengths, continuum_fluxes, continuum_errors = lff.get_range(wavelength, flux, flux_err, continuum_left, continuum_right)

            continuum_wavelength_list.append(continuum_wavelengths)
            continuum_flux_list.append(continuum_fluxes)
            continuum_error_list.append(continuum_errors)

            #Create buttons to allow more regions to be selected or fit the data that has already been chosen
            global cont_label, cont_bttn, done_bttn, happy_bttn_no
            cont_label = ttk.Label(self, text="Click 'Enter' to add another region, or click 'Done' to fit selected regions: ")
            cont_label.place(bordermode=tk.OUTSIDE, x=0, y=0) #Label will be placed in upper left corner of GUI

            cont_bttn = ttk.Button(self, text="Enter", command = lambda: fitdisplay())
            cont_bttn.place(x=0, y=30)

            done_bttn = ttk.Button(self, text="Done", command = lambda: done())
            done_bttn.place(x=100, y=30)

            #Plots clean continuum that has already been chosen
            global ax1, cont_fluxes, central_wavelength
            ax1.plot(continuum_wavelengths, continuum_fluxes, "r")
            ax1.set_ylim(np.amin(cont_fluxes), np.amax(cont_fluxes))
            ax1.set_xlim(central_wavelength - 20., central_wavelength + 20.)
            canvas.draw()

        #Prompts to select bluest/reddest wavelengths of clean continuum region
        def pickcontinuum():

            global continuumfit_label, continuumfit_bttn
            continuumfit_label.destroy()
            continuumfit_bttn.destroy()

            global guide_label_blue
            guide_label_blue = ttk.Label(self, text="Select the bluest wavelength in the clean region:")
            guide_label_blue.place(bordermode=tk.OUTSIDE, x=0, y=0) #Label will be placed in upper left corner of GUI

            global cid_blue
            cid_blue = canvas.mpl_connect("button_press_event", lambda event1: getleftside(event1)) #Connect to canvas to get user click location
            tk.Frame.wait_window(guide_label_blue)

            global guide_label_red
            guide_label_red = ttk.Label(self, text="Select the reddest wavelength in the clean region:")
            guide_label_red.place(bordermode=tk.OUTSIDE, x=0, y=0) #Label will be placed in upper left corner of GUI

            global cid_red
            cid_red = canvas.mpl_connect("button_press_event", lambda event2: getrightside(event2)) #Connect to canvas to get user click location

        #Plots +/-20 Angstrom range around central wavelength, initial prompt to navigate to "clean" continuum region
        def choose_regions():

            global start_label, start_entry, start_bttn

            global wavelength, flux, flux_err

            global central_wavelength, ax1, canvas

            central_wavelength = float(start_entry.get())

            start_label.destroy()
            start_entry.destroy()
            start_bttn.destroy()

            #Just want to use this range to set y-axis limits on plot - don't actually want to limit range of spectrum that's plotted (may need to go further than +/- 50 Angstroms to find enough clean continuum regions)
            global cont_wavelengths, cont_fluxes, cont_errors
            #cont_wavelengths, cont_fluxes, cont_errors = get_range(wavelength, flux, flux_err, central_wavelength - 20., central_wavelength + 20.)
            cont_wavelengths, cont_fluxes, cont_errors = lff.get_range(wavelength, flux, flux_err, central_wavelength - 3., central_wavelength + 3.)

            nan_indices = np.where(np.isnan(cont_fluxes) == True)[0] #Get rid of NaNs (e.g. from masking)
            cont_wavelengths = np.delete(cont_wavelengths, nan_indices)
            cont_fluxes = np.delete(cont_fluxes, nan_indices)
            cont_errors = np.delete(cont_errors, nan_indices)

            ax1.set_title("Continuum Fit: " + str(central_wavelength) + " Angstroms")
            ax1.set_xlim(central_wavelength - 10., central_wavelength + 10.)
            ax1.set_ylim(np.amin(cont_fluxes), np.amax(cont_fluxes))
            ax1.minorticks_on()
            canvas.draw()

            global continuumfit_label, continuumfit_bttn
            continuumfit_label = ttk.Label(self, text="Navigate to a nearby region free of emission lines, then click 'Enter':") #Fit continuum to put upper limit on RMS flux
            continuumfit_label.place(bordermode=tk.OUTSIDE, x=0, y=0)

            continuumfit_bttn = ttk.Button(self, text="Enter", command = lambda: pickcontinuum())
            continuumfit_bttn.place(x=0, y=30)

        #Start by reading in file with spectrum to look at/fit
        def getfilename(newtext):

            filename = E1.get() #Read in filename that user entered - should add something here to prompt another filename entry if an invalid string is entered; file should have wavelength as "col1", flux as "col2," and flux errors as "col3"
            E1.destroy() #Destroy text entry box

            label1.destroy() #Destroy text entry label

            submit.destroy() #Destroy "Enter" button for this text box

            data_table = ascii.read(filename) #Will work with this file unless the GUI is reset and a new one is input
            global wavelength, flux, flux_err
            wavelength = np.array(data_table["col1"]).astype(np.float64)
            flux = np.array(data_table["col2"]).astype(np.float64)
            flux_err = np.array(data_table["col3"]).astype(np.float64)

            #####INITIALIZE INTERACTIVE PLOT FOR SPECTRUM########################################
            global f, ax1, canvas, start_label, start_entry, start_bttn
            f = Figure()
            ax1 = f.add_subplot(111)
            ax1.plot(wavelength, flux, "k") #First, plot full spectrum - user will be able to interact with it before fitting emission lines, but should add option to go back to this
            ax1.set_xlabel("Wavelength (Angstroms)")
            ax1.set_ylabel("Flux")
            canvas = FigureCanvasTkAgg(f, self)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP) #Outputs plot to GUI window

            toolbar = NavigationToolbar2Tk(canvas, self) #Initializes matplotlib toolbar
            toolbar.update()
            canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=False)

            global continuum_wavelength_list, continuum_flux_list, continuum_error_list
        #Will store user-selected ranges for continuum fits in these lists (ranges don't need to be the same length)
            continuum_wavelength_list = []
            continuum_flux_list = []
            continuum_error_list = []

            start_label = ttk.Label(self, text="Enter a central wavelength (Angstroms):") #Set up widgets to guide user through fitting process
            start_label.pack(side=tk.TOP)
            start_label.place(bordermode=tk.OUTSIDE, x=0, y=0)

            start_entry = ttk.Entry(self)
            start_entry.pack(side=tk.TOP)
            start_entry.place(x=0, y=30)

            start_bttn = ttk.Button(self, text="Enter", command = lambda: choose_regions())
            start_bttn.pack(side=tk.TOP)
            start_bttn.place(x=0, y=60)

#############DEFINE METHODS FOR MEASURING STIS EMISSION LINES - NONE HAVE BEEN ADDED YET!!###############
class STISPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent, bg="gray")
        label = tk.Label(self, text="Measure STIS Emission Lines", font=LARGE_FONT) #Title of page
        label.pack(pady=10, padx=10)

        button1 = ttk.Button(self, text="Back to Home",
                             command=lambda: controller.show_frame(HomePage)) #Click this button to go back to the home page

        button1.pack()

        button2 = ttk.Button(self, text="Extract Continuum",
                             command=lambda: controller.show_frame(ContinuumPage)) #Click this button to extract the FUV continuum

        button2.pack()

        button3 = ttk.Button(self, text="Measure COS Emission Lines",
                             command=lambda: controller.show_frame(COSPage)) #Click this button to measure the COS emission lines

        button3.pack()

        button4 = ttk.Button(self, text="Reset",
                             command=lambda: reset()) #Click this button to reset the GUI - it will open a new instance, and the old one can be closed
        button4.pack()
        button4.place(x=1175, y=0)

        #####ENTRY BOX FOR USER TO ENTER FILE NAME OF SPECTRUM WITH LINES TO FIT#######################
        label1 = ttk.Label(self, text="Spectrum Filename: ")
        label1.pack(side=tk.TOP, pady=10)
        E1 = ttk.Entry(self)
        E1.pack(side=tk.TOP, pady=10)

        submit = ttk.Button(self, text="Enter", command=lambda: getfilename("Data File with Lines to Fit:")) #Get list of lines to fit! Should have same format as France2012_Table2.csv to work correctly
        submit.pack(side=tk.TOP, pady=10)

        def getfilename(newtext):

            filename = E1.get() #Read in filename that user entered - should add something here to prompt another filename entry if an invalid string is entered; file should have wavelength as "col1", flux as "col2," and flux errors as "col3"
            E1.destroy() #Destroy text entry box

            label1.destroy() #Destroy text entry label

            submit.destroy() #Destroy "Enter" button for this text box

            data_table = ascii.read(filename) #Will work with this file unless the GUI is reset and a new one is input
            global wavelength, flux, flux_err
            wavelength = np.array(data_table["col1"]).astype(np.float64)
            flux = np.array(data_table["col2"]).astype(np.float64)
            flux_err = np.array(data_table["col3"]).astype(np.float64)

            #####INITIALIZE INTERACTIVE PLOT FOR SPECTRUM########################################
            global f, ax1, canvas
            f = Figure()
            ax1 = f.add_subplot(111)
            ax1.plot(wavelength, flux, "k") #First, plot full spectrum - user will be able to interact with it before fitting emission lines, but should add option to go back to this
            ax1.set_xlabel("Wavelength (Angstroms)")
            ax1.set_ylabel("Flux")
            canvas = FigureCanvasTkAgg(f, self)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP) #Outputs plot to GUI window

            toolbar = NavigationToolbar2Tk(canvas, self) #Initializes matplotlib toolbar
            toolbar.update()
            canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=False)

############DEFINE METHODS FOR COS EMISSION LINE FITTING###################################################
class COSPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent, bg="gray")
        #Title of page
        title_label = ttk.Label(self, text="Measure COS Emission Lines", font=LARGE_FONT)
        title_label.pack(pady=10, padx=10)

        #Click this button to go back to the home page
        home_button = ttk.Button(self, text="Back to Home",
                                 command=lambda: controller.show_frame(HomePage))
        home_button.pack()

        #Click this button to extract the FUV continuum
        continuum_button = ttk.Button(self, text="Extract Continuum",
                                      command=lambda: controller.show_frame(ContinuumPage))
        continuum_button.pack()

        #Click this button to measure the STIS emission lines
        STIS_button = ttk.Button(self, text="Measure STIS Emission Lines",
                                 command=lambda: controller.show_frame(STISPage))
        STIS_button.pack()

        global spectrum_entry_label, spectrum_entry_box, spectrum_entry_submit

        #####ENTRY BOX FOR USER TO ENTER FILE NAME OF SPECTRUM WITH LINES TO FIT#######################
        spectrum_entry_label = ttk.Label(self, text="Choose a spectrum: ")
        spectrum_entry_label.pack(side=tk.TOP, pady=10)

        #Get list of lines to fit when user hits 'Enter'
        spectrum_entry_submit = ttk.Button(self, text="Select", command=lambda: getfilename())
        spectrum_entry_submit.pack(side=tk.TOP, pady=10)

        #Click this button to reset the GUI - it will close the old window and open a new instance
        reset_button = ttk.Button(self, text="Reset",
                                  command=lambda: reset())
        reset_button.pack()
        reset_button.place(x=1175, y=0)

        #Define lists here to append line-fitting parameters (will eventually put in a dataframe, but it's more computationally expensive to append one row at a time to a dataframe)
        global model_params, model_errors, fit_flag, chi2, DOF
        model_params = []
        model_errors = []
        fit_flag = []
        chi2 = []
        DOF = []

        ##########METHOD TO READ IN USER INPUT FILENAME - THIS/THE FOLLOWING METHODS ARE STORED UNDER THE __init__ METHOD BECAUSE I WANT THEM TO RUN EVERYTIME THIS PAGE IS OPENED ###########################################
        def getfilename():

            #Read in filename that user entered - should add something here to prompt another filename entry if an invalid string is entered
            spectrum_filename = askopenfilename()
            #Destroy label
            spectrum_entry_label.destroy()
            #Destroy "Select" button for this text box
            spectrum_entry_submit.destroy()

            global param_file_label
            #Prompt the user to choose a file with a list of lines to fit
            param_file_label = ttk.Label(self, text="Select a file with lines to fit: ")
            param_file_label.pack(side=tk.TOP, pady=10)

            global spectrum_df
            try:
                spectrum_df = pd.read_csv(spectrum_filename) #Will work with this file unless the GUI is reset and a new one is input
                spectrum_df.columns = ['wavelength', 'flux', 'fluxerr']
            except:
                #Add option to read in different types of files, since delimiter needs to be specified for pd.read_csv()
                potential_separators = [' ', '\t', ';', ':']
                separator_idx = 0
                while len(list(spectrum_df)) < 3 and separator_idx < len(potential_separators):
                    spectrum_df = pd.read_csv(spectrum_filename, sep=potential_separators[separator_idx])
                    separator_idx += 1
                try:
                    spectrum_df.columns = ['wavelength', 'flux', 'fluxerr']
                except:
                    #Print error message to GUI window if file still isn't read correctly
                    if len(list(spectrum_df)) < 3:
                        read_fail_label = ttk.Label(self, text="Couldn't read your file... hit 'Reset' to start over!")
                        read_fail_label.pack(side=tk.TOP)
                        read_fail_label.place(bordermode=tk.OUTSIDE, x=12.5, y=150)
                        read_fail_label.config(font=("Courier", 32), foreground="red")

            global param_file_submit
            param_file_submit = ttk.Button(self, text="Select", command=lambda: getlinelist())
            param_file_submit.pack(side=tk.TOP, pady=0)

            #####INITIALIZE INTERACTIVE PLOT FOR SPECTRUM########################################
            #First, plot full spectrum - user will be able to interact with it before fitting emission lines, but should add option to go back to this
            global f, ax1, canvas
            f = Figure()
            ax1 = f.add_subplot(111)
            ax1.plot(spectrum_df['wavelength'],
                     spectrum_df['flux'],
                     color="k")
            ax1.set_xlabel("Wavelength")
            ax1.set_ylabel("Flux")
            canvas = FigureCanvasTkAgg(f, self)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP) #Outputs plot to GUI window

            toolbar = NavigationToolbar2Tk(canvas, self) #Initializes matplotlib toolbar
            toolbar.update()
            canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=False)

        ###########METHOD THAT ALLOWS USER TO REFIT EMISSION LINE - DOESN'T WORK IF CLICKED BEFORE FIRST LINE IS FIT#######
        def reduce_count():

            global line_idx, ax1, test_wave

            ax1.clear()
            #"line_idx" is the number of the line in the text file - decrementing it will cause the line to pop up again when the user clicks "Fit next line" button
            line_idx = line_idx-1
            #Replace last element in list with a flag to indicate that this was a bad fit - bad fit parameters will still be output to text file at the end
            fit_flag[-1] = 1

            global redonelabel
            redonelabel = ttk.Label(self, text="Line reset! Select 'Fit next line:' to continue.") #Label appears in top left corner of GUI to guide user to click the "Fit next line" button
            redonelabel.place(bordermode=tk.OUTSIDE, x=12.5, y=150)
            redonelabel.config(font=("Courier", 32), foreground="blue")

        ##############READS IN FILE WITH LIST OF LINES TO FIT#############################################
        def getlinelist():

            param_file = askopenfilename()
            param_file_label.destroy()
            param_file_submit.destroy()

            global line_param_df
            try:
                line_param_df = pd.read_csv(param_file) #Will work with this file unless the GUI is reset and a new one is input
                line_param_df.columns = ['Line ID', 'Central Wavelength']
            except:
                #Add option to read in different types of files, since delimiter needs to be specified for pd.read_csv()
                potential_separators = [' ', '\t', ';', ':']
                separator_idx = 0
                while len(list(spectrum_df)) != 2 and separator_idx < len(potential_separators):
                    line_param_df = pd.read_csv(param_file, sep=potential_separators[separator_idx])
                    separator_idx += 1
                try:
                    line_param_df.columns = ['Line ID', 'Central Wavelength']
                except:
                    #Print error message to GUI window if file still isn't read correctly
                    if len(list(line_param_df)) != 2:
                        read_fail_label = ttk.Label(self, text="Couldn't read your file... hit 'Reset' to start over!")
                        read_fail_label.pack(side=tk.TOP)
                        read_fail_label.place(bordermode=tk.OUTSIDE, x=12.5, y=150)
                        read_fail_label.config(font=("Courier", 32), foreground="red")

            global fit_line_submit
            fit_line_submit = ttk.Button(self, text="Fit Next Line", command = lambda: fitline()) #This button will be up until the last line is fit - users can skip lines, but other labels won't be destroyed if they do
            fit_line_submit.place(bordermode=tk.OUTSIDE, x=500, y=625)

            global submitback
            submitback = ttk.Button(self, text="Redo Fit", command = lambda: reduce_count())
            submitback.place(bordermode=tk.OUTSIDE, x=700, y=625)

        ############GET RIGHT-HAND-SIDE OF FWHM FOR SIGMA INITIAL CONDITION##################################
        def onclick_FWHMright(event3):

            global FWHM_right
            #Get x- and y-coordinates (wavelength and flux) of user click called "event3"
            FWHM_right = event3.xdata
            FWHM_flux_right = event3.ydata

            global cid
            #Disconnect from canvas after user click
            canvas.mpl_disconnect(cid)

            #Calculate Gaussian sigma = FWHM / 2*sqrt(2*ln(2))
            sigma = np.absolute(FWHM_right - FWHM_left) / (2.*np.sqrt(2.*np.log(2.)))

            global a_params
            #Add Gaussian sigma initial condition to parameter array
            a_params.append(sigma)

            #Flux of "FWHM" is greater than amplitude -> absorption line, set negative amplitude
            if FWHM_flux_right > a_params[-3]:
                a_params[-3] = a_params[-3]*-1.

            global right_FWHM_label
            #Get rid of label prompting user to click right side of line
            right_FWHM_label.destroy()

        ############GET LEFT-HAND-SIDE OF FWHM FOR SIGMA INITIAL CONDITION##################################
        def onclick_FWHMleft(event2):

            global FWHM_left
            #Get x-coordinate (wavelength) of user click called "event2"
            FWHM_left = event2.xdata

            global left_FWHM_label
            #Get rid of label prompting user to click left side of line
            left_FWHM_label.destroy()

            global cid
            #Disconnect from canvas after user click
            canvas.mpl_disconnect(cid)

            global right_FWHM_label
            #Add prompt to have user click right side of line - will appear in top left corner of GUI
            right_FWHM_label = ttk.Label(self, text="Click the right intersection of the red dashed line and the data:")
            right_FWHM_label.place(bordermode=tk.OUTSIDE, x=0, y=0)

        ############GET INITIAL CONDITION FOR AMPLITUDE OF LINE###########################################
        def onclick_waveamp(event1):

            global central_wavelength, amplitude
            central_wavelength = event1.xdata #x-coordinate (wavelength) from user click called "event1"
            amplitude = event1.ydata #y-coordinate (peak flux of emission line) from user click called "event1"
            a_params.append(amplitude) #Insert initial conditions for amplitude/center of emission line into parameter array
            a_params.append(central_wavelength)

            global param_label
            param_label.destroy() #Destroy label prompting user to click peak of emission line

            global cid
            canvas.mpl_disconnect(cid) #Disconnect from canvas after user click

            #Plot red, dashed line at half of the maximum amplitude to show user where to click for FWHM guess
            ax1.hlines(amplitude / 2.0,
                       central_wavelength - (0.5*linespan),
                       central_wavelength + (0.5*linespan),
                       color = "red", linestyle = "dashed")
            canvas.draw()

            global left_FWHM_label
            left_FWHM_label = ttk.Label(self, text="Click the left intersection of the red dashed line and the data:") #Prompt user to click left side of line - will appear at top left corner of GUI
            left_FWHM_label.pack(side=tk.TOP)
            left_FWHM_label.place(bordermode=tk.OUTSIDE, x=0, y=0)

        ##########DETERMINES HOW MANY EMISSION LINES IN THE REGION OF THE CENTRAL LINE THE USER WANTS TO FIT################
        def numlines():

            global numlines_label, numlines_entry, numlines_bttn

            answer = numlines_entry.get() #Get user input for the number of lines to fit
            numlines = int(answer)

            numlines_entry.destroy() #Destroy prompt
            numlines_label.destroy()
            numlines_bttn.destroy()

            global a_params
            a_params = [] #Define list of parameters for line fitting (amplitude, central wavelength, sigma for each Gaussian as well as continuum fit parameters) - will be passed as initial conditions for fit function
            count = 1 #Get initial conditions for line 1 (must be central line)
            while count <= numlines: #Execute until initial conditions have been acquired for all lines the user wants to fit
                global param_label
                param_label = ttk.Label(self, text="Click on the peak of an emission line:") #Prompt user to select initial condition for amplitude of Gaussian - will appear in upper left corner of GUI
                param_label.pack(side=tk.TOP)
                param_label.place(bordermode=tk.OUTSIDE, x=0, y=0)

                global cid
                cid = canvas.mpl_connect("button_press_event", lambda event1: onclick_waveamp(event1)) #Connect mouse to canvas so user can click to select amplitude
                tk.Frame.wait_window(param_label) #Wait until this is destroyed in onclick_FWHMright before going to next loop
                global left_FWHM_label
                cid = canvas.mpl_connect("button_press_event", lambda event2: onclick_FWHMleft(event2)) #Connect mouse to canvas so user can click on left side of emission line
                tk.Frame.wait_window(left_FWHM_label)

                global right_FWHM_label
                cid = canvas.mpl_connect("button_press_event", lambda event3: onclick_FWHMright(event3)) #Connect mouse to canvas so user can click on right side of emission line
                tk.Frame.wait_window(right_FWHM_label)

                count += 1 #Ready to get initial conditions for the next line

            #ADD INITIAL CONDITIONS FOR FIRST ORDER POLYNOMIAL FIT TO CONTINUUM
            a_params.append(0.)
            a_params.append(1.e-14)
            #Convert list of parameters to numpy array to make it easier to index over
            a_params = np.array(a_params)

            global test_wave
            #Must put central wavelength from line list file here as initial condition for central line - otherwise, later code will extract variable number of data points for fitting depending on where the user clicked.
            a_params[1] = test_wave

            ltp = "LTP0" #most of the time
            polyorder = 1

            a_wavelengths = a_params[range(1, len(a_params)-(polyorder+1), 3)]
            if len(a_wavelengths) > 1:
                a_central = min(a_wavelengths) + 0.5*(max(a_wavelengths) - min(a_wavelengths))
            else:
                a_central = a_wavelengths[0]

            subset_df = spectrum_df.loc[(spectrum_df['wavelength'] >= a_central - 0.5*linespan) &
                                        (spectrum_df['wavelength'] <= a_central + 0.5*linespan)]
            #Use Python version of Kevin's multi-Gaussian fitting code to convolve Gaussians with COS LSFs - return 3 Angstrom range of wavelengths/corresponding fluxes/errors, model based on optimal parameters, optimal parameters, and errors on model fits
            f_model, popt, model_err = lff.fit_model(subset_df['wavelength'],
                                                     subset_df['flux'],
                                                     subset_df['fluxerr'],
                                                     a_params, numlines, polyorder, ltp)
            #Add model emission line fluxes to subset dataframe
            subset_df['Model Flux'] = f_model

            global model_params, model_errors
            model_params.append(popt) #Add optimal parameters/errors to lists - will output to text file later in code
            model_errors.append(model_err)
            print (model_params)

            amplitude_indices = list(range(0, len(popt), 3)) #Extract amplitudes of all Gaussians fit to set y-axis limit on plots
            line_amplitudes = popt[amplitude_indices]

            global chi2, DOF
            #DOF = number of points fit - number of parameters used to fit
            dof_line = len(subset_df) - len(popt)
            #Get rid of error bars < 10^-16 or so because of problem with COS reduction pipeline
            min_flux_err = 1.e-17
            #Calculate reduced chi2 statistic
            test_stat = np.sum(np.square((subset_df.loc[subset_df['fluxerr'] >= min_flux_err]['flux'] - subset_df.loc[subset_df['fluxerr'] >= min_flux_err]['Model Flux']) / subset_df.loc[subset_df['fluxerr'] >= min_flux_err]['fluxerr'])) / float(dof_line)

            DOF.append(dof_line)
            chi2.append(test_stat)
            fit_flag.append(0) #Add flag of 0 to indicate that the fit was a good one

            global lineID

            #############PLOT DATA AND MODEL FIT - WILL SAVE FILE TO DIRECTORY, BUT SHOULD GENERALIZE THAT PART###########################
            ax1.clear()
            ax1.set_title("%s: %s Angstroms" % (str(lineID), str(test_wave)))
            ax1.plot(subset_df['wavelength'],
                     subset_df['flux'],
                     color="k", lw=3.)
            ax1.plot(subset_df['wavelength'],
                     f_model,
                     color=(0., 0.5, 0.5), lw=3.)
            ax1.errorbar(subset_df['wavelength'],
                         subset_df['flux'],
                         yerr = subset_df['fluxerr'],
                         color = "red", linestyle = "none")
            for j in range(0, int(answer)): #Plot Gaussians
                z2 = (subset_df['wavelength'] - popt[1+3*j]) / popt[2+3*j]
                gauss = popt[0+3*j]*np.exp(-np.square(z2) / 2.)
                ax1.plot(subset_df['wavelength'],
                         gauss,
                         color = "purple", linestyle = "dashed")
            ax1.set_xlabel(r"Wavelength $\left[ \AA \right]$")
            increment = (max(subset_df['wavelength']) - min(subset_df['wavelength'])) / 5.
            xtick_array = np.arange(min(subset_df['wavelength']), max(subset_df['wavelength']), increment)
            ax1.set_xticks(xtick_array)
            ax1.set_xticklabels([str(tick)[0:6] for tick in xtick_array])
            ax1.set_xlim(min(subset_df['wavelength']), max(subset_df['wavelength']))
            ax1.set_ylabel("Flux \n" + r"[erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$]")
            if max(line_amplitudes) > np.amax(f_model):
                ax1.set_ylim(min(subset_df['flux']), max(subset_df['flux']) + 0.1e-14)
            else:
                ax1.set_ylim(min(subset_df['flux']), max(f_model) + 0.3e-14) #Need to keep top of error bars on plot
            canvas.draw()
            f.savefig(filename, format = "pdf")

        #################DETERMINES X-COORDINATE (WAVELENGTH) OF LEFT SIDE OF LINE AT FWHM############################################
        def getleftside(event3):

            global continuum_left
            continuum_left = event3.xdata #Save x-coordinate of button click event (event3)

            global guide_label
            guide_label.destroy() #Destroy label prompting user to click left side of line

            global cid
            canvas.mpl_disconnect(cid) #Disconnect mouse from canvas

        ##############FUNCTION TO GIVE FITTING ALGORITHM IN SCIPY.OPTIMIZE###########################################################
        def continuum_func(wavelength, slope, intercept):
            return intercept + (wavelength*slope)

        #################ALLOWS USER TO PLACE UPPER LIMIT ON EMISSION LINE IF NOT VISIBLE############################################
        def pickcontinuum():

            global continuumfit_label, continuumfit_bttn
            continuumfit_label.destroy() #Destroy label and button directing user to select region of continuum to fit
            continuumfit_bttn.destroy()

            global guide_label
            guide_label = ttk.Label(self, text="Select the bluest wavelength in the clean region:") #Prompts user to select left side of emission line at FWHM
            guide_label.place(bordermode=tk.OUTSIDE, x=0, y=0) #Label will be placed in upper left corner of GUI

            global cid
            cid = canvas.mpl_connect("button_press_event", lambda event3: getleftside(event3)) #Connect to canvas to get user click location
            tk.Frame.wait_window(guide_label)

            continuum_right = continuum_left + 3. #3 Angstrom region to fit continuum (can make size of region user input)

            continuum_subset_df = spectrum_df.loc[(spectrum_df['wavelength'] >= continuum_left) &
                                        (spectrum_df['wavelength'] <= continuum_right)]

            popt_continuum, pcov_continuum = sciopt.curve_fit(continuum_func,
                                                              continuum_subset_df['wavelength'],
                                                              continuum_subset_df['flux'],
                                                              p0=[0., 1.e-14]) #Fit continuum function to data

            global test_wave
            spectrum_subset_df = spectrum_df.loc[(spectrum_df['wavelength'] >= test_wave - 0.5*linespan) &
                                                 (spectrum_df['wavelength'] <= test_wave + 0.5*linespan)]
            spectrum_subset_df['Calculated Continuum Fluxes'] = popt_continuum[1] + (spectrum_subset_df['wavelength']*popt_continuum[0])
            #Subtract continuum before calculating RMS flux - should have option to fit other polynomial orders for continuum
            spectrum_subset_df['Continuum Subtracted Fluxes'] = spectrum_subset_df['flux'] - spectrum_subset_df['Calculated Continuum Fluxes']

            flux_RMS = np.sqrt(np.mean(np.square(spectrum_subset_df['Continuum Subtracted Fluxes'])))
            flux_average = np.average(spectrum_subset_df['Continuum Subtracted Fluxes']) #Should be ~0
            st_dev = np.std(spectrum_subset_df['Continuum Subtracted Fluxes'])

            global model_params, model_errors
            popt = np.zeros(3) + np.nan #Won't have model parameters, but still need to append to arrays for output file
            popt[0] = flux_RMS
            popt[1] = test_wave
            model_params.append(popt)
            errors = np.zeros(3) + np.nan
            errors[0] = st_dev
            model_errors.append(errors)

            global chi2, DOF
            chi2.append(np.nan) #No model to calculate chi^2 for
            DOF.append(0)
            fit_flag.append(0)

            lower_flux_lim = min(spectrum_subset_df['Continuum Subtracted Fluxes'])
            upper_flux_lim = max(spectrum_subset_df['Continuum Subtracted Fluxes'])

            #Saved figure will have blue horizontal line at value of RMS flux
            ax1.clear()
            ax1.set_title(' '.join([lineID, '\n', str(test_wave)[0:6], 'Angstroms']))
            ax1.plot(spectrum_subset_df['wavelength'],
                     spectrum_subset_df['Continuum Subtracted Fluxes'],
                     color="k", lw=3.)
            ax1.hlines(flux_RMS,
                       min(spectrum_subset_df['wavelength']),
                       max(spectrum_subset_df['wavelength']),
                       color = "blue", linestyle="solid", lw=3.)
            ax1.set_xlabel(r"Wavelength $\left[ \AA \right]$")
            increment = (max(spectrum_subset_df['wavelength']) - min(spectrum_subset_df['wavelength'])) / 5.
            xtick_array = np.arange(min(spectrum_subset_df['wavelength']), max(spectrum_subset_df['wavelength']), increment)
            ax1.set_xticks(xtick_array)
            ax1.set_xticklabels([str(tick)[0:6] for tick in xtick_array])
            ax1.set_xlim(min(spectrum_subset_df['wavelength']), max(spectrum_subset_df['wavelength']))
            ax1.set_ylim(lower_flux_lim-0.1*lower_flux_lim, upper_flux_lim+0.1*upper_flux_lim)
            ax1.set_ylabel("Flux \n" + r"[erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$]")
            canvas.draw()
            f.savefig(filename, format = "pdf")

        ###########################ASK USER WHETHER OR NOT THE EMISSION LINE IS THERE TO FIT########################################
        def tofitornot():

            global bright_entry, bright_label, bright_bttn

            answer = bright_entry.get() #User will be asked whether or not to fit the line, should answer "yes" or "no"
            bright_entry.destroy()
            bright_label.destroy()
            bright_bttn.destroy()

            if answer in ["yes", "YES", "y", "Y", "True", "T", "TRUE"]:

                global numlines_label, numlines_entry, numlines_bttn

                numlines_label = ttk.Label(self, text="Enter the number of lines to fit:") #Proceed to multi-Gaussian fit
                numlines_label.place(bordermode=tk.OUTSIDE, x=0, y=0)

                numlines_entry = ttk.Entry(self)
                numlines_entry.place(x=0, y=30)

                numlines_bttn = ttk.Button(self, text="Enter", command = lambda: numlines())
                numlines_bttn.place(x=0, y=60)

            else: #Should add something to indicate if filename is bad...

                global continuumfit_label, continuumfit_bttn
                continuumfit_label = ttk.Label(self, text="Navigate to a nearby region free of emission lines, then click 'Enter':") #Fit continuum to put upper limit on RMS flux
                continuumfit_label.place(bordermode=tk.OUTSIDE, x=0, y=0)

                continuumfit_bttn = ttk.Button(self, text="Enter", command = lambda: pickcontinuum())
                continuumfit_bttn.place(x=0, y=30)

        global line_idx #Counter to keep track of lines in user-provided list of lines to fit (e.g. France2012_Table2.csv)
        line_idx = 0

        global amplitudes, amp_err, observed_wavelengths, central_err, sigma, sigma_err #Define lists to hold Gaussian properties
        amplitudes = []
        amp_err = []
        observed_wavelengths = []
        central_err = []
        sigma = []
        sigma_err = []

        global redonelabel
        redonelabel = ttk.Label(self, text="Line reset! Select 'Fit next line:' to continue.") #Will display if user clicks "Redo Fit" (add next/previous line buttons)

        def getlinespan():

            global linespan_label, linespan_entry, linespan_bttn, linespan

            answer = linespan_entry.get() #Get user input for the number of lines to fit
            linespan = float(answer)

            linespan_entry.destroy() #Destroy prompt
            linespan_label.destroy()
            linespan_bttn.destroy()

            lower_limit = test_wave - (linespan / 2.)
            upper_limit = test_wave + (linespan / 2.)

            subset_df = spectrum_df.loc[(spectrum_df['wavelength'] >= lower_limit) &
                                        (spectrum_df['wavelength'] <= upper_limit)]

            increment = (upper_limit - lower_limit) / 5. #5 tick marks seems reasonable for these plots
            xtick_array = np.arange(lower_limit, upper_limit+increment, increment)
            xtick_labels = [str(tick)[0:6] for tick in xtick_array]

            #Redo plot to reflect different linespan
            lowerfluxlim = min(subset_df['flux'])
            upperfluxlim = max(subset_df['flux'])
            #Plot region around central wavelength of line so user can decide whether or not to fit line
            ax1.clear()
            ax1.set_title(' '.join([lineID, '\n', str(test_wave)[0:6], 'Angstroms']))
            ax1.plot(spectrum_df['wavelength'],
                     spectrum_df['flux'],
                     color="k", lw=3.) #Need to plot full spectrum so user can pan around
            ax1.set_xlabel(r"Wavelength $\left[ \AA \right]$")
            ax1.set_xticks(xtick_array)
            ax1.set_xticklabels(xtick_labels)
            ax1.set_xlim(lower_limit, upper_limit)
            ax1.set_ylabel("Flux \n" + r"[erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$]")
            ax1.set_ylim(lowerfluxlim-0.1*lowerfluxlim, upperfluxlim+0.1*upperfluxlim)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=False)

        ################INCREMENT line_idx TO RUN ALL METHODS ABOVE ON THE NEXT EMISSION LINE IN THE LIST#######################################
        def fitline():

            global line_idx

            if line_idx >= len(line_param_df['Line ID']): #Reached end of line list -> destroy all plotting/prompting widgets
                fit_line_submit.destroy()
                canvas.get_tk_widget().destroy()
                canvas._tkcanvas.destroy()
                submitback.destroy()

                global amplitudes, amp_err, observed_wavelengths, central_err, sigma, sigma_err, chi2, DOF, fit_flag, final_lab_wavelengths

                final_lab_wavelengths = []
                chi2index = 0
                for k in range(0, np.shape(model_params)[0]): #Go through each set of parameters in model_params and extract properties for calculations

                    param_set = model_params[k]
                    err_set = model_errors[k]

                    #Fix the "2" when you set polynomial order as user input
                    wavelength_indices = list(range(1, len(param_set)-2, 3))

                    if len(wavelength_indices) == 0: #Didn't fit line, calculated RMS flux instead

                        final_lab_wavelengths.append(param_set[1]) #This is just the lab wavelength
                        amplitudes.append(np.nan)
                        amp_err.append(np.nan)
                        observed_wavelengths.append(np.nan)
                        central_err.append(np.nan)
                        sigma.append(np.nan)
                        sigma_err.append(np.nan)

                    else: #May need to add multiple lines within +/-50 km/s of the central wavelength (e.g. lines with narrow/broad components)
                        line_count = 0
                        for loc in wavelength_indices:
                            wave_totest = param_set[loc]
                            line_param_df["Diff"] = np.absolute(wave_totest - line_param_df['Central Wavelength']) #Find closest matching lab wavelength in list
                            closest_match = line_param_df["Diff"].idxmin()

                            test_velocity = np.absolute(wave_totest - line_param_df['Central Wavelength'].iloc[closest_match])*C_kms / line_param_df['Central Wavelength'].iloc[closest_match] #Figure out if velocity is close to central wavelength

                            if test_velocity <= 1000.: #+/- 50 km/s for H2 lines - otherwise want to return all values
                                final_lab_wavelengths.append(line_param_df['Central Wavelength'].iloc[closest_match])

                                amplitudes.append(param_set[loc-1])
                                amp_err.append(err_set[loc-1])

                                observed_wavelengths.append(param_set[loc])
                                central_err.append(err_set[loc])

                                sigma.append(param_set[loc+1])
                                sigma_err.append(err_set[loc+1])

                                if line_count > 0:

                                    if chi2index < np.size(chi2):
                                        chi2 = np.insert(chi2, chi2index, chi2[chi2index-line_count])
                                        DOF = np.insert(DOF, chi2index, DOF[chi2index-line_count])
                                        fit_flag = np.insert(fit_flag, chi2index, fit_flag[chi2index-line_count])

                                    else:
                                        chi2 = np.append(chi2, chi2[-1])
                                        DOF = np.append(DOF, DOF[-1])
                                        fit_flag = np.append(fit_flag, fit_flag[-1])

                                line_count += 1

                            chi2index += 1

                #Convert all lists to numpy arrays/make sure values are stored as double precision floats
                fit_param_df = pd.DataFrame({'Amplitude': np.array(amplitudes).astype(np.float64),
                                              'Amplitude Uncertainty': np.array(amp_err).astype(np.float64),
                                              'Observed Wavelength': np.array(observed_wavelengths).astype(np.float64),
                                              'Observed Wavelength Uncertainty': np.array(central_err).astype(np.float64),
                                              'Sigma': np.array(sigma).astype(np.float64),
                                              'Sigma Uncertainty': np.array(sigma_err).astype(np.float64),
                                              'Chi^2': np.array(chi2).astype(np.float64),
                                              'DOF': np.array(DOF),
                                              'fit_flag': np.array(fit_flag),
                                              'Lab Wavelengths': np.array(final_lab_wavelengths).astype(np.float64).flatten()
                                              })

                #Calculate properties for all lines/fits
                fit_param_df['Radial Velocity'] = (fit_param_df['Observed Wavelength'] - fit_param_df['Lab Wavelengths'])*C_kms / fit_param_df['Lab Wavelengths']
                fit_param_df['Radial Velocity Uncertainty'] = fit_param_df['Observed Wavelength Uncertainty']*C_kms / fit_param_df['Lab Wavelengths']
                fit_param_df['FWHM'] = C_kms*2.*fit_param_df['Sigma']*np.sqrt(2*np.log(2.)) / fit_param_df['Observed Wavelength'] #km/s
                fit_param_df['FWHM Uncertainty'] = np.sqrt(np.square(fit_param_df['Sigma Uncertainty']*C_kms*2.*np.sqrt(2.*np.log(2.)) / fit_param_df['Observed Wavelength']) + np.square(fit_param_df['Observed Wavelength Uncertainty']*C_kms*2.*fit_param_df['Sigma']*np.sqrt(2.*np.log(2.)) / np.square(fit_param_df['Observed Wavelength'])))
                fit_param_df['Integrated Line Flux'] = fit_param_df['Amplitude']*fit_param_df['Sigma']*np.sqrt(2.*np.pi)
                fit_param_df['Integrated Line Flux Uncertainty'] = np.sqrt( np.square(fit_param_df['Amplitude Uncertainty']*fit_param_df['Sigma']*np.sqrt(2.*np.pi)) + np.square(fit_param_df['Sigma Uncertainty']*fit_param_df['Amplitude']*np.sqrt(2.*np.pi)))

                #Calculated integrated line flux with amplitudes for lines that were fit, so values for upper limits will be NaN - now insert RMS flux instead
                for row in range(0, len(fit_param_df)):
                    if np.isnan(fit_param_df['Integrated Line Flux'].iloc[row]) == True:
                        for col in list(fit_param_df):
                            fit_param_df[col].iloc[row] = 99.999
                        fit_param_df['Integrated Line Flux'].iloc[row] = model_params[row][0] #Append RMS flux
                        fit_param_df['Integrated Line Flux Uncertainty'].iloc[row] = model_errors[row][0]

                #Output all arrays to a text file! Should make output filename user input eventually
                fit_param_df.to_csv('SELFiE_line_fit_params.csv', sep=',', index=False)

                global bye_label

                bye_label = ttk.Label(self, text="Line fit parameters saved to file!")
                bye_label.pack(side=tk.TOP)
                bye_label.place(bordermode=tk.OUTSIDE, x=300, y=300)
                bye_label.configure(font=("Courier", 32), foreground="blue")

            else: #Keep fitting lines!

                ax1.clear()

                global redonelabel
                redonelabel.destroy()

                global lineID, test_wave, lower_limit, upper_limit
                lineID = line_param_df['Line ID'][line_idx]
                test_wave = line_param_df['Central Wavelength'][line_idx]

                #####BUTTON TO ASK FOR USER INPUT LINESPAN##########################
                global linespan_label, linespan_entry, linespan_bttn, linespan

                linespan_label = ttk.Label(self, text="Change wavelength span of line to fit:")
                linespan_label.place(bordermode=tk.OUTSIDE, x=850, y=0) #Should sit right under "Reset" button

                linespan_entry = ttk.Entry(self)
                linespan_entry.insert(0, "3") #Set default value to 3 Angstroms, will be converted to float by getlinespan() function
                linespan_entry.place(x=870, y=30)

                linespan_bttn = ttk.Button(self, text="Enter", command = lambda: getlinespan())
                linespan_bttn.place(x=915, y=60)

                #Default, if the user doesn't click the button to change the linespan
                linespan = 3.

                lower_limit = test_wave - 0.5 #0.5*linespan
                upper_limit = test_wave + 0.5 #0.5*linespan

                #Generate 5 tick mark locations and labels for plot
                increment = (upper_limit - lower_limit) / 5.
                xtick_array = np.arange(lower_limit, upper_limit+increment, increment)
                xtick_labels = [str(tick)[0:6] for tick in xtick_array]

                global filename #This shouldn't really be defined here...
                filename = ''.join([lineID.replace(" ", ""), "_", str(test_wave)[0:str(test_wave).find('.')], ".pdf"])

                #Set y-axis limits so default setting doesn't accommodate height of LyA profile
                lowerfluxlim = min(spectrum_df.loc[(spectrum_df['wavelength'] >= lower_limit) &
                                                   (spectrum_df['wavelength'] <= upper_limit)]['flux'])
                upperfluxlim = max(spectrum_df.loc[(spectrum_df['wavelength'] >= lower_limit) &
                                                   (spectrum_df['wavelength'] <= upper_limit)]['flux'])
                #Plot region around central wavelength of line so user can decide whether or not to fit line
                ax1.set_title(' '.join([lineID, '\n', str(test_wave)[0:6], "Angstroms"]))
                ax1.plot(spectrum_df['wavelength'],
                         spectrum_df['flux'],
                         color="k", lw=3.)
                ax1.set_xlabel(r"Wavelength $\left[ \AA \right]$")
                ax1.set_xticks(xtick_array)
                ax1.set_xticklabels(xtick_labels)
                ax1.set_xlim(lower_limit, upper_limit)
                ax1.set_ylabel("Flux \n" + r"[erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$]")
                ax1.set_ylim(lowerfluxlim-0.1*lowerfluxlim, upperfluxlim+0.1*upperfluxlim)
                canvas.draw()
                canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=False)

                global bright_label, bright_entry, bright_bttn

                bright_label = ttk.Label(self, text="Fit central line?") #Set up widgets to guide user through fitting process
                bright_label.pack(side=tk.TOP)
                bright_label.place(bordermode=tk.OUTSIDE, x=0, y=0)

                bright_entry = ttk.Entry(self)
                bright_entry.pack(side=tk.TOP)
                bright_entry.place(x=0, y=30)

                bright_bttn = ttk.Button(self, text="Enter", command = lambda: tofitornot())
                bright_bttn.pack(side=tk.TOP)
                bright_bttn.place(x=0, y=60)

                line_idx += 1
                redonelabel = ttk.Label(self, text="Line reset! Select 'Fit next line:' to continue.")

######ACTUALLY RUN THE APP!###################################
print ("Launching SELFiE!")
app = SELFiE() #From tkinter
app.geometry("1280x720")
app.mainloop()
