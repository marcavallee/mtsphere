
import sys
import os
import time
import numpy as np
from scipy.constants import pi, mu_0
import tkinter as tk
import matplotlib
matplotlib.use('tkAgg')
from layeredearthfunctions import PlaneWaveImpedance, apparentresistivity
from mtsphere3d import mtsphere3d, getimpedance
from mtplot import MTPlot
from tkinter import filedialog,messagebox
from waitwindow import WaitWindow
import threading

developer = 'Developed by: Marc A. Vallée'
company = 'for: Geo Data Solutions GDS Inc.'

class Layer:

    def __init__(self,nlyr):
        self.nlyr = nlyr
        self.resistivity = 0.001
        self.thickness = 10
        self.resistivityvariable = tk.DoubleVar()
        self.thicknessvariable = tk.DoubleVar()

    def displayline(self,frame,row,i):
        label = tk.Label(frame, text='Resistivity '+str(i+1))
        label.grid(row=row, sticky=tk.W)
        entry = tk.Entry(frame, textvariable=self.resistivityvariable )
        entry.grid(row=row, column=1)
        if i < self.nlyr-1:
            label = tk.Label(frame, text='Thickness '+str(i+1))
            label.grid(row=row, column=2)
            entry = tk.Entry(frame, textvariable=self.thicknessvariable )
            entry.grid(row=row, column=3)
        self.refresh()
        self.resistivityvariable.trace('w',self.updateresistivity)
        self.thicknessvariable.trace('w',self.updatethickness)

    def refresh(self):
        self.resistivityvariable.set(self.resistivity)
        self.thicknessvariable.set(self.thickness)

    def updateresistivity(self, *args):
        try:
            self.resistivity = self.resistivityvariable.get()
        except:
            return

    def updatethickness(self, *args):
        try:
            self.thickness = self.thicknessvariable.get()
        except:
            return

class pyMTSphere:

    def __init__(self):
        self.root = tk.Tk()
        self.mainframe = tk.Frame()
        self.root.title("pyMTSphere")
        self.directory = ''
        self.directoryvariable = tk.StringVar()
        self.title = 'MTSphere Model'
        self.titlevariable = tk.StringVar()
        self.outfile = None
        self.mf1file = None
        self.mf2file = None
        self.pxyd = np.array([])
        self.version = '2.2.2'
        self.releasedate = 'November 2026'
        self.htarg = {'dlf':'key_401_2009'}
        self.mtplot = MTPlot(self)
        self.diag = 0
        self.nterms = 6
        self.nlyr = 1
        self.minfreq = 1
        self.maxfreq = 1000
        self.nf = 1
        self.logarithmic = 1
        self.radius = 50
        self.depth = 100
        self.sphres = 1
        self.nx = 41
        self.ny = 1
        self.nd = 1
        self.xmin = -400
        self.xmax = 400
        self.ymin = 0
        self.ymax = 0
        self.dmin = 0
        self.dmax = 0
        self.layerlist = [Layer(self.nlyr)]
        self.ntermsvariable = tk.IntVar()
        self.minfreqvariable = tk.DoubleVar()
        self.maxfreqvariable = tk.DoubleVar()
        self.nfvariable = tk.IntVar()
        self.diagvariable = tk.IntVar()
        self.logarithmicvariable = tk.IntVar()
        self.nlyrvariable = tk.IntVar()
        self.radiusvariable = tk.DoubleVar()
        self.depthvariable = tk.DoubleVar()
        self.sphresvariable = tk.DoubleVar()
        self.nxvariable = tk.IntVar()
        self.nyvariable = tk.IntVar()
        self.ndvariable = tk.IntVar()
        self.xminvariable = tk.DoubleVar()
        self.xmaxvariable = tk.DoubleVar()
        self.yminvariable = tk.DoubleVar()
        self.ymaxvariable = tk.DoubleVar()
        self.dminvariable = tk.DoubleVar()
        self.dmaxvariable = tk.DoubleVar()

    def displaymainframe(self):
        self.mainframe.destroy()
        self.mainframe = tk.Frame(self.root)
        row = 0
        self.mainframe.grid(row=row)
        title = tk.Label(self.mainframe, text="MTSphere, developed by Marc A. Vallée, Geo Data Solutions GDS, Inc.", font=("Times", 10), relief=tk.RIDGE)
        title.grid(row=row)
        row = row + 1
        directoryframe = tk.Frame(self.mainframe)
        directoryframe.grid(row=row)
        self.displaydirectoryframe(directoryframe)
        row = row + 1
        parameterframe = tk.Frame(self.mainframe)
        parameterframe.grid(row=row)
        self.displayparameterframe(parameterframe)
        row = row + 1
        commandframe = tk.Frame(self.mainframe)
        commandframe.grid(row=row, sticky=tk.W)
        self.displaycommandframe(commandframe)

    def displaydirectoryframe(self,directoryframe):
        row = 0
        label = tk.Label(directoryframe, text='Directory')
        label.grid(row=row, sticky=tk.W)
        entry = tk.Entry(directoryframe, textvariable=self.directoryvariable, width=60)
        entry.grid(row=row, column=1)
        button = tk.Button(directoryframe, text='Select', command=self.selectdirectory)
        button.grid(row=row, column=2)

    def updatetitle(self, *args):
        try:
            self.title = self.titlevariable.get()
        except:
            return

    def updatenterms(self, *args):
        try:
            self.nterms = self.ntermsvariable.get()
        except:
            return

    def updateminfreq(self, *args):
        try:
            self.minfreq = self.minfreqvariable.get()
        except:
            return

    def updatemaxfreq(self, *args):
        try:
            self.maxfreq = self.maxfreqvariable.get()
        except:
            return

    def updatediag(self, *args):
        try:
            self.diag = self.diagvariable.get()
        except:
            return

    def updateradius(self, *args):
        try:
            self.radius = self.radiusvariable.get()
        except:
            return

    def updatedepth(self, *args):
        try:
            self.depth = self.depthvariable.get()
        except:
            return

    def updatesphres(self, *args):
        try:
            self.sphres = self.sphresvariable.get()
        except:
            return

    def updatenlyr(self, *args):
        try:
            self.nlyr = self.nlyrvariable.get()
            self.displaymainframe()
        except:
            return

    def updatenx(self, *args):
        try:
            self.nx = self.nxvariable.get()
        except:
            return

    def updateny(self, *args):
        try:
            self.ny = self.nyvariable.get()
        except:
            return

    def updatend(self, *args):
        try:
            self.nd = self.ndvariable.get()
        except:
            return

    def updatexmin(self, *args):
        try:
            self.xmin = self.xminvariable.get()
        except:
            return

    def updatexmax(self, *args):
        try:
            self.xmax = self.xmaxvariable.get()
        except:
            return

    def updateymin(self, *args):
        try:
            self.ymin = self.yminvariable.get()
        except:
            return

    def updateymax(self, *args):
        try:
            self.ymax = self.ymaxvariable.get()
        except:
            return

    def updatedmin(self, *args):
        try:
            self.dmin = self.dminvariable.get()
        except:
            return

    def updatedmax(self, *args):
        try:
            self.dmax = self.dmaxvariable.get()
        except:
            return

    def updatenf(self, *args):
        try:
            self.nf = self.nfvariable.get()
        except:
            return

    def setlogarithmic(self):
        self.logarithmic = self.logarithmicvariable.get()

    def displayparameterframe(self,parameterframe):
        row = 0
        label = tk.Label(parameterframe, text='Title')
        label.grid(row=row, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.titlevariable, width=80)
        entry.grid(row=row, column=1, columnspan=4)
        self.titlevariable.trace('w',self.updatetitle)
        row += 1
        label = tk.Label(parameterframe, text='Number of spherical degrees')
        label.grid(row=row, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.ntermsvariable)
        entry.grid(row=row, column=1)
        self.ntermsvariable.trace('w',self.updatenterms)
        label = tk.Label(parameterframe, text='Diagnostic')
        label.grid(row=row, column=2, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.diagvariable)
        entry.grid(row=row, column=3)
        self.diagvariable.trace('w',self.updatediag)
        row += 1
        label = tk.Label(parameterframe, text='Number of frequencies')
        label.grid(row=row, column=0, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.nfvariable)
        entry.grid(row=row, column=1)
        label = tk.Label(parameterframe, text='Minimum frequency')
        label.grid(row=row, column=2, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.minfreqvariable)
        entry.grid(row=row, column=3)
        label = tk.Label(parameterframe, text='Maximum frequency')
        label.grid(row=row, column=4, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.maxfreqvariable)
        entry.grid(row=row, column=5)
        self.nfvariable.trace('w',self.updatenf)
        self.minfreqvariable.trace('w',self.updateminfreq)
        self.maxfreqvariable.trace('w',self.updatemaxfreq)
        row += 1
        check = tk.Checkbutton(parameterframe, text='Logarithmic', variable=self.logarithmicvariable, command=self.setlogarithmic)
        self.logarithmicvariable.set(self.logarithmic)
        check.grid(row=row, column=0, sticky=tk.W)
        row += 1
        label = tk.Label(parameterframe, text='Radius')
        label.grid(row=row, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.radiusvariable)
        entry.grid(row=row, column=1)
        label = tk.Label(parameterframe, text='Depth')
        label.grid(row=row, column=2, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.depthvariable)
        entry.grid(row=row, column=3)
        label = tk.Label(parameterframe, text='Sphere resistivity')
        label.grid(row=row, column=4, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.sphresvariable)
        entry.grid(row=row, column=5)
        self.radiusvariable.trace('w',self.updateradius)
        self.depthvariable.trace('w',self.updatedepth)
        self.sphresvariable.trace('w',self.updatesphres)
        row += 1
        label = tk.Label(parameterframe, text='Number of layers')
        label.grid(row=row, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.nlyrvariable)
        entry.grid(row=row, column=1)
        self.nlyrvariable.trace('w',self.updatenlyr)
        row += 1
        for i in range(self.nlyr):
            self.layerlist[i].displayline(parameterframe,row,i)
            row += 1
        label = tk.Label(parameterframe, text='X number')
        label.grid(row=row, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.nxvariable)
        entry.grid(row=row, column=1)
        label = tk.Label(parameterframe, text='X minimum')
        label.grid(row=row, column=2, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.xminvariable)
        entry.grid(row=row, column=3)
        label = tk.Label(parameterframe, text='X maximum')
        label.grid(row=row, column=4, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.xmaxvariable)
        entry.grid(row=row, column=5)
        self.nxvariable.trace('w',self.updatenx)
        self.xminvariable.trace('w',self.updatexmin)
        self.xmaxvariable.trace('w',self.updatexmax)
        row += 1
        label = tk.Label(parameterframe, text='Y number')
        label.grid(row=row, column=0, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.nyvariable)
        entry.grid(row=row, column=1)
        label = tk.Label(parameterframe, text='Y minimum')
        label.grid(row=row, column=2, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.yminvariable)
        entry.grid(row=row, column=3)
        label = tk.Label(parameterframe, text='Y maximum')
        label.grid(row=row, column=4, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.ymaxvariable)
        entry.grid(row=row, column=5)
        self.nyvariable.trace('w', self.updateny)
        self.yminvariable.trace('w',self.updateymin)
        self.ymaxvariable.trace('w',self.updateymax)
        row += 1
        label = tk.Label(parameterframe, text='Depth number')
        label.grid(row=row, column=0, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.ndvariable)
        entry.grid(row=row, column=1)
        label = tk.Label(parameterframe, text='D minimum')
        label.grid(row=row, column=2, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.dminvariable)
        entry.grid(row=row, column=3)
        label = tk.Label(parameterframe, text='D maximum')
        label.grid(row=row, column=4, sticky=tk.W)
        entry = tk.Entry(parameterframe, textvariable=self.dmaxvariable)
        entry.grid(row=row, column=5)
        self.ndvariable.trace('w',self.updatend)
        self.dminvariable.trace('w',self.updatedmin)
        self.dmaxvariable.trace('w',self.updatedmax)
        self.refreshparameters()

    def displaycommandframe(self,commandframe):
        row = 0
        column = 0
        button = tk.Button(commandframe, text='Read Model', command=self.readcflfile)
        button.grid(row=0, column=column, sticky=tk.W)
        column += 1
        button = tk.Button(commandframe, text='Save Model', command=self.savecflfile)
        button.grid(row=0, column=column, sticky=tk.W)
        column += 1
        button = tk.Button(commandframe, text='Run Model', command=self.runmodel_with_wait)
        button.grid(row=0, column=column, sticky=tk.W)
        column += 1
        button = tk.Button(commandframe, text='Plot Sounding', command=self.mtplot.plotsounding)
        button.grid(row=0, column=column, sticky=tk.W)
        column += 1
        button = tk.Button(commandframe, text='Plot MT Profiles', command=self.mtplot.plotmtprofiles)
        button.grid(row=0, column=column, sticky=tk.W)
        column += 1
        button = tk.Button(commandframe, text='Plot Field Profiles', command=self.mtplot.plotfieldprofiles)
        button.grid(row=0, column=column, sticky=tk.W)
        column = column + 1
        button = tk.Button(commandframe, text='Plot Field Plane', command=self.mtplot.plotfieldplane)
        button.grid(row=0, column=column, sticky=tk.W)
        column = column + 1
        button = tk.Button(commandframe, text='Quit', command=self.quit)
        button.grid(row=0, column=column, sticky=tk.E)

    def selectdirectory(self):
        self.directory = filedialog.askdirectory(title='Working Directory')
        self.directoryvariable.set(self.directory)
        self.readcflfile()

    def savecflfile(self):
        if len(self.directory) == 0:
            messagebox.showerror('MTSphere','Directory undefined')
            return
        file = open(os.path.join(self.directory,'MTSphere.cfl'),'w')
        file.write('{}\n'.format(self.title))
        file.write(' {} {} {} {} {}! NLYR, NTERMS, NX, NY, ND, DIAG\n'.format(self.nlyr,self.nterms,self.nx,self.ny,self.nd,self.diag))
        file.write(' {} {} {} {} ! NF, MINFREQ, MAXFREQ, LOG\n'.format(self.nf,self.minfreq,self.maxfreq,self.logarithmic))
        file.write(' {} {} {} ! RADIUS, DEPTH, SPHERESISTIVITY\n'.format(self.radius,self.depth,self.sphres))
        for i in range(self.nlyr):
            file.write(' {}'.format(self.layerlist[i].resistivity))
            if i < self.nlyr-1:
                file.write(' {} ! RES, THK\n'.format(self.layerlist[i].thickness))
            else:
                file.write(' ! RES\n')
        file.write(' {} {} ! XMIN, XMAX\n'.format(self.xmin,self.xmax))
        file.write(' {} {} ! YMIN, YMAX\n'.format(self.ymin,self.ymax))
        file.write(' {} {} ! DMIN, DMAX\n'.format(self.dmin,self.dmax))
        file.close()

    def readcflfile(self):
        if len(self.directory) == 0:
            messagebox.showerror('MTSphere','Directory not set')
            return False
        filename = os.path.join(self.directory, 'MTSphere.cfl')
        if not os.path.exists(filename):
            messagebox.showwarning('MTSphere','No model file')
            return False
        cflfile = open(filename, 'r')
        self.title = cflfile.readline().rstrip()
        columns = cflfile.readline().split()
        self.nlyr = int(columns[0])
        self.layerlist = []
        for i in range(self.nlyr):
            self.layerlist.append(Layer(self.nlyr))
        self.nterms = int(columns[1])
        self.nx = int(columns[2])
        self.ny = int(columns[3])
        self.nd = int(columns[4])
        self.diag = int(columns[5])
        if self.nx < 0 or self.ny < 0:
            self.position = True
            self.nx = abs(self.nx)
            self.ny = abs(self.ny)
            self.nd = abs(self.nd)
        else:
            self.position = False
        columns = cflfile.readline().split()
        self.nf = int(columns[0])
        self.minfreq = float(columns[1])
        self.maxfreq = float(columns[2])
        self.logarithmic = int(columns[3])
        self.freq = np.zeros(self.nf)
        if self.nf == 1:
            self.freq[0] = self.minfreq
        else:
            if self.logarithmic == 0:
                for i in range(self.nf):
                    self.freq[i] = self.minfreq + i * ( self.maxfreq - self.minfreq ) / ( self.nf - 1 )
            elif self.logarithmic == 1:
                for i in range(self.nf):
                    self.freq[i] = self.minfreq * ( self.maxfreq / self.minfreq ) ** ( float(i) / float(self.nf-1) )
            else:
                columns = cflfile.readline().split()
                self.freq = [ float(x) for x in columns ]
        columns = cflfile.readline().split()
        self.radius = float(columns[0])
        self.depth = float(columns[1])
        self.sphres = float(columns[2])
        for i in range(self.nlyr):
            columns = cflfile.readline().split()
            self.layerlist[i].resistivity = float(columns[0])
            if i < self.nlyr - 1:
                self.layerlist[i].thickness = float(columns[1])
        x = np.zeros(self.nx)
        y = np.zeros(self.ny)
        d = np.zeros(self.nd)
        if self.position:
            columns = cflfile.readline().split()
            for i in range(self.nx):
                x[i] = float(columns[i])
            columns = cflfile.readline().split()
            for i in range(self.ny):
                y[i] = float(columns[i])
            columns = cflfile.readline().split()
            for i in range(self.nd):
                d[i] = float(columns[i])
        else:
            columns = cflfile.readline().split()
            self.xmin = float(columns[0])
            self.xmax = float(columns[1])
            columns = cflfile.readline().split()
            self.ymin = float(columns[0])
            self.ymax = float(columns[1])
            columns = cflfile.readline().split()
            self.dmin = float(columns[0])
            self.dmax = float(columns[1])
            if self.nx > 1:
                dx = ( self.xmax - self.xmin ) / ( self.nx - 1 )
            else:
                dx = 0
            if self.ny > 1:
                dy = ( self.ymax - self.ymin ) / ( self.ny - 1 )
            else:
                dy = 0
            if self.nd > 1:
                dd = ( self.dmax - self.dmin ) / ( self.nd - 1 )
            else:
                dd = 0
            for jx in range(self.nx):
                x[jx] = self.xmin + jx * dx
            for jy in range(self.ny):
                y[jy] = self.ymin + jy * dy
            for jd in range(self.nd):
                d[jd] = self.dmin + jd * dd
        self.pxyd = np.zeros((self.nx, self.ny, self.nd, 3))
        for jx in range(self.nx):
            for jy in range(self.ny):
                for jd in range(self.nd):
                    self.pxyd[jx,jy,jd, 0] = x[jx]
                    self.pxyd[jx,jy,jd, 1] = y[jy]
                    self.pxyd[jx,jy,jd, 2] = d[jd]
        self.xmin = x.min()
        self.xmax = x.max()
        self.ymin = y.min()
        self.ymax = y.max()
        self.dmin = d.min()
        self.dmax = d.max()
        cflfile.close()
        self.refreshparameters()
        return True

    def refreshparameters(self):
        self.titlevariable.set(self.title)
        self.nlyrvariable.set(self.nlyr)
        self.titlevariable.set(self.title)
        self.ntermsvariable.set(self.nterms)
        self.nxvariable.set(self.nx)
        self.nyvariable.set(self.ny)
        self.ndvariable.set(self.nd)
        self.diagvariable.set(self.diag)
        self.nfvariable.set(self.nf)
        self.minfreqvariable.set(self.minfreq)
        self.maxfreqvariable.set(self.maxfreq)
        self.logarithmicvariable.set(self.logarithmic)
        self.radiusvariable.set(self.radius)
        self.depthvariable.set(self.depth)
        self.sphresvariable.set(self.sphres)
        for layer in self.layerlist:
            layer.refresh()
        self.xminvariable.set(self.xmin)
        self.xmaxvariable.set(self.xmax)
        self.yminvariable.set(self.ymin)
        self.ymaxvariable.set(self.ymax)
        self.dminvariable.set(self.dmin)
        self.dmaxvariable.set(self.dmax)

    def setvectors(self):
        self.freq = np.zeros(self.nf)
        self.zhat = np.zeros((self.nf, self.nlyr))
        self.imp = np.zeros((self.nf,self.nx,self.ny,3,2))
        self.res = np.zeros(self.nlyr+1)
        self.thk = np.zeros(self.nlyr)
        if self.nf == 1:
            self.freq[0] = self.minfreq
        else:
            if self.logarithmic == 0:
                for i in range(self.nf):
                    self.freq[i] = self.minfreq + i * ( self.maxfreq - self.minfreq ) / ( self.nf - 1 )
            else:
                for i in range(self.nf):
                    self.freq[i] = self.minfreq * ( self.maxfreq / self.minfreq ) ** ( float(i) / float(self.nf-1) )
        for i in range(self.nlyr):
            self.res[i+1] = self.layerlist[i].resistivity
            if i < self.nlyr-1:
                self.thk[i+1] = self.layerlist[i].thickness

    def write0(self,message,level):
        self.outfile.write(message)
        if self.diag > level:
            sys.stdout.write(message)

    def write1(self,message,level):
        self.outfile.write(message)
        self.mf1file.write(message)
        if self.diag > level:
            sys.stdout.write(message)

    def write2(self,message,level):
        self.outfile.write(message)
        self.mf2file.write(message)
        if self.diag > level:
            sys.stdout.write(message)

    def openoutputfiles(self):
        if len(self.directory) == 0:
            messagebox.showerror('MTSphere','Directory not set')
            return
        filename = os.path.join(self.directory,'MTSphere.cfl')
        if not os.path.exists(filename):
            messagebox.showerror('MTSphere','No model file')
            return
        cflfile = open(filename,'r')
        self.outfile = open(os.path.join(self.directory,'MTSphere.out'),'w')
        self.mf1file = open(os.path.join(self.directory, 'MTSphere.mf1'), 'w')
        self.mf2file = open(os.path.join(self.directory,'MTSphere.mf2'), 'w')
        self.outfile.write('                             MTSphere - Version {}   {}\n'.\
                           format(self.version,self.releasedate))
        self.outfile.write('                             {}\n'.format(developer))
        self.outfile.write('                             {}\n\n\n\n'.format(company))
        print('MTSphere - Version ',self.version,self.releasedate)
        self.outfile.write('           INPUT DATA\n           ___________\n\n')
        while 1:
            line = cflfile.readline()
            if not line:
                break
            self.outfile.write('{}'.format(line))
        cflfile.close()
        self.readcflfile()
        self.setvectors()
        self.outfile.write('-------------------------------------------------------------------------------\n')
        self.outfile.write('\n{}\n'.format(self.title))
        self.outfile.write(' Number of layers: {}\n'.format(self.nlyr))
        self.outfile.write(' Number of frequencies: {}\n'.format(self.nf))
        for i in range(1,self.nlyr):
            self.outfile.write(' Layer {} : Resistivity:{}, Thickness: {}\n'.format(i,self.res[i],self.thk[i]))
        self.outfile.write(' Layer {} : Resistivity:{}\n'.format(self.nlyr,self.res[self.nlyr]))
        self.outfile.write(' Frequencies:')
        for i in range(self.nf):
            self.outfile.write(' {:15.7g}'.format(self.freq[i]))
        self.outfile.write('\n')

    def writeresults(self):

        for i in range(2):
            for j in range(2):
                if i == 0:
                    if j == 0:
                        self.write1('\n XX Apparent resistivity, phase and impedance\n',2)
                    else:
                        self.write1('\n XY Apparent resistivity, phase and impedance\n',2)
                else:
                    if j == 0:
                        self.write1('\n YX Apparent resistivity, phase and impedance\n',2)
                    else:
                        self.write1('\n YY Apparent resistivity, phase and impedance\n',2)
                self.write1('{:>15s}{:>10s}{:>10s}{:>10s}{:>15s}{:>15s}{:>15s}{:>15s}\n'. \
                         format('Frequency','X','Y','D','Apparent_res.','Phase',\
                                'Impedance(real)','Imp.(imag)'),2)
                for jf in range(self.nf):
                    for jd in range(self.nd):
                        for jx in range(self.nx):
                            for jy in range(self.ny):
                                self.write1('{:15.7g}{:10.1f}{:10.1f}{:10.1f}{:15.7g}{:15.7g}{:15.7g}{:15.7g}\n'.\
                                         format(self.freq[jf],\
                                        self.pxyd[jx,jy,jd,0],self.pxyd[jx,jy,jd,1],self.pxyd[jx,jy,jd,2], \
                                        self.appres[jf,jx,jy,jd,i,j], self.phase[jf,jx,jy,jd,i,j],\
                                        self.imp[jf,jx,jy,jd,i,j].real, self.imp[jf,jx,jy,jd,i,j].imag),2)

        self.write1('\n Tipper\n',2)
        self.write1('{:>15s}{:>10s}{:>10s}{:>10s}{:>15s}{:>15s}{:>15s}{:>15s}\n'. \
         format('Frequency', 'X', 'Y', 'D', 'Kzx(real)', 'Kzx(imag)', 'Kzy(real)', 'Kzy(imag'),2)
        for jf in range(self.nf):
            for jd in range(self.nd):
                for jx in range(self.nx):
                    for jy in range(self.ny):
                        self.write1('{:15.7g}{:10.1f}{:10.1f}{:10.1f}{:15.7g}{:15.7g}{:15.7g}{:15.7}\n'.\
                             format(self.freq[jf],self.pxyd[jx,jy,jd,0],self.pxyd[jx,jy,jd,1],self.pxyd[jx,jy,jd,2],
                                    self.imp[jf,jx,jy,jd,2,0].real,self.imp[jf,jx,jy,jd,2,0].imag,
                                    self.imp[jf,jx,jy,jd,2,1].real,self.imp[jf,jx,jy,jd,2,1].imag),2)

        for ic in range(2):
            if ic == 0:
                self.write2( '\n X polarization electric and magnetic fields\n',2)
            else:
                self.write2( '\n Y polarization electric and magnetic fields\n',2)
            self.write2(('{:>15s}{:>10s}{:>10s}{:>10s}'+
                      '{:>15s}{:>15s}{:>15s}{:>15s}'+
                      '{:>15s}{:>15s}{:>15s}{:>15s}{:>15s}{:>15s}\n'). \
                    format('Frequency','X','Y','D','Ex(real)','Ex(imag)','Ey(real)','Ey(imag)', 'Ez(real)', 'Ez(imag)', \
                    'Hx(real)','Hx(imag)','Hy(real)','Hy(imag)','Hz(real)','Hz(imag)'),2)
            self.write2('\nSecondary\n',2)
            for jf in range(self.nf):
                for jd in range(self.nd):
                    for jx in range(self.nx):
                        for jy in range(self.ny):
                            self.write2(("{:15.7g}{:10.1f}{:10.1f}{:10.1f}"+
                                      "{:15.7g}{:15.7g}{:15.7g}{:15.7}"+
                                      "{:15.7g}{:15.7g}{:15.7g}{:15.7}{:15.7g}{:15.7}{:15.7g}{:15.7}\n"). \
                                     format(self.freq[jf], self.pxyd[jx, jy, jd, 0], self.pxyd[jx, jy, jd, 1],
                                            self.pxyd[jx, jy, jd, 2],
                                     self.Es[jf, jx, jy, jd, 0, ic].real, self.Es[jf, jx, jy, jd, 0, ic].imag,
                                     self.Es[jf, jx, jy, jd, 1, ic].real, self.Es[jf, jx, jy, jd, 1, ic].imag,
                                     self.Es[jf, jx, jy, jd, 2, ic].real, self.Es[jf, jx, jy, jd, 2, ic].imag,
                                     self.Hs[jf, jx, jy, jd, 0, ic].real, self.Hs[jf, jx, jy, jd, 0, ic].imag,
                                     self.Hs[jf, jx, jy, jd, 1, ic].real, self.Hs[jf, jx, jy, jd, 1, ic].imag,
                                     self.Hs[jf, jx, jy, jd, 2, ic].real, self.Hs[jf, jx, jy, jd, 2, ic].imag),2)
            self.write2('\nTotal\n',2)
            for jf in range(self.nf):
                for jd in range(self.nd):
                    for jx in range(self.nx):
                        for jy in range(self.ny):
                            self.write2(("{:15.7g}{:10.1f}{:10.1f}{:10.1f}" +
                                      "{:15.7g}{:15.7g}{:15.7g}{:15.7}" +
                                      "{:15.7g}{:15.7g}{:15.7g}{:15.7}{:15.7g}{:15.7}{:15.7g}{:15.7}\n"). \
                                     format(self.freq[jf], self.pxyd[jx, jy, jd, 0], self.pxyd[jx, jy, jd, 1], self.pxyd[jx, jy, jd, 2],
                                            self.Et[jf, jx, jy, jd, 0, ic].real, self.Et[jf, jx, jy, jd, 0, ic].imag,
                                            self.Et[jf, jx, jy, jd, 1, ic].real, self.Et[jf, jx, jy, jd, 1, ic].imag,
                                            self.Et[jf, jx, jy, jd, 2, ic].real, self.Et[jf, jx, jy, jd, 2, ic].imag,
                                            self.Ht[jf, jx, jy, jd, 0, ic].real, self.Ht[jf, jx, jy, jd, 0, ic].imag,
                                            self.Ht[jf, jx, jy, jd, 1, ic].real, self.Ht[jf, jx, jy, jd, 1, ic].imag,
                                            self.Ht[jf, jx, jy, jd, 2, ic].real, self.Ht[jf, jx, jy, jd, 2, ic].imag),2)
        self.outfile.close()
        self.mf1file.close()
        self.mf2file.close()

    def runmodel_with_wait(self):
        if len(self.directory) == 0:
            messagebox.showerror('MTSphere','Directory not set')
            return
        title = 'pyMTSphere'
        wait_win = WaitWindow(self.root,title)
        def worker():
            try:
                self.runmodel()
            finally:
                elapsed = wait_win.elapsed
                wait_win.destroy()
                heures = elapsed // 3600
                minutes = (elapsed % 3600) // 60
                secondes = elapsed % 60
                messagebox.showinfo(title,
                f'Computation completed\nafter {heures:02d}:{minutes:02d}:{secondes:02d} elapsed time')
        thread = threading.Thread(target=worker, daemon=True)
        thread.start()

    def runmodel(self):
        self.openoutputfiles()
        if self.diag > 0:
            print(' Compute plane wave impedance')
        k, Z, E = PlaneWaveImpedance(self.write0, self.nlyr, self.thk, self.res, self.nf, self.freq)
        dpthl = np.zeros(self.nlyr + 1)
        for jz in range(1, self.nlyr):
            dpthl[jz + 1] = dpthl[jz] + self.thk[jz]
        if self.diag > 0:
            print(' Start computation of secondary fields')
        self.imp, self.Es, self.Hs = \
            mtsphere3d(self.write0, self.nf, self.nlyr, self.nterms, \
                   self.nx,  self.ny, self.nd, self.freq, self.pxyd, self.thk, self.res, self.depth, \
                   self.radius, self.sphres, E, dpthl, self.htarg)
        if self.diag > 0:
            print(' Compute magnetotelluric impedances')
        self.imp, self.Et, self.Ht = \
            getimpedance(self.nf, self.nx, self.ny, self.nd,
                         self.nlyr, self.nterms, self.pxyd, dpthl, k, Z, self.depth,
                         self.radius, E, self.Es, self.Hs )
        if self.diag > 0:
            print(' Compute apparent resistivity')
        self.appres, self.phase = apparentresistivity( self.nf, self.nx, self.ny, self.nd, self.freq, self.imp)
        if self.diag > 0:
            print(' Write results')
        self.writeresults()

    def retrievemodel(self):
        pass

    def quit(self):
        self.root.destroy()

if __name__ == '__main__':
    pymtsphere = pyMTSphere()
    pymtsphere.displaymainframe()
    pymtsphere.root.mainloop()
