
import os
import numpy as np
import tkinter as tk
from tkinter import filedialog,messagebox
import matplotlib
matplotlib.use('tkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

class MTPlot:

    def __init__(self,main):
        self.main = main

    def plot(self):
        if self.main.main == 0:
            self.main.setdirectory('Sounding')
            self.readcflfile()
            self.sounding()
        if self.main.main == 1:
            self.main.setdirectory('XProfile')
            self.readcflfile()
            self.plotxprofile()
        elif self.main.main == 2:
            self.main.setdirectory('YProfile')
            self.readcflfile()
            self.plotyprofile()

    def apparentresistivity(self):
        figure = Figure([8, 8])
        axes = []
        for i in range(4):
            axes.append(figure.add_subplot(2,2,i+1))

        self.readcflfile()
        mf1file = open(os.path.join(self.main.directory,'MTSphere.mf1'), 'r')
        px = np.zeros(self.nx)
        appres = np.zeros((self.nx,2,2))
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    mf1file.readline()
                for jx in range(self.nx):
                    columns = [ float(f) for f in mf1file.readline().split() ]
                    px[jx] = columns[1]
                    appres[jx,i,j] = columns[3]
        mf1file.close()

        fer12 = []
        file = open(os.path.join(self.main.directory,'xy_R12.txt'))
        while 1:
            line = file.readline()
            if not line:
                break
            columns = line.split()
            if float(columns[0]) != 25:
                fer12.append(float(columns[2]))
        file.close()
        fer12 = np.array(fer12)
        fer21 = []
        file = open(os.path.join(self.main.directory,'xy_R21.txt'))
        while 1:
            line = file.readline()
            if not line:
                break
            columns = line.split()
            if float(columns[0]) != 25:
                fer21.append(float(columns[2]))
        file.close()
        fer21 = np.array(fer21)

        axes[0].plot(px,fer12,'.b',label='FE')
        axes[0].plot(px,appres[:,0,1],'.r',label='Sphere')
        axes[0].set_xlabel('X distance (m)')
        axes[0].set_ylabel('Apparent resistivity ($\Omega$.m)')
        axes[0].set_title('XY')
        axes[0].legend(loc='center right')
        axes[1].plot(px,fer21,'.b',label='FE')
        axes[1].plot(px,appres[:,1,0],'.r',label='Sphere')
        axes[1].set_xlabel('X distance (m)')
        axes[1].set_title('YX')
        axes[1].legend(loc='center right')

        self.readcflfile()
        mf1file = open(os.path.join(self.main.directory,'MTSphere.mf1'), 'r')
        py = np.zeros(self.ny)
        appres = np.zeros((self.ny,2,2))
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    mf1file.readline()
                for jy in range(self.ny):
                    columns = [ float(f) for f in mf1file.readline().split() ]
                    py[jy] = columns[2]
                    appres[jy,i,j] = columns[3]
        mf1file.close()

        fer12 = []
        file = open(os.path.join(self.main.directory,'xy_R12.txt'))
        while 1:
            line = file.readline()
            if not line:
                break
            columns = line.split()
            if float(columns[0]) == 25:
                fer12.append(float(columns[2]))
        file.close()
        fer12 = np.array(fer12)
        fer21 = []
        file = open(os.path.join(self.main.directory,'xy_R21.txt'))
        while 1:
            line = file.readline()
            if not line:
                break
            columns = line.split()
            if float(columns[0]) == 25:
                fer21.append(float(columns[2]))
        file.close()
        fer21 = np.array(fer21)

        axes[2].plot(py,fer12,'.b',label='FE')
        axes[2].plot(py,appres[:,0,1],'.r',label='Sphere')
        axes[2].set_xlabel('Y distance (m)')
        axes[2].set_ylabel('Apparent resistivity ($\Omega$.m)')     #   axes[2].set_title('XY')
        axes[2].legend(loc='center right')
        axes[3].plot(py,fer21,'.b',label='FE')
        axes[3].plot(py,appres[:,1,0],'.r',label='Sphere')
        axes[3].set_xlabel('Y distance (m)')
        axes[3].legend(loc='center right')
        for i in range(4):
            if self.main.main == 0:
                axes[i].set_ylim(60,110)
            else:
                axes[i].set_ylim(300,400)

        master = tk.Toplevel()
        canvasplot = FigureCanvasTkAgg(figure, master=master)
        canvasplot.draw()
        canvasplot.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        canvasplot._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        button = tk.Button(master, text='Quit', command=master.destroy)
        button.pack(side=tk.RIGHT)
        button = tk.Button(master=master, text='Save JPG', command=lambda: self.savefigure(figure))
        button.pack(side=tk.RIGHT)

    def difference(self):
        figure = Figure([8, 8])
        axes = []
        for i in range(4):
            axes.append(figure.add_subplot(2,2,i+1))

        self.main.readcflfile()
        mf1file = open(os.path.join(self.main.directory,'MTSphere.mf1'), 'r')
        px = np.zeros(self.nx)
        appres = np.zeros((self.nx,2,2))
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    mf1file.readline()
                for jx in range(self.main.nx):
                    columns = [ float(f) for f in mf1file.readline().split() ]
                    px[jx] = columns[1]
                    appres[jx,i,j] = columns[3]
        mf1file.close()

        fer12 = []
        file = open(os.path.join(self.main.directory,'xy_R12.txt'))
        while 1:
            line = file.readline()
            if not line:
                break
            columns = line.split()
            if float(columns[0]) != 25:
                fer12.append(float(columns[2]))
        file.close()
        fer12 = np.array(fer12)
        fer21 = []
        file = open(os.path.join(self.main.directory,'xy_R21.txt'))
        while 1:
            line = file.readline()
            if not line:
                break
            columns = line.split()
            if float(columns[0]) != 25:
                fer21.append(float(columns[2]))
        file.close()
        fer21 = np.array(fer21)

        axes[0].plot(px,100*(np.divide(fer12,appres[:,0,1])-1),'.k')
        axes[0].set_xlabel('X distance (m)')
        axes[0].set_ylabel('Percent difference (%)')
        axes[0].set_title('XY')
        axes[0].legend(loc='center right')
        axes[1].plot(px,100*(np.divide(fer21,appres[:,1,0])-1),'.k')
        axes[1].set_xlabel('X distance (m)')
        axes[1].set_title('YX')
        axes[1].legend(loc='center right')

        self.main.readcflfile()
        mf1file = open(os.path.join(self.main.directory,'MTSphere.mf1'), 'r')
        py = np.zeros(self.ny)
        appres = np.zeros((self.ny,2,2))
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    mf1file.readline()
                for jy in range(self.ny):
                    columns = [ float(f) for f in mf1file.readline().split() ]
                    py[jy] = columns[2]
                    appres[jy,i,j] = columns[3]
        mf1file.close()

        fer12 = []
        file = open(os.path.join(self.main.directory,'xy_R12.txt'))
        while 1:
            line = file.readline()
            if not line:
                break
            columns = line.split()
            if float(columns[0]) == 25:
                fer12.append(float(columns[2]))
        file.close()
        fer12 = np.array(fer12)
        fer21 = []
        file = open(os.path.join(self.main.directory,'xy_R21.txt'))
        while 1:
            line = file.readline()
            if not line:
                break
            columns = line.split()
            if float(columns[0]) == 25:
                fer21.append(float(columns[2]))
        file.close()
        fer21 = np.array(fer21)

        axes[2].plot(px,100*(np.divide(fer12,appres[:,0,1])-1),'.k')
        axes[2].set_xlabel('Y distance (m)')
        axes[2].set_ylabel('Percent difference (%)')
        axes[3].plot(px,100*(np.divide(fer21,appres[:,1,0])-1),'.k')
        axes[3].set_xlabel('Y distance (m)')
        for i in range(4):
            axes[i].set_ylim(-1,3)

        master = tk.Toplevel()
        canvasplot = FigureCanvasTkAgg(figure, master=master)
        canvasplot.draw()
        canvasplot.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        canvasplot._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        button = tk.Button(master, text='Quit', command=master.destroy)
        button.pack(side=tk.RIGHT)
        button = tk.Button(master=master, text='Save JPG', command=lambda: self.savefigure(figure))
        button.pack(side=tk.RIGHT)

    def orthodifference(self):
        figure = Figure([8, 8])
        axes = []
        for i in range(2):
            axes.append(figure.add_subplot(2,1,i+1))

        self.main.readcflfile()
        mf1file = open(os.path.join(self.main.directory,'MTSphere.mf1'), 'r')
        px = np.zeros(self.nx)
        xappres = np.zeros((self.nx,2,2))
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    mf1file.readline()
                for jx in range(self.nx):
                    columns = [ float(f) for f in mf1file.readline().split() ]
                    px[jx] = columns[1]
                    xappres[jx,i,j] = columns[3]
        mf1file.close()

        self.main.readcflfile()
        mf1file = open(os.path.join(self.main.directory,'MTSphere.mf1'), 'r')
        py = np.zeros(self.ny)
        yappres = np.zeros((self.ny,2,2))
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    mf1file.readline()
                for jy in range(self.ny):
                    columns = [ float(f) for f in mf1file.readline().split() ]
                    py[jy] = columns[2]
                    yappres[jy,i,j] = columns[3]
        mf1file.close()

        axes[0].plot(px,100*(np.divide(yappres[:,1,0],xappres[:,0,1])-1),'.k')
        axes[0].set_xlabel('Distance (m)')
        axes[0].set_ylabel('Percent difference (%)')
        axes[0].set_title('XY')

        axes[1].plot(px,100*(np.divide(yappres[:,0,1],xappres[:,1,0])-1),'.k')
        axes[1].set_xlabel('Distance (m)')
        axes[1].set_ylabel('Percent difference (%)')
        axes[1].set_title('YX')

        master = tk.Toplevel()
        canvasplot = FigureCanvasTkAgg(figure, master=master)
        canvasplot.draw()
        canvasplot.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        canvasplot._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        button = tk.Button(master, text='Quit', command=master.destroy)
        button.pack(side=tk.RIGHT)
        button = tk.Button(master=master, text='Save JPG', command=lambda: self.savefigure(figure))
        button.pack(side=tk.RIGHT)

    def tipper(self):

        self.readcflfile()
        mf1file = open(os.path.join(self.test.directory, 'MTSphere.mf1'), 'r')
        px = np.zeros(self.nx)
        Kzx = np.zeros(self.nx, dtype=complex)
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    mf1file.readline()
                for jx in range(self.nx):
                    columns = [float(f) for f in mf1file.readline().split()]
                    px[jx] = columns[1]
        for k in range(3):
            mf1file.readline()
        for jx in range(self.nx):
            columns = [float(f) for f in mf1file.readline().split()]
            Kzx[jx] = columns[3] + 1j * columns[4]
        mf1file.close()

        Tzx = []
        file = open(os.path.join(self.test.directory,'xy_Tzx.txt'))
        while 1:
            ligne = file.readline().rstrip()
            if not ligne:
                break
            import re

            # Test du match
            pattern = r'([-\d.E+-]+)\s+([-\d.E+-]+)\s*\(\s*([-\d.E+-]+)\s*,\s*([-\d.E+-]+)\s*\)'
            match = re.match(pattern, ligne.strip())
            if match:
                a, b, c, d = map(float, match.groups())
                if a != 25:
                    Tzx.append(-c-1j*d)
            else:
                print("Échec du match. Pattern utilisé :", pattern)
        file.close()
        Tzx = np.array(Tzx, dtype=complex)

        figure = Figure([8, 8])
        axe = figure.add_subplot(1, 1, 1)
        axe.plot(px, 100 * Kzx.real, 'b', label='Sphere Real')
        axe.plot(px, 100 * Kzx.imag, 'r', label='Sphere Imaginary')
        axe.plot(px, 100 * Tzx.real, '.b', label='FE Real')
        axe.plot(px, 100 * Tzx.imag, '.r', label='FE Imaginary')
        axe.set_xlabel('X Distance (m)')
        axe.set_ylabel('X Tipper ratio (%)')
        axe.legend(loc='upper right')

        master = tk.Toplevel()
        canvasplot = FigureCanvasTkAgg(figure, master=master)
        canvasplot.draw()
        canvasplot.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        canvasplot._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        button = tk.Button(master, text='Quit', command=master.destroy)
        button.pack(side=tk.RIGHT)
        button = tk.Button(master=master, text='Save JPG', command=lambda: self.savefigure(figure))
        button.pack(side=tk.RIGHT)

    def plotsounding(self):
        self.main.readcflfile()
        if self.main.nx > 1:
            messagebox.showerror('MTSphere', 'More than one X point')
            return
        self.main.readcflfile()
        if self.main.ny > 1:
            messagebox.showerror('MTSphere', 'More than one Y point')
            return
        mf1file = open(os.path.join(self.main.directory, 'MTSphere.mf1'), 'r')
        freq = np.zeros(self.main.nf)
        appres = np.zeros((self.main.nf,2,2))
        phase = np.zeros((self.main.nf,2,2))
        imp = np.zeros((self.main.nf,2,2), dtype=complex)
        Kzx = np.zeros(self.main.nf, dtype=complex)
        Kzy = np.zeros(self.main.nf, dtype=complex)
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    mf1file.readline()
                for jf in range(self.main.nf):
                    columns = [ float(f) for f in mf1file.readline().split() ]
                    freq[jf] = columns[0]
                    appres[jf,i,j] = columns[3]
                    phase[jf,i,j] = columns[4]
                    imp[jf,i,j] = columns[5] + 1j * columns[6]
        for k in range(3):
            mf1file.readline()
        for jf in range(self.main.nf):
            columns = [ float(f) for f in mf1file.readline().split() ]
            Kzx[jf] = columns[3] + 1j * columns[4]
            Kzy[jf] = columns[5] + 1j * columns[6]
        mf1file.close()
        figure = Figure([8, 8])
        axes = []
        for i in range(4):
            axes.append(figure.add_subplot(4,1,i+1))
    #        axe1s.append(figure1.add_subplot(212))
        axes[0].plot(freq,appres[:,0,1],'b',label='XY')
        axes[0].plot(freq,appres[:,1,0],'r',label='YX')
        axes[1].plot(freq,phase[:,0,1],'b',label='XY')
        axes[1].plot(freq,phase[:,1,0]+180,'r',label='YX+180')
        axes[1].set_xlabel('Frequency (Hz)')
        axes[0].set_ylabel('Apparent resistivity ($\Omega$.m)')
        axes[1].set_ylabel('Phase ($^\circ$)')
        axes[0].legend(loc='upper right')
        axes[1].legend(loc='upper right')
        axes[0].set_xscale('log')
        axes[0].set_xlim(min(freq),max(freq))
        axes[1].set_xscale('log')
        axes[1].set_xlim(min(freq),max(freq))

#        figure2 = Figure([8, 6])
#        axe = figure2.add_subplot(111)
        axes[2].plot(freq,100*Kzx.real,'b',label='Real')
        axes[2].plot(freq,100*Kzx.imag,'r',label='Imaginary')
        axes[2].set_xlabel('Frequency (Hz)')
        axes[2].set_ylabel('X Tipper ratio (%)')
        axes[2].legend(loc='upper right')
        axes[2].set_xscale('log')
        axes[2].set_xlim(min(freq),max(freq))

#        figure3 = Figure([8, 6])
#        axe = figure3.add_subplot(111)
        axes[3].plot(freq,100*Kzy.real,'b',label='Real')
        axes[3].plot(freq,100*Kzy.imag,'r',label='Imaginary')
#        figure3.suptitle('Tipper Sounding')
        axes[3].set_xlabel('Frequency (Hz)')
        axes[3].set_ylabel('Y Tipper ratio (%)')
        axes[3].legend(loc='upper right')
        axes[3].set_ylim(-1,1)
        axes[3].set_xscale('log')
        axes[3].set_xlim(min(freq),max(freq))

        master = tk.Toplevel()
        canvasplot = FigureCanvasTkAgg(figure, master=master)
        canvasplot.draw()
        canvasplot.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        canvasplot._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        button = tk.Button(master, text='Quit', command=master.destroy)
        button.pack(side=tk.RIGHT)
        button = tk.Button(master=master, text='Save JPG', command=lambda: self.savefigure(figure))
        button.pack(side=tk.RIGHT)

#        master2 = tk.Toplevel()
#        canvasplot = FigureCanvasTkAgg(figure2, master=master2)
#        canvasplot.draw()
#        canvasplot.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
#        canvasplot._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
#        button = tk.Button(master2, text='Quit', command=master2.destroy)
#        button.pack(side=tk.RIGHT)
#        button = tk.Button(master=master2, text='Save JPG', command=lambda: self.savefigure(figure2))
#        button.pack(side=tk.RIGHT)

#        master3 = tk.Toplevel()
#        canvasplot = FigureCanvasTkAgg(figure3, master=master3)
#        canvasplot.draw()
#        canvasplot.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
#        canvasplot._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
#        button = tk.Button(master3, text='Quit', command=master3.destroy)
#        button.pack(side=tk.RIGHT)
#        button = tk.Button(master=master3, text='Save JPG', command=lambda: self.savefigure(figure3))
#        button.pack(side=tk.RIGHT)

    def plotxprofile(self):
        self.main.readcflfile()
        if self.main.ny > 1:
            messagebox.showerror('MTSphere', 'More than one Y point')
            return
        self.main.readcflfile()
        if self.main.nf > 1:
            messagebox.showerror('MTSphere', 'More than one frequency')
            return
        mf1file = open(os.path.join(self.main.directory, 'MTSphere.mf1'), 'r')
        px = np.zeros(self.main.nx)
        appres = np.zeros((self.main.nx,2,2))
        phase = np.zeros((self.main.nx,2,2))
        imp = np.zeros((self.main.nx,2,2), dtype=complex)
        Kzx = np.zeros(self.main.nx, dtype=complex)
        Kzy = np.zeros(self.main.nx, dtype=complex)
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    mf1file.readline()
                for jx in range(self.main.nx):
                    columns = [ float(f) for f in mf1file.readline().split() ]
                    px[jx] = columns[1]
                    appres[jx,i,j] = columns[3]
                    phase[jx,i,j] = columns[4]
                    imp[jx,i,j] = columns[5] + 1j * columns[6]
        for k in range(3):
            mf1file.readline()
        for jx in range(self.main.nx):
            columns = [ float(f) for f in mf1file.readline().split() ]
            Kzx[jx] = columns[3] + 1j * columns[4]
            Kzy[jx] = columns[5] + 1j * columns[6]
        mf1file.close()
        figure = Figure([8, 8])
        axes = []
        for i in range(4):
            axes.append(figure.add_subplot(4,1,i+1))
        axes[0].plot(px,appres[:,0,1],'b',label='XY')
        axes[0].plot(px,appres[:,1,0],'r',label='YX')
        axes[0].set_ylabel('Apparent resistivity ($\Omega$.m)')
        axes[0].legend(loc='upper right')
        axes[1].plot(px,phase[:,0,1],'b',label='XY')
        axes[1].plot(px,phase[:,1,0]+180,'r',label='YX+180')
        axes[1].set_xlabel('Distance (m)')
        axes[1].set_ylabel('Phase ($^\circ$)')
        axes[1].legend(loc='upper right')
        axes[2].plot(px,100*Kzx.real,'b',label='Real')
        axes[2].plot(px,100*Kzx.imag,'r',label='Imaginary')
        axes[2].set_xlabel('X Distance (m)')
        axes[2].set_ylabel('X Tipper ratio (%)')
        axes[2].legend(loc='upper right')
        axes[3].plot(px,100*Kzy.real,'b',label='Real')
        axes[3].plot(px,100*Kzy.imag,'r',label='Imaginary')
        axes[3].set_xlabel('X Distance (m)')
        axes[3].set_ylabel('Y Tipper ratio (%)')
        axes[3].legend(loc='upper right')
        axes[3].set_ylim(-1,1)

        master = tk.Toplevel()
        canvasplot = FigureCanvasTkAgg(figure, master=master)
        canvasplot.draw()
        canvasplot.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        canvasplot._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        button = tk.Button(master, text='Quit', command=master.destroy)
        button.pack(side=tk.RIGHT)
        button = tk.Button(master=master, text='Save JPG', command=lambda: self.savefigure(figure))
        button.pack(side=tk.RIGHT)

    def plotyprofile(self):
        self.main.readcflfile()
        if self.main.nx > 1:
            messagebox.showerror('MTSphere','More than one X point')
            return
        self.main.readcflfile()
        if self.main.nf > 1:
            messagebox.showerror('MTSphere','More than one frequency')
            return
        mf1file = open(os.path.join(self.main.directory, 'MTSphere.mf1'), 'r')
        py = np.zeros(self.main.ny)
        appres = np.zeros((self.main.ny,2,2))
        phase = np.zeros((self.main.ny,2,2))
        imp = np.zeros((self.main.ny,2,2), dtype=complex)
        Kzx = np.zeros(self.main.ny, dtype=complex)
        Kzy = np.zeros(self.main.ny, dtype=complex)
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    mf1file.readline()
                for jy in range(self.main.ny):
                    columns = [ float(f) for f in mf1file.readline().split() ]
                    py[jy] = columns[2]
                    appres[jy,i,j] = columns[3]
                    phase[jy,i,j] = columns[4]
                    imp[jy,i,j] = columns[5] + 1j * columns[6]
        for k in range(3):
            mf1file.readline()
        for jy in range(self.main.ny):
            columns = [ float(f) for f in mf1file.readline().split() ]
            Kzx[jy] = columns[3] + 1j * columns[4]
            Kzy[jy] = columns[5] + 1j * columns[6]
        mf1file.close()
        figure = Figure([8, 8])
        axes = []
        for i in range(4):
            axes.append(figure.add_subplot(4,1,i+1))
        axes[0].plot(py,appres[:,0,1],'b',label='XY')
        axes[0].plot(py,appres[:,1,0],'r',label='YX')
        axes[0].set_ylabel('Apparent resistivity ($\Omega$.m)')
        axes[0].legend(loc='upper right')
        axes[1].plot(py,phase[:,0,1],'b',label='XY')
        axes[1].plot(py,phase[:,1,0]+180,'r',label='YX+180')
        axes[1].set_xlabel('Y Distance (m)')
        axes[1].set_ylabel('Phase ($^\circ$)')
        axes[1].legend(loc='upper right')
        axes[2].plot(py,100*Kzx.real,'b',label='Real')
        axes[2].plot(py,100*Kzx.imag,'r',label='Imaginary')
        axes[2].set_xlabel('Y Distance (m)')
        axes[2].set_ylabel('X Tipper ratio (%)')
        axes[2].legend(loc='upper right')
        axes[3].plot(py,100*Kzy.real,'b',label='Real')
        axes[3].plot(py,100*Kzy.imag,'r',label='Imaginary')
        axes[3].set_xlabel('Y Distance (m)')
        axes[3].set_ylabel('Y Tipper ratio (%)')
        axes[3].legend(loc='upper right')

        master = tk.Toplevel()
        canvasplot = FigureCanvasTkAgg(figure, master=master)
        canvasplot.draw()
        canvasplot.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        canvasplot._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        button = tk.Button(master, text='Quit', command=master.destroy)
        button.pack(side=tk.RIGHT)
        button = tk.Button(master=master, text='Save JPG', command=lambda: self.savefigure(figure))
        button.pack(side=tk.RIGHT)

    def savefigure(self,figure):
        filename = filedialog.asksaveasfilename(title='Save JPG', filetypes=(("JPG files", "*.jpg"), ("all files", "*.*")))
        if len(filename) == 0:
            return
        figure.savefig(filename, dpi=300)

    def quit(self):
        self.root.destroy()
        quit()
