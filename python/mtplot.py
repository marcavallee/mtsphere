
import os
import numpy as np
import tkinter as tk
from tkinter import messagebox
import matplotlib
matplotlib.use('tkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
import matplotlib.pyplot as plt

class MTPlot:

    def __init__(self,main):
        self.main = main

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
                    x = columns[1]
                    y = columns[2]
                    appres[jf,i,j] = columns[4]
                    phase[jf,i,j] = columns[5]
                    imp[jf,i,j] = columns[6] + 1j * columns[7]
        for k in range(3):
            mf1file.readline()
        for jf in range(self.main.nf):
            columns = [ float(f) for f in mf1file.readline().split() ]
            Kzx[jf] = columns[4] + 1j * columns[5]
            Kzy[jf] = columns[6] + 1j * columns[7]
        mf1file.close()
        if abs(x) > 0 and abs(y) > 0:
            diagonal = True
        else:
            diagonal = False
        figure, axes = plt.subplots(4,1,figsize=(8,8))
        axes[0].plot(freq,appres[:,0,1],'b',label='XY')
        axes[0].plot(freq,appres[:,1,0],'r',label='YX')
        axes[0].set_ylabel(r'Apparent resistivity ($\Omega \cdot m$)')
        axes[0].legend(loc='upper right')
        axes[0].set_xscale('log')
        axes[0].set_xlim(min(freq),max(freq))
        axes[1].plot(freq,phase[:,0,1],'b',label='XY')
        axes[1].plot(freq,phase[:,1,0]+180,'r',label='YX+180')
        axes[1].set_xlabel('Frequency (Hz)')
        axes[1].set_ylabel(r'Phase ($^\circ$)')
        axes[1].legend(loc='upper right')
        axes[1].set_xscale('log')
        axes[1].set_xlim(min(freq),max(freq))
        if diagonal:
            axes[2].plot(freq, appres[:, 0, 0], 'b', label='XX')
            axes[2].plot(freq, appres[:, 1, 1], 'r', label='YY')
            axes[2].set_ylabel(r'Apparent resistivity ($\Omega \cdot m$)')
            axes[2].legend(loc='upper right')
            axes[2].set_xscale('log')
            axes[2].set_xlim(min(freq),max(freq))
            axes[3].plot(freq, phase[:, 0, 0], 'b', label='XX')
            axes[3].plot(freq, phase[:, 1, 1] + 180, 'r', label='YY+180')
            axes[3].set_xlabel('Frequency (Hz)')
            axes[3].set_ylabel(r'Phase ($^\circ$)')
            axes[3].legend(loc='upper right')
            axes[3].set_xscale('log')
            axes[3].set_xlim(min(freq),max(freq))
        else:
            axes[2].plot(freq,100*Kzx.real,'b',label='Real')
            axes[2].plot(freq,100*Kzx.imag,'r',label='Imaginary')
            axes[2].set_xlabel('Frequency (Hz)')
            axes[2].set_ylabel('X Tipper ratio (%)')
            axes[2].legend(loc='upper right')
            axes[2].set_xscale('log')
            axes[2].set_xlim(min(freq),max(freq))
            ylim = max(np.max(100*np.abs(Kzx)),1)
            axes[2].set_ylim(-ylim,ylim)

            axes[3].plot(freq,100*Kzy.real,'b',label='Real')
            axes[3].plot(freq,100*Kzy.imag,'r',label='Imaginary')
            axes[3].set_xlabel('Frequency (Hz)')
            axes[3].set_ylabel('Y Tipper ratio (%)')
            axes[3].legend(loc='upper right')
            ylim = max(np.max(100*np.abs(Kzy)),1)
            axes[3].set_ylim(-ylim,ylim)
            axes[3].set_xscale('log')
            axes[3].set_xlim(min(freq),max(freq))

        master = tk.Toplevel()
        canvasplot = FigureCanvasTkAgg(figure, master=master)
        canvasplot.draw()
        toolbar = NavigationToolbar2Tk(canvasplot, master)
        toolbar.update()
        canvasplot.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        canvasplot._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        toolbar.pack(side=tk.TOP, fill=tk.X)
        button = tk.Button(master, text='Quit', command=master.destroy)
        button.pack(side=tk.RIGHT)

    def plotmtprofiles(self):
        if not self.main.readcflfile():
            return
        if self.main.nf > 1:
            messagebox.showerror('MTSphere', 'More than one frequency')
            return
        mf1file = open(os.path.join(self.main.directory, 'MTSphere.mf1'), 'r')
        if self.main.nx > 1:
            if self.main.ny > 1 or self.main.nd > 1:
                messagebox.showerror('MTSphere', 'More than one distance variation')
                return
            profile = 'X'
        elif self.main.ny > 1:
            if self.main.nd > 1:
                messagebox.showerror('MTSphere', 'More than one distance variation')
                return
            profile = 'Y'
        else:
            profile = 'Z'
        if profile == 'X':
            nu = self.main.nx
            pu = np.zeros(self.main.nx)
            xlabel = 'X Distance (m)'
        elif profile == 'Y':
            nu = self.main.ny
            pu = np.zeros(self.main.ny)
            xlabel = 'Y Distance (m)'
        elif profile == 'Z':
            nu = self.main.nd
            pu = np.zeros(self.main.nd)
            xlabel = 'Depth (m)'
        appres = np.zeros((nu,2,2))
        phase = np.zeros((nu,2,2))
        imp = np.zeros((nu,2,2), dtype=complex)
        Kzx = np.zeros(nu, dtype=complex)
        Kzy = np.zeros(nu, dtype=complex)
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    mf1file.readline()
                for ju in range(nu):
                    columns = [ float(f) for f in mf1file.readline().split() ]
                    if profile == 'X':
                        pu[ju] = columns[1]
                    elif profile == 'Y':
                        pu[ju] = columns[2]
                    elif profile == 'Z':
                        pu[ju] = columns[3]
                    appres[ju,i,j] = columns[4]
                    phase[ju,i,j] = columns[5] % 360
                    imp[ju,i,j] = columns[6] + 1j * columns[7]
        for k in range(3):
            mf1file.readline()
        for ju in range(nu):
            columns = [ float(f) for f in mf1file.readline().split() ]
            Kzx[ju] = columns[4] + 1j * columns[5]
            Kzy[ju] = columns[6] + 1j * columns[7]
        mf1file.close()
        figure, axes = plt.subplots(4,1,figsize=(8,8))
        axes[0].plot(pu,appres[:,0,1],'b',label='XY')
        axes[0].plot(pu,appres[:,1,0],'r',label='YX')
        axes[0].set_ylabel(r'Apparent resistivity ($\Omega \cdot m$)')
        axes[0].legend(loc='upper left')
        ax2 = axes[0].twinx()
        ax2.plot(pu,appres[:,0,0],'b-+',label='XX')
        ax2.plot(pu,appres[:,1,1],'r-+',label='YY')
        ax2.set_ylabel(r'Apparent resistivity ($\Omega \cdot m$)')
        ax2.legend(loc='upper right')
        axes[1].plot(pu,phase[:,0,1],'b',label='XY')
        axes[1].plot(pu,(phase[:,1,0]+180) % 360,'r',label='YX+180')
        axes[1].set_ylabel(r'Phase ($^\circ$)')
        axes[1].legend(loc='upper left')
        ax2 = axes[1].twinx()
        ax2.plot(pu,phase[:,0,0],'b-+',label='XX')
        ax2.plot(pu,(phase[:,1,1]+180) % 360,'r-+',label='YY+180')
        ax2.set_ylabel(r'Phase ($^\circ$)')
        ax2.legend(loc='upper right')
        axes[2].plot(pu,100*Kzx.real,'b',label='Real')
        axes[2].plot(pu,100*Kzx.imag,'r',label='Imaginary')
        axes[2].set_ylabel('X Tipper ratio (%)')
        axes[2].legend(loc='upper right')
        axes[3].plot(pu,100*Kzy.real,'b',label='Real')
        axes[3].plot(pu,100*Kzy.imag,'r',label='Imaginary')
        axes[3].set_xlabel(xlabel)
        axes[3].set_ylabel('Y Tipper ratio (%)')
        axes[3].legend(loc='upper right')
        if abs(np.max([100*Kzx.real,100*Kzy.imag])) < 1:
            axes[3].set_ylim(-1,1)

        master = tk.Toplevel()
        canvasplot = FigureCanvasTkAgg(figure, master=master)
        canvasplot.draw()
        toolbar = NavigationToolbar2Tk(canvasplot, master)
        toolbar.update()
        canvasplot.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        canvasplot._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        toolbar.pack(side=tk.TOP, fill=tk.X)
        button = tk.Button(master, text='Quit', command=master.destroy)
        button.pack(side=tk.RIGHT)

    def plotfieldprofiles(self):
        if not self.main.readcflfile():
            return
        if sum([x>1 for x in (self.main.nx, self.main.ny, self.main.nd)])>=2:
            messagebox.showerror('MTSphere', 'More than one series longuer than 1')
            return
        mf2file = open(os.path.join(self.main.directory, 'MTSphere.mf2'), 'r')
        px = []
        py = []
        pd = []
        es = np.zeros((self.main.nx*self.main.ny*self.main.nd,3,2),dtype=complex)
        hs = np.zeros( (self.main.nx*self.main.ny*self.main.nd,3,2),dtype=complex)
        et = np.zeros((self.main.nx*self.main.ny*self.main.nd,3,2),dtype=complex)
        ht = np.zeros( (self.main.nx*self.main.ny*self.main.nd,3,2),dtype=complex)
        for k in range(5):
            mf2file.readline()
        for j in range(self.main.nx*self.main.ny*self.main.nd):
            columns = [ float(f) for f in mf2file.readline().split() ]
            px.append(columns[1])
            py.append(columns[2])
            pd.append(columns[3])
            es[j,0,0] = columns[4] + 1j * columns[5]
            es[j,1,0] = columns[6] + 1j * columns[7]
            es[j,2,0] = columns[8] + 1j * columns[9]
            hs[j,0,0] = columns[10] + 1j * columns[11]
            hs[j,1,0] = columns[12] + 1j * columns[13]
            hs[j,2,0] = columns[14] + 1j * columns[15]
        for k in range(2):
            mf2file.readline()
        for j in range(self.main.nx*self.main.ny*self.main.nd):
            columns = [ float(f) for f in mf2file.readline().split() ]
            et[j,0,0] = columns[4] + 1j * columns[5]
            et[j,1,0] = columns[6] + 1j * columns[7]
            et[j,2,0] = columns[8] + 1j * columns[9]
            ht[j,0,0] = columns[10] + 1j * columns[11]
            ht[j,1,0] = columns[12] + 1j * columns[13]
            ht[j,2,0] = columns[14] + 1j * columns[15]
        for j in range(5):
            mf2file.readline()
        for j in range(self.main.nx*self.main.ny*self.main.nd):
            columns = [ float(f) for f in mf2file.readline().split() ]
            es[j,0,1] = columns[4] + 1j * columns[5]
            es[j,1,1] = columns[6] + 1j * columns[7]
            es[j,2,1] = columns[8] + 1j * columns[9]
            hs[j,0,1] = columns[10] + 1j * columns[11]
            hs[j,1,1] = columns[12] + 1j * columns[13]
            hs[j,2,1] = columns[14] + 1j * columns[15]
        for k in range(2):
            mf2file.readline()
        for j in range(self.main.nx*self.main.ny*self.main.nd):
            columns = [ float(f) for f in mf2file.readline().split() ]
            et[j,0,1] = columns[4] + 1j * columns[5]
            et[j,1,1] = columns[6] + 1j * columns[7]
            et[j,2,1] = columns[8] + 1j * columns[9]
            ht[j,0,1] = columns[10] + 1j * columns[11]
            ht[j,1,1] = columns[12] + 1j * columns[13]
            ht[j,2,1] = columns[14] + 1j * columns[15]
        mf2file.close()
        e = et
        h = ht
        px = list(dict.fromkeys(px))
        py = list(dict.fromkeys(py))
        pd = list(dict.fromkeys(pd))
        if len(px) > 1:
            label = 'X Distance (m)'
            ps = px
        elif len(py) > 1:
            label = 'Y Distance (m)'
            ps = py
        elif len(pd) > 1:
            label = 'Depth (m)'
            ps = pd
        else:
            messagebox.showerror('No length variable')
            return

        figure, axes = plt.subplots(2, 6, figsize=(18, 10))
        for axe in axes.flatten():
            axe.tick_params(axis='y', direction='in', pad=-45)
        def plotfigures(row,column,f, flabel):

            axes[row,3*column].plot(ps,f[:,0].real,'b',label='real')
            axes[row,3*column].plot(ps,f[:,0].imag,'r',label='imag')
            axes[row,3*column].set_ylabel(flabel+'x ')
            axes[row,3*column].set_xlabel(label)
            axes[row,3*column].legend(loc='upper right')

            axes[row,1+3*column].plot(ps,f[:,1].real,'b',label='real')
            axes[row,1+3*column].plot(ps,f[:,1].imag,'r',label='imag')
            axes[row,1+3*column].set_xlabel(label)
            axes[row,1+3*column].set_ylabel(flabel+'y')
            axes[row,1+3*column].legend(loc='upper right')
            if row == 0:
                if column == 0:
                    axes[0,1].set_title('X polarization')
                else:
                    axes[0,4].set_title('Y polarization')

            axes[row,2+3*column].plot(ps, f[:, 2].real, 'b', label='real')
            axes[row,2+3*column].plot(ps, f[:, 2].imag, 'r', label='imag')
            axes[row,2+3*column].set_xlabel(label)
            axes[row,2+3*column].set_ylabel(flabel+'z')
            axes[row,2+3*column].legend(loc='upper right')

        plotfigures(0,0, e[:,:,0],'E')
        plotfigures(1,0, h[:,:,0],'H')
        plotfigures(0,1, e[:,:,1],'E')
        plotfigures(1,1, h[:,:,1],'H')

        master = tk.Toplevel()
        canvasplot = FigureCanvasTkAgg(figure, master=master)
        canvasplot.draw()
        toolbar = NavigationToolbar2Tk(canvasplot, master)
        toolbar.update()
        canvasplot.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        canvasplot._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        toolbar.pack(side=tk.TOP, fill=tk.X)
        button = tk.Button(master, text='Quit', command=master.destroy)
        button.pack(side=tk.RIGHT)

    def plotfieldplane(self):
        if not self.main.readcflfile():
            return
        sum = np.sum([x>1 for x in (self.main.nx, self.main.ny, self.main.nd)])
        if sum < 2:
            messagebox.showerror('MTSphere', 'Only one series is longuer than 1')
            return
        if sum > 2:
            messagebox.showerror('MTSphere', 'More than two series longuer than 1')
            return
        mf2file = open(os.path.join(self.main.directory, 'MTSphere.mf2'), 'r')
        px = []
        py = []
        pd = []
        if self.main.ny > self.main.nd:
            ny = self.main.ny
        else:
            ny = self.main.nd
        es = np.zeros((self.main.nx,ny,3),dtype=complex)
        hs = np.zeros( (self.main.nx,ny,3),dtype=complex)
        et = np.zeros((self.main.nx,ny,3),dtype=complex)
        ht = np.zeros( (self.main.nx,ny,3),dtype=complex)
        for k in range(5):
            mf2file.readline()
        for j in range(self.main.nx):
            for k in range(ny):
                columns = [ float(f) for f in mf2file.readline().split() ]
                px.append(columns[1])
                py.append(columns[2])
                pd.append(columns[3])
                es[j,k,0] = columns[4] + 1j * columns[5]
                es[j,k,1] = columns[6] + 1j * columns[7]
                es[j,k,2] = columns[8] + 1j * columns[9]
        for k in range(2):
            mf2file.readline()
        for j in range(self.main.nx):
            for k in range(ny):
                columns = [float(f) for f in mf2file.readline().split()]
                et[j,k,0] = columns[4] + 1j * columns[5]
                et[j,k,1] = columns[6] + 1j * columns[7]
                et[j, k, 2] = columns[8] + 1j * columns[9]
        for j in range(5):
            mf2file.readline()
        for j in range(self.main.nx):
            for k in range(ny):
                columns = [float(f) for f in mf2file.readline().split()]
                hs[j,k,0] = columns[10] + 1j * columns[11]
                hs[j,k,1] = columns[12] + 1j * columns[13]
                hs[j,k,2] = columns[14] + 1j * columns[15]
        for k in range(2):
            mf2file.readline()
        for j in range(self.main.nx):
            for k in range(ny):
                columns = [float(f) for f in mf2file.readline().split()]
                ht[j, k, 0] = columns[10] + 1j * columns[11]
                ht[j, k, 1] = columns[12] + 1j * columns[13]
                ht[j, k, 2] = columns[14] + 1j * columns[15]

        mf2file.close()
        e = et
        h = ht
        px = list(dict.fromkeys(px))
        py = list(dict.fromkeys(py))
        pd = list(dict.fromkeys(pd))

        # 1. Origines des vecteurs
        if len(pd) > len(py):
            Y = np.array(pd)
            X = np.array(px)
            xlabel = 'X Distance (m)'
            ylabel = 'Depth (m)'
        else:
            X = np.array(px)
            Y = np.array(py)
            xlabel = 'X Distance (m)'
            ylabel = 'Y Distance (m)'
        z = 0.0

        taille_fleche = X[1] - X[0]
        # 3. Calcul automatique de l'amplitude (Norme)
        E_real_amplitude = 100 * np.sqrt(e[:,:,0].real ** 2 + e[:,:,2].real ** 2)
        E_imag_amplitude = 100 * np.sqrt(e[:,:,0].imag ** 2 + e[:,:,2].imag ** 2)
        H_real_amplitude = 100 * np.sqrt(h[:,:,0].real ** 2 + h[:,:,2].real ** 2)
        H_imag_amplitude = 100 * np.sqrt(h[:,:,0].imag ** 2 + h[:,:,2].imag ** 2)

        # 1. Normalisation des composantes pour le Champ Électrique (E)
        U_E_real = np.divide(e[:, :, 0].real, E_real_amplitude, out=np.zeros_like(E_real_amplitude),
                             where=E_real_amplitude > 0) * taille_fleche
        V_E_real = np.divide(e[:, :, 2].real, E_real_amplitude, out=np.zeros_like(E_real_amplitude),
                             where=E_real_amplitude > 0) * taille_fleche
        U_E_imag = np.divide(e[:, :, 0].imag, E_imag_amplitude, out=np.zeros_like(E_imag_amplitude),
                             where=E_imag_amplitude > 0) * taille_fleche
        V_E_imag = np.divide(e[:, :, 2].imag, E_imag_amplitude, out=np.zeros_like(E_imag_amplitude),
                             where=E_imag_amplitude > 0) * taille_fleche
        color_E_real = np.full_like(E_real_amplitude, np.nan)
        color_E_imag = np.full_like(E_imag_amplitude, np.nan)
        np.log10 ( E_real_amplitude, where=(E_real_amplitude>0), out=color_E_real)
        np.log10 ( E_imag_amplitude, where=(E_imag_amplitude>0), out=color_E_imag)

        # 2. Normalisation des composantes pour le Champ Magnétique (H)
        U_H_real = np.divide(h[:, :, 0].real, H_real_amplitude, out=np.zeros_like(H_real_amplitude),
                             where=H_real_amplitude > 0) * taille_fleche
        V_H_real = np.divide(h[:, :, 2].real, H_real_amplitude, out=np.zeros_like(H_real_amplitude),
                             where=H_real_amplitude > 0) * taille_fleche
        U_H_imag = np.divide(h[:, :, 0].imag, H_imag_amplitude, out=np.zeros_like(H_imag_amplitude),
                             where=H_imag_amplitude > 0) * taille_fleche
        V_H_imag = np.divide(h[:, :, 2].imag, H_imag_amplitude, out=np.zeros_like(H_imag_amplitude),
                             where=H_imag_amplitude > 0) * taille_fleche
        color_H_real = np.full_like(H_real_amplitude, np.nan)
        color_H_imag = np.full_like(H_imag_amplitude, np.nan)
        np.log10 ( H_real_amplitude, where=(H_real_amplitude>0), out=color_H_real)
        np.log10 ( H_imag_amplitude, where=(H_imag_amplitude>0), out=color_H_imag)

        # 4. Affichage avec colormap (la couleur change selon l'amplitude)
        fig, axes = plt.subplots(2,2,figsize=(12, 12))
        im0 = axes[0, 0].pcolormesh(X, Y, color_E_real.T, cmap='jet', shading='nearest')
        axes[0,0].quiver(X, Y, U_E_real, V_E_real,
                                 angles='xy', scale_units='xy', scale=1, color='black')
        # Ajouter une barre de couleur pour l'amplitude
        fig.colorbar(im0, ax=axes[0,0], label='log10(E real)')
        axes[0,0].set_title('E field, X polarization')
        axes[0,0].set_ylabel(ylabel)

        im1 = axes[1,0].pcolormesh(X, Y, color_E_imag.T, cmap='jet', shading='nearest')
        axes[1,0].quiver(X, Y, U_E_imag , V_E_imag,
                                 angles='xy', scale_units='xy', scale=1, color='black')
        # Ajouter une barre de couleur pour l'amplitude
        fig.colorbar(im1, ax=axes[1,0], label='log10(E imag)')
        axes[1,0].set_ylabel(ylabel)
        axes[1,0].set_xlabel(xlabel)

        im2 = axes[0,1].pcolormesh(X, Y, color_H_real.T, cmap='jet', shading='nearest')
        axes[0,1].quiver(X, Y, U_H_real, V_H_real,
                                 angles='xy', scale_units='xy', scale=1, color='black')
        # Ajouter une barre de couleur pour l'amplitude
        fig.colorbar(im2, ax=axes[0,1], label='log10(H real)')
        axes[0,1].set_title('H field, Y polarization')

        im3 = axes[1, 1].pcolormesh(X, Y, color_H_imag.T, cmap='jet', shading='nearest')
        axes[1,1].quiver(X, Y, U_H_imag , V_H_imag,
                                 angles='xy', scale_units='xy', scale=1, color='black')
        # Ajouter une barre de couleur pour l'amplitude
        fig.colorbar(im3, ax=axes[1,1], label='log10(H imag)')
        axes[1,1].set_xlabel(xlabel)

        if self.main.nd > self.main.ny:
            for axe in axes.flatten():
                axe.invert_yaxis()

        master = tk.Toplevel()
        canvasplot = FigureCanvasTkAgg(fig, master=master)
        canvasplot.draw()
        toolbar = NavigationToolbar2Tk(canvasplot, master)
        toolbar.update()
        canvasplot.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        canvasplot._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        toolbar.pack(side=tk.TOP, fill=tk.X)
        button = tk.Button(master, text='Quit', command=master.destroy)
        button.pack(pady=5)
