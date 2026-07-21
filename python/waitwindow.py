import sys
import tkinter as tk
from tkinter import ttk
import time
from tkinter import messagebox

class WaitWindow:
    """Fenêtre d'attente simple et réutilisable"""
    def __init__(self, parent=None, title="Calcul en cours",
                 message="Please wait.\nComputation in progress."):
        if parent is None:
            parent = tk._default_root
        self.root = tk.Toplevel(parent)
        self.title = title
        self.root.title(self.title)
        self.root.geometry("420x180")
        self.root.resizable(False, False)
        self.root.transient()  # Reste au-dessus de la fenêtre principale
        self.root.grab_set()  # Rend la fenêtre modale

        # Centrer la fenêtre à l'écran
        self.root.update_idletasks()
        x = (self.root.winfo_screenwidth() // 2) - (self.root.winfo_width() // 2)
        y = (self.root.winfo_screenheight() // 2) - (self.root.winfo_height() // 2)
        self.root.geometry(f"+{x}+{y}")

        # Message
        tk.Label(self.root, text=message, font=("Arial", 14), pady=20, justify="center").pack()

        # Minuterie
        self.start_time = time.time()
        self.time_label = tk.Label(self.root, text="Elapsed time : 00:00:00",
                                  font=("Arial", 14, "bold"), fg="blue")
        self.time_label.pack(pady=5)
        self.elapsed = 0

#        # Barre de progression indéterminée
#        self.progress = ttk.Progressbar(self.root, mode='indeterminate', length=320)
#        self.progress.pack(pady=10)
#        self.progress.start(12)  # Vitesse d'animation

        # Petit texte en bas

        btn = tk.Button(self.root, text='Quit', command=self.quit)
        btn.pack(side='right')

        self.update_timer()

        self.root.update()

    def update_timer(self):
        """Met à jour la minuterie chaque seconde"""
        if not self.root.winfo_exists():
            return
        self.elapsed = int(time.time() - self.start_time)
        heures = self.elapsed // 3600
        minutes = ( self.elapsed % 3600 ) // 60
        secondes = self.elapsed % 60
        self.time_label.config(text=f"Elapsed time : {heures:02d}:{minutes:02d}:{secondes:02d}")

        # Rappel toutes les 1000 ms (1 seconde)
        self.root.after(1000, self.update_timer)

    def quit(self):
        """Ferme complètement le programme"""
        if messagebox.askyesno(self.title,
                               "Do you want to stop the program?",
                               icon='warning'):
            self.root.destroy()
            sys.exit(0)  # ← Fermeture complète du programme

    def destroy(self):
        """Ferme la fenêtre"""
        try:
#            self.progress.stop()
            self.root.destroy()
        except:
            pass