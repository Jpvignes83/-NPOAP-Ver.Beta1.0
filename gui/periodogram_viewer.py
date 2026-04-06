import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog, ttk
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np 

class PeriodogramViewer:
    def __init__(self, root):
        self.root = root
        if isinstance(self.root, tk.Tk) or isinstance(self.root, tk.Toplevel):
            self.root.title("Périodogramme Viewer")

        self.frame = ttk.Frame(self.root)
        self.frame.pack(fill=tk.BOTH, expand=True)

        self.markers = [] # Liste pour stocker les valeurs Périodes (floats)
        self.data = None # DataFrame pandas

        # Barre d'outils
        control_frame = ttk.Frame(self.frame)
        control_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=10)

        
        ttk.Button(control_frame, text="📍 Ajouter Marqueur (Manuel)", command=self.add_marker).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="🗑️ Vider Marqueurs", command=self.clear_markers).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="💾 Sauvegarder Figure", command=self.save_figure).pack(side=tk.LEFT, padx=5)
        

        # Matplotlib Canvas
        self.figure, self.ax = plt.subplots(figsize=(6, 4))
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Connexion du clic souris
        self.figure.canvas.mpl_connect("button_press_event", self.on_click)

    def _navigation_mode_active(self) -> bool:
        """
        Retourne True si un mode de navigation Matplotlib (zoom/pan) est actif.
        Dans ce cas on ignore les clics métier (ajout/suppression de marqueur).
        """
        toolbar = getattr(self.figure.canvas, "toolbar", None)
        if toolbar is None:
            return False
        mode = getattr(toolbar, "mode", None) or getattr(toolbar, "_active", None)
        if mode is None:
            return False
        mode_txt = str(mode).strip().upper()
        return mode_txt not in ("", "NONE")
    
    def _log_markers(self):
        """Imprime la liste triée des marqueurs de période dans la console."""
        if self.markers:
            # Tri et formatage des valeurs
            sorted_markers = sorted(self.markers)
            marker_list = [f"{p:.4f} j" for p in sorted_markers]
            print(f"MARQUEURS ACTUELS: {marker_list}")
        else:
            print("MARQUEURS ACTUELS: (Aucun)")


    def plot_periodogram(self, period=None, power=None, title="Périodogramme"):
        """
        Affiche le périodogramme et ré-affiche les marqueurs avec leurs valeurs numériques
        sur le graphique.
        """
        # Si de nouvelles données arrivent (depuis l'onglet Analyse), on écrase les anciennes
        if period is not None and power is not None:
            self.data = pd.DataFrame({"period": period, "power": power})
            # On conserve self.markers

        # Si aucune donnée n'est disponible
        if self.data is None:
            return

        self.ax.clear()
        
        # Tracé principal
        self.ax.plot(self.data["period"], self.data["power"], color="black", linewidth=1)
        
        # Cosmétique
        self.ax.set_xlabel("Période (jours)")
        self.ax.set_ylabel("Puissance")
        self.ax.set_title(title)
        self.ax.grid(True, linestyle=':', alpha=0.6)

        # Ré-affichage des marqueurs et AJOUT DE LEURS VALEURS NUMÉRIQUES SUR LE GRAPHIQUE
        
        # Récupération des limites Y et calcul de la position du texte (95% de la hauteur)
        y_lim = self.ax.get_ylim()
        y_text_pos = y_lim[1] * 0.95 
        
        for marker_period in self.markers:
            # 1. Dessin de la ligne verticale (Marqueur)
            self.ax.axvline(marker_period, color="red", linestyle="--", linewidth=2, zorder=3)
            
            # 2. Ajout du texte de la valeur de période
            self.ax.text(
                marker_period, 
                y_text_pos, 
                f'{marker_period:.4f} j', 
                color='red', 
                ha='left', 
                va='top',  
                rotation=0, # Texte vertical pour la lisibilité
                fontsize=12,
                # Boîte pour assurer la visibilité sur la courbe
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.2'), 
                zorder=10 
            )
            
        self.canvas.draw()

    
    def on_click(self, event):
        # Neutraliser les actions souris si mode zoom/pan actif
        if self._navigation_mode_active():
            return

        # 1. Vérification essentielle : le clic doit être dans les axes
        if event.inaxes != self.ax:
            return
        
        # 2. Récupérer la valeur X (Période ou Époque)
        # N'UTILISEZ JAMAIS round() SUR LES DONNÉES DE TEMPS!
        x_value = event.xdata 
        
        if x_value is None:
            return
            
        # VÉRIFICATION DE LA DISPONIBILITÉ DES DONNÉES
        if self.data is None:
            return

        current_title = self.ax.get_title()

        # Clic gauche (Ajout de marqueur) : Button 1
        if event.button == 1:
            self.markers.append(x_value)
            
            # Appel de la méthode de traçage appropriée (assurez-vous d'avoir une méthode de traçage pour le TTVViewer)
            if hasattr(self, 'plot_periodogram'):
                 self.plot_periodogram(title=current_title)
            elif hasattr(self, 'plot_ttv'):
                 # Si c'est le TTVViewer, vous pouvez appeler la méthode de traçage TTV.
                 # Sinon, appelez la méthode qui affiche l'O-C (comme plot_external_data)
                 self.plot_ttv(ttv_col='o-c', title=current_title) # Exemple
                 
            self._log_markers() 
            
        # Clic droit (Suppression de marqueur) : Button 3
        elif event.button == 3:
            if not self.markers:
                messagebox.showinfo("Suppression", "Aucun marqueur à supprimer.")
                return
            
            # Logique pour trouver le marqueur le plus proche et le supprimer
            markers_array = np.array(self.markers)
            closest_index = np.argmin(np.abs(markers_array - x_value))
            period_to_remove = markers_array[closest_index]
            
            # Tolérance de clic (1% de la plage X)
            x_range = self.ax.get_xlim()[1] - self.ax.get_xlim()[0]
            tolerance = x_range * 0.01 
            
            if np.abs(period_to_remove - x_value) < tolerance:
                self.markers.pop(closest_index)
                
                # Appel de la méthode de traçage appropriée
                if hasattr(self, 'plot_periodogram'):
                     self.plot_periodogram(title=current_title)
                elif hasattr(self, 'plot_ttv'):
                     self.plot_ttv(ttv_col='o-c', title=current_title) # Exemple
                     
                messagebox.showinfo("Suppression", f"Marqueur à P/E={period_to_remove:.4f} supprimé.")
                self._log_markers()
            else:
                 print("Clic droit ignoré : trop éloigné d'un marqueur existant.")
                 
    def add_marker(self):
        value = simpledialog.askfloat("Ajouter Marqueur", "Entrer la période à marquer :")
        if value is not None:
            self.markers.append(value)
            current_title = self.ax.get_title()
            self.plot_periodogram(title=current_title)
            self._log_markers() # <-- Appel pour imprimer la liste

    def clear_markers(self):
        """Vide la liste des marqueurs et redessine."""
        if messagebox.askyesno("Confirmer", "Voulez-vous vraiment supprimer tous les marqueurs ?"):
            self.markers = []
            current_title = self.ax.get_title()
            self.plot_periodogram(title=current_title)
            self._log_markers() # <-- Appel pour imprimer la liste (vide)

    def save_figure(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".png")
        if file_path:
            self.figure.savefig(file_path)
            messagebox.showinfo("Succès", f"Figure sauvegardée : {file_path}")