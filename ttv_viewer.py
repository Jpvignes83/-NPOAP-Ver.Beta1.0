import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter.filedialog import asksaveasfilename
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import logging
try:
    from scipy.interpolate import interp1d
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

logger = logging.getLogger(__name__)

class TTVViewer:
    def __init__(self, root):
        self.root = root
        if isinstance(self.root, tk.Tk):
            self.root.title("TTV Viewer")

        self.frame = ttk.Frame(self.root)
        self.frame.pack(fill=tk.BOTH, expand=True)

        # CRÉATION DE LA FIGURE avec 2 panneaux (données + résidus)
        self.fig = plt.Figure(figsize=(10, 7))
        gs = self.fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.3)
        self.ax = self.fig.add_subplot(gs[0])  # Panneau principal
        self.ax_res = self.fig.add_subplot(gs[1])  # Panneau résidus
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.canvas.mpl_connect("button_press_event", self.on_click)

        self.data = None
        self.markers = []

    # --- Méthodes de l'interface graphique ---
    
    def load_csv(self, file_path=None): # Rendre la méthode plus flexible
        if file_path is None:
             file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        
        if not file_path:
            return
        try:
            self.data = pd.read_csv(file_path)
            self.data.columns = [c.strip().lower() for c in self.data.columns]
            
            if 'epoch' not in self.data.columns:
                raise ValueError("Le fichier doit contenir la colonne 'epoch'.")
            
            ttv_col = next((c for c in self.data.columns if c in ['ttv', 'o-c', 'oc']), None)
            if ttv_col is None:
                 raise ValueError("Le fichier doit contenir la colonne 'TTV' ou 'O-C'.")
                 
            # NOTE : load_csv est ici seulement pour un usage standalone/test. 
            # L'application principale utilise plot_external_data
            self.plot_ttv(ttv_col) 
            messagebox.showinfo("Succès", f"Fichier chargé: {file_path}")
            
        except Exception as e:
            messagebox.showerror("Erreur Chargement", str(e))

    def plot_external_data(self, epoch_data, oc_data, yerr_data, model_x=None, model_y=None):
        """
        Trace les données O-C (points), le modèle (ligne) et les barres d'erreur.
        """
        self.ax.clear()
        
        # --- FIX MAJEUR 1: ASSURER LE TYPE float64 pour le tracé ---
        try:
            epochs = np.asarray(epoch_data, dtype=np.float64)
            oc = np.asarray(oc_data, dtype=np.float64)
            # Yerr peut être None ou un array, le caster si possible
            yerr = np.asarray(yerr_data, dtype=np.float64) if yerr_data is not None and len(yerr_data) == len(oc) else None
            
            valid_data_len = len(epochs)
        except Exception as e:
            logger.error(f"Erreur de conversion de données pour le tracé: {e}", exc_info=True)
            self.canvas.draw()
            return
        
        # 1. Tracé de la ligne de référence (Zéro) - amélioré
        self.ax.axhline(y=0, color='#d62728', linestyle='--', linewidth=2.0, alpha=0.8, zorder=1, label='Réf. Linéaire (T0)')

        # 2. Tracé du modèle sinusoïdal (si présent) - tracé en premier pour être sous les points
        if model_x is not None and model_y is not None and len(model_x) > 0:
             # Augmenter le nombre de points pour une courbe plus lisse si nécessaire
             if len(model_x) < 500 and SCIPY_AVAILABLE:
                 try:
                     # Interpolation pour courbe plus lisse
                     f_interp = interp1d(model_x, model_y, kind='cubic', bounds_error=False, fill_value='extrapolate')
                     model_x_smooth = np.linspace(model_x.min(), model_x.max(), 500)
                     model_y_smooth = f_interp(model_x_smooth)
                     self.ax.plot(model_x_smooth, model_y_smooth, color='#2ca02c', linestyle='-', linewidth=2.5, 
                                 alpha=0.9, zorder=2, label='Modèle Fit')
                 except Exception:
                     # En cas d'erreur d'interpolation, utiliser les points originaux
                     self.ax.plot(model_x, model_y, color='#2ca02c', linestyle='-', linewidth=2.5, 
                                 alpha=0.9, zorder=2, label='Modèle Fit')
             else:
                 self.ax.plot(model_x, model_y, color='#2ca02c', linestyle='-', linewidth=2.5, 
                             alpha=0.9, zorder=2, label='Modèle Fit')

        # 3. Tracé des données observées (Points avec barres d'erreur) - amélioré
        if valid_data_len > 0:
            self.ax.errorbar(
                epochs, 
                oc,     
                yerr=yerr,
                fmt='o',
                color='black',
                markersize=6,
                markerfacecolor='white',
                markeredgewidth=1.5,
                markeredgecolor='black',
                capsize=4,
                capthick=1.5,
                elinewidth=1.5,
                alpha=0.8,
                zorder=3,
                label='O-C Observés'
            )
            
            # Stockage temporaire des données
            self.data = pd.DataFrame({'epoch': epochs, 'o-c': oc})
            if yerr is not None:
                self.data['uncertainty'] = yerr

        # 4. TRACÉ DES MARQUEURS MANUELS
        for marker in self.markers:
             self.ax.axvline(marker, color='purple', linestyle=':', alpha=0.5)

        # 5. Définition des Axes
        min_epoch = None
        max_epoch = None
        if valid_data_len > 0:
            min_epoch = np.min(epochs)
            max_epoch = np.max(epochs)
            
            # AXE X: Plage correcte (du min au max)
            self.ax.set_xlim(min_epoch - 5, max_epoch + 10)
            
            # --- ADAPTATION DYNAMIQUE DE L'AXE Y BASÉE SUR L'AMPLITUDE ---
            # Si un modèle est présent, utiliser ses min/max pour définir les limites
            if model_x is not None and model_y is not None and len(model_y) > 0:
                # Utiliser l'amplitude réelle du modèle
                y_min_model = np.min(model_y)
                y_max_model = np.max(model_y)
                y_center = (y_min_model + y_max_model) / 2.0
                y_range_model = y_max_model - y_min_model
                
                # Ajouter une marge de 15% de chaque côté
                margin = y_range_model * 0.15
                if margin < 0.001: margin = 0.001  # Marge minimale
                
                y_min_lim = y_min_model - margin
                y_max_lim = y_max_model + margin
                
                self.ax.set_ylim(y_min_lim, y_max_lim)
            else:
                # Pas de modèle : utiliser les données observées
                y_max_abs = np.max(np.abs(oc))
                y_range = y_max_abs * 1.5 
                if y_range < 0.005: y_range = 0.005 # Valeur minimale pour la visibilité
                
                self.ax.set_ylim(-y_range, y_range) # Centrer l'axe Y autour de zéro

        else:
            self.ax.set_xlim(0, 100) 
            self.ax.set_ylim(-0.01, 0.01) # Plage par défaut

        # Cosmétique du graphique - améliorée
        self.ax.set_xlabel("")  # Pas de label X sur le panneau supérieur (partagé avec résidus)
        self.ax.set_ylabel("O-C (jours)", fontsize=12, fontweight='bold')
        self.ax.set_title("Diagramme TTV", fontsize=14, fontweight='bold', pad=15)
        
        # Grille améliorée
        self.ax.grid(True, linestyle='--', alpha=0.4, linewidth=0.8, color='gray')
        self.ax.set_axisbelow(True)  # Grille sous les données
        
        # Améliorer les ticks
        self.ax.tick_params(axis='both', which='major', labelsize=10, width=1.2, length=5)
        self.ax.tick_params(axis='both', which='minor', labelsize=8, width=0.8, length=3)
        
        # Afficher la légende améliorée
        if self.ax.get_legend_handles_labels()[0]:
            legend = self.ax.legend(loc='upper right', frameon=True, fancybox=True, shadow=True, 
                                   fontsize=10, framealpha=0.95, edgecolor='gray')
            legend.get_frame().set_linewidth(1.2)
        
        # --- PANEL BAS : RESIDUS ---
        self.ax_res.clear()
        if model_x is not None and model_y is not None and len(model_x) > 0 and valid_data_len > 0:
            # Calculer les résidus aux positions des données observées
            model_at_data = np.interp(epochs, model_x, model_y)
            residuals = oc - model_at_data
            
            # Tracer les résidus
            self.ax_res.scatter(epochs, residuals, c='black', s=15, alpha=0.6, zorder=3)
            
            # Ligne de référence à zéro
            self.ax_res.axhline(y=0, color='red', linestyle='--', linewidth=1.5, alpha=0.7, zorder=1)
            
            # Calculer les limites Y pour les résidus
            y_res_max_abs = np.max(np.abs(residuals))
            y_res_range = y_res_max_abs * 1.5 if y_res_max_abs > 0 else 0.01
            if y_res_range < 0.001: y_res_range = 0.001
            
            self.ax_res.set_ylim(-y_res_range, y_res_range)
            self.ax_res.set_xlim(min_epoch - 5, max_epoch + 10)
            
            # Labels et style
            self.ax_res.set_xlabel("Epoch", fontsize=11, fontweight='bold')
            self.ax_res.set_ylabel("Résidus (jours)", fontsize=11, fontweight='bold')
            self.ax_res.grid(True, linestyle='--', alpha=0.4, linewidth=0.8, color='gray')
            self.ax_res.set_axisbelow(True)
            
            # Calculer et afficher le RMS des résidus
            rms_residuals = np.std(residuals)
            self.ax_res.text(0.02, 0.95, f'RMS = {rms_residuals:.5f} jours', 
                           transform=self.ax_res.transAxes, fontsize=9,
                           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            
            # Améliorer les ticks
            self.ax_res.tick_params(axis='both', which='major', labelsize=9, width=1.0, length=4)
        else:
            # Pas de modèle : afficher un panneau vide
            self.ax_res.set_xlabel("Epoch", fontsize=11, fontweight='bold')
            self.ax_res.set_ylabel("Résidus (jours)", fontsize=11, fontweight='bold')
            self.ax_res.text(0.5, 0.5, 'Aucun modèle ajusté', 
                           transform=self.ax_res.transAxes, fontsize=10,
                           ha='center', va='center', style='italic', color='gray')
            self.ax_res.grid(True, linestyle='--', alpha=0.4, linewidth=0.8, color='gray')

        # Rafraichissement
        self.canvas.draw()
    
    # Correction de la signature pour accepter un titre
    def plot_ttv(self, ttv_col='o-c', title=None): 
        """Tracé interne utilisé par load_csv et on_click."""
        self.ax.clear()
        
        # Ligne de référence à zéro - améliorée
        self.ax.axhline(y=0, color='#d62728', linestyle='--', linewidth=2.0, alpha=0.8, zorder=1, label='Réf. Linéaire (T0)')
        
        # Tracé des données chargées
        # Ajout des barres d'erreur pour plot_ttv, si les données sont stockées
        if self.data is not None and 'epoch' in self.data.columns and ttv_col in self.data.columns:
            epochs = self.data['epoch'].values
            oc = self.data[ttv_col].values
            yerr = self.data['uncertainty'].values if 'uncertainty' in self.data.columns else None

            self.ax.errorbar(epochs, oc, yerr=yerr, fmt='o', color='black', markersize=6,
                           markerfacecolor='white', markeredgewidth=1.5, markeredgecolor='black',
                           capsize=4, capthick=1.5, elinewidth=1.5, alpha=0.8, zorder=3)
            
            # Définition de l'axe Y pour la visibilité (copier la logique de plot_external_data)
            y_max_abs = np.max(np.abs(oc))
            y_range = y_max_abs * 1.5 
            if y_range < 0.005: y_range = 0.005
            self.ax.set_ylim(-y_range, y_range)
            
            # Définition de l'axe X
            min_epoch = np.min(epochs)
            max_epoch = np.max(epochs)
            self.ax.set_xlim(min_epoch - 5, max_epoch + 10)

        # Cosmétique - améliorée
        self.ax.set_xlabel("")  # Pas de label X sur le panneau supérieur (partagé avec résidus)
        self.ax.set_ylabel("O-C (jours)", fontsize=12, fontweight='bold')
        self.ax.set_title(title if title else "Diagramme TTV", fontsize=14, fontweight='bold', pad=15)
        self.ax.grid(True, linestyle='--', alpha=0.4, linewidth=0.8, color='gray')
        self.ax.set_axisbelow(True)
        
        # Améliorer les ticks
        self.ax.tick_params(axis='both', which='major', labelsize=10, width=1.2, length=5)
        self.ax.tick_params(axis='both', which='minor', labelsize=8, width=0.8, length=3)
        
        # Affichage des marqueurs
        for marker in self.markers:
             self.ax.axvline(marker, color='purple', linestyle=':', alpha=0.5, zorder=1)
        
        # --- PANEL BAS : RESIDUS (pour plot_ttv aussi) ---
        self.ax_res.clear()
        if self.data is not None and 'epoch' in self.data.columns and ttv_col in self.data.columns:
            epochs = self.data['epoch'].values
            oc = self.data[ttv_col].values
            
            # Si on a un modèle, calculer les résidus
            # Pour plot_ttv, on n'a généralement pas de modèle, donc on affiche juste un panneau vide
            self.ax_res.set_xlabel("Epoch", fontsize=11, fontweight='bold')
            self.ax_res.set_ylabel("Résidus (jours)", fontsize=11, fontweight='bold')
            min_epoch = np.min(epochs)
            max_epoch = np.max(epochs)
            self.ax_res.set_xlim(min_epoch - 5, max_epoch + 10)
            self.ax_res.set_ylim(-0.01, 0.01)  # Plage par défaut
            self.ax_res.text(0.5, 0.5, 'Aucun modèle ajusté', 
                           transform=self.ax_res.transAxes, fontsize=10,
                           ha='center', va='center', style='italic', color='gray')
            self.ax_res.grid(True, linestyle='--', alpha=0.4, linewidth=0.8, color='gray')
            self.ax_res.set_axisbelow(True)
            self.ax_res.tick_params(axis='both', which='major', labelsize=9, width=1.0, length=4)
             
        self.canvas.draw()

    def on_click(self, event):
        # 1. Vérification essentielle : le clic doit être dans les axes
        if event.inaxes != self.ax:
            return
        
        # 2. Récupérer la valeur X (Époque)
        x_epoch = event.xdata 
        
        if x_epoch is None or self.data is None:
            return

        current_title = self.ax.get_title()

        # Clic gauche (Ajout de marqueur) : Button 1
        if event.button == 1:
            self.markers.append(x_epoch)
            
            # Re-traçage
            ttv_col = next((c for c in self.data.columns if c in ['ttv', 'o-c', 'oc']), 'o-c')
            self.plot_ttv(ttv_col, title=current_title) # Passé le titre
            # Note: _log_markers n'est pas fourni, je l'ignore ou vous devez le définir
            # self._log_markers() 
            
        # Clic droit (Suppression de marqueur) : Button 3
        elif event.button == 3:
            if not self.markers:
                messagebox.showinfo("Suppression", "Aucun marqueur à supprimer.")
                return
            
            # Logique pour trouver le marqueur le plus proche et le supprimer
            markers_array = np.array(self.markers)
            closest_index = np.argmin(np.abs(markers_array - x_epoch)) 
            period_to_remove = markers_array[closest_index]
            
            # Tolérance de clic (1% de la plage X)
            x_range = self.ax.get_xlim()[1] - self.ax.get_xlim()[0]
            tolerance = x_range * 0.01 
            
            if np.abs(period_to_remove - x_epoch) < tolerance:
                self.markers.pop(closest_index)
                
                # Re-traçage
                ttv_col = next((c for c in self.data.columns if c in ['ttv', 'o-c', 'oc']), 'o-c')
                self.plot_ttv(ttv_col, title=current_title) # Passé le titre
                         
                messagebox.showinfo("Suppression", f"Marqueur à Époque={int(period_to_remove)} supprimé.")
                # self._log_markers() # Note: _log_markers n'est pas fourni
            else:
                print("Clic droit ignoré : trop éloigné d'un marqueur existant.")
        
    def save_figure(self):
        """
        Sauvegarde simple de ce qui est actuellement affiché sur le Canvas.
        """
        filename = asksaveasfilename(
            title="Sauvegarder le graphique",
            defaultextension=".png",
            filetypes=[("Image PNG", "*.png"), ("PDF", "*.pdf")]
        )
        if filename:
            try:
                self.fig.savefig(filename, dpi=150, bbox_inches='tight')
                messagebox.showinfo("Succès", f"Graphique sauvegardé sous :\n{filename}")
            except Exception as e:
                messagebox.showerror("Erreur", f"Échec de la sauvegarde : {e}")

# --- Bloc de test Standalone ---
if __name__ == "__main__":
    root = tk.Tk()
    
    # Création du Viewer
    app = TTVViewer(root)
    
    # Création du contrôle Frame et des boutons externes (pour le test standalone)
    control_frame = ttk.Frame(root)
    control_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=10)
    
    ttk.Button(control_frame, text="📂 Charger Fichier TTV", command=app.load_csv).pack(side=tk.LEFT, padx=5)
    ttk.Button(control_frame, text="💾 Sauvegarder Figure", command=app.save_figure).pack(side=tk.LEFT, padx=5)
    
    root.mainloop()