# gui/binary_stars_viewer.py
"""
Fenêtre de visualisation 3D et animation des étoiles binaires
Inspiré du simulateur NAAP : https://astro.unl.edu/naap/ebs/animations/ebs.html
"""
import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle
import logging

logger = logging.getLogger(__name__)

class BinaryStarsViewer(tk.Toplevel):
    """
    Fenêtre de visualisation interactive des étoiles binaires
    Affiche l'animation du système binaire en orbite avec vue 3D/2D
    """
    
    def __init__(self, parent, bundle=None):
        super().__init__(parent)
        self.title("Visualisation 3D - Système Binaire")
        self.geometry("1600x1000")
        
        self.bundle = bundle
        self.animation = None
        self.is_playing = False
        
        # Paramètres par défaut (seront remplacés par le bundle si disponible)
        self.period = 1.0  # jours
        self.inclination = 90.0  # degrés
        self.r1 = 1.0  # rayon étoile 1 (R☉)
        self.r2 = 0.8  # rayon étoile 2 (R☉)
        self.sma = 10.0  # demi-grand axe (R☉)
        self.q = 0.8  # ratio de masse (M2/M1)
        
        # Charger les paramètres du bundle si disponible
        if bundle is not None:
            self.load_parameters_from_bundle()
        
        # Temps d'animation
        self.time = 0.0
        self.time_step = 0.01  # fraction de période
        
        self.create_widgets()
        
        # Démarrer l'animation
        self.start_animation()
    
    def load_parameters_from_bundle(self):
        """Charge les paramètres depuis le bundle PHOEBE2"""
        try:
            # Vérifier si PHOEBE2 est disponible (importé dans binary_stars_tab)
            try:
                import phoebe
                phoebe_available = True
            except:
                phoebe_available = False
            
            if phoebe_available and self.bundle is not None:
                self.period = float(self.bundle.get_value('period@binary', unit='d'))
                self.inclination = float(self.bundle.get_value('incl@binary', unit='deg'))
                
                # Rayons
                try:
                    self.r1 = float(self.bundle.get_value('requiv@primary', unit='solRad'))
                    self.r2 = float(self.bundle.get_value('requiv@secondary', unit='solRad'))
                except:
                    pass
                
                # Demi-grand axe
                try:
                    self.sma = float(self.bundle.get_value('sma@binary', unit='solRad'))
                except:
                    pass
                
                # Ratio de masse
                try:
                    m1 = float(self.bundle.get_value('mass@primary', unit='solMass'))
                    m2 = float(self.bundle.get_value('mass@secondary', unit='solMass'))
                    self.q = m2 / m1 if m1 > 0 else 0.8
                except:
                    pass
        except Exception as e:
            logger.warning(f"Impossible de charger les paramètres du bundle: {e}")
    
    def create_widgets(self):
        """Crée l'interface de la fenêtre"""
        
        # Frame principal
        main_frame = ttk.Frame(self, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Contrôles en haut
        control_frame = ttk.LabelFrame(main_frame, text="Contrôles", padding=10)
        control_frame.pack(fill="x", pady=(0, 10))
        
        # Boutons de contrôle
        btn_frame = ttk.Frame(control_frame)
        btn_frame.pack(fill="x")
        
        self.play_pause_btn = ttk.Button(
            btn_frame,
            text="⏸ Pause",
            command=self.toggle_animation
        )
        self.play_pause_btn.pack(side="left", padx=5)
        
        ttk.Button(
            btn_frame,
            text="⏮ Reset",
            command=self.reset_animation
        ).pack(side="left", padx=5)
        
        ttk.Button(
            btn_frame,
            text="🔄 Actualiser",
            command=self.update_visualization
        ).pack(side="left", padx=5)
        
        # Séparateur
        ttk.Separator(btn_frame, orient=tk.VERTICAL).pack(side="left", fill="y", padx=10)
        
        # Vitesse d'animation
        ttk.Label(btn_frame, text="Vitesse:").pack(side="left", padx=5)
        self.speed_var = tk.DoubleVar(value=1.0)
        speed_scale = ttk.Scale(
            btn_frame,
            from_=0.1,
            to=5.0,
            variable=self.speed_var,
            orient=tk.HORIZONTAL,
            length=150
        )
        speed_scale.pack(side="left", padx=5)
        self.speed_label = ttk.Label(btn_frame, text="1.0x")
        self.speed_label.pack(side="left", padx=5)
        speed_scale.configure(command=self.update_speed_label)
        
        # Zone de visualisation
        viz_frame = ttk.Frame(main_frame)
        viz_frame.pack(fill=tk.BOTH, expand=True)
        
        # Créer la figure matplotlib
        self.fig = plt.figure(figsize=(12, 8))
        self.fig.suptitle("Visualisation du Système Binaire", fontsize=14, fontweight='bold')
        
        # 3 sous-graphiques : Vue de face, Vue de côté, Vue de dessus
        self.ax_face = self.fig.add_subplot(2, 2, 1)
        self.ax_side = self.fig.add_subplot(2, 2, 2)
        self.ax_top = self.fig.add_subplot(2, 2, 3)
        self.ax_info = self.fig.add_subplot(2, 2, 4)
        self.ax_info.axis('off')
        
        # Canvas
        self.canvas = FigureCanvasTkAgg(self.fig, master=viz_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Toolbar
        toolbar = NavigationToolbar2Tk(self.canvas, viz_frame)
        toolbar.update()
        
        # Configuration des axes
        self.setup_axes()
    
    def setup_axes(self):
        """Configure les axes des graphiques"""
        self.ax_face.set_title("Vue de Face (X-Y)")
        self.ax_face.set_xlabel("X (R☉)")
        self.ax_face.set_ylabel("Y (R☉)")
        self.ax_face.set_aspect('equal')
        self.ax_face.grid(True, alpha=0.3)
        
        self.ax_side.set_title("Vue de Côté (X-Z)")
        self.ax_side.set_xlabel("X (R☉)")
        self.ax_side.set_ylabel("Z (R☉)")
        self.ax_side.set_aspect('equal')
        self.ax_side.grid(True, alpha=0.3)
        
        self.ax_top.set_title("Vue de Dessus (Y-Z)")
        self.ax_top.set_xlabel("Y (R☉)")
        self.ax_top.set_ylabel("Z (R☉)")
        self.ax_top.set_aspect('equal')
        self.ax_top.grid(True, alpha=0.3)
    
    def calculate_positions(self, phase):
        """
        Calcule les positions des deux étoiles à une phase donnée
        phase: fraction de la période (0.0 à 1.0)
        """
        # Angle orbital
        theta = 2 * np.pi * phase
        
        # Distance du centre de masse
        # Pour un système binaire : r1 = a * q/(1+q), r2 = a * 1/(1+q)
        a1 = self.sma * self.q / (1 + self.q)
        a2 = self.sma / (1 + self.q)
        
        # Positions dans le plan orbital (avant inclinaison)
        x1_orb = -a1 * np.cos(theta)
        y1_orb = -a1 * np.sin(theta)
        z1_orb = 0
        
        x2_orb = a2 * np.cos(theta)
        y2_orb = a2 * np.sin(theta)
        z2_orb = 0
        
        # Application de l'inclinaison (rotation autour de l'axe X)
        incl_rad = np.radians(self.inclination)
        cos_i = np.cos(incl_rad)
        sin_i = np.sin(incl_rad)
        
        # Position étoile 1
        x1 = x1_orb
        y1 = y1_orb * cos_i
        z1 = y1_orb * sin_i
        
        # Position étoile 2
        x2 = x2_orb
        y2 = y2_orb * cos_i
        z2 = y2_orb * sin_i
        
        return (x1, y1, z1), (x2, y2, z2)
    
    def update_visualization(self):
        """Met à jour la visualisation"""
        phase = (self.time % 1.0)
        pos1, pos2 = self.calculate_positions(phase)
        
        # Effacer les axes
        self.ax_face.clear()
        self.ax_side.clear()
        self.ax_top.clear()
        self.ax_info.clear()
        self.ax_info.axis('off')
        
        # Configuration des axes
        self.setup_axes()
        
        # Limites des graphiques (basées sur le demi-grand axe)
        limit = self.sma * 1.5
        self.ax_face.set_xlim(-limit, limit)
        self.ax_face.set_ylim(-limit, limit)
        self.ax_side.set_xlim(-limit, limit)
        self.ax_side.set_ylim(-limit, limit)
        self.ax_top.set_xlim(-limit, limit)
        self.ax_top.set_ylim(-limit, limit)
        
        # Dessiner les étoiles
        # Vue de face (X-Y)
        circle1_face = Circle((pos1[0], pos1[1]), self.r1, color='red', alpha=0.7, label='Étoile 1')
        circle2_face = Circle((pos2[0], pos2[1]), self.r2, color='blue', alpha=0.7, label='Étoile 2')
        self.ax_face.add_patch(circle1_face)
        self.ax_face.add_patch(circle2_face)
        
        # Vue de côté (X-Z)
        circle1_side = Circle((pos1[0], pos1[2]), self.r1, color='red', alpha=0.7)
        circle2_side = Circle((pos2[0], pos2[2]), self.r2, color='blue', alpha=0.7)
        self.ax_side.add_patch(circle1_side)
        self.ax_side.add_patch(circle2_side)
        
        # Vue de dessus (Y-Z)
        circle1_top = Circle((pos1[1], pos1[2]), self.r1, color='red', alpha=0.7)
        circle2_top = Circle((pos2[1], pos2[2]), self.r2, color='blue', alpha=0.7)
        self.ax_top.add_patch(circle1_top)
        self.ax_top.add_patch(circle2_top)
        
        # Légende
        self.ax_face.legend(loc='upper right')
        
        # Informations texte
        info_text = f"""
Paramètres du Système:
━━━━━━━━━━━━━━━━━━━━━━━━━━
Période: {self.period:.3f} jours
Inclinaison: {self.inclination:.1f}°
Demi-grand axe: {self.sma:.2f} R☉

Rayons:
  Étoile 1: {self.r1:.2f} R☉
  Étoile 2: {self.r2:.2f} R☉
  Ratio: {self.r2/self.r1:.2f}

Ratio de masse (q): {self.q:.3f}

État Actuel:
━━━━━━━━━━━━━━━━━━━━━━━━━━
Phase: {phase:.3f}
Temps: {self.time:.3f} périodes
        """
        self.ax_info.text(0.1, 0.5, info_text, fontsize=10, 
                         verticalalignment='center', fontfamily='monospace')
        
        # Trajectoires (optionnel, afficher quelques orbites précédentes)
        self.draw_orbits(phase)
        
        self.canvas.draw()
    
    def draw_orbits(self, current_phase):
        """Dessine les trajectoires orbitales"""
        n_points = 100
        phases = np.linspace(0, 1, n_points)
        
        x1_orbit = []
        y1_orbit = []
        z1_orbit = []
        x2_orbit = []
        y2_orbit = []
        z2_orbit = []
        
        for p in phases:
            pos1, pos2 = self.calculate_positions(p)
            x1_orbit.append(pos1[0])
            y1_orbit.append(pos1[1])
            z1_orbit.append(pos1[2])
            x2_orbit.append(pos2[0])
            y2_orbit.append(pos2[1])
            z2_orbit.append(pos2[2])
        
        # Vue de face
        self.ax_face.plot(x1_orbit, y1_orbit, 'r--', alpha=0.3, linewidth=1)
        self.ax_face.plot(x2_orbit, y2_orbit, 'b--', alpha=0.3, linewidth=1)
        
        # Vue de côté
        self.ax_side.plot(x1_orbit, z1_orbit, 'r--', alpha=0.3, linewidth=1)
        self.ax_side.plot(x2_orbit, z2_orbit, 'b--', alpha=0.3, linewidth=1)
        
        # Vue de dessus
        self.ax_top.plot(y1_orbit, z1_orbit, 'r--', alpha=0.3, linewidth=1)
        self.ax_top.plot(y2_orbit, z2_orbit, 'b--', alpha=0.3, linewidth=1)
        
        # Centre de masse
        self.ax_face.plot(0, 0, 'k+', markersize=10, markeredgewidth=2)
        self.ax_side.plot(0, 0, 'k+', markersize=10, markeredgewidth=2)
        self.ax_top.plot(0, 0, 'k+', markersize=10, markeredgewidth=2)
    
    def animate(self, frame):
        """Fonction d'animation"""
        if self.is_playing:
            self.time += self.time_step * self.speed_var.get()
            self.update_visualization()
        return []
    
    def start_animation(self):
        """Démarre l'animation"""
        self.is_playing = True
        self.animation = FuncAnimation(
            self.fig,
            self.animate,
            interval=50,  # ms
            blit=False,
            cache_frame_data=False
        )
        self.canvas.draw()
    
    def toggle_animation(self):
        """Met en pause/reprend l'animation"""
        self.is_playing = not self.is_playing
        if self.is_playing:
            self.play_pause_btn.config(text="⏸ Pause")
        else:
            self.play_pause_btn.config(text="▶ Play")
    
    def reset_animation(self):
        """Remet l'animation à zéro"""
        self.time = 0.0
        self.update_visualization()
    
    def update_speed_label(self, value):
        """Met à jour le label de vitesse"""
        speed = float(value)
        self.speed_label.config(text=f"{speed:.1f}x")
    
    def on_close(self):
        """Gère la fermeture de la fenêtre"""
        if self.animation:
            self.animation.event_source.stop()
        self.destroy()

# Import conditionnel pour PHOEBE_AVAILABLE
try:
    import sys
    if sys.platform == 'win32' and 'readline' not in sys.modules:
        class MockReadline:
            def add_history(self, *args, **kwargs): pass
            def read_history_file(self, *args, **kwargs): pass
            def write_history_file(self, *args, **kwargs): pass
            def set_history_length(self, *args, **kwargs): pass
            def get_history_length(self, *args, **kwargs): return 0
            def set_completer(self, *args, **kwargs): pass
            def set_completer_delims(self, *args, **kwargs): pass
            def set_completion_display_matches_hook(self, *args, **kwargs): pass
            def parse_and_bind(self, *args, **kwargs): pass
        sys.modules['readline'] = MockReadline()
    import phoebe
    PHOEBE_AVAILABLE = True
except:
    PHOEBE_AVAILABLE = False

