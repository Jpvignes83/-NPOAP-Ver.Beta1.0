# gui/main_window.py

import tkinter as tk
from tkinter import ttk

from gui.home_tab import HomeTab
from gui.planetarium_tab import PlanetariumTab
try:
    from gui.night_observation_tab import NightObservationTab
    NIGHT_OBS_AVAILABLE = True
except ImportError as e:
    NIGHT_OBS_AVAILABLE = False
    NightObservationTab = None
    print(f"Warning: NightObservationTab non disponible: {e}")
from gui.data_reduction_tab import CCDProcGUI
from gui.photometry_exoplanets_tab import PhotometryExoplanetsTab
from gui.asteroid_photometry_tab import AsteroidPhotometryTab
from gui.transient_photometry_tab import TransientPhotometryTab
from gui.data_analysis_tab import DataAnalysisTab
from gui.binary_stars_tab import BinaryStarsTab
try:
    from gui.easy_lucky_imaging_tab import EasyLuckyImagingTab
    EASY_LUCKY_AVAILABLE = True
except ImportError as e:
    EASY_LUCKY_AVAILABLE = False
    EasyLuckyImagingTab = None
    print(f"Warning: EasyLuckyImagingTab non disponible: {e}")
from gui.cluster_analysis_tab import ClusterAnalysisTab
from gui.spectroscopy_tab import SpectroscopyTab
try:
    from gui.catalogues_tab import CataloguesTab
    CATALOGUES_AVAILABLE = True
except ImportError as e:
    CATALOGUES_AVAILABLE = False
    CataloguesTab = None
    print(f"Warning: CataloguesTab non disponible: {e}")


class MainWindow:
    def __init__(self, root):
        self.root = root
        self.root.title("NPOAP")
        self.root.geometry("1600x1000")

        # Configuration du style des onglets : police plus petite pour une meilleure lisibilité
        style = ttk.Style()
        style.configure(
            "TNotebook.Tab",
            font=("Arial", 8, "bold"),
            padding=(10, 12),
        )

        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # ---------------------------
        # Initialisation des onglets
        # ---------------------------

        self.home_tab = HomeTab(self.notebook)
        if NIGHT_OBS_AVAILABLE:
            self.night_observation_tab = NightObservationTab(self.notebook)
        else:
            self.night_observation_tab = None
        
        self.data_reduction_tab = CCDProcGUI(self.notebook)
        self.photometry_exoplanets_tab = PhotometryExoplanetsTab(self.notebook, base_dir=None)
        self.asteroid_photometry_tab = AsteroidPhotometryTab(self.notebook)
        self.transient_photometry_tab = TransientPhotometryTab(self.notebook)
        self.data_analysis_tab = DataAnalysisTab(self.notebook)
        self.binary_stars_tab = BinaryStarsTab(self.notebook)
        self.planetarium_tab = PlanetariumTab(self.notebook)
        if EASY_LUCKY_AVAILABLE:
            self.easy_lucky_imaging_tab = EasyLuckyImagingTab(self.notebook)
        else:
            self.easy_lucky_imaging_tab = None
        self.cluster_analysis_tab = ClusterAnalysisTab(self.notebook)
        self.spectroscopy_tab = SpectroscopyTab(self.notebook)
        if CATALOGUES_AVAILABLE:
            self.catalogues_tab = CataloguesTab(self.notebook)
        else:
            self.catalogues_tab = None

        # ---------------------------
        # Ajout dans le Notebook
        # ---------------------------

        self.notebook.add(self.home_tab, text="🏠 Accueil")
        if self.night_observation_tab is not None:
            # Lien NightObservationTab → PlanetariumTab (C2A)
            try:
                self.night_observation_tab.set_c2a_visualizer(self.planetarium_tab.update_targets_from_night_tab)
            except Exception:
                pass
            self.notebook.add(self.night_observation_tab, text="🌙 Observation de la Nuit")
        # Onglet Planétarium immédiatement après « Observation de la Nuit »
        self.notebook.add(self.planetarium_tab, text="🪐 Planétarium (C2A)")
        # Puis réduction et le reste de la chaîne
        self.notebook.add(self.data_reduction_tab.frame, text="🛠️ Réduction de Données")
        self.notebook.add(self.photometry_exoplanets_tab, text="🔭 Photométrie Exoplanètes")
        self.notebook.add(self.asteroid_photometry_tab.frame, text="🛰️ Photométrie Astéroïdes")
        self.notebook.add(self.transient_photometry_tab, text="💥 Photométrie Transitoires")
        self.notebook.add(self.data_analysis_tab.main_frame, text="📈 Analyse des Données")
        self.notebook.add(self.binary_stars_tab, text="⭐ Étoiles Binaires")
        if self.easy_lucky_imaging_tab is not None:
            self.notebook.add(self.easy_lucky_imaging_tab, text="✨ Easy Lucky Imaging")
        self.notebook.add(self.cluster_analysis_tab, text="📊 Analyse d'amas")
        self.notebook.add(self.spectroscopy_tab, text="🔬 Spectroscopie")
        if self.catalogues_tab is not None:
            self.notebook.add(self.catalogues_tab, text="📚 Catalogues")

        # Gestion fermeture propre
        self.root.protocol("WM_DELETE_WINDOW", self.on_quit)

    def on_quit(self):
        self.root.destroy()


if __name__ == "__main__":
    root = tk.Tk()
    app = MainWindow(root)
    root.mainloop()
