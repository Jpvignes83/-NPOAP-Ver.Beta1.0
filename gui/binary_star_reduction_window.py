# gui/binary_star_reduction_window.py
"""
Fenêtre dédiée pour les techniques de réduction d'images d'étoiles binaires
inspirées de REDUC (http://www.astrosurf.com/hfosaf/reduc/tutoriel.htm)
"""
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import logging
from pathlib import Path
import threading
import json

logger = logging.getLogger(__name__)

# Import du module de réduction d'images (techniques REDUC)
try:
    from core.binary_star_reduction import BinaryStarReduction
    REDUCTION_AVAILABLE = True
except ImportError as e:
    REDUCTION_AVAILABLE = False
    logger.warning(f"Module de réduction non disponible: {e}")


class BinaryStarReductionWindow(tk.Toplevel):
    """
    Fenêtre pour le traitement d'images d'étoiles binaires avec techniques REDUC
    """
    
    def __init__(self, parent, base_dir=None):
        super().__init__(parent)
        
        self.title("Réduction d'Images - Techniques REDUC")
        self.geometry("1200x900")
        
        if base_dir is None:
            self.base_dir = Path.home()
        else:
            self.base_dir = Path(base_dir)
        
        # Processeur de réduction
        if REDUCTION_AVAILABLE:
            self.reducer = BinaryStarReduction()
        else:
            self.reducer = None
            messagebox.showerror(
                "Erreur",
                "Module de réduction non disponible!\n\n"
                "Le module core/binary_star_reduction.py doit être disponible."
            )
        
        # Liste de travail (images retenues)
        self.work_list = []  # Liste de tuples (path, metrics)
        self.work_list_file = None  # Chemin du fichier de liste sauvegardé
        
        self.create_widgets()
        
        # Barre de progression globale
        self.progress = ttk.Progressbar(self, mode="indeterminate")
        self.progress.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=5)
    
    def create_widgets(self):
        """Crée l'interface utilisateur"""
        
        # En-tête
        header_frame = ttk.Frame(self, padding=10)
        header_frame.pack(fill="x")
        
        title_label = ttk.Label(
            header_frame,
            text="Techniques de Réduction REDUC",
            font=("Helvetica", 14, "bold")
        )
        title_label.pack()
        
        info_label = ttk.Label(
            header_frame,
            text="Basé sur REDUC (http://www.astrosurf.com/hfosaf/reduc/tutoriel.htm)",
            font=("Helvetica", 8),
            foreground="gray"
        )
        info_label.pack()
        
        # Notebook pour les différentes fonctionnalités
        notebook = ttk.Notebook(self)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Onglet 1: Lucky Imaging (BestOf + ELI fusionnés)
        lucky_frame = ttk.Frame(notebook, padding=10)
        notebook.add(lucky_frame, text="Lucky Imaging (BestOf + ELI)")
        self.create_lucky_imaging_tab(lucky_frame)
        
        # Onglet 2: Mesure de séparation
        measure_frame = ttk.Frame(notebook, padding=10)
        notebook.add(measure_frame, text="Mesure de séparation")
        self.create_measure_tab(measure_frame)
    
    def create_lucky_imaging_tab(self, parent):
        """Crée l'onglet Lucky Imaging (BestOf + ELI fusionnés)"""
        
        # Chargement d'images
        load_frame = ttk.LabelFrame(parent, text="1. Chargement des images", padding=10)
        load_frame.pack(fill="x", pady=5)
        
        ttk.Label(load_frame, text="Dossier d'images:").pack(anchor="w", pady=2)
        
        self.lucky_dir_var = tk.StringVar()
        dir_frame = ttk.Frame(load_frame)
        dir_frame.pack(fill="x", pady=5)
        ttk.Entry(dir_frame, textvariable=self.lucky_dir_var, width=50).pack(side="left", fill="x", expand=True)
        ttk.Button(dir_frame, text="📁 Parcourir", command=lambda: self.browse_folder(self.lucky_dir_var)).pack(side="left", padx=(5, 0))
        
        # Section 1: Analyse et tri (BestOf)
        analysis_frame = ttk.LabelFrame(parent, text="2. Analyse et tri des images (BestOf)", padding=10)
        analysis_frame.pack(fill="x", pady=5)
        
        ttk.Label(analysis_frame, text="Pourcentage d'images à analyser/conserver:").pack(anchor="w", pady=2)
        
        self.lucky_percent_var = tk.StringVar(value="10")
        scale_frame = ttk.Frame(analysis_frame)
        scale_frame.pack(fill="x", pady=5)
        scale = ttk.Scale(scale_frame, from_=5, to=50, orient=tk.HORIZONTAL,
                         variable=self.lucky_percent_var, length=300,
                         command=lambda v: self.lucky_percent_label.config(text=f"{float(v):.0f}%"))
        scale.pack(side="left", fill="x", expand=True)
        self.lucky_percent_label = ttk.Label(scale_frame, text="10%", width=6)
        self.lucky_percent_label.pack(side="left", padx=(10, 0))
        
        ttk.Label(analysis_frame, 
                 text="Les images seront triées par qualité (FWHM, SNR, contraste) et seules les meilleures seront utilisées pour le stacking.",
                 font=("Helvetica", 8),
                 foreground="gray",
                 wraplength=600).pack(anchor="w", pady=(5, 0))
        
        # Bouton d'analyse dans la section Analyse
        analysis_button_frame = ttk.Frame(analysis_frame)
        analysis_button_frame.pack(fill="x", pady=(10, 0))
        
        ttk.Button(
            analysis_button_frame,
            text="🔍 Analyser les images (crée la liste de travail)",
            command=self.run_analysis_only
        ).pack(pady=5)
        
        # Section 3: Résultats de l'analyse et liste de travail
        result_frame = ttk.LabelFrame(parent, text="3. Résultats de l'analyse et liste de travail", padding=10)
        result_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        
        # Zone de texte avec scrollbar pour les résultats
        text_frame = ttk.Frame(result_frame)
        text_frame.pack(fill=tk.BOTH, expand=True)
        
        scrollbar = ttk.Scrollbar(text_frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.lucky_result_text = tk.Text(text_frame, height=10, yscrollcommand=scrollbar.set)
        self.lucky_result_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.lucky_result_text.yview)
        
        # Informations sur la liste de travail
        work_list_frame = ttk.Frame(result_frame)
        work_list_frame.pack(fill="x", pady=(5, 0))
        
        self.work_list_info_label = ttk.Label(
            work_list_frame,
            text="Aucune liste de travail active",
            font=("Helvetica", 9),
            foreground="gray"
        )
        self.work_list_info_label.pack(side="left")
        
        # Boutons de gestion de la liste
        list_buttons_frame = ttk.Frame(result_frame)
        list_buttons_frame.pack(fill="x", pady=5)
        
        ttk.Button(
            list_buttons_frame,
            text="💾 Sauvegarder la liste",
            command=self.save_work_list,
            state="disabled"
        ).pack(side="left", padx=2)
        self.save_list_button = list_buttons_frame.winfo_children()[-1]
        
        ttk.Button(
            list_buttons_frame,
            text="📂 Charger une liste",
            command=self.load_work_list
        ).pack(side="left", padx=2)
        
        ttk.Button(
            list_buttons_frame,
            text="📋 Voir la liste",
            command=self.view_work_list
        ).pack(side="left", padx=2)
        
        # Section 4: Stacking (ELI)
        stacking_frame = ttk.LabelFrame(parent, text="4. Stacking (ELI - Easy Lucky Imaging)", padding=10)
        stacking_frame.pack(fill="x", pady=5)
        
        ttk.Label(stacking_frame, text="Position de référence pour l'alignement sub-pixel (x, y):").pack(anchor="w", pady=2)
        ref_pos_frame = ttk.Frame(stacking_frame)
        ref_pos_frame.pack(fill="x", pady=5)
        self.ref_x_var = tk.StringVar(value="")
        self.ref_y_var = tk.StringVar(value="")
        ttk.Entry(ref_pos_frame, textvariable=self.ref_x_var, width=15).pack(side="left", padx=(0, 5))
        ttk.Label(ref_pos_frame, text=", ").pack(side="left")
        ttk.Entry(ref_pos_frame, textvariable=self.ref_y_var, width=15).pack(side="left")
        ttk.Label(ref_pos_frame, text="(en pixels - position d'une étoile brillante)", 
                 font=("Helvetica", 8), foreground="gray").pack(side="left", padx=(10, 0))
        
        ttk.Label(stacking_frame, text="Méthode de stacking:").pack(anchor="w", pady=(10, 2))
        self.lucky_method_var = tk.StringVar(value="median")
        method_frame = ttk.Frame(stacking_frame)
        method_frame.pack(fill="x", pady=5)
        for method, label in [("median", "Médiane (recommandé)"), ("mean", "Moyenne"), ("sigma_clip", "Sigma-clip")]:
            ttk.Radiobutton(method_frame, text=label, variable=self.lucky_method_var, value=method).pack(side="left", padx=(0, 20))
        
        # Option de sauvegarde automatique
        save_options_frame = ttk.Frame(stacking_frame)
        save_options_frame.pack(fill="x", pady=(10, 0))
        
        self.auto_save_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            save_options_frame,
            text="Enregistrement automatique dans le dossier source",
            variable=self.auto_save_var
        ).pack(side="left")
        
        ttk.Label(
            save_options_frame,
            text="(sinon, demande de choisir l'emplacement)",
            font=("Helvetica", 8),
            foreground="gray"
        ).pack(side="left", padx=(10, 0))
        
        # Bouton de stacking dans la section ELI
        stacking_button_frame = ttk.Frame(stacking_frame)
        stacking_button_frame.pack(fill="x", pady=(15, 5))
        
        ttk.Button(
            stacking_button_frame,
            text="✨ Créer image empilée (utilise la liste de travail)",
            command=self.run_stacking_from_work_list
        ).pack(pady=5)
        
        ttk.Label(
            stacking_button_frame,
            text="⚠️ Une liste de travail doit être créée ou chargée avant de lancer le stacking",
            font=("Helvetica", 8),
            foreground="orange"
        ).pack()
    
    def create_measure_tab(self, parent):
        """Crée l'onglet de mesure de séparation"""
        
        # Chargement d'image
        load_frame = ttk.LabelFrame(parent, text="1. Chargement de l'image", padding=10)
        load_frame.pack(fill="x", pady=5)
        
        self.measure_image_var = tk.StringVar()
        img_frame = ttk.Frame(load_frame)
        img_frame.pack(fill="x", pady=5)
        ttk.Entry(img_frame, textvariable=self.measure_image_var, width=50).pack(side="left", fill="x", expand=True)
        ttk.Button(img_frame, text="📁 Parcourir", command=self.browse_image_file).pack(side="left", padx=(5, 0))
        
        # Positions des étoiles
        positions_frame = ttk.LabelFrame(parent, text="2. Positions approximatives des étoiles", padding=10)
        positions_frame.pack(fill="x", pady=5)
        
        ttk.Label(positions_frame, text="Étoile 1 (x, y):").pack(anchor="w", pady=2)
        star1_frame = ttk.Frame(positions_frame)
        star1_frame.pack(fill="x", pady=5)
        self.star1_x_var = tk.StringVar(value="")
        self.star1_y_var = tk.StringVar(value="")
        ttk.Entry(star1_frame, textvariable=self.star1_x_var, width=15).pack(side="left", padx=(0, 5))
        ttk.Label(star1_frame, text=", ").pack(side="left")
        ttk.Entry(star1_frame, textvariable=self.star1_y_var, width=15).pack(side="left")
        ttk.Label(star1_frame, text="(en pixels)", font=("Helvetica", 8), foreground="gray").pack(side="left", padx=(10, 0))
        
        ttk.Label(positions_frame, text="Étoile 2 (x, y):").pack(anchor="w", pady=(10, 2))
        star2_frame = ttk.Frame(positions_frame)
        star2_frame.pack(fill="x", pady=5)
        self.star2_x_var = tk.StringVar(value="")
        self.star2_y_var = tk.StringVar(value="")
        ttk.Entry(star2_frame, textvariable=self.star2_x_var, width=15).pack(side="left", padx=(0, 5))
        ttk.Label(star2_frame, text=", ").pack(side="left")
        ttk.Entry(star2_frame, textvariable=self.star2_y_var, width=15).pack(side="left")
        ttk.Label(star2_frame, text="(en pixels)", font=("Helvetica", 8), foreground="gray").pack(side="left", padx=(10, 0))
        
        # Échelle pixel
        ttk.Label(positions_frame, text="Échelle pixel (arcsec/pixel):").pack(anchor="w", pady=(10, 2))
        self.pixel_scale_var = tk.StringVar(value="1.0")
        ttk.Entry(positions_frame, textvariable=self.pixel_scale_var, width=15).pack(anchor="w", pady=5)
        
        # Résultats
        result_frame = ttk.LabelFrame(parent, text="3. Résultats", padding=10)
        result_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        
        self.measure_result_text = tk.Text(result_frame, height=10)
        self.measure_result_text.pack(fill=tk.BOTH, expand=True)
        
        # Bouton d'exécution
        button_frame = ttk.Frame(parent)
        button_frame.pack(fill="x", pady=10)
        
        ttk.Button(
            button_frame,
            text="📏 Mesurer séparation",
            command=self.measure_separation
        ).pack(pady=5)
    
    def browse_folder(self, var):
        """Ouvre un dialogue pour sélectionner un dossier"""
        folder = filedialog.askdirectory(initialdir=self.base_dir)
        if folder:
            var.set(folder)
    
    def browse_image_file(self):
        """Ouvre un dialogue pour sélectionner un fichier image"""
        path = filedialog.askopenfilename(
            initialdir=self.base_dir,
            filetypes=[("FITS", "*.fits"), ("Tous les fichiers", "*.*")]
        )
        if path:
            self.measure_image_var.set(path)
    
    def run_analysis_only(self):
        """Exécute uniquement l'analyse BestOf sans créer l'image empilée"""
        if not REDUCTION_AVAILABLE or self.reducer is None:
            messagebox.showerror("Erreur", "Module de réduction non disponible")
            return
        
        image_dir = self.lucky_dir_var.get()
        if not image_dir:
            messagebox.showwarning("Attention", "Sélectionnez d'abord un dossier d'images!")
            return
        
        folder = Path(image_dir)
        # Rechercher les fichiers FITS (insensible à la casse)
        fits_files = list(folder.glob("*.fits")) + list(folder.glob("*.FITS"))
        # Dédupliquer en utilisant le chemin absolu normalisé (pour Windows)
        # Utiliser resolve() pour normaliser les chemins et gérer la casse
        seen_paths = set()
        unique_files = []
        for f in fits_files:
            resolved = f.resolve()
            if resolved not in seen_paths:
                seen_paths.add(resolved)
                unique_files.append(f)
        image_files = sorted(unique_files, key=str)
        
        if not image_files:
            messagebox.showwarning("Attention", "Aucune image FITS trouvée dans le dossier!")
            return
        
        def analysis_task():
            try:
                self.progress.start()
                percent = float(self.lucky_percent_var.get()) / 100.0
                
                results = self.reducer.bestof_sort(image_files, top_percent=percent)
                
                self.progress.stop()
                
                # Afficher les résultats dans la zone de texte
                self.lucky_result_text.delete(1.0, tk.END)
                result_text = f"Analyse terminée!\n"
                result_text += f"{'='*60}\n\n"
                result_text += f"Total images analysées: {len(image_files)}\n"
                result_text += f"Meilleures images sélectionnées: {len(results)} ({percent*100:.0f}%)\n\n"
                result_text += f"Top {min(20, len(results))} meilleures images:\n"
                result_text += f"{'-'*60}\n"
                
                for i, (path, metrics) in enumerate(results[:20], 1):
                    score = metrics.get('score', 0.0)
                    fwhm = metrics.get('fwhm', 0.0)
                    snr = metrics.get('snr', 0.0)
                    contrast = metrics.get('contrast', 0.0)
                    n_stars = metrics.get('n_stars', 0)
                    
                    result_text += f"\n{i:3d}. {path.name}\n"
                    result_text += f"      Score: {score:.2f}\n"
                    result_text += f"      FWHM: {fwhm:.2f}\" | SNR: {snr:.2f} | Contraste: {contrast:.3f} | Étoiles: {n_stars}\n"
                
                result_text += f"\n{'='*60}\n"
                result_text += f"Pour créer l'image empilée, remplissez les paramètres de stacking ci-dessous\n"
                result_text += f"et cliquez sur 'Analyser et créer image empilée'.\n"
                
                self.lucky_result_text.insert(1.0, result_text)
                
                # Sauvegarder la liste de travail
                self.work_list = results
                self.work_list_file = None  # Réinitialiser le fichier (nouvelle analyse)
                
                # Mettre à jour l'interface
                self.update_work_list_info()
                
                logger.info(f"Analyse: {len(results)} meilleures images sur {len(image_files)}")
                
            except Exception as e:
                self.progress.stop()
                messagebox.showerror("Erreur", f"Erreur lors de l'analyse: {e}")
                logger.error(f"Erreur analyse: {e}", exc_info=True)
        
        threading.Thread(target=analysis_task, daemon=True).start()
    
    def update_work_list_info(self):
        """Met à jour l'affichage des informations sur la liste de travail"""
        if self.work_list:
            n_images = len(self.work_list)
            if self.work_list_file:
                file_name = Path(self.work_list_file).name
                self.work_list_info_label.config(
                    text=f"✅ Liste de travail active: {n_images} images (depuis {file_name})",
                    foreground="green"
                )
            else:
                self.work_list_info_label.config(
                    text=f"✅ Liste de travail active: {n_images} images (analysées)",
                    foreground="green"
                )
            self.save_list_button.config(state="normal")
        else:
            self.work_list_info_label.config(
                text="Aucune liste de travail active",
                foreground="gray"
            )
            self.save_list_button.config(state="disabled")
    
    def save_work_list(self):
        """Sauvegarde la liste de travail dans un fichier JSON"""
        if not self.work_list:
            messagebox.showwarning("Attention", "Aucune liste de travail à sauvegarder!")
            return
        
        # Demander le fichier de sauvegarde
        save_path = filedialog.asksaveasfilename(
            initialdir=self.base_dir,
            defaultextension=".json",
            filetypes=[("JSON", "*.json"), ("Tous les fichiers", "*.*")],
            title="Sauvegarder la liste de travail"
        )
        
        if not save_path:
            return
        
        try:
            # Préparer les données à sauvegarder
            data = {
                'total_images': len(self.work_list),
                'source_directory': str(self.lucky_dir_var.get()),
                'percent_used': float(self.lucky_percent_var.get()),
                'images': []
            }
            
            for path, metrics in self.work_list:
                data['images'].append({
                    'path': str(path),
                    'filename': path.name,
                    'metrics': {
                        'score': float(metrics.get('score', 0.0)),
                        'fwhm': float(metrics.get('fwhm', 0.0)),
                        'snr': float(metrics.get('snr', 0.0)),
                        'contrast': float(metrics.get('contrast', 0.0)),
                        'n_stars': int(metrics.get('n_stars', 0))
                    }
                })
            
            # Sauvegarder en JSON
            with open(save_path, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            
            self.work_list_file = save_path
            self.update_work_list_info()
            
            messagebox.showinfo("Succès", f"Liste de travail sauvegardée!\n\n{save_path}\n\n{len(self.work_list)} images")
            logger.info(f"Liste de travail sauvegardée: {save_path} ({len(self.work_list)} images)")
            
        except Exception as e:
            messagebox.showerror("Erreur", f"Erreur lors de la sauvegarde: {e}")
            logger.error(f"Erreur sauvegarde liste: {e}", exc_info=True)
    
    def load_work_list(self):
        """Charge une liste de travail depuis un fichier JSON"""
        load_path = filedialog.askopenfilename(
            initialdir=self.base_dir,
            filetypes=[("JSON", "*.json"), ("Tous les fichiers", "*.*")],
            title="Charger une liste de travail"
        )
        
        if not load_path:
            return
        
        try:
            with open(load_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            # Vérifier la structure
            if 'images' not in data:
                raise ValueError("Format de fichier invalide: 'images' manquant")
            
            # Reconstruire la liste de travail (en évitant les doublons)
            self.work_list = []
            seen_paths = set()
            for item in data['images']:
                path = Path(item['path'])
                # Vérifier que le fichier existe
                if not path.exists():
                    logger.warning(f"Image non trouvée: {path}")
                    continue
                
                # Dédupliquer en utilisant le chemin résolu
                resolved = path.resolve()
                if resolved in seen_paths:
                    logger.debug(f"Fichier dupliqué ignoré: {path}")
                    continue
                seen_paths.add(resolved)
                
                metrics = item.get('metrics', {})
                self.work_list.append((path, metrics))
            
            if not self.work_list:
                messagebox.showwarning("Attention", "Aucune image valide trouvée dans la liste!")
                return
            
            self.work_list_file = load_path
            
            # Restaurer les paramètres si disponibles
            if 'source_directory' in data:
                self.lucky_dir_var.set(data['source_directory'])
            if 'percent_used' in data:
                self.lucky_percent_var.set(str(data['percent_used']))
                self.lucky_percent_label.config(text=f"{data['percent_used']:.0f}%")
            
            # Afficher les informations
            self.display_work_list_info(data)
            
            # Mettre à jour l'interface
            self.update_work_list_info()
            
            messagebox.showinfo("Succès", f"Liste de travail chargée!\n\n{len(self.work_list)} images")
            logger.info(f"Liste de travail chargée: {load_path} ({len(self.work_list)} images)")
            
        except Exception as e:
            messagebox.showerror("Erreur", f"Erreur lors du chargement: {e}")
            logger.error(f"Erreur chargement liste: {e}", exc_info=True)
    
    def display_work_list_info(self, data=None):
        """Affiche les informations de la liste de travail dans la zone de texte"""
        if not self.work_list:
            return
        
        if data is None:
            # Reconstruire les infos depuis la liste
            result_text = f"Liste de travail chargée\n"
            result_text += f"{'='*60}\n\n"
            result_text += f"Nombre d'images: {len(self.work_list)}\n"
            if self.work_list_file:
                result_text += f"Source: {Path(self.work_list_file).name}\n"
            result_text += f"\n{'='*60}\n\n"
            result_text += f"Top {min(20, len(self.work_list))} images:\n"
            result_text += f"{'-'*60}\n"
            
            for i, (path, metrics) in enumerate(self.work_list[:20], 1):
                score = metrics.get('score', 0.0)
                fwhm = metrics.get('fwhm', 0.0)
                snr = metrics.get('snr', 0.0)
                contrast = metrics.get('contrast', 0.0)
                
                result_text += f"\n{i:3d}. {path.name}\n"
                result_text += f"      Score: {score:.2f} | FWHM: {fwhm:.2f}\" | SNR: {snr:.2f}\n"
        else:
            # Utiliser les données du fichier
            result_text = f"Liste de travail chargée\n"
            result_text += f"{'='*60}\n\n"
            result_text += f"Nombre d'images: {len(self.work_list)}\n"
            if 'source_directory' in data:
                result_text += f"Dossier source: {data['source_directory']}\n"
            if 'percent_used' in data:
                result_text += f"Pourcentage utilisé: {data['percent_used']:.0f}%\n"
            result_text += f"Fichier: {Path(self.work_list_file).name}\n"
            result_text += f"\n{'='*60}\n\n"
            result_text += f"Top {min(20, len(self.work_list))} images:\n"
            result_text += f"{'-'*60}\n"
            
            for i, item in enumerate(data['images'][:20], 1):
                filename = item.get('filename', item['path'])
                metrics = item.get('metrics', {})
                score = metrics.get('score', 0.0)
                fwhm = metrics.get('fwhm', 0.0)
                snr = metrics.get('snr', 0.0)
                
                result_text += f"\n{i:3d}. {filename}\n"
                result_text += f"      Score: {score:.2f} | FWHM: {fwhm:.2f}\" | SNR: {snr:.2f}\n"
        
        result_text += f"\n{'='*60}\n"
        result_text += f"Vous pouvez maintenant créer l'image empilée en utilisant cette liste.\n"
        
        self.lucky_result_text.delete(1.0, tk.END)
        self.lucky_result_text.insert(1.0, result_text)
    
    def view_work_list(self):
        """Affiche la liste de travail dans la zone de texte"""
        if not self.work_list:
            messagebox.showwarning("Attention", "Aucune liste de travail active!")
            return
        
        self.display_work_list_info()
    
    def run_stacking_from_work_list(self):
        """Exécute le stacking en utilisant la liste de travail"""
        if not REDUCTION_AVAILABLE or self.reducer is None:
            messagebox.showerror("Erreur", "Module de réduction non disponible")
            return
        
        if not self.work_list:
            messagebox.showwarning(
                "Attention",
                "Aucune liste de travail active!\n\n"
                "Veuillez d'abord analyser les images ou charger une liste sauvegardée."
            )
            return
        
        ref_x_str = self.ref_x_var.get()
        ref_y_str = self.ref_y_var.get()
        if not ref_x_str or not ref_y_str:
            messagebox.showwarning("Attention", "Entrez la position de référence (x, y) pour l'alignement!")
            return
        
        try:
            ref_x = float(ref_x_str)
            ref_y = float(ref_y_str)
        except ValueError:
            messagebox.showerror("Erreur", "Les positions doivent être des nombres!")
            return
        
        # Déterminer le chemin de sauvegarde
        if self.auto_save_var.get():
            # Sauvegarde automatique dans le dossier source
            if self.work_list:
                source_dir = self.work_list[0][0].parent
            elif self.lucky_dir_var.get():
                source_dir = Path(self.lucky_dir_var.get())
            else:
                source_dir = self.base_dir
            
            # Générer un nom de fichier avec timestamp
            from datetime import datetime
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            method_name = self.lucky_method_var.get()
            n_images = len(self.work_list) if self.work_list else 0
            output_filename = f"ELI_{method_name}_{n_images}imgs_{timestamp}.fits"
            output_path = source_dir / output_filename
        else:
            # Demander à l'utilisateur de choisir l'emplacement
            default_dir = self.base_dir
            if self.work_list:
                default_dir = self.work_list[0][0].parent
            elif self.lucky_dir_var.get():
                default_dir = Path(self.lucky_dir_var.get())
            
            # Suggérer un nom de fichier
            from datetime import datetime
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            method_name = self.lucky_method_var.get()
            n_images = len(self.work_list) if self.work_list else 0
            default_filename = f"ELI_{method_name}_{n_images}imgs_{timestamp}.fits"
            
            output_path = filedialog.asksaveasfilename(
                initialdir=default_dir,
                initialfile=default_filename,
                defaultextension=".fits",
                filetypes=[("FITS", "*.fits"), ("Tous les fichiers", "*.*")]
            )
            
            if not output_path:
                return
            
            output_path = Path(output_path)
        
        def stacking_task():
            try:
                self.progress.start()
                
                method = self.lucky_method_var.get()
                best_image_files = [path for path, _ in self.work_list]
                
                # Afficher le statut
                self.lucky_result_text.delete(1.0, tk.END)
                result_text = f"Stacking en cours...\n"
                result_text += f"{'='*60}\n\n"
                result_text += f"Utilisation de la liste de travail: {len(best_image_files)} images\n"
                if self.work_list_file:
                    result_text += f"Source: {Path(self.work_list_file).name}\n"
                result_text += f"Méthode: {method}\n"
                result_text += f"Position référence: ({ref_x}, {ref_y})\n"
                self.lucky_result_text.insert(1.0, result_text)
                self.lucky_result_text.update()
                
                # Stacking (ELI)
                success = self.reducer.eli_lucky_imaging(
                    best_image_files,
                    output_path,
                    (ref_x, ref_y),
                    top_percent=1.0,  # Toutes les images sont déjà sélectionnées
                    method=method
                )
                
                self.progress.stop()
                
                if success:
                    result_text += f"\n{'='*60}\n"
                    result_text += f"✅ Image empilée créée avec succès!\n"
                    result_text += f"Fichier: {output_path.name}\n"
                    result_text += f"Chemin complet: {output_path}\n"
                    
                    self.lucky_result_text.delete(1.0, tk.END)
                    self.lucky_result_text.insert(1.0, result_text)
                    
                    messagebox.showinfo("Succès", f"Image empilée créée avec succès!\n\n{output_path}")
                    logger.info(f"Stacking terminé: {output_path}")
                else:
                    messagebox.showerror("Erreur", "Échec de la création de l'image empilée")
                    
            except Exception as e:
                self.progress.stop()
                messagebox.showerror("Erreur", f"Erreur lors du stacking: {e}")
                logger.error(f"Erreur stacking: {e}", exc_info=True)
        
        threading.Thread(target=stacking_task, daemon=True).start()
    
    def measure_separation(self):
        """Mesure la séparation d'un système binaire"""
        if not REDUCTION_AVAILABLE or self.reducer is None:
            messagebox.showerror("Erreur", "Module de réduction non disponible")
            return
        
        image_path = self.measure_image_var.get()
        if not image_path:
            messagebox.showwarning("Attention", "Sélectionnez d'abord une image!")
            return
        
        star1_x_str = self.star1_x_var.get()
        star1_y_str = self.star1_y_var.get()
        star2_x_str = self.star2_x_var.get()
        star2_y_str = self.star2_y_var.get()
        pixel_scale_str = self.pixel_scale_var.get()
        
        if not all([star1_x_str, star1_y_str, star2_x_str, star2_y_str]):
            messagebox.showwarning("Attention", "Entrez les positions des deux étoiles!")
            return
        
        try:
            star1_pos = (float(star1_x_str), float(star1_y_str))
            star2_pos = (float(star2_x_str), float(star2_y_str))
            pixel_scale = float(pixel_scale_str)
        except ValueError:
            messagebox.showerror("Erreur", "Les valeurs doivent être des nombres!")
            return
        
        def measure_task():
            try:
                self.progress.start()
                
                result = self.reducer.measure_binary_separation(
                    Path(image_path),
                    star1_pos,
                    star2_pos,
                    pixel_scale
                )
                
                self.progress.stop()
                
                # Afficher les résultats
                self.measure_result_text.delete(1.0, tk.END)
                result_text = "Résultats de la mesure:\n"
                result_text += f"{'='*60}\n\n"
                result_text += f"Étoile 1 (centroïde précis):\n"
                result_text += f"  x = {result['x1']:.3f} pixels\n"
                result_text += f"  y = {result['y1']:.3f} pixels\n\n"
                result_text += f"Étoile 2 (centroïde précis):\n"
                result_text += f"  x = {result['x2']:.3f} pixels\n"
                result_text += f"  y = {result['y2']:.3f} pixels\n\n"
                result_text += f"{'-'*60}\n\n"
                result_text += f"Séparation: {result['separation_pix']:.3f} pixels\n"
                result_text += f"Séparation: {result['separation_arcsec']:.3f} arcsec\n"
                result_text += f"Angle de position: {result['position_angle']:.2f}°\n\n"
                if result.get('sigma_theta', 0) > 0 or result.get('sigma_rho', 0) > 0:
                    result_text += f"Incertitudes: σθ = {result.get('sigma_theta', 0):.4f}°, σρ = {result.get('sigma_rho', 0):.4f}\"\n\n"
                result_text += f"(Angle mesuré depuis le nord vers l'est)\n"
                
                self.measure_result_text.insert(1.0, result_text)
                logger.info(f"Mesure séparation: {result['separation_arcsec']:.3f}\", PA: {result['position_angle']:.2f}°")
                
            except Exception as e:
                self.progress.stop()
                messagebox.showerror("Erreur", f"Erreur lors de la mesure: {e}")
                logger.error(f"Erreur mesure séparation: {e}", exc_info=True)
        
        threading.Thread(target=measure_task, daemon=True).start()

