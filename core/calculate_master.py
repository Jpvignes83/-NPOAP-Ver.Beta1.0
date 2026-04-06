import logging
import time
import numpy as np
from astropy.io import fits
from astropy.nddata import CCDData
from pathlib import Path
import tkinter as tk
from tkinter import ttk

# Configuration du logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

class CalculateMaster:
    def __init__(self, output_dir, progress_bar=None, progress_label=None):
        """Initialisation avec le répertoire de sortie et widgets de progression."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.progress_bar = progress_bar
        self.progress_label = progress_label

    def create_master(self, files, master_name, combine_function=np.median, normalize=False):
        """
        Crée un fichier maître (bias, dark, flat) à partir des fichiers fournis.
        
        Parameters:
          files : liste de chemins vers les fichiers FITS à combiner.
          master_name : nom du fichier maître à enregistrer.
          combine_function : fonction de combinaison (ex : np.median).
          normalize : booléen indiquant si une normalisation doit être appliquée (pour les flats).

        Returns:
          Le chemin du fichier maître créé, ou None en cas d'erreur.
        """
        if not files:
            logging.warning(f"⚠️ Aucun fichier fourni pour {master_name}.")
            return None

        try:
            num_files = len(files)
            data_list = []

            # Chargement des données FITS et conversion en float32 dès la lecture
            for i, file in enumerate(files):
                ccd = CCDData.read(file, unit='adu')
                data_list.append(ccd.data.astype(np.float32))
                self.update_progress(i + 1, num_files, f"📂 Lecture de {file}")

            # Combinaison des données sur l'axe 0
            master_data = combine_function(data_list, axis=0).astype(np.float32)

            # Normalisation pour les flats si demandé
            if normalize:
                median_value = np.median(master_data)
                if median_value == 0:
                    logging.warning("⚠️ Master Flat a une médiane de 0, impossible de normaliser.")
                else:
                    master_data /= median_value
                    logging.info("📏 Normalisation appliquée au Master Flat.")

            # Correction pour éviter la division par zéro dans le cas des Flats
            if "flat" in master_name.lower() and np.any(master_data == 0):
                logging.warning("⚠️ Master Flat contient des pixels à zéro, risque de division par zéro.")
                master_data[master_data == 0] = 1

            # Sauvegarde du fichier maître en forçant le type float32 et le byteorder
            master_path = self.output_dir / master_name
            hdu = fits.PrimaryHDU(master_data.newbyteorder('<'))
            hdu.writeto(master_path, overwrite=True)
            logging.info(f"✅ Fichier maître enregistré : {master_path}")

            self.update_progress(num_files, num_files, "🎉 Création terminée !")

            # Validation de la bonne création du fichier maître
            self.validate_master_file(master_path, master_name)
            return master_path

        except Exception as e:
            logging.error(f"❌ Erreur lors de la création de {master_name} : {e}")
            return None

    def validate_master_file(self, file_path, master_name):
        """Vérifie si un fichier maître a bien été créé et est en float32."""
        try:
            with fits.open(file_path) as hdul:
                data = hdul[0].data
                dtype = data.dtype
                mean_value = np.mean(data)
            if dtype != np.float32:
                logging.warning(f"⚠️ {master_name} enregistré en {dtype} au lieu de float32 !")
            if "bias" in master_name.lower():
                logging.info(f"✅ Master Bias vérifié : Moyenne = {mean_value}, Type = {dtype}")
            elif "dark" in master_name.lower():
                logging.info(f"✅ Master Dark vérifié : Moyenne = {mean_value}, Type = {dtype}")
            elif "flat" in master_name.lower():
                logging.info(f"✅ Master Flat vérifié : Moyenne = {mean_value}, Type = {dtype}")
        except Exception as e:
            logging.error(f"❌ Erreur lors de la vérification de {master_name} : {e}")

    def update_progress(self, value, total, message):
        """Mise à jour de la barre de progression et du texte (si les widgets sont définis)."""
        if self.progress_bar and self.progress_label:
            progress = (value / total) * 100
            self.progress_bar['value'] = progress
            self.progress_label.config(text=f"{message} ({progress:.1f}%)")
            self.progress_bar.update_idletasks()
            time.sleep(0.1)
