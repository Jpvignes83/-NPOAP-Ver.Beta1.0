import os
import logging
from astropy.io import fits


def get_fits_files(directory):
    """Récupère la liste des fichiers FITS dans un dossier."""
    return [os.path.join(directory, f) for f in os.listdir(directory) if f.lower().endswith(".fits")]


def load_fits(image_file):
    """Charge un fichier FITS et retourne ses données et son en-tête."""
    try:
        with fits.open(image_file) as hdul:
            data = hdul[0].data
            header = hdul[0].header

            if data is None:
                logging.error(f"L’image FITS {image_file} ne contient pas de données valides.")
                return None, None

            return data, header
    except Exception as e:
        logging.error(f"Erreur lors de l'ouverture du FITS {image_file}: {e}")
        return None, None


def save_fits(data, header, output_path):
    """Enregistre les données et le header dans un nouveau fichier FITS."""
    try:
        fits.writeto(output_path, data, header, overwrite=True)
        logging.info(f"Fichier FITS enregistré: {output_path}")
    except Exception as e:
        logging.error(f"Erreur lors de l'écriture du FITS {output_path}: {e}")


def save_calibrated_fits(raw_path, master_bias_path, master_dark_path, master_flat_path, output_path):
    """
    Calibre une image brute (RAW) en appliquant bias, dark et flat,
    puis enregistre le FITS calibré avec mise à jour de l'en-tête.
    """
    # Chargement des données et header original
    data, header = load_fits(raw_path)
    if data is None or header is None:
        logging.error(f"Abandon de la calibration pour {raw_path}: données non valides.")
        return

    # Chargement des masters
    bias_data, _ = load_fits(master_bias_path)
    dark_data, _ = load_fits(master_dark_path)
    flat_data, _ = load_fits(master_flat_path)
    if None in (bias_data, dark_data, flat_data):
        logging.error(f"Masters non valides pour {raw_path}. Calibration interrompue.")
        return

    # Conversion en float pour éviter les under/overflow
    calibrated = (data.astype(float) - bias_data.astype(float) - dark_data.astype(float))
    calibrated /= flat_data.astype(float)

    # Mise à jour de l'en-tête
    header['BIASCORR'] = ('YES', 'Bias soustrait')
    header['BIASFILE'] = (os.path.basename(master_bias_path), 'Master bias utilisé')
    header['DARKCORR'] = ('YES', 'Dark soustrait')
    header['DARKFILE'] = (os.path.basename(master_dark_path), 'Master dark utilisé')
    header['FLATCORR'] = ('YES', 'Flat appliqué')
    header['FLATFILE'] = (os.path.basename(master_flat_path), 'Master flat utilisé')
    header.add_history(
        f"Calibration appliquée: Bias={os.path.basename(master_bias_path)}, "
        f"Dark={os.path.basename(master_dark_path)}, "
        f"Flat={os.path.basename(master_flat_path)}"
    )

    # Enregistrement
    save_fits(calibrated, header, output_path)


class FileHandler:
    """
    Wrapper orienté objet autour des fonctions utilitaires de ce module.

    Permet d'utiliser:
        from utils import FileHandler
        fh = FileHandler()
        fh.get_fits_files(...)
    """

    @staticmethod
    def get_fits_files(directory):
        return get_fits_files(directory)

    @staticmethod
    def load_fits(image_file):
        return load_fits(image_file)

    @staticmethod
    def save_fits(data, header, output_path):
        return save_fits(data, header, output_path)

    @staticmethod
    def save_calibrated_fits(raw_path, master_bias_path, master_dark_path, master_flat_path, output_path):
        return save_calibrated_fits(raw_path, master_bias_path, master_dark_path, master_flat_path, output_path)
