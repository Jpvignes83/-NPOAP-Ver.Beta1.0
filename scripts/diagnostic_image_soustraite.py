#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Diagnostic rapide d'une image soustraite (ZOGY ou autre).
Usage: python diagnostic_image_soustraite.py "chemin/vers/subtracted_xxx.fits"
"""

import sys
import numpy as np

def main():
    if len(sys.argv) < 2:
        print("Usage: python diagnostic_image_soustraite.py <fichier_fits_soustrait>")
        print("Exemple: python diagnostic_image_soustraite.py \"d:\\AT 2026fbz\\...\\subtracted_xxx.fits\"")
        sys.exit(1)
    path = sys.argv[1]
    try:
        from astropy.io import fits
    except ImportError:
        print("Erreur: astropy requis (pip install astropy)")
        sys.exit(1)
    try:
        with fits.open(path) as hdul:
            data = hdul[0].data.astype(float)
            header = hdul[0].header
    except Exception as e:
        print(f"Erreur lecture FITS: {e}")
        sys.exit(1)
    # Stats
    valid = np.isfinite(data)
    if not np.any(valid):
        print("Aucune valeur finie dans l'image.")
        sys.exit(1)
    d = data[valid]
    n = d.size
    pos = np.sum(d > 0)
    neg = np.sum(d < 0)
    zero = np.sum(d == 0)
    print("=" * 60)
    print("DIAGNOSTIC IMAGE SOUSTRAITE")
    print("=" * 60)
    print(f"Fichier    : {path}")
    print(f"Dimensions : {data.shape}")
    print("")
    print("Statistiques (pixels valides):")
    print(f"  Min      : {np.nanmin(data):.4g}")
    print(f"  Max      : {np.nanmax(data):.4g}")
    print(f"  Moyenne  : {np.nanmean(data):.4g}")
    print(f"  Médiane  : {np.nanmedian(data):.4g}")
    print(f"  Écart-type : {np.nanstd(data):.4g}")
    print("")
    print("Répartition du signe:")
    print(f"  Pixels > 0 : {pos:12d} ({100*pos/n:.1f} %)")
    print(f"  Pixels < 0 : {neg:12d} ({100*neg/n:.1f} %)")
    print(f"  Pixels = 0 : {zero:12d} ({100*zero/n:.1f} %)")
    print("")
    # Interprétation
    print("Interprétation attendue pour une soustraction correcte:")
    print("  - Médiane proche de 0 (fond soustrait).")
    print("  - Répartition ~50% positif / ~50% négatif (bruit centré).")
    print("  - Un transitoire (nouvelle source) apparaît en PEAK POSITIF.")
    print("  - Une source absente en science mais en référence → négatif.")
    print("")
    med = np.nanmedian(data)
    if abs(med) > 0.1 * np.nanstd(data):
        print(f"  [Info] Médiane = {med:.4g} : écart au zéro détecté (normal si fond variable).")
    if pos > 0.6 * n or neg > 0.6 * n:
        print("  [Attention] Fort déséquilibre positif/négatif : vérifier alignement ou méthode.")
    else:
        print("  [OK] Répartition des signes cohérente avec une soustraction.")
    print("=" * 60)

if __name__ == "__main__":
    main()
