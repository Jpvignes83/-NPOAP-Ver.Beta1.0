import os
import csv
import numpy as np

# Notez l'ajout de period_input=None ici
def extraire_et_sauvegarder(dossier_source, fichier_sortie="mid-time.csv", period_input=None):
    """
    1. Extrait mid-time, uncertainty.
    2. Utilise period_input (Littérature) SI fournie, sinon calcule la médiane.
    3. Calcule Epochs et O-C.
    4. Sauvegarde le CSV.
    """
    if not dossier_source or not os.path.exists(dossier_source):
        print("Erreur : Dossier invalide.")
        return

    donnees_brutes = []

    # 1. EXTRACTION
    fichiers = [f for f in os.listdir(dossier_source) 
                if os.path.isfile(os.path.join(dossier_source, f)) 
                and f != fichier_sortie]

    for nom_fichier in fichiers:
        chemin_complet = os.path.join(dossier_source, nom_fichier)
        mid_time = None
        uncertainty = None
        period_local = None

        try:
            with open(chemin_complet, 'r', encoding='utf-8', errors='ignore') as f:
                lignes = f.readlines()
                for ligne in lignes:
                    ligne_clean = ligne.strip()
                    # Extraction Mid-time
                    if "Final mid-time (BJD_TDB)" in ligne_clean:
                        parts = ligne_clean.replace('=', ' ').split()
                        if parts: mid_time = float(parts[-1])
                    # Extraction Uncertainty
                    if "Final mid-time uncertainty" in ligne_clean:
                        parts = ligne_clean.replace('=', ' ').split()
                        if parts: uncertainty = float(parts[-1])
                    # Extraction Période locale (si dispo)
                    if "P (days):" in ligne_clean:
                        parts = ligne_clean.split(':')
                        if len(parts) > 1:
                            try: period_local = float(parts[1].strip())
                            except ValueError: period_local = None

            if mid_time is not None:
                if uncertainty is None: uncertainty = 0.0
                donnees_brutes.append({
                    'fichier': nom_fichier, 
                    'time': mid_time, 
                    'err': uncertainty, 
                    'p_raw': period_local
                })
        except Exception as e:
            print(f"Erreur lecture {nom_fichier}: {e}")

    if not donnees_brutes:
        print("Aucune donnée extraite.")
        return

    # 2. CALCUL DES PARAMÈTRES (T0, P)
    donnees_brutes.sort(key=lambda x: x['time'])
    times = np.array([d['time'] for d in donnees_brutes])

    # A) T0 = Première valeur de mid-time
    T0 = times[0]

    # B) CHOIX DE LA PÉRIODE (CORRIGÉ)
    if period_input is not None and period_input > 0:
        # Cas 1 : L'utilisateur a donné une valeur manuelle (Littérature)
        P_final = float(period_input)
        source_p = "Manuel (Littérature)"
    else:
        # Cas 2 : Calcul de la Période par Régression Linéaire
        
        # 1. Calcul des epochs initiales (basées sur une estimation grossière)
        if len(times) > 1:
            diffs = np.diff(times)
            P_guess = np.median(diffs) if len(diffs) > 0 else 1.0
        else:
            P_guess = 1.0

        # Ces epochs servent juste de variable X pour la régression
        epochs_ref = np.round((times - T0) / P_guess)

        # 2. Régression linéaire : Time = T0_new + P_new * Epoch
        try:
            # np.polyfit(X, Y, degré) -> polyfit(epochs, times, 1)
            # poly[0] = Pente (P_new), poly[1] = Ordonnée (T0_new)
            poly = np.polyfit(epochs_ref, times, 1)
            P_final = poly[0]
            T0 = poly[1] # OPTIONNEL: On peut même utiliser le T0 ajusté du fit
            source_p = "Auto (Régression Linéaire)"
        except Exception as e:
            print(f"Erreur de régression : {e}. Fallback vers médiane.")
            # Fallback si la régression échoue
            all_periods = [d['p_raw'] for d in donnees_brutes if d['p_raw'] is not None]
            P_final = np.median(all_periods) if all_periods else 1.0
            source_p = "Auto (Fallback Médiane)"

    print(f"Paramètres utilisés -> T0: {T0:.6f}, P: {P_final:.8f} [{source_p}]")

    # 3. CALCUL DES EPOCHS ET O-C
    # Utilisation de la nouvelle T0 et P_final
    if P_final != 0:
        # Le calcul de l'Epoch est refait avec la P_final trouvée
        epochs = np.round((times - T0) / P_final)
    else:
        epochs = np.zeros(len(times))
    
    # C = T0 + Epoch * P_final
    calculated = T0 + epochs * P_final
    
    # (O-C) = Observé - Calculé
    oc_values = times - calculated

    # 4. ÉCRITURE DU CSV
    chemin_sortie = os.path.join(dossier_source, fichier_sortie)
    try:
        with open(chemin_sortie, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['Epoch', 'O-C', 'Mid-time', 'Uncertainty', 'Fichier'])
            for i, d in enumerate(donnees_brutes):
                writer.writerow([
                    int(epochs[i]),
                    f"{oc_values[i]:.6f}",
                    f"{d['time']:.6f}",
                    f"{d['err']:.6f}",
                    d['fichier']
                ])
        print(f"Fichier créé : {chemin_sortie}")
    except Exception as e:
        print(f"Erreur écriture CSV : {e}")