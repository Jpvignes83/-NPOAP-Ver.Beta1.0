# gui/target_selector.py
import logging
import numpy as np
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
from astropy.io import fits
from astropy.visualization import ZScaleInterval
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.vizier import Vizier
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
from matplotlib.backend_bases import MouseButton 
import pandas as pd
from astropy.wcs.utils import proj_plane_pixel_scales

# Import des outils d'analyse depuis le pipeline
from core.photometry_pipeline import (
    launch_photometry_aperture, 
    refine_centroid, 
    estimate_fwhm_marginal, 
    show_diagnostics_windows
)

# --- CONFIGURATION VIZIER ---
Vizier.VIZIER_SERVER = "vizier.cfa.harvard.edu" 
Vizier.ROW_LIMIT = 5000 # On demande beaucoup d'étoiles pour filtrer géométriquement après
# ----------------------------

# Stockage pour fermer les fenêtres de diagnostic précédentes
_current_diag_windows = []

def close_previous_diagnostics():
    """Ferme les fenêtres de diagnostic ouvertes précédemment."""
    global _current_diag_windows
    for win in _current_diag_windows:
        try:
            if win.winfo_exists():
                win.destroy()
        except:
            pass
    _current_diag_windows = []

def launch_target_selection(
    fits_path,
    on_selection_done=launch_photometry_aperture,
    is_asteroid=False,
    ephemeris_data=None,
    header_target_coord=None
):
    """
    Lance la fenêtre de sélection de T1 et des étoiles de comparaison.
    
    Parameters
    ----------
    fits_path : str or Path
        Chemin vers l'image FITS de référence
    on_selection_done : callable, optional
        Fonction callback appelée après la sélection (default: launch_photometry_aperture)
    is_asteroid : bool, optional
        Si True, T1 est un astéroïde (ne pas chercher dans Gaia, utiliser magnitude éphémérides)
    ephemeris_data : Table, optional
        Table d'éphémérides avec colonnes 'datetime_jd' et 'V' pour les astéroïdes
    header_target_coord : SkyCoord, optional
        Coordonnées de la cible lues depuis l'entête FITS pour affichage.
    """
    fits_path = Path(fits_path)

    # 1. Chargement FITS
    try:
        with fits.open(fits_path) as hdul:
            data = hdul[0].data.astype(float) 
            header = hdul[0].header
            wcs = WCS(header)
    except Exception as e:
        messagebox.showerror("Erreur FITS", f"Impossible d'ouvrir l'image :\n{e}")
        return

    # 2. ZScale
    try:
        interval = ZScaleInterval()
        vmin, vmax = interval.get_limits(data)
    except:
        vmin, vmax = np.percentile(data, [10, 90])

    # 3. Champ de recherche (Rayon fixe de 15 arcmin)
    search_radius = 15 * u.arcmin

    # 4. Interface
    root = tk.Toplevel()
    root.title(f"Sélection Cibles - {fits_path.name}")
    root.geometry("1450x850") 

    # --- Zone Gauche (Liste) ---
    left_frame = tk.Frame(root, width=500, bg="#f0f0f0") 
    left_frame.pack(side=tk.LEFT, fill=tk.Y)
    left_frame.pack_propagate(False)

    tk.Label(left_frame, text=" Étoiles Sélectionnées (Coordonnées Gaia) ", bg="#ddd", font=("Arial", 10, "bold"), pady=5).pack(fill=tk.X)
    
    # En-têtes colonnes
    h_frame = tk.Frame(left_frame, bg="#ccc")
    h_frame.pack(fill="x", padx=2, pady=2)
    tk.Label(h_frame, text="", width=5, bg="#ccc").pack(side="left") 
    tk.Label(h_frame, text="ID", width=4, bg="#ccc", font=("Arial",9,"bold")).pack(side="left")
    tk.Label(h_frame, text="Mag", width=6, bg="#ccc", font=("Arial",9,"bold")).pack(side="left")
    tk.Label(h_frame, text="BP-RP", width=6, bg="#ccc", font=("Arial",9,"bold")).pack(side="left")
    tk.Label(h_frame, text="RA (deg)", width=12, bg="#ccc", font=("Arial",9,"bold")).pack(side="left")
    tk.Label(h_frame, text="Dec (deg)", width=12, bg="#ccc", font=("Arial",9,"bold")).pack(side="left")

    list_container = tk.Frame(left_frame, bg="#f0f0f0")
    list_container.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=5, pady=5)

    btn_frame = tk.Frame(left_frame, bg="#e0e0e0", pady=10)
    btn_frame.pack(side=tk.BOTTOM, fill=tk.X)
    
    # --- Zone Droite (Graphique) ---
    right_frame = tk.Frame(root, bg="black")
    right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

    status_label = tk.Label(right_frame, text="Mode: SÉLECTION (Cliquez gauche T1, Droit Diagnostic)", 
                            bg="black", fg="#00ff00", font=("Consolas", 11, "bold"), pady=5)
    status_label.pack(side=tk.TOP, fill=tk.X)

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': wcs})
    ax.imshow(data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
    
    def format_coord(x, y):
        try:
            c = wcs.pixel_to_world(x, y)
            return f"RA={c.ra.deg:.5f} Dec={c.dec.deg:.5f}"
        except: return ""
    ax.format_coord = format_coord

    def _parse_header_coord(header):
        ra_h = header.get("OBJCTRA") or header.get("RA")
        dec_h = header.get("OBJCTDEC") or header.get("DEC")
        if not ra_h or not dec_h:
            return None
        
        try:
            ra_val = float(ra_h)
            dec_val = float(dec_h)
            if ra_val > 24.0:
                return SkyCoord(ra_val, dec_val, unit=(u.deg, u.deg))
        except Exception:
            pass
        
        try:
            return SkyCoord(ra_h, dec_h, unit=(u.hourangle, u.deg))
        except Exception:
            try:
                return SkyCoord(float(ra_h), float(dec_h), unit=(u.deg, u.deg))
            except Exception:
                return None

    # --- Header T1 ---
    try:
        target_coord_h = header_target_coord
        if target_coord_h is None:
            target_coord_h = _parse_header_coord(header)
        
        if target_coord_h is not None:
            tx, ty = wcs.world_to_pixel(target_coord_h)
            if np.isfinite(tx) and np.isfinite(ty):
                c = Circle(
                    (tx, ty),
                    radius=30,
                    edgecolor='yellow',
                    ls='--',
                    lw=1.5,
                    fill=False,
                    label='header_t1'
                )
                ax.add_patch(c)
    except Exception as e:
        logging.warning(f"Impossible d'afficher T1 du header: {e}")

    canvas = FigureCanvasTkAgg(fig, master=right_frame)
    canvas.draw()

    toolbar = NavigationToolbar2Tk(canvas, right_frame)
    toolbar.update()
    toolbar.pack(side=tk.BOTTOM, fill=tk.X)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    # Changer le curseur pour le rendre plus visible en mode zoom
    # Utiliser "pencil" qui est plus visible
    def update_cursor(event):
        mode = toolbar.mode if hasattr(toolbar, 'mode') else ""
        if mode and mode != "" and 'zoom' in mode.lower():
            canvas_widget.config(cursor="pencil")
        else:
            canvas_widget.config(cursor="")
    
    canvas.mpl_connect('motion_notify_event', update_cursor)

    # Variables d'état
    state = {"comps": [], "vars": [], "mag_t1": None, "bp_rp_t1": None, "target": None, "fwhm_t1": None}
    
    def update_list():
        for w in list_container.winfo_children(): w.destroy()
        
        # --- Ligne T1 ---
        if state["target"]:
            row = tk.Frame(list_container, bg="#f0f0f0")
            row.pack(fill="x", pady=2)
            
            t1_var = tk.BooleanVar(value=True)
            row.t1_var = t1_var 
            chk = tk.Checkbutton(row, variable=t1_var, state="disabled", bg="#f0f0f0", width=2)
            chk.pack(side="left")
            
            tk.Label(row, text="T1", width=4, fg="red", font=("Arial",9,"bold"), bg="#f0f0f0").pack(side="left")
            
            mag = state.get("mag_t1")
            mag_val = f"{mag:.2f}" if mag else "?"
            tk.Label(row, text=mag_val, width=6, bg="#f0f0f0").pack(side="left")
            
            bprp = state.get("bp_rp_t1")
            bprp_val = f"{bprp:.2f}" if bprp is not None else "-"
            tk.Label(row, text=bprp_val, width=6, bg="#f0f0f0").pack(side="left")
            
            t = state["target"]
            tk.Label(row, text=f"{t.ra.deg:.5f}", width=12, bg="#eef").pack(side="left", padx=1)
            tk.Label(row, text=f"{t.dec.deg:.5f}", width=12, bg="#eef").pack(side="left", padx=1)
        
        # --- Lignes Comparateurs ---
        if state["comps"]:
            for i, comp_data in enumerate(state["comps"]):
                bp_rp = None
                if isinstance(comp_data, dict):
                    coord = comp_data["coord"]
                    mag = comp_data["mag"]
                    bp_rp = comp_data.get("bp_rp")
                elif len(comp_data) == 3:
                    coord, mag, bp_rp = comp_data
                else:
                    coord, mag = comp_data 

                f = tk.Frame(list_container, bg="#f0f0f0")
                f.pack(fill="x", pady=1)
                
                v = state["vars"][i]
                
                tk.Checkbutton(f, variable=v, bg="#f0f0f0", width=2).pack(side="left") 
                tk.Label(f, text=f"C{i+1}", width=4, fg="blue", bg="#f0f0f0").pack(side="left")
                
                m_txt = f"{mag:.2f}" if mag else "?"
                tk.Label(f, text=m_txt, width=6, bg="#f0f0f0").pack(side="left")
                
                br_txt = f"{bp_rp:.2f}" if bp_rp is not None else "-"
                tk.Label(f, text=br_txt, width=6, bg="#f0f0f0").pack(side="left")
                
                tk.Label(f, text=f"{coord.ra.deg:.5f}", width=12, bg="#fff").pack(side="left", padx=1)
                tk.Label(f, text=f"{coord.dec.deg:.5f}", width=12, bg="#fff").pack(side="left", padx=1)

    def on_click(event):
        # Le clic droit doit fonctionner même en mode zoom/pan
        is_right_click = event.button == MouseButton.RIGHT
        
        # Bloquons seulement le clic gauche si le mode toolbar est actif
        if not is_right_click and toolbar.mode != "" and toolbar.mode is not None:
            status_label.config(text=f"⚠️ MODE {str(toolbar.mode).upper()} ACTIF ! Désactivez la loupe.", fg="white", bg="red")
            return
        
        if not event.inaxes: return

        x_raw, y_raw = event.xdata, event.ydata
        if not np.isfinite(x_raw) or not np.isfinite(y_raw): return
            
        try:
            x_c, y_c = refine_centroid(data, x_raw, y_raw)
            if not np.isfinite(x_c) or not np.isfinite(y_c): raise ValueError("Centroidage échoué.")
            
            fwhm_info = estimate_fwhm_marginal(data, x_c, y_c)
            fwhm_val = fwhm_info[0] if fwhm_info and fwhm_info[0] else 0.0
            coord = wcs.pixel_to_world(x_c, y_c)
            
        except Exception as e:
            logging.error(f"Erreur de centroidage/FWHM: {e}")
            return
            
        if event.button == MouseButton.RIGHT:
            try:
                logging.info(f"Clic droit détecté à ({x_c:.1f}, {y_c:.1f}), fwhm_info={fwhm_info}")
                
                # Stocker le FWHM si c'est sur T1
                if state["target"] and fwhm_info and fwhm_info[0] is not None:
                    # Vérifier si le clic est proche de T1 (dans un rayon de 50 pixels)
                    tx, ty = wcs.world_to_pixel(state["target"])
                    dist = np.sqrt((x_c - tx)**2 + (y_c - ty)**2)
                    if dist < 50:  # Si le clic droit est sur T1
                        state["fwhm_t1"] = fwhm_info[0]
                        logging.info(f"FWHM T1 mesuré par clic droit : {state['fwhm_t1']:.2f} px")
                
                close_previous_diagnostics() 
                new_windows = show_diagnostics_windows(data, x_c, y_c, fwhm_info, parent=root) 
                logging.info(f"show_diagnostics_windows a retourné {len(new_windows) if new_windows else 0} fenêtres")
                if new_windows:
                    _current_diag_windows.extend(new_windows)
                    status_label.config(text=f"Diagnostic @ {x_c:.1f},{y_c:.1f}", bg="#ccffcc")
                else:
                    status_label.config(text=f"❌ Impossible d'afficher le diagnostic (FWHM invalide?)", bg="orange")
            except Exception as e:
                logging.error(f"Erreur lors de l'affichage du diagnostic: {e}", exc_info=True)
                status_label.config(text=f"❌ Erreur diagnostic: {e}", bg="red")
            
        elif event.button == MouseButton.LEFT:
            status_label.config(text="Interrogation Gaia...", fg="black", bg="orange")
            root.config(cursor="watch")
            root.update()

            try:
                state["comps"] = []; state["vars"] = []
                state["mag_t1"] = None; state["bp_rp_t1"] = None
                
                to_remove = [p for p in ax.patches if p.get_label() != 'header_t1']
                for p in to_remove: p.remove()
                [t.remove() for t in ax.texts]

                state["target"] = coord
                ax.add_patch(Circle((x_c, y_c), 20, edgecolor='red', lw=2, fill=False))
                ax.annotate("T1", (x_c+25, y_c+25), color='red', fontweight='bold')
                canvas.draw()

                # --- TRAITEMENT T1 ---
                if is_asteroid:
                    # Pour les astéroïdes : utiliser directement les coordonnées cliquées (pas de recherche Gaia)
                    state["bp_rp_t1"] = None
                    
                    # Obtenir la magnitude depuis les éphémérides Horizons
                    if ephemeris_data is not None:
                        try:
                            # Récupérer la date d'observation depuis le header
                            from astropy.time import Time
                            date_obs = header.get("DATE-OBS") or header.get("DATE")
                            if date_obs:
                                t_obs = Time(date_obs, scale='utc')
                                jd_obs = t_obs.jd
                                
                                # Interpoler la magnitude V depuis les éphémérides
                                if 'datetime_jd' in ephemeris_data.colnames and 'V' in ephemeris_data.colnames:
                                    eph_jds = ephemeris_data['datetime_jd']
                                    mag_v = ephemeris_data['V']
                                    state["mag_t1"] = float(np.interp(jd_obs, eph_jds, mag_v))
                                    logging.info(f"T1 (astéroïde) : Magnitude V={state['mag_t1']:.2f} depuis éphémérides Horizons (JD={jd_obs:.6f})")
                                else:
                                    logging.warning("Colonnes datetime_jd ou V manquantes dans éphémérides")
                                    state["mag_t1"] = None
                            else:
                                logging.warning("DATE-OBS manquant dans header, impossible d'interpoler magnitude")
                                state["mag_t1"] = None
                        except Exception as e:
                            logging.error(f"Erreur interpolation magnitude T1 depuis éphémérides: {e}", exc_info=True)
                            state["mag_t1"] = None
                    else:
                        logging.warning("Éphémérides non disponibles, magnitude T1 non disponible")
                        state["mag_t1"] = None
                else:
                    # Pour les étoiles (exoplanètes) : recherche dans Gaia
                    # --- REQUETE GAIA ---
                    # Inclure phot_variable_flag pour exclure les variables
                    v = Vizier(columns=["RA_ICRS", "DE_ICRS", "Gmag", "bp_rp", "phot_variable_flag", "**"], row_limit=5000)
                    
                    # 1. T1 - Recherche de l'étoile Gaia la plus proche
                    # Rayon de 2.5 arcsec (approprié pour un FWHM typique de 1.2 arcsec)
                    res_t1 = v.query_region(state["target"], radius=2.5*u.arcsec, catalog="I/355/gaiadr3")
                    
                    if res_t1 and len(res_t1) > 0 and len(res_t1[0]) > 0:
                        # Calculer la distance de chaque étoile Gaia par rapport à la position cliquée
                        target_coord = state["target"]  # Coordonnée de la position cliquée
                        gaia_table = res_t1[0]
                        
                        # Créer des SkyCoord pour toutes les étoiles Gaia trouvées
                        gaia_coords = SkyCoord(ra=gaia_table['RA_ICRS'], dec=gaia_table['DE_ICRS'], unit=(u.deg, u.deg))
                        
                        # Calculer les séparations angulaires
                        separations = target_coord.separation(gaia_coords)
                        
                        # Trouver l'étoile la plus proche
                        idx_closest = np.argmin(separations)
                        t1_row = gaia_table[idx_closest]
                        separation_arcsec = separations[idx_closest].arcsec
                        
                        logging.info(f"T1: {len(gaia_table)} étoiles Gaia trouvées dans 2.5\". Étoile la plus proche: {separation_arcsec:.2f}\"")
                        
                        if "Gmag" in t1_row.colnames: state["mag_t1"] = float(t1_row["Gmag"])
                        
                        if "bp_rp" in t1_row.colnames and not np.ma.is_masked(t1_row["bp_rp"]):
                             state["bp_rp_t1"] = float(t1_row["bp_rp"])
                        elif "BP-RP" in t1_row.colnames and not np.ma.is_masked(t1_row["BP-RP"]):
                             state["bp_rp_t1"] = float(t1_row["BP-RP"])
                        
                        if state["mag_t1"] is None:
                            logging.error("T1 sélectionnée sans magnitude Gmag.")
                            messagebox.showerror("Données Manquantes", "Cette étoile n'a pas de magnitude (Gmag).\nVeuillez choisir une autre cible.")
                            root.config(cursor="")
                            return
                        
                        # Avertir si la séparation est importante
                        if separation_arcsec > 2.0:
                            logging.warning(f"T1: Étoile Gaia la plus proche à {separation_arcsec:.2f}\" de la position cliquée")
                    else:
                        messagebox.showerror("Erreur", "Aucune étoile Gaia trouvée pour T1 à cette position.")
                        root.config(cursor="")
                        return

                # --- 2. COMPARATEURS ---
                # Requête Gaia pour les comparateurs (toujours nécessaire)
                v = Vizier(columns=["RA_ICRS", "DE_ICRS", "Gmag", "bp_rp", "phot_variable_flag", "**"], row_limit=5000)
                res_c = v.query_region(state["target"], radius=search_radius, catalog="I/355/gaiadr3")
                
                if res_c and len(res_c[0]) > 0:
                    cand = []
                    t1_mag = state["mag_t1"]  # Peut être None pour astéroïde si éphémérides non disponibles
                    
                    # --- FILTRE 1 : DANS L'IMAGE ? ---
                    h_img, w_img = data.shape
                    margin = 20 # Marge de sécurité en pixels (bord)

                    for r in res_c[0]:
                            try:
                                if "Gmag" not in r.colnames or np.ma.is_masked(r["Gmag"]): continue
                                m = float(r["Gmag"])
                                c_sky = SkyCoord(ra=r["RA_ICRS"]*u.deg, dec=r["DE_ICRS"]*u.deg)
                                
                                # Conversion en pixels
                                cx, cy = wcs.world_to_pixel(c_sky)
                                
                                # VERIFICATION GEOMETRIQUE STRICTE
                                # Si l'étoile est hors du rectangle de l'image (moins la marge), on l'oublie
                                if not (margin <= cx <= w_img - margin and margin <= cy <= h_img - margin):
                                    continue
                                
                                # Si on est ici, l'étoile est SUR le capteur.
                                
                                bprp = None
                                if "bp_rp" in r.colnames and not np.ma.is_masked(r["bp_rp"]):
                                    bprp = float(r["bp_rp"])
                                elif "BP-RP" in r.colnames and not np.ma.is_masked(r["BP-RP"]):
                                    bprp = float(r["BP-RP"])
                                
                                sep = c_sky.separation(state["target"]).arcsec
                                
                                # Exclure les étoiles variables
                                is_variable = False
                                if "phot_variable_flag" in r.colnames and not np.ma.is_masked(r["phot_variable_flag"]):
                                    var_flag = r["phot_variable_flag"]
                                    if var_flag is not None and var_flag != '':
                                        var_str = str(var_flag).strip().upper()
                                        if var_str in ['VARIABLE', 'TRUE', '1', 'Y', 'YES']:
                                            is_variable = True
                                
                                # Filtre distance T1, magnitude proche (si disponible), et exclusion des variables
                                if sep > 10 and not is_variable:
                                    if t1_mag is not None:
                                        # Filtre par magnitude si disponible
                                        if abs(m - t1_mag) <= 2.5:
                                            cand.append((c_sky, m, bprp))
                                    else:
                                        # Pas de filtre magnitude si T1 mag non disponible
                                        cand.append((c_sky, m, bprp))
                            except: continue
                        
                    # --- FILTRE 2 : TRIER PAR MAGNITUDE (Les 15 meilleurs) ---
                    if t1_mag is not None:
                        cand = sorted(cand, key=lambda x: abs(x[1]-t1_mag))[:15]
                    else:
                        # Si pas de magnitude T1, trier par magnitude absolue (les plus brillantes)
                        cand = sorted(cand, key=lambda x: x[1])[:15]
                    
                    for i, (cc, mm, val_bp) in enumerate(cand):
                        state["comps"].append((cc, mm, val_bp))
                        state["vars"].append(tk.BooleanVar(value=True))
                        
                        cx, cy = wcs.world_to_pixel(cc)
                        ax.add_patch(Rectangle((cx-15, cy-15), 30, 30, edgecolor='cyan', fill=False))
                        ax.annotate(f"C{i+1}", (cx+18, cy), color='cyan', fontsize=9)
                else:
                    logging.warning("Aucune étoile trouvée dans la zone large.")
                
                update_list()
                msg = f"T1 (Mag: {state['mag_t1'] or '?'}). {len(state['comps'])} comparateurs."
                status_label.config(text=msg, fg="#00ff00", bg="black")
                
            except Exception as e:
                logging.error(f"Erreur Gaia: {e}")
                status_label.config(text=f"Erreur: {e}", fg="white", bg="red")
            finally:
                root.config(cursor="")
                canvas.draw()

    fig.canvas.mpl_connect("button_press_event", on_click)

    def validate():
        if state["target"]:
            final_comps = [c[0] for i, c in enumerate(state["comps"]) if state["vars"][i].get()]
            close_previous_diagnostics()
            
            # Créer target_data avec coord et FWHM (si mesuré)
            target_data_dict = {"coord": state["target"]}
            if state["fwhm_t1"] is not None:
                target_data_dict["fwhm"] = state["fwhm_t1"]
                logging.info(f"FWHM T1 transmis à launch_photometry_aperture : {state['fwhm_t1']:.2f} px")
            
            on_selection_done(fits_path, target_data_dict, final_comps)
            root.destroy()
        else:
            messagebox.showwarning("Stop", "Définissez T1 avant de valider.")

    ttk.Button(btn_frame, text="✅ VALIDER LA SÉLECTION", command=validate).pack(fill=tk.X, ipady=5, padx=5)
    root.mainloop()

class TargetSelector:
    def __init__(self, fits_path=None, on_selection_done=launch_photometry_aperture):
        self.fits_path = fits_path
        self.on_selection_done = on_selection_done
    def run(self, fits_path=None, on_selection_done=None):
        launch_target_selection(fits_path or self.fits_path, on_selection_done or self.on_selection_done)
    @staticmethod
    def open(fits_path, on_selection_done=launch_photometry_aperture):
        launch_target_selection(fits_path, on_selection_done)