#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script KBMOD a executer sous WSL (Windows Subsystem for Linux).
Lit les FITS d'un dossier (chemin WSL), lance la detection Synthetic Tracking,
et ecrit les candidats dans kbmod_candidates.csv dans ce meme dossier.

Usage (depuis WSL ou via: wsl python3 /mnt/c/.../scripts/kbmod_wsl_detect.py ...):
  python3 kbmod_wsl_detect.py <dossier_fits_wsl> [--vmin 0.1] [--vmax 10.0] [--scale 1.0] [--max-results 200] [--gpu]
      [--min-lh 10] [--static-mask-sigma 8] [--static-mask-fraction 0.4] [--no-static-mask]

Le CSV produit a les colonnes: ra_deg,dec_deg,x0,y0,vx_px_per_day,vy_px_per_day,likelihood,jd_ref
"""
import argparse
import csv
import sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="KBMOD detection (run under WSL)")
    parser.add_argument("fits_dir", help="WSL path to folder containing FITS files")
    parser.add_argument("--vmin", type=float, default=0.02, help="Min velocity arcsec/min (0.02 ~ 76 px/j @ 0.38\"/px)")
    parser.add_argument("--vmax", type=float, default=60.0, help="Max velocity arcsec/min (60 = astéroïdes ceinture)")
    parser.add_argument("--scale", type=float, default=1.0, help="Scale arcsec/px")
    parser.add_argument("--max-results", type=int, default=200, help="Max candidates")
    parser.add_argument("--gpu", action="store_true", help="Use GPU for search (requires CUDA in WSL)")
    parser.add_argument("--min-lh", type=float, default=10.0, help="Minimum likelihood retained in output")
    parser.add_argument("--static-mask-sigma", type=float, default=8.0, help="Robust sigma threshold for static mask")
    parser.add_argument("--static-mask-fraction", type=float, default=0.4, help="Min frame fraction to mark a pixel as static")
    parser.add_argument("--no-static-mask", action="store_true", help="Disable robust static mask (recommended for very slow targets)")
    parser.add_argument("--dedup-pos-px", type=float, default=3.0, help="Duplicate grouping radius in starting position (pixels)")
    parser.add_argument("--dedup-vel-pxday", type=float, default=5.0, help="Duplicate grouping radius in velocity (px/day)")
    args = parser.parse_args()

    fits_dir = Path(args.fits_dir)
    if not fits_dir.is_dir():
        print(f"ERROR: Not a directory: {fits_dir}", file=sys.stderr)
        sys.exit(1)

    try:
        from kbmod.core.psf import PSF
    except ImportError as e:
        print(f"ERROR: KBMOD not installed in this environment: {e}", file=sys.stderr)
        print("See docs/INSTALL_KBMOD_WSL.md for installation under WSL.", file=sys.stderr)
        sys.exit(1)

    import numpy as np
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.time import Time

    import kbmod.search as kbmod_search
    StackSearch = getattr(kbmod_search, "StackSearch", None)
    fill_psi_phi = getattr(kbmod_search, "fill_psi_phi_array_from_image_arrays", None)
    SearchParameters = getattr(kbmod_search, "SearchParameters", None)
    search_cpu_only = getattr(kbmod_search, "search_cpu_only", None)
    search_gpu = getattr(kbmod_search, "search", None) or getattr(kbmod_search, "search_gpu", None)
    PsiPhiArray = getattr(kbmod_search, "PsiPhiArray", None)
    TrajectoryList = getattr(kbmod_search, "TrajectoryList", None)

    # API v2: pas d'ImageStack, utilisation de fill_psi_phi + search_cpu_only(psi_phi, params, out1, out2)
    use_api_v2 = (
        StackSearch is not None
        and fill_psi_phi is not None
        and search_cpu_only is not None
        and TrajectoryList is not None
    )
    use_gpu = args.gpu and search_gpu is not None and callable(search_gpu)

    ImageStack = getattr(kbmod_search, "ImageStack", None)
    LayeredImage = getattr(kbmod_search, "LayeredImage", None)
    RawImage = getattr(kbmod_search, "RawImage", None)
    if not use_api_v2:
        if ImageStack is None or StackSearch is None:
            avail = [x for x in dir(kbmod_search) if not x.startswith("_")]
            print("ERROR: KBMOD API: ImageStack absent et API v2 (fill_psi_phi_array_from_image_arrays, PsiPhiArray) incomplete.", file=sys.stderr)
            print("Symboles kbmod.search: %s" % avail, file=sys.stderr)
            sys.exit(1)

    # Liste des FITS (tri par nom)
    fits_paths = sorted(fits_dir.glob("*.fits")) + sorted(fits_dir.glob("*.fit"))
    if not fits_paths:
        print("ERROR: No FITS files in directory", file=sys.stderr)
        sys.exit(1)

    science_list = []
    mask_list = []
    variance_list = []
    times_used = []
    wcs_ref = None

    for fpath in fits_paths:
        try:
            with fits.open(str(fpath)) as hdul:
                data = hdul[0].data.astype(np.float32)
                if data is None or not data.size:
                    continue
                header = hdul[0].header
                if wcs_ref is None:
                    wcs_ref = WCS(header)
                jd = header.get("JD-UTC", 0.0) or header.get("JD", 0.0)
                if jd == 0.0:
                    date_obs = header.get("DATE-OBS", header.get("DATE"))
                    if date_obs:
                        exptime = float(header.get("EXPTIME", header.get("EXPOSURE", 0.0)))
                        t = Time(date_obs, scale="utc")
                        jd = float(t.jd) + exptime / (2.0 * 86400.0)
                if jd == 0.0:
                    continue
                gain = float(header.get("GAIN", header.get("EGAIN", 1.0)))
                rn = float(header.get("RDNOISE", header.get("READNOIS", 5.0)))
                var = np.abs(data) / np.clip(gain, 1e-6, None) + rn * rn
                var = np.clip(var, 1e-6, None).astype(np.float32)
                mask = np.zeros_like(data, dtype=np.float32)
                science_list.append(data)
                mask_list.append(mask)
                variance_list.append(var)
                times_used.append(jd)
        except Exception as e:
            print(f"Skip {fpath}: {e}", file=sys.stderr)
            continue

    if len(science_list) == 0:
        print("ERROR: No valid images for KBMOD", file=sys.stderr)
        sys.exit(1)

    # Masquage statique robuste: réduit les faux positifs sur sources fixes/artéfacts.
    # Pour objets très lents, il peut supprimer le signal utile -> option de désactivation.
    if not args.no_static_mask:
        try:
            cube = np.stack(science_list, axis=0)
            med = np.median(cube, axis=0)
            mad = np.median(np.abs(cube - med[None, :, :]), axis=0)
            robust_sigma = 1.4826 * mad + 1e-6
            thresh = med + max(1.0, float(args.static_mask_sigma)) * robust_sigma
            frac = np.mean(cube > thresh[None, :, :], axis=0)
            persistent = frac >= float(np.clip(args.static_mask_fraction, 0.05, 0.95))
            # Petite dilatation 1 pixel (8-connexité) sans dépendance externe.
            persistent_dilated = persistent.copy()
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    if dx == 0 and dy == 0:
                        continue
                    persistent_dilated |= np.roll(np.roll(persistent, dy, axis=0), dx, axis=1)
            static_mask = persistent_dilated.astype(np.float32)
            mask_list = [static_mask.copy() for _ in mask_list]
            static_frac = 100.0 * float(np.mean(static_mask > 0))
            print("KBMOD: static mask active (%.2f%% pixels masked)" % static_frac, file=sys.stderr)
        except Exception as e:
            print("KBMOD: static mask skipped (%s)" % e, file=sys.stderr)
    else:
        print("KBMOD: static mask disabled (--no-static-mask)", file=sys.stderr)

    scale = max(1e-6, args.scale)
    v_min_px_day = max(0.01, args.vmin / scale * 60.0 * 24.0)
    v_max_px_day = args.vmax / scale * 60.0 * 24.0
    jd_ref = times_used[0]

    raw_results = []
    if use_api_v2:
        # API v2: KBMOD 2.x utilise PsiPhiArray() sans arguments; fill_psi_phi peut creer ou remplir l'array
        try:
            PsiPhiArray = getattr(kbmod_search, "PsiPhiArray", None)
            if PsiPhiArray is None:
                raise RuntimeError("PsiPhiArray not found in kbmod.search")
            ny, nx = science_list[0].shape
            n_im = len(science_list)
            psi_phi = None
            # Essai 1: fill_psi_phi cree et retourne le PsiPhiArray (signature: n, science, mask, var, times)
            try:
                psi_phi = fill_psi_phi(n_im, science_list, mask_list, variance_list, times_used)
            except TypeError:
                pass
            # Essai 1b: fill_psi_phi(science, mask, var, times) -> PsiPhiArray
            if psi_phi is None:
                try:
                    psi_phi = fill_psi_phi(science_list, mask_list, variance_list, times_used)
                except TypeError:
                    pass
            # Essai 1c: fill_psi_phi(n_im, ny, nx, science, mask, var, times) -> PsiPhiArray
            if psi_phi is None:
                try:
                    psi_phi = fill_psi_phi(n_im, ny, nx, science_list, mask_list, variance_list, times_used)
                except TypeError:
                    pass
            # Essai 2: fill_psi_phi(psi_phi, n, ...) avec PsiPhiArray() vide + resize ou setters
            if psi_phi is None:
                psi_phi = PsiPhiArray()
                resized = False
                for method_name in ("resize", "set_dimensions", "resize_from_dims", "init"):
                    meth = getattr(psi_phi, method_name, None)
                    if meth is not None and callable(meth):
                        try:
                            import inspect
                            sig = inspect.signature(meth)
                            nparams = len([p for p in sig.parameters if p != "self"])
                            if nparams >= 3:
                                meth(n_im, ny, nx)
                            elif nparams == 2:
                                meth(ny, nx)
                            else:
                                meth(n_im)
                            resized = True
                            break
                        except (TypeError, ValueError):
                            pass
                if not resized:
                    for a, b, c in [("set_num_times", "set_height", "set_width"), ("set_num_images", "set_height", "set_width")]:
                        sa, sb, sc = getattr(psi_phi, a, None), getattr(psi_phi, b, None), getattr(psi_phi, c, None)
                        if callable(sa) and callable(sb) and callable(sc):
                            try:
                                sa(n_im)
                                sb(ny)
                                sc(nx)
                                resized = True
                                break
                            except (TypeError, ValueError):
                                pass
                # API KBMOD: (psi_phi, num_bytes, science_list, mask_list, variance_list, times_used) — pas de n_im
                fill_psi_phi(psi_phi, -1, science_list, mask_list, variance_list, times_used)
            if psi_phi is None:
                raise RuntimeError("Could not create or fill PsiPhiArray")
            # Lancer la recherche (GPU si --gpu et search disponible, sinon CPU)
            search_fn = search_gpu if use_gpu else search_cpu_only
            if use_gpu:
                print("KBMOD: using GPU search", file=sys.stderr)
            if SearchParameters is not None and TrajectoryList is not None:
                params = SearchParameters()
                # Cette version de SearchParameters n'a pas min_velocity/max_velocity (vitesse définie ailleurs)
                for attr, val in [
                    ("total_results", args.max_results),
                    ("results_per_pixel", 8),
                    ("min_lh", float(args.min_lh)),
                    ("min_observations", max(2, n_im // 2)),
                    ("x_start_min", 0),
                    ("x_start_max", nx - 1 if nx else 9999),
                    ("y_start_min", 0),
                    ("y_start_max", ny - 1 if ny else 9999),
                ]:
                    if hasattr(params, attr):
                        try:
                            setattr(params, attr, val)
                        except Exception:
                            pass
                print("KBMOD: v_min=%.2f v_max=%.1f \"/min -> %.0f .. %.0f px/jour (scale=%.4f \"/px)" % (args.vmin, args.vmax, v_min_px_day, v_max_px_day, scale), file=sys.stderr)
                results_traj = TrajectoryList(args.max_results)
                stamps_traj = TrajectoryList(args.max_results)
                # Appel standard (psi_phi, params, out1, out2) puis essai avec v_min, v_max en plus
                try:
                    search_cpu_only(psi_phi, params, results_traj, stamps_traj)
                except TypeError:
                    try:
                        search_cpu_only(psi_phi, params, results_traj, stamps_traj, v_min_px_day, v_max_px_day)
                    except TypeError:
                        pass
                except Exception as gpu_err:
                    if use_gpu and search_gpu is not None:
                        print("KBMOD GPU failed, falling back to CPU: %s" % gpu_err, file=sys.stderr)
                        try:
                            search_cpu_only(psi_phi, params, results_traj, stamps_traj)
                        except TypeError:
                            search_cpu_only(psi_phi, params, results_traj, stamps_traj, v_min_px_day, v_max_px_day)
                    else:
                        raise
                # Extraire les trajectoires : TrajectoryList (C++) n'est pas iterable en Python, utiliser .size() + .get(i) ou .at(i)
                raw_results = []
                try:
                    size_fn = getattr(results_traj, "size", None)
                    n = size_fn() if callable(size_fn) else (len(results_traj) if hasattr(results_traj, "__len__") else 0)
                    get_fn = getattr(results_traj, "get", None) or getattr(results_traj, "at", None)
                    if callable(get_fn):
                        for i in range(min(n, args.max_results)):
                            raw_results.append(get_fn(i))
                except (TypeError, AttributeError) as e:
                    print("KBMOD: could not read TrajectoryList (%s)" % e, file=sys.stderr)
                    raw_results = []
                if len(raw_results) == 0:
                    try:
                        attrs = [x for x in dir(params) if not x.startswith("_")]
                        print("KBMOD: 0 candidats. Attributs SearchParameters: %s" % (attrs,), file=sys.stderr)
                        syms = [s for s in dir(kbmod_search) if not s.startswith("_") and ("veloc" in s.lower() or "generator" in s.lower() or "config" in s.lower() or "search" in s.lower())]
                        if syms:
                            print("KBMOD: symboles kbmod.search (veloc/generator/config/search): %s" % syms[:20], file=sys.stderr)
                    except Exception:
                        pass
                    # Repli API classique (StackSearch avec v_min/v_max) car l'API v2 ne permet pas de fixer la vitesse
                    if ImageStack is not None and LayeredImage is not None and StackSearch is not None:
                        print("KBMOD: repli vers API classique (ImageStack + StackSearch) avec v=%.0f..%.0f px/jour" % (v_min_px_day, v_max_px_day), file=sys.stderr)
                        try:
                            psf = PSF(1.0)
                            kernel = getattr(psf, "kernel", np.ones((3, 3)) / 9.0)
                            if hasattr(psf, "set_kernel"):
                                psf.set_kernel(kernel)
                            imstack = ImageStack()
                            for i in range(len(science_list)):
                                data = science_list[i]
                                jd = times_used[i]
                                var = np.abs(data) / 1.0 + 25.0
                                var = np.clip(var, 1e-6, None).astype(np.float32)
                                mask = np.zeros_like(data, dtype=np.float32)
                                try:
                                    if RawImage is not None:
                                        sci = RawImage(data)
                                        msk = RawImage(mask)
                                        var_img = RawImage(var)
                                        limg = LayeredImage(sci, msk, var_img, jd, psf)
                                    else:
                                        limg = LayeredImage(data, mask, var, jd, psf)
                                    imstack.add_image(limg)
                                except Exception as e:
                                    print("KBMOD classic fallback LayeredImage: %s" % e, file=sys.stderr)
                                    continue
                            if imstack.img_count() > 0:
                                search = StackSearch(imstack)
                                if hasattr(search, "search"):
                                    search.search(v_min_px_day, v_max_px_day, -v_max_px_day, v_max_px_day, args.max_results)
                                raw_results = search.get_results(args.max_results) if hasattr(search, "get_results") else []
                                if raw_results:
                                    raw_results = list(raw_results)[: args.max_results] if not isinstance(raw_results, list) else raw_results[: args.max_results]
                                    print("KBMOD: API classique a retourné %d candidats" % len(raw_results), file=sys.stderr)
                        except Exception as e:
                            import traceback
                            print("KBMOD: échec repli API classique: %s" % e, file=sys.stderr)
                            traceback.print_exc(file=sys.stderr)
            elif SearchParameters is None:
                # Fallback: ancienne signature (psi_phi, v_min, v_max, max_results) si disponible
                try:
                    traj_list = search_cpu_only(psi_phi, v_min_px_day, v_max_px_day, args.max_results)
                    if traj_list is not None:
                        raw_results = list(traj_list)[: args.max_results]
                    else:
                        raw_results = []
                except TypeError:
                    raw_results = []
            else:
                raw_results = []
        except Exception as e:
            import traceback
            print("ERROR: KBMOD API v2: %s" % e, file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
            sys.exit(1)
    else:
        # API classique: ImageStack + LayeredImage + StackSearch
        psf = PSF(1.0)
        kernel = getattr(psf, "kernel", np.ones((3, 3)) / 9.0)
        if hasattr(psf, "set_kernel"):
            psf.set_kernel(kernel)
        imstack = ImageStack()
        LayeredImage = getattr(kbmod_search, "LayeredImage", None)
        RawImage = getattr(kbmod_search, "RawImage", None)
        for i in range(len(science_list)):
            data = science_list[i]
            jd = times_used[i]
            var = np.abs(data) / 1.0 + 25.0
            var = np.clip(var, 1e-6, None).astype(np.float32)
            mask = np.zeros_like(data, dtype=np.float32)
            try:
                if RawImage is not None:
                    sci = RawImage(data)
                    msk = RawImage(mask)
                    var_img = RawImage(var)
                    limg = LayeredImage(sci, msk, var_img, jd, psf)
                else:
                    limg = LayeredImage(data, mask, var, jd, psf)
                imstack.add_image(limg)
            except Exception as e:
                print("LayeredImage image %s: %s" % (i, e), file=sys.stderr)
                continue
        if imstack.img_count() == 0:
            print("ERROR: No valid images for KBMOD", file=sys.stderr)
            sys.exit(1)
        try:
            search = StackSearch(imstack)
            if hasattr(search, "search"):
                search.search(v_min_px_day, v_max_px_day, -v_max_px_day, v_max_px_day, args.max_results)
            raw_results = search.get_results(args.max_results) if hasattr(search, "get_results") else []
        except Exception as e:
            print("ERROR: StackSearch: %s" % e, file=sys.stderr)
            sys.exit(1)

    jd_ref = times_used[0] if times_used else 0.0
    all_results = []
    for r in raw_results:
        try:
            if hasattr(r, "x") and hasattr(r, "y"):
                x0, y0 = float(r.x), float(r.y)
            elif hasattr(r, "trajectory"):
                tr = r.trajectory
                x0 = float(getattr(tr, "x", 0))
                y0 = float(getattr(tr, "y", 0))
            else:
                x0, y0 = 0.0, 0.0
            if hasattr(r, "vx") and hasattr(r, "vy"):
                vx, vy = float(r.vx), float(r.vy)
            elif hasattr(r, "trajectory"):
                tr = r.trajectory
                vx = float(getattr(tr, "vx", 0))
                vy = float(getattr(tr, "vy", 0))
            else:
                vx = vy = 0.0
            lh = float(getattr(r, "lh", getattr(r, "likelihood", 0.0)))
            try:
                sky = wcs_ref.pixel_to_world(x0, y0)
                ra_deg = float(sky.ra.deg)
                dec_deg = float(sky.dec.deg)
            except Exception:
                ra_deg = dec_deg = float("nan")
            all_results.append({
                "ra_deg": ra_deg,
                "dec_deg": dec_deg,
                "x0": x0,
                "y0": y0,
                "vx_px_per_day": vx,
                "vy_px_per_day": vy,
                "likelihood": lh,
                "jd_ref": jd_ref,
            })
        except Exception as e:
            print(f"Skip result: {e}", file=sys.stderr)
            continue

    # Filtrage + déduplication (inspiré des étapes de post-traitement KBMOD publiées).
    min_lh = float(args.min_lh)
    candidates = []
    for r in all_results:
        lh = float(r.get("likelihood", 0.0))
        if not np.isfinite(lh) or lh < min_lh:
            continue
        x0 = float(r.get("x0", np.nan))
        y0 = float(r.get("y0", np.nan))
        vx = float(r.get("vx_px_per_day", np.nan))
        vy = float(r.get("vy_px_per_day", np.nan))
        if not (np.isfinite(x0) and np.isfinite(y0) and np.isfinite(vx) and np.isfinite(vy)):
            continue
        candidates.append(r)

    candidates.sort(key=lambda row: float(row.get("likelihood", 0.0)), reverse=True)
    dedup_pos = max(0.5, float(args.dedup_pos_px))
    dedup_vel = max(0.5, float(args.dedup_vel_pxday))
    kept = []
    for cand in candidates:
        x0 = float(cand["x0"])
        y0 = float(cand["y0"])
        vx = float(cand["vx_px_per_day"])
        vy = float(cand["vy_px_per_day"])
        duplicate = False
        for sel in kept:
            dx = x0 - float(sel["x0"])
            dy = y0 - float(sel["y0"])
            dvx = vx - float(sel["vx_px_per_day"])
            dvy = vy - float(sel["vy_px_per_day"])
            if (dx * dx + dy * dy) ** 0.5 <= dedup_pos and (dvx * dvx + dvy * dvy) ** 0.5 <= dedup_vel:
                duplicate = True
                break
        if not duplicate:
            kept.append(cand)
        if len(kept) >= args.max_results:
            break
    results = kept
    print(
        "KBMOD: raw=%d, after_lh=%d, after_dedup=%d (min_lh=%.2f)"
        % (len(all_results), len(candidates), len(results), min_lh),
        file=sys.stderr,
    )

    out_csv = fits_dir / "kbmod_candidates.csv"
    fieldnames = ["ra_deg", "dec_deg", "x0", "y0", "vx_px_per_day", "vy_px_per_day", "likelihood", "jd_ref"]
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    print(f"Wrote {len(results)} candidates to {out_csv}")
    if len(results) == 0:
        print("Aucun candidat trouvé. Vérifiez: plage vitesses (vmin/vmax \"/min), échelle (\"/px), nombre d'images.", file=sys.stderr)
        print("Pour un astéroïde de ceinture (ex. 17 Téthys): typiquement 20-60 \"/min → augmentez --vmax si besoin.", file=sys.stderr)


if __name__ == "__main__":
    main()
