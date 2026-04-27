#!/usr/bin/env python3
"""
Worker WSL pour NPOAP Prospector.

Ce worker valide l'environnement Prospector+FSPS côté WSL et renvoie un
résultat JSON exploitable par l'UI NPOAP.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: prospector_wsl_worker.py <input_json> <output_json>", file=sys.stderr)
        return 2

    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])

    payload = json.loads(input_path.read_text(encoding="utf-8"))
    sed_data = payload.get("sed_data", {}) or {}
    _ = payload.get("n_walkers", 100)
    _ = payload.get("n_steps", 1000)

    try:
        import numpy as np
        import prospect
        import fsps
    except Exception as e:
        output = {
            "success": False,
            "message": f"Backend WSL indisponible: {e}",
            "parameters": {},
            "uncertainties": {},
            "backend": "wsl",
        }
        output_path.write_text(json.dumps(output), encoding="utf-8")
        return 0

    # Validation minimale FSPS réelle côté WSL.
    try:
        sps = fsps.StellarPopulation()
        wave, spec = sps.get_spectrum(tage=5.0, peraa=True)
        fsps_points = int(len(wave))
    except Exception as e:
        output = {
            "success": False,
            "message": f"FSPS WSL non exploitable: {e}",
            "parameters": {},
            "uncertainties": {},
            "backend": "wsl",
        }
        output_path.write_text(json.dumps(output), encoding="utf-8")
        return 0

    wavelengths = np.asarray(sed_data.get("wavelength", []), dtype=float)
    fluxes = np.asarray(sed_data.get("flux", []), dtype=float)

    # Résumé rapide et explicite (pas d'échantillonnage MCMC ici).
    msg = (
        "Mode Prospector via WSL actif. "
        "Prospector+FSPS sont validés côté Linux "
        f"(prospect={getattr(prospect, '__version__', '?')}, fsps={getattr(fsps, '__version__', '?')}, "
        f"points_fsps={fsps_points}). "
        "Ce mode exécute actuellement une validation backend et un résumé SED, "
        "sans fit MCMC complet."
    )

    params = {}
    if wavelengths.size and fluxes.size and wavelengths.size == fluxes.size:
        idx = int(np.nanargmax(fluxes))
        params["lambda_peak_A"] = float(wavelengths[idx])
        params["flux_peak"] = float(fluxes[idx])
        params["n_points_sed"] = int(wavelengths.size)

    output = {
        "success": True,
        "message": msg,
        "parameters": params,
        "uncertainties": {},
        "backend": "wsl",
    }
    output_path.write_text(json.dumps(output), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
