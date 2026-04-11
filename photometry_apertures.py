"""
Rayons d'ouverture photométriques (astéroïdes) — source unique pour SET-UP et batch.

Formule alignée sur l'onglet photométrie astéroïdes : r_ap ≈ 2 × FWHM (avec plafonds),
anneau de ciel cohérent avec les mêmes coefficients que le GUI.
"""
from __future__ import annotations

import logging
from typing import Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

DEFAULT_FWHM_PIX = 4.0

# Plage FWHM acceptée en entrée pour le batch (voiring plus large que le GUI interactif)
BATCH_MAX_INPUT_FWHM_PIX = 14.0


def aperture_radii_from_fwhm_pixels(
    fwhm: Optional[float],
    *,
    default_fwhm: float = DEFAULT_FWHM_PIX,
    max_input_fwhm: float = 8.0,
) -> Tuple[float, float, float]:
    """
    Calcule (r_ap, r_in, r_out) à partir d'une estimation de FWHM en pixels.

    Parameters
    ----------
    fwhm :
        FWHM estimé. Si None, non fini, ou hors [2, max_input_fwhm], on utilise default_fwhm.
    default_fwhm :
        Valeur de repli (typ. 4 px).
    max_input_fwhm :
        Borne haute inclusive pour faire confiance au FWHM (8 = comportement GUI historique ;
        14 = batch après mesure PSF sur comps).
    """
    if fwhm is None or not np.isfinite(fwhm) or fwhm < 2.0 or fwhm > max_input_fwhm:
        logger.debug(
            "FWHM hors plage ou invalide (%s), utilisation défaut %.1f px",
            fwhm,
            default_fwhm,
        )
        fwhm = float(default_fwhm)

    r_ap = max(2.0, min(round(2.0 * fwhm, 1), 20.0))
    r_in = max(4.0, min(round(2.5 * fwhm, 1), 40.0))
    r_out = max(6.0, min(round(3.5 * fwhm, 1), 50.0))

    if r_in <= r_ap:
        r_in = r_ap + 2.0
    if r_out <= r_in:
        r_out = r_in + 2.0

    return float(r_ap), float(r_in), float(r_out)


def stabilize_fwhm_sequential(
    raw: Optional[float],
    last_valid: Optional[float],
    *,
    max_step_ratio: float = 1.35,
    blend_with_last: float = 0.35,
) -> Optional[float]:
    """
    Évite les sauts brutaux de FWHM d'une image à l'autre sans « figer » la valeur précédente.

    Si ``raw`` est dans [last/max_step_ratio, last*max_step_ratio], on garde ``raw``.
    Sinon on ramène d'abord ``raw`` sur le bord de la bande, puis on mélange avec ``last_valid`` :
    ``(1-blend)*borne + blend*last`` pour laisser tout de même évoluer le PSF dans le temps.

    Retourne None si raw est None ou non positif.
    """
    if raw is None or not np.isfinite(raw) or raw <= 0:
        return None
    raw = float(raw)
    if last_valid is None or not np.isfinite(last_valid) or last_valid <= 0:
        return raw
    last_valid = float(last_valid)
    lo = last_valid / max_step_ratio
    hi = last_valid * max_step_ratio
    if lo <= raw <= hi:
        return raw
    edge = float(np.clip(raw, lo, hi))
    blended = (1.0 - blend_with_last) * edge + blend_with_last * last_valid
    logger.warning(
        "[BATCH] FWHM PSF brut=%.2f px hors bande [%.2f, %.2f] (précédent=%.2f) -> lissage %.2f px",
        raw,
        lo,
        hi,
        last_valid,
        blended,
    )
    return blended
