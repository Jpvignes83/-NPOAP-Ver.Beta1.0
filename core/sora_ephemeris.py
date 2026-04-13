# core/sora_ephemeris.py
"""SORA : import, prédiction d’occultations (sora.prediction.prediction)."""

from __future__ import annotations

from typing import Any, Callable, Optional, Tuple, Union


def try_import_sora() -> Tuple[bool, str]:
    """
    Vérifie que le paquet sora-astro (requirements.txt) est importable.

    Retourne (True, "") ou (False, message d’erreur).
    """
    try:
        import sora  # noqa: F401
        from sora.prediction import prediction  # noqa: F401

        return True, ""
    except ImportError as e:
        return False, str(e)


def run_sora_import_check(callback: Callable[[bool, str], None]) -> None:
    """Pour thread UI : callback(ok, message)."""
    ok, msg = try_import_sora()
    callback(ok, msg)


def ucac4_vizier_catalogue() -> Any:
    """
    Catalogue UCAC4 (VizieR ``I/322A/out``) pour ``prediction(..., catalogue=...)``,
    aligné sur les désignations du type ``563-021605`` (colonne ``UCAC4``).

    Réf. : Zacharias+ 2012 ; voir doc SORA « Using different Vizier Catalogues ».
    """
    from astropy.time import Time
    from sora.star.catalog import VizierCatalogue

    return VizierCatalogue(
        name="UCAC4",
        cat_path="I/322A/out",
        code="UCAC4",
        ra="RAJ2000",
        dec="DEJ2000",
        pmra="pmRA",
        pmdec="pmDE",
        epoch=Time("J2000"),
        band={"V": "Vmag", "f": "f.mag", "J": "Jmag", "H": "Hmag", "K": "Kmag"},
        errors=["e_RAJ2000", "e_DEJ2000", "e_pmRA", "e_pmDE", None, None],
    )


def _resolve_prediction_catalogue(catalogue: Union[str, Any]) -> Any:
    """Convertit le libellé UI ``ucac4`` en objet ``VizierCatalogue`` SORA."""
    if isinstance(catalogue, str) and catalogue.strip().lower() == "ucac4":
        return ucac4_vizier_catalogue()
    return catalogue


def run_occultation_prediction(
    body: str,
    time_beg: str,
    time_end: str,
    *,
    mag_lim: Optional[Union[float, dict]] = None,
    step: float = 60.0,
    divs: int = 1,
    catalogue: Union[str, Any] = "gaiadr3",
    use_geocenter: bool = False,
    observer_name: str = "Site",
    lon_deg: float = 0.0,
    lat_deg: float = 0.0,
    height_m: float = 0.0,
    verbose: bool = False,
) -> Any:
    """
    Appelle ``sora.prediction.prediction`` (téléchargement éphémérides + étoiles Gaia ou Vizier).

    ``catalogue`` peut être ``gaiadr2`` / ``gaiaedr3`` / ``gaiadr3``, ou la chaîne ``ucac4``
    (catalogue UCAC4 via VizieR, pour rapprochement d’OccultWatcher).

    Parameters
    ----------
    body
        Nom ou désignation résolvable par SORA / Small Body Database (ex. ``\"136199\"``, ``\"Europa\"``).
    time_beg, time_end
        ISO ou tout format accepté par ``astropy.time.Time``.
    use_geocenter
        Si True, ``reference_center='geocenter'``. Sinon observateur topocentrique (lon, lat, h).
        Par défaut False (site réel).
    """
    from astropy import units as u
    from astropy.time import Time
    from sora.observer import Observer
    from sora.prediction import prediction

    body = (body or "").strip()
    if not body:
        raise ValueError("Indiquez un corps (nom ou numéro).")

    tb = Time(time_beg)
    te = Time(time_end)
    if te <= tb:
        raise ValueError("L’instant de fin doit être postérieur au début.")

    if use_geocenter:
        ref: Any = "geocenter"
    else:
        ref = Observer(
            name=observer_name.strip() or "Site",
            lon=lon_deg * u.deg,
            lat=lat_deg * u.deg,
            height=height_m * u.m,
        )

    cat_obj = _resolve_prediction_catalogue(catalogue)

    return prediction(
        tb,
        te,
        body=body,
        mag_lim=mag_lim,
        catalogue=cat_obj,
        step=float(step),
        divs=int(divs),
        verbose=verbose,
        reference_center=ref,
    )


def format_prediction_table(table: Any, max_lines: int = 80) -> str:
    """Texte lisible pour un ``PredictionTable`` SORA / table Astropy."""
    if table is None:
        return ""
    try:
        n = len(table)
    except TypeError:
        return str(table)
    if n == 0:
        return "Aucune occultation stellaire trouvée sur l’intervalle et les paramètres donnés."
    try:
        lines = table.pformat_all(max_width=120, max_lines=max_lines)
        return "\n".join(lines)
    except Exception:
        return str(table)
