#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Télécharge le catalogue GLADE+ (ASCII) et le convertit en HDF5 (Gladeplus.h5).

Source catalogue : http://elysium.elte.hu/~dalyag/GLADE+.txt (cf. https://glade.elte.hu/)
Conversion : ``scripts/vendor/convert_glade_catalog.py`` (adapté du convertisseur communautaire GLADE).

Sortie par défaut : ``%USERPROFILE%\\.npoap\\catalogues\\`` avec ``GLADE+.txt`` et ``Gladeplus.h5``.

Prérequis conversion (PyTables, hors requirements principaux) :

  pip install -r requirements_glade_build.txt

Usage :
  python scripts/download_and_convert_glade_catalog.py
  python scripts/download_and_convert_glade_catalog.py --skip-download
  python scripts/download_and_convert_glade_catalog.py --output-dir D:\\Donnees\\catalogues
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path

GLADE_PLUS_URL = "http://elysium.elte.hu/~dalyag/GLADE+.txt"


def _default_output_dir() -> Path:
    return Path.home() / ".npoap" / "catalogues"


def _vendor_converter_path() -> Path:
    here = Path(__file__).resolve().parent
    return here / "vendor" / "convert_glade_catalog.py"


def download_glade_txt(dest_txt: Path, chunk_bytes: int = 8 * 1024 * 1024) -> None:
    """Télécharge GLADE+.txt par blocs (fichier ~6 Go)."""
    import urllib.request

    dest_txt.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Téléchargement depuis %s vers %s …", GLADE_PLUS_URL, dest_txt)
    req = urllib.request.Request(
        GLADE_PLUS_URL,
        headers={"User-Agent": "NPOAP/GLADE-downloader"},
        method="GET",
    )
    total = 0
    milestone = 50 * 1024 * 1024
    next_report = milestone
    with urllib.request.urlopen(req, timeout=120) as resp, open(dest_txt, "wb") as out:
        while True:
            buf = resp.read(chunk_bytes)
            if not buf:
                break
            out.write(buf)
            total += len(buf)
            if total >= next_report:
                logging.info("  … %d Mo téléchargés", total // (1024 * 1024))
                next_report += milestone
    logging.info("Téléchargement terminé : %d octets (~%.2f Go)", total, total / (1024**3))


def run_conversion(
    py_exe: Path,
    converter: Path,
    txt_path: Path,
    h5_path: Path,
    max_luminosity_distance_mpc: float,
) -> None:
    cmd = [
        str(py_exe),
        str(converter),
        "-i",
        str(txt_path),
        "-o",
        str(h5_path),
        "-md",
        str(max_luminosity_distance_mpc),
    ]
    logging.info("Conversion HDF5 : %s", " ".join(cmd))
    r = subprocess.run(cmd)
    if r.returncode != 0:
        raise RuntimeError("Conversion GLADE → HDF5 a échoué (code %s)" % r.returncode)


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    p = argparse.ArgumentParser(
        description="Télécharge GLADE+ et produit Gladeplus.h5 pour NPOAP (tuiles alertes / NINA)."
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=_default_output_dir(),
        help="Répertoire de sortie (défaut : ~/.npoap/catalogues)",
    )
    p.add_argument(
        "--txt-name",
        default="GLADE+.txt",
        help="Nom du fichier ASCII téléchargé",
    )
    p.add_argument(
        "--h5-name",
        default="Gladeplus.h5",
        help="Nom du fichier HDF5 produit",
    )
    p.add_argument(
        "--max-luminosity-distance",
        type=float,
        default=500.0,
        help="Coupure distance luminosité (Mpc), défaut 500",
    )
    p.add_argument(
        "--skip-download",
        action="store_true",
        help="Ne pas télécharger ; utilise --txt-name déjà présent dans output-dir",
    )
    p.add_argument(
        "--skip-convert",
        action="store_true",
        help="Télécharger seulement, ne pas lancer la conversion",
    )
    args = p.parse_args(argv)

    out_dir: Path = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    txt_path = out_dir / args.txt_name
    h5_path = out_dir / args.h5_name
    converter = _vendor_converter_path()
    if not converter.is_file():
        logging.error("Script de conversion introuvable : %s", converter)
        return 1

    if not args.skip_download:
        if txt_path.exists():
            logging.info("Fichier existant (réutilisation) : %s", txt_path)
        else:
            try:
                download_glade_txt(txt_path)
            except Exception as e:
                logging.error("Échec du téléchargement : %s", e, exc_info=True)
                logging.error(
                    "Téléchargez manuellement depuis %s et placez le fichier en :\n  %s",
                    GLADE_PLUS_URL,
                    txt_path,
                )
                return 1
    else:
        if not txt_path.is_file():
            logging.error("Fichier manquant (--skip-download) : %s", txt_path)
            return 1

    if args.skip_convert:
        logging.info("Conversion ignorée (--skip-convert). Texte : %s", txt_path)
        return 0

    py_exe = Path(sys.executable)
    try:
        run_conversion(
            py_exe,
            converter,
            txt_path,
            h5_path,
            args.max_luminosity_distance,
        )
    except Exception as e:
        logging.error("%s", e, exc_info=True)
        return 1

    logging.info("Terminé. HDF5 : %s", h5_path)
    logging.info("Référence config : config.GLADE_CATALOG_H5_PATH")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
