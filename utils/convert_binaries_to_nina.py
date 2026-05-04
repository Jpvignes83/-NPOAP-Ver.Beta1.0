import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, Any


def ra_deg_to_hms(ra_deg: float):
    """Convertit une ascension droite en degrés vers (heures, minutes, secondes)."""
    ra_hours_total = ra_deg / 15.0
    h = int(ra_hours_total)
    m_total = (ra_hours_total - h) * 60.0
    m = int(m_total)
    s = (m_total - m) * 60.0
    return h, m, s


def dec_deg_to_dms(dec_deg: float):
    """Convertit une déclinaison en degrés vers (signe, degrés, minutes, secondes)."""
    negative = dec_deg < 0
    abs_deg = abs(dec_deg)
    d = int(abs_deg)
    m_total = (abs_deg - d) * 60.0
    m = int(m_total)
    s = (m_total - m) * 60.0
    return negative, d, m, s


def load_template(template_path: Path) -> Dict[str, Any]:
    with template_path.open("r", encoding="utf-8") as f:
        return json.load(f)


def update_template_from_row(template: Dict[str, Any], row: Dict[str, str], index: int) -> Dict[str, Any]:
    """
    Remplit un template NINA à partir d'une ligne du CSV laurent_all_pairs.csv.

    Colonnes utilisées :
      - ra1, dec1 : coordonnées de la première étoile (centre de la cible)
      - source_id_1, source_id_2 : pour construire un nom lisible si disponibles
    """
    ra_deg = float(row["ra1"])
    dec_deg = float(row["dec1"])
    ra_h, ra_m, ra_s = ra_deg_to_hms(ra_deg)
    neg_dec, dec_d, dec_m, dec_s = dec_deg_to_dms(dec_deg)

    # Construire un nom de cible
    src1 = row.get("source_id_1", "").strip()
    src2 = row.get("source_id_2", "").strip()
    if src1 and src2:
        target_name = f"Laurent_{src1}_{src2}"
    else:
        target_name = f"Laurent_{ra_deg:.5f}_{dec_deg:.5f}"

    # Cloner le template pour éviter les effets de bord
    data = json.loads(json.dumps(template))

    # Accès au bloc "Target" / "InputCoordinates" dans le template fourni
    target = data.get("Target", {})
    coords = target.get("InputCoordinates", {})

    target["TargetName"] = target_name
    data["Name"] = target_name

    coords["RAHours"] = ra_h
    coords["RAMinutes"] = ra_m
    coords["RASeconds"] = ra_s

    coords["NegativeDec"] = bool(neg_dec)
    coords["DecDegrees"] = dec_d
    coords["DecMinutes"] = dec_m
    coords["DecSeconds"] = dec_s

    target["InputCoordinates"] = coords
    data["Target"] = target

    return data


def convert_csv_to_nina(csv_path: Path, output_dir: Path, template_path: Path) -> int:
    """
    Convertit un fichier CSV laurent_all_pairs.csv en plusieurs fichiers JSON NINA.

    Retourne le nombre de fichiers créés.
    """
    template = load_template(template_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    count = 0
    with csv_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for idx, row in enumerate(reader, start=1):
            try:
                nina_obj = update_template_from_row(template, row, idx)
            except Exception as e:
                # En cas de problème sur une ligne, on la saute mais on continue
                print(f"[WARN] Ligne {idx} ignorée ({e})")
                continue

            name = nina_obj.get("Name") or f"Laurent_{idx}"
            # Nettoyer un peu le nom pour l'utiliser comme nom de fichier
            safe_name = "".join(c if c.isalnum() or c in ("_", "-", " ") else "_" for c in name).strip()
            if not safe_name:
                safe_name = f"Laurent_{idx}"

            out_file = output_dir / f"{safe_name}.json"
            with out_file.open("w", encoding="utf-8") as out_f:
                json.dump(nina_obj, out_f, ensure_ascii=False, indent=2)
            count += 1

    return count


def main():
    parser = argparse.ArgumentParser(
        description="Convertir un CSV laurent_all_pairs.csv en fichiers JSON NINA."
    )
    parser.add_argument("csv_file", type=str, help="Chemin vers le fichier CSV laurent_all_pairs.csv")
    parser.add_argument("output_dir", type=str, help="Répertoire de sortie pour les fichiers JSON NINA")
    parser.add_argument(
        "--template",
        type=str,
        default=str(Path(__file__).parent.parent / "templates" / "nina_target_template.json"),
        help="Chemin du template JSON NINA à utiliser",
    )
    args = parser.parse_args()

    csv_path = Path(args.csv_file)
    output_dir = Path(args.output_dir)
    template_path = Path(args.template)

    if not csv_path.exists():
        print(f"Erreur : fichier CSV introuvable : {csv_path}")
        raise SystemExit(1)
    if not template_path.exists():
        print(f"Erreur : template JSON NINA introuvable : {template_path}")
        raise SystemExit(1)

    count = convert_csv_to_nina(csv_path, output_dir, template_path)
    print(f"Conversion terminée : {count} fichiers JSON créés dans {output_dir}")


if __name__ == "__main__":
    main()

