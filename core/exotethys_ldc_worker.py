"""
Processus auxiliaire : ExoTETHyS appelle ``exit()`` en cas d'erreur.
Exécution en sous-processus pour ne pas terminer l'application NPOAP.

Usage : python -m core.exotethys_ldc_worker <chemin_cfg.txt>
"""
from __future__ import annotations

import sys


def main() -> None:
    if len(sys.argv) < 2:
        sys.stderr.write("usage: exotethys_ldc_worker <config.txt>\n")
        sys.exit(2)
    from exotethys.sail import ldc_calculate

    ldc_calculate(sys.argv[1])


if __name__ == "__main__":
    main()
