# Avant tout import chargeant CuPy : certaines entrées site-packages ont metadata=None,
# ce qui fait planter cupy._detect_duplicate_installation (AttributeError sur .get).
import importlib.metadata as _importlib_metadata

_orig_distributions = _importlib_metadata.distributions


def _distributions_skip_null_metadata():
    for dist in _orig_distributions():
        if getattr(dist, "metadata", None) is not None:
            yield dist


_importlib_metadata.distributions = _distributions_skip_null_metadata  # type: ignore[assignment]

import tkinter as tk
from gui.main_window import MainWindow
import utils.logging_handler as logging_handler
from utils.package_check import check_pylightcurve_update
import logging


def _install_global_toplevel_front_policy(root: tk.Tk) -> None:
    """
    Force toutes les futures fenêtres Toplevel à rester au premier plan.
    S'applique à l'ensemble des onglets NPOAP.
    """
    if getattr(tk.Toplevel, "_npoap_front_policy_installed", False):
        return

    original_init = tk.Toplevel.__init__

    def patched_toplevel_init(self, *args, **kwargs):
        original_init(self, *args, **kwargs)
        try:
            # Rattacher à la fenêtre principale pour éviter qu'une Toplevel
            # se perde derrière l'application.
            self.transient(root)
        except Exception:
            pass
        try:
            self.wm_attributes("-topmost", 1)
            self.lift()
        except Exception:
            pass
        try:
            self.after_idle(self.focus_force)
        except Exception:
            pass

    tk.Toplevel.__init__ = patched_toplevel_init
    tk.Toplevel._npoap_front_policy_installed = True


if __name__ == "__main__":
    logging_handler.setup_logging()
    logger = logging.getLogger(__name__)
    
    # Vérifier la mise à jour de pylightcurve au démarrage
    try:
        update_available, installed_ver, latest_ver = check_pylightcurve_update()
        if update_available:
            logger.info(
                f"⚠ Mise à jour de pylightcurve disponible: "
                f"version {installed_ver} → {latest_ver}. "
                f"Pour mettre à jour: pip install --upgrade pylightcurve"
            )
        elif latest_ver:
            logger.debug(f"pylightcurve est à jour (version {installed_ver})")
    except Exception as e:
        logger.debug(f"Erreur lors de la vérification de mise à jour de pylightcurve: {e}")
    
    root = tk.Tk()
    _install_global_toplevel_front_policy(root)
    app = MainWindow(root)
    root.mainloop()
    
