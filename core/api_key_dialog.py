import os
import tkinter as tk
from tkinter import simpledialog, messagebox
from pathlib import Path

API_KEY_PATH = Path.home() / ".astrometry_api_key"


def ask_api_key():
    """Ouvre une boîte de dialogue pour saisir la clé API Astrometry.net."""
    root = tk.Tk()
    root.withdraw()  # Ne pas afficher la fenêtre principale
    api_key = simpledialog.askstring(
        "Clé API Astrometry.net",
        "Veuillez entrer votre clé API Astrometry.net :"
    )
    root.destroy()
    return api_key


def get_astrometry_api_key():
    """Récupère la clé API : fichier local sinon demande à l'utilisateur."""
    api_key = None

    # 1. Lire depuis le fichier local s'il existe
    if API_KEY_PATH.exists():
        try:
            with open(API_KEY_PATH, "r", encoding="utf-8") as f:
                api_key = f.read().strip()
                if api_key:
                    return api_key
        except Exception as e:
            print(f"Erreur lecture clé API locale : {e}")

    # 2. Sinon demander via boîte de dialogue
    api_key = ask_api_key()

    # 3. Sauvegarder la clé si saisie valide
    if api_key:
        try:
            with open(API_KEY_PATH, "w", encoding="utf-8") as f:
                f.write(api_key.strip())
            messagebox.showinfo(
                "Clé API enregistrée",
                f"La clé a été sauvegardée dans :\n{API_KEY_PATH}"
            )
        except Exception as e:
            messagebox.showerror("Erreur", f"Impossible d'enregistrer la clé : {e}")
    else:
        messagebox.showwarning("Clé manquante", "Aucune clé API saisie.")

    return api_key


class APIKeyDialog(simpledialog.Dialog):
    """
    Boîte de dialogue modale pour saisir/modifier la clé API Astrometry.net.

    Utilisation typique depuis la GUI :
        dialog = APIKeyDialog(parent)
        if dialog.api_key:
            # nouvelle clé disponible
    """

    def __init__(self, parent, title="Clé API Astrometry.net"):
        self.api_key = None
        super().__init__(parent, title)

    def body(self, master):
        tk.Label(master, text="Entrez votre clé API Astrometry.net :").grid(
            row=0, column=0, padx=10, pady=(10, 5), sticky="w"
        )

        self.entry = tk.Entry(master, show="*", width=40)
        self.entry.grid(row=1, column=0, padx=10, pady=(0, 10))

        # Pré-remplit avec la clé existante si disponible
        if API_KEY_PATH.exists():
            try:
                with open(API_KEY_PATH, "r", encoding="utf-8") as f:
                    key = f.read().strip()
                if key:
                    self.entry.insert(0, key)
            except Exception:
                pass

        return self.entry  # focus initial

    def apply(self):
        key = self.entry.get().strip()
        self.api_key = key

        if not key:
            messagebox.showwarning(
                "Clé API",
                "Aucune clé API saisie."
            )
            return

        # Sauvegarde dans le fichier API_KEY_PATH
        try:
            with open(API_KEY_PATH, "w", encoding="utf-8") as f:
                f.write(key)
        except Exception as e:
            messagebox.showerror(
                "Clé API",
                f"Impossible d'enregistrer la clé API : {e}"
            )
