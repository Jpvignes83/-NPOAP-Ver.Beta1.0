import logging
import os
from logging import FileHandler, Formatter
from tkinter import Text
from pathlib import Path

class TextHandler(logging.Handler):
    def __init__(self, text_widget: Text):
        super().__init__()
        self.text_widget = text_widget
        self.setFormatter(Formatter("%(asctime)s - %(levelname)s - %(message)s"))

    def emit(self, record):
        msg = self.format(record)
        def append():
            self.text_widget.configure(state="normal")
            self.text_widget.insert("end", msg + "\n")
            self.text_widget.see("end")
            self.text_widget.configure(state="disabled")
        self.text_widget.after(0, append)

def setup_logging(log_dir: str | Path = "logs", level=logging.INFO):
    log_dir = Path(log_dir)
    if not log_dir.is_absolute():
        base_dir = Path(__file__).resolve().parent.parent
        log_dir = base_dir / log_dir
    log_dir.mkdir(exist_ok=True)
    log_file = log_dir / f"npoap_{os.getpid()}.log"

    # Éviter doublons si appelé plusieurs fois
    root = logging.getLogger()
    if root.handlers:
        root.handlers.clear()

    # Fichier + console
    file_handler = FileHandler(log_file, encoding="utf-8")
    file_handler.setFormatter(Formatter("%(asctime)s - %(levelname)s - %(name)s - %(message)s"))

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(Formatter("%(levelname)s - %(message)s"))

    root.setLevel(level)
    root.addHandler(file_handler)
    root.addHandler(console_handler)

    logging.info(f"Logging initialisé → {log_file}")