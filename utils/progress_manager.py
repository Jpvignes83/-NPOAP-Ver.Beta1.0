# utils/progress_manager.py
import tkinter.ttk as ttk
from tkinter import Tk

class ProgressManager(ttk.Progressbar):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, mode="determinate", **kwargs)
        self.max_value = 0
        self.root = master if isinstance(master, Tk) else master.winfo_toplevel()

    def start(self, max_value: int):
        self.max_value = max_value
        self["maximum"] = max_value
        self["value"] = 0
        self.update_idletasks()

    def step(self, increment: int = 1):
        if self.max_value > 0:
            new_value = min(self["value"] + increment, self.max_value)
            self["value"] = new_value
            self.update_idletasks()

    def update_progress(self, value: int):
        if self.max_value > 0:
            self["value"] = min(value, self.max_value)
            self.update_idletasks()

    def finish(self):
        self["value"] = self.max_value
        self.update_idletasks()

    def start_indeterminate(self):
        self.configure(mode="indeterminate")
        self.start()

    def stop_indeterminate(self):
        self.stop()
        self.configure(mode="determinate")