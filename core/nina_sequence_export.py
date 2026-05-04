# core/nina_sequence_export.py
"""
Export JSON de conteneur cible NINA (DeepSkyObjectContainer), format attendu par NINA 2/3.
"""
from __future__ import annotations

from typing import Any, Dict


def ra_dec_to_nina_format(ra_deg: float, dec_deg: float) -> Dict[str, Any]:
    """RA/Dec en degrés → champs NINA InputCoordinates (heures pour RA, DMS pour Dec)."""
    ra_hours = ra_deg / 15.0
    ra_h = int(ra_hours)
    ra_minutes_float = (ra_hours - ra_h) * 60.0
    ra_m = int(ra_minutes_float)
    ra_s = (ra_minutes_float - ra_m) * 60.0

    dec_abs = abs(dec_deg)
    dec_negative = dec_deg < 0
    dec_d = int(dec_abs)
    dec_minutes_float = (dec_abs - dec_d) * 60.0
    dec_m = int(dec_minutes_float)
    dec_s = (dec_minutes_float - dec_m) * 60.0

    return {
        "RAHours": ra_h,
        "RAMinutes": ra_m,
        "RASeconds": ra_s,
        "NegativeDec": dec_negative,
        "DecDegrees": dec_d,
        "DecMinutes": dec_m,
        "DecSeconds": dec_s,
    }


def build_nina_deep_sky_container_dict(target_name: str, ra_deg: float, dec_deg: float) -> dict:
    """Structure identique à un export NINA (fichier .json de séquence cible)."""
    coords = ra_dec_to_nina_format(ra_deg, dec_deg)
    return {
        "$id": "1",
        "$type": "NINA.Sequencer.Container.DeepSkyObjectContainer, NINA.Sequencer",
        "Target": {
            "$id": "2",
            "$type": "NINA.Astrometry.InputTarget, NINA.Astrometry",
            "Expanded": True,
            "TargetName": target_name,
            "PositionAngle": 0.0,
            "InputCoordinates": {
                "$id": "3",
                "$type": "NINA.Astrometry.InputCoordinates, NINA.Astrometry",
                "RAHours": coords["RAHours"],
                "RAMinutes": coords["RAMinutes"],
                "RASeconds": coords["RASeconds"],
                "NegativeDec": coords["NegativeDec"],
                "DecDegrees": coords["DecDegrees"],
                "DecMinutes": coords["DecMinutes"],
                "DecSeconds": coords["DecSeconds"],
            },
        },
        "ExposureInfoListExpanded": False,
        "ExposureInfoList": {
            "$id": "4",
            "$type": "NINA.Core.Utility.AsyncObservableCollection`1[[NINA.Sequencer.Utility.ExposureInfo, NINA.Sequencer]], NINA.Core",
            "$values": [],
        },
        "Strategy": {
            "$type": "NINA.Sequencer.Container.ExecutionStrategy.SequentialStrategy, NINA.Sequencer"
        },
        "Name": target_name,
        "Conditions": {
            "$id": "5",
            "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Conditions.ISequenceCondition, NINA.Sequencer]], System.Collections.ObjectModel",
            "$values": [],
        },
        "IsExpanded": True,
        "Items": {
            "$id": "6",
            "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.SequenceItem.ISequenceItem, NINA.Sequencer]], System.Collections.ObjectModel",
            "$values": [],
        },
        "Triggers": {
            "$id": "7",
            "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Trigger.ISequenceTrigger, NINA.Sequencer]], System.Collections.ObjectModel",
            "$values": [],
        },
        "Parent": None,
        "ErrorBehavior": 0,
        "Attempts": 1,
    }


def sanitize_nina_filename_component(name: str, max_len: int = 80) -> str:
    """Nom de fichier sûr pour Windows."""
    s = "".join(c if c.isalnum() or c in "-_" else "_" for c in str(name))
    return s[:max_len] if s else "target"
