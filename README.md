# -NPOAP-Ver.Beta1.0
première version en BêtaTest
NPOAP — description du projet

NPOAP (Nouvelle Plateforme d’Observation et d’Analyse Photométrique)

C'est une application desktop sous Python et Tkinter, pensée pour enchaîner réduction d’images, photométrie et analyse sur des données d’amateur ou semi-professionnel. Elle s’appuie sur l’écosystème scientifique Python (Astropy, photutils, etc.), un environnement Conda astroenv (Python 3.11), et peut s’interfacer avec WSL pour des outils Linux (astrométrie locale, KBMOD, chaînes lourdes).

Résumé des fonctions (par onglet)

Accueil
Paramètres observatoire, clés / services (ex. Astrometry.net), calculateur d’échelle, point d’entrée configuration.
Observation de la nuit
Préparation de nuit : éphémérides, visibilité, lien avec le planétarium / export (ex. NINA) selon la version.
Planétarium (C2A)
Visualisation / cible : cartes du ciel et intégration avec la planification d’observation.
Réduction de données
Chaîne CCD : bias, darks, flats, lights, astrométrie (local ou en ligne), alignement WCS, empilements.
Photométrie exoplanètes
Workflow type HOPS : données → photométrie → ajustement de transits (outil intégré ou lancé depuis NPOAP).
Photométrie astéroïdes
Courbes de lumière et astrométrie d’astéroïdes ; prise en charge d’astrométrie type zero-aperture (CuPy optionnel).
Photométrie transitoires
Mesures sur événements transitoires (supernovae, etc.), avec appuis catalogues / TNS selon les modules.
Analyse des données
Périodogrammes, TTV, systèmes multiples, simulations N-body / ajustements avancés sur courbes de lumière.
Étoiles binaires
Modélisation et analyse de binaires à éclipses (ex. PHOEBE).
Easy Lucky Imaging
Chaîne simplifiée lucky imaging et mesure de séparation pour couples serrés.
Analyse d’amas
Diagrammes couleur-magnitude, âge / distance (ex. isochrones PARSEC via ezpadova si installé).
Spectroscopie
Chargement et analyse de spectres (specutils, etc.) ; installation optionnelle de Prospector (SED).
Catalogues
Extraction et usage de catalogues (Gaia, astéroïdes, exoplanètes, etc.) pour alimenter les autres onglets.
