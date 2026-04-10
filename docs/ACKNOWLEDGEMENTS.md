# Remerciements (Acknowledgments)

NPOAP utilise de nombreuses bibliothèques et outils open-source de la communauté scientifique. Nous tenons à remercier les développeurs et contributeurs suivants :

## Crédits NPOAP / HOPS-modified

- **Responsable de l'intégration HOPS-modified** : J.P Vignes
- **Contact** : jeanpascal.vignes@gmail.com

## Bibliothèques Python principales

### Calcul scientifique
- **NumPy** : Bibliothèque fondamentale pour le calcul numérique (https://numpy.org/)
- **SciPy** : Bibliothèque scientifique pour Python, incluant optimisation, statistiques, traitement du signal (https://scipy.org/)
- **Pandas** : Manipulation et analyse de données structurées (https://pandas.pydata.org/)

### Astronomie
- **Astropy** : Bibliothèque Python pour l'astronomie (https://www.astropy.org/)
  - Utilisée pour la manipulation des images FITS, les coordonnées célestes, les transformations WCS, et les calculs astronomiques
  - Développée par la communauté Astropy

- **photutils** : Outils de photométrie astronomique (https://photutils.readthedocs.io/)
  - Utilisée pour la détection d'étoiles, la photométrie par ouverture, et l'estimation du fond de ciel
  - Développée dans le cadre du projet Astropy

- **astroquery** : Interface Python pour les services astronomiques en ligne (https://astroquery.readthedocs.io/)
  - Utilisée pour interroger les catalogues Gaia (Vizier), VSX, SB9, DEBCat et les éphémérides JPL Horizons
  - Développée dans le cadre du projet Astropy

- **reproject** : Reprojection d'images astronomiques (https://reproject.readthedocs.io/)
  - Utilisée pour l’**alignement WCS** des images en réduction (onglet Réduction de Données) : reprojection de chaque image sur le WCS de référence pour supprimer les décalages avant empilement
  - Utilisée pour l'alignement et la réorientation d'images lors de la soustraction d'images dans l'onglet Photométrie Transitoires
  - Permet la mise à l'échelle, le repositionnement et la correction d'orientation (notamment pour les images de l'hémisphère sud)
  - Développée dans le cadre du projet Astropy

### Visualisation
- **Matplotlib** : Bibliothèque de visualisation (https://matplotlib.org/)
  - Utilisée pour l'affichage graphique des images, courbes de lumière, et périodogrammes

### Étoiles binaires
- **PHOEBE2** : Bibliothèque Python pour la modélisation d'étoiles binaires à éclipses (https://phoebe-project.org/)
  - Utilisée pour la modélisation et l'analyse des systèmes binaires dans l'onglet "Étoiles Binaires"
  - Développée par le projet PHOEBE

### Analyse de données
- **emcee** : Bibliothèque MCMC pour l'ajustement de modèles (https://emcee.readthedocs.io/)
  - Utilisée pour l'analyse des variations de temps de transit (TTV)

- **statsmodels** : Bibliothèque d'analyse statistique (https://www.statsmodels.org/)
  - Utilisée pour les analyses statistiques avancées des courbes de lumière

- **pylightcurve** : Bibliothèque Python pour la modélisation et l'analyse de courbes de lumière d'exoplanètes (https://github.com/ucl-exoplanets/pylightcurve)
  - Développée par l'équipe UCL Exoplanets (University College London)
  - Utilisée pour la modélisation des transits d'exoplanètes, le calcul des coefficients d'assombrissement du limbe, et l'ajustement flexible de courbes de lumière multi-époques
  - Licence MIT

- **Wōtan (wotan)** : détrendage robuste de séries temporelles pour la photométrie de transits et, plus généralement, les courbes de lumière (https://github.com/hippke/wotan)
  - Développé par **Michael Hippke** et collaborateurs ; interface Python unifiée pour de nombreux estimateurs (biweight, médiane, splines, etc.)
  - Utilisé dans NPOAP en **option** dans le viewer d’affinage du temps de transit central (Analyse de données → période → courbe avec marqueurs) : `flatten` avec **masque transit** dérivé du Tc courant et de la largeur T14 saisie, puis normalisation locale hors transit
  - Utilisé également dans l’onglet **Photométrie astéroïdes** (analyse de courbe de lumière, cadre **Détrendage**) : méthodes **Wōtan biweight** et **Wōtan median** sur toute la série, via le module `core/asteroid_lc_detrend.py` (fenêtre temporelle en jours, automatique ou saisie ; les options Wōtan n’apparaissent que si le paquet `wotan` est installé)
  - **Citation** : Hippke, David, Mulders & Heller (2019), *AJ*, 158, 143 ; DOI [10.3847/1538-3881/ab3984](https://doi.org/10.3847/1538-3881/ab3984) ; arXiv:1906.00966
  - Licence MIT

- **HOPS (HOlomon Photometric Software)** : logiciel de photométrie exoplanètes intégré dans l'onglet Photométrie Exoplanètes
  - Auteurs : **Angelos Tsiaras** et **Konstantinos Karpouzas**
  - Copyright (c) 2017 Angelos Tsiaras and Konstantinos Karpouzas
  - Licence MIT
  - Les écrans d’export vers **ExoClock**, **ETD** et **AAVSO** renvoient notamment vers la base **Exoplanet Transit Database (ETD)** hébergée sur **VarAstro** (https://var.astro.cz/fr/Exoplanets)

- **exoclock** : bibliothèque Python liée au projet **ExoClock** (éphémérides et catalogue d’exoplanètes en transit pour la mission ARIEL et le suivi au sol)
  - Utilisée dans la chaîne **HOPS** et, dans NPOAP, comme source complémentaire dans l’onglet **Observation de la nuit** (agrégation avec la NASA Exoplanet Archive)
  - Site du projet : https://www.exoclock.space/

### Simulation N-body
- **rebound** : Bibliothèque de simulation N-body pour systèmes planétaires (https://rebound.readthedocs.io/)
  - Utilisée pour les simulations gravitationnelles dans l'analyse de systèmes multiples (onglet Analyse de Données)
  - Permet de modéliser les interactions gravitationnelles entre planètes
  - **Citation** : Pour obtenir la citation exacte de rebound, utilisez :
    ```python
    import rebound
    sim = rebound.Simulation()
    sim.cite()
    ```
  - Référence principale : Rein & Spiegel (2015), "IAS15: a fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine precision over a billion orbits", *MNRAS*, 446, 1424

- **ultranest** : Bibliothèque d'échantillonnage bayésien avec nested sampling (https://johannesbuchner.github.io/UltraNest/)
  - Utilisée pour le fitting bayésien des modèles N-body aux observations TTV
  - Alternative robuste aux méthodes MCMC pour l'estimation de paramètres
  - Utilise l'algorithme MLFriends (Buchner, 2014; 2019) pour le nested sampling Monte Carlo
  - Permet de dériver les distributions de probabilité a posteriori et l'évidence bayésienne
  - Référence : Buchner (2021), "UltraNest - a robust, general purpose Bayesian inference engine"

### Photométrie de transitoires
- **STDPipe** : Simple Transient Detection Pipeline (https://stdpipe.readthedocs.io/)
  - Auteur principal : Sergey Karpov
  - Utilisée pour la photométrie des transitoires : astrométrie automatique, soustraction d'images, détection de transitoires, photométrie calibrée
  - Fournit des méthodes avancées de détection (segmentation, DAOStarFinder, IRAFStarFinder)
  - Permet le téléchargement d'images de référence depuis Pan-STARRS, SDSS, DES
  - Licence MIT

### Autres bibliothèques
- **Pillow (PIL)** : Bibliothèque de traitement d'images (https://python-pillow.org/)
  - Utilisée pour le traitement et la manipulation d'images

- **reportlab** : Bibliothèque de génération de PDF (https://www.reportlab.com/)
  - Utilisée pour la génération de documentation PDF

- **requests** : Bibliothèque HTTP pour Python (https://requests.readthedocs.io/)
  - Utilisée pour les requêtes HTTP vers les services en ligne

- **specutils** : Bibliothèque Python pour l'analyse de données spectroscopiques (https://specutils.readthedocs.io/)
  - Utilisée pour la représentation, le chargement, la manipulation et l'analyse de spectres d'étoiles
  - Développée dans le cadre du projet Astropy

- **synphot** : Bibliothèque de photométrie synthétique (https://synphot.readthedocs.io/)
  - Utilisée pour les calculs de photométrie synthétique avec filtres U, B, V, R, I
  - Développée dans le cadre du projet Astropy

### Analyse de populations stellaires et spectroscopie avancée

- **Prospector** : Code Python pour inférer les propriétés des populations stellaires à partir de données photométriques et spectroscopiques (https://github.com/bd-j/prospector)
  - Utilisé pour l'analyse de spectres de galaxies et l'inférence de propriétés stellaires (âge, métallicité, extinction, redshift)
  - Utilise l'inférence bayésienne avec des modèles de populations stellaires simples (SSP - Simple Stellar Population)
  - Développé par l'équipe BD-J (Benjamin Johnson et collaborateurs)
  - **Important** : Le package `prospector` sur PyPI est un autre outil. La version astronomique doit être installée depuis GitHub : `pip install git+https://github.com/bd-j/prospector.git`

- **sedpy** : Bibliothèque Python pour la manipulation de SED (Spectral Energy Distribution) (https://github.com/bd-j/sedpy)
  - Utilisée par Prospector pour créer et manipuler des SED
  - **Important** : Doit être installée depuis GitHub (pas PyPI) pour avoir le module `observate` requis par Prospector
  - Développée par l'équipe BD-J

- **FSPS (Flexible Stellar Population Synthesis)** : Code Fortran avec bindings Python pour générer des modèles de populations stellaires (https://github.com/cconroy20/fsps)
  - Utilisé par Prospector comme moteur pour générer les templates SSP
  - Nécessite CMake et gfortran pour la compilation
  - Peut être installé depuis GitHub : `git clone https://github.com/dfm/python-fsps.git` puis `pip install .`
  - Développé par Charlie Conroy et collaborateurs

## Services et catalogues externes

### Catalogues astrométriques et photométriques
- **Gaia DR3** : Catalogue astrométrique de l'Agence Spatiale Européenne (ESA)
  - Utilisé comme référence astrométrique et photométrique (https://www.cosmos.esa.int/web/gaia)
  - Fournit positions précises, magnitudes, couleurs pour des millions d'étoiles

### Catalogues d'étoiles binaires
- **VSX (AAVSO Variable Star Index)** : Catalogue des étoiles variables (https://www.aavso.org/vsx/)
  - Utilisé pour la recherche de paramètres de systèmes binaires variables

- **SB9 (Ninth Catalog of Spectroscopic Binary Orbits)** : Catalogue d'orbites de binaires spectroscopiques (https://sb9.astro.ulb.ac.be/)
  - Utilisé pour récupérer les paramètres orbitaux de systèmes binaires spectroscopiques

- **DEBCat** : Database of Eclipsing Binaries Catalog (http://www.astro.keele.ac.uk/jkt/debcat/)
  - Utilisé pour les paramètres détaillés de systèmes binaires à éclipses

- **Washington Double Star Catalog (WDS)** : Catalogue des étoiles doubles de l'US Naval Observatory (https://crf.usno.navy.mil/wds/)
  - Format de référence pour les rapports de mesure de séparation binaire
  - Utilisé pour générer les rapports au format standard WDS dans l'onglet Easy Lucky Imaging

- **Minor Planet Center (MPC)** : Centre de données pour les astéroïdes et comètes (https://www.minorplanetcenter.net/)
  - **NEA.txt** : Catalogue des astéroïdes géocroiseurs (Near-Earth Asteroids)
  - **AllCometEls.txt** : Catalogue des éléments orbitaux de toutes les comètes
  - Utilisés dans l'onglet "Catalogues" → "Astéroïdes & Comètes" pour la sélection et le tri d'objets
  - Utilisés dans l’onglet **Observation de la nuit** ; les éphémérides détaillées et la magnitude apparente affichée (**Mag vis**) proviennent du service en ligne MPC (**MPES**, `cgi.minorplanetcenter.net`)

### Services d'astrométrie
- **Astrometry.net / NOVA** : Service d'astrométrie automatique (http://astrometry.net/, https://nova.astrometry.net/)
  - Utilisé pour la résolution de champs stellaires (plate solving) via NOVA (service en ligne) et solve-field (installation locale)
  - Développé par l'équipe Astrometry.net

### Services d'éphémérides
- **JPL Horizons** : Service d'éphémérides du Jet Propulsion Laboratory (NASA)
  - Utilisé pour obtenir les éphémérides des astéroïdes et planètes (https://ssd.jpl.nasa.gov/horizons/)

### Services de catalogues
- **Vizier** : Service de catalogues astronomiques du Centre de données astronomiques de Strasbourg (CDS)
  - Utilisé pour interroger les catalogues Gaia et autres catalogues astronomiques (https://vizier.cds.unistra.fr/)

- **TESS EBS (Eclipsing Binary Stars)** : Catalogue des étoiles binaires à éclipses observées par TESS (https://tessebs.villanova.edu/)
  - Utilisé dans l'onglet "Catalogues" pour l'extraction de données d'étoiles binaires

- **Exoplanet.eu** : Base de données d'exoplanètes (https://exoplanet.eu/)
  - Utilisé dans l'onglet "Catalogues" pour l'extraction de données d'exoplanètes

### Images de référence
- **Pan-STARRS** : Panoramic Survey Telescope and Rapid Response System
  - Utilisé pour télécharger des images de référence profondes pour la soustraction d'images (https://panstarrs.stsci.edu/)

- **SDSS** : Sloan Digital Sky Survey
  - Utilisé pour télécharger des images de référence pour la soustraction d'images (https://www.sdss.org/)

- **DES** : Dark Energy Survey
  - Utilisé pour télécharger des images de référence pour les champs de l'hémisphère sud (https://www.darkenergysurvey.org/)

### Bases de données d'exoplanètes
- **NASA Exoplanet Archive** : Base de données d'exoplanètes (https://exoplanetarchive.ipac.caltech.edu/)
  - Utilisée pour les éphémérides et paramètres des exoplanètes (dont requêtes **TAP** sur la table `pscomppars` pour les planètes en transit dans l’onglet **Observation de la nuit**)

- **ExoClock** : projet communautaire de suivi des éphémérides d’exoplanètes en transit (https://www.exoclock.space/)
  - Catalogue et outils utilisés via le module Python **exoclock** dans NPOAP (observation de la nuit)

- **VarAstro — Exoplanet Transit Database (ETD)** : portail de la section « étoiles variables et exoplanètes » de la Société astronomique tchèque (https://var.astro.cz/fr/Exoplanets)
  - L’interface web et l’API catalogue sont des services tiers ; l’accès API typique requiert une authentification — NPOAP désactive la source ETD tant qu’aucune session valide n’est disponible


## Outils et frameworks

- **Python** : Langage de programmation (https://www.python.org/)
  - Langage principal de développement

- **Conda/Miniconda** : Système de gestion d'environnements et de paquets (https://docs.conda.io/)
  - Utilisé pour la gestion des environnements Python et des dépendances

- **Tkinter** : Bibliothèque d'interface graphique (incluse dans Python)
  - Utilisée pour l'interface utilisateur graphique

- **WSL (Windows Subsystem for Linux)** : Système pour exécuter Linux sur Windows (Microsoft)
  - Utilisé pour l'installation locale d'Astrometry.net sur Windows

## Techniques et méthodologies

- **Zero-Aperture-Astrometry** (Ben Sharkey) : Outil et méthodologie pour dériver des positions astrométriques zero-aperture à partir d’observations (https://github.com/bensharkey/Zero-Aperture-Astrometry)
  - La méthode zero-aperture de NPOAP (onglet Photométrie Astéroïdes) s’en inspire : fit linéaire pondéré des positions (ou des résidus moyens RA/Dec) en fonction de l’aperture photométrique, extrapolation à aperture = 0, et correction du WCS
  - Utilisé pour les observations au format ADES (PSV/XML) dans le projet original ; NPOAP adapte l’approche aux images FITS (plusieurs apertures en pixels, résidus moyens par aperture, puis extrapolation)
  - Licence GPL-3.0

- **REDUC** : Techniques de réduction d'images inspirées du tutoriel REDUC (http://www.astrosurf.com/hfosaf/reduc/tutoriel.htm)
  - **BestOf** : Tri et sélection des meilleures images selon leur qualité (FWHM, SNR, contraste)
  - **ELI (Easy Lucky Imaging)** : Empilement d'images avec alignement sub-pixel pour améliorer la résolution
  - **Centroiding** : Mesure précise des positions d'étoiles par centroïde
  - Utilisées dans l'onglet Easy Lucky Imaging pour le traitement d'images d'étoiles binaires

## Références bibliographiques

### Astrométrie
- **Farnocchia et al. (2022)** : "International Asteroid Warning Network Timing Campaign: 2019 XS", *Planetary Science Journal*, 3:156
  - Guide de référence pour les améliorations astrométriques implémentées dans NPOAP
  - DOI: https://doi.org/10.3847/PSJ/ac7224

### Étoiles binaires
- **Prša (2018)** : "Modeling and Analysis of Eclipsing Binary Stars: The Theory and Design Principles of PHOEBE"
  - Référence principale pour l'intégration de PHOEBE2
  - DOI: 10.1088/978-0-7503-1287-5

- **Laurent (2022)** : "La séparation linéaire utilisée pour déterminer le caractère physique d'une étoile double visuelle"
  - Article publié dans Étoiles Doubles - n°05 (Décembre 2022)
  - Auteur : Philippe Laurent (SAF - Commission des Étoiles Doubles, Président de l'Association Astronomie en Provence)
  - Méthode utilisée pour identifier les couples d'étoiles binaires physiques à partir des données Gaia DR3
  - Implémentée dans l'onglet "Catalogues" → "Étoiles Binaires" pour la création de catalogues de binaires physiques

### Simulation N-body
- **Rein & Spiegel (2015)** : "IAS15: a fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine precision over a billion orbits", *Monthly Notices of the Royal Astronomical Society*, 446, 1424
  - Intégrateur utilisé par rebound pour les simulations N-body
  - DOI: https://doi.org/10.1093/mnras/stu2254

- **Rein & Liu (2012)** : "REBOUND: an open-source multi-purpose N-body code for collisional dynamics", *Astronomy & Astrophysics*, 537, A128
  - Article principal sur la bibliothèque rebound
  - DOI: https://doi.org/10.1051/0004-6361/201118085

### Échantillonnage bayésien
- **Buchner (2021)** : "UltraNest - a robust, general purpose Bayesian inference engine"
  - Article principal sur la bibliothèque UltraNest
  - URL: https://johannesbuchner.github.io/UltraNest/

- **Buchner (2014, 2019)** : MLFriends - Nested sampling Monte Carlo algorithm
  - Algorithme utilisé par UltraNest pour le nested sampling
  - Permet de dériver les distributions de probabilité a posteriori et l'évidence bayésienne

- **dynesty** : Bibliothèque Python pour le nested sampling dynamique (https://dynesty.readthedocs.io/)
  - Utilisée par Prospector pour l'inférence bayésienne (alternative à emcee)
  - Permet l'estimation de l'évidence bayésienne et des distributions a posteriori
  - Référence : Speagle (2020), "Dynesty: a dynamic nested sampling package for estimating Bayesian posteriors and evidences", *MNRAS*, 493, 3132

**Notes** : 
- Pour obtenir la citation complète et à jour de rebound avec toutes les références, exécutez dans Python :
  ```python
  import rebound
  sim = rebound.Simulation()
  sim.cite()
  ```
- UltraNest utilise l'algorithme MLFriends (Buchner, 2014; 2019) pour le nested sampling Monte Carlo, permettant de dériver les distributions de probabilité a posteriori et l'évidence bayésienne.

## Licences

NPOAP est développé en utilisant des bibliothèques sous diverses licences open-source (principalement BSD, MIT, et Apache 2.0). Nous respectons les licences de toutes les bibliothèques utilisées.

Pour plus d'informations sur les licences spécifiques, consultez :
- Les fichiers LICENSE ou COPYRIGHT de chaque bibliothèque
- Les sites web officiels des projets mentionnés ci-dessus

## Remerciements spéciaux

Nous remercions particulièrement :
- **Ben Sharkey** pour le projet Zero-Aperture-Astrometry et la méthodologie d’extrapolation zero-aperture
- Les équipes du **NASA Exoplanet Archive** / **IPAC** et du **Minor Planet Center** pour l’accès aux catalogues et services d’éphémérides
- Le projet **ExoClock** (groupe de travail éphémérides ARIEL) et les contributeurs du module **exoclock**
- La **Section étoiles variables et exoplanètes** de la **Société astronomique tchèque** pour **VarAstro** et la base **ETD**
- La communauté **Astropy** pour leurs outils fondamentaux
- L'équipe **STDPipe** (Sergey Karpov) pour leur pipeline de détection de transitoires
- L'équipe **PHOEBE** pour leur bibliothèque de modélisation d'étoiles binaires
- L'équipe **BD-J** (Benjamin Johnson et collaborateurs) pour Prospector et sedpy, outils essentiels pour l'analyse spectroscopique avancée
- L'équipe **FSPS** (Charlie Conroy et collaborateurs) pour leur code de synthèse de populations stellaires
- L'équipe **UCL Exoplanets** pour pylightcurve
- **Michael Hippke** et collaborateurs pour **Wōtan** (wotan), outil de référence pour le détrendage de courbes de lumière de transits
- La communauté astronomique open-source pour ses contributions continues et le partage de connaissances qui rendent des outils comme NPOAP possibles

---

**NPOAP - Acknowledgments**

*Merci à tous les développeurs et contributeurs de ces projets open-source qui rendent NPOAP possible.*

