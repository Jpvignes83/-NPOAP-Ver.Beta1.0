call "C:\ProgramData\Miniconda3\Scripts\activate.bat" astroenv
cd /d "C:\Users\docke\OneDrive\Documents\NPOAP"
set "PATH=C:\ProgramData\mingw64\mingw64\bin;%PATH%"
set "FC=C:\ProgramData\mingw64\mingw64\bin\gfortran.exe"
where python
python --version
gfortran --version
python -m pip install --no-cache-dir fsps