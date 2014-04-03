set HOME=C:\Users\Owner\
cd C:\eclipse\workspace\HG_AN_EXAMPLE_PYPI_PROJECT
C:\Python24\python.exe setup.py bdist_egg upload --identity="Oren Livne" --sign --quiet
C:\Python25\python.exe setup.py bdist_egg upload --identity="Oren Livne" --sign --quiet
C:\Python26\python.exe setup.py bdist_egg upload --identity="Oren Livne" --sign --quiet
C:\Python24\python.exe setup.py bdist_wininst --target-version=2.4 register upload --identity="Oren Livne" --sign --quiet
C:\Python25\python.exe setup.py bdist_wininst --target-version=2.5 register upload --identity="Oren Livne" --sign --quiet
C:\Python26\python.exe setup.py bdist_wininst --target-version=2.6 register upload --identity="Oren Livne" --sign --quiet
C:\Python26\python.exe setup.py sdist upload --identity="Oren Livne" --sign
pause