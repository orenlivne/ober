@echo off
set PYTHONPATH=C:/ober/impute;C:/ober/hera;C:/ober/network;C:/ober/util;C:/ober/famplot
set TEST_DATA_DIR=C:/ober/testdata
set OBER=C:/ober
set PATH=C:\cygwin\home\oren\plink-1.07-dos;%PATH%

rem python -c "import sys; print sys.path"
python %*
