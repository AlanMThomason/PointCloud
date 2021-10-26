python -m numpy.f2py -c --compiler=msvc FT_20170106a.f90 -m FT_20170106a
copy .\FT_20170106a\.libs\*.dll .\