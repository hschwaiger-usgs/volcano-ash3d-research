echo -n "      character(len=40),parameter,public :: Ash3d_GitComID ='" > version.h
git log -n 1 | grep commit | cut -f 2 -d' ' | tr -d $'\n' >> version.h
echo -n "'" >> version.h
