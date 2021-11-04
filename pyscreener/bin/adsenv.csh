#!/bin/csh
######
## Set some environment variables.
setenv ADS_ROOT /home/mburlage/pyscreener/ADFRsuite-1.0/ADFRsuite_x86_64Linux_1.0/pyscreener

########
## plaform we run on
##
setenv ADS_ARCHOSV `$ADS_ROOT/bin/archosv`

#######
## path to the extralibs directory.
##
setenv ADS_EXTRALIBS $ADS_ROOT/lib

#######
## path to the extrainclude directory
setenv ADS_EXTRAINCLUDE $ADS_ROOT/include

########
## add the path to the directory holding the python interpreter to your path
##
set path=($ADS_ROOT/bin:$path)
# Open Babel formats, plugins directory:
setenv BABEL_LIBDIR $ADS_ROOT/lib/openbabel/2.4.1
setenv BABEL_DATADIR $ADS_ROOT/share/openbabel/2.4.1

#REDUCE
setenv REDUCE_HET_DICT $ADS_ROOT/bin


# set the LD_LIBRARY PATH for each platform

if (`uname -s` == Darwin) then
    setenv DISPLAY :0.0
    set isdefined=`printenv DYLD_LIBRARY_PATH`
    if ( $#isdefined ) then
	setenv DYLD_LIBRARY_PATH $ADS_ROOT/lib:$DYLD_LIBRARY_PATH
    else
	setenv DYLD_LIBRARY_PATH $ADS_ROOT/lib
    endif

else
    set isdefined=`printenv LD_LIBRARY_PATH`
    if ( $#isdefined ) then
	setenv LD_LIBRARY_PATH $ADS_ROOT/lib:$LD_LIBRARY_PATH
    else
	setenv LD_LIBRARY_PATH $ADS_ROOT/lib
    endif

endif

# use python interpreter that comes with the tools

unset PYTHONHOME
setenv PYTHONHOME $ADS_ROOT
setenv PYTHONPATH $ADS_ROOT/CCSBpckgs
setenv python $ADS_ROOT/bin/python

