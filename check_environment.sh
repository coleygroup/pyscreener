#!/bin/bash

software=$1

if [ "$software" = "vina" ] \
   || [ "$software" == "smina"  ] \
   || [ "$software" == "psovina" ] \
   || [ "$software" == "qvina" ]; then
    if ! which prepare_receptor >/dev/null ; then
        echo prepare_receptor is not on PATH!
        false
    fi
    if ! which $software >/dev/null ; then
        echo $software is not on PATH!
        false
    fi
elif [ "$software" == "dock6" ]; then
    VDW_DEFN_FILE=$DOCK6/parameters/vdw_AMBER_parm99.defn
    FLEX_DEFN_FILE=$DOCK6/parameters/flex.defn
    FLEX_DRIVE_FILE=$DOCK6/parameters/flex_drive.tbl
    DOCK=$DOCK6/bin/dock6

    if [ ! -d $DOCK6 ]; then
        echo DOCK6 environment variable is not set or does not exist!
        false
    fi
    if [ ! -f $VDW_DEFN_FILE ]; then
        echo DOCK6 directory is not configured properly! no vdw definition file at $VDW_DEFN_FILE!
        false
    fi
    if [ ! -f $FLEX_DEFN_FILE ]; then
        echo DOCK6 directory is not configured properly! no flex definition file at $FLEX_DEFN_FILE!
        false
    fi
    if [ ! -f $FLEX_DRIVE_FILE ]; then
        echo DOCK6 directory is not configured properly! no flex drive file at $FLEX_DRIVE_FILE!
        false
    fi
    if [ ! -f $DOCK ]; then
        echo DOCK6 directory is not configured properly! no dock6 executable at $DOCK!
        false
    fi
elif [ -z "$software" ]; then
    echo "no software provided!"
    exit 1
else
    echo $software is not a supported psycreener software!
    exit 1
fi

echo environment is set up correctly to use $software!