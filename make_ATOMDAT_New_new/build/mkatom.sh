#!/bin/bash

#   define catalog path-variables ---------------------------
  export ATOM_SRC=../Atomic-Models-Source
  export EKUR_DIR=$ABDIR/EmKurFlux-Bin-Source
  export SNHD_DIR=$ABDIR/SNIa-Hydro
  export ATOM_DIR=$ABDIR/Atomic-Models-ASCII
  export LINE_DIR=$ABDIR/Line-List-Bin-Source
  export HTML_DIR=html
  export LOGO=wdst.ps


#   print variables and warning -----------------------------
  echo "
            **************************************
            *       CREATE ATOMIC DATABASE       *
            **************************************

This batch script creates atomic data files for use by WM-basic in
\"${ATOM_DIR}\"

from sources in
\"${ATOM_SRC}\"

The corresponding grotrian diagrams and the HTML catalogue is written to
\"${ATOM_SRC}/${HTML_DIR}\"

Remember to set WMD_DIR correctly!
WMBDIR: \"${WMBDIR}\"

            **************************************
            *  EXISTING DATA WILL BE DELETED !!  *
            *                                    *
            *     Press ENTER to continue,       *
            *           CTRL-C to exit.          *
            **************************************
"
  read dummy

#   make directories --------------------------------------
  echo "
           Deleting old/preparing new data...
"
# rm -r ${ATOM_SRC}/${HTML_DIR}
  mkdir -p ${ATOM_SRC}/${HTML_DIR}
  mkdir -p ${ATOM_SRC}/${HTML_DIR}/cat
  mkdir -p ${ATOM_SRC}/${HTML_DIR}/plot
  mkdir -p ${ATOM_SRC}/${HTML_DIR}/data

  mkdir -p ${ABDIR}
  mkdir -p ${ATOM_DIR}
  mkdir -p ${LINE_DIR}

  cp -pr ${EKUR_DIR} ${ABDIR}
  cp -pr ${SNHD_DIR} ${ABDIR}

  cp ${ATOM_SRC}/_atdiestwnel ${ATOM_DIR}/atdiestwnel
  cp ${ATOM_SRC}/_ATOMIC      ${ATOM_DIR}/ATOMIC


#   running mkatom.exe ------------------------------------
  echo "
           Running mkatom...
"
  ${WMBDIR}/build/mkatom


#   running idl -------------------------------------------
  echo "
           Creating PS-files...
"
  FILE=IDL_grot.pro
  echo "!quiet=1
    device,retain=2,TRANSLATION=TRANSARR,decompose=0
    .r grot
    exit" >${FILE}
  IDL_STARTUP=${FILE} export IDL_STARTUP

  if test -n "$COMSPEC"
  then
    cygwrap "$IDLDE"
  else
    idl
  fi

  rm ${FILE}


#   running ghostscript -----------------------------------
  echo "
           Creating PDF- and PNG-files...
"
  cp ${LOGO} ${ATOM_SRC}/${HTML_DIR}/plot/logo.ps
  cd ${ATOM_SRC}/${HTML_DIR}/plot
  sed -e '
    /@echo off/d
    s/@echo/echo/
    s/gswin32c/gs/
  ' mkplot.bat >mkplot.sh
  sh mkplot.sh
  cd ..


#   deleting files ----------------------------------------
  rm data/ener.*
  rm data/tran.*
  rm data/_list.txt
  rm plot/mkplot.*


  echo "
            **************************************
            *                                    *
            *                DONE!               *
            *                                    *
            **************************************
"
