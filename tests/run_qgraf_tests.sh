#!/bin/bash

QGRAF_VERSION="3.6.7"
QGRAF_DIR="qgraf-${QGRAF_VERSION}"
SITE="http://qgraf.tecnico.ulisboa.pt"
FORT=gfortran

if [[ ! -e "${QGRAF_DIR}" ]]; then
    echo "Missing qgraf dependency, downloading."
    wget --quiet --user=anonymous --password=anonymous -O ./qgraf.tgz \
         "${SITE}/links/qgraf-3.6.7.tgz" || exit 1
    mkdir "${QGRAF_DIR}" || exit 1
    tar -C "${QGRAF_DIR}" -xvf qgraf.tgz || exit 1
    (
        cd "${QGRAF_DIR}" &&
            ${FORT} "qgraf-${QGRAF_VERSION}.f08" -o qgraf
    )
fi

