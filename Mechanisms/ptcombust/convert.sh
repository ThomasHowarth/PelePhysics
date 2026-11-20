#!/usr/bin/env bash

MECH_HOME="$(pwd)"
MECH_FILE="${MECH_HOME}/mechanism.yaml"
bash ../converter.sh -f "${MECH_FILE}" -- -c heterogeneous --gas_name gas --interface_name Pt_surf
