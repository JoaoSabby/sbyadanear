#!/bin/bash
set -e
source /opt/intel/oneapi/setvars.sh >/dev/null
exec "$@"
