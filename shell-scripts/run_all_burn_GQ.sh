#!/bin/bash
#
# execute all burn model generated-quantities shell scripts
#

shell-scripts/run_g1_GQ.sh
sleep 1
shell-scripts/run_g2_GQ.sh
sleep 1
shell-scripts/run_lognorm_GQ.sh
