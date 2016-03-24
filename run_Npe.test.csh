#!/bin/bash

starver SL12d

source cons.sh
rm -rf Npe_test.root

root4star -l -b StRoot/macros/run_StNpeMaker.C\(\"Npefile.list\"\,\"Npe_test.root\"\) 
