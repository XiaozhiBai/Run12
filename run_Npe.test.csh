#!/bin/bash

starver SL12d

cons
rm -rf Npe_test.root

root4star -l -b StRoot/macros/run_StNpeMaker.C\(\"Npefile.list\"\,\"Npe_test.root\"\) 