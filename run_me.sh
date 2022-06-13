#!/bin/bash

starver SL21c
cons
#root4star -q -l -b analyzePico.C\(\"./st_physics_20057003_raw_1500005.picoDst.root\",\"test\",\"config_19p6GeV.txt\"\)
root4star -b -q -l analyzeTree.C\(3,1,\"0E801ACD6F61FBAC58201E930AD9B7F7_1421_ME.root\",\"test\",0,1000000000,0\)
echo 'Good morning/afternoon/evening/night! I love you, Ding Chen!'
