#!/bin/bash
date
starver SL21c

rm -rf sched*
rm *.dataset 
rm *.session.xml

rm -rf Output_SE
rm -rf Err_SE
rm -rf Log_SE
rm -rf Script_SE
rm -rf List_SE
rm -rf flow_Phi_SE
mkdir Output_SE
mkdir Err_SE
mkdir Log_SE
mkdir Script_SE
mkdir List_SE
mkdir flow_Phi_SE
rm -rf Output_ME
rm -rf Err_ME
rm -rf Log_ME
rm -rf Script_ME
rm -rf List_ME
rm -rf flow_Phi_ME
mkdir Output_ME
mkdir Err_ME
mkdir Log_ME
mkdir Script_ME
mkdir List_ME
mkdir flow_Phi_ME


chmod -R 777 Output_ME
chmod -R 777 Err_ME
chmod -R 777 Log_ME
chmod -R 777 Script_ME
chmod -R 777 List_ME
chmod -R 777 flow_Phi_ME
chmod -R 777 Output_SE
chmod -R 777 Err_SE
chmod -R 777 Log_SE
chmod -R 777 Script_SE
chmod -R 777 List_SE
chmod -R 777 flow_Phi_SE
