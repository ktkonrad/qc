#!/bin/sh

cp nodal_domain_driver_no_main.c nodal_domain_driver_no_main.c.backup
cp nodal_domain_driver.c nodal_domain_driver_no_main.c
sed -i -e 's/int main(/int count_main(/' nodal_domain_driver_no_main.c
sed -i -e 's/processArgs/count_processArgs/' nodal_domain_driver_no_main.c
sed -i -e 's/int verb//' nodal_domain_driver_no_main.c