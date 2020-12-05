#!/usr/bin/sh
#mysql -h focal.usc.edu -u cybershk_ro -pCyberShake2007 -DCyberShake <query_rotd50.sql > Northridge_sites.txt

mysql -h focal.usc.edu -u cybershk_ro -pCyberShake2007 -DCyberShake <query_events.sql > Northridge_events.txt
