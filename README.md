# OIM_Crystal_Python
Preparation for this code:
1. provide a file(ang) with ferrite and austenite
2. export a file(txt) containing inly austenite OR information

Crystal_V1.1 can do following tasks:
1. calculate average austenite orientation, and export the data as filename.ang
2. calculate the NW and KS OR according to the average austenite orientation, and export the data as nw_oim_or and ks_oim_or
3. calculate the OR between austenite and ferrite according to the raw data, and export data as oim_or and oim_or_centered (by using 100 010 001 austenite orientation)
4. calculate OR in three different forms of pole figures, that is 100, 110, and 111 pole figures, and export the file as origin_or_100 and so forth.
5. calculate deviation map between ferrite111 direction and austenite110 direction, and export the map as dev_110111_map
6. calculate deviation angle from average, and export the data as _deviation_map
7. calcualte the schmid factors on different slip systems (only for ausformed bainite)
8. calculate distribution of variants
9. according to the calculated schmid factors, predict the assisted bainite variants in ausformed bainite.
