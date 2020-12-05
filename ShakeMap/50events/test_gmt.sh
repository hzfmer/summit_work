gmt begin taiwan_relief pdf
gmt grdimage @earth_relief_01s -JM15c -R121.9/122/25.9/26 -Baf -BWSen -I+d
# gmt colorbar -DJMR+w10c+o1.5c/0c+ml -Bxa1000f -By+l"m"
#gmt grdimage @earth_relief_30s -JM15c -R118/125/20/26 -Baf -BWSen -I+d
#gmt colorbar -DJMR+w10c+o1.5c/0c+ml -Bxa10f -By+l"m"
gmt end
