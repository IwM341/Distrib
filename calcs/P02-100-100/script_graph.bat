conversions.exe -config config_rd.txt -grid capt/l_grid.txt  -values capt/o_hl.bmat -LE capt/LE.txt -func -o graphics/L-.txt
conversions.exe -config config_rd.txt -grid capt/h_grid.txt -values capt/o_lh.bmat -LE capt/LE.txt -func -o graphics/H-.txt

conversions.exe -config config_rd.txt -grid capt/l_grid.txt -values capt/cl.bmat -LE capt/LE.txt -o graphics/CL_INIT.txt
conversions.exe -config config_rd.txt -grid capt/h_grid.txt -values capt/ch.bmat -LE capt/LE.txt -o graphics/CH_INIT.txt

r_distrib.exe  -config config_rd.txt -distrib.grid capt/l_grid.txt -distrib.values capt/cl.bmat -out graphics/CL_R_INIT.txt
r_distrib.exe  -config config_rd.txt -distrib.grid capt/h_grid.txt -distrib.values capt/ch.bmat -out graphics/CH_R_INIT.txt
