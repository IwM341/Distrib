conversions.exe -config config_rd.txt -grid capt/l_grid.txt -values evolute/cl_o.bmat -LE capt/LE.txt -o graphics/CL_FINAL.txt
conversions.exe -config config_rd.txt -grid capt/h_grid.txt -values evolute/ch_o.bmat  -LE capt/LE.txt -o graphics/CH_FINAL.txt

conversions.exe -config config_rd.txt -grid capt/l_grid.txt -values evolute/dl_o.bmat -LE capt/LE.txt -o graphics/NL_FINAL.txt
conversions.exe -config config_rd.txt -grid capt/h_grid.txt -values evolute/dh_o.bmat -LE capt/LE.txt -o graphics/NH_FINAL.txt


r_distrib.exe  -config config_rd.txt -distrib.grid capt/l_grid.txt -distrib.values evolute/cl_o.bmat -out graphics/CL_R_FINAL.txt
r_distrib.exe  -config config_rd.txt -distrib.grid capt/h_grid.txt -distrib.values evolute/ch_o.bmat -out graphics/CH_R_FINAL.txt
r_distrib.exe  -config config_rd.txt -distrib.grid capt/l_grid.txt -distrib.values evolute/dl_o.bmat -out graphics/NL_R_FINAL.txt
r_distrib.exe  -config config_rd.txt -distrib.grid capt/h_grid.txt -distrib.values evolute/dh_o.bmat -out graphics/NH_R_FINAL.txt

annihilation_rd.exe -config config_rd.txt -l_d.grid capt/l_grid.txt -l_d.values evolute/dl_o.bmat -h_d.grid capt/h_grid.txt -h_d.values evolute/dh_o.bmat -out graphics/ANN_HL_R.txt
annihilation_rd.exe -config config_rd.txt -l_d.grid capt/l_grid.txt -l_d.values evolute/dl_o.bmat -h_d.grid capt/l_grid.txt -h_d.values evolute/dl_o.bmat -out graphics/ANN_LL_R.txt
annihilation_rd_vector -config config_rd.txt distrib.grid capt/l_grid.txt distrib.values evolute/dl_o.bmat -out graphics/ANN_L.txt