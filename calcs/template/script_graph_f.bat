conversions -config config_rd.txt -grid capt/l_grid.txt -values evolute/cl_o.bmat -LE capt/LE.txt -o graphics/CL_FINAL.txt
conversions -config config_rd.txt -grid capt/h_grid.txt -values evolute/ch_o.bmat  -LE capt/LE.txt -o graphics/CH_FINAL.txt

conversions -config config_rd.txt -grid capt/l_grid.txt -values evolute/dl_o.bmat -LE capt/LE.txt -o graphics/NL_FINAL.txt
conversions -config config_rd.txt -grid capt/h_grid.txt -values evolute/dh_o.bmat -LE capt/LE.txt -o graphics/NH_FINAL.txt


r_distrib  -config config_rd.txt -distrib.grid capt/l_grid.txt -distrib.values evolute/cl_o.bmat -o graphics/CL_R_FINAL.txt
r_distrib  -config config_rd.txt -distrib.grid capt/h_grid.txt -distrib.values evolute/ch_o.bmat -o graphics/CH_R_FINAL.txt
r_distrib  -config config_rd.txt -distrib.grid capt/l_grid.txt -distrib.values evolute/dl_o.bmat -o graphics/NL_R_FINAL.txt
r_distrib  -config config_rd.txt -distrib.grid capt/h_grid.txt -distrib.values evolute/dh_o.bmat -o graphics/NH_R_FINAL.txt

annihilation_rd -config config_rd.txt -l_d.grid capt/l_grid.txt -l_d.values evolute/dl_o.bmat -h_d.grid capt/h_grid.txt -h_d.values evolute/dh_o.bmat -o graphics/ANN_HL_R.txt
annihilation_rd -config config_rd.txt -l_d.grid capt/l_grid.txt -l_d.values evolute/dl_o.bmat -h_d.grid capt/l_grid.txt -h_d.values evolute/dl_o.bmat -o graphics/ANN_LL_R.txt
annihilation_rd_vector -config config_rd.txt distrib.grid capt/l_grid.txt distrib.values evolute/dl_o.bmat -o graphics/ANN_L.txt