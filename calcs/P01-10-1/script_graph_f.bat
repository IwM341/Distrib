..\..\release\conversions.exe -config config_rd.txt -grid capt/l_grid.txt -values evolute/cl_o.bmat -LE capt/LE.txt -o graphics/CL_FINAL.txt
..\..\release\conversions.exe -config config_rd.txt -grid capt/h_grid.txt -values evolute/ch_o.bmat  -LE capt/LE.txt -o graphics/CH_FINAL.txt

..\..\release\conversions.exe -config config_rd.txt -grid capt/l_grid.txt -values evolute/dl_o.bmat -LE capt/LE.txt -o graphics/NL_FINAL.txt
..\..\release\conversions.exe -config config_rd.txt -grid capt/h_grid.txt -values evolute/dh_o.bmat -LE capt/LE.txt -o graphics/NH_FINAL.txt


..\..\release\r_distrib.exe  -config config_rd.txt -distrib.grid capt/l_grid.txt -distrib.values evolute/cl_o.bmat -o graphics/CL_R_FINAL.txt
..\..\release\r_distrib.exe  -config config_rd.txt -distrib.grid capt/h_grid.txt -distrib.values evolute/ch_o.bmat -o graphics/CH_R_FINAL.txt
..\..\release\r_distrib.exe  -config config_rd.txt -distrib.grid capt/l_grid.txt -distrib.values evolute/dl_o.bmat -o graphics/NL_R_FINAL.txt
..\..\release\r_distrib.exe  -config config_rd.txt -distrib.grid capt/h_grid.txt -distrib.values evolute/dh_o.bmat -o graphics/NH_R_FINAL.txt

..\..\release\annihilation_rd.exe -config config_rd.txt -l_d.grid capt/l_grid.txt -l_d.values evolute/dl_o.bmat -h_d.grid capt/h_grid.txt -h_d.values evolute/dh_o.bmat -o graphics/ANN_HL_R.txt
..\..\release\annihilation_rd.exe -config config_rd.txt -l_d.grid capt/l_grid.txt -l_d.values evolute/dl_o.bmat -h_d.grid capt/l_grid.txt -h_d.values evolute/dl_o.bmat -o graphics/ANN_LL_R.txt
..\..\release\annihilation_rd_vector -config config_rd.txt distrib.grid capt/l_grid.txt distrib.values evolute/dl_o.bmat -o graphics/ANN_L.txt