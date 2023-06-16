..\..\release\conversions.exe -config config_rd.txt -grid capt/l_grid -values evolute/c_o.bmat -o graphics/C_FINAL.txt

..\..\release\conversions.exe -config config_rd.txt -grid capt/l_grid -values evolute/d_o.bmat -o graphics/N_FINAL.txt


..\..\release\r_distrib.exe  -config config_rd.txt -distrib.grid capt/l_grid -distrib.values evolute/c_o.bmat -o graphics/C_R_FINAL.txt

..\..\release\r_distrib.exe  -config config_rd.txt -distrib.grid capt/l_grid -distrib.values evolute/d_o.bmat -o graphics/N_R_FINAL.txt

..\..\release\annihilation_rd.exe -config config_rd.txt -l_d.grid capt/l_grid.txt -l_d.values evolute/d_o.bmat -h_d.grid capt/l_grid.txt -h_d.values evolute/d_o.bmat -out graphics/ANN_R.txt
