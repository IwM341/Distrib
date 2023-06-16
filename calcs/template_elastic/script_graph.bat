..\..\release\conversions.exe -config config_rd.txt -grid capt/l_grid -values capt/o_hl.bmat -func -o graphics/N-.txt

..\..\release\conversions.exe -config config_rd.txt -grid capt/l_grid -values capt/cl.bmat -o graphics/CL_INIT.txt

..\..\release\r_distrib.exe  -config config_rd.txt -distrib.grid capt/l_grid -distrib.values capt/cl.bmat -o graphics/CL_R_INIT.txt
