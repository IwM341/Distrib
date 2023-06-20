..\..\release\conversions.exe -config config_rd.txt -grid capt/l_grid.txt -values capt/o_hl.bmat -func -o graphics/N-.txt

..\..\release\conversions.exe -config config_rd.txt -grid capt/l_grid.txt -values capt/cl.bmat -o graphics/CL_INIT.txt

..\..\release\r_distrib.exe  -config config_rd.txt -distrib.grid capt/l_grid.txt -distrib.values capt/cl.bmat -out graphics/CL_R_INIT.txt
