..\..\release\annihilation.exe -config config_rd.txt -l_grid capt/l_grid.txt -h_grid capt/h_grid.txt -ann ann/ann_hl.bmat
..\..\release\annihilation.exe -config config_rd.txt -l_grid capt/l_grid.txt -h_grid capt/l_grid.txt -ann ann/ann_ll.bmat
..\..\release\annihilation.exe -config config_rd.txt -l_grid capt/h_grid.txt -h_grid capt/h_grid.txt -ann ann/ann_hh.bmat

..\..\release\annihilation_vector.exe -config config_rd.txt -grid capt/l_grid.txt -ann ann/ann_l.bmat
..\..\release\annihilation_vector.exe -config config_rd.txt -grid capt/h_grid.txt -ann ann/ann_h.bmat