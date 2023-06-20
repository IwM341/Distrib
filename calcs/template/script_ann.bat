annihilation -config config_rd.txt -l_grid capt/l_grid.txt -h_grid capt/h_grid.txt -o ann/ann_hl.bmat
annihilation -config config_rd.txt -l_grid capt/l_grid.txt -h_grid capt/l_grid.txt -o ann/ann_ll.bmat
annihilation -config config_rd.txt -l_grid capt/h_grid.txt -h_grid capt/h_grid.txt -o ann/ann_hh.bmat

annihilation_vector -config config_rd.txt -grid capt/l_grid.txt -o ann/ann_el.bmat
annihilation_vector -config config_rd.txt -grid capt/h_grid.txt -o ann/ann_eh.bmat