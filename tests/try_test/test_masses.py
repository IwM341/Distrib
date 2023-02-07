import test_functions as tf
import pandas as pd

'''python script for testing programm'''
'''masses - array of masses being tested'''
mk = [10,30,50,70,100]
'''dmk - delta mass'''
dmk = 0.0

'''remark: to test dmk with fixed mass, 
        just swap types (i.e. write mk =  mass, dmk = [list of dmk])'''


'''tf.test_parametr_array_mks runs programm, 
to specify path to programm, write in args prog_name = path/to/prog.
to test dmk use function test_parametr_array_dmks'''
#tr = tf.test_parametr_array_mks(mk,dmk,parallel = False)

'''tf.load_list_mk loads all fraction files. Use it if programm has worked and previous line is commented.
if files are not ready, uncomment previous line and comment this line'''
tr = tf.load_list_mk(mk,dmk)

'''save results to datatable'''
df = pd.DataFrame([ tf.dict_concat(el['H'],{'mk':el['mk'],'dmk':el['dmk']}) for el in tr])
df.to_csv('df.txt', header=None, index=None, sep = '\t')

'''this line prints gnuplot command to plot distributions'''
print(tf.plotting_command(df.columns,'df.txt','mk'))
