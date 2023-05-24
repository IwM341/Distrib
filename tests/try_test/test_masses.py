import test_functions as tf
import pandas as pd

'''python script for testing programm'''
'''masses - array of masses being tested'''
#mk = [10,15,20,30,40,50,60,70,85,100]
mk = 10
'''dmk - delta mass'''
dmk = [(i-4)*0.5e-4 for i in range(9)]
'''remark: to test dmk with fixed mass, 
        just swap types (i.e. write mk =  mass, dmk = [list of dmk])'''


'''tf.test_parametr_array_mks runs programm, 
to specify path to programm, write in args prog_name = path/to/prog.
to test dmk use function test_parametr_array_dmks'''


'''tf.load_list_mk loads all fraction files. Use it if programm has worked and previous line is commented.
if files are not ready, uncomment previous line and comment this line'''
#tr = tf.load_list_mk(mk,dmk)


df_fname = 'df_dmk10'

calculate = False

'''save results to datatable'''
if(calculate):
        tr = tf.test_parametr_array_dmks(mk,dmk,parallel = True)
        df = pd.DataFrame([ tf.dict_concat(el['H'],{'mk':el['mk'],'dmk':el['dmk']}) for el in tr])
        df.to_csv(df_fname + "_pure.txt", header=None, index=None, sep = '\t')
        df.to_csv(df_fname + ".txt", index=None, sep = '\t')
else:
        df = pd.read_csv(df_fname + ".txt", sep = '\t')


'''this line prints gnuplot command to plot distributions'''
print(tf.plotting_command(df.columns,df_fname + "_pure.txt",'dmk'))
