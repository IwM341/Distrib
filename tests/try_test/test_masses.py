import test_functions as tf
import pandas as pd
masses = [10,30,50,70,100]
dmk = 0.0
#tr = tf.test_parametr_array_mks(masses,dmk,parallel = False)
tr = tf.load_list_mk(masses,dmk)
df = pd.DataFrame([ tf.dict_concat(el['H'],{'mk':el['mk'],'dmk':el['dmk']}) for el in tr])
df.to_csv('df.txt', header=None, index=None, sep = '\t')
print(tf.plotting_command(df.columns,'df.txt','mk'))