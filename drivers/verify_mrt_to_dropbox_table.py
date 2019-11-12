import pandas as pd, numpy as np
from astropy.table import Table

df1 = pd.read_csv('../data/string_table1.csv')
df2 = pd.read_csv('../data/string_table2.csv')
df3 = pd.read_csv('../data/string_table3.csv')

t1 = Table.read('../data/ajab339at1_mrt.txt', format='ascii.cds')
t2 = Table.read('../data/ajab339at2_mrt.txt', format='ascii.cds')
t3 = Table.read('../data/ajab339at3_mrt.txt', format='ascii.cds')

np.testing.assert_array_equal(
    list(map(str,df1.source_id)),
    list(map(str,t1['Gaia'])),
)

print('verified csv and MRT tables have equal source_id')

np.testing.assert_allclose(
    list(map(float,df2.age)),
    list(map(float,t2['logage'])),
)

np.testing.assert_allclose(
    list(map(float,df3.l)),
    list(map(float,t3['GLON'])),
)

print('verified table2 and table3 params close to rtol=1e-7')
