import numpy as np
import pandas
import sys

psa = pandas.read_csv('/gpfs/alpine/geo112/scratch/hzfmer/LA/gp_rupture_test/LA/data_for_Zhifeng/study_15_4_psa.csv')

site_names6 = [line.rstrip('\n') for line in open('/gpfs/alpine/geo112/scratch/hzfmer/LA/gp_rupture_test/LA/data_for_Zhifeng/site_names6_35.txt')]
site_names7 = [line.rstrip('\n') for line in open('/gpfs/alpine/geo112/scratch/hzfmer/LA/gp_rupture_test/LA/data_for_Zhifeng/site_names7_35.txt')]
site_names8 = [line.rstrip('\n') for line in open('/gpfs/alpine/geo112/scratch/hzfmer/LA/gp_rupture_test/LA/data_for_Zhifeng/site_names8_45.txt')]

nsite = {6.35:len(site_names6), 7.35:len(site_names7), 8.45:len(site_names8)}
site_names = {6.35:site_names6, 7.35:site_names7, 8.45:site_names8}
src_id = {6.35:76, 7.35:128, 8.45:68}
for M in [6.35, 7.35, 8.45]:
    for T in [1, 2, 5]:
        pd = psa[(psa['Source_ID'] == src_id[M]) & (psa['IM_Type_Component'] == 'RotD50') & (psa['IM_Type_Value'] == T)]
        pd.reset_index(inplace=True, drop=True)
        nrows, _ = pd.shape
        print(f'M={M}, T={T}, nrows={nrows}')
        sa = np.zeros((nsite[M], 1), dtype='f')
        for i in range(nrows):
            try:
                idx = site_names[M].index(pd.loc[i, 'CS_Short_Name'])
                sa[idx] = pd.loc[i, 'IM_Value'] / 100 / 9.8 # convert from cm/s^2 to g
            except IndexError:
                print(f'i={i}, nrows={nrows}')
                print(pd.loc[i, :])
                sys.exit(-1)
            except KeyError:
                print(f"KeyError: pd.loc[i, 'CS_Short_Name']")
                continue
        print(f'Max PSA = {sa.max()}')
        sa.tofile(f"../sa_cb_site_M{str(M).replace('.', '_')}_T{T}.bin")

