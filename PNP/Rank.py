import scipy.io as sio
import pandas as pd 
from tqdm import tqdm
EnergyMat1=sio.loadmat("BB1.mat")["BB1"]
EnergyMat2=sio.loadmat("BB2.mat")["BB2"]
GUP=sio.loadmat("Gup.mat")["Gup"]
GDN=sio.loadmat("Gdn.mat")["Gdn"]
UP=pd.DataFrame({'Gene':[str(i[0][0]) for i in GUP],'Energy':[i[0] for i in EnergyMat2]})
DOWN=pd.DataFrame({'Gene':[str(i[0][0]) for i in GDN],'Energy':[i[0] for i in EnergyMat1]})
ALL_RANK=UP.append(DOWN)
Final_RAnk=pd.DataFrame({'Gene_Name':list(set(ALL_RANK.Gene))})
Final_RAnk['Energy Mean']=[ALL_RANK[ALL_RANK.Gene==i]['Energy'].mean() for i in tqdm(Final_RAnk.Gene_Name)]
Final_RAnk=Final_RAnk.sort_values(by='Energy Mean')
#%%
Final_RAnk.to_csv('Rank.csv',index=False)