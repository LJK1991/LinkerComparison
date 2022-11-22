import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from upsetplot import UpSet
from upsetplot import from_memberships

#load files with all the barcodes
def file_to_list(file):
    temp = []
    with open(file) as f:
        for line in f:
            temp.append(line.strip("\n"))
    return temp

#numerize everything
#read the SPLiT-seqbarcode set
SPLIT = file_to_list("/home/lucas/Documents/Work/2022/10-2022/SPLTsq_Round1BC.txt")
SPLTBC = pd.DataFrame({'BC':SPLIT,'num':list(range(1,97))})

#seperat strings make loop and assingn them numbers

LRsplitpipe = file_to_list("/home/lucas/Documents/Work/2022/10-2022/LRsplitpipe_cell_barcodes.txt")
paulranum = file_to_list("/home/lucas/Documents/Work/2022/10-2022/paulranum_cell_barcodes.txt")
SCSit = file_to_list("/home/lucas/Documents/Work/2022/10-2022/SCSit_cell_barcodes.txt")
#STARsolo = file_to_list("/home/lucas/Documents/Work/2022/10-2022/STARsolo_cell_barcodes.txt")
STARsolo = file_to_list("/home/lucas/Documents/Work/2022/10-2022/STARsolo_new_cell_barcodes.txt")
yzhang = file_to_list("/home/lucas/Documents/Work/2022/10-2022/yzhang_cell_barcodes.txt")
zUMI = file_to_list("/home/lucas/Documents/Work/2022/10-2022/zUMI_cell_barcodes.txt")

#LRsplitpipe
alist = []
blist = []
clist = []
for cell in LRsplitpipe:
    a = int(SPLTBC['num'].loc[SPLTBC['BC'] == cell[:8]].values)
    b = int(SPLTBC['num'].loc[SPLTBC['BC'] == cell[8:16]].values)
    c = int(SPLTBC['num'].loc[SPLTBC['BC'] == cell[16:]].values)
    alist.append(a)
    blist.append(b)
    clist.append(c)
LRdf = pd.DataFrame({'1':alist,'2':blist,'3':clist})
LRdf['BC'] = LRdf['1'].astype(str) + "_" + LRdf['2'].astype(str) + "_" + LRdf['3'].astype(str)
LRdf.to_csv("/home/lucas/Documents/Work/2022/10-2022/LRbcdist.tsv", sep="\t")

#paulranum
alist = []
blist = []
clist = []
for cell in paulranum:
    a = int(SPLTBC['num'].loc[SPLTBC['BC'] == cell[:16].replace('AGCATTCG','')].values)
    b = int(SPLTBC['num'].loc[SPLTBC['BC'] == cell[16:29].replace('ATCCA','')].values)
    c = int(SPLTBC['num'].loc[SPLTBC['BC'] == cell[29:].replace('GTGGCC','')].values)
    alist.append(a)
    blist.append(b)
    clist.append(c)
pauldf = pd.DataFrame({'1':alist,'2':blist,'3':clist})
pauldf['BC'] = pauldf['1'].astype(str) + "_" + pauldf['2'].astype(str) + "_" + pauldf['3'].astype(str)
pauldf.to_csv("/home/lucas/Documents/Work/2022/10-2022/paulbcdist.tsv", sep="\t")

#SCSit
#did not do barcode conversion as they said they would
alist = []
blist = []
clist = []
for cell in SCSit:
    cell = cell.strip("X")
    a = int(cell.split("_")[0])
    b = int(cell.split("_")[1])
    c = int(cell.split("_")[2])
    if a > 48:
        a = a - 48
    alist.append(a)
    blist.append(b)
    clist.append(c)
SCSdf = pd.DataFrame({'1':alist,'2':blist,'3':clist})
SCSdf['BC'] = SCSdf['1'].astype(str) + "_" + SCSdf['2'].astype(str) + "_" + SCSdf['3'].astype(str)
SCSdf.to_csv("/home/lucas/Documents/Work/2022/10-2022/SCSbcdist.tsv", sep="\t")

#STARsolo
alist = []
blist = []
clist = []
for cell in STARsolo:
    a = int(SPLTBC['num'].loc[SPLTBC['BC'] == cell[:8]].values)
    b = int(SPLTBC['num'].loc[SPLTBC['BC'] == cell[8:16]].values)
    c = int(SPLTBC['num'].loc[SPLTBC['BC'] == cell[16:]].values)
    alist.append(a)
    blist.append(b)
    clist.append(c)
STARdf = pd.DataFrame({'1':alist,'2':blist,'3':clist})
STARdf['BC'] = STARdf['1'].astype(str) + "_" + STARdf['2'].astype(str) + "_" + STARdf['3'].astype(str)
STARdf.to_csv("/home/lucas/Documents/Work/2022/10-2022/STARbcdist.tsv", sep="\t")

#yzhang
alist = []
blist = []
clist = []
for cell in yzhang:
    a = int(cell.split('_')[1])
    cell = cell.split('_')[0]
    b = int(SPLTBC['num'].loc[SPLTBC['BC'] == cell[8:16]].values)
    c = int(SPLTBC['num'].loc[SPLTBC['BC'] == cell[:8]].values)
    alist.append(a)
    blist.append(b)
    clist.append(c)
yzhangdf = pd.DataFrame({'1':alist,'2':blist,'3':clist})
yzhangdf['BC'] = yzhangdf['1'].astype(str) + "_" + yzhangdf['2'].astype(str) + "_" + yzhangdf['3'].astype(str)
yzhangdf.to_csv("/home/lucas/Documents/Work/2022/10-2022/yzhangbcdist.tsv", sep="\t")


#zUMI
alist = []
blist = []
clist = []
for cell in zUMI:
    a = int(SPLTBC['num'].loc[SPLTBC['BC'] == cell[:8]].values)
    b = int(SPLTBC['num'].loc[SPLTBC['BC'] == cell[8:16]].values)
    c = int(SPLTBC['num'].loc[SPLTBC['BC'] == cell[16:]].values)
    if c > 48:
        c = c - 48
    alist.append(c)
    blist.append(b)
    clist.append(a)
zUMIdf = pd.DataFrame({'1':alist,'2':blist,'3':clist})
zUMIdf['BC'] = zUMIdf['1'].astype(str) + "_" + zUMIdf['2'].astype(str) + "_" + zUMIdf['3'].astype(str)
zUMIdf.to_csv("/home/lucas/Documents/Work/2022/10-2022/zUMIbcdist.tsv", sep="\t")

'''
total_BC = list(LRdf['BC'])+list(pauldf['BC'])+list(SCSdf['BC'])+list(STARdf['BC'])+list(yzhangdf['BC'])+list(zUMIdf['BC'])
total_BC = list(set(total_BC))
print(total_BC)
indicator_df = pd.DataFrame({'BC':total_BC})
indicator_df['LR-splitpipe'] = False
indicator_df['SPLiT-seq demultiplex'] = False
indicator_df['SCSit'] = False
indicator_df['STARsolo'] = False
indicator_df['splitseq pipeline'] = False
indicator_df['zUMI'] = False
print(indicator_df.head())
print(indicator_df.columns)
for i in indicator_df.index:
    if indicator_df.BC[i] in LRdf['BC']:
        indicator_df.loc[i, 'LR-splitpipe'] = True
    if indicator_df.BC[i] in pauldf['BC']:
        indicator_df.loc[i, 'SPLiT-seq demultiplex'] = True
    if indicator_df.BC[i] in SCSdf['BC']:
        indicator_df.loc[i, 'SCSit'] = True
    if indicator_df.BC[i] in STARdf['BC']:
        indicator_df.loc[i, 'STARsolo'] = True
    if indicator_df.BC[i] in yzhangdf['BC']:
        indicator_df.loc[i, 'splitseq pipeline'] = True
    if indicator_df.BC[i] in zUMIdf['BC']:
        indicator_df.loc[i, 'zUMI'] = True

indicator_df.set_index('BC', inplace=True)

#indicator_df = indicator_df.T
#print(indicator_df)
#indicator_df.pop('BC')

UpSet(from_memberships(indicator_df),indicator_df)

plt.show()
'''
from upsetplot import from_contents

setlist = {'LR-splitpipe':set(LRdf['BC']), 'SPLiT-seq demultiplex':set(pauldf['BC']), 'SCSit':set(SCSdf['BC']), 'STARsolo':set(STARdf['BC']), 'split-seq pipeline':set(yzhangdf['BC']),'zUMI':set(zUMIdf['BC'])}
x = from_contents(setlist)
UpSet(x).plot()
plt.show()
#labels = ['LR-splitpipe', 'paulranum SPLiTseq demultiplex', 'SCSit', 'STARsolo', 'yjzhang', 'zUMI']
#supervenn(setlist, labels)
#plt.savefig('/home/lucas/Documents/Work/2022/10-2022/supervenn.png')