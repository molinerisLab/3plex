#!/usr/bin/env python3

# ./distMedian.py < *_lunp.lastcol > *_modif_zscore

import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# https://www.geeksforgeeks.org/plot-multiple-lines-in-matplotlib/
def moving_average(x, w):
	return np.convolve(x, np.ones(w), 'valid') / w


# Lista di score riportati in ultima colonna (prob unpairedness di intera finestra)
score_lst = [ float(val.strip()) for val in sys.stdin if val != "NA\n" ]

# Converto in nparray
score_nparr = np.array(score_lst)

# Calcolo mediana e MAD (Median Absolute Deviation from the median) 
median = np.median(score_nparr)
mad = np.median(np.absolute(score_nparr - median))

# Calcolo modified z-score (senza correzione per dati Normali) per ogni elemento
# Asse x parte dal centro della prima finestra di 10 nt
modif_zscore = (score_nparr - median) / mad
x = np.arange(5, len(score_nparr)+5)

# Media Mobile ordine 9
movAv_modif_zscore = moving_average(modif_zscore,9)
x2 = np.arange(5+4, len(score_nparr)+5-4)

# Unpairedness profile
plt.figure(figsize=(12, 4), dpi=150)
plt.bar(x, modif_zscore)
plt.title('SingleStrand Prob Profile')
plt.xlabel('sequence (nt)')
plt.ylabel('modified z-score')
plt.savefig('ss_prof.png')

# Unpairedness profile - Media Mobile ordine 10
sns.set_theme()
plt.figure(figsize=(12, 5), dpi=150)
plt.plot(x2, movAv_modif_zscore,linewidth=0.8)
plt.title('SingleStrand Prob Profile')
plt.xlabel('sequence (nt)')
plt.ylabel('modified z-score')
plt.savefig('ss_prof.mediaMobile.png')

# Stdout
print('# Median: {}'.format(median))
print('# Fine range: {}'.format(x[-1]))
print('# Max modified z-score position: {}'.format(np.where(modif_zscore==np.amax(modif_zscore))))
for i in modif_zscore:
	print("{:.4f}".format(i))
