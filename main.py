# -*- coding: utf-8 -*-

import pandas as pd
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
from lifelines.statistics import logrank_test


df = pd.read_csv('survivalInput.txt', sep='\t')

expression_median = df['expression'].median()
append = []
for i, exp in enumerate(df['expression']):
    if df.iloc[i]['expression'] > expression_median:
        append.append('high expression')
    else:
        append.append('low expression')
df['gene'] = append

f1 = df.gene == 'high expression'
T1 = df[f1]['futime']
C1 = df[f1]['fustat']

f2 = df.gene == 'low expression'
T2 = df[f2]['futime']
C2 = df[f2]['fustat']

kmf = KaplanMeierFitter()
ax = plt.subplot(111)
kmf.fit(T1, event_observed=C1, label=['MLK4 low expression'])
kmf.plot(ax=ax, linestyle='--', dashes=(5, 2), ci_show=False)
kmf.fit(T2, event_observed=C2, label=['MLK4 high expression'])
kmf.plot(ax=ax, linestyle='--', dashes=(5, 2), ci_show=False)
kmf2 = plt.gcf()
results = logrank_test(T1, T2, C1, C2, alpha=.999)

plt.title('survival curve (p = ' + str(results.p_value) + ')')
plt.xlabel('time (day)')
plt.ylabel('survival rate')
plt.show()











