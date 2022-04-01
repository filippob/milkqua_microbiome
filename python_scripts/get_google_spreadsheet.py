#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 15:43:11 2022

@author: filippo
"""

#%% import libraries
import os
import gspread
import pandas as pd
from oauth2client.service_account import ServiceAccountCredentials

#%% set params
print(" - setting parameters")
scope = ['https://www.googleapis.com/auth/spreadsheets',
        'https://www.googleapis.com/auth/drive']
docid = "12qfJKcg3hHhG3210pim878cTvMSGFsmR0YET_qFs55E"
home = '/home/users/filippo.biscarini.est/MILKQUA'
credential_path = 'Config/microbiomes-6719db38cddb.json'
outdir = 'Config'
label = 'milkqua_stools_swabs.csv'

#%%
fname = os.path.join(home, credential_path)
credentials = ServiceAccountCredentials.from_json_keyfile_name(fname, scope)
client = gspread.authorize(credentials)
google_sh = client.open_by_key(docid)

#%%
sheet1 = google_sh.get_worksheet(0)

#%%
df = pd.DataFrame(data=sheet1.get_all_records())

#%% write out
fname = os.path.join(home, outdir, label)
print(" - writing to file {}".format(fname))
df.to_csv(fname, index=False)


print("DONE!")
