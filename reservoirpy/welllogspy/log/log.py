# -*- coding: utf-8 -*-
import pandas as pd
from lasio import LASFile
import os

class log(LASFile):
    def __init__(self, file_ref,find_mnemonics=False,**kwargs):
        LASFile.__init__(self, file_ref=file_ref,
                         autodetect_encoding = True, **kwargs)
        
        if find_mnemonics == True:
            file_dir = os.path.dirname(__file__)
            mnemonics_path = os.path.join(file_dir,'mnemonics.csv')
            mnemonics = pd.read_csv(mnemonics_path, header=0)
            
            for curve in self.curves:
                for col in mnemonics.columns:
                    if curve['mnemonic'] in mnemonics[col].tolist():
                        print('Mnemonic: ',curve['mnemonic'],' => ', col)
                        curve.type = col 
                        break

    
