# -*- coding: utf-8 -*-
import pandas as pd
from lasio import LASFile
import os

class log(LASFile):

    def __init__(self, file_ref,find_mnemonics=False,**kwargs):
        """__init__ [summary]

        Parameters
        ----------
        file_ref : [type]
            [description]
        find_mnemonics : bool, optional
            [description], by default False
        """
        mnemonics_df = kwargs.pop('mnemonics', None)
        LASFile.__init__(self, file_ref=file_ref,
                         autodetect_encoding = True, **kwargs)
        
        if find_mnemonics == True:
            file_dir = os.path.dirname(__file__)
            mnemonics_path = os.path.join(file_dir,'mnemonics.csv')
            if mnemonics_df is None:
                try:
                    mnemonics = pd.read_csv(mnemonics_path, header=0)
                    for curve in self.curves:
                        for col in mnemonics.columns:
                            if curve['mnemonic'] in mnemonics[col].tolist():
                                print('Mnemonic: ',curve['mnemonic'],' => ',  col + '_' + curve['mnemonic'] )
                                curve.mnemonic = col + '_' + curve['mnemonic']
                                break
                except:
                    print("file mnemoics.csv not found. Provide one")
                    pass
            else:
                for curve in self.curves:
                    for col in mnemonics_df.columns:
                        if curve['mnemonic'] in mnemonics_df[col].tolist():
                            print('Mnemonic: ',curve['mnemonic'],' => ',  col + '_' + curve['mnemonic'] )
                            curve.mnemonic = col + '_' + curve['mnemonic']
                            break

            


    
