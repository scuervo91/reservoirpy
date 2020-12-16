import numpy as np 
import pandas as pd 
import os 
from datetime import date

def schedule_writer(keywords_dict, keywords, start_date=None, end_date=None):
    string = ""
    string += "RPTRST\n 'BASIC=1' /\n"
    string += "RPTSCHED\n 'FIP=1' 'WELLS=1' 'WELLSPECS' /\n"
    dates_list = []
    
    #Extract all dates from all keywords
    for d in keywords_dict:
        _dates_df = keywords_dict[d]['date']
        dates_list.append(_dates_df)
    
    dates_df = pd.concat(dates_list, axis=0).sort_values()

    # Flag to decide to write new DATE keyword
    write_date_keyword = True 
    for date in dates_df.unique():
        if start_date is not None:
            if pd.Timestamp(date) < start_date:
                continue
        if end_date is not None:
            if pd.Timestamp(date) > end_date:
                continue
        # If other keywords have been written in previous date 
        # The keyword Dates must be written
        if write_date_keyword:
            string += 'DATES\n'
        
        date_str = pd.Timestamp(date).strftime("%d '%b' %Y").upper()
        string += date_str + ' /\n'

        #Check if keyword is in the allowed keywords to write
        write_date_keyword = False
        #Iterate over dict of keywords
        c=0
        for key in keywords_dict:
            if key=='DATES':
                write_date_keyword = True if write_date_keyword else False
                continue               
            elif key not in keywords:
                write_date_keyword = True if write_date_keyword else False
                continue
            
            key_date = keywords_dict[key].loc[keywords_dict[key]['date']==pd.Timestamp(date)]
            if key_date.empty:

                write_date_keyword = True if write_date_keyword else False
                continue
            c += 1
            if c==1:
                string += '/\n'
                
            string += key + '\n'
            cols = [i for i in key_date.columns if i != 'date']
            for index, row in key_date.iterrows():
                string += key_date.loc[[index],cols].to_string(index=False, header=False) + '/\n'
            string += '/\n' 

            write_date_keyword = True
            
    return string
   
        
      
        
            
            
            
            
    
    
    
    
    
    