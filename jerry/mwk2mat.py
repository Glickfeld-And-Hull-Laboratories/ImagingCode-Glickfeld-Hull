"""
Convert .mwk2 file to .mat file
"""

import sqlite3
import zlib
import msgpack
import pandas as pd
import numpy as np
import json 
from scipy.io import savemat



##################
##################

# Class for extracting mworks data from .mwk2 files (sqlite database)

##################
##################

class MWK2Reader:

    _compressed_text_type_code = 1
    _compressed_msgpack_stream_type_code = 2

    def __init__(self, filename):
        self._conn = sqlite3.connect(filename)

    def close(self):
        self._conn.close()

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        self.close()

    @staticmethod
    def _decompress(data):
        return zlib.decompress(data, -15)

    def __iter__(self):
        unpacker = msgpack.Unpacker(strict_map_key=False)
        for code, time, data in self._conn.execute('SELECT * FROM events'):
            if not isinstance(data, bytes):
                yield (code, time, data)
            else:
                unpacker.feed(data)
                for obj in unpacker:
                    if isinstance(obj, msgpack.ExtType):
                        if obj.code == self._compressed_text_type_code:
                            obj = self._decompress(obj.data).decode('utf-8')
                        elif (obj.code ==
                              self._compressed_msgpack_stream_type_code):
                            unpacker.feed(self._decompress(obj.data))
                            continue
                    yield (code, time, obj)



####################     
####################


# Get data from mwk2 file

# Code notes (source https://mworks.github.io/documentation/0.9/reference/sysvars.html#system-event-codes)

# -- 0 RESERVED_CODEC_CODE
# -- 1 RESERVED_SYSTEM_EVENT_CODE
# -- 2 RESERVED_COMPONENT_CODEC_CODE
# -- 3 RESERVED_TERMINATION_CODE

####################
####################

def get_codec_events_subject(file) -> dict:
    """
    Gets codec from each .mwk2 files

    Params:
    file: file path
    """ 
    codes = []
    times = []
    data = []
    with MWK2Reader(file) as event_file:
        for code, time, data_pack in event_file:
            codes.append(code)
            times.append(time)  
            data.append(data_pack)  

    #Create df #Create data frame
    df = pd.DataFrame({'event_code': codes, 'times':times, 'data':data})

    #Grab codec
    codec = df.query("event_code == 0")['data'][0]
    codec_formatted = {codec.get(x).get('tagname'):x for x in codec.keys()}

    #formate codec to data frame
    codec_df = pd.DataFrame({'event_name': codec_formatted.keys(), 'event_code':codec_formatted.values()})

    #get subjectNum code id
    subnum_event_code = codec_df.query("event_name == 'subjectNum'")['event_code'].reset_index(drop = True).iloc[0]

    subnum = df.query("event_code == @subnum_event_code")['data'].unique()

    if len(subnum) > 1:
        res = 'Error - more than 1 subject id found!!!'
        return res
    else:
        subnum_str = "i" + str(subnum[0])

    #event data
    event_df = df.query("event_code != 0 and event_code != 1 and event_code != 2")

    return {'codec': codec_df, 'event':event_df, 'subject': subnum_str, 'file':file}




####################

# Convert to .mat

####################
#Path of.mwk2
mwk2_path = "G:\\Behavior\\Data\\hullglick-VisStimRet-20240801-091221.mwk2"

#Extract data
res = get_codec_events_subject(mwk2_path)

codec = res.get('codec')

events = res.get('event')

#Join codec with events
event_codec = (
    events
        .merge(codec,
        how = 'left', 
        left_on = 'event_code',
        right_on='event_code'
        )
    .filter(items=['event_name', 'times', 'data'])
    .query("~event_name.str.contains('#')",  engine='python')
) 


#Perform conversion
mat_data = {}
for col in event_codec.columns:
    mat_data[col] = event_codec[col].to_numpy()

#Add subject
mat_data['subject'] = res.get('subject')

#add file name
mat_data['file'] = res.get('file')

#Path output
opath = res.get('file')[:-1-4] + '_mwk2mat' + '.mat'

#Save file
savemat(opath, mat_data)
