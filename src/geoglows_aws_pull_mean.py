import numpy as np
import pandas as pd
import s3fs
import xarray
import datetime
import geoglows
# You may set anon to False if you have a credential file stored on your system but it is not necessary for this demonstration

def pull_mean (reach_id):
    mean = geoglows.data.annual_averages(river_id=reach_id)
    return(mean)