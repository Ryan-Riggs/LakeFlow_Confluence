import numpy as np
import pandas as pd
import s3fs
import xarray
import datetime

# You may set anon to False if you have a credential file stored on your system but it is not necessary for this demonstration
# dropna was added 10/18/2024

def pull_tributary (reach_id, start_date):
    start_date = datetime.datetime.strptime(start_date, '%m-%d-%Y').date()
    end_date = datetime.date.today()
    sim_begin = datetime.date(1940,1,1)
    first_index=start_date-sim_begin
    second_index = end_date-sim_begin
    bucket_uri = 's3://geoglows-v2-retrospective/retrospective.zarr'
    region_name = 'us-west-2'
    s3 = s3fs.S3FileSystem(anon=True, client_kwargs=dict(region_name=region_name))
    s3store = s3fs.S3Map(root=bucket_uri, s3=s3, check=False)
    ds = xarray.open_zarr(s3store)
    df = ds['Qout'].sel(rivid=reach_id)[first_index.days:second_index.days,].to_dataframe()
    df = df.reset_index().set_index('time').pivot(columns='rivid', values='Qout')
    return(df)