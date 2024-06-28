import numpy as np
import pandas as pd
import s3fs
import xarray

# You may set anon to False if you have a credential file stored on your system but it is not necessary for this demonstration

def pull_tributary (reach_id):
    bucket_uri = 's3://geoglows-v2-retrospective/retrospective.zarr'
    region_name = 'us-west-2'
    s3 = s3fs.S3FileSystem(anon=True, client_kwargs=dict(region_name=region_name))
    s3store = s3fs.S3Map(root=bucket_uri, s3=s3, check=False)
    ds = xarray.open_zarr(s3store)
    df = ds['Qout'].sel(rivid=reach_id).to_dataframe()
    df = df.reset_index().set_index('time').pivot(columns='rivid', values='Qout')
    return(df)