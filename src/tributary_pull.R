library(httr)

reach = '710825777'
pass_start_date = '20231215'

download_tributary = function(reach_id){
  url = 'https://geoglows.ecmwf.int/api/v2/retrospective/'
  website = paste0(url, reach_id, '?format=csv&start_date=', pass_start_date)
  data=content(GET(website))
  return(data)
}

testing = download_tributary(reach)
