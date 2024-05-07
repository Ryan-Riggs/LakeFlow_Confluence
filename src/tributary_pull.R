library(httr)
library(data.table)

reaches = c('710825777', '710722628')
pass_start_date = '20240401'

download_tributary = function(reach_id){
  url = 'https://geoglows.ecmwf.int/api/v2/retrospective/'
  website = paste0(url, reach_id, '?format=csv&start_date=', pass_start_date)
  data=content(GET(website))
  return(data.table(data))
}

tributary_flow = lapply(reaches, download_tributary)

tributary_flow_merge = Reduce(function(...) merge(..., all = TRUE), tributary_flow)

tributary_aggregated = tributary_flow_merge[,tributary_total:=rowSums(.SD), .SDcols=-1]
