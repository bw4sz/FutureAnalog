fnFAPlot <- function(arbthresh, comm.type, scenario.type, stat, high.col, low.col){
  out.rasters <- list.files(paste(out_path, arbthresh, sep="/"), 
                            pattern = "NonAnalogRasters", 
                            full.names = TRUE
  )
  
  NonAnalogRasters <- lapply(out.rasters, 
                             function(x) {
                               load(x)
                               return(results)
                             }
  )
  
  names(NonAnalogRasters) <- substr(out.rasters, 
                                    nchar(out.rasters) - 11, 
                                    nchar(out.rasters) - 8)
  
  mean_res <- list()
  mean_scenario <- list()
  
  if(scenario.type == "RCP"){
    scenario <- unique(substr(names(NonAnalogRasters), 3, 4))
  }
  if(scenario.type == "GCM"){
    scenario <- unique(substr(names(NonAnalogRasters), 1, 2))
  } 
  
  for(s in scenario) {
    rasters <- grep(s, names(NonAnalogRasters))
    rasters <- NonAnalogRasters[rasters]
    for(nam in names(rasters[[1]])){
      if(grepl(comm.type, nam) == TRUE) {
        res <- lapply(rasters, function(x) x[[nam]])
        res <- stack(res)
        mean_res[[nam]] <- calc(res, get(stat))
      }
    }
    mean_scenario[[s]] <- stack(mean_res)
  }
  
  plot_dat <- stack(mean_scenario)
  repWithScenario <- function(y){
    x <- substr(y, nchar(y), nchar(y))
    gsub(x, scenario[as.numeric(x)], y)
  }
  names(plot_dat) <- lapply(names(plot_dat), function(x) repWithScenario(x))
  
  gplot(plot_dat) + geom_tile(aes(fill = value)) +
    facet_wrap(~variable, nrow = 3) +
    scale_fill_gradient(low = low.col, high = high.col) +
    coord_equal() + 
    labs(x = "longitude", y = "latitude") + 
    scale_x_continuous(breaks = seq(-81, -75, 2)) + 
    theme_classic() + 
    theme(strip.background=element_blank(), panel.margin = unit(1, "lines"))
  
  filenam <- paste(arbthresh, comm.type, scenario.type, stat, sep="_", ".pdf")
  ggsave(paste(out_path, filenam, sep="/"), width = 17, height = 9)
}
