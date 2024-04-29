library(here)

experiments <- list.dirs(path = here("data", "temp"), full.names = F, recursive = F)

for(experiment in experiments){
    load(here("data", "temp", experiment, "trial_num.RData"))
    
    complete_num <- length(list.files(path = here("data", "temp", experiment),
                                      pattern="^\\d+?\\.RData$"))
    running_num <- length(list.files(path = here("data", "temp", experiment, "progress"),
                                     pattern="^\\d+?$"))
    
    status <- paste0(experiment, ": ", complete_num, "/", trial_num, " completed")
    if(complete_num != trial_num){
        status <- paste0(status, ", ", running_num, " are running.")
    }
    
    print(status)
}






