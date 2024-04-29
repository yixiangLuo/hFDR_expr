library(here)

get_expr_index <- function(trial_id, index_data){
    expr_chunk_size <- index_data$sample_len
    
    expr_index <- (trial_id-1) %/% expr_chunk_size + 1
    
    sample_index <- (trial_id-1) %% expr_chunk_size + 1
    
    return(list(expr_index = expr_index, sample_index = sample_index,
                trial_id = trial_id))
}

get_data_indeces <- function(job_id, index_data){
    trial_num <- index_data$n_expr * index_data$sample_len
    trials_per_job <- ceiling(trial_num / index_data$n_jobs)
    
    job_trials <- ((job_id-1) * trials_per_job + 1) : (job_id * trials_per_job)
    job_trials <- job_trials[job_trials <= trial_num]
    
    job_trial_indeces <- lapply(job_trials, function(trial_id){
        get_expr_index(trial_id, index_data)
    })
    
    return(job_trial_indeces)
}

