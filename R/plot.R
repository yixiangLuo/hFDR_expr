library(here)
library(tidyverse)
library(latex2exp)
library(ggpubr)
library(magick)
library(grid)
# library(cowplot)

source(here("R", "utils.R"))


# hFDR (mean and quantile) v.s. true FDR
draw_hFDR <- function(experiment, fig_var, band_mass = 0.95,
                      method_names, method_colors, method_shapes){
  method_names <- c("hFDR", "FDP")
  method_colors <- unname(method_colors[method_names])
  method_shapes <- unname(method_shapes[method_names])
  sum_stat <- function(data){
    digits <- 3
    paste(round(mean(data), digits),
          round(quantile(data, probs = (1-band_mass)/2, names = F), digits),
          round(quantile(data, probs = (1+band_mass)/2, names = F), digits))
  }
  extract_value <- function(data, index){
    sapply(data, function(dat){
      as.numeric(str_split(dat, " ")[[1]][index])
    })
  }
  factor_to_num <- function(factors){
    sapply(factors, function(factor){
      as.numeric(toString(factor))
    })
  }
  
  load(here("data", paste0(experiment, ".RData")))
  
  plot_list <- lapply(1:length(results), function(iter){
    result <- results[[iter]]
    
    result$hFDR <- result$hFDR %>%
      mutate(across(all_of(method_names), factor_to_num))
    
    hFDR <- result$hFDR %>%
      group_by(tune_index) %>%
      summarise(across(all_of(method_names), sum_stat)) %>%
      pivot_longer(-tune_index, names_to = "method", values_to = "summary") %>%
      mutate(mean = extract_value(summary, 1),
             low = extract_value(summary, 2),
             up = extract_value(summary, 3),
             log_lambda = sqrt(result$tune_seq[tune_index])) %>%
      select(-summary)
    
    typeIIerror <- result$hFDR %>%
      group_by(tune_index) %>%
      summarise(typeII_error = mean(typeII_error)) %>%
      mutate(log_lambda = sqrt(result$tune_seq[tune_index]),
             method = "typeII_error")
    
    y_range <- c(0, max(hFDR$up))
    bar.width <- 0.05 * abs(diff(range(hFDR$log_lambda)))
    plot <- ggplot(hFDR, aes(x = log_lambda, y = mean, color = method)) +
      geom_line(data = typeIIerror, aes(x = log_lambda, y = typeII_error), linetype = "dashed", color = "grey") +
      geom_line() +
      geom_errorbar(aes(x = log_lambda, y = mean, ymin = low, ymax = up,
                        color = method), width = bar.width,
                    position = position_dodge(width = bar.width)) +
      # geom_point(aes(shape = method), size = 2) +
      scale_color_manual(values = method_colors, labels = method_names, breaks = method_names, name = "Estimator") +
      scale_shape_manual(values = method_shapes, labels = method_names, breaks = method_names, name = "Estimator") +
      ylim(y_range) +
      theme_bw() +
      theme(aspect.ratio = 1,
            panel.grid = element_blank(),
            strip.text = element_text(size = 15),
            axis.title = element_text(size = 13),
            axis.text = element_text(size = 10),
            legend.position = "right",
            legend.title=element_text(size=9),
            legend.text=element_text(size=9)) +
      labs(x = "log(lambda)", y = "estimated/true FDR"
           # , title = fig_var$value[iter]
      )
    # ggsave(filename = here("figs", paste0(experiment, "-", fig_var$value[iter], ".pdf")),
    #        plot, width = 3.5, height = 3)
    
    return(plot)
  })
  
  n_figs <- length(plot_list)
  fig <- ggarrange(plotlist = plot_list, 
                   ncol = min(2, n_figs), nrow = ceiling(n_figs/2),
                   labels = fig_var$value,
                   common.legend = T, legend = "right",
                   align = "hv")
  fig <- annotate_figure(fig, top = text_grob(paste0("subfigures: ", fig_var$name),
                                              face = "bold", size = 14))
  
  ggsave(filename = here("figs", paste0(experiment, ".pdf")),
         fig, width = 7, height = 3*(ceiling(n_figs/2))+1)
  
  return(fig)
}

draw_aggregate <- function(aggregate, band_mass = 0.9){
  method_colors <- c("#e41a1c", "black", "#010101", "grey", "darkorange", "dodgerblue4")
  method_lines <- c("solid", "solid", "solid", "dashed", "solid", "solid")
  method_names <- c("hFDR", "FDP", "FDR", "FPR", "std", "std_true")
  show_names <- c(TeX("$\\widehat{FDR}$"), "FDP", "FDR", "FPR", TeX("$\\widehat{s.e.}(\\widehat{FDR})$"), TeX("$s.e.(\\widehat{FDR})$"))
  subfig_name <- function(fig_vars){
    sapply(fig_vars, function(fig_var){
      if(fig_var == "Inv_Sparse") "Sparse Precision"
      else if(fig_var == "X_AR") "X-AR"
      else if(fig_var == "Coef_AR") "Coef-AR"
      else str_replace(fig_var, "_", " ")
    })
  }
  sel_method_names <- c("Lasso", "FS", "Logistic", "Graphical Lasso")
  names(sel_method_names) <- c("lasso", "FS", "logistic", "glasso")
  
  verticle <- sum(sapply(aggregate, length)) > 1
  
  sum_stat <- function(data){
    digits <- 3
    mean_val <- mean(data)
    std <- sqrt(var(data))
    low_quantile <- quantile(data, probs = (1-band_mass)/2, names = F)
    up_quantile <- quantile(data, probs = (1+band_mass)/2, names = F)
    low_outlier <- mean(data[data < low_quantile])
    low_outlier <- if(!is.na(low_outlier) & abs(low_outlier-mean_val) > 3*std) low_outlier else NA
    up_outlier <- mean(data[data > up_quantile])
    up_outlier <- if(!is.na(up_outlier) & abs(up_outlier-mean_val) > 3*std) up_outlier else NA
    
    paste(round(c(mean_val, low_quantile, up_quantile, low_outlier, up_outlier), digits), collapse = " ")
  }
  extract_value <- function(data, index){
    sapply(data, function(dat){
      char <- str_split(dat, " ")[[1]][index]
      if(char == "NA") -1 else unname(as.numeric(char))
    })
  }
  factor_to_num <- function(factors){
    sapply(factors, function(factor){
      as.numeric(toString(factor))
    })
  }
  
  
  
  figs <- lapply(names(aggregate), function(experiment){
    source(here("R", "settings", paste0(experiment, ".R")))
    n <- n_seq[[1]]
    
    words <- strsplit(experiment, "-")[[1]]
    model <- words[1]
    if(model == "glasso"){
      model <- "Gaussian graphic"
      sel_method <- "glasso"
    } else{
      sel_method <- words[2]
    }
    # robust <- words[length(words)] == "robust"
    subfigs <- aggregate[[experiment]]
    
    load(here("data", paste0(experiment, ".RData")))
    
    data <- lapply(subfigs, function(subfig){
      result <- results[[subfig]]
      # n <- result$side_info[[1]]$n
      
      result$hFDR <- result$hFDR %>%
        mutate(across(all_of(c("hFDR", "FDP")), factor_to_num))
      
      hFDR <- result$hFDR %>%
        group_by(tune_index) %>%
        summarise(across(all_of(c("hFDR", "FDP")), sum_stat)) %>%
        pivot_longer(-tune_index, names_to = "method", values_to = "summary") %>%
        mutate(mean = extract_value(summary, 1),
               low = extract_value(summary, 2),
               up = extract_value(summary, 3),
               low_out = extract_value(summary, 4),
               up_out = extract_value(summary, 5),
               tuning_para = result$tune_seq[tune_index],
               subfig = subfig_name(subfig),
               type = "FDR") %>%
        select(-summary)
      hFDR$tuning_para <- log(hFDR$tuning_para * ifelse(sel_method %in% c("FS", "glasso"), 1, n)) * ifelse(sel_method == "FS", 1, -1)
      hFDR$tuning_para <- hFDR$tuning_para * (1-2*(hFDR$subfig == "FS"))
      
      typeIIerror <- result$hFDR %>%
        group_by(tune_index) %>%
        summarise(typeII_error = mean(typeII_error)) %>%
        mutate(tuning_para = result$tune_seq[tune_index],
               method = "FPR",
               subfig = subfig_name(subfig),
               type = "FDR")
      typeIIerror$tuning_para <- log(typeIIerror$tuning_para * ifelse(sel_method %in% c("FS", "glasso"), 1, n)) * ifelse(sel_method == "FS", 1, -1)
      typeIIerror$tuning_para <- typeIIerror$tuning_para * (1-2*(typeIIerror$subfig == "FS"))
      
      n_tune <- max(result$hFDR$tune_index)
      sample_size <- max(result$hFDR$sample_id)
      
      hFDR.samples <- matrix(NA, n_tune, sample_size)
      FDP.samples <- matrix(NA, n_tune, sample_size)
      
      for(sample_i in 1:sample_size){
        hFDR.samples[, sample_i] <- factor_to_num((result$hFDR %>% filter(sample_id == sample_i))$hFDR)
        FDP.samples[, sample_i] <- factor_to_num((result$hFDR %>% filter(sample_id == sample_i))$FDP)
      }
      FDR.samples <- rowMeans(FDP.samples)
      deviate.true <- calc_deviate(hFDR.samples, FDR.samples, "std")[[2]]
      
      deviate.est <- result$hFDR %>%
        mutate(std = extract_dev(std, index = 2)) %>%
        select(-hFDR, -FDP) %>%
        group_by(tune_index) %>%
        summarise(std = sum_stat(std)) %>%
        pivot_longer(-tune_index, names_to = "method", values_to = "summary") %>%
        mutate(mean = extract_value(summary, 1),
               low = extract_value(summary, 2),
               up = extract_value(summary, 3),
               low_out = extract_value(summary, 4),
               up_out = extract_value(summary, 5),
               tuning_para = result$tune_seq[tune_index],
               subfig = subfig_name(subfig),
               type = "s.e.") %>%
        select(-summary)
      deviate.est$tuning_para <- log(deviate.est$tuning_para * ifelse(sel_method %in% c("FS", "glasso"), 1, n)) * ifelse(sel_method == "FS", 1, -1)
      deviate.est$tuning_para <- deviate.est$tuning_para * (1-2*(deviate.est$subfig == "FS"))
      
      deviate.true <- data.frame(tuning_para = result$tune_seq,
                                 dev = deviate.true,
                                 method = "std_true",
                                 subfig = subfig_name(subfig),
                                 type = "s.e.")
      deviate.true$tuning_para <- log(deviate.true$tuning_para * ifelse(sel_method %in% c("FS", "glasso"), 1, n)) * ifelse(sel_method == "FS", 1, -1)
      deviate.true$tuning_para <- deviate.true$tuning_para * (1-2*(deviate.true$subfig == "FS"))
      
      return(list(hFDR = rbind(hFDR, deviate.est), typeIIerror = typeIIerror, deviate.true = deviate.true))
    })
    
    hFDR <- do.call(rbind, lapply(data, function(df_list){df_list$hFDR}))
    FDR <- hFDR %>% filter(method == "FDP") %>% mutate(method = "FDR") %>% select(-c("low", "up"))
    typeIIerror <- do.call(rbind, lapply(data, function(df_list){df_list$typeIIerror}))
    deviate.true <- do.call(rbind, lapply(data, function(df_list){df_list$deviate.true}))
    
    
    y_range <- c(0, max(hFDR$up))
    hFDR <- hFDR %>% 
      group_by(subfig, type) %>%
      mutate(bar_width = 0.05 * abs(diff(range(tuning_para))))
    # https://stackoverflow.com/questions/53490654/adding-the-errorbar-icon-in-legend-in-ggplot
    GeomErrorbar$draw_key <- function(data, params, size){
      if(!(data$colour %in% c("grey", "dodgerblue4", "#010101"))){
        # lwd = data$size * .pt
        data$linetype[is.na(data$linetype)] <- 0
        segmentsGrob(c(0.4, 0.4, 0.5), c(0.2, 0.8, 0.2), c(0.6, 0.6, 0.5), c(0.2, 0.8, 0.8),
                           gp = gpar(col = alpha(data$colour, data$alpha),
                                           lwd = 0.8, lty = data$linetype, lineend = "butt"),
                           arrow = params$arrow) 
      }
    }
    GeomLine$draw_key <- function(data, params, size){
      if(!(data$colour %in% c("black"))){
        if (is.null(data$linetype)) { data$linetype <- 0 }
        else { data$linetype[is.na(data$linetype)] <- 0 }
        segmentsGrob(0.1, 0.5, 0.9, 0.5,
                           gp = gpar(col = alpha(data$colour %||% data$fill %||% "black", data$alpha),
                                           fill = alpha(params$arrow.fill %||% data$colour %||% data$fill %||% "black", data$alpha), 
                                           lwd = (data$linewidth %||% 0.5) * .pt, lty = data$linetype %||% 1,
                                           lineend = params$lineend %||% "butt"), arrow = params$arrow)
      }
    }
    GeomPoint$draw_key <- function(data, params, size){
      if(!(data$colour %in% c("grey", "dodgerblue4", "#010101"))){
        pointsGrob(c(0.5, 0.5), c(0.05, 0.95), gp = gpar(col = alpha(data$colour %||% "black", data$alpha), fill = alpha(data$fill %||% "black", data$alpha), fontsize = 1, lwd = 1))
      }
    }
    plot <- ggplot(hFDR, aes(x = tuning_para, y = mean, color = method, linetype = method)) +
      geom_line(data = typeIIerror, aes(x = tuning_para, y = typeII_error)) +
      geom_line() +
      geom_line(data = FDR, aes(x = tuning_para, y = mean)) +
      geom_errorbar(aes(x = tuning_para, y = mean, ymin = low, ymax = up,
                        color = method, width = bar_width), 
                    size = 0.5, position = "dodge") +
      geom_line(data = deviate.true, aes(y = dev, color = method)) +
      # geom_point(aes(x = tuning_para, y = low_out, color = method), size = 0.2) +
      # geom_point(aes(x = tuning_para, y = up_out, color = method), size = 0.2) +
      scale_color_manual(values = method_colors, labels = show_names, breaks = method_names, name = "value") +
      scale_linetype_manual(values = method_lines, labels = show_names, breaks = method_names, name = "value") +
      ylim(c(0, NA)) +
      theme_bw() +
      theme(aspect.ratio = 1,
            panel.grid = element_blank(),
            strip.text = element_text(size = 10),
            strip.text.y = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 10),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 9, margin = margin(r = 10, unit = "pt")),
            legend.text.align = 0) +
      guides(colour = guide_legend(nrow = 1)) +
      labs(x = ifelse(sel_method == "FS", "log(# steps)", TeX("$-log(\\lambda)$")), y = "") +
      ggtitle(sel_method_names[sel_method])
    
    if(verticle){
      plot <- plot +
        facet_grid(vars(factor(type, levels = c("FDR", "s.e."))),
                   vars(factor(subfig, levels = subfig_name(subfigs))), scales = "free") 
    } else{
      plot <- plot +
        ggh4x::facet_grid2(cols = vars(factor(type, levels = c("FDR", "s.e."))), scales = "free", independent = "y")
    }
    
    return(plot)
  })
  
  fig <- ggarrange(plotlist = figs, ncol = length(figs), nrow = 1,
                   common.legend = T, legend = "bottom") +
    theme(plot.margin = margin(0, 0, 0, 0, "cm")) 
  fig_cols <- sum(sapply(aggregate, length))
  
  if(verticle){
    fig <- annotate_figure(fig, left = "          s.e.                                FDR")
    fig_width <- max(6, 2 * fig_cols + 2.5)
    fig_height <- 5
  } else{
    fig_width <- 6
    fig_height <- 3
  }
  
  models <- sapply(names(aggregate), function(experiment){
    strsplit(experiment, "-")[[1]][1]
  })
  
  if("MCC" %in% do.call(c, aggregate)){
    fig_name <- "unfavorable"
  } else{
    aggregate_names <- paste0(names(aggregate), collapse = "-")
    aggregate_names <- unique(strsplit(aggregate_names, "-")[[1]])
    fig_name <- paste0(aggregate_names, collapse = "-")
  }
  
  ggsave(filename = here("figs", paste0(fig_name, ".pdf")),
         fig, width = fig_width, height = fig_height)
  
return(fig)
}

draw_hFDR_aggre <- function(experiment, sel_methods, X_types,
                            band_mass = 0.9, method_names, method_colors, method_shapes){
  method_names <- c("hFDR", "FDP")
  method_colors <- unname(method_colors[method_names])
  method_shapes <- unname(method_shapes[method_names])
  sum_stat <- function(data){
    digits <- 3
    paste(round(mean(data), digits),
          round(quantile(data, probs = (1-band_mass)/2, names = F), digits),
          round(quantile(data, probs = (1+band_mass)/2, names = F), digits))
  }
  extract_value <- function(data, index){
    sapply(data, function(dat){
      as.numeric(str_split(dat, " ")[[1]][index])
    })
  }
  factor_to_num <- function(factors){
    sapply(factors, function(factor){
      as.numeric(toString(factor))
    })
  }
  
  data <- lapply(sel_methods, function(sel_method){
    load(here("data", paste0(experiment, sel_method, ".RData")))
    
    data.sel <- lapply(X_types, function(X_type){
      result <- results[[X_type]]
      
      result$hFDR <- result$hFDR %>%
        mutate(across(all_of(method_names), factor_to_num))
      
      hFDR <- result$hFDR %>%
        group_by(tune_index) %>%
        summarise(across(all_of(method_names), sum_stat)) %>%
        pivot_longer(-tune_index, names_to = "method", values_to = "summary") %>%
        mutate(mean = extract_value(summary, 1),
               low = extract_value(summary, 2),
               up = extract_value(summary, 3),
               tuning_para = result$tune_seq[tune_index],
               X_type = X_type) %>%
        select(-summary)
      # if(sel_method != "FS") hFDR$tuning_para <- -log(hFDR$tuning_para)
      # else hFDR$tuning_para <- sqrt(hFDR$tuning_para)
      hFDR$tuning_para <- log(hFDR$tuning_para) * ifelse(sel_method == "FS", 1, -1)
      
      typeIIerror <- result$hFDR %>%
        group_by(tune_index) %>%
        summarise(typeII_error = mean(typeII_error)) %>%
        mutate(tuning_para = result$tune_seq[tune_index],
               method = "typeII_error",
               X_type = X_type)
      typeIIerror$tuning_para <- log(typeIIerror$tuning_para) * ifelse(sel_method == "FS", 1, -1)
      # if(sel_method != "FS") typeIIerror$tuning_para <- -log(typeIIerror$tuning_para)
      # else typeIIerror$tuning_para <- sqrt(typeIIerror$tuning_para)
      
      return(list(hFDR = hFDR, typeIIerror = typeIIerror))
    })
    
    hFDR <- do.call(rbind, lapply(data.sel, function(df_list){df_list$hFDR})) %>%
      mutate(sel_method = sel_method)
    typeIIerror <- do.call(rbind, lapply(data.sel, function(df_list){df_list$typeIIerror})) %>%
      mutate(sel_method = sel_method)
    
    return(list(hFDR = hFDR, typeIIerror = typeIIerror))
  })
  
  hFDR <- do.call(rbind, lapply(data, function(df_list){df_list$hFDR}))
  typeIIerror <- do.call(rbind, lapply(data, function(df_list){df_list$typeIIerror}))
  
  y_range <- c(0, max(hFDR$up))
  hFDR <- hFDR %>% 
    group_by(X_type, sel_method) %>%
    mutate(bar_width = 0.05 * abs(diff(range(tuning_para))))
  plot <- ggplot(hFDR, aes(x = tuning_para, y = mean, color = method)) +
    geom_line(data = typeIIerror, aes(x = tuning_para, y = typeII_error), linetype = "dashed", color = "grey") +
    geom_line() +
    geom_errorbar(aes(x = tuning_para, y = mean, ymin = low, ymax = up,
                      color = method, width = bar_width), 
                  position = "dodge") +
    # geom_point(aes(shape = method), size = 2) +
    scale_color_manual(values = method_colors, labels = method_names, breaks = method_names, name = "Estimator") +
    scale_shape_manual(values = method_shapes, labels = method_names, breaks = method_names, name = "Estimator") +
    ylim(y_range) +
    theme_bw() +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          strip.text = element_text(size = 15),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 10),
          legend.position = "right",
          legend.title=element_text(size=9),
          legend.text=element_text(size=9)) +
    labs(x = "tuning parameter", y = "estimated/true FDR")
  plot <- plot +
    if(length(unique(hFDR$X_type)) == 1 && length(unique(hFDR$sel_method)) > 1){
      fig_width <- 3*(ceiling(length(sel_methods)))+1
      fig_height <- 3*(ceiling(length(X_types)))
      ggh4x::facet_grid2(vars(factor(X_type, levels = X_types)),
                         vars(factor(sel_method, levels = sel_methods)), scales = "free")
    } else{
      fig_width <- 3*(ceiling(length(X_types)))+1
      fig_height <- 3*(ceiling(length(sel_methods)))
      ggh4x::facet_grid2(vars(factor(sel_method, levels = sel_methods)),
                         vars(factor(X_type, levels = X_types)), scales = "free", independent = "x")
    }
  
  
  ggsave(filename = here("figs", paste0(experiment, "hFDR-", paste(X_types, collapse = "_"), ".pdf")),
         plot, width = fig_width, height = fig_height)
  
  return(plot)
}

draw_hFDR_asymp <- function(experiment, sel_methods, X_types, p_seq, at_tune_index = 1){
  method_names <- c("hFDR", "FDP")
  sum_stat <- function(data){
    digits <- 3
    paste(round(mean(data), digits),
          round(sqrt(var(data)), digits))
  }
  extract_value <- function(data, index){
    sapply(data, function(dat){
      as.numeric(str_split(dat, " ")[[1]][index])
    })
  }
  factor_to_num <- function(factors){
    sapply(factors, function(factor){
      as.numeric(toString(factor))
    })
  }
  
  data <- lapply(sel_methods, function(sel_method){
    load(here("data", paste0(experiment, "-", sel_method, "-asymp.RData")))
    
    data.sel <- lapply(X_types, function(X_type){
      data.X <- lapply(p_seq, function(p){
        result <- results[[paste0(X_type, "-", p)]]
        
        result$hFDR <- result$hFDR %>%
          mutate(across(all_of(method_names), factor_to_num))
        
        hFDR <- result$hFDR %>%
          filter(tune_index == at_tune_index) %>%
          mutate(bias = mean(hFDR) - mean(FDP),
                 std = sqrt(var(hFDR)),
                 X_type = X_type,
                 p = p) %>%
          select(bias, std, X_type, p)
        
        return(hFDR)
      })
      hFDR <- do.call(rbind, data.X)
      return(hFDR)
    })
    
    hFDR <- do.call(rbind, data.sel) %>% mutate(sel_method = sel_method)
    return(hFDR)
  })
  
  hFDR <- do.call(rbind, data) %>%
    pivot_longer(c("bias", "std"), names_to = "moment", values_to = "value")
  
  plot <- ggplot(hFDR, aes(x = p, y = value, color = moment)) +
    geom_line() +
    facet_grid(vars(factor(sel_method, levels = sel_methods)),
               vars(factor(X_type, levels = X_types))) +
    theme_bw() +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          strip.text = element_text(size = 15),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 10),
          legend.position = "right",
          legend.title=element_text(size=9),
          legend.text=element_text(size=9)) +
    labs(x = "number of variables", y = "Bias and std of hFDR")
  
  ggsave(filename = here("figs", paste0(experiment, "hFDR", "-asymp", ".pdf")),
         plot, width = 7, height = 3*(ceiling(length(sel_methods)))+1)
  
  return(plot)
}


# hFDR - FDP(mean and quantile) v.s. FDP - FDP = 0
draw_FDP_diff <- function(experiment, fig_var, band_mass = 0.95,
                          method_names, method_colors, method_shapes){
  method_names <- c(method_names, "FDP")
  method_colors <- unname(method_colors[method_names])
  method_shapes <- unname(method_shapes[method_names])
  sum_stat <- function(data){
    digits <- 3
    paste(round(mean(data), digits),
          round(quantile(data, probs = (1-band_mass)/2, names = F), digits),
          round(quantile(data, probs = (1+band_mass)/2, names = F), digits))
  }
  extract_value <- function(data, index){
    sapply(data, function(dat){
      as.numeric(str_split(dat, " ")[[1]][index])
    })
  }
  
  load(here("data", paste0(experiment, ".RData")))
  
  plot_list <- lapply(results, function(result){
    oracle_FDP <- result$hFDR$FDP
    hFDR <- result$hFDR %>%
      mutate_at(all_of(method_names), ~ (.x - oracle_FDP))
    hFDR <- hFDR %>%
      group_by(tune_index) %>%
      summarise(across(all_of(method_names), sum_stat)) %>%
      pivot_longer(-tune_index, names_to = "method", values_to = "summary") %>%
      mutate(mean = extract_value(summary, 1),
             low = extract_value(summary, 2),
             up = extract_value(summary, 3),
             log_lambda = log(result$tune_seq[tune_index])) %>%
      select(-summary)
    
    plot <- ggplot(hFDR, aes(x = log_lambda, y = mean, color = method)) +
      geom_line() +
      geom_errorbar(aes(x = log_lambda, y = mean, ymin = low, ymax = up,
                        color = method), width = 0.05,
                    position = position_dodge(width = 0.05)) +
      # geom_point(aes(shape = method), size = 2) +
      scale_color_manual(values = method_colors, labels = method_names, breaks = method_names) +
      scale_shape_manual(values = method_shapes, labels = method_names, breaks = method_names) +
      theme_bw() +
      theme(aspect.ratio = 1,
            panel.grid = element_blank(),
            strip.text = element_text(size = 15),
            axis.title = element_text(size = 13),
            axis.text = element_text(size = 10),
            legend.position = "right",
            legend.title=element_text(size=9),
            legend.text=element_text(size=9)) +
      labs(x = "log(lambda)", y = "Difference with true FDP")
    
    return(plot)
  })
  
  n_figs <- length(plot_list)
  fig <- ggarrange(plotlist = plot_list, 
                   ncol = min(2, n_figs), nrow = ceiling(n_figs/2),
                   labels = fig_var$value,
                   common.legend = T, legend = "right",
                   align = "hv")
  fig <- annotate_figure(fig, top = text_grob(paste0("subfigures: ", fig_var$name),
                                              face = "bold", size = 14))
  
  ggsave(filename = here("figs", paste0(experiment, "-fdp-diff.pdf")),
         fig, width = 7, height = 3*(ceiling(n_figs/2))+1)
  
  return(fig)
}

# samples of hFDR v.s. true FDR
draw_est <- function(experiment, fig_var, sample = 1,
                     method_names, method_colors, method_shapes){
  method_colors <- c(unname(method_colors[method_names]), "black")
  method_shapes <- c(unname(method_shapes[method_names]), "1")
  method_names <- c(method_names, "true_FDR")
  # method_colors <- c(unname(method_colors[method_names]))
  # method_shapes <- c(unname(method_shapes[method_names]))
  # method_names <- c(method_names, "FDP")
  
  
  load(here("data", paste0(experiment, ".RData")))
  
  plot_list <- lapply(results, function(result){
    FDR <- result$hFDR %>%
      group_by(tune_index) %>%
      summarise(true_FDR = mean(FDP))
    
    hFDR <- result$hFDR %>%
      # filter(sample_id == sample) %>%
      cbind(FDR) %>%
      pivot_longer(all_of(method_names), names_to = "method", values_to = "hfdr") %>%
      mutate(log_lambda = log(result$tune_seq[tune_index]))
    
    plot <- ggplot(hFDR, aes(x = log_lambda, y = hfdr, color = method, alpha = method)) +
      geom_line(aes(group = interaction(sample_id, method))) +
      # geom_point() +
      # geom_point(aes(shape = method), size = 2) +
      scale_color_manual(values = method_colors, labels = method_names, breaks = method_names) +
      scale_alpha_manual(values = c(0.1, 1), labels = method_names, breaks = method_names) +
      theme_bw() +
      theme(aspect.ratio = 1,
            panel.grid = element_blank(),
            strip.text = element_text(size = 15),
            axis.title = element_text(size = 13),
            axis.text = element_text(size = 10),
            legend.position = "right",
            legend.title=element_text(size=9),
            legend.text=element_text(size=9)) +
      labs(x = "log(lambda)", y = "estimated/true hFDR")
    
    return(plot)
  })
  
  n_figs <- length(plot_list)
  fig <- ggarrange(plotlist = plot_list, 
                   ncol = min(2, n_figs), nrow = ceiling(n_figs/2),
                   labels = fig_var$value,
                   common.legend = T, legend = "right",
                   align = "hv")
  fig <- annotate_figure(fig, top = text_grob(paste0("subfigures: ", fig_var$name),
                                              face = "bold", size = 14))
  
  ggsave(filename = here("figs", paste0(experiment, "-est.pdf")),
         fig, width = 7, height = 3*(ceiling(n_figs/2))+1)
  
  return(fig)
}

# estimated std of hFDR (mean and quantile) v.s. true std of hFDR
draw_dev <- function(experiment, fig_var, band_mass = 0.95,
                     method_names, method_colors){
  
  est.method_names <- method_names[method_names != "hFDR"]
  method_colors <- c("black", method_colors[method_names != "hFDR"])
  method_names <- c("std_true", method_names[method_names != "hFDR"])
  
  sum_stat <- function(data){
    digits <- 3
    paste(round(mean(data), digits),
          round(quantile(data, probs = (1-band_mass)/2, names = F), digits),
          round(quantile(data, probs = (1+band_mass)/2, names = F), digits))
  }
  extract_value <- function(data, index){
    sapply(data, function(dat){
      as.numeric(str_split(dat, " ")[[1]][index])
    })
  }
  factor_to_num <- function(factors){
    sapply(factors, function(factor){
      as.numeric(toString(factor))
    })
  }
  
  load(here("data", paste0(experiment, ".RData")))
  
  nlambda <- max(results[[1]]$hFDR$tune_index)
  sample_size <- max(results[[1]]$hFDR$sample_id)
  
  for(index in 1:2){
    plot_list <- lapply(results, function(result){
      hFDR.samples <- matrix(NA, nlambda, sample_size)
      FDP.samples <- matrix(NA, nlambda, sample_size)
      
      for(sample_i in 1:sample_size){
        hFDR.samples[, sample_i] <- factor_to_num((result$hFDR %>% filter(sample_id == sample_i))$hFDR)
        FDP.samples[, sample_i] <- factor_to_num((result$hFDR %>% filter(sample_id == sample_i))$FDP)
      }
      FDR.samples <- rowMeans(FDP.samples)
      deviate.true <- calc_deviate(hFDR.samples, FDR.samples, "std")[[index]]
      # deviate.true <- rowMeans(hFDR.samples * calc_deviate(hFDR.samples, "logstd")[[index]])
      y.range <- c(0, min(5, max(deviate.true)*1.5))
      
      data <- result$hFDR %>%
        mutate(across(all_of(est.method_names), ~ extract_dev(.x, index = index)))
      # mutate(across(all_of(est.method_names), ~ .x * factor_to_num(hFDR)))
      
      deviate.est <- data %>%
        select(-hFDR, -FDP) %>%
        group_by(tune_index) %>%
        summarise(across(all_of(est.method_names), sum_stat)) %>%
        pivot_longer(-tune_index, names_to = "method", values_to = "summary") %>%
        mutate(mean = extract_value(summary, 1),
               low = extract_value(summary, 2),
               up = extract_value(summary, 3),
               log_lambda = log(result$tune_seq[tune_index])) %>%
        select(-summary)
      
      deviate.true <- data.frame(log_lambda = log(result$tune_seq),
                                 dev = deviate.true,
                                 method = "std_true")
      bar.width <- 0.05 * abs(diff(range(log(result$tune_seq))))
      
      plot <- ggplot(deviate.est, aes(x = log_lambda)) +
        geom_line(data = deviate.true, aes(y = dev, color = method), size = 1) +
        geom_line(aes(y = mean, color = method)) +
        geom_errorbar(aes(y = mean, ymin = low, ymax = up, color = method),
                      width = bar.width, position = position_dodge(width = bar.width)) +
        scale_color_manual(values = method_colors, labels = method_names, breaks = method_names) +
        lims(y = y.range) +
        theme_bw() +
        theme(aspect.ratio = 1,
              panel.grid = element_blank(),
              strip.text = element_text(size = 15),
              axis.title = element_text(size = 13),
              axis.text = element_text(size = 10),
              legend.position = "right",
              legend.title=element_text(size=9),
              legend.text=element_text(size=9)) +
        labs(x = "log(lambda)", y = "estimated/true deviate")
      
      return(plot)
    })
    
    n_figs <- length(plot_list)
    fig <- ggarrange(plotlist = plot_list, 
                     ncol = min(2, n_figs), nrow = ceiling(n_figs/2),
                     labels = fig_var$value,
                     common.legend = T, legend = "right",
                     align = "hv")
    fig <- annotate_figure(fig, top = text_grob(paste0("subfigures: ", fig_var$name),
                                                face = "bold", size = 14))
    
    ggsave(filename = here("figs", paste0(experiment, "-dev", "-", c("low", "up")[index], ".pdf")),
           fig, width = 7, height = 3*(ceiling(n_figs/2))+1)
  }
  
  return(fig)
}

draw_dev_aggre <- function(experiment, sel_methods, X_types,
                           band_mass = 0.95, method_names, method_colors, method_shapes){
  est.method_names <- c("std")
  method_colors <- c("black", method_colors[est.method_names])
  method_names <- c("std_true", est.method_names)
  names(method_colors) <- method_names
  sum_stat <- function(data){
    digits <- 3
    paste(round(mean(data), digits),
          round(quantile(data, probs = (1-band_mass)/2, names = F), digits),
          round(quantile(data, probs = (1+band_mass)/2, names = F), digits))
  }
  extract_value <- function(data, index){
    sapply(data, function(dat){
      as.numeric(str_split(dat, " ")[[1]][index])
    })
  }
  factor_to_num <- function(factors){
    sapply(factors, function(factor){
      as.numeric(toString(factor))
    })
  }
  
  data <- lapply(sel_methods, function(sel_method){
    load(here("data", paste0(experiment, sel_method, ".RData")))
    
    n_tune <- max(results[[1]]$hFDR$tune_index)
    sample_size <- max(results[[1]]$hFDR$sample_id)
    index <- 2
    
    data.sel <- lapply(X_types, function(X_type){
      result <- results[[X_type]]
      
      hFDR.samples <- matrix(NA, n_tune, sample_size)
      FDP.samples <- matrix(NA, n_tune, sample_size)
      
      for(sample_i in 1:sample_size){
        hFDR.samples[, sample_i] <- factor_to_num((result$hFDR %>% filter(sample_id == sample_i))$hFDR)
        FDP.samples[, sample_i] <- factor_to_num((result$hFDR %>% filter(sample_id == sample_i))$FDP)
      }
      FDR.samples <- rowMeans(FDP.samples)
      deviate.true <- calc_deviate(hFDR.samples, FDR.samples, "std")[[index]]
      
      data <- result$hFDR %>%
        mutate(across(all_of(est.method_names), ~ extract_dev(.x, index = index)))
      
      deviate.est <- data %>%
        select(-hFDR, -FDP) %>%
        group_by(tune_index) %>%
        summarise(across(all_of(est.method_names), sum_stat)) %>%
        pivot_longer(-tune_index, names_to = "method", values_to = "summary") %>%
        mutate(mean = extract_value(summary, 1),
               low = extract_value(summary, 2),
               up = extract_value(summary, 3),
               tuning_para = result$tune_seq[tune_index],
               X_type = X_type) %>%
        select(-summary)
      # if(sel_method != "FS") deviate.est$tuning_para <- -log(deviate.est$tuning_para)
      # else deviate.est$tuning_para <- sqrt(deviate.est$tuning_para)
      deviate.est$tuning_para <- log(deviate.est$tuning_para) * ifelse(sel_method == "FS", 1, -1)
      
      deviate.true <- data.frame(tuning_para = result$tune_seq,
                                 dev = deviate.true,
                                 method = "std_true",
                                 X_type = X_type)
      deviate.true$tuning_para <- log(deviate.true$tuning_para) * ifelse(sel_method == "FS", 1, -1)
      # if(sel_method != "FS") deviate.true$tuning_para <- -log(deviate.true$tuning_para)
      # else deviate.true$tuning_para <- sqrt(deviate.true$tuning_para)
      
      return(list(deviate.est = deviate.est, deviate.true = deviate.true))
    })
    
    deviate.est <- do.call(rbind, lapply(data.sel, function(df_list){df_list$deviate.est})) %>%
      mutate(sel_method = sel_method)
    deviate.true <- do.call(rbind, lapply(data.sel, function(df_list){df_list$deviate.true})) %>%
      mutate(sel_method = sel_method)
    
    return(list(deviate.est = deviate.est, deviate.true = deviate.true))
  })
  
  deviate.est <- do.call(rbind, lapply(data, function(df_list){df_list$deviate.est}))
  deviate.true <- do.call(rbind, lapply(data, function(df_list){df_list$deviate.true}))
  
  deviate.est <- deviate.est %>% 
    group_by(X_type, sel_method) %>%
    mutate(bar_width = 0.05 * abs(diff(range(tuning_para))))
  
  plot <- ggplot(deviate.est, aes(x = tuning_para)) +
    geom_line(data = deviate.true, aes(y = dev, color = method), size = 1) +
    geom_line(aes(y = mean, color = method)) +
    geom_errorbar(aes(y = mean, ymin = low, ymax = up, color = method, width = bar_width),
                  position = "dodge") +
    scale_color_manual(values = method_colors, labels = c("true std", "estimation"), breaks = method_names) +
    theme_bw() +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          strip.text = element_text(size = 15),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 10),
          legend.position = "right",
          legend.title=element_text(size=9),
          legend.text=element_text(size=9)) +
    labs(x = "tuning parameter", y = "estimated/true deviate")
  
  plot <- plot +
    if(length(unique(deviate.est$X_type)) == 1 && length(unique(deviate.est$sel_method)) > 1){
      fig_width <- 3*(ceiling(length(sel_methods)))+1
      fig_height <- 3*(ceiling(length(X_types)))
      ggh4x::facet_grid2(vars(factor(X_type, levels = X_types)),
                         vars(factor(sel_method, levels = sel_methods)), scales = "free")
    } else{
      fig_width <- 3*(ceiling(length(X_types)))+1
      fig_height <- 3*(ceiling(length(sel_methods)))
      ggh4x::facet_grid2(vars(factor(sel_method, levels = sel_methods)),
                         vars(factor(X_type, levels = X_types)), scales = "free", independent = "x")
    }
  
  
  ggsave(filename = here("figs", paste0(experiment, "std-", paste(X_types, collapse = "_"), ".pdf")),
         plot, width = fig_width, height = fig_height)
  
  return(plot)
}


# samples of estimated std of hFDR v.s. true std of hFDR
draw_dev_est <- function(experiment, fig_var, lambda.index = 5, exaggerate = 1){
  
  load(here("data", paste0(experiment, ".RData")))
  
  plot_list <- lapply(results, function(result){
    hFDR.samples <- as.numeric((result$hFDR %>% filter(tune_index == lambda.index))$hFDR)
    FDP.samples <- as.numeric((result$hFDR %>% filter(tune_index == lambda.index))$FDP)
    
    deviate.true <- calc_sample_deviate(hFDR.samples, FDP.samples)
    deviate.est <- (result$hFDR %>% filter(tune_index == lambda.index))$std
    deviate.est <- extract_dev(deviate.est, 2)
    
    plot.data <- data.frame(deviate_true = deviate.true, deviate_est = exaggerate * deviate.est)
    
    plot <- ggplot(plot.data, aes(x = deviate_true, y = deviate_est)) +
      geom_point() +
      geom_abline(slope = 1, intercept = 0, color = "black") +
      # coord_fixed(ratio = 1) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            strip.text = element_text(size = 15),
            axis.title = element_text(size = 13),
            axis.text = element_text(size = 10),
            # aspect.ratio = 1,
            legend.position = "right",
            legend.title=element_text(size=9),
            legend.text=element_text(size=9)) +
      labs(x = "true deviate", y = "estimated deviate")
    
    return(plot)
  })
  
  n_figs <- length(plot_list)
  fig <- ggarrange(plotlist = plot_list, 
                   ncol = min(2, n_figs), nrow = ceiling(n_figs/2),
                   labels = fig_var$value,
                   common.legend = T, legend = "right",
                   align = "hv")
  fig <- annotate_figure(fig, top = text_grob(paste0("subfigures: ", fig_var$name),
                                              face = "bold", size = 14))
  
  ggsave(filename = here("figs", paste0(experiment, "-dev_pred.pdf")),
         fig, width = 7, height = 3*(ceiling(n_figs/2))+1)
  
  return(fig)
}

# hFDR +- std v.s. true FDP for individual problems
draw_predict <- function(experiment, fig_var, sample = 1, exaggerate = 1, save_fig = T){
  
  load(here("data", paste0(experiment, ".RData")))
  
  colors <- c("true_FDR" = "black", "hFDR" = "#e41a1c")
  
  plot_list <- lapply(results, function(result){
    
    hFDR <- result$hFDR %>%
      filter(sample_id == sample) %>%
      mutate(log_lambda = log(result$tune_seq[tune_index]),
             std = extract_dev(std, 2))
    
    FDR <- result$hFDR %>%
      group_by(tune_index) %>%
      summarise(true_FDR = mean(FDP)) %>%
      mutate(log_lambda = log(result$tune_seq[tune_index]))
    
    plot <- ggplot(hFDR, aes(x = log_lambda)) +
      geom_errorbar(aes(y = hFDR,
                        ymin = hFDR - exaggerate * std,
                        ymax = hFDR + exaggerate * std,
                        color = "hFDR"),
                    width = 0.05) +
      geom_line(aes(y = hFDR, color = "hFDR")) +
      geom_line(data = FDR, aes(y = true_FDR, color = "true_FDR")) +
      scale_color_manual(values = colors) + 
      theme_bw() +
      theme(aspect.ratio = 1,
            panel.grid = element_blank(),
            strip.text = element_text(size = 15),
            axis.title = element_text(size = 13),
            axis.text = element_text(size = 10)) +
      labs(x = "log(lambda)", y = "estimated/true FDR",
           color = "value")
    
    return(plot)
  })
  
  n_figs <- length(plot_list)
  fig <- ggarrange(plotlist = plot_list, 
                   ncol = min(2, n_figs), nrow = ceiling(n_figs/2),
                   labels = fig_var$value,
                   common.legend = T, legend = "right",
                   align = "hv")
  fig <- annotate_figure(fig, top = text_grob(paste0("subfigures: ", fig_var$name, "\n",
                                                     "sample.id: ", sample),
                                              face = "bold", size = 14))
  if(save_fig){
    ggsave(filename = here("figs", paste0(experiment, "-one_pred.pdf")),
           fig, width = 7, height = 3*(ceiling(n_figs/2))+1)
  }
  
  return(fig)
}

# gif of draw_predict for each individual problem
draw_predict_all <- function(experiment, fig_var, sample_size, exaggerate = 1){
  for(sample_id in 1:sample_size){
    print(sample_id)
    fig <- draw_predict(experiment, fig_var, sample = sample_id, exaggerate, save_fig = F)
    ggsave(filename = here("figs", "temp", paste0(experiment, "-one_pred-",
                                                  int_to_str(sample_id, nchar(sample_size)), ".png")),
           fig, width = 800, height = 800, units = "px", dpi = 150, bg = "white")
  }
  ### list all png files
  png_files <- list.files(here("figs", "temp"),
                          pattern = paste0(experiment, "-one_pred-"),
                          recursive = FALSE,
                          all.files = FALSE,
                          full.names = TRUE)
  
  ### create a GIF file from all the plots
  png_files %>%
    map(image_read) %>% # reads each path file
    image_join() %>% # joins image
    image_animate(fps = 1) %>% # animates
    image_write(here("figs", paste0(experiment, "-all_pred.gif")))
  
  for(sample_id in 1:sample_size){
    file.remove(here("figs", "temp", paste0(experiment, "-one_pred-",
                                            int_to_str(sample_id, nchar(sample_size)), ".png")))
  }
  
  
  # img <- image_graph(800, 800, res = 150)
  # for(sample_id in 1:sample_size){
  #     print(sample_id)
  #     fig <- draw_predict(experiment, fig_var, sample = sample_id, exaggerate, save_fig = F)
  #     print(fig)
  # }
  # dev.off()
  # animation <- image_animate(img, fps = 1, optimize = TRUE)
  # browser()
  # image_write(animation, here("figs", paste0(experiment, "-all_pred.gif")))
  
  return(NULL)
}


## for developing, not used now

draw_std_samples  <- function(experiment, fig_var, lambda.index = 1){
  
  beta_proj <- function(beta, beta_true){
    nonnull.bs <- abs(beta) >= max(beta_true)*0.5
    nonnull <- beta_true != 0
    
    # sum(nonnull.bs)
    # sum(abs(beta))
    # sum(abs(beta - beta_true))
    sum(abs(beta[nonnull.bs]))
    # sum(abs(beta[nonnull] - beta_true[nonnull]))
  }
  
  load(here("data", paste0(experiment, ".RData")))
  
  plot_list <- lapply(results, function(result){
    std_true_beta <- sqrt(var(filter(result$hFDR, tune_index == lambda.index)$hFDR))
    std_bs <- filter(result$hFDR, tune_index == lambda.index)$dev_est
    
    beta_index_true <- beta_proj(result$side_info[[1]]$beta, result$side_info[[1]]$beta)
    beta_index_bs <- sapply(result$side_info, function(side_info){
      beta_proj(side_info$beta_hat, side_info$beta)
    })
    data.bs <- data.frame(beta = beta_index_bs, std = std_bs)
    data.true <- data.frame(beta = beta_index_true, std = std_true_beta)
    
    sample_id = 1
    beta = result$side_info[[sample_id]]$beta
    print(max(beta))
    beta_hat = result$side_info[[sample_id]]$beta_hat
    beta_hat.ols = result$side_info[[sample_id]]$beta_hat.ols
    beta_hat.lasso = result$side_info[[sample_id]]$beta_hat.lasso
    compare = data.frame(beta = beta, beta_hat = beta_hat,
                         beta_hat.ols = beta_hat.ols, beta_hat.lasso = beta_hat.lasso)
    # browser()
    
    plot <- ggplot(data.bs, aes(x = beta, y = std)) +
      geom_point() +
      geom_point(data = data.true, color = "red")
    theme_bw() +
      theme(aspect.ratio = 1,
            panel.grid = element_blank(),
            strip.text = element_text(size = 15),
            axis.title = element_text(size = 13),
            axis.text = element_text(size = 10),
            legend.position = "right",
            legend.title=element_text(size=9),
            legend.text=element_text(size=9)) +
      labs(x = "beta L1 norm", y = "std")
    
    return(plot)
  })
  
  n_figs <- length(plot_list)
  fig <- ggarrange(plotlist = plot_list, 
                   ncol = min(2, n_figs), nrow = ceiling(n_figs/2),
                   labels = fig_var$value,
                   common.legend = T, legend = "right",
                   align = "hv")
  fig <- annotate_figure(fig, top = text_grob(paste0("subfigures: ", fig_var$name),
                                              face = "bold", size = 14))
  
  ggsave(filename = here("figs", paste0(experiment, "-std_samples.pdf")),
         fig, width = 7, height = 3*(ceiling(n_figs/2))+1)
  
  return(fig)
}


draw_beta_sample <- function(experiment, fig_var, lambda.index = 5, fig_var_index = 1){
  
  load(here("data", paste0(experiment, ".RData")))
  
  plot.names <- c("std_hFDR", "FDR", "FDP", "hFDR")
  
  plot_list <- lapply(plot.names, function(plot.name){
    result <- results[[fig_var_index]]
    
    plot.data <- result$hFDR %>% 
      filter(tune_index == lambda.index, sample_id > 2)
    
    cur_value <- filter(result$hFDR, tune_index == lambda.index, sample_id == 2)[[plot.name]]
    
    plot <- ggplot(plot.data, aes(x = beta_dist, y = !!as.name(plot.name))) +
      geom_point() +
      geom_hline(yintercept = cur_value, color = "blue") + 
      theme_bw() +
      theme(panel.grid = element_blank(),
            strip.text = element_text(size = 15),
            axis.title = element_text(size = 13),
            axis.text = element_text(size = 10),
            aspect.ratio = 1,
            legend.position = "right",
            legend.title=element_text(size=9),
            legend.text=element_text(size=9)) +
      labs(x = "mah_dis", y = plot.name) # P(mah_dis <= x)
    
    true_value <- filter(result$hFDR, tune_index == lambda.index, sample_id == 1)
    if(true_value$beta_dist > 1){
      data.true <- data.frame(beta_dist = true_value$beta_dist, value = true_value[[plot.name]])
      plot <- plot + geom_point(data = data.true, aes(x = beta_dist, y = value), color = "red")
    }
    
    return(plot)
  })
  
  n_figs <- length(plot_list)
  fig <- ggarrange(plotlist = plot_list, 
                   ncol = min(2, n_figs), nrow = ceiling(n_figs/2),
                   # labels = plot.names,
                   common.legend = T, legend = "right",
                   align = "hv")
  fig <- annotate_figure(fig, top = text_grob(paste(fig_var$name, fig_var$value[fig_var_index]),
                                              face = "bold", size = 14))
  
  ggsave(filename = here("figs", paste0(experiment, "-beta.pdf")),
         fig, width = 7, height = 3*(ceiling(n_figs/2))+1)
  
  return(fig)
}
