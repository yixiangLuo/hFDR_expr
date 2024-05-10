library(here)
library(tidyverse)
library(latex2exp)
library(ggpubr)
library(magick)
library(grid)
# library(cowplot)

source(here("R", "utils.R"))


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
  sel_method_names <- c("Lasso", "fs", "Logistic", "Graphical Lasso")
  names(sel_method_names) <- c("lasso", "fs", "logistic", "glasso")
  
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
    if(model == "test"){
      model <- "test"
      sel_method <- "test"
    } else if(model == "glasso"){
      model <- "Gaussian graphic"
      sel_method <- "glasso"
    } else {
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
      hFDR$tuning_para <- log(hFDR$tuning_para * ifelse(sel_method %in% c("fs", "glasso"), 1, n)) * ifelse(sel_method == "fs", 1, -1)
      hFDR$tuning_para <- hFDR$tuning_para * (1-2*(hFDR$subfig == "fs"))
      
      typeIIerror <- result$hFDR %>%
        group_by(tune_index) %>%
        summarise(typeII_error = mean(typeII_error)) %>%
        mutate(tuning_para = result$tune_seq[tune_index],
               method = "FPR",
               subfig = subfig_name(subfig),
               type = "FDR")
      typeIIerror$tuning_para <- log(typeIIerror$tuning_para * ifelse(sel_method %in% c("fs", "glasso"), 1, n)) * ifelse(sel_method == "fs", 1, -1)
      typeIIerror$tuning_para <- typeIIerror$tuning_para * (1-2*(typeIIerror$subfig == "fs"))
      
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
      deviate.est$tuning_para <- log(deviate.est$tuning_para * ifelse(sel_method %in% c("fs", "glasso"), 1, n)) * ifelse(sel_method == "fs", 1, -1)
      deviate.est$tuning_para <- deviate.est$tuning_para * (1-2*(deviate.est$subfig == "fs"))
      
      deviate.true <- data.frame(tuning_para = result$tune_seq,
                                 dev = deviate.true,
                                 method = "std_true",
                                 subfig = subfig_name(subfig),
                                 type = "s.e.")
      deviate.true$tuning_para <- log(deviate.true$tuning_para * ifelse(sel_method %in% c("fs", "glasso"), 1, n)) * ifelse(sel_method == "fs", 1, -1)
      deviate.true$tuning_para <- deviate.true$tuning_para * (1-2*(deviate.true$subfig == "fs"))
      
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
      labs(x = ifelse(sel_method == "fs", "log(# steps)", TeX("$-log(\\lambda)$")), y = "") +
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

plot.hFDR <- function(hFDR.obj, cv.obj, errors = NA, sign.tune = -1, xlab = "-log(lambda)", ylab = NA,
                      show_cv = T, show_hFDR = T, show_FPR = F, show_FDR = T, show_FDP = show_FDR,
                      log_cv = T, log_x = T, show_legend = T, show_axes = T, show_extra_axes = show_axes, 
                      FDP.leg = "FDP", leg_bty="o", main = NULL){
  if(is.na(ylab)) ylab <- expression("False Discovery Rate")
  
  trans_x <- if(log_x) log else identity
  
  plot.range <- range(c(hFDR.obj$hFDR-hFDR.obj$hFDR.se, hFDR.obj$hFDR+hFDR.obj$hFDR.se))
  plot.range[1] <- max(plot.range[1], 0)
  plot.range[2] <- min(plot.range[2], 1)
  plot.args <- list(x = sign.tune * trans_x(hFDR.obj$lambda),
                    y = hFDR.obj$hFDR,
                    xlab = "", ylab = "",
                    xlim = range(sign.tune * trans_x(cv.obj$lambda)),
                    ylim = plot.range,
                    type = "n")
  
  # "estimated FDR and scaled CV MSE"
  # new.args <- list(...)
  # if(length(new.args)) plot.args[names(new.args)] <- new.args
  
  # use log scale for MSE
  ep <- 1e-8
  mse_scale <- if(log_cv) function(x) log(pmax(x,ep)) else identity
  mse.mean <- mse_scale(cv.obj$cvm)
  mse.low <- mse_scale(cv.obj$cvlo)
  mse.up <- mse_scale(cv.obj$cvup)
  cv.range <- range(c(mse.low, mse.up))
  cv.transfer <- function(mse){
    (mse - min(cv.range)) / abs(diff(cv.range)) * abs(diff(plot.range)) + min(plot.range)
  }
  
  leg <- c()
  lty <- c()
  lwd <- c()
  pch <- c()
  col <- c()
  
  do.call("plot", plot.args)
  
  if(show_cv){
    lines(x = sign.tune * trans_x(cv.obj$lambda),
          y = cv.transfer(mse.mean),
          col = "dodgerblue3")
    error.bars(sign.tune * trans_x(cv.obj$lambda),
               cv.transfer(mse.up), cv.transfer(mse.low),
               width = 0.01, col = alpha("dodgerblue3", 0.1))
    
    abline(v = sign.tune * trans_x(cv.obj$lambda.min), lty = 1, col = "dodgerblue3")
    abline(v = sign.tune * trans_x(cv.obj$lambda.1se), lty = 2, col = "dodgerblue3")
    
    # leg <- c(leg, "-loglikelihood")
    leg <- c(leg, "CV MSE")
    lty <- c(lty, 1)
    lwd <- c(lwd, 1)
    pch <- c(pch, 26)
    col <- c(col, "dodgerblue3")
  }
  
  if(!any(is.na(errors))){
    if(show_FPR){
      lines(x = sign.tune * trans_x(cv.obj$lambda),
            y = errors$FPP,
            col = "orange3", lty = 2)
      lines(x = sign.tune * trans_x(cv.obj$lambda),
            y = errors$FPR,
            col = "orange3")
      
      leg <- c(leg, "Type II error", "Type II error rate")
      lty <- c(lty, 2, 1)
      lwd <- c(lwd, 1, 1)
      pch <- c(pch, 26, 26)
      col <- c(col, "orange3", "orange3")
    }
    
    if(!is.null(errors$FDP) && show_FDP){
      lines(x = sign.tune * trans_x(cv.obj$lambda),
            y = errors$FDP,
            col = "black", lty = 2)
      
      leg <- c(leg, FDP.leg)
      lty <- c(lty, 2)
      lwd <- c(lwd, 1)
      pch <- c(pch, 26)
      col <- c(col, "black")
    }
    if(!is.null(errors$FDR) && show_FDR){
      lines(x = sign.tune * trans_x(cv.obj$lambda),
            y = errors$FDR,
            col = "black")
      
      leg <- c(leg, "FDR")
      lty <- c(lty, 1)
      lwd <- c(lwd, 1)
      pch <- c(pch, 26)
      col <- c(col, "black")
    }
  }
  
  
  if(show_hFDR){
    error.bars(sign.tune * trans_x(hFDR.obj$lambda),
               hFDR.obj$hFDR+hFDR.obj$hFDR.se, hFDR.obj$hFDR-hFDR.obj$hFDR.se,
               width = 0.01, alpha("red", 0.3))
    points(sign.tune*trans_x(hFDR.obj$lambda), hFDR.obj$hFDR,
           pch = 20, col = "red")
    
    leg <- c(leg, TeX("$\\widehat{FDR}$"))
    lty <- c(lty, 1)
    lwd <- c(lwd, 0)
    pch <- c(pch, 19)
    col <- c(col, "red")
  }
  
  if(show_extra_axes) {
    axis(side = 3, at = sign.tune*trans_x(cv.obj$lambda),
         labels = paste(cv.obj$nzero), tick = FALSE, line = -0.5)
    mtext("# selections", side = 3, line = 2)
    
    cv.ticks <- if(log_cv){ seq(exp(cv.range[1]), exp(cv.range[2]), length.out = 6) }
    else { seq(cv.range[1], cv.range[2], length.out = 6) }
    axis(side = 4, at = cv.transfer(mse_scale(cv.ticks)),
         labels = formatC(cv.ticks, format = "g", digits = 2), tick = T, line = 0)
    mtext("CV MSE", side = 4, line = 2.5)  # "-log likelihood (cv)"
  }
  
  legend.pos <- legend("bottomright", inset = 0.05,
                       legend = leg, bg = "white",
                       lty = lty, lwd = lwd,
                       pch = pch, col = col, plot = F, y.intersp=1.2)[[1]]
  cap <- min((hFDR.obj$hFDR-hFDR.obj$hFDR.se)[sign.tune * trans_x(hFDR.obj$lambda) >= legend.pos$left])
  if(!any(is.na(errors))) cap <- min(cap, errors$FDP[sign.tune * trans_x(cv.obj$lambda) >= legend.pos$left])
  floor <- max(cv.transfer(mse.up[sign.tune * trans_x(cv.obj$lambda) >= legend.pos$left]))
  top <- if(cap > floor) (legend.pos$top + cap) / 2 else legend.pos$top
  if(show_legend) {
    legend(x = legend.pos$left, y = top, inset = 0.05,
           legend = leg, bg = "white",
           lty = lty, lwd = lwd,
           pch = pch, col = col, y.intersp=1.2, bty=leg_bty)
  }
  
  #par(mar = c(4,4,3,4), cex = 1, cex.axis = 1, cex.lab = 1)
  if(show_extra_axes) {
    title(xlab = xlab, ylab = ylab, line = 2.5, cex.lab = 1)
  } else {
    title(xlab = xlab, ylab = ylab, line = 2.5, cex.lab = 1)
    title(main = main, line = 1, cex.lab = 1)
  }
  invisible()
}