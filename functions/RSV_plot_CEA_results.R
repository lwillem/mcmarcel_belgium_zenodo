############################################################################ #
# This file is part of the RSV modelling project.
# 
# => VISUAL PRESENTATION OF THE RESULTS
#
#  Copyright 2025, CHERMID, UNIVERSITY OF ANTWERP
############################################################################ #

plot_CEA_results <- function(sim_output_filename)
{
  
  sim_output <- readRDS(paste0(sim_output_filename, '.rds'))
  dim(sim_output)
  
  # set pdf dimensions
  pdf_w_h <- c(2,2) 
  
  # get plot output directory
  plot_output_folder <- file.path(dirname(sim_output_filename), '/')
  
  # get all scenario-country-year combinations
  if(length(unique(sim_output$population_code)) > 1 || length(unique(sim_output$year)) > 1){
    sim_output$scen_country_year <- paste(sim_output$scenario,
                                          sim_output$population_code,
                                          sim_output$year,
                                          sep='_')
  } else{
    sim_output$scen_country_year <- sim_output$scenario
  }
  
  # get unique combinations
  scen_country_year_opt <- unique(sim_output$scen_country_year)
  
  # safety check on dimensions
  expected_num_records <- max(sim_output$scenario_id) * unique(sim_output$num_sim)
  if(nrow(sim_output) > expected_num_records){
    stop('safety stop in RSV_plot_CEA_results.R::plot_CEA_results() \n ==>>  The number of rows in "sim_output" is higher than expected based on the scenario IDs and number of simulations. Please verify that the model output contains only one resuts per stochastic realiation per scenario.')
  }
  
  # run all scenario-country_year combinations
  scen_country_year <- scen_country_year_opt[1]
  for(scen_country_year in scen_country_year_opt) 
    {
      # select subset 
      country_data <- sim_output[sim_output$scen_country_year == scen_country_year,]
      print(paste('Create country CEA plots for:',scen_country_year))
      
      # set scenario colours and legend
      country_data$intervention_factor <- as.factor(country_data$intervention)
      
      # select data and rename columns
      sel_country_data <- country_data[,c('incremental_cost_disc_program', 'total_QALY_disc_averted', 'intervention')]
      names(sel_country_data) <- c('Cost','Effect','Strategy')
      
      # align names with CEAC/CEAF graphs
      ce_legend      <- get_intervention_legend(sel_country_data$Strategy)      
      ce_legend      <- ce_legend[order(ce_legend$name, decreasing = TRUE),] 
      sel_country_data$Strategy <- replace_for(sel_country_data$Strategy, ce_legend$name, ce_legend$name_legend)
      
      # open pdf stream
      pdf(paste0(plot_output_folder, 'CE_PLANE_', scen_country_year, '.pdf'), width = 4.8 * pdf_w_h[1], height = 2.8 * pdf_w_h[2])
      
      # only average values
      print(plot.icers.edit(sel_country_data, 
                            plot_uncertainty='none',
                            xlab ='Discounted Effects',
                            ylab ='Discounted Costs',
                            currency = 'EUR') + 
              theme(legend.position = "bottom"))
      
      # including uncertainty for the frontier
      print(plot.icers.edit(sel_country_data, 
                            plot_uncertainty = 'frontier',
                            xlab ='Discounted Effects',
                            ylab ='Discounted Costs',
                            currency = 'EUR'))
      
      # including all uncertainty
      print(plot.icers.edit(sel_country_data, 
                            plot_uncertainty = 'all',
                            xlab ='Discounted Effects',
                            ylab ='Discounted Costs',
                            currency = 'EUR'))
      
      # close pdf stream
      dev.off()
      
    } -> dummy
  
  
}

get_intervention_legend <- function(factor_codes){
  
  factor_codes   <- as.factor(factor_codes)
  
  factor_legend  <- data.frame(name=levels(factor_codes),
                               color=seq(1,nlevels(factor_codes))+1,
                               stringsAsFactors = F)
  
  # fix color-scheme of the basecase interventions for original "price threshold code"
  basecase_scheme <- data.frame(name  = c("comparator","mAb","mAb_OctMar","mAb_OctMar_CB","maternal","maternal_SepMar"),
                                color = c("#8DA0CB","#E5C494","#FC8D62","#FFD92F","#A6D854","#B3B3B3"),
                                bool_basecase = TRUE,
                                stringsAsFactors = F)
  
  # identify names that are part of the baseline_scheme, and re-use these colors
  factor_legend$bool_basecase <- factor_legend$name %in% basecase_scheme$name
  factor_legend <- rbind(basecase_scheme[basecase_scheme$name %in% levels(factor_codes),],
                         factor_legend[!factor_legend$bool_basecase,])
  
  # include intervention type
  factor_legend$interventiontype <- rep(NA,nrow(factor_legend))
  factor_legend$interventiontype[grep("mAb", factor_legend$name)] <- 'mAb'
  factor_legend$interventiontype[grep("mat", factor_legend$name)] <- 'MV'
  factor_legend$interventiontype[grepl("mAb", factor_legend$name) & grepl("mat", factor_legend$name)] <- 'mAb+MV'
  factor_legend$interventiontype[grep("comparator", factor_legend$name)] <- 'comparator'
  
  # include lty and lwd
  factor_legend$lty <- 1
  factor_legend$lwd <- 3
  
  # set other colors (if any left)
  num_non_basecase <- sum(!factor_legend$bool_basecase)
  if(num_non_basecase>0){
    #use the most distinct colours for each intervention type
    qual_col_pals = brewer.pal.info
    col_vectorname='Set2'
    if(nrow(factor_legend)>8) { col_vectorname="Dark2" }
    col_vector=as.vector(unlist(mapply(brewer.pal,qual_col_pals$maxcolors[rownames(qual_col_pals)==col_vectorname],col_vectorname)))
    
    # remove the colors from the baseline (if any)
    col_vector <- col_vector[!col_vector%in% basecase_scheme$color]
    
    # manually remove indistinguishable color from Set2 and Dark2
    col_vector <- col_vector[!col_vector%in% '#1B9E77']
    col_vector <- rev(col_vector)
    
    # make sure the color vector is long enough, if not, duplicate
    if(num_non_basecase>length(col_vector)){
      incr_factor <- ceiling(nrow(factor_legend)/length(col_vector))
      lty_vector = rep(1:incr_factor,each=length(col_vector))
      col_vector = rep(col_vector,incr_factor)
      print("WARNING: DUPLICATED COLORS IN PLOT(S) DUE TO THE HIGH NUMBER OF INTERVENTIONS")
    } else{
      lty_vector <- factor_legend$lty
    }
    factor_legend$color[!factor_legend$bool_basecase] = col_vector[1:num_non_basecase]
    factor_legend$lty[!factor_legend$bool_basecase]   = lty_vector[1:num_non_basecase]
    
    # fix for light yellow color, make darker
    factor_legend$color <- gsub("#FFFFB3","#FFCC00",factor_legend$color)
  }
  
  # adjust legend names
  factor_legend$name_legend <- factor_legend$name
  factor_legend$name_legend[factor_legend$name == 'comparator']  <- 'no intervention'
  factor_legend$name_legend <- gsub('_CB', ' + catch-up',factor_legend$name_legend)
  factor_legend$name_legend <- gsub('_', ': ',factor_legend$name_legend)
  factor_legend$name_legend <- gsub('maternal', 'MV',factor_legend$name_legend)
  
  
  # update notation if seasonal
  for(i_row in which(grepl(':',factor_legend$name_legend))){
    months <- gsub('.*: ', '',factor_legend$name_legend[i_row])
    months_adjusted <- paste0(substr(months,0,3), '-',substr(months,4,nchar(months)))
    factor_legend$name_legend[i_row] <- paste0(factor_legend$interventiontype[i_row], ': ',
                                               months_adjusted)
  }
  factor_legend
  
  return(factor_legend)
}

# adjusted function from dampack package
plot.icers.edit <- function(x_all,
                       txtsize = 12,
                       currency = "$",
                       effect_units = "QALYs",
                       label = c("all", "frontier", "none"),  # changed order
                       label_max_char = NULL,
                       plot_frontier_only = FALSE,
                       alpha = 1,
                       n_x_ticks = 6,
                       n_y_ticks = 6,
                       xbreaks = NULL,
                       ybreaks = NULL,
                       xlim = NULL,
                       ylim = NULL,
                       xlab = NULL, # new
                       ylab = NULL, # new
                       xexpand = expansion(0.1),
                       yexpand = expansion(0.1),
                       max.iter = 20000,
                       plot_uncertainty = c("none", "frontier", "all"), # new
                       ...) {
  
  # make sure the function variables are not pre-used yet
  Cost <- Effect <- Status <- Strategy <- Inc_Cost <- Inc_Effect <- ICER <- NULL
  
  # get average cost and effect
  x_mean  <- aggregate(. ~ Strategy, data = x_all,mean)
  
  # include no intervention
  levels(x_mean$Strategy) <- c(x_mean$Strategy,'no intervention')
  x_mean <- x_mean[c(1, 1:nrow(x_mean)), ]
  x_mean$Strategy[1] <- 'no intervention'
  x_mean[1,2:3] <- 0
  
  # calculate ICERS
  x <- calculate_icers(cost = x_mean$Cost, 
                       effect = x_mean$Effect, 
                       strategies = x_mean$Strategy)
  
  
  
  if (ncol(x) > 7) {
    # reformat icers class object if uncertainty bounds are present
    x <- x %>%
      select(Strategy, Cost, Effect,
             Inc_Cost, Inc_Effect,
             ICER, Status)
  }
  
  # type checking
  label <- match.arg(label)
  plot_uncertainty <- match.arg(plot_uncertainty)
  
  # this is so non-dominated strategies are plotted last (on top)
  x <- arrange(x, Status)
  
  # change status text in data frame for plotting
  d_name <- "Dominated"
  ed_name <- "Weakly Dominated"
  nd_name <- "Efficient Frontier"
  
  status_expand <- c("D" = d_name, "ED" = ed_name,
                     "ND" = nd_name, "ref" = nd_name)
  x$Status <- factor(status_expand[x$Status], ordered = FALSE,
                     levels = c(d_name, ed_name, nd_name))
  
  # linetype
  plot_lines <- c("Dominated" = "blank",
                  "Weakly Dominated" = "blank",
                  "Efficient Frontier" = "solid")
  
  # names to refer to in aes
  stat_name <- "Status"
  strat_name <- "Strategy"
  eff_name <- "Effect"
  cost_name <- "Cost"
  
  # frontier only
  if (plot_frontier_only) {
    plt_data <- x[x$Status == nd_name, ]
  } else {
    plt_data <- x
  }
  
  # make plot
  icer_plot <- ggplot(plt_data, aes(x = !!sym(eff_name), y = !!sym(cost_name),
                                    shape = !!sym(stat_name)))
  
  # add all data?
  if(plot_uncertainty != 'none'){
    
    x_uncertainty <- x_all %>% left_join(x[,c('Strategy', 'Status')], by='Strategy')
    
    if(plot_uncertainty == 'frontier'){
      x_uncertainty <- x_uncertainty %>% filter(Status == 'Efficient Frontier')
    }
    
    icer_plot <- icer_plot + geom_point(data = x_uncertainty,
                           aes(y = Cost,
                               x = Effect,
                               color = Strategy),
                           shape = 4,
                           alpha = 0.2,
                           inherit.aes = FALSE) + # to prevent inheriting from icer_plot
      guides(color = guide_legend(override.aes = list(alpha = 1))) # use alpha 1 in the legend
  }
  
  # explicitly add axis
  icer_plot <-  icer_plot + geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
                            geom_hline(yintercept = 0, linetype = "dashed", color = "grey")
  
  # add mean
  icer_plot <-  icer_plot + geom_point(alpha = alpha, size = 2) +
    geom_line(aes(linetype = !!sym(stat_name), group = !!sym(stat_name))) +
    scale_linetype_manual(name = NULL, values = plot_lines) +
    scale_shape_discrete(name = NULL)
  
  # specify axis
  if(is.null(xlab)) { xlab <- "Effect" }
  if(is.null(ylab)) { ylab <- "Cost" }
  icer_plot <- icer_plot + labs(x = paste0(xlab, " (", effect_units, ")"),
                                y = paste0(ylab, " (", currency, ")"))
  
  # required to use the adjusted plot.icers function
  source('functions/lib/add_common_aes.R')  
  
  # set other details
  icer_plot <- add_common_aes(icer_plot, txtsize, col = "none",
                              continuous = c("x", "y"),
                              n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                              xbreaks = xbreaks, ybreaks = ybreaks,
                              xlim = xlim, ylim = ylim,
                              xexpand = xexpand, yexpand = yexpand)
  
  # labeling
  if (label != "none") {
    if (!is.null(label_max_char)) {
      plt_data[, strat_name] <- str_sub(plt_data[, strat_name],
                                        start = 1L, end = label_max_char)
    }
    if (label == "all") {
      lab_data <- plt_data
    }
    if (label == "frontier") {
      lab_data <- plt_data[plt_data$Status == nd_name, ]
    }
    
    icer_plot <- icer_plot +
      ggrepel::geom_label_repel(data = lab_data,
                       aes(label = !!sym(strat_name)),
                       size = 3,
                       show.legend = FALSE,
                       max.iter = max.iter,
                       direction = "both")
  }
  return(icer_plot)
}

# help function to run a gsub operation for all elements of a vector or list
replace_for <- function(x, pattern, replacement){
    if(length(pattern) == length(replacement)){
      for(i in 1:length(pattern)){
        if(any(x == pattern[i])){
          x[x == pattern[i]] <- replacement[i]
        }
      }
    } else{
      warning('replace_for not possible with different dimensions of the pattern and replacement')
    }
    
    return(x)
  }
  
