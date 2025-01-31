## Code to conduct FR-BCA

## library(dplyr)
## library(magrittr)

## Inputs (list)
## - Model name(s)
## - Analysis parameters {delta, T, losses}
##    - Baseline
##    - Ranges for sensitivity analysis
## - Cost data for each model
## - Data from performance assessment for each model
##    - Format?

## Outputs
## - Table with baseline BCR + high/low values, for each parameter included in
## the sensitvity analysis
## - Columns:
##    - model/variant
##    - sensitivity parameter name
##    - bcr
##    - label (base, low, high)
## - Ready to print tables and plots?


## for each parameter in inputs$parameters$sensitivity:
## - update "base" list with {low, high} values for that parameter
## - compute the BCRs and add a column indicating the parameter that is varied

## So there will be a single function to compute pv_cost and pv_loss/pv_benefit
## The inputs will either be base "as is" or base with a sensitivity parameter varied
## There will be a wrapper function that calls each function iteratively for each sensitivity

## The other thing to note is this will be done one model at a time
## so the input tables need to be separated, so that there is only one status quo
## and all other rows are compared to that status quo
## There will be a wrapper that separates the models in the background,
## does all the computations for each model, and then assembles back together

###
## Purpose:
## Format strings with underscores for pretty printing figure labels
###
#' Format axis tick labels without hard-coding
#'
#' @description
#' Wrapper that conducts BCA for each set of models in list
#'
#' @importFrom stringr str_to_title
#' @export
#'
label_format <- function(s) {
  return(stringr::str_to_title(gsub("_", " ", factor(s))))
}


#' @export
preprocess_model <- function(eal, cost, p) {
    ## Purpose:
    ## join eal and cost tables
    ## compute total floor area (floor area * num_stories)
    ## separate models by height
    ##
    join_cols = c('system', 'model', 'num_stories', 'intervention', 'design_s', 'design_ns')
    names_structural = c("RC III", "RC IV", "backup-frame")
    models <- list()
    dat <- eal |>
      dplyr::left_join(cost, by=join_cols) |>
      dplyr::mutate(total_area=p$floor_area*num_stories) |>
      dplyr::mutate(design=paste(design_s, design_ns, sep="-")) |>
      dplyr::mutate(design=case_when(
                      design %in% "baseline-baseline" ~ "baseline",
                      design %in% "baseline-nsfr" ~ "nonstructural",
                      design %in% paste(names_structural, "baseline", sep="-") ~ "structural",
                      design %in% paste(names_structural, "nsfr", sep="-") ~ "full"))
    systems <- dat |> dplyr::distinct(system) |> pull()
    for (j in systems) {
      dat_j = dat |> dplyr::filter(system == j)
      stories <- dat_j |> dplyr::distinct(num_stories) |> pull()
      for (i in stories) {
        models[[j]][[paste(i, "story", sep="-")]] <- dat_j |>
          dplyr::filter(num_stories == i)
      }
    }
    return(models)
}


#' @export
preprocess_cost <- function(model) {
    ## Purpose:
    ## ignore rows without nonstructural costs
    return(model |> dplyr::filter(!is.na(c_ns)))
}


## PV(total cost) = structural cost + PV(nonstructural cost)
###
## Purpose:
## Calculate present value s, ns, and total costs
###
#' Present Value Cost
#'
#' @description
#' Compute the present value (PV) structural, nonstructural, and total construction cost
#'
#' @importFrom dplyr mutate
#'
#'
#' @param model A model
#' @param params A list of parameters
#'
#' @return Updated model table including PV(Cost)
#' @export
#'
pv_cost <- function(model, params) {
    p = params$parameters$base
    return(
        model |>
        ## filter out missing NS cost
        ## preprocess_cost() |>
        dplyr::mutate(
                 pv_res = construction - c_ns,
                 ## pv_s=c_s,
                 pv_ns = c_ns
                 ## pv_ns=c_ns*(1 + 1/(1+p$delta)^(p$T/2))
               ) |>
         dplyr::mutate(pv_total=pv_res+pv_ns)
    )
}

## compute cost delta
## pv cost formula:
## pv_cost = pv_s + pv_ns + (pv_ns / (1+delta))
## pv_s = structural cost
## pv_ns = nonstructural cost * (1-delta)^(2018-2011)
#' @export
pv_dcost <- function(model, params) {
    ## Purpose:
    ## Calculate present value cost deltas, relative to status quo
    ## NB: assumes cost table filtered to only include rows with ns cost
    return(
        model |>
        pv_cost(params) |>
        dplyr::mutate(
                   cost_diff=pv_total - pv_total[intervention == 0],
                   cost_delta=(pv_total/pv_total[intervention == 0]) - 1
                      )
    )
}

#' @export
#'
compute_loss <- function(loss_name, p, occ_t, fr_t) {
  ## function to compute loss based on loss_name and parameter values
  ## NOTE: For now, to allow including new loss categories, just have to add computation to this function!
  loss = p[['loss']][[loss_name]]
  ## loss_val = loss[['value']] * (loss[['time']]*occ_t + (1-loss[['time']])*fr_t)
  if (loss_name == 'loss_displacement') {
    ## loss_val = loss_val * p[['tenant']]
    loss_val = loss * p[['tenant']] * fr_t
  } else if (grepl('(business_income|value_added)', loss_name)) {
    ## loss_val = loss_val * (1-p[['recapture']])
    loss_val = loss * (1-p[['recapture']]) * fr_t
  } else if (loss_name == 'loss_rental_income') {
    loss_val = loss * occ_t
  } else {
    loss_val = NA
  }
  return(loss_val)
}

#' @export
#'
compute_supply_chain_loss <- function(p, business_income) {
  ## NOTE: if supply chain losses need to be reformulated, just have to redefine this function!
  if (length(p[['loss']][['loss_supply_chain']]) > 0) {
    loss_val = p[['loss']][['loss_supply_chain']] * business_income
  } else {
    loss_val = NA
  }
  return(loss_val)
}

#' @export
## loss formulas:
## total_area = (floor area) * number of stories
## displacement = displace_per_area * tenant_per_area * reocc_days * total_area
## bi = (1 - recapture) * bi_per_area * fr_days * total_area
## ri = ri_per_area * fr_days * total_area
pv_loss <- function(model, p) {
    ## Purpose:
    ## Calculate losses for each model
    ## NB: assumes total_area column has been created
    loss_names = names(p$loss)
  for (loss_name in loss_names[!grepl('supply_chain', loss_names)]) {
    model = model |>
      dplyr::rowwise() |>
      dplyr::mutate(UQ(loss_name):=compute_loss(loss_name, p, re_occupancy_time, functional_recovery_time)*total_area) |>
      dplyr::ungroup()
  }
    return(
        model |>
        dplyr::mutate(loss_supply_chain=compute_supply_chain_loss(p, loss_business_income))
        )
}

#' @importFrom stats uniroot
#'
#' @export
#'
f_npv <- function(t, cf, i) {
  sum( cf / (1 + i)^t )
}

#' @export
#'
f_irr <- function(t, cf) {
  ## TODO: introduce error handling
  tryCatch(
    uniroot(f = f_npv, interval = c(0, 1), t = t, cf = cf)$root,
    error = function(e) {
      ## suppress error message, as IRR cannot be computed for all rows (eg, if cost_delta == 0)
      return(NA)
    },
    warning = function(w) {
      message("A warning occurred:\n", w)
    }
  )
}


#' @export
## pv benefits formula:
## (avoided losses) * ( (1 - (1+delta)^(-T)) / delta)
pv_benefit <- function(model, params, label='base') {
    ## Purpose:
    ## Calculate present value avoided losses, relative to status quo
    join_cols = c('system', 'model', 'num_stories', 'intervention', 'design_s', 'design_ns')
    ## loss_cols = c('repair_cost', 'displacement', 'business_income', 'rental_income', 'sc')
    p = params$parameters$base
    loss_cols = c('repair_cost', names(p$loss))
    m = model |>
        ## NB: assumes total_area column has been created
        pv_loss(p) |>
        dplyr::select(all_of(c(join_cols, loss_cols))) |>
        dplyr::rowwise() |>
        dplyr::mutate(loss_total=sum(across(all_of(loss_cols)), na.rm=TRUE)) |>
        dplyr::ungroup() |>
        dplyr::mutate(eal_diff=loss_total[intervention == 0] - loss_total) |>
        dplyr::mutate(eal_delta=eal_diff/loss_total[intervention == 0]) |>
        dplyr::rowwise() |>
        dplyr::mutate(benefit=f_npv(t=seq(from=0, to=p$T), cf=rep(eal_diff, length=p$T+1), i=p$delta)) |>
        dplyr::ungroup()
    return(
        model |>
        dplyr::select(!repair_cost) |>
        dplyr::left_join(m, by=join_cols)
    )
}


#' @export
bcr <- function(model, params, label='base') {
    ## Purpose:
    ## Compute BCR and NPV
    ## calls pv_benefit and pv_dcost
    ## TODO: add column with label {base, high, low}
    p = params$parameters$base
    model <- pv_benefit(model, params)
    model <- pv_dcost(model, params)
    return(model |>
           dplyr::rowwise() |>
           dplyr::mutate(
                      bcr=benefit/cost_diff,
                      npv=benefit-cost_diff,
                      irr=f_irr(t=seq(from=0, to=p$T), cf=c(eal_diff-cost_diff, rep(eal_diff, length=p$T))),
                    label=label) |>
           dplyr::mutate(aroi=(npv/cost_diff)/p$T) |>
           dplyr::ungroup()
           )
}


#' @export
set_params <- function(params, param, bound='low') {
    ## Purpose:
    ## Reset baseline parameter to one of {low, high}
    p <- params
    if (grepl('^loss', param)) {
    ## loss_names = names(p[['parameters']][['sensitivity']][['loss']])
    ## for (loss_name in loss_names) {
      p[['parameters']][['base']][['loss']][[param]] <- p[['parameters']][['sensitivity']][['loss']][[param]][[bound]]
    ## }
  } else {
    p[['parameters']][['base']][[param]] <- p[['parameters']][['sensitivity']][[param]][[bound]]
  }
    return(p)
}


#' @export
sensitivity <- function(model, params) {
    ## Purpose:
    ## Compute BCR and NPV, using low/high values of parameters
    ## --- ##
    ## get list of sensitivity parameters
    s <- params[['parameters']][['sensitivity']]
    ## store calculations
    m <- list()
    ## iterate low/high
    for (hi_low in c('low', 'high')) {
        ## iterate over parameters
        for (n in names(s)) {
          if (n == 'loss') {
            ss <- s[['loss']]
            for (nn in names(ss)) {
              p <- set_params(params, nn, bound=hi_low)
              m[[paste(nn, hi_low, sep='-')]] <- bcr(model, p, label=hi_low) |>
                dplyr::mutate(parameter=nn)
            }
          } else {
            p <- set_params(params, n, bound=hi_low)
            m[[paste(n, hi_low, sep='-')]] <- bcr(model, p, label=hi_low) |>
              dplyr::mutate(parameter=n)
          }
        }
    }
    return(dplyr::bind_rows(m))
}


#' @export
bca <- function(model, params) {
    ## Purpose:
    ## Wrapper that computes (1) baseline BCA and (2) sensitivity analysis
    ## --- ##
    ## first pass baseline bca
    m_b <- bcr(model, params)
    ## second pass sensitivity
    m_s <- sensitivity(model, params)
    ## TODO: pivot_wider on bcr -> bcr_base, bcr_low, bcr_high
    return(dplyr::bind_rows(m_b, m_s))
}


###
## Purpose:
## Wrapper that conducts BCA for each set of models in list
###
#' Conduct FR-BCA with sensitivity analysis
#'
#' @description
#' Wrapper that conducts BCA for each set of models in list
#'
#' @importFrom dplyr bind_rows
#'
#'
#' @param eal Table of Expected Annualized Losses (EAL)
#' @param cost Table of Structural, Nonstructural, and Total Construction Costs
#' @param params List of analysis parameters
#'
#' @return Table including EAL, Costs, Benefits, FR-BCA Outputs (BCR and NPV)
#' @export
#'
frbca <- function(eal, cost, params) {
  output = list()
  models <- preprocess_model(eal, cost, params[['parameters']][['base']])
  systems <- names(models)
  for (i in systems) {
    ## o_i = models[[i]]
    ## stories = names(o_i)
    stories = names(models[[i]])
    o_i = lapply(stories, FUN=function(x) { bca(models[[i]][[x]], params) })
    ## for (j in stories) {
    ##   o_i[[j]] <- bca(models[[i]][[j]], params)
    ## }
    output[[i]] = dplyr::bind_rows(o_i)
  }
  ## TODO: filter out NaN as base case?
  ## TODO: filter out NA for missing cost?
  ## output <- dplyr::bind_rows(output)
  return(output)
}

#' @export
#'
#' @importFrom forcats fct_rev
#' @importFrom tidyr pivot_wider
#'
postprocess_bcr <- function(output, systems="RCMF", designs="nonstructural", stories=4, out_base=FALSE) {
  ## function to postprocess output for plotting sensitivity
  ## drop baseline-baseline, if it exists
  designs = designs[!grepl('baseline', designs)]
  ## create data frame for table/plot
  plot_df <- output |>
    dplyr::filter((system %in% systems) &
                  (design %in% designs) &
                  (num_stories %in% stories)) |>
    dplyr::select(model, system, num_stories, design, bcr, label, parameter)
  base <- plot_df |>
    dplyr::filter(label == 'base') |>
    dplyr::select(!c(label, parameter))
  if (out_base) {
    return(base)
  } else {
    sen <- plot_df |>
      dplyr::filter(label != 'base') |>
      tidyr::pivot_wider(names_from=label, values_from=bcr) |>
      dplyr::left_join(base, by=c('model', 'system', 'num_stories', 'design')) |>
      dplyr::rename(bcr_low=low, bcr_high=high) |>
      dplyr::mutate(
               design=factor(design, levels=designs),
               parameter=forcats::fct_rev(parameter))
    return(sen)
  }
}

###
## Purpose:
## Post-process data and generate BCR plots for multiple systems
###
#' Plot FR-BCA Outputs for Baseline Parameters
#'
#' @description
#' Plot for FR-BCA outputs, for fixed story height and structural system
#'
#' @importFrom dplyr filter select left_join rename
#' @importFrom tidyr pivot_wider
#' @import ggplot2
#'
#'
#' @param output Output from `frbca()`
#' @param n_floors Number of stories (for figure title)
#' @param system Name of structural system being plotted (default: "RCMF")
#'
#' @return Updated model table including PV(Cost)
#' @export
#'
plot_bcr <- function(output, systems="RCMF", designs="nonstructural", stories=4) {
  ## generate plot
  ## TODO: parameterize label (ie, title)
  label <- "BCRs for baseline and recovery-based designs"
  plot.base <- postprocess_bcr(output, systems, designs, stories, out_base=TRUE) |>
    dplyr::mutate(
           system=factor(system, levels=systems),
           num_stories=factor(num_stories, levels=stories),
           design=factor(design, levels=designs)
         ) |>
  ggplot(aes(x = design, y = bcr, fill = num_stories)) +
  geom_bar(
    stat='identity',
    width = 0.5,
    position = 'dodge') +
  geom_hline(yintercept=1, colour="red", linetype="dashed") +
  geom_hline(yintercept=0, colour="black", linetype="solid", linewidth=0.5) +
  facet_wrap(~system, scales="free") +
    labs(
      ## title = label,
      x = "Stories",
      y = "BCR") +
  ggthemes::theme_few() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 13)
    ) +
    labs(fill = "")
return(plot.base)
}

###
## Purpose:
## Post-process data and generate plot for sensitivity analysis
###
#' Plot FR-BCA Outputs with Sensitivity Analysis
#'
#' @description
#' Sensitivity analysis plot for FR-BCA outputs, for fixed story height and structural system
#'
#' @importFrom dplyr filter select left_join rename
#' @importFrom tidyr pivot_wider
#' @import ggplot2
#' @importFrom ggthemes scale_fill_colorblind theme_few
#'
#'
#' @param output Output from `frbca()`
#' @param n_floors Number of stories (for figure title)
#' @param system Name of structural system being plotted (default: "RCMF")
#'
#' @return Updated model table including PV(Cost)
#' @export
#'
plot_bcr_sensitivity <- function(output, systems="RCMF", designs="nonstructural", stories=4) {
  ## generate plot
  label_begin <- 'Sensitivity Analysis: Benefit-cost ratios for'
  label_end <- 'archetypes, relative to baseline ASCE 7-16 design.'
  ## label_param <- c("Delta", "Business Income", "Displacement", "T")
  label_param <- c("T", "L_{DC}", "L_{BI}", "\\delta")
  plot.sen <- postprocess_bcr(output, systems, designs, stories) |>
    dplyr::mutate(
           system=factor(system, levels=systems),
           num_stories=factor(num_stories, levels=stories),
           design=factor(design, levels=designs)
         ) |>
    ggplot2::ggplot() +
    ggplot2::geom_segment(aes(x=parameter, xend=parameter, y=bcr_low, yend=bcr_high),
                 linewidth = 5, colour = "red", alpha = 0.6) +
    ggplot2::geom_segment(aes(x=parameter, xend=parameter, y=bcr-0.001, yend=bcr+0.001),
                 linewidth = 5, colour = "black") +
    ggplot2::geom_hline(yintercept=1, colour="red", linetype="dashed") +
    ## ggplot2::scale_x_discrete(labels=label_param) +
    ## ggplot2::scale_x_discrete(labels=rev(label_param)) +
    ## ggplot2::scale_x_discrete(labels=lapply(sprintf(r'($%s$)', label_param), TeX)) +
    ggplot2::scale_x_discrete(labels=c(expression("T"), expression(L[DC]), expression(L[BI]), expression(delta))) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(system ~ design, ncol=2) +
    ggplot2::theme_light() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 13)
    ) +
    ggplot2::labs(
      ## title=paste(label_begin, paste0(stories, '-story'), paste(systems, collapse=", "), label_end),
      x='Parameter',
      y='Benefit-cost ratio')
return(plot.sen)
}

#' @export
#'
#' @importFrom forcats fct_rev
#' @importFrom dplyr filter select rename mutate
#' @importFrom tidyr pivot_wider
#'
postprocess_eal <- function(output, systems="RCMF", designs="nonstructural", stories=4) {
  return(
    output |>
    dplyr::filter((system %in% systems) &
                  (design %in% designs) &
                  num_stories %in% stories) |>
    dplyr::select(model, system, num_stories, design, label, total_area, project, repair_cost, starts_with('loss')) |>
    dplyr::mutate(project=project/50) |>
    dplyr::select(!loss_ratio) |>
    dplyr::rename(loss_repair_cost=repair_cost) |>
    dplyr::filter(label == 'base') |>
    dplyr::mutate(design=factor(design, levels=rev(designs))) |>
    tidyr::pivot_longer(cols=!c('model', 'system', 'design', 'num_stories', 'label', 'total_area', 'project'),
                        names_to='loss_category', values_to='loss') |>
    dplyr::mutate(loss_category=forcats::fct_relevel(
                                           forcats::fct_rev(loss_category),
                                           'loss_total', after=Inf))
  )
}

#' Plot EALs disaggregated by loss category
#' @export
#'
#' @import ggplot2
#' @importFrom scales label_dollar
#' @importFrom ggthemes scale_fill_colorblind
#'
plot_eal_by_loss <- function(output, systems="RCMF", designs="nonstructural", stories=4) {
  ## TODO: add facet_wrap to plot for multiple systems
  ## TODO: do not hardcode axis tick labels, apply formatting
  labels_x = c("Total", "Business Income", "Displacement", "Rental Income", "Repair Cost", "Value Added")
  plot.eal <- postprocess_eal(output, systems, designs, stories) |>
    dplyr::mutate(
             system=factor(system, levels=systems),
             design=factor(design, levels=rev(designs)),
             num_stories=factor(num_stories)) |>
    ggplot(aes(x=loss_category, y=loss, fill=design, pattern=design)) +
    geom_col(position='dodge', width=0.75) +
    facet_wrap(~system) +
    ## ggplot2::theme_light() +
    ## ggplot2::scale_x_discrete(labels=label_format()) +
    ggplot2::scale_x_discrete(labels=rev(labels_x)) +
    ggplot2::scale_y_continuous(labels = scales::label_dollar(scale_cut=scales::cut_short_scale())) +
    ## TODO: Add geom_text labels for dollar amounts
    ggplot2::coord_flip() +
    ggthemes::theme_few(base_size=10) +
    ggplot2::scale_fill_manual(breaks=rev, values=rev(ggthemes::colorblind_pal()(length(systems)))) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    ) +
    ggplot2::labs(
               fill = "",
               x = "Loss",
               y = "Loss Category")
  return(plot.eal)
}

#' Plot Total EALs across multiple systems and stories
#' @export
#'
#' @import ggplot2
#' @importFrom scales dollar
#' @importFrom ggthemes theme_few
#'
plot_eal <- function(output, systems="RCMF", designs="nonstructural", stories=4, w=FALSE) {
  plot.eal <- postprocess_eal(output, systems, designs, stories) |>
    dplyr::filter(loss_category %in% "loss_total") |>
    ## dplyr::mutate(loss=ifelse(w==TRUE, loss/project, loss)) |>
    dplyr::mutate(
             system=factor(system, levels=systems),
             design=factor(design, levels=designs),
             num_stories=factor(num_stories)) |>
    ggplot(aes(x = num_stories, y = loss, fill = design)) +
    geom_bar(
      stat='identity',
      width = 0.5,
      position = 'dodge') +
    facet_wrap(~system) +
    ggplot2::scale_y_continuous(labels = scales::label_dollar(scale_cut=scales::cut_short_scale())) +
    labs(
      ## title = "EALs for baseline and recovery-based designs",
      x = "Stories",
      y = "EAL") +
    ggthemes::theme_few() +
    ## theme(axis.text.y = element_text(angle = 0,  hjust = 1, size = 15)) +
    ## coord_flip() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 13)
    ) +
    labs(fill = "")
  return(plot.eal)
}


#' @export
#'
#' @importFrom dplyr filter select mutate
#' @import ggplot2
#' @importFrom scales percent
#' @importFrom ggthemes scale_fill_colorblind theme_few
#'
plot_cost_delta <- function(output, systems="RCMF", designs="nonstructural", stories=4) {
  plot.cost.delta <- output |>
    dplyr::filter((system %in% systems) &
                  (design %in% designs) &
                  (num_stories %in% stories)) |>
    dplyr::select(system, design, num_stories, cost_delta) |>
    dplyr::mutate(
             system=factor(system, levels=systems),
             num_stories=factor(num_stories, levels=stories),
             design=factor(design, levels=designs)) |>
    ggplot(aes(x = design, y = cost_delta, fill = num_stories)) +
    geom_bar(
      stat='identity',
      width = 0.5,
      position = 'dodge') +
    facet_wrap(~system) +
    labs(
      ## title = "Cost deltas for recovery-based design interventions",
      x = "System",
      y = "Cost delta") +
    ggthemes::theme_few(base_size=12) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 13)
    ) +
    ggplot2::labs(fill = "Stories")
  return(plot.cost.delta)
}
