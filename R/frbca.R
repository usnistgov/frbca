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


#' @export
preprocess_model <- function(eal, cost, p) {
    ## Purpose:
    ## join eal and cost tables
    ## compute total floor area (floor area * num_stories)
    ## separate models by height
    ##
    join_cols = c('system', 'model', 'num_stories', 'intervention', 'design_s', 'design_ns')
    models <- list()
    dat <- eal |>
        dplyr::left_join(cost, by=join_cols) |>
        dplyr::mutate(total_area=p$floor_area*num_stories)
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
        dplyr::mutate(delta_loss=loss_total[intervention == 0] - loss_total) |>
        dplyr::rowwise() |>
        dplyr::mutate(benefit=f_npv(t=seq(from=0, to=p$T), cf=rep(delta_loss, length=p$T+1), i=p$delta)) |>
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
                      irr=f_irr(t=seq(from=0, to=p$T), cf=c(delta_loss-cost_diff, rep(delta_loss, length=p$T))),
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
    o_i = models[[i]]
    stories = names(o_i)
    for (j in stories) {
      o_i[[j]] <- bca(models[[i]][[j]], params)
    }
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
postprocess_bcr <- function(output, model_list=c('RCMF-4-baseline-nsfr'), out_base=FALSE) {
  ## function to postprocess output for plotting sensitivity
  ## drop baseline-baseline, if it exists
  model_list = model_list[!grepl('baseline-baseline', model_list)]
  ## create data frame for table/plot
  plot_df <- output |>
    dplyr::filter(model %in% model_list) |>
    dplyr::select(model, bcr, label, parameter)
  base <- plot_df |>
    dplyr::filter(label == 'base') |>
    dplyr::select(!c(label, parameter))
  sen <- plot_df |>
    dplyr::filter(label != 'base') |>
    tidyr::pivot_wider(names_from=label, values_from=bcr) |>
    dplyr::left_join(base, by='model') |>
    dplyr::rename(bcr_low=low, bcr_high=high) |>
    dplyr::mutate(
             model=factor(model, levels=model_list),
             parameter=forcats::fct_rev(parameter))
  if (out_base) {
    return(base)
  } else {
    return(sen)
  }
}

###
## Purpose:
## Post-process data and generate plot for sensitivity analysis
###
#' Plot FR-BCA Outputs
#'
#' @description
#' Sensitivity analysis plot for FR-BCA outputs, for fixed story height and structural system
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
plot_bcr <- function(output, model_list) {
  ## generate plot
  system <- output |> dplyr::distinct(system) |> pull()
  n_floors <- output |> dplyr::distinct(num_stories) |> pull()
  label_begin <- 'Sensitivity Analysis: Benefit-cost ratios for'
  label_end <- 'archetypes, relative to baseline ASCE 7-16 design.'
  plot.sen <- postprocess_bcr(output, model_list) |>
    ggplot2::ggplot() +
    ggplot2::geom_segment(aes(x=parameter, xend=parameter, y=bcr_low, yend=bcr_high),
                 linewidth = 5, colour = "red", alpha = 0.6) +
    ggplot2::geom_segment(aes(x=parameter, xend=parameter, y=bcr-0.001, yend=bcr+0.001),
                 linewidth = 5, colour = "black") +
    ggplot2::geom_hline(yintercept=1, colour='red') +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~model, ncol=1) +
    ## geom_hline(data=rcmf, aes(yintercept=bcr)) +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position='bottom') +
    ggplot2::labs(
      title=paste(label_begin, paste0(n_floors, '-story'), system, label_end),
      x='Parameter',
      y='Benefit-cost ratio')
return(plot.sen)
}

#' @export
#'
#' @importFrom forcats fct_rev
#' @importFrom tidyr pivot_wider
#'
postprocess_eal <- function(output, model_list=c('RCMF-4-baseline-baseline')) {
  return(
    output |>
    dplyr::filter(model %in% model_list) |>
    dplyr::select(model, label, repair_cost, starts_with('loss')) |>
    dplyr::select(!loss_ratio) |>
    dplyr::rename(loss_repair_cost=repair_cost) |>
    dplyr::filter(label == 'base') |>
    dplyr::mutate(model=factor(model, levels=rev(model_list))) |>
    tidyr::pivot_longer(cols=!c('model', 'label'), names_to='loss_category', values_to='loss') |>
    dplyr::mutate(loss_category=forcats::fct_relevel(
                                           forcats::fct_rev(loss_category),
                                           'loss_total', after=Inf))
  )
}

#' @export
#'
#' @importFrom scales label_dollar
#' @importFrom ggthemes scale_fill_colorblind
#'
plot_eal <- function(output, model_list) {
  ## PLACEHOLDER FOR PLOTTING EALs
  plot.eal <- postprocess_eal(output, model_list) |>
    ggplot(aes(x=loss_category, y=loss, fill=model, pattern=model)) +
    geom_col(position='dodge', width=0.5) +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position='bottom') +
    ggplot2::scale_y_continuous(labels = scales::label_dollar()) +
    ## TODO: Add geom_text labels for dollar amounts
    ggplot2::coord_flip() +
    ggthemes::scale_fill_colorblind()
  return(plot.eal)
}
