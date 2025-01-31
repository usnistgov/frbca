
<!-- README.md is generated from README.Rmd. Please edit that file -->

# frbca

<!-- badges: start -->

<!-- badges: end -->

The goal of frbca is to conduct benefit-cost analysis (BCA) for
buildings archetypes designed for functional recovery. The inputs
include:

- Expected Annualized Losses (EALs) for (i) a baseline (code-conforming)
  design and (ii) recovery-based design interventions
- Structural and nonstructural costs for both (i) and (ii)
- A list of analysis parameters, including discount rate, time horizon,
  and economic analysis parameters such as business income

For more details, see

Fung, JF, Cook, DT, Elsibaie, S, Zhang, Y, Sattar, S, Morris, P,
Johnson, KJ, and Burton, HV (2025). \`\`Economic Evaluation at the
Design Phase for Functional Recovery: Integrated Design, Assessment, and
Economic Evaluation for New Buildings,’’ In Review.

## Installation

You can install the development version of frbca like so:

You can install frbca from [GitHub](https://github.com) with:

``` r
devtools::install_github("juanfung/frbca")
```

## Expected inputs

The package includes the data files used for the analysis (Fung et
al. 2025). To run the analysis, provide two inputs in tabular (eg,
`data.frame`) form:

### Expected annualized losses (EALs) from a performance assessment

Each row is an archetype model. Expected (minimum) columns:

    | variable                 | type  | description                        |
    |--------------------------|-------|------------------------------------|
    | system                   | <chr> | structural system type             |
    | num_stories              | <dbl> | number of stories                  |
    | design_s                 | <chr> | structural design intervention     |
    | design_ns                | <chr> | nonstructural design intervention  |
    | model                    | <chr> | unique model name or id            |
    | intervention             | {0,1} | baseline design = 0; otherwise = 1 |
    | loss_ratio               | <dbl> | loss ratio                         |
    | repair_costs             | <dbl> | repair costs, in dollars           |
    | re_occupancy_time        | <dbl> | re-occupancy time, in days         |
    | functional_recovery_time | <dbl> | functional recovery time, in days  |

### Construction costs

Each row is an archetype model. Expected (minimum) columns:

    | variable        | type  | description                                |
    |-----------------|-------|--------------------------------------------|
    | system          | <chr> | structural system type                     |
    | num_stories     | <dbl> | number of stories                          |
    | design_s        | <chr> | structural design intervention             |
    | design_ns       | <chr> | nonstructural design intervention          |
    | model           | <chr> | unique model name or id                    |
    | intervention    | {0,1} | baseline design = 0; otherwise = 1         |
    | c_s             | <dbl> | structural construction costs (dollars)    |
    | c_ns            | <dbl> | nonstructural construction costs (dollars) |
    | construction    | <dbl> | total construction costs (dollars)         |
    | project         | <dbl> | total project costs (dollars)              |

In addition, several parameters are required in list form. The following
is an example of the base parameters used in our analysis
`input_param$parameters$base`:

    | parameter    | type    | description                                        |
    |--------------|---------|----------------------------------------------------|
    | floor_area   | <dbl>   | total square footage per story                     |
    | delta        | <delta> | discount rate                                      |
    | T            | <dbl>   | planning horizon, in years                         |
    | loss         | <dbl>   | economic losses (list)                             |
    | tenant       | <dbl>   | number of tenants per square foot                  |
    | recapture    | <dbl>   | recapture rate for rental income                   |

In addition, parameter values for sensitivity analysis should be
provided under the sub-list `input_param$parameters$sensitivity`. See
documentation and included data for more details.

## Example

This is a basic example which shows you how to conduct benefit-cost
analysis for 4-story RCMF:

``` r
## NOT RUN
library(frbca)

## load package data
input_cost = frbca::input_cost ## or: data("input_cost", package = "frbca")
input_eal = frbca::input_eal
input_param = fbca::input_param

## run analysis (returns results as list)
output <- frbca::frbca(input_eal, input_cost, input_param)

## plot baseline BCR for 4-story RCMF
output |>
  dplyr::bind_rows(output) |>
  frbca::plot_bcr(n_floors=4, system='RCMF')

## plot BCR with sensitivity analysis for 4-story RCMF
output |>
  dplyr::bind_rows(output) |>
  frbca::plot_bcr_sensitivity(n_floors=4, system='RCMF')
```
