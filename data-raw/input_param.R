## code to prepare `input_param` dataset goes here

## Set analysis parameters
## NEED TO DOCUMENT
## TODO: json file or other more programmatic way to do this (eg, some defaults)?

input_model_name = 'all'
days = 365.25
cpi_bi = 307.1/257 # analysis year = 2023 / base year = 2019
cpi_ri = 307.1/240 # analysis year = 2023 / base year = 2016
bi_low = 0.75 ## 0.76 * cpi_bi
bi_high = 5.42 ## 5.53 * cpi_bi
rent = 22.06 / days * cpi_ri
rent_low = 16.15 / days * cpi_ri
rent_high = 31.73 / days * cpi_ri
rho = 0.87
va = 7.97 ## 7.25 * cpi_bi
delta_va = 0.035
delta_va_high = 0.12
tenant_per_area = 0.0197113

## Parameters
## floor_area = total building floor area, per story
## total_floors = total stories (NB: currently obtained from building name...not robust!)
## delta = discount rate
## T = time horizon
## bi = business income
## ri = rental income
## displacement = losses due to population displacement
## tenant = ... Need to look at notes
## recapture = business income recapture rate
## sc = supply chain losses multiplier
## TODO: read/write parameters from/to csv or json?
## Need to collect:
## {total_floors, bi_low, bi_high, ri, displacement, tenant}
## Provided defaults:
## {delta, T, recapture, sc}
bca_inputs <- list(
    model = input_model_name,
    parameters = list(
        base=list(
            floor_area=120*120,
            delta=0.03,
            T=50,
          loss=list(
            loss_business_income=(bi_low+bi_high)/2,
            loss_rental_income=rent,
            loss_displacement=112,
            ##loss_displacement=list(
            ##  recurring=112,
            ##  fixed=112),
            ## tenant=0.021,
            ## recapture=rho,
            ## loss_supply_chain=4,
            loss_value_added=va * (delta_va)
          ),
            ## bi=(bi_low+bi_high)/2,
            ## ri=rent,
            ## displacement=112,
            ## sc=4,
            tenant=tenant_per_area,
            recapture=rho
            ),
        sensitivity=list(
          loss=list(
            loss_business_income=list(
              low=bi_low,
              high=bi_high),
            ## loss_supply_chain=list(
            ##   low=2,
            ##   high=10),
            ## loss_displacement=list(
            ##   low=list(recurring=53, fixed=53),
            ##   high=list(recurring=275, fixed=275)),
            loss_rental_income=list(
              low=rent_low,
              high=rent_high),
            loss_value_added=list(
              low=va * (delta_va),
              high= va * (delta_va_high))
            ),
          ## displacement=list(
          ##   low=53,
          ##   high=275),
          ## sc=list(
          ##   low=2,
          ##   high=10),
          delta=list(
                low=0.07,
                high=0.02),
            T=list(
                low=30,
                high=75)
            ## bi=list(
            ##     low=bi_low,
            ##     high=bi_high)
        )
    )
)


usethis::use_data(input_param, overwrite = TRUE)
