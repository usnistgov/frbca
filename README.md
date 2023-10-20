
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

## Installation

You can install the development version of frbca like so:

You can install frbca from [GitHub](https://github.com) with:

``` r
devtools::install_github("juanfung/frbca")
#> Downloading GitHub repo juanfung/frbca@HEAD
#> vctrs (0.6.3 -> 0.6.4) [CRAN]
#> withr (2.5.0 -> 2.5.1) [CRAN]
#> fansi (1.0.4 -> 1.0.5) [CRAN]
#> dplyr (1.1.2 -> 1.1.3) [CRAN]
#> Installing 4 packages: vctrs, withr, fansi, dplyr
#> Installing packages into '/private/var/folders/zc/54_9ngk97j71v6mcb67622tc001r_w/T/RtmpZ8bsvw/temp_libpathcc8f7839be96'
#> (as 'lib' is unspecified)
#> 
#> The downloaded binary packages are in
#>  /var/folders/zc/54_9ngk97j71v6mcb67622tc001r_w/T//RtmpfcDAdO/downloaded_packages
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/private/var/folders/zc/54_9ngk97j71v6mcb67622tc001r_w/T/RtmpfcDAdO/remotes10a6ba211608/juanfung-frbca-3469e64/DESCRIPTION’ ... OK
#> * preparing ‘frbca’:
#> * checking DESCRIPTION meta-information ... OK
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘frbca_0.0.0.9000.tar.gz’
#> Warning: invalid uid value replaced by that for user 'nobody'
#> Warning: invalid gid value replaced by that for user 'nobody'
#> Installing package into '/private/var/folders/zc/54_9ngk97j71v6mcb67622tc001r_w/T/RtmpZ8bsvw/temp_libpathcc8f7839be96'
#> (as 'lib' is unspecified)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(frbca)
## basic example code
```
