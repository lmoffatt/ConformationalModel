
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ConformationalModel

<!-- badges: start -->
<!-- badges: end -->

The goal of ConformationalModel is to build kinetic models of ion
channel receptors based on the interactions of conformational domains.

## Installation

You can install the development version of ConformationalModel from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lmoffatt/ConformationalModel")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(ConformationalModel)
```

lets define the conformational changes and its organization in domains

``` r
changes <- conformationalChanges(
  conformational_changes = c("Binding", "Rocking"),
  agonists = c("ATP", "")
)

domains <- conformationaldomains(c(1, 2, 1, 2, 1, 2))
```

now define the interactions using the conformational changes ids and the
domain ids:

``` r



interactions <- conformationalInteractions(
  labels = c("RB", "BR", "RBR"),
  players_ids = list(c(2, 1), c(1, 2), c(2, 1, 2)),
  domains_ids =list(
    list(c(6, 1), c(2, 3), c(4, 5)), 
    list(c(1, 2), c(3, 4), c(5, 6)), 
    list(c(6, 1, 2), c(2, 3, 4), c(4, 5, 6)
         )
    )
)
```

we define the conductance interactions in the same way and the function
that translate the conductance interactions to a defined interaction:

``` r

conductances <- conductanceInteractions(
  labels = c("Rocking_Current_factor"),
  players_ids = list(c(2)),
  domains_ids = list(list(c(2), c(4), c(6)))
)

conductance_info <- conductanceInteractionInfo(
  id_conductance_interaction = 1,
  kind = "equilibrium",
  leakeage_label = "Leakeage_current_ratio"
)
```

finally we define the standard states for some conformational changes in
terms of interactions of those standard states

``` r

standard_states <- standardStates(
  id_conformational_change = c(2, 2, 2, 2),
  id_interaction = c(1, 2, 3, 3) ,
  pos_within_interaction = c(1, 2, 1, 3) ,
  count = c(1, 1, 1, 1)
)
```

with all this information we define a conformational model:

``` r
model = ConformationalModel(changes,
                            domains,
                            interactions,
                            conductances,
                            standard_states,
                            conductance_info)
```

now we define each one of the model parameters and we calculate the
rates for those parameters

``` r
parameters = c(
  Binding_on = 10.0,
  Binding_off = 10000,
  Rocking_on = 10000,
  Rocking_off = 100,
  RB = 10,
  RB_0 = 1,
  RB_1 = 1,
  BR = 10,
  BR_0 = 1,
  BR_1 = 1,
  RBR = 1,
  RBR_0 = 1,
  RBR_1 = 1,
  RBR_2 = 1,
  Rocking_Current_factor = 10,
  Leakeage_current_ratio = 1e-6
)


m = make_model(model, parameters)
m
#> $Q0
#>        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
#>  [1,] 0e+00 0e+00 3e+04 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00
#>  [2,] 1e+06 0e+00 0e+00 1e+04 0e+00 1e+04 0e+00 1e+04 0e+00 0e+00 0e+00 0e+00
#>  [3,] 1e+04 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 2e+04 0e+00 0e+00 0e+00
#>  [4,] 0e+00 1e+04 1e+06 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 1e+04 0e+00 1e+04
#>  [5,] 0e+00 2e+06 0e+00 0e+00 0e+00 0e+00 1e+04 0e+00 0e+00 0e+00 1e+04 0e+00
#>  [6,] 0e+00 1e+04 1e+06 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 1e+04
#>  [7,] 0e+00 0e+00 0e+00 1e+06 1e+04 1e+06 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00
#>  [8,] 0e+00 1e+04 1e+06 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 1e+04 0e+00 0e+00
#>  [9,] 0e+00 0e+00 2e+04 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00
#> [10,] 0e+00 0e+00 0e+00 1e+04 0e+00 0e+00 0e+00 1e+04 1e+06 0e+00 0e+00 0e+00
#> [11,] 0e+00 0e+00 0e+00 1e+06 1e+04 0e+00 0e+00 1e+06 0e+00 0e+00 0e+00 0e+00
#> [12,] 0e+00 0e+00 0e+00 1e+04 0e+00 1e+04 0e+00 0e+00 1e+06 0e+00 0e+00 0e+00
#> [13,] 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 1e+04 0e+00 0e+00 1e+06 1e+04 1e+06
#> [14,] 0e+00 0e+00 0e+00 0e+00 3e+06 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00
#> [15,] 0e+00 0e+00 0e+00 0e+00 1e+04 1e+06 0e+00 1e+06 0e+00 0e+00 0e+00 0e+00
#> [16,] 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 1e+06 0e+00 0e+00 0e+00 1e+06 0e+00
#> [17,] 0e+00 0e+00 0e+00 0e+00 0e+00 1e+04 0e+00 1e+04 1e+06 0e+00 0e+00 0e+00
#> [18,] 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 1e+06 1e+04 0e+00
#> [19,] 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 1e+04 0e+00 0e+00 0e+00 0e+00 1e+06
#> [20,] 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00
#> [21,] 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 3e+04 0e+00 0e+00 0e+00
#> [22,] 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 1e+04 0e+00 1e+04
#> [23,] 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00
#> [24,] 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00
#>       [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
#>  [1,] 0e+00     0 0e+00     0 0e+00 0e+00 0e+00     0 0e+00 0e+00 0e+00     0
#>  [2,] 0e+00     0 0e+00     0 0e+00 0e+00 0e+00     0 0e+00 0e+00 0e+00     0
#>  [3,] 0e+00     0 0e+00     0 0e+00 0e+00 0e+00     0 0e+00 0e+00 0e+00     0
#>  [4,] 0e+00     0 0e+00     0 0e+00 0e+00 0e+00     0 0e+00 0e+00 0e+00     0
#>  [5,] 0e+00     0 1e+04     0 0e+00 0e+00 0e+00     0 0e+00 0e+00 0e+00     0
#>  [6,] 0e+00     0 0e+00     0 1e+04 0e+00 0e+00     0 0e+00 0e+00 0e+00     0
#>  [7,] 1e+04     0 0e+00     0 0e+00 0e+00 1e+04     0 0e+00 0e+00 0e+00     0
#>  [8,] 0e+00     0 0e+00     0 1e+04 0e+00 0e+00     0 0e+00 0e+00 0e+00     0
#>  [9,] 0e+00     0 0e+00     0 0e+00 0e+00 0e+00     0 1e+04 0e+00 0e+00     0
#> [10,] 0e+00     0 0e+00     0 0e+00 0e+00 0e+00     0 0e+00 1e+04 0e+00     0
#> [11,] 1e+04     0 0e+00     0 0e+00 1e+04 0e+00     0 0e+00 0e+00 0e+00     0
#> [12,] 0e+00     0 0e+00     0 0e+00 0e+00 0e+00     0 0e+00 1e+04 0e+00     0
#> [13,] 0e+00     0 0e+00     0 0e+00 0e+00 0e+00     0 0e+00 0e+00 1e+04     0
#> [14,] 0e+00     0 0e+00 30000 0e+00 0e+00 0e+00     0 0e+00 0e+00 0e+00     0
#> [15,] 0e+00     0 0e+00     0 0e+00 1e+04 1e+04     0 0e+00 0e+00 0e+00     0
#> [16,] 0e+00 10000 1e+06     0 0e+00 0e+00 0e+00 20000 0e+00 0e+00 0e+00     0
#> [17,] 0e+00     0 0e+00     0 0e+00 0e+00 0e+00     0 0e+00 1e+04 0e+00     0
#> [18,] 0e+00     0 1e+04     0 1e+06 0e+00 0e+00     0 0e+00 0e+00 1e+04     0
#> [19,] 0e+00     0 1e+04     0 1e+06 0e+00 0e+00     0 0e+00 0e+00 1e+04     0
#> [20,] 1e+06     0 0e+00 20000 0e+00 1e+06 1e+06     0 0e+00 0e+00 0e+00 10000
#> [21,] 0e+00     0 0e+00     0 0e+00 0e+00 0e+00     0 0e+00 0e+00 0e+00     0
#> [22,] 0e+00     0 0e+00     0 1e+04 0e+00 0e+00     0 1e+06 0e+00 0e+00     0
#> [23,] 1e+04     0 0e+00     0 0e+00 1e+04 1e+04     0 0e+00 2e+06 0e+00     0
#> [24,] 0e+00     0 0e+00     0 0e+00 0e+00 0e+00 30000 0e+00 0e+00 3e+06     0
#> 
#> $Qa
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#>  [1,]    0   30    0    0    0    0    0    0    0     0     0     0     0
#>  [2,]    0    0    0    0   20    0    0    0    0     0     0     0     0
#>  [3,]    0    0    0   10    0   10    0   10    0     0     0     0     0
#>  [4,]    0    0    0    0    0    0   10    0    0     0    10     0     0
#>  [5,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>  [6,]    0    0    0    0    0    0   10    0    0     0     0     0     0
#>  [7,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>  [8,]    0    0    0    0    0    0    0    0    0     0    10     0     0
#>  [9,]    0    0    0    0    0    0    0    0    0    10     0    10     0
#> [10,]    0    0    0    0    0    0    0    0    0     0     0     0    10
#> [11,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [12,]    0    0    0    0    0    0    0    0    0     0     0     0    10
#> [13,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [14,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [15,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [16,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [17,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [18,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [19,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [20,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [21,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [22,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [23,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [24,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>       [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
#>  [1,]     0     0     0     0     0     0     0     0     0     0     0
#>  [2,]     0     0     0     0     0     0     0     0     0     0     0
#>  [3,]     0     0     0     0     0     0     0     0     0     0     0
#>  [4,]     0     0     0     0     0     0     0     0     0     0     0
#>  [5,]    10     0     0     0     0     0     0     0     0     0     0
#>  [6,]     0    10     0     0     0     0     0     0     0     0     0
#>  [7,]     0     0    10     0     0     0     0     0     0     0     0
#>  [8,]     0    10     0     0     0     0     0     0     0     0     0
#>  [9,]     0     0     0    10     0     0     0     0     0     0     0
#> [10,]     0     0     0     0    10     0     0     0     0     0     0
#> [11,]     0     0    10     0     0     0     0     0     0     0     0
#> [12,]     0     0     0     0     0    10     0     0     0     0     0
#> [13,]     0     0     0     0     0     0    10     0     0     0     0
#> [14,]     0     0     0     0     0     0     0     0     0     0     0
#> [15,]     0     0    10     0     0     0     0     0     0     0     0
#> [16,]     0     0     0     0     0     0     0     0     0     0     0
#> [17,]     0     0     0     0    10    10     0     0     0     0     0
#> [18,]     0     0     0     0     0     0    10     0     0     0     0
#> [19,]     0     0     0     0     0     0    10     0     0     0     0
#> [20,]     0     0     0     0     0     0     0     0     0     0     0
#> [21,]     0     0     0     0     0     0     0     0    30     0     0
#> [22,]     0     0     0     0     0     0     0     0     0    20     0
#> [23,]     0     0     0     0     0     0     0     0     0     0    10
#> [24,]     0     0     0     0     0     0     0     0     0     0     0
#> 
#> $g
#>              [,1]
#>  [1,] 9.99999e-07
#>  [2,] 9.99999e-07
#>  [3,] 9.99990e-06
#>  [4,] 9.99990e-06
#>  [5,] 9.99999e-07
#>  [6,] 9.99990e-06
#>  [7,] 9.99990e-06
#>  [8,] 9.99990e-06
#>  [9,] 9.99900e-05
#> [10,] 9.99900e-05
#> [11,] 9.99990e-06
#> [12,] 9.99900e-05
#> [13,] 9.99900e-05
#> [14,] 9.99999e-07
#> [15,] 9.99990e-06
#> [16,] 9.99990e-06
#> [17,] 9.99900e-05
#> [18,] 9.99900e-05
#> [19,] 9.99900e-05
#> [20,] 9.99900e-05
#> [21,] 9.99001e-04
#> [22,] 9.99001e-04
#> [23,] 9.99001e-04
#> [24,] 9.99001e-04
```
