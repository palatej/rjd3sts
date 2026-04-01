# Creates a white noise

Creates a white noise

## Usage

``` r
.noise(var = 1)
```

## Arguments

- var:

  Variance of the noise

## Value

A raw java state block

## Examples

``` r
sb<-.noise(.01)
.ssf_T(sb, 0)
#>      [,1]
#> [1,]    0
```
