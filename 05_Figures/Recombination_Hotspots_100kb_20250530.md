Untitled
================

# Recombination Hot- and Coldspots

------------------------------------------------------------------------

## Prepare environment

``` r
## Clean environment
rm(list = ls())
#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))

# Load libraries
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.2     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.4     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
#library(dplyr)
#library(scales)
#library(ggplot2)
library(patchwork)
library(cowplot)
```

    ## 
    ## Attaching package: 'cowplot'
    ## 
    ## The following object is masked from 'package:patchwork':
    ## 
    ##     align_plots
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     stamp

``` r
#library(grid)
#library(gridExtra)
library(zoo)
```

    ## 
    ## Attaching package: 'zoo'
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     as.Date, as.Date.numeric

``` r
# Set directory
setwd("~/sciebo/RecombinationLandscape_CRIP/02_Revision_2/BMC_Genomics_CripRecLand/New_Results/")
getwd()
```

    ## [1] "/home/alle/sciebo/RecombinationLandscape_CRIP/02_Revision_2/BMC_Genomics_CripRecLand/New_Results"

## Load recombination data files

``` r
#-------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------
ind <- c("MF1", "MF2", "MF3", "MF4", "MG2", "MG3", "MG4", "MG5", "NMF1", "NMF2", "NMF3", "NMF4", 
         "SI1", "SI2", "SI3", "SI4", "SS1", "SS2", "SS3", "SS4")
chr <- c("Chr1", "Chr2", "Chr3", "Chr4")

common_path = "~/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/iSMC/"

# Chr1
files_to_read_1 = list.files(
  path = common_path,        # directory to search within
  pattern = "*Chr1_ismc.rho.100kb.bedgraph", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
data_lst_1 = lapply(files_to_read_1, read.table, header = TRUE)  # read all the matching files

# Chr2
files_to_read_2 = list.files(
  path = common_path,        # directory to search within
  pattern = "*Chr2_ismc.rho.100kb.bedgraph", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
data_lst_2 = lapply(files_to_read_2, read.table, header = TRUE)  # read all the matching files

# Chr3
files_to_read_3 = list.files(
  path = common_path,        # directory to search within
  pattern = "*Chr3_ismc.rho.100kb.bedgraph", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
data_lst_3 = lapply(files_to_read_3, read.table, header = TRUE)  # read all the matching files

# Chr4
files_to_read_4 = list.files(
  path = common_path,        # directory to search within
  pattern = "*Chr4_ismc.rho.100kb.bedgraph", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE,          # return the full path
)
data_lst_4 = lapply(files_to_read_4, read.table, header = TRUE)  # read all the matching files

#-------------------------------------------------------------------
# Rename columns
#-------------------------------------------------------------------
# Rename data list names
names(data_lst_1) <- ind
names(data_lst_2) <- ind
names(data_lst_3) <- ind
names(data_lst_4) <- ind

# Check if renaming worked fine
#Make empty list, otherwise loop won't work
check_rename <- list()

for (i in 1:20){
  check_rename[i] <- identical(names(data_lst_1[i]), colnames(data_lst_1[[i]][4]))
}
check_rename
```

    ## [[1]]
    ## [1] TRUE
    ## 
    ## [[2]]
    ## [1] TRUE
    ## 
    ## [[3]]
    ## [1] TRUE
    ## 
    ## [[4]]
    ## [1] TRUE
    ## 
    ## [[5]]
    ## [1] TRUE
    ## 
    ## [[6]]
    ## [1] TRUE
    ## 
    ## [[7]]
    ## [1] TRUE
    ## 
    ## [[8]]
    ## [1] TRUE
    ## 
    ## [[9]]
    ## [1] TRUE
    ## 
    ## [[10]]
    ## [1] TRUE
    ## 
    ## [[11]]
    ## [1] TRUE
    ## 
    ## [[12]]
    ## [1] TRUE
    ## 
    ## [[13]]
    ## [1] TRUE
    ## 
    ## [[14]]
    ## [1] TRUE
    ## 
    ## [[15]]
    ## [1] TRUE
    ## 
    ## [[16]]
    ## [1] TRUE
    ## 
    ## [[17]]
    ## [1] TRUE
    ## 
    ## [[18]]
    ## [1] TRUE
    ## 
    ## [[19]]
    ## [1] TRUE
    ## 
    ## [[20]]
    ## [1] TRUE

``` r
for (i in 1:20){
  check_rename[i] <- identical(names(data_lst_2[i]), colnames(data_lst_2[[i]][4]))
}
check_rename
```

    ## [[1]]
    ## [1] TRUE
    ## 
    ## [[2]]
    ## [1] TRUE
    ## 
    ## [[3]]
    ## [1] TRUE
    ## 
    ## [[4]]
    ## [1] TRUE
    ## 
    ## [[5]]
    ## [1] TRUE
    ## 
    ## [[6]]
    ## [1] TRUE
    ## 
    ## [[7]]
    ## [1] TRUE
    ## 
    ## [[8]]
    ## [1] TRUE
    ## 
    ## [[9]]
    ## [1] TRUE
    ## 
    ## [[10]]
    ## [1] TRUE
    ## 
    ## [[11]]
    ## [1] TRUE
    ## 
    ## [[12]]
    ## [1] TRUE
    ## 
    ## [[13]]
    ## [1] TRUE
    ## 
    ## [[14]]
    ## [1] TRUE
    ## 
    ## [[15]]
    ## [1] TRUE
    ## 
    ## [[16]]
    ## [1] TRUE
    ## 
    ## [[17]]
    ## [1] TRUE
    ## 
    ## [[18]]
    ## [1] TRUE
    ## 
    ## [[19]]
    ## [1] TRUE
    ## 
    ## [[20]]
    ## [1] TRUE

``` r
for (i in 1:20){
  check_rename[i] <- identical(names(data_lst_3[i]), colnames(data_lst_3[[i]][4]))
}
check_rename
```

    ## [[1]]
    ## [1] TRUE
    ## 
    ## [[2]]
    ## [1] TRUE
    ## 
    ## [[3]]
    ## [1] TRUE
    ## 
    ## [[4]]
    ## [1] TRUE
    ## 
    ## [[5]]
    ## [1] TRUE
    ## 
    ## [[6]]
    ## [1] TRUE
    ## 
    ## [[7]]
    ## [1] TRUE
    ## 
    ## [[8]]
    ## [1] TRUE
    ## 
    ## [[9]]
    ## [1] TRUE
    ## 
    ## [[10]]
    ## [1] TRUE
    ## 
    ## [[11]]
    ## [1] TRUE
    ## 
    ## [[12]]
    ## [1] TRUE
    ## 
    ## [[13]]
    ## [1] TRUE
    ## 
    ## [[14]]
    ## [1] TRUE
    ## 
    ## [[15]]
    ## [1] TRUE
    ## 
    ## [[16]]
    ## [1] TRUE
    ## 
    ## [[17]]
    ## [1] TRUE
    ## 
    ## [[18]]
    ## [1] TRUE
    ## 
    ## [[19]]
    ## [1] TRUE
    ## 
    ## [[20]]
    ## [1] TRUE

``` r
for (i in 1:20){
  check_rename[i] <- identical(names(data_lst_4[i]), colnames(data_lst_4[[i]][4]))
}
check_rename
```

    ## [[1]]
    ## [1] TRUE
    ## 
    ## [[2]]
    ## [1] TRUE
    ## 
    ## [[3]]
    ## [1] TRUE
    ## 
    ## [[4]]
    ## [1] TRUE
    ## 
    ## [[5]]
    ## [1] TRUE
    ## 
    ## [[6]]
    ## [1] TRUE
    ## 
    ## [[7]]
    ## [1] TRUE
    ## 
    ## [[8]]
    ## [1] TRUE
    ## 
    ## [[9]]
    ## [1] TRUE
    ## 
    ## [[10]]
    ## [1] TRUE
    ## 
    ## [[11]]
    ## [1] TRUE
    ## 
    ## [[12]]
    ## [1] TRUE
    ## 
    ## [[13]]
    ## [1] TRUE
    ## 
    ## [[14]]
    ## [1] TRUE
    ## 
    ## [[15]]
    ## [1] TRUE
    ## 
    ## [[16]]
    ## [1] TRUE
    ## 
    ## [[17]]
    ## [1] TRUE
    ## 
    ## [[18]]
    ## [1] TRUE
    ## 
    ## [[19]]
    ## [1] TRUE
    ## 
    ## [[20]]
    ## [1] TRUE

``` r
# Rename columns in every data frame
# https://stackoverflow.com/questions/28648513/changing-column-names-in-a-list-of-data-frames-in-r
new_colnames <- c("chromosome", "start", "end", "rho")

data_lst_1 <- lapply(data_lst_1, rename_with, ~ new_colnames)
data_lst_2 <- lapply(data_lst_2, rename_with, ~ new_colnames)
data_lst_3 <- lapply(data_lst_3, rename_with, ~ new_colnames)
data_lst_4 <- lapply(data_lst_4, rename_with, ~ new_colnames)

# Change names in chromosome column to the correct chromosome
for (i in 1:20){
  data_lst_1[[i]][1] <- "Chr1"
}

for (i in 1:20){
  data_lst_2[[i]][1] <- "Chr2"
}

for (i in 1:20){
  data_lst_3[[i]][1] <- "Chr3"
}

for (i in 1:20){
  data_lst_4[[i]][1] <- "Chr4"
}

# Add column with individual name
for (i in 1:20){
  data_lst_1[[i]]$individual <- names(data_lst_1[i])
}

for (i in 1:20){
  data_lst_2[[i]]$individual <- names(data_lst_2[i])
}

for (i in 1:20){
  data_lst_3[[i]]$individual <- names(data_lst_3[i])
}

for (i in 1:20){
  data_lst_4[[i]]$individual <- names(data_lst_4[i])
}

# Merge dataframes into one big dataframe (df1 = chr1, df2 = chr2, ...)
df1 <- bind_rows(data_lst_1)
df2 <- bind_rows(data_lst_2)
df3 <- bind_rows(data_lst_3)
df4 <- bind_rows(data_lst_4)

df <-  rbind(df1,df2,df3,df4)

# Remove incomplete windows
df <- df %>%
  mutate(wnd_size = end-start)

unique(df$wnd_size) # window size should be only 10,000
```

    ## [1] 100000  41515  72064  47284  16328

``` r
window_size = 100e3

df <- df %>% 
  filter(wnd_size == window_size)

# also for the separate df
df1 <- df1 %>%
  mutate(wnd_size = end-start) %>% 
  filter(wnd_size == window_size)

df2 <- df2 %>%
  mutate(wnd_size = end-start) %>% 
  filter(wnd_size == window_size)

df3 <- df3 %>%
  mutate(wnd_size = end-start) %>% 
  filter(wnd_size == window_size)

df4 <- df4 %>%
  mutate(wnd_size = end-start) %>% 
  filter(wnd_size == window_size)

unique(df4$wnd_size)
```

    ## [1] 100000

``` r
#----------------------------------------------------------
# Calculate mean and confidence interval
## Chr1
head(df1)
```

    ##   chromosome  start    end          rho individual wnd_size
    ## 1       Chr1    508 100508 7.693025e-05        MF1   100000
    ## 2       Chr1 100508 200508 7.461620e-05        MF1   100000
    ## 3       Chr1 200508 300508 7.744472e-05        MF1   100000
    ## 4       Chr1 300508 400508 9.542069e-05        MF1   100000
    ## 5       Chr1 400508 500508 1.318317e-04        MF1   100000
    ## 6       Chr1 500508 600508 2.103774e-04        MF1   100000

``` r
mf_mean1 <- df1 %>% 
  filter(individual == "MF1" | individual == "MF2" | 
           individual == "MF3" | individual == "MF4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(mf_mean1)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho   SE.low  SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>    <dbl>    <dbl>
    ## 1 Chr1          508 100508     4 0.000153 0.000103 0.000203
    ## 2 Chr1       100508 200508     4 0.000545 0.000175 0.000915
    ## 3 Chr1       200508 300508     4 0.00243  0.000992 0.00387 
    ## 4 Chr1       300508 400508     4 0.00145  0.000165 0.00274 
    ## 5 Chr1       400508 500508     4 0.00210  0.00143  0.00277 
    ## 6 Chr1       500508 600508     4 0.00807  0.00467  0.0115

``` r
mg_mean1 <- df1 %>% 
  filter(individual == "MG2" | individual == "MG3" | 
           individual == "MG4" | individual == "MG5") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(mg_mean1)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho   SE.low  SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>    <dbl>    <dbl>
    ## 1 Chr1          508 100508     4 0.000235 0.000157 0.000313
    ## 2 Chr1       100508 200508     4 0.000369 0.000209 0.000528
    ## 3 Chr1       200508 300508     4 0.000695 0.000240 0.00115 
    ## 4 Chr1       300508 400508     4 0.000584 0.000230 0.000938
    ## 5 Chr1       400508 500508     4 0.00209  0.00110  0.00309 
    ## 6 Chr1       500508 600508     4 0.00713  0.00543  0.00883

``` r
nmf_mean1 <- df1 %>% 
  filter(individual == "NMF1" | individual == "NMF2" | 
           individual == "NMF3" | individual == "NMF4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(nmf_mean1)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho   SE.low SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>    <dbl>   <dbl>
    ## 1 Chr1          508 100508     4 0.000846 0.000402 0.00129
    ## 2 Chr1       100508 200508     4 0.00797  0.00465  0.0113 
    ## 3 Chr1       200508 300508     4 0.00801  0.00335  0.0127 
    ## 4 Chr1       300508 400508     4 0.00567  0.00314  0.00820
    ## 5 Chr1       400508 500508     4 0.00780  0.00345  0.0121 
    ## 6 Chr1       500508 600508     4 0.0228   0.0116   0.0340

``` r
si_mean1 <- df1 %>% 
  filter(individual == "SI1" | individual == "SI2" | 
           individual == "SI3" | individual == "SI4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(si_mean1)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho   SE.low  SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>    <dbl>    <dbl>
    ## 1 Chr1          508 100508     4 0.000346 0.000288 0.000403
    ## 2 Chr1       100508 200508     4 0.00425  0.00155  0.00695 
    ## 3 Chr1       200508 300508     4 0.00598  0.00233  0.00963 
    ## 4 Chr1       300508 400508     4 0.00908  0.00399  0.0142  
    ## 5 Chr1       400508 500508     4 0.00902  0.00397  0.0141  
    ## 6 Chr1       500508 600508     4 0.0277   0.0135   0.0420

``` r
ss_mean1 <- df1 %>% 
  filter(individual == "SS1" | individual == "SS2" | 
           individual == "SS3" | individual == "SS4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(ss_mean1)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho   SE.low  SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>    <dbl>    <dbl>
    ## 1 Chr1          508 100508     4 0.000276 0.000223 0.000330
    ## 2 Chr1       100508 200508     4 0.000304 0.000237 0.000372
    ## 3 Chr1       200508 300508     4 0.000579 0.000309 0.000849
    ## 4 Chr1       300508 400508     4 0.00105  0.000503 0.00160 
    ## 5 Chr1       400508 500508     4 0.00170  0.00132  0.00208 
    ## 6 Chr1       500508 600508     4 0.0118   0.00877  0.0147

``` r
## Chr2
head(df2)
```

    ##   chromosome  start    end         rho individual wnd_size
    ## 1       Chr2  34119 134119 0.001633976        MF1   100000
    ## 2       Chr2 134119 234119 0.001678886        MF1   100000
    ## 3       Chr2 234119 334119 0.001851184        MF1   100000
    ## 4       Chr2 334119 434119 0.004022659        MF1   100000
    ## 5       Chr2 434119 534119 0.005846321        MF1   100000
    ## 6       Chr2 534119 634119 0.008752257        MF1   100000

``` r
mf_mean2 <- df2 %>% 
  filter(individual == "MF1" | individual == "MF2" | 
           individual == "MF3" | individual == "MF4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(mf_mean2)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho  SE.low SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr2        34119 134119     4  0.00129 0.00107 0.00152
    ## 2 Chr2       134119 234119     4  0.00159 0.00138 0.00181
    ## 3 Chr2       234119 334119     4  0.00171 0.00155 0.00188
    ## 4 Chr2       334119 434119     4  0.00214 0.00148 0.00281
    ## 5 Chr2       434119 534119     4  0.00290 0.00183 0.00397
    ## 6 Chr2       534119 634119     4  0.00626 0.00542 0.00710

``` r
mg_mean2 <- df2 %>% 
  filter(individual == "MG2" | individual == "MG3" | 
           individual == "MG4" | individual == "MG5") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(mg_mean2)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho   SE.low  SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>    <dbl>    <dbl>
    ## 1 Chr2        34119 134119     4 0.00105  0.000866 0.00123 
    ## 2 Chr2       134119 234119     4 0.00105  0.000830 0.00127 
    ## 3 Chr2       234119 334119     4 0.00107  0.000830 0.00132 
    ## 4 Chr2       334119 434119     4 0.000746 0.000658 0.000834
    ## 5 Chr2       434119 534119     4 0.000715 0.000568 0.000862
    ## 6 Chr2       534119 634119     4 0.00446  0.00263  0.00629

``` r
nmf_mean2 <- df2 %>% 
  filter(individual == "NMF1" | individual == "NMF2" | 
           individual == "NMF3" | individual == "NMF4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(nmf_mean2)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho   SE.low SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>    <dbl>   <dbl>
    ## 1 Chr2        34119 134119     4  0.00113 0.000877 0.00138
    ## 2 Chr2       134119 234119     4  0.00147 0.00113  0.00180
    ## 3 Chr2       234119 334119     4  0.00185 0.00162  0.00208
    ## 4 Chr2       334119 434119     4  0.00251 0.00177  0.00326
    ## 5 Chr2       434119 534119     4  0.00252 0.00184  0.00321
    ## 6 Chr2       534119 634119     4  0.00616 0.00421  0.00811

``` r
si_mean2 <- df2 %>% 
  filter(individual == "SI1" | individual == "SI2" | 
           individual == "SI3" | individual == "SI4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(si_mean2)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho   SE.low SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>    <dbl>   <dbl>
    ## 1 Chr2        34119 134119     4 0.000755 0.000480 0.00103
    ## 2 Chr2       134119 234119     4 0.000987 0.000577 0.00140
    ## 3 Chr2       234119 334119     4 0.00244  0.00158  0.00331
    ## 4 Chr2       334119 434119     4 0.00544  0.00183  0.00906
    ## 5 Chr2       434119 534119     4 0.00218  0.00148  0.00288
    ## 6 Chr2       534119 634119     4 0.00588  0.00351  0.00825

``` r
ss_mean2 <- df2 %>% 
  filter(individual == "SS1" | individual == "SS2" | 
           individual == "SS3" | individual == "SS4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(ss_mean2)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho   SE.low SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>    <dbl>   <dbl>
    ## 1 Chr2        34119 134119     4  0.00174 0.00164  0.00184
    ## 2 Chr2       134119 234119     4  0.00195 0.00179  0.00212
    ## 3 Chr2       234119 334119     4  0.00195 0.00165  0.00226
    ## 4 Chr2       334119 434119     4  0.00143 0.000928 0.00193
    ## 5 Chr2       434119 534119     4  0.00125 0.000704 0.00181
    ## 6 Chr2       534119 634119     4  0.00330 0.00312  0.00349

``` r
## Chr3
head(df3)
```

    ##   chromosome  start    end        rho individual wnd_size
    ## 1       Chr3  15720 115720 0.02849305        MF1   100000
    ## 2       Chr3 115720 215720 0.03097246        MF1   100000
    ## 3       Chr3 215720 315720 0.02531519        MF1   100000
    ## 4       Chr3 315720 415720 0.03671909        MF1   100000
    ## 5       Chr3 415720 515720 0.11595510        MF1   100000
    ## 6       Chr3 515720 615720 0.14063106        MF1   100000

``` r
mf_mean3 <- df3 %>% 
  filter(individual == "MF1" | individual == "MF2" | 
           individual == "MF3" | individual == "MF4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(mf_mean3)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho  SE.low SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr3        15720 115720     4  0.00940 0.00295  0.0159
    ## 2 Chr3       115720 215720     4  0.0114  0.00474  0.0182
    ## 3 Chr3       215720 315720     4  0.0105  0.00524  0.0158
    ## 4 Chr3       315720 415720     4  0.0137  0.00580  0.0217
    ## 5 Chr3       415720 515720     4  0.0339  0.00644  0.0613
    ## 6 Chr3       515720 615720     4  0.0406  0.00720  0.0741

``` r
mg_mean3 <- df3 %>% 
  filter(individual == "MG2" | individual == "MG3" | 
           individual == "MG4" | individual == "MG5") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(mg_mean3)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho  SE.low SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr3        15720 115720     4  0.00640 0.00426 0.00855
    ## 2 Chr3       115720 215720     4  0.00632 0.00420 0.00844
    ## 3 Chr3       215720 315720     4  0.00604 0.00389 0.00818
    ## 4 Chr3       315720 415720     4  0.00555 0.00374 0.00737
    ## 5 Chr3       415720 515720     4  0.00552 0.00364 0.00740
    ## 6 Chr3       515720 615720     4  0.00740 0.00494 0.00985

``` r
nmf_mean3 <- df3 %>% 
  filter(individual == "NMF1" | individual == "NMF2" | 
           individual == "NMF3" | individual == "NMF4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(nmf_mean3)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho  SE.low SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr3        15720 115720     4  0.00809 0.00526  0.0109
    ## 2 Chr3       115720 215720     4  0.0104  0.00740  0.0133
    ## 3 Chr3       215720 315720     4  0.0130  0.0104   0.0156
    ## 4 Chr3       315720 415720     4  0.0125  0.00967  0.0154
    ## 5 Chr3       415720 515720     4  0.0148  0.0116   0.0179
    ## 6 Chr3       515720 615720     4  0.0148  0.0143   0.0153

``` r
si_mean3 <- df3 %>% 
  filter(individual == "SI1" | individual == "SI2" | 
           individual == "SI3" | individual == "SI4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(si_mean3)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho  SE.low SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr3        15720 115720     4  0.00260 0.00142 0.00378
    ## 2 Chr3       115720 215720     4  0.00825 0.00540 0.0111 
    ## 3 Chr3       215720 315720     4  0.00945 0.00620 0.0127 
    ## 4 Chr3       315720 415720     4  0.00953 0.00658 0.0125 
    ## 5 Chr3       415720 515720     4  0.00997 0.00718 0.0128 
    ## 6 Chr3       515720 615720     4  0.00982 0.00656 0.0131

``` r
ss_mean3 <- df3 %>% 
  filter(individual == "SS1" | individual == "SS2" | 
           individual == "SS3" | individual == "SS4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(ss_mean3)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho  SE.low SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr3        15720 115720     4  0.00500 0.00339 0.00661
    ## 2 Chr3       115720 215720     4  0.00663 0.00496 0.00830
    ## 3 Chr3       215720 315720     4  0.00781 0.00613 0.00948
    ## 4 Chr3       315720 415720     4  0.00764 0.00648 0.00879
    ## 5 Chr3       415720 515720     4  0.00805 0.00731 0.00879
    ## 6 Chr3       515720 615720     4  0.0101  0.00916 0.0111

``` r
## Chr4
head(df4)
```

    ##   chromosome  start    end         rho individual wnd_size
    ## 1       Chr4    123 100123 0.005758207        MF1   100000
    ## 2       Chr4 100123 200123 0.008033111        MF1   100000
    ## 3       Chr4 200123 300123 0.008062696        MF1   100000
    ## 4       Chr4 300123 400123 0.008062085        MF1   100000
    ## 5       Chr4 400123 500123 0.008082399        MF1   100000
    ## 6       Chr4 500123 600123 0.008082864        MF1   100000

``` r
mf_mean4 <- df4 %>% 
  filter(individual == "MF1" | individual == "MF2" | 
           individual == "MF3" | individual == "MF4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(mf_mean4)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho  SE.low SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr4          123 100123     4  0.00705 0.00549 0.00861
    ## 2 Chr4       100123 200123     4  0.00859 0.00690 0.0103 
    ## 3 Chr4       200123 300123     4  0.0104  0.00928 0.0116 
    ## 4 Chr4       300123 400123     4  0.0110  0.00960 0.0123 
    ## 5 Chr4       400123 500123     4  0.0110  0.00962 0.0124 
    ## 6 Chr4       500123 600123     4  0.0112  0.00970 0.0126

``` r
mg_mean4 <- df4 %>% 
  filter(individual == "MG2" | individual == "MG3" | 
           individual == "MG4" | individual == "MG5") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(mg_mean4)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho  SE.low SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr4          123 100123     4   0.0101 0.00678  0.0135
    ## 2 Chr4       100123 200123     4   0.0102 0.00683  0.0136
    ## 3 Chr4       200123 300123     4   0.0101 0.00679  0.0135
    ## 4 Chr4       300123 400123     4   0.0100 0.00675  0.0133
    ## 5 Chr4       400123 500123     4   0.0102 0.00684  0.0136
    ## 6 Chr4       500123 600123     4   0.0102 0.00684  0.0136

``` r
nmf_mean4 <- df4 %>% 
  filter(individual == "NMF1" | individual == "NMF2" | 
           individual == "NMF3" | individual == "NMF4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(nmf_mean4)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho  SE.low SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr4          123 100123     4   0.0264 0.00955  0.0433
    ## 2 Chr4       100123 200123     4   0.0523 0.0118   0.0928
    ## 3 Chr4       200123 300123     4   0.0523 0.0120   0.0925
    ## 4 Chr4       300123 400123     4   0.0252 0.0117   0.0386
    ## 5 Chr4       400123 500123     4   0.0455 0.0126   0.0784
    ## 6 Chr4       500123 600123     4   0.0604 0.0153   0.105

``` r
si_mean4 <- df4 %>% 
  filter(individual == "SI1" | individual == "SI2" | 
           individual == "SI3" | individual == "SI4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(si_mean4)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho  SE.low SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr4          123 100123     4  0.00573 0.00238 0.00909
    ## 2 Chr4       100123 200123     4  0.00609 0.00250 0.00967
    ## 3 Chr4       200123 300123     4  0.00657 0.00283 0.0103 
    ## 4 Chr4       300123 400123     4  0.00708 0.00310 0.0111 
    ## 5 Chr4       400123 500123     4  0.00758 0.00333 0.0118 
    ## 6 Chr4       500123 600123     4  0.00964 0.00601 0.0133

``` r
ss_mean4 <- df4 %>% 
  filter(individual == "SS1" | individual == "SS2" | 
           individual == "SS3" | individual == "SS4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
head(ss_mean4)
```

    ## # A tibble: 6 × 7
    ## # Groups:   chromosome, start [6]
    ##   chromosome  start    end     N mean_rho  SE.low SE.high
    ##   <chr>       <int>  <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr4          123 100123     4  0.00457 0.00164 0.00750
    ## 2 Chr4       100123 200123     4  0.00588 0.00281 0.00894
    ## 3 Chr4       200123 300123     4  0.00695 0.00368 0.0102 
    ## 4 Chr4       300123 400123     4  0.00723 0.00389 0.0106 
    ## 5 Chr4       400123 500123     4  0.00751 0.00390 0.0111 
    ## 6 Chr4       500123 600123     4  0.00704 0.00369 0.0104

``` r
#----------------------------------------------------------
```

## Define themes

``` r
pop <-  c("MF", 
          "NMF", 
          "MG", 
          "SI", 
          "SS")

#my_palette <- c("#355BA4",
#                "#8f4e8a",
#                "#287f97",
#                 "#b26712",
#                "#AB3232")

my_palette <- c("#4C83EB",
                "#CD70C6",
                "#3AB6D8",
                "#FF941A",
                "#AB3232")

# Define theme
mytheme <- theme(axis.text.x = element_text(size = 14, color = "black"),
                 axis.text.y = element_text(size = 14, color = "black"),
                 axis.title.y = element_text(size = 14,face = "bold", color = "black"),
                 axis.title.x = element_text(size = 14,face = "bold", color = "black"),
                 title = element_text(size = 12, color = "black"),
                 panel.background = element_rect(fill="white"),
                 panel.grid.major = element_line("lightgrey"),
                 panel.grid.minor = element_line("white"),
                 panel.border = element_rect("black", fill = NA),
                 plot.background = element_rect(fill="white"),
                 legend.background = element_rect(fill="white"),
                 legend.position="bottom")
```

## Define hot- and colgdspots

Mean + 2xSD = hot

Mean - 2x SD = cold

for each chromosome

``` r
# Calculate mean, SD and confidence interval per chromosome
head(df)
```

    ##   chromosome  start    end          rho individual wnd_size
    ## 1       Chr1    508 100508 7.693025e-05        MF1   100000
    ## 2       Chr1 100508 200508 7.461620e-05        MF1   100000
    ## 3       Chr1 200508 300508 7.744472e-05        MF1   100000
    ## 4       Chr1 300508 400508 9.542069e-05        MF1   100000
    ## 5       Chr1 400508 500508 1.318317e-04        MF1   100000
    ## 6       Chr1 500508 600508 2.103774e-04        MF1   100000

``` r
# calculate global mean and sd per chromosome
rho_chr_mean <- df %>%
  group_by(chromosome) %>% 
  summarise(N_global = n(),
            mean_rho_global = mean(rho),
            sd.rho_global = sd(rho),
            SE.low_global = mean_rho_global - (sd(rho)/sqrt(N_global)),
            SE.high_global = mean_rho_global + (sd(rho)/sqrt(N_global)))
  

# now define hotspots
# mean per population per postition
df_mean <- df %>% 
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            sd.rho = sd(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
```

    ## `summarise()` has grouped output by 'chromosome', 'start'. You can override
    ## using the `.groups` argument.

``` r
# join global mean with mean per populations
df_mean_joined <- df_mean %>%
  left_join(rho_chr_mean, by = "chromosome") %>% 
  mutate(spot = case_when(
    mean_rho >= (mean_rho_global + 2*sd.rho_global) ~ "hotspot", 
    mean_rho <= (mean_rho_global - 2*sd.rho_global) ~ "coldspot",
    TRUE ~ "expected_range")
    )
write.csv(df_mean_joined, "chromosome_mean_rho_and_global_mean_with_hotspots_100kb.csv")

unique(df_mean_joined$spot)
```

    ## [1] "expected_range"

``` r
count(df_mean_joined, spot)
```

    ## # A tibble: 1,902 × 4
    ## # Groups:   chromosome, start [1,902]
    ##    chromosome  start spot               n
    ##    <chr>       <int> <chr>          <int>
    ##  1 Chr1          508 expected_range     1
    ##  2 Chr1       100508 expected_range     1
    ##  3 Chr1       200508 expected_range     1
    ##  4 Chr1       300508 expected_range     1
    ##  5 Chr1       400508 expected_range     1
    ##  6 Chr1       500508 expected_range     1
    ##  7 Chr1       600508 expected_range     1
    ##  8 Chr1       700508 expected_range     1
    ##  9 Chr1       800508 expected_range     1
    ## 10 Chr1       900508 expected_range     1
    ## # ℹ 1,892 more rows

``` r
# no hotspot found

# join global mean with full data set
df_joined <- df %>%
  left_join(rho_chr_mean, by = "chromosome") %>% 
  mutate(spot = case_when(
    rho >= (mean_rho_global + 2*sd.rho_global) ~ "hotspot", 
    rho <= (mean_rho_global - 2*sd.rho_global) ~ "coldspot",
    TRUE ~ "expected_range")
    )

unique(df_joined$spot)
```

    ## [1] "expected_range" "hotspot"

``` r
count(df_joined, spot)
```

    ##             spot     n
    ## 1 expected_range 36808
    ## 2        hotspot  1232

``` r
# plot it
ggplot(df_joined, aes(x = rho)) +
  geom_histogram(bins = 50, fill = "grey") +
  geom_vline(aes(xintercept = mean_rho_global +  2*sd.rho_global), color = "red") +
  geom_vline(aes(xintercept = mean_rho_global -  2*sd.rho_global), color = "blue") +  facet_wrap(~chromosome)
```

![](Recombination_Hotspots_100kb_20250530_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

no coldspots found

## Plot heatmap

### Load chromosome info

``` r
#read-in chromosome and scaffold data 
chr.data <- read.table("~/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/input-files_figures/crip4.0_genome_chr_length.csv", sep = ";", header = TRUE)

class(chr.data)
```

    ## [1] "data.frame"

``` r
head(chr.data)
```

    ##        Name Start      End
    ## 1      Chr1     1 61357614
    ## 2      Chr2     1 58906861
    ## 3      Chr3     1 53163979
    ## 4      Chr4     1 17018963
    ## 5 Scaffold1     1   645517
    ## 6 Scaffold3     1   204234

``` r
chr.data <- chr.data %>% 
  dplyr::rename("chr" = "Name",
         "chr_start" = "Start",
         "chr_end" = "End")

str(chr.data)
```

    ## 'data.frame':    14 obs. of  3 variables:
    ##  $ chr      : chr  "Chr1" "Chr2" "Chr3" "Chr4" ...
    ##  $ chr_start: int  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ chr_end  : int  61357614 58906861 53163979 17018963 645517 204234 175650 101174 90000 75203 ...

``` r
chr.data
```

    ##           chr chr_start  chr_end
    ## 1        Chr1         1 61357614
    ## 2        Chr2         1 58906861
    ## 3        Chr3         1 53163979
    ## 4        Chr4         1 17018963
    ## 5   Scaffold1         1   645517
    ## 6   Scaffold3         1   204234
    ## 7   Scaffold4         1   175650
    ## 8   Scaffold5         1   101174
    ## 9   Scaffold6         1    90000
    ## 10  Scaffold7         1    75203
    ## 11  Scaffold8         1    28000
    ## 12  Scaffold9         1    24901
    ## 13 Scaffold10         1    23999
    ## 14 Scaffold11         1    21354

``` r
# define chromosomes of interest
target <- c("Chr1", "Chr2", "Chr3", "Chr4")

# modify chromosome data frame to filter and start at 0 and get end in Mb
chromosomes <- chr.data %>%
  filter(chr %in% target) %>%
  select(chr, chr_end) %>%
  mutate(chr_start = 0) %>%
  mutate(chr_end = chr_end/1e6)

chromosomes
```

    ##    chr  chr_end chr_start
    ## 1 Chr1 61.35761         0
    ## 2 Chr2 58.90686         0
    ## 3 Chr3 53.16398         0
    ## 4 Chr4 17.01896         0

``` r
# add centromere ranges
# Read in centromere ranges
cenrange <- read.csv("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/CentromereRange.csv", sep = ",", header = TRUE)
head(cenrange)
```

    ##   censtart   cenend Chromosome   length
    ## 1 24007501 35892501       Chr1 11885000
    ## 2 27232501 36202501       Chr2  8970000
    ## 3 20902501 27762501       Chr3  6860000
    ## 4  7722501  9572501       Chr4  1850000

``` r
cenrange <- cenrange %>%
  mutate(
    censtart = censtart/1e6,
    cenend = cenend/1e6
  ) %>%
  rename("chr" = "Chromosome")
```

## Load Cla-element data

``` r
cla_colnames = c("chr", "start", "end", "repeatname")

#cla <- read.table("/Users/laupe/Downloads/Cla1.master_AllPop.bedgraph" , header = FALSE, na.string = "NA", sep = "\t", col.names = cla_colnames)
cla <- read.table("/home/alle/recombination-map/MELT/MELT-run2/Cla1/AllPop/2-GroupAnalysis//Cla1.master_AllPop.bedgraph" , header = FALSE, na.string = "NA", sep = "\t", col.names = cla_colnames)

cla <- cla %>%
  mutate(length = end - start) %>%
  filter(chr %in% c("Chr1", "Chr2", "Chr3", "Chr4"))

cla$chr <- as.factor(cla$chr)
```

### Plot it

``` r
df_mean <-  df_mean %>%
  rename(chr = chromosome)

hotspots <- df_joined %>%
  filter(spot == "hotspot") %>%
  rename(chr = chromosome)

hotspots <- hotspots %>%
  group_by(chr, start, end) %>%
  summarise(
    n_ind = n_distinct(individual),
    individuals = paste(unique(individual), collapse = ",")
  ) %>%
  ungroup
```

    ## `summarise()` has grouped output by 'chr', 'start'. You can override using the
    ## `.groups` argument.

``` r
write.csv(df_joined, "all_rho_and_global_mean_with_hotspots_100kb.csv")
write.csv(hotspots, "hotspots_rho_100kb.csv")

mytheme <- theme(axis.text.x = element_text(size = 14, color = "black"),
                 #axis.text.y = element_text(size = 16, color = "black"),
                 #axis.title.y = element_text(size = 18,face = "bold", color = "black"),
                 axis.title.x = element_text(size = 18,face = "bold", color = "black"),
                 title = element_text(size = 17, color = "black"),
                 text = element_text(size=17, color = "black"),
                 panel.background = element_rect(fill="white"),
                 panel.grid.major = element_line("white"),
                 panel.grid.minor = element_line("white"),
                 panel.grid.major.x = element_line("lightgrey", linetype = 3),
                 panel.border = element_rect("black", fill = NA),
                 plot.background = element_rect(fill="white"),
                 legend.background = element_rect(fill="white"),
                 legend.position="bottom", 
                 axis.title.y=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank()) 

p <- ggplot(chromosomes) +
  geom_rect(mapping = aes(xmin = chr_start, xmax = chr_end, ymin = 0, ymax = 1), fill = "white") +
  # centromere ranges
  geom_segment(data = cenrange,
             aes(x = censtart, xend = cenend, y = 1.05, yend = 1.05),
             color = "grey25", size = 1.5) +
  # mark hotspots
  #geom_point(data = cla, aes(x = start / 1e6, y = -0.05), color = "black", size = 1.5) +
  geom_rect(data = hotspots,
          aes(xmin = start / 10^6, xmax = end / 10^6, ymin = -0.25, ymax = -0.05),
          color = "grey50", alpha = 0.002) +
  # add cla-elements
  geom_rect(data = cla,
          aes(xmin = start / 10^6, xmax = end / 10^6, ymin = -0.50, ymax = -0.25),
          color = "darkred", alpha = 0.002) +
  # add mean recombination
  geom_rect(data = df_mean, 
            mapping = aes(xmin = start / 10^6, xmax = end / 10^6, ymin = 0, ymax = 1 ,  color = mean_rho)) +
  mytheme + 
  labs(x = "Position (Mb)") +
  scale_x_continuous(labels = scales::comma, n.breaks = 6, position = "top") +
  scale_y_continuous(limits = c(-0.50, 1.05)) +
  #scale_color_gradient2(low = muted("blue"), mid = "orange", high = muted("red"), midpoint = 3, trans = "sqrt") +
  scale_color_viridis_c(option = "F", direction = 1, name = "Mean \u03c1") +
    
  facet_grid(chr ~ .) +
  theme(axis.text.x = element_text(size = 14, color = "black", angle = 0),
        legend.position="bottom",
        legend.key.width = unit(1, "cm"))  # Adjusting text size and angle 
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
print(p)
```

![](Recombination_Hotspots_100kb_20250530_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
cowplot::save_plot("~/sciebo/RecombinationLandscape_CRIP/02_Revision_2/BMC_Genomics_CripRecLand/New_Results/Heatmap_rho_100kb_centromere_cla.svg", p,nrow = 2, ncol = 6, base_asp = .9, dpi = 600, bg = "white", scale = 0.7)
cowplot::save_plot("~/sciebo/RecombinationLandscape_CRIP/02_Revision_2/BMC_Genomics_CripRecLand/New_Results/Heatmap_rho_100kb_centromere_cla.png", p, nrow = 2, ncol = 6, base_asp = .9, dpi = 600, bg = "white", scale = 0.7)
```
