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

common_path = "~/recombination-map/ismc-output/final_bedgraphs/"

# Chr1
files_to_read_1 = list.files(
  path = common_path,        # directory to search within
  pattern = "*Chr1_ismc.rho.1kb.bedgraph", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
data_lst_1 = lapply(files_to_read_1, read.table, header = TRUE)  # read all the matching files

# Chr2
files_to_read_2 = list.files(
  path = common_path,        # directory to search within
  pattern = "*Chr2_ismc.rho.1kb.bedgraph", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
data_lst_2 = lapply(files_to_read_2, read.table, header = TRUE)  # read all the matching files

# Chr3
files_to_read_3 = list.files(
  path = common_path,        # directory to search within
  pattern = "*Chr3_ismc.rho.1kb.bedgraph", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
data_lst_3 = lapply(files_to_read_3, read.table, header = TRUE)  # read all the matching files

# Chr4
files_to_read_4 = list.files(
  path = common_path,        # directory to search within
  pattern = "*Chr4_ismc.rho.1kb.bedgraph", # regex pattern, some explanation below
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

    ## [1] 1000  515   64  284  328

``` r
window_size = 1e3

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

    ## [1] 1000

``` r
#----------------------------------------------------------
# Calculate mean and confidence interval
## Chr1
head(df1)
```

    ##   chromosome start  end          rho individual wnd_size
    ## 1       Chr1   508 1508 9.463359e-05        MF1     1000
    ## 2       Chr1  1508 2508 9.074583e-05        MF1     1000
    ## 3       Chr1  2508 3508 8.812839e-05        MF1     1000
    ## 4       Chr1  3508 4508 8.630482e-05        MF1     1000
    ## 5       Chr1  4508 5508 8.504117e-05        MF1     1000
    ## 6       Chr1  5508 6508 8.408031e-05        MF1     1000

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
    ##   chromosome start   end     N mean_rho   SE.low  SE.high
    ##   <chr>      <int> <int> <int>    <dbl>    <dbl>    <dbl>
    ## 1 Chr1         508  1508     4 0.000197 0.000132 0.000263
    ## 2 Chr1        1508  2508     4 0.000172 0.000120 0.000224
    ## 3 Chr1        2508  3508     4 0.000159 0.000112 0.000205
    ## 4 Chr1        3508  4508     4 0.000151 0.000107 0.000195
    ## 5 Chr1        4508  5508     4 0.000147 0.000104 0.000189
    ## 6 Chr1        5508  6508     4 0.000143 0.000101 0.000186

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
    ##   chromosome start   end     N mean_rho   SE.low  SE.high
    ##   <chr>      <int> <int> <int>    <dbl>    <dbl>    <dbl>
    ## 1 Chr1         508  1508     4 0.000278 0.000200 0.000355
    ## 2 Chr1        1508  2508     4 0.000266 0.000190 0.000341
    ## 3 Chr1        2508  3508     4 0.000257 0.000182 0.000332
    ## 4 Chr1        3508  4508     4 0.000250 0.000175 0.000325
    ## 5 Chr1        4508  5508     4 0.000245 0.000170 0.000321
    ## 6 Chr1        5508  6508     4 0.000242 0.000167 0.000318

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
    ##   chromosome start   end     N mean_rho   SE.low  SE.high
    ##   <chr>      <int> <int> <int>    <dbl>    <dbl>    <dbl>
    ## 1 Chr1         508  1508     4 0.000565 0.000300 0.000831
    ## 2 Chr1        1508  2508     4 0.000489 0.000272 0.000706
    ## 3 Chr1        2508  3508     4 0.000439 0.000252 0.000625
    ## 4 Chr1        3508  4508     4 0.000414 0.000242 0.000586
    ## 5 Chr1        4508  5508     4 0.000402 0.000236 0.000568
    ## 6 Chr1        5508  6508     4 0.000394 0.000232 0.000556

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
    ##   chromosome start   end     N mean_rho   SE.low  SE.high
    ##   <chr>      <int> <int> <int>    <dbl>    <dbl>    <dbl>
    ## 1 Chr1         508  1508     4 0.00188  0.000904 0.00285 
    ## 2 Chr1        1508  2508     4 0.00113  0.000631 0.00164 
    ## 3 Chr1        2508  3508     4 0.000766 0.000480 0.00105 
    ## 4 Chr1        3508  4508     4 0.000598 0.000406 0.000789
    ## 5 Chr1        4508  5508     4 0.000517 0.000370 0.000664
    ## 6 Chr1        5508  6508     4 0.000459 0.000344 0.000574

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
    ##   chromosome start   end     N mean_rho   SE.low  SE.high
    ##   <chr>      <int> <int> <int>    <dbl>    <dbl>    <dbl>
    ## 1 Chr1         508  1508     4 0.000338 0.000261 0.000415
    ## 2 Chr1        1508  2508     4 0.000316 0.000250 0.000382
    ## 3 Chr1        2508  3508     4 0.000304 0.000243 0.000365
    ## 4 Chr1        3508  4508     4 0.000296 0.000238 0.000355
    ## 5 Chr1        4508  5508     4 0.000292 0.000235 0.000348
    ## 6 Chr1        5508  6508     4 0.000288 0.000232 0.000344

``` r
## Chr2
head(df2)
```

    ##   chromosome start   end         rho individual wnd_size
    ## 1       Chr2 34119 35119 0.001679298        MF1     1000
    ## 2       Chr2 35119 36119 0.001668090        MF1     1000
    ## 3       Chr2 36119 37119 0.001658504        MF1     1000
    ## 4       Chr2 37119 38119 0.001651847        MF1     1000
    ## 5       Chr2 38119 39119 0.001647120        MF1     1000
    ## 6       Chr2 39119 40119 0.001643633        MF1     1000

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
    ##   chromosome start   end     N mean_rho  SE.low SE.high
    ##   <chr>      <int> <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr2       34119 35119     4  0.00147 0.00130 0.00164
    ## 2 Chr2       35119 36119     4  0.00143 0.00125 0.00161
    ## 3 Chr2       36119 37119     4  0.00139 0.00120 0.00158
    ## 4 Chr2       37119 38119     4  0.00137 0.00117 0.00157
    ## 5 Chr2       38119 39119     4  0.00135 0.00115 0.00156
    ## 6 Chr2       39119 40119     4  0.00134 0.00114 0.00155

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
    ##   chromosome start   end     N mean_rho   SE.low SE.high
    ##   <chr>      <int> <int> <int>    <dbl>    <dbl>   <dbl>
    ## 1 Chr2       34119 35119     4  0.00116 0.000967 0.00135
    ## 2 Chr2       35119 36119     4  0.00114 0.000949 0.00133
    ## 3 Chr2       36119 37119     4  0.00112 0.000933 0.00132
    ## 4 Chr2       37119 38119     4  0.00111 0.000923 0.00131
    ## 5 Chr2       38119 39119     4  0.00111 0.000915 0.00130
    ## 6 Chr2       39119 40119     4  0.00110 0.000910 0.00129

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
    ##   chromosome start   end     N mean_rho   SE.low SE.high
    ##   <chr>      <int> <int> <int>    <dbl>    <dbl>   <dbl>
    ## 1 Chr2       34119 35119     4  0.00143 0.00126  0.00160
    ## 2 Chr2       35119 36119     4  0.00131 0.00111  0.00152
    ## 3 Chr2       36119 37119     4  0.00125 0.00102  0.00148
    ## 4 Chr2       37119 38119     4  0.00122 0.000979 0.00146
    ## 5 Chr2       38119 39119     4  0.00120 0.000957 0.00144
    ## 6 Chr2       39119 40119     4  0.00119 0.000943 0.00143

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
    ##   chromosome start   end     N mean_rho   SE.low SE.high
    ##   <chr>      <int> <int> <int>    <dbl>    <dbl>   <dbl>
    ## 1 Chr2       34119 35119     4 0.00155  0.000893 0.00221
    ## 2 Chr2       35119 36119     4 0.00118  0.000765 0.00159
    ## 3 Chr2       36119 37119     4 0.00102  0.000675 0.00136
    ## 4 Chr2       37119 38119     4 0.000944 0.000627 0.00126
    ## 5 Chr2       38119 39119     4 0.000903 0.000599 0.00121
    ## 6 Chr2       39119 40119     4 0.000877 0.000580 0.00117

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
    ##   chromosome start   end     N mean_rho  SE.low SE.high
    ##   <chr>      <int> <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr2       34119 35119     4  0.00178 0.00167 0.00189
    ## 2 Chr2       35119 36119     4  0.00177 0.00166 0.00187
    ## 3 Chr2       36119 37119     4  0.00176 0.00165 0.00186
    ## 4 Chr2       37119 38119     4  0.00175 0.00165 0.00185
    ## 5 Chr2       38119 39119     4  0.00175 0.00164 0.00185
    ## 6 Chr2       39119 40119     4  0.00174 0.00164 0.00185

``` r
## Chr3
head(df3)
```

    ##   chromosome start   end         rho individual wnd_size
    ## 1       Chr3 15720 16720 0.008413159        MF1     1000
    ## 2       Chr3 16720 17720 0.008381837        MF1     1000
    ## 3       Chr3 17720 18720 0.008387939        MF1     1000
    ## 4       Chr3 18720 19720 0.007947215        MF1     1000
    ## 5       Chr3 19720 20720 0.007294738        MF1     1000
    ## 6       Chr3 20720 21720 0.006712812        MF1     1000

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
    ##   chromosome start   end     N mean_rho  SE.low SE.high
    ##   <chr>      <int> <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr3       15720 16720     4  0.00492 0.00341 0.00644
    ## 2 Chr3       16720 17720     4  0.00491 0.00340 0.00642
    ## 3 Chr3       17720 18720     4  0.00490 0.00339 0.00642
    ## 4 Chr3       18720 19720     4  0.00479 0.00335 0.00622
    ## 5 Chr3       19720 20720     4  0.00461 0.00329 0.00594
    ## 6 Chr3       20720 21720     4  0.00446 0.00323 0.00568

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
    ##   chromosome start   end     N mean_rho  SE.low SE.high
    ##   <chr>      <int> <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr3       15720 16720     4  0.00618 0.00411 0.00825
    ## 2 Chr3       16720 17720     4  0.00617 0.00411 0.00824
    ## 3 Chr3       17720 18720     4  0.00617 0.00410 0.00824
    ## 4 Chr3       18720 19720     4  0.00617 0.00410 0.00824
    ## 5 Chr3       19720 20720     4  0.00617 0.00410 0.00824
    ## 6 Chr3       20720 21720     4  0.00617 0.00410 0.00824

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
    ##   chromosome start   end     N mean_rho  SE.low SE.high
    ##   <chr>      <int> <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr3       15720 16720     4  0.00552 0.00280 0.00824
    ## 2 Chr3       16720 17720     4  0.00544 0.00269 0.00819
    ## 3 Chr3       17720 18720     4  0.00540 0.00263 0.00816
    ## 4 Chr3       18720 19720     4  0.00537 0.00260 0.00814
    ## 5 Chr3       19720 20720     4  0.00531 0.00253 0.00810
    ## 6 Chr3       20720 21720     4  0.00521 0.00241 0.00800

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
    ##   chromosome start   end     N mean_rho  SE.low SE.high
    ##   <chr>      <int> <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr3       15720 16720     4  0.00219 0.00125 0.00313
    ## 2 Chr3       16720 17720     4  0.00219 0.00124 0.00313
    ## 3 Chr3       17720 18720     4  0.00219 0.00124 0.00314
    ## 4 Chr3       18720 19720     4  0.00219 0.00123 0.00314
    ## 5 Chr3       19720 20720     4  0.00219 0.00123 0.00315
    ## 6 Chr3       20720 21720     4  0.00218 0.00122 0.00315

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
    ##   chromosome start   end     N mean_rho  SE.low SE.high
    ##   <chr>      <int> <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr3       15720 16720     4  0.00378 0.00285 0.00471
    ## 2 Chr3       16720 17720     4  0.00377 0.00284 0.00471
    ## 3 Chr3       17720 18720     4  0.00377 0.00283 0.00471
    ## 4 Chr3       18720 19720     4  0.00378 0.00283 0.00472
    ## 5 Chr3       19720 20720     4  0.00379 0.00283 0.00474
    ## 6 Chr3       20720 21720     4  0.00380 0.00283 0.00477

``` r
## Chr4
head(df4)
```

    ##   chromosome start  end         rho individual wnd_size
    ## 1       Chr4   123 1123 0.005037760        MF1     1000
    ## 2       Chr4  1123 2123 0.005040148        MF1     1000
    ## 3       Chr4  2123 3123 0.005043371        MF1     1000
    ## 4       Chr4  3123 4123 0.005042901        MF1     1000
    ## 5       Chr4  4123 5123 0.005043123        MF1     1000
    ## 6       Chr4  5123 6123 0.005044316        MF1     1000

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
    ##   chromosome start   end     N mean_rho  SE.low SE.high
    ##   <chr>      <int> <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr4         123  1123     4  0.00540 0.00433 0.00647
    ## 2 Chr4        1123  2123     4  0.00540 0.00433 0.00647
    ## 3 Chr4        2123  3123     4  0.00540 0.00433 0.00648
    ## 4 Chr4        3123  4123     4  0.00540 0.00433 0.00648
    ## 5 Chr4        4123  5123     4  0.00541 0.00433 0.00649
    ## 6 Chr4        5123  6123     4  0.00542 0.00434 0.00649

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
    ##   chromosome start   end     N mean_rho  SE.low SE.high
    ##   <chr>      <int> <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr4         123  1123     4  0.00999 0.00672  0.0133
    ## 2 Chr4        1123  2123     4  0.00999 0.00672  0.0133
    ## 3 Chr4        2123  3123     4  0.0100  0.00672  0.0133
    ## 4 Chr4        3123  4123     4  0.0100  0.00672  0.0133
    ## 5 Chr4        4123  5123     4  0.0100  0.00672  0.0133
    ## 6 Chr4        5123  6123     4  0.0100  0.00672  0.0133

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
    ##   chromosome start   end     N mean_rho  SE.low SE.high
    ##   <chr>      <int> <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr4         123  1123     4   0.0183 0.00802  0.0286
    ## 2 Chr4        1123  2123     4   0.0181 0.00803  0.0281
    ## 3 Chr4        2123  3123     4   0.0178 0.00804  0.0276
    ## 4 Chr4        3123  4123     4   0.0173 0.00802  0.0266
    ## 5 Chr4        4123  5123     4   0.0171 0.00802  0.0263
    ## 6 Chr4        5123  6123     4   0.0171 0.00804  0.0262

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
    ##   chromosome start   end     N mean_rho  SE.low SE.high
    ##   <chr>      <int> <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr4         123  1123     4  0.00522 0.00226 0.00818
    ## 2 Chr4        1123  2123     4  0.00522 0.00226 0.00818
    ## 3 Chr4        2123  3123     4  0.00522 0.00226 0.00818
    ## 4 Chr4        3123  4123     4  0.00522 0.00226 0.00819
    ## 5 Chr4        4123  5123     4  0.00522 0.00226 0.00819
    ## 6 Chr4        5123  6123     4  0.00523 0.00226 0.00821

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
    ##   chromosome start   end     N mean_rho  SE.low SE.high
    ##   <chr>      <int> <int> <int>    <dbl>   <dbl>   <dbl>
    ## 1 Chr4         123  1123     4  0.00464 0.00185 0.00743
    ## 2 Chr4        1123  2123     4  0.00461 0.00181 0.00741
    ## 3 Chr4        2123  3123     4  0.00458 0.00177 0.00739
    ## 4 Chr4        3123  4123     4  0.00455 0.00173 0.00737
    ## 5 Chr4        4123  5123     4  0.00453 0.00171 0.00735
    ## 6 Chr4        5123  6123     4  0.00452 0.00169 0.00735

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

    ##   chromosome start  end          rho individual wnd_size
    ## 1       Chr1   508 1508 9.463359e-05        MF1     1000
    ## 2       Chr1  1508 2508 9.074583e-05        MF1     1000
    ## 3       Chr1  2508 3508 8.812839e-05        MF1     1000
    ## 4       Chr1  3508 4508 8.630482e-05        MF1     1000
    ## 5       Chr1  4508 5508 8.504117e-05        MF1     1000
    ## 6       Chr1  5508 6508 8.408031e-05        MF1     1000

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
write.csv(df_mean_joined, "chromosome_mean_rho_and_global_mean_with_hotspots_1kb.csv")

unique(df_mean_joined$spot)
```

    ## [1] "expected_range"

``` r
count(df_mean_joined, spot)
```

    ## # A tibble: 190,376 × 4
    ## # Groups:   chromosome, start [190,376]
    ##    chromosome start spot               n
    ##    <chr>      <int> <chr>          <int>
    ##  1 Chr1         508 expected_range     1
    ##  2 Chr1        1508 expected_range     1
    ##  3 Chr1        2508 expected_range     1
    ##  4 Chr1        3508 expected_range     1
    ##  5 Chr1        4508 expected_range     1
    ##  6 Chr1        5508 expected_range     1
    ##  7 Chr1        6508 expected_range     1
    ##  8 Chr1        7508 expected_range     1
    ##  9 Chr1        8508 expected_range     1
    ## 10 Chr1        9508 expected_range     1
    ## # ℹ 190,366 more rows

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

    ##             spot       n
    ## 1 expected_range 3688810
    ## 2        hotspot  118710

``` r
# plot it
ggplot(df_joined, aes(x = rho)) +
  geom_histogram(bins = 50, fill = "grey") +
  geom_vline(aes(xintercept = mean_rho_global +  2*sd.rho_global), color = "red") +
  geom_vline(aes(xintercept = mean_rho_global -  2*sd.rho_global), color = "blue") +  facet_wrap(~chromosome)
```

![](Recombination_Hotspots_1kb_20250530_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

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
write.csv(df_joined, "all_rho_and_global_mean_with_hotspots_1kb.csv")
write.csv(hotspots, "hotspots_rho_1kb.csv")

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

![](Recombination_Hotspots_1kb_20250530_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
cowplot::save_plot("~/sciebo/RecombinationLandscape_CRIP/02_Revision_2/BMC_Genomics_CripRecLand/New_Results/Heatmap_rho_1kb_centromere_cla.svg", p,nrow = 2, ncol = 6, base_asp = .9, dpi = 600, bg = "white", scale = 0.7)
cowplot::save_plot("~/sciebo/RecombinationLandscape_CRIP/02_Revision_2/BMC_Genomics_CripRecLand/New_Results/Heatmap_rho_1kb_centromere_cla.png", p, nrow = 2, ncol = 6, base_asp = .9, dpi = 600, bg = "white", scale = 0.7)
```
