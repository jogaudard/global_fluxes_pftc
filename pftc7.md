# Reprocessing PFTC7 fluxes
Joseph Gaudard
2025-08-05

# Importing and reading the files

``` r
my_packages <- c(
  "dataDownloader",
  "tidyverse",
  "licoread",
  "fluxible",
  "here",
  "ggpmisc"
)

lapply(my_packages, library, character.only = TRUE)

sites <- c(1:5)

files <- paste0("LI7500_Site ", sites, ".zip")

lapply(files,
       get_file,
       node = "hk2cy",
       remote_path = "raw_data/x_raw_ecosystem_fluxes/LI7500",
       path = "pftc7_raw_data")

files_zip <- paste0("pftc7_raw_data/", files)

lapply(files_zip,
       unzip,
       exdir = "pftc7_raw_data/LI7500")

# ok so those files are full of things we don't want...
not_use <- c("not used", "old")
files_remove <- files |>
  str_remove(".zip") |>
  rep(2) |>
  paste0("pftc7_raw_data/LI7500/", . = _, "/", not_use) |>
  append(
    c("pftc7_raw_data/LI7500/__MACOSX",
      "pftc7_raw_data/LI7500/LI7500_Site 4/4_2600_west_2_night_resp-2023-12-16T200310.txt") # this file is empty and makes import crash
  )

unlink(here(files_remove), recursive = TRUE)
```

``` r

pftc7_raw_data <- import7500("pftc7_raw_data/LI7500", version = "post2023")

pftc7_data <- pftc7_raw_data |>
  mutate(
  f_fluxid_temp = str_replace_all(f_fluxid, c(" " = "_", "__" = "_")), # name convention...
  replicate = case_when(
    str_detect(f_fluxid_temp, "_redo") ~ 2,
    .default = 1
  ),
  f_fluxid_temp = str_remove_all(f_fluxid_temp, "_redo")
  ) |>
  separate_wider_delim(f_fluxid_temp, "_", names = c(
    "site_id",
    "elevation_m_asl",
    "aspect",
    "plot_id",
    "day_night",
    "flux_type"
    ),
    too_many = "drop"
  ) |>
  mutate(
    flux_type = str_extract(flux_type, "^(\\w+)(?=-)|^(\\w+)(?=\\.txt$)|^\\w"), # getting rid of the rest
    flux_type = str_replace_all(flux_type, "r$", "resp"),
    site_id = as.double(site_id),
    elevation_m_asl = as.double(elevation_m_asl),
    plot_id = as.double(plot_id)
  )
```

# Processing

## CO<sub>2</sub>

Wet air correction:

``` r
pftc7_data <- flux_drygas(pftc7_data, `CO2 (umol/mol)`, `H2O (mmol/mol)`)
```

### Exponential model

Fitting the model presented in Zhao, Hammerle, Zeeman, & Wohlfahrt
(2018) .

``` r
# just to save time

start_cut_pftc7 <- 20
end_cut_pftc7 <- 80

pftc7_fits_exp_co2 <- flux_fitting(pftc7_data,
                                   f_conc = `CO2 (umol/mol)_dry`,
                                   fit_type = "exp_zhao18",
                                   start_cut = start_cut_pftc7,
                                   end_cut = end_cut_pftc7,
                                   cut_direction = "from_start")
```

Using `fluxible::flux_quality` to assess the quality of the dataset.

``` r
pftc7_flags_exp_co2 <- flux_quality(pftc7_fits_exp_co2,
                                    f_conc = `CO2 (umol/mol)_dry`,
                                    force_discard = c(
                                      "4_2600_east_4_day_resp.txt", # messed up
                                      "4_2600_east_3_day_photo.txt", # same
                                      "3_2400_east_1_day_resp.txt", # strange
                                      "2_2200_west_1_night_resp.txt", # missing part of measurement
                                      "3_2400_west_4_day_resp−2023−12−14T142617.txt", # non sense
                                      "3_2400_east_4_day_photo.txt" # peak
                                    ),
                                    force_zero = c(
                                      "4_2600_east_2_day_a.txt", # just noise
                                      "4_2600_east_3_night_a−2023−12−16T204255.txt" # just noise
                                    ),
                                    force_lm = c(
                                      "5_2800_east_3_night_a.txt", # exp fit messed up
                                      "1_2000_west_5_day_a−2023−12−14T094709.txt" # outliers are messing up the exp
                                    ),
                                    force_ok = c(
                                      "1_2000_east_1_night_resp−2023−12−11T220436.txt" # not sure why it was discarded
                                    ))
#> 
#>  Total number of measurements: 258
#> 
#>  ok   216     84 %
#>  zero     28      11 %
#>  discard      7   3 %
#>  force_discard    5   2 %
#>  force_zero   1   0 %
#>  force_lm     1   0 %
#>  start_error      0   0 %
#>  no_data      0   0 %
#>  force_ok     0   0 %
#>  no_slope     0   0 %
```

``` r
pftc7_flags_exp_co2 |>
  flux_plot(f_conc = `CO2 (umol/mol)_dry`,
            print_plot = FALSE,
            output = "longpdf",
            f_ylim_upper = 500,
            f_ylim_lower = 350,
            y_text_position = 430,
            f_plotname = "pftc7_exp_co2")
```

Now let’s calculate the fluxes with `fluxible::flux_calc`.

``` r
pftc7_fluxes_exp_co2 <- flux_calc(pftc7_flags_exp_co2,
                                  slope_col = f_slope_corr,
                                  temp_air_col = `Temperature (C)`,
                                  setup_volume = 2197,
                                  atm_pressure = pressure_atm,
                                  plot_area = 1.44,
                                  conc_unit = "ppm",
                                  flux_unit = "umol/m2/s",
                                  cols_keep = c(
                                    "Date", "Time", "f_quality_flag",
                                    "site_id", "elevation_m_asl", "aspect",
                                    "plot_id", "day_night", "flux_type",
                                    "replicate"
                                  ))
```

### Linear model

``` r
pftc7_fits_lin_co2 <- flux_fitting(pftc7_data,
                                   f_conc = `CO2 (umol/mol)_dry`,
                                   fit_type = "linear",
                                   start_cut = start_cut_pftc7,
                                   end_cut = end_cut_pftc7,
                                   cut_direction = "from_start")
```

Using `fluxible::flux_quality` to assess the quality of the dataset.

``` r
pftc7_flags_lin_co2 <- flux_quality(pftc7_fits_lin_co2,
                                    f_conc = `CO2 (umol/mol)_dry`,
                                    rsquared_threshold = 0.5,
                                    force_discard = c(
                                      "4_2600_east_4_day_resp.txt", # messed up
                                      "2_2200_west_1_night_resp.txt", # missing part of measurement
                                      "3_2400_east_4_day_photo.txt", # peak
                                      "3_2400_west_4_day_resp−2023−12−14T142617.txt" # non sense
                                    ),
                                    force_ok = c(
                                      "1_2000_east_1_night_resp−2023−12−11T220436.txt" # not sure why it was discarded
                                    )
                                    )
#> 
#>  Total number of measurements: 258
#> 
#>  ok   151     59 %
#>  zero     85      33 %
#>  discard      19      7 %
#>  force_discard    3   1 %
#>  start_error      0   0 %
#>  no_data      0   0 %
#>  force_ok     0   0 %
#>  force_zero   0   0 %
#>  force_lm     0   0 %
#>  no_slope     0   0 %
```

``` r
pftc7_flags_lin_co2 |>
  flux_plot(f_conc = `CO2 (umol/mol)_dry`,
            print_plot = FALSE,
            output = "longpdf",
            f_ylim_upper = 500,
            f_ylim_lower = 350,
            y_text_position = 430,
            f_plotname = "pftc7_lin_co2")
```

Now let’s calculate the fluxes with `fluxible::flux_calc`.

``` r
pftc7_fluxes_lin_co2 <- flux_calc(pftc7_flags_lin_co2,
                                  slope_col = f_slope_corr,
                                  temp_air_col = `Temperature (C)`,
                                  setup_volume = 2197,
                                  atm_pressure = pressure_atm,
                                  plot_area = 1.44,
                                  conc_unit = "ppm",
                                  flux_unit = "umol/m2/s",
                                  cols_keep = c(
                                    "Date", "Time", "f_quality_flag",
                                    "site_id", "elevation_m_asl", "aspect",
                                    "plot_id", "day_night", "flux_type",
                                    "replicate"
                                  ))
```

## H<sub>2</sub>O

Wet air correction:

``` r
pftc7_data <- flux_drygas(pftc7_data, `H2O (mmol/mol)`, `H2O (mmol/mol)`)
```

### Exponential model

Fitting the model presented in Zhao et al. (2018) .

``` r
pftc7_fits_exp_h2o <- flux_fitting(pftc7_data,
                                   f_conc = `H2O (mmol/mol)`,
                                   fit_type = "exp_zhao18",
                                   start_cut = start_cut_pftc7,
                                   end_cut = end_cut_pftc7,
                                   cut_direction = "from_start")
```

Using `fluxible::flux_quality` to assess the quality of the dataset.

``` r
pftc7_flags_exp_h2o <- flux_quality(pftc7_fits_exp_h2o,
                                    f_conc = `H2O (mmol/mol)`,
                                    rsquared_threshold = 0.5,
                                    ambient_conc = 20,
                                    error = 10)
#> 
#>  Total number of measurements: 258
#> 
#>  ok   237     92 %
#>  zero     16      6 %
#>  start_error      3   1 %
#>  discard      2   1 %
#>  force_discard    0   0 %
#>  no_data      0   0 %
#>  force_ok     0   0 %
#>  force_zero   0   0 %
#>  force_lm     0   0 %
#>  no_slope     0   0 %
```

``` r
pftc7_flags_exp_h2o |>
  flux_plot(f_conc = `H2O (mmol/mol)`,
            print_plot = FALSE,
            output = "longpdf",
            f_ylim_upper = 40,
            f_ylim_lower = 0,
            y_text_position = 15,
            f_plotname = "pftc7_exp_h2o")
```

Now let’s calculate the fluxes with `fluxible::flux_calc`.

``` r
pftc7_fluxes_exp_h2o <- flux_calc(pftc7_flags_exp_h2o,
                                  slope_col = f_slope_corr,
                                  temp_air_col = `Temperature (C)`,
                                  setup_volume = 2197,
                                  atm_pressure = pressure_atm,
                                  plot_area = 1.44,
                                  conc_unit = "mmol/mol",
                                  flux_unit = "mmol/m2/s",
                                  cols_keep = c(
                                    "Date", "Time", "f_quality_flag",
                                    "site_id", "elevation_m_asl", "aspect",
                                    "plot_id", "day_night", "flux_type",
                                    "replicate"
                                  ))
```

### Linear model

``` r
pftc7_fits_lin_h2o <- flux_fitting(pftc7_data,
                                   f_conc = `H2O (mmol/mol)`,
                                   fit_type = "linear",
                                   start_cut = start_cut_pftc7,
                                   end_cut = end_cut_pftc7,
                                   cut_direction = "from_start")
```

Using `fluxible::flux_quality` to assess the quality of the dataset.

``` r
pftc7_flags_lin_h2o <- flux_quality(pftc7_fits_lin_h2o,
                                    f_conc = `H2O (mmol/mol)`,
                                    rsquared_threshold = 0.5,
                                    ambient_conc = 20,
                                    error = 10,
                                    force_discard = c(
                                      "1_2000_east_1_day_photo−2023−12−14T105344.txt" # clearly not zero
                                    ))
#> 
#>  Total number of measurements: 258
#> 
#>  ok   169     66 %
#>  zero     83      32 %
#>  discard      3   1 %
#>  start_error      3   1 %
#>  force_discard    0   0 %
#>  no_data      0   0 %
#>  force_ok     0   0 %
#>  force_zero   0   0 %
#>  force_lm     0   0 %
#>  no_slope     0   0 %
```

``` r
pftc7_flags_lin_h2o |>
  flux_plot(f_conc = `H2O (mmol/mol)`,
            print_plot = FALSE,
            output = "longpdf",
            f_ylim_upper = 40,
            f_ylim_lower = 0,
            y_text_position = 15,
            f_plotname = "pftc7_lin_h2o")
```

Now let’s calculate the fluxes with `fluxible::flux_calc`.

``` r
pftc7_fluxes_lin_h2o <- flux_calc(pftc7_flags_lin_h2o,
                                  slope_col = f_slope_corr,
                                  temp_air_col = `Temperature (C)`,
                                  setup_volume = 2197,
                                  atm_pressure = pressure_atm,
                                  plot_area = 1.44,
                                  conc_unit = "mmol/mol",
                                  flux_unit = "mmol/m2/s",
                                  cols_keep = c(
                                    "Date", "Time", "f_quality_flag",
                                    "site_id", "elevation_m_asl", "aspect",
                                    "plot_id", "day_night", "flux_type",
                                    "replicate"
                                  ))
```

# Comparison

Adjusting to PFTC7 naming convention

``` r
pftc7_fluxible_fluxes_co2 <- pftc7_fluxes_exp_co2 |>
  bind_rows(pftc7_fluxes_lin_co2) |>
  mutate(
    flux_category = "Carbon",
    flux_type = case_when(
      flux_type == "resp" & day_night == "day" ~ "resp_day",
      flux_type == "resp" & day_night == "night" ~ "resp_night",
      flux_type == "photo" ~ "nee",
      .default = flux_type
    )
  )

pftc7_fluxible_fluxes_h2o <- pftc7_fluxes_exp_h2o |>
  bind_rows(pftc7_fluxes_lin_h2o) |>
  mutate(
    flux_category = "Water",
    flux_type = case_when(
      flux_type == "resp" & day_night == "day" ~ "evap_day",
      flux_type == "resp" & day_night == "night" ~ "evap_night",
      flux_type == "photo" ~ "evapotrans",
      .default = flux_type
    )
  )

# keeping the redo when there is a redo
pftc7_fluxible_fluxes <- bind_rows(pftc7_fluxible_fluxes_co2, pftc7_fluxible_fluxes_h2o) |>
  arrange(desc(replicate)) |>
  distinct(site_id, elevation_m_asl, aspect, plot_id, day_night, flux_type, flux_category, f_model, .keep_all = TRUE)
```

Downloading and adjusting published data

``` r
get_file("hk2cy",
         "x_ecosystem_fluxes",
         "x_PFTC7_clean_ecosystem_fluxes_2023.csv",
         "pftc7_raw_data")

pftc7_published_fluxes <- read_csv("pftc7_raw_data/x_PFTC7_clean_ecosystem_fluxes_2023.csv")

pftc7_published_fluxes <- pftc7_published_fluxes |>
  filter(device == "LI-7500")
```

Gathering the data

``` r
pfct7_comparison <- pftc7_published_fluxes |>
  left_join(pftc7_fluxible_fluxes, by = join_by(site_id, elevation_m_asl, aspect, plot_id, day_night, flux_type, flux_category)) |>
  drop_na(f_model) # missing file in raw data
```

![Comparison of fluxes calculated with fluxible and published PFTC7 data
paper.](pftc7_files/figure-commonmark/figure_h2o-1.png)

And we clean after ourselves :)

``` r
unlink(here("pftc7_raw_data"),
       recursive = TRUE)
```

### References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0" line-spacing="2">

<div id="ref-zhaoCalculationDaytimeCO22018" class="csl-entry">

Zhao, P., Hammerle, A., Zeeman, M., & Wohlfahrt, G. (2018). On the
calculation of daytime CO2 fluxes measured by automated closed
transparent chambers. *Agricultural and Forest Meteorology*, *263*,
267–275. doi:
[10.1016/j.agrformet.2018.08.022](https://doi.org/10.1016/j.agrformet.2018.08.022)

</div>

</div>
