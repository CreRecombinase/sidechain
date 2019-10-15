library(purrr)
library(hdf5r)
library(dplyr)
library(tibble)


Rcpp::sourceCpp("src/mz.cpp")


gen_index <- function(offsets, sizes){
  map2(as.numeric(offsets), as.numeric(sizes), ~seq(from = .x + 1, length.out = .y))
}

get_spectra <- function(h5, scans, cols = get_names(h5[["spectra"]])){
  scan_tibble <- tibble::tibble(ScanNumber = scans)
  spec_table <- read_df_h5(h5[["spectra_table"]], c("ScanNumber", "offset", "N")) %>%
    semi_join(scan_tibble, by = "ScanNumber") %>%
    slice(order(scans))
  stopifnot(all.equal(sort(spec_table$ScanNumber), sort(scans)))
  index <- gen_index(spec_table$offset, spec_table$N)
  names(index) <- spec_table$ScanNumber
  map_dfr(index, ~read_df_h5(h5[["spectra"]], cols, index = .x), .id = "scan")
}




get_xic <- function(h5,
                    mz,
                    ppm = 100,
                    cols = get_names(h5[["trace"]])){
  delta <- ppm * 0.000001
  table <- read_df_h5(h5[["trace_table"]])
  id_range <- table[filter_pred_range(table$MZLow, table$MZHigh, mz, delta = delta),]
  index <- gen_index(id_range$offset, id_range$N)
  cols <- unique(c(cols, "MZ"))
  names(index) <- id_range$ID
  map_dfr(index,function(x){
    tdf <- read_df_h5(h5[["trace"]], cols, index = x)
    tdf[filter_pred(tdf$MZ, mz, delta),]
  })
}


sample_mz <- function(h5, seed = 123, N_spectra = 25, n_elem = 100){
  spectra_grp <- h5[["spectra_table/ScanNumber"]]
  set.seed(seed)
  num_spectra <- spectra_grp$dims
  scans <- sample(1:num_spectra, N_spectra, replace = FALSE)
  get_spectra(h5, scans = scans, cols = c("Intensity", "MZ")) %>%
    arrange(desc(Intensity)) %>%
    pull(MZ)  %>%
    magrittr::extract(seq_len(n_elem))
}

xic_test <- function(h5, mzmin = 400, mzmax = 1000, ppm = 100, N = 100, seed = 123){
    set.seed(seed)
    mzs <- sample_mz(h5, seed = seed, N_spectra = N)
    return(map_df(mzs, ~get_xic(h5, .x, ppm)))
}

get_xic_from_spectrum <- function(h5, mz, ppm=100){
  scans <- read_df_h5(h5[["spectra_table"]], c("ScanNumber", "MSOrder", "offset", "N")) %>%
    dplyr::filter(MSOrder == 1)
  indexes <- gen_index(scans$offset, sizes = scans$N)
  map(indexes, function(x){
      tx <- h5[["spectra/MZ"]][x]
      return(filter_vec(target = tx, query = mz, delta = ppm * 0.000001))
  }) %>% flatten_dbl()
}




xic_spectrum_test <- function(h5,
                              mzmin = 400,
                              mzmax = 1000,
                              ppm = 100,
                              N = 100,
                              seed = 123){

  set.seed(seed)

  mzs <- sample_mz(h5, seed = seed, N_spectra = N)

  len  <- length(get_xic_from_spectrum(h5, mzs, ppm = ppm))
  return(len)
}


benchmark <- function(h5,N=100,seed=123){
  #spectrum_test
  spec_df <- read_df_h5(h5[["spectra_table"]])
  max_scan <- nrow(spec_df)
  set.seed(seed)
  scans <- sample(seq_len(max_scan), N, replace = TRUE)
  walk(scans, ~get_spectra(h5, .x))
  #spec_test_end
  xic_spectrum_test(h5)


  }
