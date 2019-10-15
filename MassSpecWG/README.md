The Data
========

We’ll focus on the file, `blacktea_hcd_7500_top3_1.h5`
(`md5sum: 67307cc6d651019f774a79c35d3b28f6`)

Packages
========

There 2 main packages for working with HDF5 data. `hdf5r` and `rhdf5`.
`hdf5r` is in my opinion a little more ergonomic for R users, but it
does require the user to have a pre-existing installation of HDF5 on
their machine if they are a Mac or Linux user (Windows users don’t have
to worry as a windows-compatible binary of HDF5 comes with the package).

I’m a `tidyverse` devotee so I’ll also be loading the `tidyverse`
meta-package as well. In particular I’ll be using the functional
programming library `purrr` for dataframe/list manipulation. While the
`purrr` library is convenient, there is a performance penalty inherent
to my using an immutable data structure/functional programming style, a
penalty that is inherent neither to `HDF5` nor to the `hdf5r` package.

I also wrote some helper functions for getting data out of the
`pytables`/`pandas` HDF5 format and a few more for going from the 1 data
cube per spectrum format to a 1 datacube format. You can find all of
that in `R/functions.R`.

Generating the new H5 file
--------------------------

You can run the R code below to generate the new and improved file
assuming you’ve put the `blacktea_hcd_7500_top3_1.h5` h5 file in the
`data` folder.

``` r
library(hdf5r)
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.2.1           ✔ purrr   0.3.2.9000 
    ## ✔ tibble  2.1.3           ✔ dplyr   0.8.3      
    ## ✔ tidyr   0.8.99.9000     ✔ stringr 1.4.0      
    ## ✔ readr   1.3.1           ✔ forcats 0.4.0

    ## ── Conflicts ────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter()     masks stats::filter()
    ## ✖ purrr::flatten_df() masks hdf5r::flatten_df()
    ## ✖ dplyr::lag()        masks stats::lag()

``` r
source("R/functions.R")
#source("R/mz_fun.R")

file_h5 <- H5File$new("data/blacktea_hcd_7500_top3_1.h5",mode="r")

spectra_grp <- pytable_ls(file_h5[["spectra"]],id_pattern="S([0-9]+)")
trace_grp <- pytable_ls(file_h5[["trace"]],id_pattern="T([0-9]+)")  %>% select(-N)
spec_table_df <- read_pytables(grp = file_h5[["spectra_table"]]) %>%
  left_join(spectra_grp,by=c("ScanNumber"="group"))
trace_table_df <- read_pytables(grp = file_h5[["trace_table"]]) %>%
  left_join(trace_grp,by=c("ID"="group"))  %>%
  mutate(offset=lag(cumsum(N),default=0))

t_1 <- trace_grp$group[1]
proto_trace <- read_pytables(file_h5[[paste0("trace/T",t_1)]])
proto_spec <-  read_pytables(file_h5[[paste0("spectra/S1")]])
```

``` r
new_file <- "new_data/blacktea_hcd_7500_top3_1.h5"

file.remove(new_file[file.exists(new_file)]) #Remove the file if it already exists
```

    ## [1] TRUE

``` r
alt_file1 <- make_file("new_data/blacktea_hcd_7500_top3_1.h5",mode="w")

write_df_h5(spec_table_df,make_group(alt_file1,"spectra_table"))
write_df_h5(trace_table_df,make_group(alt_file1,"trace_table"))

ctd1_trace <- create_transposed_dataset(alt_file1,"trace",proto_trace,trace_table_df)
ctd1_spec <- create_transposed_dataset(alt_file1,"spectra",proto_spec,spec_table_df)

sub_trace_df <- filter(trace_table_df,N>0) %>% arrange(ID)

trace_t1 <- system.time(
  copy_pytables_h5(
    source_grp = file_h5[["trace"]],
    source_names = paste0("T",sub_trace_df$ID),
    dest_offsets = sub_trace_df$offset,
    dsl= ctd1_trace))

spec_t1 <- system.time(
  copy_pytables_h5(
    source_grp = file_h5[["spectra"]],
    source_names = paste0("S",spec_table_df$ScanNumber),
    dest_offsets = spec_table_df$offset,
    dsl= ctd1_spec))

file_h5$close_all()
alt_file1$close_all()
```

Benchmark
=========

To run the benchmark you’ll want to make sure you have `numpy`,`pandas`,
and `h5py` installed. I’ve tweaked the original benchmark code and put
it in `python/original.py` and the benchmark for the refactored data
model is in `python/refactor.py`. You can run the benchmark by doing
something like this:

``` shell
cd ./python && python3 benchmark.py
```


You can interact with the benchmark data using [snakeviz](https://jiffyclub.github.io/snakeviz/).  Once you've installed it you can run it on the command like like this:

``` shell
snakeviz profiles/refactor.stat
```


Kita
=========

To use Kita you'll need to sign up for an account on [Kita Lab](https://hdflab.hdfgroup.org/). Once you've logged on to kita lab you can locally install [h5pyd](https://github.com/HDFGroup/h5pyd) (or use the copy of h5pyd on kita lab).  If you are running on kita lab copy the `new_data/blacktea_hcd_7500_top3_1.h5` file onto the server using the icon that looks like this: ![Up Arrow Icon](profiles/UpArrow.png).  If you are running locally, you'll want to copy the file `~/.hscfg` from kita lab to your local machine (Go to `File` then `Open from path...` then type in `.hscfg` and copy the text (Shift + Right click -> Copy) and paste it in to a file at `~/.hscfg` (I have no idea if this works on Windows). Either way you can then issue the following command to load the HDF5 file to the HSDS/Kita service

``` shell
hsload new_data/blacktea_hcd_7500_top3_1.h5 /home/your_kita_username/blacktea_hcd_7500_top3_1.h5
```

The benchmark code might then look something like this. 

``` python
import kitaMZ
kita_file = "/home/your_kita_username/blacktea_hcd_7500_top3_1.h5"

def Benchmark_kita(fn, seed=123, N=100, scans=None):
    with h5py.File(fn, "r") as h5:
        print("Testing " + str(h5))
        kitaMZ.spectrumTest(h5=h5, N=N, seed=seed, scans=scans)
        kitaMZ.XICTest(h5, seed=seed, N=N, scans=scans)
        kitaMZ.XICFromSpectrumTest(h5, N=N, seed=seed)
    return True

pr_kita = cProfile.Profile()
pr_kita.enable()
Benchmark_kita(kita_file)
pr_kita.disable()
pr_kita.dump_stats('../profiles/kita.stat')
```


 You'll of course want to make sure you've cloned this repo onto Kita Lab if you're running on Kita Lab and also make sure that your working directory (or pythonpath) has `kitaMZ.py` in it.
