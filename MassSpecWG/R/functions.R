library(hdf5r)
library(stringr)
library(purrr)

ls_list <- function(node, recurse = TRUE) {
  node_ls <- node$ls(detailed = FALSE, recursive = FALSE)
  names <- node_ls$name
  types <- node_ls$obj_type
  if (!recurse | (all(types != "H5I_GROUP"))) {
    return(map(names, ~ node[[.x]]))
  }
  map2(names, types, function(name, type) {
    if (type == "H5I_GROUP") {
      return(ls_list(node[[name]], recurse = TRUE))
    }
    return(node[[name]])
  }) %>% set_names(names)
}


pytable_ls <- function(grp, id_pattern = "T([0-9]+)") {
  all_ds <- grp$ls(recursive = T)

  ds_df <- as_tibble(all_ds) %>%
    dplyr::select(
      -link.type,
      -group.nlinks,
      -group.mounted,
      -dataset.space_class
    ) %>%
    filter(obj_type == "H5I_DATASET") %>%
    dplyr::select(
      -obj_type,
      -committed_type,
      -num_attrs,
      -dataset.maxdims,
      -dataset.type_class,
      -dataset.rank
    ) %>%
    tidyr::separate(col = name, into = c("group", "name"), sep = "/") %>%
    dplyr::mutate(group = as.integer(str_replace(group, id_pattern, "\\1"))) %>%
    tidyr::spread(
      key = "name",
      value = "dataset.dims", convert = TRUE
    ) %>%
    arrange(group) %>%
    dplyr::select(group, N = axis1) %>%
    mutate(offset = lag(cumsum(N), default = 0))
  return(ds_df)
}

create_tbl_dataset <- function(dest_h5, name, prototype_df, table_df, chunksize = as.integer(mean(table_df$N)), gzip_level = NULL) {
  out_h5 <- make_file(dest_h5, mode = "r+")
  ct <- hdf5r:::guess_dtype(prototype_df)
  tot_size <- sum(table_df$N)
  sp <- H5S$new(type = "simple", dims = tot_size)
  dest_h5$create_dataset(
    name = name,
    dtype = ct,
    space = sp,
    chunk_dims = chunksize,
    gzip_level = gzip_level
  )
}

create_transposed_dataset <- function(dest_h5, name, prototype_df, table_df, chunksize = as.integer(mean(table_df$N)), gzip_level = NULL) {
  out_h5 <- make_file(dest_h5, mode = "r+")
  out_grp <- out_h5$create_group(name)
  ct <- map(prototype_df, hdf5r:::guess_dtype)
  tot_size <- sum(table_df$N)
  sp <- H5S$new(type = "simple", dims = tot_size)
  imap(ct, ~ out_grp$create_dataset(
    name = .y,
    dtype = .x,
    space = sp,
    chunk_dims = chunksize,
    gzip_level = gzip_level
  ))
}


get_names <- function(grp) {
  if (inherits(grp, "H5D")) {
    return(grp$get_type()$describe()$labels)
  }
  return(grp$ls()$name)
}

read_block_df <- function(items, values, extract_cols = NULL) {
  if (!is.null(extract_cols)) {
    cls <- items$read()
    cls_keep <- which(cls %in% extract_cols)
    if (length(cls_keep) == 0) {
      return(NULL)
    }
    dmt <- values$dims
    ret_df <- tibble::as_tibble(
      set_colnames(
        t(values$read(
          args = list(cls_keep, seq_len(dmt[2])),
          drop = FALSE
        )),
        cls[cls_keep]
      )
    )
    return(ret_df)
  }
  tibble::as_tibble(
    magrittr::set_colnames(
      t(values$read(drop = FALSE)),
      items$read()
    )
  )
}

copy_pytables_h5 <- function(source_grp, source_names, dest_offsets, dsl) {
  stopifnot(length(source_names) == length(dest_offsets))
  input_list <- list(
    source = source_names,
    offset = dest_offsets
  )

  pb <- progress_estimated(length(source_names))
  pwalk(input_list, function(source, offset, source_grp, dsl) {
    tdf <- read_pytables(source_grp[[source]])
    N <- nrow(tdf)
    pb$tick()$print()
    walk2(tdf, dsl, function(vec, ds, offset, N) {
      sl <- seq(from = offset + 1, length.out = N)
      ds[sl] <- vec
    }, offset = offset, N = N)
  }, source_grp = source_grp, dsl = dsl)
}



copy_pytables_tbl_h5 <- function(source_grp, source_names, dest_offsets, ids = seq_along(source_names), idname, ds) {
  input_list <- list(
    source = source_names,
    offset = dest_offsets,
    id = ids
  )
  pb <- progress_estimated(length(source_names))
  for (i in seq_along(source_names)) {
    source <- source_names[i]
    offset <- dest_offsets[i]
    id <- ids[i]
    pb$tick()$print()
    tdf <- read_pytables(source_grp[[source]]) %>%
      mutate({
        {
          idname
        }
      } := id)
    stopifnot(nrow(tdf) > 0)
    N <- nrow(tdf)
    sl <- seq(from = offset + 1, length.out = N)
    ds[sl] <- tdf
  }
}



make_group <- function(h5, name) {
  if (!h5$exists(name)) {
    return(h5$create_group(name))
  }
  grp <- h5[[name]]
  stopifnot(any("H5Group" %in% class(grp)))
  return(grp)
}


make_file <- function(path, mode = "r+") {
  if (inherits(path, "H5File")) {
    return(path)
  }
  if (!fs::file_exists(path)) {
    if (str_detect(mode, "r")) {
      warning(path, " does not exist, changing mode from ", mode, " to 'w'")
      mode <- "w"
    }
  }
  return(H5File$new(path, mode = mode))
}




write_df_h5 <- function(df, grp, append = FALSE) {
  ct <- map(df, typeof) %>% map(function(x) {
    if (x == "integer") {
      return(h5types$H5T_NATIVE_INT)
    } else {
      return(h5types[[x]])
    }
  })
  iwalk(df, function(x, y) {
    exst <- grp$exists(y)
    if (exst) {
      ds <- grp[[y]]
    } else {
      ds <- grp$create_dataset(name = y, robj = x, gzip_level = 1L)
    }
    if (append && exst) {
      odims <- ds$dims
      newdims <- odims + length(x)
      ds$set_extent(newdims)
      ds[(odims + 1):newdims] <- x
    } else {
      ds[1:length(x)] <- x
    }
  })
}

read_pytables <- function(grp, colnames = NULL) {
  grp_names <- names(grp)
  items_names <- map(grp_names[str_detect(grp_names, "_items")], ~ grp[[.x]])
  all_items <- map(items_names, ~ .x$read()) %>% flatten_chr()
  if (is.null(colnames)) {
    stopifnot(all(colnames %in% all_items))
  }
  values_names <- map(grp_names[str_detect(grp_names, "_values")], ~ grp[[.x]])
  ret_df <- list(items = items_names, values = values_names) %>%
    transpose() %>%
    purrr::discard(~ length(.x$values$dims) == 1) %>%
    map(~ read_block_df(.x$items, .x$values)) %>%
    compact() %>%
    bind_cols()
  return(ret_df)
}

read_df_h5 <- function(grp, names = get_names(grp), index = NULL) {
  if (inherits(grp, "H5D")) {
    if (is.null(index)) {
      return(dplyr::select(grp[], !!!names))
    }
    return(dplyr::select(grp[index], !!!names))
  }

  if (is.null(index)) {
    return(magrittr::set_colnames(purrr::map_dfc(names, ~ grp[[.x]][]), names))
  } else {
    return(magrittr::set_colnames(purrr::map_dfc(names, ~ grp[[.x]][index]), names))
  }
}
