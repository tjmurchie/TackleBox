#!/usr/bin/env Rscript
# MetaMerge heatmap / stacked-bar report renderer
# Version 1.1.0

suppressPackageStartupMessages({
  pkgs <- c("readr", "dplyr", "tidyr", "ggplot2", "patchwork",
            "stringr", "forcats", "scales", "tibble", "ggtext")
  to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(to_install) > 0) {
    message("Installing missing packages: ", paste(to_install, collapse = ", "))
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[idx + 1]
}

input_dir   <- get_arg("--input-dir")
out_dir     <- get_arg("--outdir", "metamerge_reports")
top_n       <- as.integer(get_arg("--top-n", "0"))
bar_top_n   <- as.integer(get_arg("--bar-top-n", "40"))
page_rows   <- as.integer(get_arg("--page-rows", "28"))
order_file  <- get_arg("--sample-order")
plaus_file  <- get_arg("--plausibility")

plausibility_tbl <- NULL
if (!is.null(plaus_file) && file.exists(plaus_file)) {
  plausibility_tbl <- readr::read_csv(plaus_file, show_col_types = FALSE) %>%
    dplyr::select(scientific_name, plausibility)
  message("Plausibility table loaded: ", nrow(plausibility_tbl), " taxa")
}

if (is.null(input_dir)) {
  stop("Usage: Rscript metamerge_heatmap_report.R --input-dir path/to/report_inputs --outdir path/to/reports [--top-n 0] [--bar-top-n 20] [--page-rows 28] [--sample-order sample_order.csv]")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load sample order file (if supplied)
# ---------------------------------------------------------------------------
# order_map: long-format table with columns plot_group_new, plot_sample_id,
# sort_order.  Built by pivoting the wide CSV (columns = groups, rows = IDs).
order_map <- NULL
if (!is.null(order_file)) {
  order_tbl <- readr::read_csv(order_file,
                               col_types = readr::cols(.default = "c"),
                               show_col_types = FALSE)
  order_map <- order_tbl %>%
    tidyr::pivot_longer(dplyr::everything(),
                        names_to  = "plot_group_new",
                        values_to = "plot_sample_id") %>%
    dplyr::filter(!is.na(plot_sample_id),
                  nchar(trimws(plot_sample_id)) > 0) %>%
    dplyr::mutate(plot_sample_id = trimws(plot_sample_id)) %>%
    dplyr::group_by(plot_group_new) %>%
    dplyr::mutate(sort_order = dplyr::row_number()) %>%
    dplyr::ungroup()
  # Samples may intentionally appear in multiple columns (duplication across
  # groups is supported).  Do NOT deduplicate here.
  message("Sample order file loaded: ", nrow(order_map), " entries across ",
          length(unique(order_map$plot_group_new)), " groups (",
          length(unique(order_map$plot_sample_id)), " unique samples)")
}

# ---------------------------------------------------------------------------
# Bins, colours, shape maps
# ---------------------------------------------------------------------------
count_breaks <- c(-0.1, 0.5, 5, 10, 50, 100, 500, 1000, 5000, Inf)
count_labels <- c("0", "1-4", "5-9", "10-49", "50-99", "100-499", "500-999", "1,000-4,999", ">=5,000")
count_palette <- c(
  "0"           = "#f4efcf",
  "1-4"         = "#dbe8a4",
  "5-9"         = "#b8df8a",
  "10-49"       = "#7fcdbb",
  "50-99"       = "#41b6c4",
  "100-499"     = "#2c7fb8",
  "500-999"     = "#253494",
  "1,000-4,999" = "#1a237e",
  ">=5,000"     = "#0d164d"
)
dark_bins <- c("100-499", "500-999", "1,000-4,999", ">=5,000")

status_priority_map <- c(
  "Very high confidence" = 0, "High confidence" = 1, "Supported" = 2,
  "Tentative" = 3, "Weak support" = 4, "Blank-associated" = 5
)
tax_rank_order_map <- c(
  "family" = 0, "subfamily" = 1, "tribe" = 2, "subtribe" = 3,
  "genus" = 4, "subgenus" = 5, "species group" = 6, "species" = 7,
  "subspecies" = 8
)
lib_support_levels <- c(
  "Very high confidence", "High confidence", "Supported",
  "Tentative", "Weak support", "Blank-associated",
  "Environmental-control"
)
# Unicode shapes with decreasing polygon-face count for visual confidence hierarchy.
# Rendered via cairo_pdf + DejaVu Sans (all glyphs confirmed present in that font).
lib_shape_map <- c(
  "Very high confidence"  = "\u2B22",  # ⬢ black hexagon        (6 faces)
  "High confidence"       = "\u2B1F",  # ⬟ black pentagon       (5 faces)
  "Supported"             = "\u25A0",  # ■ black square         (4 faces)
  "Tentative"             = "\u25C6",  # ◆ black diamond        (4 rotated)
  "Weak support"          = "\u25B2",  # ▲ black triangle up    (3 faces)
  "Blank-associated"      = "\u2715",  # ✕ multiplication X     (control/blank)
  "Environmental-control" = "\u25CB"   # ○ open circle          (env. sample)
)

# ---------------------------------------------------------------------------
# Ecological plausibility (supplied via --plausibility CSV)
# ---------------------------------------------------------------------------
# pch codes: 16=filled circle, 18=filled diamond, 2=open triangle up, 13=circle-X.
# Chosen to conceptually parallel aDNA support (more filled = more likely) while remaining
# visually distinct: open △ (pch=2) differs from filled ▲ (Weak support); ⊗ (pch=13) differs from ✕ (Blank-associated).
plaus_cat_levels <- c("Very likely present", "Plausible", "Unlikely", "Likely false positive")
plaus_cat_shapes <- stats::setNames(c(16L, 18L, 2L, 13L),                         plaus_cat_levels)
plaus_cat_colors <- stats::setNames(c("#2ca02c", "#1f77b4", "#e6a500", "#CC0000"), plaus_cat_levels)

add_plausibility <- function(label_df) {
  if (is.null(plausibility_tbl)) return(label_df)
  label_df %>%
    dplyr::left_join(
      plausibility_tbl %>% dplyr::mutate(
        plaus_cat = dplyr::case_when(
          plausibility == "very_likely"           ~ "Very likely present",
          plausibility == "plausible"             ~ "Plausible",
          plausibility == "unlikely"              ~ "Unlikely",
          plausibility == "likely_false_positive" ~ "Likely false positive",
          TRUE ~ NA_character_
        )
      ) %>% dplyr::select(scientific_name, plaus_cat),
      by = "scientific_name"
    ) %>%
    dplyr::mutate(plaus_cat = factor(plaus_cat, levels = plaus_cat_levels))
}

make_plaus_panel <- function(label_df, boundaries, shared_y, base_theme) {
  if (is.null(plausibility_tbl) || !"plaus_cat" %in% names(label_df)) return(NULL)

  # Mirror the aDNA support legend approach exactly:
  #   real points  → colour via identity column, show.legend = FALSE
  #   legend dummy → shape only (no colour in aes), show.legend = TRUE
  #   colour forced via override.aes on the single shape guide
  # This avoids ggplot2 4.x merged-colour+shape guide deduplication which can
  # silently drop plausibility levels absent from the current page's data.
  label_plot <- label_df %>%
    dplyr::mutate(plaus_color = dplyr::coalesce(
      unname(plaus_cat_colors[as.character(plaus_cat)]), NA_character_
    ))

  # One row per plausibility level — shape only, alpha=0, show.legend=TRUE.
  # Same anchor pattern as the aDNA support legend_df.
  plaus_legend_df <- data.frame(
    taxon_key = rep(label_df$taxon_key[[1]], 4),
    x         = rep(0, 4),
    plaus_cat = factor(plaus_cat_levels, levels = plaus_cat_levels)
  )

  ggplot(label_plot, aes(x = 0, y = taxon_key)) +
    geom_hline(data = boundaries, aes(yintercept = yint),
               inherit.aes = FALSE, linewidth = 0.3, colour = "grey30") +
    geom_point(aes(shape = plaus_cat, colour = plaus_color),
               size = 2.5, na.rm = TRUE, show.legend = FALSE) +
    geom_point(data = plaus_legend_df,
               aes(x = x, y = taxon_key, shape = plaus_cat),
               inherit.aes = FALSE, alpha = 0, size = 3.0, show.legend = TRUE) +
    scale_colour_identity() +
    scale_shape_manual(
      values = plaus_cat_shapes, breaks = plaus_cat_levels,
      limits = plaus_cat_levels, drop = FALSE,
      name   = "Ecological Plausibility\n(E. Ontario ~43.5°N)"
    ) +
    guides(
      shape = guide_legend(order = 3, override.aes = list(
        alpha  = 1,
        size   = 3,
        colour = unname(plaus_cat_colors)
      ))
    ) +
    shared_y +
    coord_cartesian(xlim = c(-0.1, 0.1), clip = "off") +
    labs(title = " ", subtitle = " ") +
    base_theme +
    theme(legend.position = "right")
}

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

save_plot_safe <- function(plot_obj, out_file, width, height, dpi = 300) {
  grDevices::cairo_pdf(out_file, width = min(width, 24), height = min(height, 18),
                       family = "DejaVu Sans")
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot_obj)
}

save_plot_pages <- function(plot_list, out_file, width, height) {
  grDevices::cairo_pdf(out_file, width = min(width, 24), height = min(height, 18),
                       family = "DejaVu Sans", onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (p in plot_list) print(p)
}

pretty_common_name <- function(x) {
  x <- ifelse(is.na(x) | x == "", " ", x)
  x <- stringr::str_squish(x)
  x <- ifelse(x == " ", x, stringr::str_to_title(tolower(x)))
  x
}

make_taxon_key <- function(scientific_name, tax_rank, display_group) {
  paste(scientific_name %||% "", tax_rank %||% "", display_group %||% "", sep = "|||")
}

normalize_support_col <- function(dat) {
  if (!"aDNA_support_status" %in% names(dat) && "DNA_support_status" %in% names(dat))
    dat <- dat %>% rename(aDNA_support_status = DNA_support_status)
  dat
}

subset_label <- function(group_name) {
  group_name <- as.character(group_name %||% "")
  if (group_name == "" || is.na(group_name)) return("subset")
  stringr::str_replace_all(group_name, "[^A-Za-z0-9_-]", "_")
}

# ---------------------------------------------------------------------------
# resolve_x_order: build x-axis levels + display label map for one plot group
# ---------------------------------------------------------------------------
# grp_dat    — data already filtered (or filterable) to a single plot_group
# subset_name — the plot_group name being rendered
#
# Returns a named character vector: names = plot_sample_id (factor levels),
# values = display labels (with "(xN)" suffix when N > 1 merged libraries).
# Order follows order_map when supplied, otherwise the default sort.
resolve_x_order <- function(grp_dat, subset_name) {
  present_ids <- unique(grp_dat$plot_sample_id)

  # Per-sample replicate count for label suffix.
  if ("n_technical_libraries" %in% names(grp_dat)) {
    n_tech_df <- grp_dat %>%
      dplyr::group_by(plot_sample_id) %>%
      dplyr::summarise(
        n_tech         = suppressWarnings(max(as.integer(n_technical_libraries), na.rm = TRUE)),
        plot_lib_label = dplyr::first(stats::na.omit(as.character(plot_library_label))),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        n_tech = dplyr::if_else(is.na(n_tech) | is.infinite(n_tech), 1L, n_tech)
      )
  } else {
    n_tech_df <- grp_dat %>%
      dplyr::distinct(plot_sample_id, plot_library_label) %>%
      dplyr::rename(plot_lib_label = plot_library_label) %>%
      dplyr::mutate(n_tech = 1L)
  }

  n_tech_df <- n_tech_df %>%
    dplyr::mutate(
      display_label = dplyr::if_else(
        n_tech > 1L,
        paste0(dplyr::coalesce(plot_lib_label, plot_sample_id), " (x", n_tech, ")"),
        dplyr::coalesce(plot_lib_label, plot_sample_id)
      )
    )

  # Determine ordered x_levels.
  if (!is.null(order_map)) {
    # Use ALL samples listed in order_map for this group — including zero-detection
    # samples that are absent from grp_dat.  Those will appear as empty columns.
    ord <- order_map %>%
      dplyr::filter(plot_group_new == subset_name) %>%
      dplyr::arrange(sort_order)
    x_levels <- as.character(ord$plot_sample_id)
  } else {
    level_df <- grp_dat %>%
      dplyr::distinct(plot_sample_id, plot_library_label, is_negative_control,
                      is_environmental_control, sample_type) %>%
      dplyr::mutate(
        is_neg = dplyr::if_else(is.na(is_negative_control), FALSE,
                                as.logical(is_negative_control)),
        is_pos = grepl("positive", tolower(as.character(sample_type))),
        is_env = dplyr::if_else(is.na(is_environmental_control), FALSE,
                                as.logical(is_environmental_control))
      ) %>%
      dplyr::arrange(is_neg, is_pos, is_env, plot_sample_id, plot_library_label)
    x_levels <- as.character(level_df$plot_sample_id)
  }

  label_map <- stats::setNames(
    as.character(n_tech_df$display_label[match(x_levels, n_tech_df$plot_sample_id)]),
    x_levels
  )
  # Fall back to plot_sample_id for any unmatched entries.
  missing <- is.na(label_map)
  label_map[missing] <- names(label_map)[missing]
  label_map
}

bin_counts <- function(x) {
  factor(cut(x, breaks = count_breaks, labels = count_labels,
             include.lowest = TRUE, right = FALSE), levels = count_labels)
}

label_contrast_color  <- function(cb) { if (!is.na(cb) && as.character(cb) %in% dark_bins) "white" else "black" }
symbol_contrast_color <- function(cb) { if (!is.na(cb) && as.character(cb) %in% dark_bins) "white" else "black" }

map_library_support <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x == "Damage-supported (>=100 reads)" ~ "Very high confidence",
    x == "Damage-supported"               ~ "High confidence",
    x == "Lineage-supported"              ~ "Supported",
    x == "Count-only"                     ~ "Tentative",
    x == "Blank-library"                  ~ "Blank-associated",
    x == "Environmental-control"          ~ "Environmental-control",
    TRUE ~ NA_character_
  )
}

# Panel width: sized tightly to the longest text on the current page.
panel_width <- function(x, min_width = 0.35, max_width = 3.0, char_scale = 0.055) {
  x <- x[!is.na(x) & nchar(as.character(x)) > 0]
  if (length(x) == 0) return(min_width)
  mx <- max(nchar(as.character(x)), na.rm = TRUE)
  max(min_width, min(max_width, 0.10 + mx * char_scale))
}

parse_tax_path_tbl <- function(path) {
  if (is.na(path) || path == "") return(tibble::tibble(name = character(), rank = character()))
  parts <- unlist(strsplit(path, "\t"))
  nm <- stringr::str_match(parts, "^[0-9]+:(.*?):\"[^\"]+\"$")[, 2]
  rk <- tolower(stringr::str_match(parts, "^[0-9]+:.*?:\"([^\"]+)\"$")[, 2])
  tibble::tibble(name = nm, rank = rk)
}

# Determine a human-readable grouping label from a taxon's rank and tax_path.
# Key fix v0.5.7: subfamilies and tribes are grouped under their parent
# FAMILY (not under themselves), so e.g. Bovinae appears within Bovidae.
infer_display_group_single <- function(scientific_name, tax_rank, tax_path, display_group) {
  tr  <- tolower(as.character(tax_rank  %||% ""))
  sci <- as.character(scientific_name %||% "")

  # Family-level taxa are their own group label.
  if (tr == "family") return(sci)

  # Subfamily / tribe / subtribe: look up parent family from path.
  # This ensures e.g. Bovinae groups with Bison under Bovidae.
  if (tr %in% c("subfamily", "tribe", "subtribe")) {
    tbl <- parse_tax_path_tbl(tax_path)
    if (nrow(tbl) > 0) {
      fam <- tbl %>% filter(rank == "family") %>% slice_head(n = 1) %>% pull(name)
      if (length(fam) == 1 && !is.na(fam) && fam != "") return(fam)
    }
    return(sci)
  }

  # For genus / species and below: search path for family → subfamily → tribe → order.
  tbl <- parse_tax_path_tbl(tax_path)
  if (nrow(tbl) > 0) {
    for (rk in c("family", "subfamily", "tribe", "order")) {
      hit <- tbl %>% filter(rank == rk) %>% slice_head(n = 1) %>% pull(name)
      if (length(hit) == 1 && !is.na(hit) && hit != "") return(hit)
    }
  }

  if (!is.na(display_group) && display_group != "") return(as.character(display_group))
  sci
}

with_inferred_group <- function(dat) {
  dat %>%
    rowwise() %>%
    mutate(display_group = infer_display_group_single(scientific_name, tax_rank, tax_path, display_group)) %>%
    ungroup()
}

prepare_subset <- function(dat, subset_name, top_n = 0, page_rows = 28) {
  dat <- normalize_support_col(dat) %>%
    filter(plot_group == subset_name) %>%
    mutate(
      display_group       = ifelse(is.na(display_group) | display_group == "", scientific_name, display_group),
      aDNA_support_status = factor(aDNA_support_status, levels = names(status_priority_map)),
      common_name         = pretty_common_name(common_name),
      count               = as.numeric(count),
      taxon_key           = make_taxon_key(scientific_name, tax_rank, display_group)
    ) %>%
    filter(!is.na(count), count > 0)
  if (nrow(dat) == 0) return(NULL)

  dat <- dat %>%
    mutate(
      support_priority = dplyr::coalesce(.data$support_priority,
                           unname(status_priority_map[as.character(aDNA_support_status)])),
      tax_rank_order   = dplyr::coalesce(.data$tax_rank_order,
                           unname(tax_rank_order_map[tolower(as.character(tax_rank))]), 999)
    ) %>%
    group_by(plot_group, taxon_key) %>%
    mutate(
      subset_taxon_sum = dplyr::coalesce(.data$subset_taxon_sum, sum(count, na.rm = TRUE)),
      subset_taxon_max = dplyr::coalesce(.data$subset_taxon_max, max(count, na.rm = TRUE))
    ) %>%
    ungroup() %>%
    group_by(plot_group, display_group) %>%
    mutate(
      display_group_sum              = dplyr::coalesce(.data$display_group_sum, sum(count, na.rm = TRUE)),
      display_group_max              = dplyr::coalesce(.data$display_group_max, max(count, na.rm = TRUE)),
      display_group_status_priority  = dplyr::coalesce(.data$display_group_status_priority,
                                         min(support_priority, na.rm = TRUE))
    ) %>%
    ungroup()

  taxon_tbl <- dat %>%
    group_by(taxon_key) %>%
    summarise(
      scientific_name                = dplyr::first(scientific_name),
      display_group                  = dplyr::first(display_group),
      common_name                    = dplyr::first(common_name),
      tax_rank                       = dplyr::first(tax_rank),
      support_priority               = suppressWarnings(min(support_priority,              na.rm = TRUE)),
      subset_taxon_sum               = suppressWarnings(max(subset_taxon_sum,              na.rm = TRUE)),
      subset_taxon_max               = suppressWarnings(max(subset_taxon_max,              na.rm = TRUE)),
      display_group_sum              = suppressWarnings(max(display_group_sum,             na.rm = TRUE)),
      display_group_max              = suppressWarnings(max(display_group_max,             na.rm = TRUE)),
      display_group_status_priority  = suppressWarnings(min(display_group_status_priority, na.rm = TRUE)),
      tax_rank_order                 = suppressWarnings(min(tax_rank_order,                na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      support_priority              = ifelse(is.infinite(support_priority),             999, support_priority),
      display_group_status_priority = ifelse(is.infinite(display_group_status_priority),999, display_group_status_priority),
      tax_rank_order                = ifelse(is.infinite(tax_rank_order),               999, tax_rank_order)
    ) %>%
    arrange(display_group_status_priority, desc(display_group_sum), desc(display_group_max),
            display_group, support_priority, desc(subset_taxon_sum), desc(subset_taxon_max),
            tax_rank_order, scientific_name)

  if (top_n > 0 && nrow(taxon_tbl) > top_n) taxon_tbl <- taxon_tbl %>% slice_head(n = top_n)
  taxon_tbl <- taxon_tbl %>% mutate(row_id = row_number())
  split(taxon_tbl$taxon_key, ceiling(taxon_tbl$row_id / page_rows))
}

# ---------------------------------------------------------------------------
# make_heatmap
# ---------------------------------------------------------------------------
make_heatmap <- function(dat, subset_name, out_file, top_n = 0, page_rows = 28) {
  dat <- normalize_support_col(dat)
  dat <- with_inferred_group(dat) %>%
    mutate(
      aDNA_support_status = factor(aDNA_support_status, levels = names(status_priority_map)),
      common_name         = pretty_common_name(common_name),
      count               = as.numeric(count),
      taxon_key           = make_taxon_key(scientific_name, tax_rank, display_group)
    )

  if (identical(subset_name, "__ALL__")) {
    dat <- dat %>% mutate(plot_group = "all_libraries")
    subset_name <- "all_libraries"
  }

  pages <- prepare_subset(dat, subset_name, top_n = top_n, page_rows = page_rows)
  if (is.null(pages)) return(invisible(NULL))

  if (!"is_environmental_control" %in% names(dat)) dat$is_environmental_control <- FALSE

  x_label_map <- resolve_x_order(dat %>% filter(plot_group == subset_name), subset_name)
  x_levels    <- names(x_label_map)
  x_labels    <- function(x) x_label_map[x]

  if (length(x_levels) == 0) return(invisible(NULL))

  # Compute max page size before loop so shared_y can pad shorter pages to
  # match — this keeps the physical tile height identical on every page.
  max_rows <- max(vapply(pages, length, numeric(1)))

  plots <- list()
  idx   <- 1
  for (page_taxa in pages) {
    raw_page_dat <- dat %>%
      filter(plot_group == subset_name, taxon_key %in% page_taxa) %>%
      mutate(
        common_name = pretty_common_name(common_name),
        count       = as.numeric(count),
        taxon_key   = make_taxon_key(scientific_name, tax_rank, display_group)
      ) %>%
      filter(!is.na(count), count > 0)

    label_df <- raw_page_dat %>%
      group_by(taxon_key) %>%
      summarise(
        scientific_name               = dplyr::first(scientific_name),
        display_group                 = dplyr::first(display_group),
        common_name                   = dplyr::first(common_name),
        display_group_status_priority = suppressWarnings(min(display_group_status_priority, na.rm = TRUE)),
        display_group_sum             = suppressWarnings(max(display_group_sum,             na.rm = TRUE)),
        display_group_max             = suppressWarnings(max(display_group_max,             na.rm = TRUE)),
        support_priority              = suppressWarnings(min(support_priority,              na.rm = TRUE)),
        subset_taxon_sum              = suppressWarnings(max(subset_taxon_sum,              na.rm = TRUE)),
        subset_taxon_max              = suppressWarnings(max(subset_taxon_max,              na.rm = TRUE)),
        tax_rank_order                = suppressWarnings(min(tax_rank_order,                na.rm = TRUE)),
        tax_rank                      = dplyr::first(tax_rank),
        .groups = "drop"
      ) %>%
      mutate(
        display_group_status_priority = ifelse(is.infinite(display_group_status_priority), 999, display_group_status_priority),
        display_group_sum  = ifelse(is.infinite(display_group_sum),  0, display_group_sum),
        display_group_max  = ifelse(is.infinite(display_group_max),  0, display_group_max),
        support_priority   = ifelse(is.infinite(support_priority),   999, support_priority),
        subset_taxon_sum   = ifelse(is.infinite(subset_taxon_sum),   0, subset_taxon_sum),
        subset_taxon_max   = ifelse(is.infinite(subset_taxon_max),   0, subset_taxon_max),
        tax_rank_order     = ifelse(is.infinite(tax_rank_order),     999, tax_rank_order)
      ) %>%
      arrange(display_group_status_priority, desc(display_group_sum), desc(display_group_max),
              display_group, support_priority, desc(subset_taxon_sum), desc(subset_taxon_max),
              tax_rank_order, scientific_name) %>%
      mutate(
        row_id               = row_number(),
        display_group_label  = ifelse(duplicated(display_group), "", display_group)
      )

    label_df <- add_plausibility(label_df)
    y_levels <- label_df$taxon_key

    page_dat <- tidyr::expand_grid(plot_sample_id = x_levels, taxon_key = y_levels) %>%
      left_join(raw_page_dat, by = c("plot_sample_id", "taxon_key"), relationship = "many-to-many") %>%
      group_by(plot_sample_id, taxon_key) %>%
      summarise(
        scientific_name          = dplyr::first(stats::na.omit(scientific_name)),
        display_group            = dplyr::first(stats::na.omit(display_group)),
        common_name              = dplyr::first(stats::na.omit(common_name)),
        tax_rank                 = dplyr::first(stats::na.omit(tax_rank)),
        count                    = sum(count, na.rm = TRUE),
        raw_library_adna_support = dplyr::first(stats::na.omit(as.character(library_adna_support))),
        is_negative_control      = dplyr::first(stats::na.omit(is_negative_control)),
        sample_type              = dplyr::first(stats::na.omit(sample_type)),
        .groups = "drop"
      ) %>%
      mutate(
        count      = replace(count, is.na(count), 0),
        # Explicit factor with ALL 9 levels prevents ggplot2 from silently
        # dropping unused bins from the legend even with drop=FALSE.
        count_bin  = factor(bin_counts(count), levels = count_labels),
        count_text = ifelse(count > 0, scales::comma(round(count, 1)), ""),
        count_text_color         = vapply(as.character(count_bin), label_contrast_color,  character(1)),
        library_support_display  = map_library_support(raw_library_adna_support),
        symbol_color             = ifelse(
          map_library_support(raw_library_adna_support) == "Blank-associated",
          "#CC0000",   # red for blank-associated — makes contamination signal visible
          vapply(as.character(count_bin), symbol_contrast_color, character(1))
        )
      ) %>%
      left_join(label_df %>% select(taxon_key, row_id), by = "taxon_key") %>%
      mutate(
        taxon_key               = factor(taxon_key,               levels = rev(y_levels)),
        plot_sample_id          = factor(plot_sample_id,          levels = x_levels),
        library_support_display = factor(library_support_display, levels = lib_support_levels),
        count_bin               = factor(count_bin,               levels = count_labels)
      )

    # Dummy data containing one row for EVERY count bin level.  This is the
    # most reliable way to force all 9 bins into the fill legend regardless of
    # which bins happen to be present on this specific page.
    fill_dummy_df <- tibble::tibble(
      plot_sample_id = factor(rep(x_levels[1], 9), levels = x_levels),
      taxon_key      = factor(rep(rev(y_levels)[1], 9), levels = rev(y_levels)),
      count_bin      = factor(count_labels, levels = count_labels)
    )

    n_rows     <- nrow(label_df)
    boundaries <- label_df %>%
      group_by(display_group) %>%
      summarise(group_end = max(row_id), .groups = "drop") %>%
      filter(group_end < n_rows) %>%
      mutate(yint = n_rows - group_end + 0.5)

    symbol_dat <- page_dat %>% filter(!is.na(library_support_display))

    # Per-page column widths — sized to longest text on THIS page only.
    group_width  <- panel_width(label_df$display_group_label, min_width = 0.35, max_width = 1.8, char_scale = 0.055)
    sci_width    <- panel_width(label_df$scientific_name,      min_width = 0.65, max_width = 2.8, char_scale = 0.060)
    common_width <- panel_width(label_df$common_name,          min_width = 0.40, max_width = 2.4, char_scale = 0.053)
    # Heat panel width: proportional to library count; no excessive minimum.
    n_libs       <- length(x_levels)
    heat_width   <- max(0.65, n_libs * 0.52)

    # Shared y-scale — identical limits and expand across all panels ensure
    # perfect alignment of the geom_hline group-boundary lines.
    # extra_expand pads shorter pages so every page has the same physical
    # row height and tiles look consistent across all pages.
    n_page_rows  <- length(page_taxa)
    extra_expand <- max(0, (max_rows - n_page_rows)) / 2
    shared_y <- scale_y_discrete(limits = rev(y_levels), expand = expansion(add = 0.4 + extra_expand))

    # Consistent title line on ALL panels so patchwork aligns data regions.
    # Label panels: blank title + column name as subtitle.
    # Heat panel: subset name as title + "Library" as subtitle.
    # Equal title overhead = equal vertical offset before data region starts.
    title_size    <- 8.0   # pt — same for all panel titles
    subtitle_size <- 7.0   # pt — same for all panel subtitles

    # Theme shared by the three left-hand text columns.
    # axis.text.x reserves the same space as the heat panel's rotated labels
    # (colour=NA makes it invisible) so patchwork aligns data areas correctly.
    label_theme <- theme(
      axis.text.y      = element_blank(),
      axis.ticks       = element_blank(),
      axis.title       = element_blank(),
      # Invisible x-axis text: same size/angle as heat panel → equal bottom margin.
      axis.text.x      = element_text(colour = NA, size = 7, angle = 45, hjust = 1, vjust = 1),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid       = element_blank(),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.title       = element_text(size = title_size,    face = "plain", hjust = 0,
                                      margin = margin(b = 1)),
      plot.subtitle    = element_text(size = subtitle_size, face = "bold",  hjust = 0,
                                      margin = margin(b = 2)),
      plot.margin      = margin(5.5, 2, 5.5, 2)
    )

    # NOTE: y aesthetic uses taxon_key as plain character (not factor) to avoid
    # the lazy-quosure page-ordering bug described in previous versions.
    group_panel <- ggplot(label_df, aes(x = 0, y = taxon_key)) +
      geom_hline(data = boundaries, aes(yintercept = yint),
                 inherit.aes = FALSE, linewidth = 0.3, colour = "grey30") +
      geom_text(aes(label = display_group_label), hjust = 0, size = 2.7) +
      shared_y +
      coord_cartesian(xlim = c(0, 1), clip = "off") +
      labs(title = " ", subtitle = "Family / Group") +
      label_theme

    sci_panel <- ggplot(label_df, aes(x = 0, y = taxon_key)) +
      geom_hline(data = boundaries, aes(yintercept = yint),
                 inherit.aes = FALSE, linewidth = 0.3, colour = "grey30") +
      geom_text(aes(label = scientific_name), hjust = 0, size = 2.65, fontface = "italic") +
      shared_y +
      coord_cartesian(xlim = c(0, 1), clip = "off") +
      labs(title = " ", subtitle = "Scientific Name") +
      label_theme

    common_panel <- ggplot(label_df, aes(x = 0, y = taxon_key)) +
      geom_hline(data = boundaries, aes(yintercept = yint),
                 inherit.aes = FALSE, linewidth = 0.3, colour = "grey30") +
      geom_text(aes(label = common_name), hjust = 0, size = 2.55) +
      shared_y +
      coord_cartesian(xlim = c(0, 1), clip = "off") +
      labs(title = " ", subtitle = "Common Name") +
      label_theme

    plaus_panel <- make_plaus_panel(label_df, boundaries, shared_y, label_theme)

    legend_df <- tibble::tibble(
      plot_sample_id          = factor(x_levels[1], levels = x_levels),
      taxon_key               = factor(rev(y_levels)[1], levels = rev(y_levels)),
      library_support_display = factor(lib_support_levels, levels = lib_support_levels)
    )

    heat_panel <- ggplot(page_dat, aes(x = plot_sample_id, y = taxon_key)) +
      geom_hline(data = boundaries, aes(yintercept = yint),
                 inherit.aes = FALSE, linewidth = 0.3, colour = "grey30") +
      # Invisible dummy tiles to force ALL 9 count bins into the legend.
      geom_tile(data = fill_dummy_df,
                aes(x = plot_sample_id, y = taxon_key, fill = count_bin),
                alpha = 0, width = 0, height = 0, colour = NA,
                inherit.aes = FALSE, show.legend = NA, na.rm = TRUE) +
      # Real tiles — rectangular (height ~1/3 the column width).
      geom_tile(aes(fill = count_bin),
                colour = "grey84", linewidth = 0.18,
                width = 0.96, height = 0.88, na.rm = TRUE) +
      # Support symbols — nudged to upper-right; size reduced for short tiles.
      geom_point(data = symbol_dat,
                 aes(shape = library_support_display, colour = symbol_color),
                 position = position_nudge(x = 0.37, y = 0.18),
                 size = 1.6, stroke = 0.2, show.legend = FALSE, na.rm = TRUE) +
      # Invisible dummy points for shape legend.
      geom_point(data = legend_df,
                 aes(x = plot_sample_id, y = taxon_key, shape = library_support_display),
                 inherit.aes = FALSE, alpha = 0, size = 3.0, show.legend = TRUE) +
      # Count text labels.
      geom_text(aes(x = plot_sample_id, y = taxon_key,
                    label = count_text, colour = count_text_color),
                size = 2.1, fontface = "plain", na.rm = TRUE,
                inherit.aes = FALSE, show.legend = FALSE) +
      scale_fill_manual(
        values   = count_palette,
        breaks   = count_labels,
        limits   = count_labels,
        labels   = count_labels,
        drop     = FALSE,
        na.value = count_palette[["0"]],
        name     = "MEGAN Count of Unique Reads\n(stepped square-root bins)"
      ) +
      scale_shape_manual(
        values = lib_shape_map, breaks = lib_support_levels,
        limits = lib_support_levels, drop = FALSE,
        name   = "Per-library aDNA support"
      ) +
      scale_colour_identity() +
      guides(
        fill   = guide_legend(order = 1, override.aes = list(alpha = 1, size = 4, shape = NA, colour = NA)),
        shape  = guide_legend(order = 2, override.aes = list(alpha = 1, size = 2.8,
                             colour = c("#000000","#000000","#000000","#000000","#000000","#CC0000","#000000"))),
        colour = "none"
      ) +
      shared_y +
      scale_x_discrete(labels = x_labels, drop = FALSE) +
      labs(
        x        = NULL, y = NULL,
        title    = paste0(subset_name, if (length(pages) > 1) paste0(" (page ", idx, "/", length(pages), ")") else ""),
        subtitle = "Library"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid      = element_blank(),
        axis.text.y     = element_blank(),
        axis.title      = element_blank(),
        axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
        legend.position = "right",
        # Single outer border only — no inner per-section boxes.
        legend.background     = element_blank(),
        legend.box.background = element_rect(colour = "grey40", fill = "white", linewidth = 0.5),
        legend.box.margin     = margin(4, 4, 4, 4),
        legend.spacing.y      = unit(0.1, "cm"),
        legend.key            = element_rect(fill = "white", colour = NA),
        legend.key.height     = unit(0.50, "cm"),
        legend.key.width      = unit(0.50, "cm"),
        legend.text           = element_text(size = 7),
        legend.title          = element_text(size = 7.5, face = "bold"),
        plot.title    = element_text(face = "plain", size = title_size,    margin = margin(b = 1)),
        plot.subtitle = element_text(face = "bold",  size = subtitle_size, colour = "grey35",
                                     margin = margin(b = 2)),
        plot.margin   = margin(5.5, 5.5, 5.5, 2)
      )

    if (!is.null(plaus_panel)) {
      plots[[idx]] <- group_panel + sci_panel + plaus_panel + common_panel + heat_panel +
        patchwork::plot_layout(widths  = c(group_width, sci_width, 0.12, common_width, heat_width),
                               guides  = "collect")
    } else {
      plots[[idx]] <- group_panel + sci_panel + common_panel + heat_panel +
        patchwork::plot_layout(widths = c(group_width, sci_width, common_width, heat_width))
    }
    idx <- idx + 1
  }

  plaus_w <- if (!is.null(plausibility_tbl)) 0.12 else 0.0
  width  <- max(5.5, group_width + sci_width + common_width + plaus_w + heat_width + 2.8)
  height <- max(4.0, 1.4 + max_rows * 0.16)
  save_plot_pages(plots, out_file, width = width, height = height)
}

# ---------------------------------------------------------------------------
# make_stacked_bars
# ---------------------------------------------------------------------------
make_stacked_bars <- function(dat, subset_name, out_file, bar_top_n = 20) {
  dat <- normalize_support_col(dat)
  dat <- with_inferred_group(dat) %>%
    filter(plot_group == subset_name) %>%
    mutate(count = as.numeric(count)) %>%
    filter(!is.na(count), count > 0)
  if (nrow(dat) == 0) return(invisible(NULL))

  if (!"support_priority" %in% names(dat))
    dat <- dat %>% mutate(support_priority = unname(status_priority_map[as.character(aDNA_support_status)]))

  dat <- dat %>%
    group_by(plot_group, scientific_name) %>%
    mutate(
      subset_taxon_sum = dplyr::coalesce(.data$subset_taxon_sum, sum(count, na.rm = TRUE)),
      subset_taxon_max = dplyr::coalesce(.data$subset_taxon_max, max(count, na.rm = TRUE))
    ) %>%
    ungroup()

  # Meaningful ranks for stacked bars: exclude uninformative high-rank nodes
  # (clades, kingdoms, phyla) that aggregate reads from many organisms.
  bar_meaningful_ranks <- c("species", "genus", "family", "subfamily",
                            "tribe", "subtribe", "order", "suborder",
                            "superfamily", "infraorder")

  taxa <- dat %>%
    filter(tolower(as.character(tax_rank)) %in% bar_meaningful_ranks) %>%
    group_by(scientific_name, display_group) %>%
    summarise(
      support_priority = min(support_priority, na.rm = TRUE),
      subset_taxon_sum = max(subset_taxon_sum, na.rm = TRUE),
      subset_taxon_max = max(subset_taxon_max, na.rm = TRUE),
      tax_rank         = dplyr::first(tax_rank),
      .groups = "drop"
    ) %>%
    mutate(
      tax_rank_order = dplyr::coalesce(
        unname(tax_rank_order_map[tolower(as.character(tax_rank))]), 999L)
    ) %>%
    arrange(support_priority, desc(subset_taxon_sum), desc(subset_taxon_max), scientific_name) %>%
    slice_head(n = bar_top_n)

  bar_dat <- dat %>%
    semi_join(taxa, by = c("scientific_name", "display_group")) %>%
    group_by(plot_sample_id, scientific_name, display_group) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

  if (nrow(bar_dat) == 0) return(invisible(NULL))

  total_by_lib <- bar_dat %>% group_by(plot_sample_id) %>%
    summarise(total = sum(count), .groups = "drop")
  bar_dat <- bar_dat %>% left_join(total_by_lib, by = "plot_sample_id") %>%
    mutate(prop = ifelse(total > 0, count / total, 0))

  if (!"is_environmental_control" %in% names(dat)) dat$is_environmental_control <- FALSE
  x_label_map <- resolve_x_order(dat, subset_name)
  lib_order   <- names(x_label_map)
  x_labels    <- function(x) x_label_map[x]

  # Order groups by total abundance so dominant groups appear first / at bottom.
  grp_totals <- taxa %>%
    group_by(display_group) %>%
    summarise(grp_total = sum(subset_taxon_sum), .groups = "drop") %>%
    arrange(desc(grp_total))
  grp_levels <- grp_totals$display_group

  # Evenly-spaced hues, one per group.
  hue_vals <- seq(15, 375, length.out = length(grp_levels) + 1)[seq_along(grp_levels)]
  names(hue_vals) <- grp_levels

  # Build fill map and legend entry order: transparent group-name headers then
  # taxa within each group (highest rank = darkest colour), with spacers between groups.
  fill_map       <- c()
  header_labels  <- c()
  spacer_labels  <- c()
  legend_entries <- c()   # headers + taxa + spacers (for legend factor levels)
  taxa_order     <- c()   # taxa only (for bar_dat stacking factor levels)

  for (i in seq_along(grp_levels)) {
    grp <- grp_levels[i]
    # Sort within group: highest taxonomic rank first (family darkest), then by count.
    grp_taxa <- taxa %>% filter(display_group == grp) %>%
      arrange(tax_rank_order, desc(subset_taxon_sum), scientific_name) %>% pull(scientific_name)
    n <- length(grp_taxa)
    if (n == 0) next
    # l = 35 (darkest) for highest-rank entry, up to 70 (lightest) for lowest-rank.
    cols <- if (n == 1) grDevices::hcl(h = hue_vals[[grp]], c = 65, l = 55) else
            grDevices::hcl(h = hue_vals[[grp]], c = 65, l = seq(35, 70, length.out = n))
    names(cols) <- grp_taxa
    fill_map   <- c(fill_map, cols)
    taxa_order <- c(taxa_order, grp_taxa)
    # Always add a bold group header when the group has >1 taxon.
    # When the group name IS itself one of the plotted taxa (e.g. "Bovidae"
    # classified at family rank), the header uses a bold HTML label that maps
    # to the same fill key as the transparent header entries.
    if (length(grp_taxa) == 1 && grp %in% grp_taxa) {
      # Single-taxon group whose name is the taxon itself: no extra header needed.
      legend_entries <- c(legend_entries, grp_taxa)
    } else {
      hdr_key <- paste0("__hdr__", grp)   # unique key not in taxa names
      header_labels  <- c(header_labels,  stats::setNames(grp, hdr_key))
      legend_entries <- c(legend_entries, hdr_key, grp_taxa)
    }
    # Add a small blank spacer between groups (not after the last group).
    if (i < length(grp_levels)) {
      sp <- strrep("\u00a0", i)   # unique run of non-breaking spaces per group
      spacer_labels  <- c(spacer_labels,  sp)
      legend_entries <- c(legend_entries, sp)
    }
  }

  # Transparent fill for group headers and spacers: only text shows in the legend.
  # header_labels is a named vector: name = __hdr__GroupName, value = display label.
  # Build bold HTML labels for headers (rendered via ggtext::element_markdown).
  hdr_bold_labels <- if (length(header_labels) > 0)
    stats::setNames(paste0("<b>", header_labels, "</b>"), names(header_labels))
  else
    character(0)

  fill_map_full <- c(
    fill_map,
    stats::setNames(rep("transparent", length(header_labels)), names(header_labels)),
    stats::setNames(rep("transparent", length(spacer_labels)),  spacer_labels)
  )

  # Legend label map: replace __hdr__X keys with bold display names.
  legend_label_map <- stats::setNames(legend_entries, legend_entries)
  for (k in names(hdr_bold_labels)) legend_label_map[[k]] <- hdr_bold_labels[[k]]
  for (sp in spacer_labels) legend_label_map[[sp]] <- sp

  bar_dat$scientific_name <- factor(bar_dat$scientific_name, levels = taxa_order)
  bar_dat$plot_sample_id  <- factor(bar_dat$plot_sample_id,  levels = lib_order)

  bar_theme <- theme_minimal(base_size = 11) +
    theme(
      axis.text.x           = element_text(angle = 45, hjust = 1, size = 7),
      legend.position       = "right",
      legend.background     = element_blank(),
      legend.box.background = element_rect(colour = "grey40", fill = "white", linewidth = 0.5),
      legend.box.margin     = margin(4, 4, 4, 4),
      legend.key            = element_blank(),
      legend.text           = ggtext::element_markdown(size = 8),
      plot.title            = element_text(face = "plain", size = 10)
    )

  p_count <- ggplot(bar_dat, aes(x = plot_sample_id, y = count, fill = scientific_name)) +
    geom_col(width = 0.88, colour = "grey25", linewidth = 0.15) +
    scale_fill_manual(values = fill_map_full, breaks = legend_entries,
                      limits = legend_entries, labels = legend_label_map,
                      drop = FALSE, na.value = "#aaaaaa", name = "Taxon") +
    scale_x_discrete(labels = x_labels, drop = FALSE) +
    labs(title = paste0(subset_name, " — MEGAN counts"),
         x = NULL, y = "MEGAN count of unique reads") +
    bar_theme + theme(legend.position = "none")

  p_prop <- ggplot(bar_dat, aes(x = plot_sample_id, y = prop, fill = scientific_name)) +
    geom_col(width = 0.88, colour = "grey25", linewidth = 0.15) +
    scale_fill_manual(values = fill_map_full, breaks = legend_entries,
                      limits = legend_entries, labels = legend_label_map,
                      drop = FALSE, na.value = "#aaaaaa", name = "Taxon") +
    scale_x_discrete(labels = x_labels, drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = paste0(subset_name, " — proportions"),
         x = NULL, y = "Proportion of unique reads assigned") +
    bar_theme

  # Collect legend from p_prop (top panel) so it sits centred against the full figure.
  combined <- p_prop / p_count +
    patchwork::plot_layout(heights = c(1, 1), guides = "collect")
  width <- max(9.5, length(lib_order) * 0.55 + 4)
  save_plot_safe(combined, out_file, width = width, height = 8.8)
}

# ---------------------------------------------------------------------------
# make_damage_heatmap
# Companion to make_heatmap: tiles are coloured by Holi damage (C→T rate)
# rather than MEGAN read count.  Only rendered when plot_damage data exist.
# ---------------------------------------------------------------------------
damage_breaks  <- c(-Inf, 0.001, 0.02, 0.05, 0.10, 0.20, Inf)
damage_labels  <- c("No data", "<2%", "2-5%", "5-10%", "10-20%", ">20%")
damage_palette <- c(
  "No data" = "#eeeeee",
  "<2%"     = "#fff7bc",
  "2-5%"    = "#fec44f",
  "5-10%"   = "#fe9929",
  "10-20%"  = "#cc4c02",
  ">20%"    = "#662506"
)
damage_dark_bins <- c("10-20%", ">20%")

bin_damage <- function(x) {
  factor(cut(as.numeric(x), breaks = damage_breaks, labels = damage_labels,
             include.lowest = TRUE, right = FALSE), levels = damage_labels)
}

make_damage_heatmap <- function(dat, subset_name, out_file, top_n = 0, page_rows = 28) {
  if (!"plot_damage" %in% names(dat)) return(invisible(NULL))
  dat <- normalize_support_col(dat)
  dat <- with_inferred_group(dat) %>%
    mutate(
      plot_damage         = as.numeric(plot_damage),
      aDNA_support_status = factor(aDNA_support_status, levels = names(status_priority_map)),
      common_name         = pretty_common_name(common_name),
      count               = as.numeric(count),
      taxon_key           = make_taxon_key(scientific_name, tax_rank, display_group)
    )

  if (identical(subset_name, "__ALL__")) {
    dat <- dat %>% mutate(plot_group = "all_libraries")
    subset_name <- "all_libraries"
  }

  if (!"is_environmental_control" %in% names(dat)) dat$is_environmental_control <- FALSE

  # Only render if at least one real damage value exists for this group.
  grp_dat <- dat %>% filter(plot_group == subset_name)
  if (nrow(grp_dat) == 0) return(invisible(NULL))
  if (all(is.na(grp_dat$plot_damage))) return(invisible(NULL))

  pages <- prepare_subset(dat, subset_name, top_n = top_n, page_rows = page_rows)
  if (is.null(pages)) return(invisible(NULL))

  x_label_map <- resolve_x_order(grp_dat, subset_name)
  x_levels    <- names(x_label_map)
  x_labels    <- function(x) x_label_map[x]

  if (length(x_levels) == 0) return(invisible(NULL))

  max_rows <- max(vapply(pages, length, numeric(1)))

  plots <- list()
  idx   <- 1
  for (page_taxa in pages) {
    raw_page_dat <- dat %>%
      filter(plot_group == subset_name, taxon_key %in% page_taxa) %>%
      mutate(
        common_name = pretty_common_name(common_name),
        count       = as.numeric(count),
        plot_damage = as.numeric(plot_damage),
        taxon_key   = make_taxon_key(scientific_name, tax_rank, display_group)
      )

    label_df <- raw_page_dat %>%
      group_by(taxon_key) %>%
      summarise(
        scientific_name               = dplyr::first(scientific_name),
        display_group                 = dplyr::first(display_group),
        common_name                   = dplyr::first(common_name),
        display_group_status_priority = suppressWarnings(min(display_group_status_priority, na.rm = TRUE)),
        display_group_sum             = suppressWarnings(max(display_group_sum,             na.rm = TRUE)),
        display_group_max             = suppressWarnings(max(display_group_max,             na.rm = TRUE)),
        support_priority              = suppressWarnings(min(support_priority,              na.rm = TRUE)),
        subset_taxon_sum              = suppressWarnings(max(subset_taxon_sum,              na.rm = TRUE)),
        subset_taxon_max              = suppressWarnings(max(subset_taxon_max,              na.rm = TRUE)),
        tax_rank_order                = suppressWarnings(min(tax_rank_order,                na.rm = TRUE)),
        tax_rank                      = dplyr::first(tax_rank),
        .groups = "drop"
      ) %>%
      mutate(
        display_group_status_priority = ifelse(is.infinite(display_group_status_priority), 999, display_group_status_priority),
        display_group_sum  = ifelse(is.infinite(display_group_sum),  0, display_group_sum),
        display_group_max  = ifelse(is.infinite(display_group_max),  0, display_group_max),
        support_priority   = ifelse(is.infinite(support_priority),   999, support_priority),
        subset_taxon_sum   = ifelse(is.infinite(subset_taxon_sum),   0, subset_taxon_sum),
        subset_taxon_max   = ifelse(is.infinite(subset_taxon_max),   0, subset_taxon_max),
        tax_rank_order     = ifelse(is.infinite(tax_rank_order),     999, tax_rank_order)
      ) %>%
      arrange(display_group_status_priority, desc(display_group_sum), desc(display_group_max),
              display_group, support_priority, desc(subset_taxon_sum), desc(subset_taxon_max),
              tax_rank_order, scientific_name) %>%
      mutate(
        row_id               = row_number(),
        display_group_label  = ifelse(duplicated(display_group), "", display_group)
      )

    label_df <- add_plausibility(label_df)
    y_levels <- label_df$taxon_key

    page_dat <- tidyr::expand_grid(plot_sample_id = x_levels, taxon_key = y_levels) %>%
      left_join(raw_page_dat, by = c("plot_sample_id", "taxon_key"), relationship = "many-to-many") %>%
      group_by(plot_sample_id, taxon_key) %>%
      summarise(
        scientific_name  = dplyr::first(stats::na.omit(scientific_name)),
        display_group    = dplyr::first(stats::na.omit(display_group)),
        common_name      = dplyr::first(stats::na.omit(common_name)),
        tax_rank         = dplyr::first(stats::na.omit(tax_rank)),
        count            = sum(count, na.rm = TRUE),
        plot_damage      = suppressWarnings(max(plot_damage, na.rm = TRUE)),
        raw_library_adna_support = dplyr::first(stats::na.omit(as.character(library_adna_support))),
        .groups = "drop"
      ) %>%
      mutate(
        count       = replace(count, is.na(count) | (count == 0), 0),
        plot_damage = ifelse(is.infinite(plot_damage), NA_real_, plot_damage),
        # Tiles with count>0 but no damage data show "No data"; count=0 tiles also "No data"
        dmg_for_bin = ifelse(count > 0 & !is.na(plot_damage), plot_damage, NA_real_),
        damage_bin  = factor(ifelse(count > 0,
                               as.character(bin_damage(dmg_for_bin)),
                               "No data"),
                             levels = damage_labels),
        damage_text = ifelse(count > 0 & !is.na(plot_damage) & plot_damage > 0.001,
                             paste0(round(plot_damage * 100, 1), "%"), ""),
        damage_text_color = ifelse(as.character(damage_bin) %in% damage_dark_bins, "white", "black")
      ) %>%
      left_join(label_df %>% select(taxon_key, row_id), by = "taxon_key") %>%
      mutate(
        taxon_key      = factor(taxon_key,      levels = rev(y_levels)),
        plot_sample_id = factor(plot_sample_id, levels = x_levels),
        damage_bin     = factor(damage_bin,     levels = damage_labels)
      )

    # Dummy data — forces all damage bins into legend on every page.
    fill_dummy_dmg <- tibble::tibble(
      plot_sample_id = factor(rep(x_levels[1], length(damage_labels)), levels = x_levels),
      taxon_key      = factor(rep(rev(y_levels)[1], length(damage_labels)), levels = rev(y_levels)),
      damage_bin     = factor(damage_labels, levels = damage_labels)
    )

    n_rows     <- nrow(label_df)
    boundaries <- label_df %>%
      group_by(display_group) %>%
      summarise(group_end = max(row_id), .groups = "drop") %>%
      filter(group_end < n_rows) %>%
      mutate(yint = n_rows - group_end + 0.5)

    group_width  <- panel_width(label_df$display_group_label, min_width = 0.35, max_width = 1.8, char_scale = 0.055)
    sci_width    <- panel_width(label_df$scientific_name,      min_width = 0.65, max_width = 2.8, char_scale = 0.060)
    common_width <- panel_width(label_df$common_name,          min_width = 0.40, max_width = 2.4, char_scale = 0.053)
    n_libs       <- length(x_levels)
    heat_width   <- max(0.65, n_libs * 0.52)

    n_page_rows  <- length(page_taxa)
    extra_expand <- max(0, (max_rows - n_page_rows)) / 2
    shared_y    <- scale_y_discrete(limits = rev(y_levels), expand = expansion(add = 0.4 + extra_expand))
    title_size    <- 8.0
    subtitle_size <- 7.0

    label_theme <- theme(
      axis.text.y      = element_blank(),
      axis.ticks       = element_blank(),
      axis.title       = element_blank(),
      axis.text.x      = element_text(colour = NA, size = 7, angle = 45, hjust = 1, vjust = 1),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid       = element_blank(),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.title       = element_text(size = title_size,    face = "plain", hjust = 0,
                                      margin = margin(b = 1)),
      plot.subtitle    = element_text(size = subtitle_size, face = "bold",  hjust = 0,
                                      margin = margin(b = 2)),
      plot.margin      = margin(5.5, 2, 5.5, 2)
    )

    group_panel <- ggplot(label_df, aes(x = 0, y = taxon_key)) +
      geom_hline(data = boundaries, aes(yintercept = yint),
                 inherit.aes = FALSE, linewidth = 0.3, colour = "grey30") +
      geom_text(aes(label = display_group_label), hjust = 0, size = 2.7) +
      shared_y +
      coord_cartesian(xlim = c(0, 1), clip = "off") +
      labs(title = " ", subtitle = "Family / Group") +
      label_theme

    sci_panel <- ggplot(label_df, aes(x = 0, y = taxon_key)) +
      geom_hline(data = boundaries, aes(yintercept = yint),
                 inherit.aes = FALSE, linewidth = 0.3, colour = "grey30") +
      geom_text(aes(label = scientific_name), hjust = 0, size = 2.65, fontface = "italic") +
      shared_y +
      coord_cartesian(xlim = c(0, 1), clip = "off") +
      labs(title = " ", subtitle = "Scientific Name") +
      label_theme

    common_panel <- ggplot(label_df, aes(x = 0, y = taxon_key)) +
      geom_hline(data = boundaries, aes(yintercept = yint),
                 inherit.aes = FALSE, linewidth = 0.3, colour = "grey30") +
      geom_text(aes(label = common_name), hjust = 0, size = 2.55) +
      shared_y +
      coord_cartesian(xlim = c(0, 1), clip = "off") +
      labs(title = " ", subtitle = "Common Name") +
      label_theme

    plaus_panel <- make_plaus_panel(label_df, boundaries, shared_y, label_theme)

    heat_panel <- ggplot(page_dat, aes(x = plot_sample_id, y = taxon_key)) +
      geom_hline(data = boundaries, aes(yintercept = yint),
                 inherit.aes = FALSE, linewidth = 0.3, colour = "grey30") +
      geom_tile(data = fill_dummy_dmg,
                aes(x = plot_sample_id, y = taxon_key, fill = damage_bin),
                alpha = 0, width = 0, height = 0, colour = NA,
                inherit.aes = FALSE, show.legend = TRUE, na.rm = TRUE) +
      geom_tile(aes(fill = damage_bin),
                colour = "grey84", linewidth = 0.18,
                width = 0.96, height = 0.88, na.rm = TRUE) +
      geom_text(aes(x = plot_sample_id, y = taxon_key,
                    label = damage_text, colour = damage_text_color),
                size = 2.1, fontface = "plain", na.rm = TRUE,
                inherit.aes = FALSE, show.legend = FALSE) +
      scale_fill_manual(
        values   = damage_palette,
        breaks   = damage_labels,
        limits   = damage_labels,
        labels   = damage_labels,
        drop     = FALSE,
        na.value = damage_palette[["No data"]],
        name     = "Holi Damage (C\u2192T rate)"
      ) +
      scale_colour_identity() +
      guides(fill = guide_legend(order = 1, override.aes = list(alpha = 1, size = 4))) +
      shared_y +
      scale_x_discrete(labels = x_labels, drop = FALSE) +
      labs(
        x        = NULL, y = NULL,
        title    = paste0(subset_name, if (length(pages) > 1) paste0(" (page ", idx, "/", length(pages), ")") else ""),
        subtitle = "Library"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid      = element_blank(),
        axis.text.y     = element_blank(),
        axis.title      = element_blank(),
        axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
        legend.position = "right",
        legend.background     = element_blank(),
        legend.box.background = element_rect(colour = "grey40", fill = "white", linewidth = 0.5),
        legend.box.margin     = margin(4, 4, 4, 4),
        plot.title    = element_text(face = "plain", size = title_size,    margin = margin(b = 1)),
        plot.subtitle = element_text(face = "bold",  size = subtitle_size, colour = "grey35",
                                     margin = margin(b = 2)),
        plot.margin   = margin(5.5, 5.5, 5.5, 2)
      )

    if (!is.null(plaus_panel)) {
      plots[[idx]] <- group_panel + sci_panel + plaus_panel + common_panel + heat_panel +
        patchwork::plot_layout(widths  = c(group_width, sci_width, 0.12, common_width, heat_width),
                               guides  = "collect")
    } else {
      plots[[idx]] <- group_panel + sci_panel + common_panel + heat_panel +
        patchwork::plot_layout(widths = c(group_width, sci_width, common_width, heat_width))
    }
    idx <- idx + 1
  }

  plaus_w <- if (!is.null(plausibility_tbl)) 0.12 else 0.0
  width  <- max(5.5, group_width + sci_width + common_width + plaus_w + heat_width + 2.8)
  height <- max(4.0, 1.4 + max_rows * 0.16)
  save_plot_pages(plots, out_file, width = width, height = height)
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
files <- list.files(input_dir, pattern = "_heatmap_input\\.tsv$", full.names = TRUE)
if (length(files) == 0)
  stop("No *_heatmap_input.tsv files found in input directory: ", input_dir)

run_lines <- c()
for (f in files) {
  dat   <- readr::read_tsv(f, show_col_types = FALSE)
  dat   <- normalize_support_col(dat)
  broad <- sub("_heatmap_input\\.tsv$", "", basename(f))
  message("Rendering ", basename(f))
  run_lines <- c(run_lines, paste("Rendering", basename(f)))

  # When a sample order file is supplied, reassign plot_group from its column
  # headers and exclude any samples not listed in the file.
  if (!is.null(order_map)) {
    dat <- dat %>%
      dplyr::left_join(
        order_map %>%
          dplyr::select(plot_sample_id, plot_group_new) %>%
          dplyr::distinct(),
        by = "plot_sample_id",
        relationship = "many-to-many"
      ) %>%
      dplyr::mutate(
        plot_group = dplyr::if_else(!is.na(plot_group_new), plot_group_new, NA_character_)
      ) %>%
      dplyr::filter(!is.na(plot_group)) %>%
      dplyr::select(-plot_group_new)
  }

  make_heatmap(dat, "__ALL__",
               file.path(out_dir, paste0(broad, "_all_libraries_heatmap.pdf")),
               top_n = top_n, page_rows = page_rows)

  for (grp in unique(dat$plot_group)) {
    grp_lab <- subset_label(grp)
    make_heatmap(dat, grp,
                 file.path(out_dir, paste0(broad, "_", grp_lab, "_heatmap.pdf")),
                 top_n = top_n, page_rows = page_rows)
    make_damage_heatmap(dat, grp,
                        file.path(out_dir, paste0(broad, "_", grp_lab, "_damage_heatmap.pdf")),
                        top_n = top_n, page_rows = page_rows)
    make_stacked_bars(dat, grp,
                      file.path(out_dir, paste0(broad, "_", grp_lab, "_stacked_bars.pdf")),
                      bar_top_n = bar_top_n)
  }
}
writeLines(c(run_lines, paste("Reports written to:", out_dir)),
           con = file.path(out_dir, "run_report.txt"))
message("Reports written to: ", out_dir)
