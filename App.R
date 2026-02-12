# =============================================================================
# SeqDedupe v0.1.0
# Lightweight sequence deduplication & preprocessing for wet-lab biologists
#
# Supports: FASTA / FASTQ
# Pipeline: Filter (length, GC, ambiguous bases) → Deduplicate (hash-based)
# Outputs:  Clean sequences, audit trail, cluster report, run report
#
# No Bioconductor dependency. Pure base R + digest + ggplot2.
# =============================================================================

APP_VERSION <- "0.1.0"

# --- Dependencies ------------------------------------------------------------
for (pkg in c("shiny", "DT", "digest", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(shiny)
library(DT)
library(digest)
library(ggplot2)

# =============================================================================
# CORE ENGINE
# =============================================================================

parse_fasta <- function(filepath, source_name = "file") {
  lines <- readLines(filepath, warn = FALSE)
  lines <- lines[nchar(trimws(lines)) > 0]
  
  header_idx <- which(startsWith(lines, ">"))
  if (length(header_idx) == 0) return(NULL)
  
  n <- length(header_idx)
  headers <- sub("^>", "", lines[header_idx])
  sequences <- character(n)
  
  for (i in seq_len(n)) {
    start <- header_idx[i] + 1L
    end <- if (i < n) header_idx[i + 1L] - 1L else length(lines)
    if (start <= end) {
      sequences[i] <- paste0(lines[start:end], collapse = "")
    }
  }
  
  data.frame(header = headers, sequence = sequences,
             file_source = source_name, stringsAsFactors = FALSE)
}

parse_fastq <- function(filepath, source_name = "file") {
  lines <- readLines(filepath, warn = FALSE)
  n_lines <- length(lines)
  if (n_lines < 4) return(NULL)
  
  idx <- seq(1, n_lines, by = 4)
  valid <- startsWith(lines[idx], "@")
  idx <- idx[valid]
  
  data.frame(
    header = sub("^@", "", lines[idx]),
    sequence = lines[idx + 1L],
    quality = lines[idx + 3L],
    file_source = source_name,
    stringsAsFactors = FALSE
  )
}

detect_format <- function(filename) {
  ext <- tolower(tools::file_ext(filename))
  if (ext %in% c("fastq", "fq")) return("fastq")
  return("fasta")
}

#' Auto-detect whether sequences are nucleotide or protein
#' Samples up to 100 sequences; if >5% of non-gap characters are
#' amino-acid-only (not ATCGU), classify as protein.
detect_seq_type <- function(sequences, n_sample = 100) {
  sample_seqs <- toupper(utils::head(sequences, n_sample))
  combined <- paste(sample_seqs, collapse = "")
  clean <- gsub("[\\-\\.\\*NnXx]", "", combined)
  total <- nchar(clean)
  if (total == 0) return("nucleotide")
  protein_only <- nchar(gsub("[ATCGU]", "", clean))
  if (protein_only / total > 0.05) return("protein")
  return("nucleotide")
}

compute_metrics <- function(df, seq_type = "nucleotide") {
  df$seq_length <- nchar(df$sequence)
  
  if (seq_type == "nucleotide") {
    gc_count <- nchar(gsub("[^GCgc]", "", df$sequence))
    df$gc_content <- ifelse(df$seq_length > 0,
                            round(gc_count / df$seq_length * 100, 2), 0)
    # N = unknown nucleotide
    df$n_count <- nchar(gsub("[^Nn]", "", df$sequence))
    # Ambiguous = anything not standard ATCGU
    ambig_count <- nchar(gsub("[ATCGUatcgu]", "", df$sequence))
  } else {
    # Protein: GC content not applicable
    df$gc_content <- NA_real_
    # X = unknown amino acid residue
    df$n_count <- nchar(gsub("[^Xx]", "", df$sequence))
    # Ambiguous = anything not one of the 20 standard amino acids
    ambig_count <- nchar(gsub("[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]", "", df$sequence))
  }
  
  df$n_ambiguous <- ambig_count
  df$pct_ambiguous <- ifelse(df$seq_length > 0,
                             round(ambig_count / df$seq_length * 100, 2), 0)
  df$seq_hash <- vapply(df$sequence, function(s) digest(s, algo = "md5"),
                        character(1), USE.NAMES = FALSE)
  df
}

apply_filters <- function(df, min_len = 0, max_len = Inf,
                          min_gc = 0, max_gc = 100,
                          max_n_count = Inf, max_ambig_pct = 100,
                          seq_type = "nucleotide") {
  keep_len <- df$seq_length >= min_len & df$seq_length <= max_len
  # GC filter only applies to nucleotide sequences
  if (seq_type == "nucleotide") {
    keep_gc <- df$gc_content >= min_gc & df$gc_content <= max_gc
  } else {
    keep_gc <- rep(TRUE, nrow(df))
  }
  keep_n   <- df$n_count <= max_n_count
  keep_amb <- df$pct_ambiguous <= max_ambig_pct
  
  # Tag filter reason (sequential logic: first failing filter wins)
  unknown_label <- if (seq_type == "nucleotide") "N count" else "X count"
  df$filter_reason <- ""
  df$filter_reason[!keep_len] <- "length"
  df$filter_reason[df$filter_reason == "" & !keep_gc] <- "GC content"
  df$filter_reason[df$filter_reason == "" & !keep_n] <- unknown_label
  df$filter_reason[df$filter_reason == "" & !keep_amb] <- "ambiguous residues"
  
  keep <- keep_len & keep_gc & keep_n & keep_amb
  
  list(
    passed = df[keep, ],
    filtered_out = df[!keep, ],
    n_removed_length = sum(!keep_len),
    n_removed_gc = sum(keep_len & !keep_gc),
    n_removed_n = sum(keep_len & keep_gc & !keep_n),
    n_removed_ambig = sum(keep_len & keep_gc & keep_n & !keep_amb),
    n_total_removed = sum(!keep)
  )
}

deduplicate_seqs <- function(df, method = "sequence", keep = "first",
                             case_insensitive = TRUE) {
  n <- nrow(df)
  
  if (case_insensitive) {
    seq_for_key <- toupper(df$sequence)
    hash_key <- vapply(seq_for_key, function(s) digest(s, algo = "md5"),
                       character(1), USE.NAMES = FALSE)
  } else {
    hash_key <- df$seq_hash
  }
  
  dedup_key <- switch(method,
    "sequence" = hash_key,
    "header"   = df$header,
    "both"     = paste(df$header, hash_key, sep = "|||")
  )
  
  key_to_cluster <- match(dedup_key, unique(dedup_key))
  df$cluster_id <- key_to_cluster
  cluster_sizes <- table(key_to_cluster)
  df$cluster_size <- as.integer(cluster_sizes[as.character(key_to_cluster)])
  
  if (keep == "first") {
    is_dup <- duplicated(dedup_key)
  } else {
    is_dup <- duplicated(dedup_key, fromLast = TRUE)
  }
  
  df$status <- ifelse(is_dup, "removed_duplicate", "retained")
  
  deduped <- df[!is_dup, ]
  removed <- df[is_dup, ]
  
  dup_cluster_ids <- unique(df$cluster_id[df$cluster_size > 1])
  cluster_summary <- if (length(dup_cluster_ids) > 0) {
    do.call(rbind, lapply(dup_cluster_ids, function(cid) {
      members <- df[df$cluster_id == cid, ]
      retained <- members[members$status == "retained", ]
      data.frame(
        cluster_id = cid,
        cluster_size = nrow(members),
        retained_header = retained$header[1],
        retained_length = retained$seq_length[1],
        removed_headers = paste(
          members$header[members$status == "removed_duplicate"], collapse = "; "),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(cluster_id = integer(0), cluster_size = integer(0),
               retained_header = character(0), retained_length = integer(0),
               removed_headers = character(0), stringsAsFactors = FALSE)
  }
  
  list(
    all_annotated   = df,
    deduped         = deduped,
    removed         = removed,
    cluster_summary = cluster_summary,
    n_original      = n,
    n_kept          = nrow(deduped),
    n_removed       = nrow(removed),
    n_dup_clusters  = length(dup_cluster_ids)
  )
}

write_fasta_out <- function(df, filepath, line_width = 80) {
  con <- file(filepath, "w")
  on.exit(close(con))
  for (i in seq_len(nrow(df))) {
    writeLines(paste0(">", df$header[i]), con)
    s <- df$sequence[i]
    starts <- seq(1, nchar(s), by = line_width)
    for (st in starts) {
      writeLines(substr(s, st, min(st + line_width - 1, nchar(s))), con)
    }
  }
}

write_fastq_out <- function(df, filepath) {
  con <- file(filepath, "w")
  on.exit(close(con))
  for (i in seq_len(nrow(df))) {
    writeLines(paste0("@", df$header[i]), con)
    writeLines(df$sequence[i], con)
    writeLines("+", con)
    qual <- if ("quality" %in% names(df) && !is.na(df$quality[i])) {
      df$quality[i]
    } else {
      paste(rep("I", nchar(df$sequence[i])), collapse = "")
    }
    writeLines(qual, con)
  }
}

generate_run_report <- function(params, results) {
  divider <- "==========================================================="
  is_nuc <- (params$seq_type == "nucleotide")
  unknown_label <- if (is_nuc) "Max N count" else "Max X count"
  unit <- if (is_nuc) "bp" else "aa"
  
  lines <- c(
    divider,
    sprintf("  SeqDedupe v%s  -  Run Report", APP_VERSION),
    divider,
    sprintf("  Date:       %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("  R version:  %s", R.version.string),
    sprintf("  Platform:   %s", R.version$platform),
    "-----------------------------------------------------------",
    "  INPUT",
    sprintf("    Files:          %s", params$files),
    sprintf("    Format:         %s", toupper(params$format)),
    sprintf("    Sequence type:  %s", if (is_nuc) "Nucleotide" else "Protein"),
    sprintf("    Sequences:      %s", format(results$n_input, big.mark = ",")),
    "-----------------------------------------------------------",
    "  FILTER PARAMETERS",
    sprintf("    Min length:     %s %s", params$min_len, unit),
    sprintf("    Max length:     %s", if (is.na(params$max_len)) "none"
            else paste(format(params$max_len, big.mark = ","), unit))
  )
  
  if (is_nuc) {
    lines <- c(lines,
      sprintf("    GC range:       %s - %s%%", params$min_gc, params$max_gc))
  }
  
  lines <- c(lines,
    sprintf("    %-16s%s", paste0(unknown_label, ":"),
            if (is.na(params$max_n)) "none" else params$max_n),
    sprintf("    Max ambig %%:    %s%%", params$max_ambig),
    "-----------------------------------------------------------",
    "  DEDUPLICATION PARAMETERS",
    sprintf("    Method:         by %s", params$dedup_method),
    sprintf("    Keep:           %s occurrence", params$keep_which),
    sprintf("    Case-insensitive: %s", params$case_insensitive),
    "-----------------------------------------------------------",
    "  PIPELINE RESULTS",
    ""
  )
  
  st <- results$step_table
  for (i in seq_len(nrow(st))) {
    removed_str <- if (is.na(st$Removed[i]) || st$Removed[i] == 0) ""
                   else sprintf("  (-%s)", format(st$Removed[i], big.mark = ","))
    lines <- c(lines, sprintf("    %-28s %8s%s",
                               st$Step[i],
                               format(st$Sequences[i], big.mark = ","),
                               removed_str))
  }
  
  lines <- c(lines,
    "",
    "-----------------------------------------------------------",
    sprintf("  FINAL OUTPUT:     %s sequences", format(results$dedup$n_kept, big.mark = ",")),
    sprintf("  DUPLICATE GROUPS: %s", results$dedup$n_dup_clusters),
    divider,
    "",
    "Citation:",
    sprintf("  SeqDedupe v%s", APP_VERSION),
    "  https://github.com/mbaffour/seqdedupe",
    ""
  )
  
  paste(lines, collapse = "\n")
}

# =============================================================================
# UI
# =============================================================================

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: 'Segoe UI', system-ui, -apple-system, sans-serif; background: #f5f6fa; }
    
    .app-header {
      background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
      color: white; padding: 22px 30px; margin: -15px -15px 20px;
      border-bottom: 3px solid #e94560;
    }
    .app-header h2 { margin: 0 0 4px; font-weight: 700; letter-spacing: -0.5px; }
    .app-header .subtitle { opacity: 0.75; font-size: 13px; }
    .app-header .version-badge { 
      float: right; background: rgba(233,69,96,0.85); 
      padding: 4px 14px; border-radius: 12px; font-size: 12px; 
      font-weight: 600; margin-top: 4px;
    }
    
    .stat-card {
      background: white; border-radius: 10px; padding: 14px 16px;
      text-align: center; box-shadow: 0 2px 8px rgba(0,0,0,0.06);
      border-left: 4px solid #ccc; margin-bottom: 12px;
    }
    .stat-card .val { font-size: 24px; font-weight: 700; }
    .stat-card .lab { 
      font-size: 10px; color: #7f8c8d; text-transform: uppercase; 
      letter-spacing: 0.5px; margin-top: 2px;
    }
    .s-input  { border-left-color: #3498db; } .s-input .val  { color: #2c3e50; }
    .s-filt   { border-left-color: #f39c12; } .s-filt .val   { color: #f39c12; }
    .s-dup    { border-left-color: #e74c3c; } .s-dup .val    { color: #e74c3c; }
    .s-kept   { border-left-color: #27ae60; } .s-kept .val   { color: #27ae60; }
    .s-clust  { border-left-color: #9b59b6; } .s-clust .val  { color: #9b59b6; }
    .s-pct    { border-left-color: #e74c3c; } .s-pct .val    { color: #e74c3c; }
    
    .card { 
      background: white; border-radius: 10px; padding: 20px;
      box-shadow: 0 2px 8px rgba(0,0,0,0.06); margin-bottom: 15px;
    }
    .card h4 { 
      margin-top: 0; color: #1a1a2e; 
      border-bottom: 2px solid #ecf0f1; padding-bottom: 10px; 
    }
    
    .step-table { width: 100%; border-collapse: collapse; font-size: 14px; }
    .step-table th { 
      text-align: left; padding: 10px 14px; background: #f8f9fa; 
      border-bottom: 2px solid #dee2e6; color: #495057; font-weight: 600;
    }
    .step-table td { padding: 10px 14px; border-bottom: 1px solid #eee; }
    .step-table tr:last-child td { 
      font-weight: 700; background: #eafaf1; border-bottom: 2px solid #27ae60; 
    }
    .step-delta { color: #e74c3c; font-size: 12px; }
    
    .section-label { 
      font-weight: 600; color: #1a1a2e; font-size: 12px; 
      text-transform: uppercase; letter-spacing: 0.5px;
      margin: 14px 0 6px; border-bottom: 1px solid #eee; padding-bottom: 5px;
    }
    
    .btn-run {
      background: linear-gradient(135deg, #e94560, #c0392b); color: white;
      border: none; font-weight: 700; font-size: 15px; padding: 10px;
      width: 100%; border-radius: 6px; margin-top: 5px; cursor: pointer;
    }
    .btn-run:hover { background: linear-gradient(135deg, #c0392b, #e94560); color: white; }
    
    .btn-dl {
      background: #27ae60; color: white; border: none; font-weight: 600;
      width: 100%; border-radius: 6px; padding: 8px; margin-top: 4px;
      font-size: 12px;
    }
    .btn-dl:hover { background: #219a52; color: white; }
    
    .file-warning { 
      background: #fef9e7; border-left: 4px solid #f39c12; 
      padding: 10px 14px; border-radius: 4px; margin: 8px 0; font-size: 12px;
    }
    
    .repro-stamp {
      background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 6px;
      padding: 14px 18px; font-family: 'Consolas', 'Fira Code', monospace;
      font-size: 12px; color: #495057; white-space: pre-wrap;
      max-height: 500px; overflow-y: auto; line-height: 1.5;
    }
    
    .nav-tabs > li > a { color: #1a1a2e; font-weight: 500; }
    hr { border-color: #ecf0f1; margin: 10px 0; }
  "))),
  
  div(class = "app-header",
    span(class = "version-badge", paste0("v", APP_VERSION)),
    h2("SeqDedupe"),
    div(class = "subtitle",
        "Deduplicate, filter & QC your FASTA/FASTQ sequences before downstream analysis")
  ),
  
  sidebarLayout(
    sidebarPanel(width = 3,
      div(class = "card",
        h4("Upload"),
        fileInput("seq_files", "FASTA / FASTQ files",
                  multiple = TRUE,
                  accept = c(".fasta", ".fa", ".fna", ".faa", ".fas",
                             ".fastq", ".fq", ".txt")),
        uiOutput("file_warning_ui"),
        
        div(class = "section-label", "Sequence Type"),
        selectInput("seq_type", NULL,
          c("Auto-detect" = "auto", "Nucleotide (DNA/RNA)" = "nucleotide",
            "Protein" = "protein")),
        
        div(class = "section-label", "Length Filter"),
        fluidRow(
          column(6, numericInput("min_len", "Min", value = 0, min = 0)),
          column(6, numericInput("max_len", "Max", value = NA, min = 1))
        ),
        
        div(class = "section-label", "Composition Filter"),
        uiOutput("composition_filter_ui"),
        
        div(class = "section-label", "Deduplication"),
        selectInput("dedup_method", NULL,
          c("By sequence" = "sequence",
            "By header" = "header",
            "By both" = "both")
        ),
        fluidRow(
          column(6, radioButtons("keep_which", NULL,
            c("Keep first" = "first", "Keep last" = "last"), inline = FALSE)),
          column(6, checkboxInput("case_insensitive", "Case-insensitive", TRUE))
        ),
        
        div(class = "section-label", "Output"),
        fluidRow(
          column(6, selectInput("output_format", NULL,
            c("Auto" = "auto", "FASTA" = "fasta", "FASTQ" = "fastq"))),
          column(6, numericInput("line_width", "Width", 80, min = 10, max = 200))
        ),
        
        actionButton("run_btn", "Run Pipeline", class = "btn-run",
                     icon = icon("play")),
        br(), br(),
        uiOutput("download_buttons_ui")
      )
    ),
    
    mainPanel(width = 9,
      uiOutput("stats_ui"),
      
      tabsetPanel(id = "tabs", type = "tabs",
        tabPanel("Pipeline Summary",
          div(class = "card", style = "margin-top: 15px;",
            h4("Stepwise Processing Summary"),
            uiOutput("step_table_ui")
          ),
          fluidRow(
            column(6,
              div(class = "card",
                h4("Input Summary"),
                verbatimTextOutput("input_summary")
              )
            ),
            column(6,
              div(class = "card",
                h4("Length Distribution"),
                plotOutput("length_hist", height = "280px")
              )
            )
          )
        ),
        tabPanel("Duplicate Clusters",
          div(class = "card", style = "margin-top: 15px;",
            h4("Duplicate Cluster Detail"),
            p("Each row groups identical sequences. ",
              tags$b("Retained Header"), " = kept. ",
              tags$b("Removed Headers"), " = dropped."),
            DTOutput("cluster_table")
          )
        ),
        tabPanel("Filtered Out",
          div(class = "card", style = "margin-top: 15px;",
            h4("Sequences Removed by Filters"),
            DTOutput("filter_table")
          )
        ),
        tabPanel("Full Audit",
          div(class = "card", style = "margin-top: 15px;",
            h4("Complete Audit Trail"),
            p("Every input sequence tagged with its final disposition."),
            DTOutput("audit_table")
          )
        ),
        tabPanel("Run Report",
          div(class = "card", style = "margin-top: 15px;",
            h4("Reproducibility Report"),
            p("Parameters, versions, and results. Include in supplementary materials or lab notebooks."),
            div(class = "repro-stamp", textOutput("run_report_text"))
          )
        ),
        tabPanel("Preview",
          div(class = "card", style = "margin-top: 15px;",
            h4("Output Preview (first 100)"),
            verbatimTextOutput("output_preview")
          )
        )
      )
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {
  
  raw_data      <- reactiveVal(NULL)
  input_format  <- reactiveVal("fasta")
  pipeline_out  <- reactiveVal(NULL)
  run_params    <- reactiveVal(NULL)
  
  # Resolved sequence type: auto-detect or manual override
  seq_type_resolved <- reactive({
    if (input$seq_type != "auto") return(input$seq_type)
    rd <- raw_data()
    if (is.null(rd)) return("nucleotide")
    detect_seq_type(rd$sequence)
  })
  
  # --- Dynamic composition filter UI ----------------------------------------
  output$composition_filter_ui <- renderUI({
    st <- seq_type_resolved()
    is_nuc <- (st == "nucleotide")
    unknown_label <- if (is_nuc) "Max N count" else "Max X count"
    
    tagList(
      if (is_nuc) {
        fluidRow(
          column(6, numericInput("min_gc", "Min GC%", value = 0, min = 0, max = 100)),
          column(6, numericInput("max_gc", "Max GC%", value = 100, min = 0, max = 100))
        )
      },
      fluidRow(
        column(6, numericInput("max_n_count", unknown_label, value = NA, min = 0)),
        column(6, numericInput("max_ambig", "Max ambig %", value = 100, min = 0, max = 100))
      )
    )
  })
  
  # --- Upload ----------------------------------------------------------------
  observeEvent(input$seq_files, {
    req(input$seq_files)
    pipeline_out(NULL)
    
    all_dfs <- list()
    detected_formats <- character()
    
    withProgress(message = "Reading files...", value = 0, {
      n_files <- nrow(input$seq_files)
      for (i in seq_len(n_files)) {
        fp <- input$seq_files$datapath[i]
        fname <- input$seq_files$name[i]
        fmt <- detect_format(fname)
        detected_formats <- c(detected_formats, fmt)
        incProgress(1 / n_files, detail = fname)
        
        tryCatch({
          df <- if (fmt == "fastq") parse_fastq(fp, fname) else parse_fasta(fp, fname)
          if (!is.null(df) && nrow(df) > 0) {
            all_dfs[[fname]] <- df
          } else {
            showNotification(paste("No sequences in", fname), type = "warning")
          }
        }, error = function(e) {
          showNotification(paste("Error reading", fname, ":", e$message),
                           type = "error", duration = 8)
        })
      }
    })
    
    if (length(all_dfs) > 0) {
      has_quality <- all(sapply(all_dfs, function(d) "quality" %in% names(d)))
      if (!has_quality) {
        all_dfs <- lapply(all_dfs, function(d) {
          if (!"quality" %in% names(d)) d$quality <- NA_character_
          d
        })
      }
      combined <- do.call(rbind, all_dfs)
      rownames(combined) <- NULL
      combined$original_index <- seq_len(nrow(combined))
      
      fmt_table <- table(detected_formats)
      input_format(names(which.max(fmt_table)))
      raw_data(combined)
      
      showNotification(
        sprintf("Loaded %s sequences from %d file(s)",
                format(nrow(combined), big.mark = ","), length(all_dfs)),
        type = "message", duration = 4)
    }
  })
  
  # --- File size warning -----------------------------------------------------
  output$file_warning_ui <- renderUI({
    req(input$seq_files)
    total_mb <- sum(input$seq_files$size) / 1e6
    if (total_mb > 50) {
      div(class = "file-warning", icon("triangle-exclamation"),
          sprintf(" Large upload (%.0f MB). For >500 MB consider seqkit rmdup.", total_mb))
    } else if (total_mb > 10) {
      div(class = "file-warning", icon("circle-info"),
          sprintf(" %.0f MB uploaded. Should process fine.", total_mb))
    }
  })
  
  # --- Input summary ---------------------------------------------------------
  output$input_summary <- renderText({
    rd <- raw_data()
    if (is.null(rd)) return("Upload files to begin.")
    st <- seq_type_resolved()
    unit <- if (st == "nucleotide") "bp" else "aa"
    sl <- nchar(rd$sequence)
    sources <- table(rd$file_source)
    file_lines <- paste(sprintf("  %s : %s seqs", names(sources),
                                format(sources, big.mark = ",")), collapse = "\n")
    paste(sep = "\n",
      sprintf("Format: %s", toupper(input_format())),
      sprintf("Type:   %s", if (st == "nucleotide") "Nucleotide" else "Protein"),
      sprintf("Files:  %d", length(sources)),
      file_lines, "",
      sprintf("Total: %s sequences", format(nrow(rd), big.mark = ",")),
      sprintf("Unique headers: %s", format(length(unique(rd$header)), big.mark = ",")),
      "",
      sprintf("Length: %s - %s %s (median %s)",
              format(min(sl), big.mark = ","),
              format(max(sl), big.mark = ","),
              unit,
              format(median(sl), big.mark = ",")))
  })
  
  # --- Run pipeline ----------------------------------------------------------
  observeEvent(input$run_btn, {
    rd <- raw_data()
    req(rd)
    
    params <- list(
      files = paste(input$seq_files$name, collapse = ", "),
      format = input_format(),
      seq_type = seq_type_resolved(),
      min_len = input$min_len,
      max_len = input$max_len,
      min_gc = if (!is.null(input$min_gc)) input$min_gc else 0,
      max_gc = if (!is.null(input$max_gc)) input$max_gc else 100,
      max_n = input$max_n_count,
      max_ambig = input$max_ambig,
      dedup_method = input$dedup_method,
      keep_which = input$keep_which,
      case_insensitive = input$case_insensitive
    )
    run_params(params)
    
    results <- list()
    st <- seq_type_resolved()
    is_nuc <- (st == "nucleotide")
    
    withProgress(message = "Running pipeline...", value = 0, {
      
      incProgress(0.1, detail = "Computing sequence metrics...")
      df <- compute_metrics(rd, seq_type = st)
      results$n_input <- nrow(df)
      results$input_lengths <- df$seq_length
      
      incProgress(0.25, detail = "Applying filters...")
      max_len <- if (is.na(input$max_len)) Inf else input$max_len
      max_n   <- if (is.na(input$max_n_count)) Inf else input$max_n_count
      min_gc  <- if (is_nuc && !is.null(input$min_gc)) input$min_gc else 0
      max_gc  <- if (is_nuc && !is.null(input$max_gc)) input$max_gc else 100
      
      filt <- apply_filters(df,
        min_len = input$min_len, max_len = max_len,
        min_gc = min_gc, max_gc = max_gc,
        max_n_count = max_n, max_ambig_pct = input$max_ambig,
        seq_type = st)
      
      results$filter <- filt
      n_after_filter <- nrow(filt$passed)
      
      if (n_after_filter == 0) {
        showNotification("All sequences filtered out! Adjust settings.",
                         type = "error")
        pipeline_out(NULL)
        return()
      }
      
      incProgress(0.5, detail = sprintf("Deduplicating %s sequences...",
                                        format(n_after_filter, big.mark = ",")))
      dedup <- deduplicate_seqs(filt$passed,
        method = input$dedup_method, keep = input$keep_which,
        case_insensitive = input$case_insensitive)
      
      results$dedup <- dedup
      results$output_lengths <- dedup$deduped$seq_length
      
      incProgress(0.85, detail = "Building reports...")
      
      # Stepwise summary table
      n_after_len <- results$n_input - filt$n_removed_length
      n_after_gc  <- n_after_len - filt$n_removed_gc
      n_after_n   <- n_after_gc - filt$n_removed_n - filt$n_removed_ambig
      
      unknown_label <- if (is_nuc) "N" else "X"
      
      if (is_nuc) {
        results$step_table <- data.frame(
          Step = c("Input",
                   "After length filter",
                   "After GC filter",
                   sprintf("After %s/ambiguous filter", unknown_label),
                   "After deduplication"),
          Sequences = c(results$n_input, n_after_len, n_after_gc,
                        n_after_n, dedup$n_kept),
          Removed = c(NA, filt$n_removed_length, filt$n_removed_gc,
                      filt$n_removed_n + filt$n_removed_ambig, dedup$n_removed),
          stringsAsFactors = FALSE
        )
      } else {
        results$step_table <- data.frame(
          Step = c("Input",
                   "After length filter",
                   sprintf("After %s/ambiguous filter", unknown_label),
                   "After deduplication"),
          Sequences = c(results$n_input, n_after_len,
                        n_after_n, dedup$n_kept),
          Removed = c(NA, filt$n_removed_length,
                      filt$n_removed_n + filt$n_removed_ambig, dedup$n_removed),
          stringsAsFactors = FALSE
        )
      }
      
      # Audit trail
      audit_cols <- c("original_index", "header", "file_source",
                      "seq_length", "gc_content", "n_count", "pct_ambiguous")
      
      dedup_audit <- dedup$all_annotated
      filt_out <- filt$filtered_out
      
      if (nrow(filt_out) > 0) {
        filt_out$status <- "filtered"
        filt_out$cluster_id <- NA_integer_
        filt_out$cluster_size <- NA_integer_
      }
      
      for (col in c(audit_cols, "status")) {
        if (!col %in% names(dedup_audit)) dedup_audit[[col]] <- NA
        if (nrow(filt_out) > 0 && !col %in% names(filt_out)) filt_out[[col]] <- NA
      }
      
      results$audit <- rbind(
        dedup_audit[, c(audit_cols, "status"), drop = FALSE],
        if (nrow(filt_out) > 0) filt_out[, c(audit_cols, "status"), drop = FALSE]
      )
      
      incProgress(1, detail = "Done!")
    })
    
    pipeline_out(results)
    
    showNotification(
      sprintf("Done: %s input -> %s output",
              format(results$n_input, big.mark = ","),
              format(results$dedup$n_kept, big.mark = ",")),
      type = "message", duration = 5)
  })
  
  # --- Stats cards -----------------------------------------------------------
  output$stats_ui <- renderUI({
    res <- pipeline_out()
    if (is.null(res)) return(NULL)
    
    n_filt <- res$filter$n_total_removed
    n_dup <- res$dedup$n_removed
    pct <- if (res$n_input > 0) round((n_filt + n_dup) / res$n_input * 100, 1) else 0
    
    fluidRow(
      column(2, div(class = "stat-card s-input",
        div(class = "val", format(res$n_input, big.mark = ",")),
        div(class = "lab", "Input"))),
      column(2, div(class = "stat-card s-filt",
        div(class = "val", format(n_filt, big.mark = ",")),
        div(class = "lab", "Filtered"))),
      column(2, div(class = "stat-card s-dup",
        div(class = "val", format(n_dup, big.mark = ",")),
        div(class = "lab", "Duplicates"))),
      column(2, div(class = "stat-card s-kept",
        div(class = "val", format(res$dedup$n_kept, big.mark = ",")),
        div(class = "lab", "Output"))),
      column(2, div(class = "stat-card s-clust",
        div(class = "val", format(res$dedup$n_dup_clusters, big.mark = ",")),
        div(class = "lab", "Dup Clusters"))),
      column(2, div(class = "stat-card s-pct",
        div(class = "val", paste0(pct, "%")),
        div(class = "lab", "Removed")))
    )
  })
  
  # --- Stepwise summary table ------------------------------------------------
  output$step_table_ui <- renderUI({
    res <- pipeline_out()
    if (is.null(res)) return(p("Run the pipeline to see results."))
    
    st <- res$step_table
    rows <- lapply(seq_len(nrow(st)), function(i) {
      delta <- if (is.na(st$Removed[i]) || st$Removed[i] == 0) ""
               else sprintf('<span class="step-delta">\u2212%s</span>',
                            format(st$Removed[i], big.mark = ","))
      tags$tr(
        tags$td(st$Step[i]),
        tags$td(style = "text-align:right;", format(st$Sequences[i], big.mark = ",")),
        tags$td(style = "text-align:right;", HTML(delta))
      )
    })
    
    tags$table(class = "step-table",
      tags$thead(tags$tr(
        tags$th("Pipeline Step"),
        tags$th(style = "text-align:right;", "Sequences"),
        tags$th(style = "text-align:right;", "Removed")
      )),
      tags$tbody(rows)
    )
  })
  
  # --- Before/After histogram ------------------------------------------------
  output$length_hist <- renderPlot({
    res <- pipeline_out()
    rd <- raw_data()
    if (is.null(rd)) return(NULL)
    
    st <- seq_type_resolved()
    x_label <- if (st == "nucleotide") "Sequence Length (bp)" else "Sequence Length (aa)"
    
    if (!is.null(res)) {
      df_plot <- rbind(
        data.frame(length = res$input_lengths, stage = "Before"),
        data.frame(length = res$output_lengths, stage = "After")
      )
      df_plot$stage <- factor(df_plot$stage, levels = c("Before", "After"))
      fill_vals <- c("Before" = "#bdc3c7", "After" = "#27ae60")
    } else {
      df_plot <- data.frame(length = nchar(rd$sequence), stage = "Before")
      df_plot$stage <- factor(df_plot$stage, levels = "Before")
      fill_vals <- c("Before" = "#3498db")
    }
    
    n_bins <- min(50, length(unique(df_plot$length)))
    
    ggplot(df_plot, aes(x = length, fill = stage)) +
      geom_histogram(bins = n_bins, alpha = 0.65, position = "identity",
                     color = "white", linewidth = 0.3) +
      scale_fill_manual(values = fill_vals) +
      theme_minimal(base_size = 13) +
      labs(x = x_label, y = "Count", fill = NULL) +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = "#dee2e6",
                                         linewidth = 0.3)
      )
  })
  
  # --- Cluster table ---------------------------------------------------------
  output$cluster_table <- renderDT({
    res <- pipeline_out()
    if (is.null(res) || nrow(res$dedup$cluster_summary) == 0) {
      return(datatable(data.frame(Message = "No duplicate clusters found."),
                       options = list(dom = "t"), rownames = FALSE))
    }
    cs <- res$dedup$cluster_summary
    st <- seq_type_resolved()
    len_unit <- if (st == "nucleotide") "bp" else "aa"
    names(cs) <- c("Cluster", "Size", "Retained Header",
                    paste0("Retained Length (", len_unit, ")"), "Removed Headers")
    datatable(cs, options = list(pageLength = 20, scrollX = TRUE,
                                  order = list(list(1, "desc"))),
              rownames = FALSE, filter = "top")
  })
  
  # --- Filter table ----------------------------------------------------------
  output$filter_table <- renderDT({
    res <- pipeline_out()
    if (is.null(res) || nrow(res$filter$filtered_out) == 0) {
      return(datatable(data.frame(Message = "No sequences filtered out."),
                       options = list(dom = "t"), rownames = FALSE))
    }
    fo <- res$filter$filtered_out[, c("original_index", "header", "file_source",
                                       "seq_length", "gc_content", "n_count",
                                       "pct_ambiguous", "filter_reason")]
    st <- seq_type_resolved()
    unknown_col <- if (st == "nucleotide") "N Count" else "X Count"
    gc_col <- if (st == "nucleotide") "GC%" else "GC% (n/a)"
    names(fo) <- c("Index", "Header", "Source", "Length", gc_col,
                    unknown_col, "Ambig%", "Reason")
    datatable(fo, options = list(pageLength = 20, scrollX = TRUE),
              rownames = FALSE, filter = "top")
  })
  
  # --- Audit table -----------------------------------------------------------
  output$audit_table <- renderDT({
    res <- pipeline_out()
    if (is.null(res)) {
      return(datatable(data.frame(Message = "Run the pipeline first."),
                       options = list(dom = "t"), rownames = FALSE))
    }
    a <- res$audit
    st <- seq_type_resolved()
    unknown_col <- if (st == "nucleotide") "N Count" else "X Count"
    gc_col <- if (st == "nucleotide") "GC%" else "GC% (n/a)"
    names(a) <- c("Index", "Header", "Source", "Length", gc_col,
                   unknown_col, "Ambig%", "Status")
    datatable(a, options = list(pageLength = 25, scrollX = TRUE),
              rownames = FALSE, filter = "top") |>
      formatStyle("Status", backgroundColor = styleEqual(
        c("retained", "removed_duplicate", "filtered"),
        c("#eafaf1", "#fdedec", "#fef9e7")
      ), fontWeight = styleEqual("retained", "bold"))
  })
  
  # --- Run report ------------------------------------------------------------
  output$run_report_text <- renderText({
    res <- pipeline_out()
    params <- run_params()
    if (is.null(res) || is.null(params)) return("Run the pipeline to generate a report.")
    generate_run_report(params, res)
  })
  
  # --- Output preview --------------------------------------------------------
  output$output_preview <- renderText({
    res <- pipeline_out()
    if (is.null(res)) return("Run the pipeline to preview output.")
    deduped <- res$dedup$deduped
    n_show <- min(100, nrow(deduped))
    lines <- character()
    for (i in seq_len(n_show)) {
      lines <- c(lines, paste0(">", deduped$header[i]))
      s <- deduped$sequence[i]
      if (nchar(s) > 200) {
        lines <- c(lines, paste0(substr(s, 1, 200), "... [", nchar(s), " bp]"))
      } else {
        lines <- c(lines, s)
      }
    }
    if (nrow(deduped) > n_show) {
      lines <- c(lines, "", sprintf("... and %s more",
                                    format(nrow(deduped) - n_show, big.mark = ",")))
    }
    paste(lines, collapse = "\n")
  })
  
  # --- Downloads -------------------------------------------------------------
  output$download_buttons_ui <- renderUI({
    res <- pipeline_out()
    if (is.null(res)) return(NULL)
    tagList(
      div(class = "section-label", "Downloads"),
      downloadButton("dl_seqs", "Sequences", class = "btn-dl"),
      br(), br(),
      downloadButton("dl_audit", "Audit CSV", class = "btn-dl"),
      br(), br(),
      downloadButton("dl_clusters", "Clusters CSV", class = "btn-dl"),
      br(), br(),
      downloadButton("dl_report", "Run Report", class = "btn-dl")
    )
  })
  
  output$dl_seqs <- downloadHandler(
    filename = function() {
      fmt <- if (input$output_format == "auto") input_format() else input$output_format
      ext <- if (fmt == "fastq") "fastq" else "fasta"
      base <- tools::file_path_sans_ext(input$seq_files$name[1])
      paste0(base, "_deduped.", ext)
    },
    content = function(file) {
      res <- pipeline_out(); req(res)
      fmt <- if (input$output_format == "auto") input_format() else input$output_format
      if (fmt == "fastq") write_fastq_out(res$dedup$deduped, file)
      else write_fasta_out(res$dedup$deduped, file, line_width = input$line_width)
    }
  )
  
  output$dl_audit <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$seq_files$name[1]), "_audit.csv")
    },
    content = function(file) {
      res <- pipeline_out(); req(res)
      write.csv(res$audit, file, row.names = FALSE)
    }
  )
  
  output$dl_clusters <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$seq_files$name[1]), "_clusters.csv")
    },
    content = function(file) {
      res <- pipeline_out(); req(res)
      cs <- res$dedup$cluster_summary
      if (nrow(cs) > 0) write.csv(cs, file, row.names = FALSE)
      else write.csv(data.frame(Message = "No duplicates"), file, row.names = FALSE)
    }
  )
  
  output$dl_report <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$seq_files$name[1]),
             "_seqdedupe_report.txt")
    },
    content = function(file) {
      res <- pipeline_out(); params <- run_params(); req(res, params)
      writeLines(generate_run_report(params, res), file)
    }
  )
}

# =============================================================================
shinyApp(ui = ui, server = server)
