# app.R
library(shiny)
library(bslib)
library(fontawesome)
library(reticulate)
library(data.table)
library(DT)
library(plotly)
library(base64enc)

# Define the two sets of population codes
default_pops <- c("EUR", "AFR", "AMR", "EAS", "SAS")
hapmap_pops  <- c(
  "YRI",  # Yoruba in Ibadan, Nigeria
  "CHB",  # Han Chinese in Beijing, China
  "JPT",  # Japanese in Tokyo, Japan
  "CEU",  # CEPH/Utah (CEU)
  "MKK",  # Maasai in Kinyawa, Kenya
  "LWK",  # Luhya in Webuye, Kenya
  "CHD",  # Chinese in Denver, USA
  "GIH",  # Gujarati Indians in Houston, USA
  "TSI",  # Toscani in Italia
  "MXL",  # Mexican Ancestry in LA, USA
  "ASW"   # African Ancestry in SW USA
)

ui <- navbarPage(
  title = div(style = "display: flex; align-items: center;"),
  theme = bs_theme(
    version = 5,
    bootswatch = "united",
    bg = "#f8f9fa",
    fg = "#343a40",
    primary = "#a00000",
    secondary = "#343a40",
    base_font = font_google("Roboto"),
    heading_font = font_google("Lato")
  ) %>%
    bs_add_rules("
      /* Navbar gradient background */
      .navbar {
        background: linear-gradient(135deg, #a00000 0%, #343a40 100%) !important;
      }
      .navbar-brand, .nav-link {
        color: #fff !important;
        transition: background-color .3s;
      }
      .nav-link:hover, .nav-link.active {
        background-color: rgba(255,255,255,0.15) !important;
        border-radius: .25rem;
      }
      .card {
        border: none;
        border-radius: 1rem;
        box-shadow: 0 4px 12px rgba(0,0,0,0.1);
        margin-bottom: 1.5rem;
        transition: transform .2s, box-shadow .2s;
      }
      .card:hover {
        transform: translateY(-5px);
        box-shadow: 0 8px 20px rgba(0,0,0,0.15);
      }
      .card-header {
        background-color: #0a566e !important;
        color: #fff !important;
        border-top-left-radius: 1rem;
        border-top-right-radius: 1rem;
      }
      .form-control, .selectize-control.single .selectize-input {
        border-radius: .5rem !important;
      }
      .btn-danger {
        background: linear-gradient(135deg, #a00000 0%, #0a566e 100%) !important;
        border-color: #0a566e !important;
        box-shadow: 0 2px 6px rgba(0,0,0,0.2);
        transition: background .3s, box-shadow .3s;
      }
      .btn-danger:hover {
        background: linear-gradient(135deg, #0a566e 0%, #a00000 100%) !important;
        box-shadow: 0 4px 10px rgba(0,0,0,0.3);
      }
    "),
  fluid = TRUE,
  
  # ---------- Home tab ----------
  tabPanel(
    "Home",
    fluidPage(
      fluidRow(
        column(
          12, align = "center",
          div(
            class = "logo-container",
            tags$img(
              src = base64enc::dataURI(file = "logo.png", mime = "image/png"),
              alt = "LDSeeker Logo",
              style = "max-width: 180px; margin: 1rem auto;"
            )
          )
        )
      ),
      fluidRow(
        column(
          12,
          card(
            title = "Settings",
            status = "primary",
            height = "280px",
            width = NULL,
            style = "padding: 1.5rem;",
            fluidRow(
              column(2, radioButtons("pairwise", "Pairwise LD",
                                     choices = c("YES", "NO"),
                                     selected = "NO", inline = TRUE)),
              ## Changed fileInput from width=3 to width=2
              column(2, fileInput("input_file", "Upload Variant File",
                                  accept = c(".txt", ".tsv", ".csv"))),
              column(2, selectInput("ref_panel", "Reference Panel",
                                    choices = c(
                                      "Pheno Scanner"       = "Pheno_Scanner",
                                      "1000 Genomes (hg38)" = "1000G_hg38",
                                      "TOP-LD"              = "TOP_LD",
                                      "HapMap"              = "Hap_Map"
                                    ),
                                    selected = "1000G_hg38")),
              column(1, numericInput("r2threshold", "R²",
                                     value = 0.8, min = 0, max = 1,
                                     step = 0.01)),
              column(2, selectInput("population", "Population",
                                    choices = default_pops,
                                    selected = default_pops[1])),
              column(2, numericInput("maf", "MAF",
                                     value = 0.01, min = 0, max = 1,
                                     step = 0.01)),
              column(
                1,
                style = "display: flex; align-items: center;",
                actionButton("run_ld", "Run",
                             class = "btn btn-danger",
                             style = "width:100%;")
              )
            )
          )
        )
      ),
      fluidRow(
        column(
          12,
          card(
            title = "LDSeeker Results",
            status = "info",
            width = NULL,
            DTOutput("ld_table")
          )
        )
      )
    )
  ), 
  # ---------- LD Heatmap tab ----------
  tabPanel(
    "LD Heatmap",
    fluidPage(
      titlePanel("LD Heatmap Viewer (R² Matrix)"),
      
      sidebarLayout(
        sidebarPanel(
          fileInput("ld_heat_file", "Upload LD File (.txt)", 
                    accept = c(".txt", ".tsv")),
          sliderInput("r2_range", "Adjust R² Color Scale", 
                      min = 0, max = 1, value = c(0, 1), step = 0.01),
          checkboxInput("show_values", "Show R² values on tiles", value = FALSE),
          selectInput("color_scale", "Select Colorscale",
                      choices = c("Viridis" = "Viridis",
                                  "Cividis" = "Cividis",
                                  "Plasma" = "Plasma",
                                  "Inferno" = "Inferno",
                                  "Magma" = "Magma",
                                  "Jet" = "Jet",
                                  "Hot" = "Hot",
                                  "Cool" = "Cool",
                                  "RdBu" = "RdBu",
                                  "YlGnBu" = "YlGnBu",
                                  "Greys" = "Greys"),
                      selected = "RdBu")
        ),
        
        mainPanel(
          plotlyOutput("heatmap", height = "800px")
        )
      )
    )
  ),
  # ---------- LD Heatmap tab ----------
  
  
  
  
  # ---------- Help tab ----------
  tabPanel(
    "Help",
    fluidPage(
      h3("Help & Documentation"),
      h5("Overview"),
      h5("Use LDSeeker to explore linkage disequilibrium metrics for your variant dataset. Configure settings on the Home tab, click **Run**, then review results."),
     h5(" Settings Controls"),
      tags$ul(
        tags$li(strong("Pairwise LD:"), " YES = compute pairwise r² for all SNP pairs; NO = Finds LD among all the SNPs."),
        tags$li(strong("Upload rsIDs file:"), " Supported: .txt, .tsv, .csv. Must include SNP and CHR columns."),
        tags$li(strong("Reference Panel:"), " Choose population reference: Pheno_Scanner, 1000G_hg38, TOP_LD, Hap_Map."),
        tags$li(strong("R² Threshold:"), " Numeric input (0–1). Filters LD output by minimum r² (default = 0.8)."),
        tags$li(strong("Population:"), " Select from EUR, AFR, AMR, EAS, SAS for allele frequencies (or HapMap symbols if HapMap is chosen)."),
        tags$li(strong("MAF:"), " Minor allele frequency filter between 0–1 (default = 0.01)."),
        tags$li(strong("Run Button:"), " Executes the  LDSeeker with the selected parameters.")
      ),
      h5("Results Table"),
      tags$ul(
        tags$li("Displays LD metrics using DataTables. Use the buttons above to copy, download (CSV/Excel/PDF)."),
        tags$li("Columns include: SNP, CHR, BP, A1, A2, R² (or R, Dprime), and any summary stats.")
      )
    )
  ),
  

  # ---------- About tab ----------
  tabPanel(
    "About",
    # wrap in fluidPage so we can inject CSS
    fluidPage(
      # CSS to style our logos
      tags$style(HTML("
        .about-logo {
          max-height: 150px;   /* limit height to 150px */
          width: auto;         /* auto width to keep aspect ratio */
          margin: 1rem;        /* optional spacing */
        }
        .about-logos-row {
          display: flex;
          justify-content: center;
          align-items: center;
          flex-wrap: wrap;
          gap: 2rem;
          margin-top: 1rem;
        }
      ")),
      fluidRow(
        column(
          width = 12,
          h2("Welcome to LDSeeker!"),
         h4("LDSeeker: An open source tool for seeking LD"),
         h4("Correspondence: Dr. Pantelis G. Bagos (email: pbagos@compgen.org)"),
         h4("Download the source code of LDSeeker from ",
            tags$a(href = "https://github.com/gmanios/LDSeeker", "GitHub", target = "_blank")
          ),
         h4(tags$a(href = "https://sites.google.com/compgen.org/lab/home?authuser=0",
                   "Computational Genetics Group", target = "_blank")),
         h4("Main developer:"),
          tags$ul(
            tags$li("Georgios A. Manios - Department of Computer Science and Biomedical Informatics,
                      University of Thessaly (email: gmanios@uth.gr)")
          ),
          # logos in a flex container to center them
          div(class = "about-logos-row",
              tags$img(
                src = base64enc::dataURI(file = "uth.png", mime = "image/png"),
                alt = "UTH Logo",
                class = "about-logo"
              ),
              tags$img(
                src = base64enc::dataURI(file = "cgg.png", mime = "image/png"),
                alt = "CGG Logo",
                class = "about-logo"
              )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # ========== LD Heatmap Logic ===================
  # ====== LD Heatmap logic ======
  
  ld_heat_data <- reactive({
    req(input$ld_heat_file)
    df <- fread(input$ld_heat_file$datapath, header = TRUE, sep = "\t", data.table = FALSE)
    
    # Normalize column names
    if ("SNP_A" %in% names(df)) names(df)[names(df) == "SNP_A"] <- "rsID1"
    if ("SNP_B" %in% names(df)) names(df)[names(df) == "rsID2"]
    
    validate(
      need(all(c("rsID1", "rsID2", "R2") %in% names(df)),
           "File must contain rsID1, rsID2, and R2 columns (or SNP_A, SNP_B).")
    )
    
    df
  })
  
  heatmap_data <- reactive({
    df <- ld_heat_data()
    
    snps <- sort(unique(c(df$rsID1, df$rsID2)))
    mat <- matrix(NA, nrow = length(snps), ncol = length(snps),
                  dimnames = list(snps, snps))
    
    for (i in 1:nrow(df)) {
      a <- df$rsID1[i]
      b <- df$rsID2[i]
      r2 <- df$R2[i]
      mat[a, b] <- r2
      mat[b, a] <- r2
    }
    
    mat[is.na(mat)] <- 0
    diag(mat) <- 1.00
    
    for (i in 1:nrow(mat)) {
      for (j in 1:ncol(mat)) {
        if (j > i) {
          mat[i, j] <- NA
        }
      }
    }
    
    mat
  })
  
  get_font_color <- function(value, r2_min, r2_max) {
    normalized <- (value - r2_min) / (r2_max - r2_min)
    if (normalized > 0.5) return("black")
    else return("white")
  }
  
  output$heatmap <- renderPlotly({
    mat <- heatmap_data()
    r2_min <- input$r2_range[1]
    r2_max <- input$r2_range[2]
    selected_colorscale <- input$color_scale
    
    snps <- rownames(mat)
    
    hover_text <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
    for (i in seq_along(snps)) {
      for (j in seq_along(snps)) {
        val <- mat[i, j]
        if (!is.na(val)) {
          hover_text[i, j] <- paste0("rsID1: ", snps[i], "<br>rsID2: ", snps[j], "<br>R²: ", round(val, 3))
        }
      }
    }
    
    p <- plot_ly(
      x = snps, y = snps, z = mat,
      type = "heatmap",
      colorscale = selected_colorscale,
      zmin = r2_min,
      zmax = r2_max,
      text = hover_text,
      hoverinfo = "text",
      showscale = TRUE
    )
    
    if (input$show_values) {
      annots <- list()
      for (i in seq_along(snps)) {
        for (j in seq_along(snps)) {
          val <- mat[i, j]
          if (!is.na(val)) {
            text_col <- get_font_color(val, r2_min, r2_max)
            annots[[length(annots) + 1]] <- list(
              x = snps[j],
              y = snps[i],
              text = sprintf("%.2f", val),
              showarrow = FALSE,
              font = list(color = text_col, size = 10)
            )
          }
        }
      }
      p <- layout(p, annotations = annots)
    }
    
    layout(p,
           xaxis = list(title = "SNP", tickangle = -90, scaleanchor = "y"),
           yaxis = list(title = "SNP", autorange = "reversed"),
           margin = list(l = 120, b = 120),
           title = "Interactive LD Heatmap (R²)"
    )
  })
  #============================================================================
  rv <- reactiveValues(ld = NULL)
  
  # ========== Dynamic population picker ==========
  observeEvent(input$ref_panel, {
    if (input$ref_panel == "Hap_Map") {
      updateSelectInput(
        session, "population",
        label   = "Population (HapMap)",
        choices = hapmap_pops,
        selected = hapmap_pops[1]
      )
    } else {
      updateSelectInput(
        session, "population",
        label   = "Population",
        choices = default_pops,
        selected = default_pops[1]
      )
    }
  })
  
  # ========== Run LDSeeker ==========
  observeEvent(input$run_ld, {
    req(input$input_file)
    withProgress(message = "Running LDSeeker...", value = 0, {
      incProgress(0.1, "Preparing inputs...")
      input_path <- input$input_file$datapath
      
      incProgress(0.1, "Building command...")
      cmd <- paste(
        "python LDSeeker.py",
        "--file-path", shQuote(input_path),
        "--r2threshold", input$r2threshold,
        "--pop", input$population,
        "--maf", input$maf,
        "--ref", input$ref_panel,
        "--pairwise", input$pairwise
      )
      
      incProgress(0.1, "Executing script...")
      result_log <- tryCatch(
        paste(system(cmd, intern = TRUE), collapse = "\n"),
        error = function(e) e$message
      )
      incProgress(0.4, "Completed")
      output$result_output <- renderText(result_log)
      
      file_name <- if (input$pairwise == "YES")
        "LD_info_chr_all_pairwise.txt" else "LD_info_chr_all.txt"
      txt_path <- file.path(getwd(), file_name)
      incProgress(0.1, paste("Reading", file_name))
      df <- if (file.exists(txt_path)) {
        tryCatch(fread(txt_path), error = function(e) data.table(Error = e$message))
      } else {
        data.table(Message = paste(file_name, "not found."))
      }
      rv$ld <- df
      
        # ← HERE: set server = FALSE so all rows are sent to client for Buttons export
      output$ld_table <- renderDT(
        datatable(
          rv$ld,
          extensions = 'Buttons',
          options = list(
            pageLength = 10,
            dom        = 'Bfrtip',
            buttons    = list(
              list(extend='copy',  text='Copy',
                   exportOptions=list(modifier=list(page='all'))),
              list(extend='csvHtml5', text='CSV',  filename='LDSeeker_results',
                   exportOptions=list(modifier=list(page='all'))),
              list(extend='excelHtml5', text='Excel',filename='LDSeeker_results',
                   exportOptions=list(modifier=list(page='all'))),
              list(extend='csvHtml5', text='TXT',   filename='LDSeeker_results',
                   extension='.txt', fieldSeparator='\t',
                   exportOptions=list(modifier=list(page='all')))
            )
          )
        ),
        server = FALSE
      )
    })
  })
}

shinyApp(ui, server)
