library(shinythemes)
library(GenomicRanges)
library(Gviz)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm39)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(data.table)

# Mouse genome and gene annotation

genome <- BSgenome.Mmusculus.UCSC.mm39
txdb <- TxDb.Mmusculus.UCSC.mm39.refGene

axis_track <- GenomeAxisTrack()

gene_track_mouse <- GeneRegionTrack(
  txdb, transcriptAnnotation= "geneid", 
  stacking = "squish",
  name="Gene model",background.title="white",col.title="black",col.axis="black")

options(ucscChromosomeNames=FALSE)

# Human genome and gene annotation

genome_human <- BSgenome.Hsapiens.UCSC.hg38
txdb_human <- TxDb.Hsapiens.UCSC.hg38.refGene

gene_track_human <- GeneRegionTrack(
  txdb_human, genome=genome_human, showId=TRUE,transcriptAnnotation= "geneid", 
  stacking = "dense",
  name="Gene model",background.title="white",col.title="black",col.axis="black")
options(ucscChromosomeNames=FALSE)


plot_genome <- function(location,y_lim,y_lim_CLIP, dataset_ctr, dataset_cko, sashimi, sashimi_num, clip_dataset) {
  # Load data, at a reasonable level of detail for width(location)
  n <- min(width(location), 1000)
  
  tdp43_CLIP <- AlignmentsTrack(clip_dataset,
                                isPaired = FALSE, 
                                ylim=c(0,as.numeric(y_lim_CLIP)),
                                fill = "coral4",
                                #col="coral4",
                                name="TDP43 CLIP",
                                background.title="white",col.title="black",col.axis="black")
  
  ctr <- AlignmentsTrack(dataset_ctr,
                         isPaired = TRUE, #Paired-End Sequencing
                         ylim=c(0,as.numeric(y_lim)),
                         fill = "lightgrey",
                         #col="lightgrey",
                         name="Controls",
                         background.title="white",col.title="black",col.axis="black")
  
  cko <- AlignmentsTrack(dataset_cko,
                         isPaired = TRUE, #Paired-End Sequencing
                         ylim=c(0,as.numeric(y_lim)),
                         fill = "#a7c7fa",
                         #col="lightgrey",
                         name="TDP43 cKOs",
                         background.title="white",col.title="black",col.axis="black")
  
  plotTracks(
    list(axis_track,tdp43_CLIP, ctr, cko, gene_track_mouse),sizes=c(1.5,2,2,2,2),
    chromosome=as.character(seqnames(location)),
    from=start(location), to=end(location),type = c("coverage","sashimi"),
    sashimiNumbers = as.character(sashimi_num),
    sashimiScore = as.numeric(sashimi),
    add35=TRUE)
}

# ==================================================================


browser_ui <- function(id) {
  ns <- NS(id)
  
  div(
    h2("Select Datasets"),
    fluidRow(
      column(2,
             selectInput(ns('dataset_ctr'),
                         label = 'Choose a Controls Dataset:',
                         choices = c("Nestin-Cre Cortex e14"="/mnt/gtklab01/linglab/tdp43/STAR/tdp43_nestin_ctx_e14/tdp43_nestin_ctx_e14_ctr_merged.bam",
                                     "Nestin-Cre Thalamus e14"="/mnt/gtklab01/linglab/tdp43/STAR/tdp43_nestin_thlm_e14/tdp43_nestin_thlm_e14_ctr_merged.bam",
                                     "Camk2a-Cre Hippocampus p90"="/mnt/gtklab01/linglab/external_datasets/tdp43_camk2a_hip_p90/STAR/CAMK2A_HIP_P90_ctr_merged.bam",
                                     "CNP-Cre Sciatic nerve p21"="/mnt/gtklab01/linglab/tdp43/STAR/tdp43_cnp_scn_p21/tdp43_cnp_scn_p21_ctr_merged.bam",
                                     "CNP-Cre Spinal cord p60"="/mnt/gtklab01/linglab/tdp43/STAR/tdp43_cnp_spc_p60/tdp43_cnp_spc_p60_ctr.bam"),
                         selected = "Nestin-Cre Cortex e14")),
      column(2,
             selectInput(ns('dataset_cko'),
                         label = 'Choose a TDP43 cKOs Dataset:',
                         choices = c("Nestin-Cre Cortex e14"="/mnt/gtklab01/linglab/tdp43/STAR/tdp43_nestin_ctx_e14/tdp43_nestin_ctx_e14_cko_merged.bam",
                                     "Nestin-Cre Thalamus e14"="/mnt/gtklab01/linglab/tdp43/STAR/tdp43_nestin_thlm_e14/tdp43_nestin_thlm_e14_cko_merged.bam",
                                     "Camk2a-Cre Hippocampus p90"="/mnt/gtklab01/linglab/tdp43/STAR/CAMK2A_HIP_P90/CAMK2A_HIP_P90_cko_merged.bam",
                                     "CNP-Cre Sciatic nerve p21"="/mnt/gtklab01/linglab/tdp43/STAR/tdp43_cnp_scn_p21/tdp43_cnp_scn_p21_cko_merged.bam",
                                     "CNP-Cre Spinal cord p60"="/mnt/gtklab01/linglab/tdp43/STAR/tdp43_cnp_spc_p60/tdp43_cnp_spc_p60_cko_merged.bam"),
                         selected = "Nestin-Cre Cortex e14")),
      column(2,
             selectInput(ns('clip_dataset'),
                         label = 'Choose a TDP43 CLIP Dataset:',
                         choices = c("Mouse"="/mnt/gtklab01/linglab/tdp43/clipseq_tdp43/mouse/tdp43_clip_merged.bam",
                                     "Human"="/mnt/gtklab01/linglab/tdp43/clipseq_tdp43/human/Tollervey_2011_human_tdp43_clipseq.bam"),
                         selected = "Mouse"))),
    h2("Configure plot"),
    fluidRow(
      column(2,
             textInput(ns("location_str"), "Genomic location", "chr1:132532449-132539104")),
      column(2,
             textInput(ns("y_limit"), "Scale (Reads)", 1500)),
      column(2,
             textInput(ns("y_limit_CLIP"), "Scale (TDP43 CLIP)", 3000)),
      column(2,
             textInput(ns("min_sashimi"), "Minimum splice junction count", 100)),
      column(2,checkboxInput(ns("sashimi_num"), "Show splice junction counts", value = FALSE, width = NULL))),
    h2("Alignment"),
    plotOutput(ns("genome_plot"),width = "100%"),
    h2("Adjust plot"),
    fluidRow(
      column(6,
             actionButton(ns("go_left"), "<<"),
             actionButton(ns("go_right"), ">>"),
             actionButton(ns("zoom_in"), "Zoom in"),
             actionButton(ns("zoom_out"), "Zoom out"),
             downloadButton(ns("genome_alignment_plot"), "Download plot"))))
  
}

browser_server <- function(input, output, session) {
  
  location <- reactive({
    GRanges(input$location_str, seqinfo=seqinfo(genome))
  })
  
  y_lim <- reactive({
    input$y_limit
  })
  
  y_lim_CLIP <- reactive({
    input$y_limit_CLIP
  })
  
  sashimi <- reactive({
    input$min_sashimi
  })
  
  sashimi_num <- reactive({
    input$sashimi_num
  })
  
  dataset_ctr <- reactive({
    input$dataset_ctr
  })
  
  dataset_cko <- reactive({
    input$dataset_cko
  })
  
  clip_dataset <- reactive({
    input$clip_dataset
  })
  
  observeEvent(input$go_left, {
    amount <- max(1, width(location()) %/% 4)
    new_location <- shift(location(), -amount)
    new_location <- trim(new_location)
    updateTextInput(session, "location_str", value=as.character(new_location))
  })
  
  observeEvent(input$go_right, {
    amount <- max(1, width(location()) %/% 4)
    new_location <- shift(location(), amount)
    new_location <- trim(new_location)
    updateTextInput(session, "location_str", value=as.character(new_location))
  })
  
  observeEvent(input$zoom_in, {
    amount <- width(location()) %/% 4
    new_location <- location()
    start(new_location) <- start(new_location) + amount
    end(new_location) <- end(new_location) - amount
    updateTextInput(session, "location_str", value=as.character(new_location))
  })
  
  observeEvent(input$zoom_out, {
    amount <- width(location()) %/% 2
    new_location <- location()
    start(new_location) <- start(new_location) - amount
    end(new_location) <- end(new_location) + amount
    new_location <- trim(new_location)
    updateTextInput(session, "location_str", value=as.character(new_location))
  })
  
  
  output$genome_plot <- renderPlot({
    plot_genome( location(),y_lim(),y_lim_CLIP(), dataset_ctr(), dataset_cko(), sashimi(), sashimi_num(), clip_dataset())
  })
  
  output$genome_alignment_plot <- downloadHandler(
    filename = "genome_alignment_plot.png",
    content = function(file){
      png(file = file,width = 1000, height = 600)
      plot_genome( location(),y_lim(),y_lim_CLIP(),dataset_ctr(), dataset_cko(), sashimi(), sashimi_num(), clip_dataset())
      dev.off()
    }
  )
  
}
