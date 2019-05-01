library(shiny)
library(kernlab)
library(dplyr)
# Setup constant values
completedanalysis = FALSE
splitnames = FALSE
batf2mean = -6.48989899
batf2stdev = 0.792471586
sweeneymean = -2.071616162
sweeneystdev = 0.412494173
sulimanmean = -3.901313131
sulimanstdev = 1.38555054
roemean = -49.99192115
roestdev = 1.644894428
zakmean = -2.919968884
zakstdev = 0.815010259

# Setup redcap lookup table
fullredcapnames = c('baseline_arm_1', 'longitudinal_m3_arm_1', 'longitudinal_m6_arm_1', 'longitudinal_m9_arm_1',
                    'longitudinal_m12_arm_1', 'longitudinal_m15_arm_1', 'longitudinal_m18_arm_1', 'longitudinal_m21_arm_1',
                    'longitudinal_m24_arm_1', 'treatment_week_4_arm_1', 'treatment_week_8_arm_1', 'treatment_week_12_arm_1', 'treatment_week_24_arm_1')
lookups = c('base', 'longm3', 'longm6', 'longm9', 'longm12', 'longm15', 'longm18', 'longm21', 'longm24', 'treatw4', 'treatw8', 'treatw12', 'treatw24')
lookupdf = data.frame(Type = lookups, redcap = fullredcapnames)
# Load SVM Model
load("Roe3.model") #roe3model
load("Zak16.model") #zak16model

runanalysis <- function(inputfile){

  normgapdh <- function(x) {
    return(x - test_data$GAPDH)
  }
  
    completedanalysis <<- FALSE
  # Try to import the data sheet
  datasheet = read.csv(inputfile, 1, stringsAsFactors = FALSE, as.is = TRUE)
  
  # Check if names are in 3-part format
  if (lengths(regmatches(datasheet[1,9], gregexpr("_", datasheet[1,9]))) == 2){
    splitnames = TRUE
  }

  # Get IDs of rows with genes of interest
  idbatf2 = which(datasheet$Probe.Name == "BATF2")
  idgbp5 = which(datasheet$Probe.Name == "GBP5")
  iddusp3 = which(datasheet$Probe.Name == "DUSP3")
  idklf2 = which(datasheet$Probe.Name == "KLF2")
  idgas6 = which(datasheet$Probe.Name == "GAS6")
  idsept4 = which(datasheet$Probe.Name == "SEPT4.")
  idcd1c = which(datasheet$Probe.Name == "CD1C")
  idblk = which(datasheet$Probe.Name == "BLK")
  idgapdh = which(datasheet$Probe.Name == "GAPDH")
  
  # Setup data container
  databuffer = data.frame()
  
  # Cycle through samples
  for (i in seq(9, ncol(datasheet))){
    # Obtain correct row IDs for each gene
    batf2 = as.numeric(as.character(datasheet[idbatf2,i]))
    gapdh = as.numeric(as.character(datasheet[idgapdh,i]))
    gbp5 = as.numeric(as.character(datasheet[idgbp5,i]))
    dusp3 = as.numeric(as.character(datasheet[iddusp3,i]))
    klf2 = as.numeric(as.character(datasheet[idklf2,i]))
    gas6 = as.numeric(as.character(datasheet[idgas6,i]))
    sept4 = as.numeric(as.character(datasheet[idsept4,i]))
    cd1c = as.numeric(as.character(datasheet[idcd1c,i]))
    blk = as.numeric(as.character(datasheet[idblk,i]))
    
    # Normalise to GAPDH
    batf2 = batf2 - gapdh
    gbp5 = gbp5 - gapdh
    dusp3 = dusp3 - gapdh
    klf2 = klf2 - gapdh
    gas6 = gas6 - gapdh
    sept4 = sept4 - gapdh
    cd1c = cd1c - gapdh
    blk = blk - gapdh
    
    # Perform Calculations
    batf2score = batf2
    sweeney3 = ((gbp5 + dusp3)/2) - klf2
    suliman4 = (gas6 + sept4) - (cd1c + blk)
    
    # Calculate Z-Scores
    batf2z = (batf2 - batf2mean) / batf2stdev
    sweeneyz = (sweeney3 - sweeneymean) / sweeneystdev
    sulimanz = (suliman4 - sulimanmean) / sulimanstdev
    
    
    # Store results, parsing sample data if it's present in the sample name
    if (splitnames == TRUE) {
      samplename = as.character(datasheet[1,i])
      if (lengths(regmatches(samplename, gregexpr("_", samplename))) == 2) {
      split = unlist(strsplit(samplename,split="_"))
      name = split[1]
      type = split[2]
      date = split[3]
      } else {
        name = samplename
        type = "undefined"
        date = "undefined"
      }
      # Replace type with lookup
      if (tolower(type) %in% lookupdf$Type) {
        type = lookupdf$redcap[match(tolower(type), lookupdf$Type)]
      }
      sampleresults = data.frame("uin" = name, "redcap_event_name" = type, "sample_date" = date, "batf2_value" = batf2score, "batf2_zscore" = batf2z,
                                 "sweeney3_value" = sweeney3, "sweeney3_zscore" = sweeneyz, "suliman4_value" = suliman4, "suliman4_zscore" = sulimanz)
    } else{
      name = as.character(datasheet[1,i])
      sampleresults = data.frame("uin" = name, "batf2_value" = batf2score, "batf2_zscore" = batf2z,
                                 "sweeney3_value" = sweeney3, "sweeney3_zscore" = sweeneyz,"suliman4_value" = suliman4, "suliman4_zscore" = sulimanz)
    }
    
    databuffer = rbind(databuffer, sampleresults)
  }
  
  # Fetch SVM results
  # Transpose and comvert input file to a more useful format 
  inputcleaned <- datasheet[-c(2), -c(2:8)]
  transposed = t(inputcleaned)
  rownames(transposed) <- c()
  colnames(transposed) <- transposed[1, ]
  colnames(transposed)[1] = "patient"
  transposed <- transposed[-1, ]
  transposed = transposed[, c(1:37)]
  nsdatasheet=as.data.frame(transposed)
  nsdatasheet_c = as.data.frame(apply(nsdatasheet, 2, as.character))
  nsdatasheet_n = as.data.frame(apply(nsdatasheet_c[,-1], 2, as.numeric))
  test_df = cbind("patient" = nsdatasheet_c[,1], nsdatasheet_n)
  testdata <- test_df[,-c(1)] # remove the sample ID from the test data to run it through the model
  numobs_test <- dim(testdata)[1]
  
  # Convert names of columns to correct format for SVM
  names(testdata)[names(testdata)=="FCGR1A (CD64)"] <- "FCGR1A"
  names(testdata)[names(testdata)=="SEPT4."] <- "SEPT4"

  
  # Cycle through samples and get SVM prediction
  svmresults <- rbind() # create the output list
  for (i in 1:numobs_test){
    test_data <- testdata[i,]
    test_data[] <- lapply(test_data, normgapdh)
    # Predict with SVM model
    roe_value <- predict(roe3model, test_data, type="decision")
    roe_value <- roe_value * -1 # Correct for inverted SVM value issue
    roe_zscore <- (roe_value - roemean) / roestdev

    zak_value <- predict(zak16model, test_data, type="decision")
    #zak_value <- zak_value * -1 # No need to correct for inverted SVM value issue
	  zak_zscore <- (zak_value - zakmean) / zakstdev
	  
    resultrow = cbind(roe_value, roe_zscore, zak_value, zak_zscore)
    svmresults <- rbind(svmresults, resultrow) 
  }
  # Combine other test data with SVM results for output
  decisions <<- as.data.frame(svmresults)
  fulldata <<- as.data.frame(cbind(databuffer, "roe3_value" = decisions$V1, "roe3_zscore" = decisions$V2, "zak16_value" = decisions$V3, "zak16_zscore" = decisions$V4))
  completedanalysis <<- TRUE
  return(fulldata)
}

# Shiny App UI
ui <- fluidPage(
  titlePanel('TBAnalyseR'),
  sidebarLayout(
    sidebarPanel(
      fileInput("inputfile", "Upload Nanostring data file",  multiple = FALSE, accept = c("text/csv", "text/comma-separated-values", "text/plain", ".csv")),
      uiOutput("download")
    ),
    mainPanel(
      p(strong("Update 1.2"), "- Zak-16 scores added. (01/05/19)"),
      p("This app is designed to process raw Nanostring output CSV files (log2-transformed) and perform diagnostic tests for TB based on gene expression."),
      tableOutput('table')
    )
  )
)

# Shiny App Server
server <- function(input, output) {
  generateResults <- reactive({
    file = input$inputfile
    try(df <- runanalysis(file$datapath))
    validate(
      need(typeof(df) != "closure", "Analysis failed. Check that input file is an unaltered Nanostring output file")
    )
    return(df)
  })
  
  output$table <- renderTable({
    req(input$inputfile)
    generateResults()
  })

  output$downloadData <- downloadHandler(
      filename = function() {
      origfile = input$inputfile
      paste(substr(origfile$name, 1, nchar(origfile$name)-4),"_results.csv", sep = "")
    },
    content = function(file) {
      # Write to a file specified by the 'file' argument
      write.csv(generateResults(), file, row.names = FALSE)
    }
  )

  output$download <- renderUI({
    if(!is.null(input$inputfile) & (completedanalysis == TRUE)) {
      downloadButton('downloadData', 'Download Results')
    }
  })
  }

shinyApp(ui = ui, server = server)

