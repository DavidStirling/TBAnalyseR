library(shiny)
library(kernlab)
library(dplyr)
completedanalysis = FALSE
# Setup constant values
splitnames = FALSE
batf2mean = -6.48989899
batf2stdev = 0.792471586
sweeneymean = -2.071616162
sweeneystdev = 0.412494173
sulimanmean = -3.901313131
sulimanstdev = 1.38555054
# Setup redcap lookup table
fullredcapnames = c('baseline_arm_1', 'longitudinal_m3_arm_1', 'longitudinal_m6_arm_1', 'longitudinal_m9_arm_1',
                    'longitudinal_m12_arm_1', 'longitudinal_m15_arm_1', 'longitudinal_m18_arm_1', 'longitudinal_m21_arm_1',
                    'longitudinal_m24_arm_1', 'treatment_week_4_arm_1', 'treatment_week_8_arm_1', 'treatment_week_12_arm_1', 'treatment_week_24_arm_1')
lookups = c('base', 'longm3', 'longm6', 'longm9', 'longm12', 'longm15', 'longm18', 'longm21', 'longm24', 'treatw4', 'treatw8', 'treatw12', 'treatw24')
lookupdf = data.frame(Type = lookups, redcap = fullredcapnames)
# Load SVM Model
load("svmmodelroe3")

runanalysis <- function(inputfile){
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
      split = unlist(strsplit(samplename,split="_"))
      name = split[1]
      type = split[2]
      date = split[3]
      # Replace type with lookup
      if (tolower(type) %in% lookupdf$Type) {
        type = lookupdf$redcap[match(tolower(type), lookupdf$Type)]
      }
      sampleresults = data.frame("UIN" = name, "Redcap_Event_Name" = type, "Sample_Date" = date, "BATF2_Value" = batf2score, "BATF2_Zscore" = batf2z,
                                 "Sweeney3_Value" = sweeney3, "Sweeney3_Zscore" = sweeneyz, "Suliman4_Value" = suliman4, "Suliman4_Zscore" = sulimanz)
    } else{
      name = as.character(datasheet[1,i])
      sampleresults = data.frame("UIN" = name, "BATF2_Value" = batf2score, "BATF2_Zscore" = batf2z,
                                 "Sweeney3_Value" = sweeney3, "Sweeney3_Zscore" = sweeneyz,"Suliman4_Value" = suliman4, "Suliman4_Zscore" = sulimanz)
    }
    
    databuffer = rbind(databuffer, sampleresults)
  }
  
  # Fetch SVM results
  # Transpose and comvert input file to a more useful format 
  inputcleaned <- datasheet[-c(2), -c(2:8)]
  transposed = t(inputcleaned)
  rownames(transposed) <- c()
  colnames(transposed) <- transposed[1, ]
  colnames(transposed)[1] = "Patient"
  transposed <- transposed[-1, ]
  transposed = transposed[, c(1:37)]
  nsdatasheet=as.data.frame(transposed)
  nsdatasheet_c = as.data.frame(apply(nsdatasheet, 2, as.character))
  nsdatasheet_n = as.data.frame(apply(nsdatasheet_c[,-1], 2, as.numeric))
  test_df = cbind("Patient" = nsdatasheet_c[,1], nsdatasheet_n)
  
  testdata <- test_df[,-c(1,2)] # remove the sample ID from the test data to run it through the model
  numobs_test <- dim(testdata)[1]
  
  # Cycle through samples and get SVM prediction
  pred.value_out <- rbind() # create the output list
  for (i in 1:numobs_test){
    test_data <- testdata[i,]
    # Predict with SVM model
    pred.value <- predict(svmmodel, test_data, type="decision")
    pred.value <- pred.value * -1 # Correct for inverted SVM value issue
    pred.value_out <- rbind(pred.value_out, pred.value) 
  }
  # Combine other test data with SVM results for output
  decisions <<- as.data.frame(pred.value_out)
  fulldata <<- as.data.frame(cbind(databuffer, "Roe3_Score"=decisions$V1))
  completedanalysis <<- TRUE
  return(fulldata)
}

ui <- fluidPage(
  titlePanel('TBAnalyseR'),
  sidebarLayout(
    sidebarPanel(
      fileInput("inputfile", "Upload Nanostring data file",  multiple = FALSE, accept = c("text/csv", "text/comma-separated-values", "text/plain", ".csv")),
      uiOutput("download")
      #downloadButton('downloadData', 'Download Results')
    ),
    mainPanel(
      p("This app is designed to process raw Nanostring output CSV files (log2-transformed) and perform diagnostic tests for TB based on gene expression."),
      tableOutput('table')
    )
  )
)

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
  
  # downloadHandler() takes two arguments, both functions.
  # The content function is passed a filename as an argument, and
  #   it should write out data to that filename.
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
