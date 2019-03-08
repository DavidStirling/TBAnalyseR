# Nanostring Data Analyser by David Stirling, 2018-2019

# This is a script to process nanostring data in R instead of using excel.
# Takes in an unaltered csv file produced by the Nanostring (log2 transformed by Nanostring software)
# Puts out a CSV file containing analysis scores for BATF2, Sweeney3, Suliman4 and Roe5 TB tests

#####

library(kernlab)

# Setup constant values
splitnames = FALSE
batf2mean = -6.48989899
batf2stdev = 0.792471586
sweeneymean = -2.071616162
sweeneystdev = 0.412494173
sulimanmean = -3.901313131
sulimanstdev = 1.38555054

load("C:/Users/daves/Desktop/Test Stage 2/Nanostring/SVM/svmmodelroe3")

# Prompt user for input and output files
input = choose.files(caption = "Select csv file to import...", multi = FALSE, filters = c("Comma Delimited Files (.csv)","*.csv"))
outputname = choose.files(caption = "Choose file to save data to...", multi = FALSE, filters = c("Comma Delimited Files (.csv)","*.csv"))


# Sanity check input/output files
if (identical(input,character(0))){
  winDialog(type = c("ok"),"Input file not set, aborting run")
  stop("Input file not set, aborting run")
  }
if (identical(outputname,character(0))){
  winDialog(type = c("ok"),"Output file not set, aborting run")
  stop("Output file not set, aborting run")
  }
if (!file.exists(input)){
  winDialog(type = c("ok"),"Specified input file does not exist, aborting run")
  stop("Specified input file does not exist, aborting run")
  } 
print("Parameters set.")

# Try to import the data sheet
datasheet = read.csv(input, 1, stringsAsFactors = FALSE, as.is = TRUE)
print("Data imported.")

# Check if names are in 3-part format
if (lengths(regmatches(datasheet[1,9], gregexpr("_", datasheet[1,9]))) == 2){
  splitnames = TRUE
  print("Detected splittable names")
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
    sampleresults = data.frame("UIN" = name, "Redcap_Event_Name" = type, "Sample_Date" = date, "BATF2_Value" = batf2score, "BATF2_Zscore" = batf2z, "Sweeney3_Value" = sweeney3, "Sweeney3_Zscore" = sweeneyz, "Suliman4_Value" = suliman4, "Suliman4_Zscore" = sulimanz)
  } else{
    name = as.character(datasheet[1,i])
    sampleresults = data.frame("UIN" = name, "BATF2_Value" = batf2score, "BATF2_Zscore" = batf2z, "Sweeney3_Value" = sweeney3, "Sweeney3_Zscore" = sweeneyz, "Suliman4_Value" = suliman4, "Suliman4_Zscore" = sulimanz)
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
databuffer <<- as.data.frame(cbind(databuffer, "Roe3_Score"=decisions$V1))

# Write Results File
write.csv(databuffer, outputname)
print("Analysis Complete.")
print(paste("Data saved to", outputname))
winDialog(type = c("ok"),"Analysis Complete")
