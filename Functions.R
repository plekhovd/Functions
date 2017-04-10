###Function for connecting to Pleiades site and downloading .csv files
pleiadesFunction = function(file){
  con = gzcon(url(paste("http://atlantides.org/downloads/pleiades/dumps/", file, sep="")))
  txt = readLines(con)
  dat = read.csv(textConnection(txt))
}