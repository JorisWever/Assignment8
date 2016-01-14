importFromUrl <- function(url){
	# create data folder (if it doesnt exist already)
	dir.create('data', showWarnings = FALSE)
	
	# create name variables
	fName <- basename(url)
	outName <- paste('data/',fName)
	
	if(!file.exists(outName)){
  	# download and unzip the data
  	download.file(url = url, destfile = outName, method = 'wget')
  }
	unzip(outName, exdir = 'data')
}