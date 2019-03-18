con <- file("C:/Users/zty20/Documents/R/Protein_AD/inbiomap.txt", "r")
line=readLines(con,n=1)
while( length(line) != 0 ) {
	 if (substr(line, 3, 7) =='<edge'){
		protein1=substr(line, 25, 30) 
		protein2=substr(line, 13, 18) 
		line=readLines(con,n=1)
		score=substr(line, 41, 45)
		result=data.frame(protein1,protein2,score)
		write.table(file='net_result.txt',result,append=T,col.names=F,row.names=F,quote = FALSE)
	}
	else {
	line=readLines(con,n=1)
	}
	 
}
close(con)