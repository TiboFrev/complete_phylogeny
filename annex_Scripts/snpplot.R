vcfC <- readLines("C:\\Users\\frevillet\\Desktop\\SNPpositions.txt")
library(stringr)

segment <- c()
nombres <- c()

##62 genomes 3seq
portions <- c('135051-135065-143895-144102', '116174-118590-133754-134405', '121690-121728-136855-136875', '119365-119403-138310-138364', '4589-4778-135051-135065', '17306-17401-146483-146683', '91-2017-150562-150562')
frequences <- c(16, 33, 1, 2, 1, 1, 1)

##66 genomes 3seq
#portions <- c('135051-135065-143895-144102', '116174-118590-133754-134405', '2425-2425-148371-148612', '167-169-148371-148612', '1904-2425-148827-148915', '119365-119403-138310-138364', '4589-4778-135051-135065', '17306-17401-146483-146683', '91-2017-150562-150562')
#frequences <- c(18, 33, 1, 1, 1, 2, 1, 1, 1)

##62 genomes geneconv
#portions <- c('134818-139288', '135052-139288', '135917-141125', '136195-138364', '18514-56346', '136126-138364', '2297-17401', '2297-8642', '2297-8157', '18760-56346', '18912-56346', '66825-84054', '66825-84334', '4590-19554', '4590-20056', '136992-140558', '68715-83591', '67410-83428', '111048-119700', '119366-131542', '137079-140558', '138311-148393', '68715-83428', '123717-131542', '111048-119520', '67410-84352', '138311-148235', '111048-119584', '2523-7010', '122535-135065', '135917-139288', '134221-138077', '121691-127633', '132820-136875', '4590-19555', '127716-131838', '4590-7304', '121727-127633', '5159-8208', '133031-136875', '68792-83428', '85113-95398', '3594-7010', '122080-127633', '11762-17401', '133119-136875', '133068-136875', '140565-148331', '12242-17401', '122098-127633', '138424-142031', '27549-45675', '143896-148393', '85599-95398', '103652-119403', '27962-45675', '103363-119403')
#frequences <- c(168, 231, 19, 24, 7, 56, 10, 7, 1, 11, 1, 7, 140, 372, 14, 16, 7, 173, 6, 7, 2, 34, 11, 31, 11, 79, 2, 1, 20, 378, 21, 2, 7, 8, 4, 18, 6, 8, 1, 8, 1, 14, 1, 3, 2, 2, 1, 7, 1, 1, 1, 11, 16, 4, 2, 3, 1)




for (i in 3:length(vcfC)-1){
  segment <- c(segment, as.integer(str_extract(vcfC[i], "\\d+")))
  nombres <- c(nombres, str_remove(str_extract(vcfC[i], ":\\s+(\\d+)"), ":\\s+"))
}

get_bounds <- function(portion) {
  as.numeric(unlist(strsplit(portion, "-")))
}


colors = rainbow(length(frequences))
png(file = "Dplot_3seq_66.png", width = 1900, height = 1000, units = "px")
par(fig=c(0,1,0.2,1))
plot(nombres, type = "s", xaxt = "n", main = "Density of the SNPs across the LSDV genome", cex.main = 2, xlab = "", ylab = "number of SNP",cex.lab = 1.5)  
axis(side = 1, at = seq_along(segment[1:length(segment)]), labels = segment[1:length(segment)], las = 2)
abline(v = 13.851, col = "red")
text(14, 60, "beginning of core genome", col = "red", pos=4)
abline(v = 106.910, col = "red")
text(107.5, 60, "ending of core genome", col = "red",pos=4)

# Créer un nouvel espace graphique pour les segments
par(fig=c(0,1,0,0.2), new=TRUE)
plot.window(xlim = c(min(segment)/1000, max(segment)/1000), ylim = c(0, sum(frequences)),ylab = "number of recombination events")
for (i in 1:length(portions)) {
  bounds <- get_bounds(portions[i])
  color = colors[i]
  for (j in 1:frequences[i]) {
    segments(bounds[1]/1000, 3*j, bounds[2]/1000, 3*j,lwd=1,col = color)
    segments(bounds[3]/1000, 3*j, bounds[4]/1000, 3*j,lwd=1,col = color)
  }
}
dev.off()

######################

colors <- rainbow(length(frequences))
png(file = "test_3seq_62.png", width = 1900, height = 1000, units = "px")
par(fig = c(0, 1, 0.1, 1))
plot(nombres, type = "s", xaxt = "n", main = "Density of the SNPs across the LSDV genome", cex.main = 2, xlab = "", ylab = "number of SNP", cex.lab = 1.5)
abline(v = 13.851, col = "red")
text(14, 60, "beginning of core genome", col = "red", pos=4)
abline(v = 106.910, col = "red")
text(107.5, 60, "ending of core genome", col = "red",pos=4)

axis(side = 1, at = seq_along(segment[1:length(segment)]), labels = segment[1:length(segment)], las = 2)


# Créer un nouvel espace graphique pour les segments
par(fig = c(0, 1, 0, 0.1), new = TRUE)
plot.window(xlim = c(min(segment)/1000, max(segment)/1000), ylim = c(0, sum(frequences)), ylab = "number of recombination events")
for (i in 1:length(portions)) {
  bounds <- get_bounds(portions[i])
  color <- colors[i]
    rect(bounds[1]/1000, 1 , bounds[2]/1000, 100, col = color, border = NA)
    rect(bounds[3]/1000, 1, bounds[4]/1000, 100, col = color, border = NA)
}
dev.off()





