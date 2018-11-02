windnetplot <-
function () 
{

#
# Blank plot for wordlayout to work with
#

plot(x=NA, xlim=c(160,635), ylim=c(25,425), axes=TRUE, asp=1, xlab="", ylab="")

#
# Work out label layout and coordinates
#

labellay <- wordlayout(x=GNAR::vswindcoords[,1], y=GNAR::vswindcoords[,2],
                       words=GNAR::vswindnames)
labelpos <- cbind(c(labellay[,1]+0.5*labellay[,3]), c(labellay[,2]+0.5*labellay[,4]))

#
# Plot the network, and distances and names of wind stations
#

plot(GNAR::vswindnetD, layout=labelpos, vertex.label=GNAR::vswindnames,
     vertex.size=0, vertex.label.color="black", edge.lty=2, edge.color="red",
     edge.label=round(1/E(GNARtoigraph(GNAR::vswindnetD))$weight,digits=2),
     edge.label.cex=0.6,
     rescale=FALSE, xlim=c(160,635), ylim=c(25,425), axes=TRUE, asp=1, add=TRUE)

}
