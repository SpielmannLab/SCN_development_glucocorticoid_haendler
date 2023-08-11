library(ggplot2)
colours <- c("#FDE725", 
"#B4DE2C", 
"#6DCD59", 
"#35B779", 
"#1F9E89", 
"#26828E", 
"#31688E", 
"#3E4A89", 
"#482878", 
"#440154")

colfunc <- colorRamp(c("#ff5733", "#838383", "#355edf"))
cols <- colfunc(seq(0,1,len=11))
rgb(cols[,1], cols[,2], cols[,3], maxColorValue = 255) %>% paste(collapse = ", ")

# Use the output of this for Sankeyplot.

# Sankeyplot color examples at:
https://r-graph-gallery.com/322-custom-colours-in-sankey-diagram.html


colrs2 <- c("#ff5733", "#838383", "#355edf")
colors2 <- c("#FF5733", "#E65F43", "#CD6852", "#B47163", "#9B7A72", "#838383", "#737B95", "#6374A7", "#546CBA", "#4465CC", "#355EDF")