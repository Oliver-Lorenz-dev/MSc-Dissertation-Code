# create venn diagram seen in figure 17
# numbers obtained from python code (see python_code/venn_diagram folder)

library(VennDiagram)

draw.triple.venn(area1 = 4803, area2 = 16476, area3 = 5781,
                 n12 = 2993, n13 = 718, n23 = 3683, n123 = 479,
                 alpha = rep(0.3, 3), fill = c("yellow", "red", "blue"),
                 category = c("Normal v Embryonic"
                              ,"BPH v Normal", "Cancer v BPH"))