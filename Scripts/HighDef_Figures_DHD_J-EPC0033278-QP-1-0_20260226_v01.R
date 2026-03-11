# Script creating High definition images
# By David Demaree
# Physical Scientist, Edison NJ
# 3/10/2026

# We changed the publication that we were using and were using water research instead
# of water research X.  Water research requires high definition images of anything
# made out side of microsoft office.  So this remakes the figures in the proper scale.

# Before you run this script you have to run the processing script first.
# N_Am_Meta_DHD_J-EPC0033278-QP-1-0_20260311_v04.R
# Or if you ran the script previously you can run it from the workspace.

##############################################################################

# Forest plot


# Removes the USA subcategory from the Forestplot list.
ForestSubSort1 <- ForestSubSort %>% filter(Animal!="USA")


# Creates a dataframe from the ForestSubSort list.
anum <- nrow(ForestSubSort1)+1



ForestDat1 <- structure(list(mean  = c(NA, ForestSubSort1$Random.effects.mean), 
                            lower = c(NA, ForestSubSort1$Random.effects.95low),
                            upper = c(NA, ForestSubSort1$Random.effects.95high)),
                       .Names = c("mean", "lower", "upper"), 
                       row.names = c(NA, -anum),
                       class = "data.frame")
Foresttext1 <- cbind(c("Pathogen", ForestSubSort1$Pathogen),
                    c("Animal", ForestSubSort1$Animal),
                    c("Report#", ForestSubSort1$Report.),
                    c("Positive", ForestSubSort1$Positive),
                    c("Sample(n)", ForestSubSort1$Sample.n.),
                    c("GRADE", ForestSubSort1$GRADE))

# Sets the height of the graph
fheight <- 40*anum+264

#Multiplies ForestDat by 100
ForestDat100 <- ForestDat1*100

# Creates the Forest plot 
tiff(paste('HighDef/Figure4_FecalPathDHD_20260226v1.tiff'),
    # These adjust the size of the figure.
    width = 2244, height = fheight , units = "px") 
print(ForestDat100 |> 
        forestplot(labeltext = Foresttext1, 
                   is.summary = c(rep(TRUE, 1), rep(FALSE, anum - 1)),
                   line.margin = unit(-5, "cm"),
                   boxsize = 0.4,
                   xticks = c(0,20,40,60,80,100),
                   xlab = "% Prevalence",
                   # title = "Prevalence of pathogens in animal feces",

                   graphwidth = unit (20, "cm"),) |>
        fp_add_lines(h_2 = gpar(lty = 1),
                     h_6 = gpar(lty = 1),
                     h_10 = gpar(lty = 1),
                     h_13 = gpar(lty = 1),
                     h_17 = gpar(lty = 1),
                     h_20 = gpar(lty = 1),
                     h_24 = gpar(lty = 1),
                     h_30 = gpar(lty = 1))|>
        fp_set_style(box = "black",
                     line = "darkgrey",
                     summary = "black",
                     txt_gp = fpTxtGp(  cex = 2.5,
                                      ticks = gpar(fontfamily = "",  cex = 3),
                                      xlab = gpar(fontfamily = "", cex =3))
        )           

)

dev.off()

##############################################################################




# script that removes studies with USA data.
EnumboxDF1 <- EnumboxDF %>% filter(AniboxC!="USA")


# Creates a boxplot of the generated enumerated data in tiff format
# This boxplot is created off of generated data so if you regenerate the data
# in the parent script N_Am_Meta_DHD_J-EPC0033278-QP-1-0_20260311_v04.R you will
# get slightly different results.  But the results should be pretty close due
# to the sheer number points that were generated.

tiff(paste('HighDef/Figure5_FecalPathDHD_20260226v1.tiff'),
     width = 2244, height = 2244 , units = "px")
print(Enumboxplot <- ggplot(EnumboxDF1, aes(x=PathboxC, y=EnumRandC, fill=AniboxC)) + 
        geom_boxplot(position = position_dodge2(preserve = "single"),
                     outlier.size = 3) +
        geom_hline(yintercept=1)+
        facet_wrap(~PathboxC, scale="free") + 
        scale_y_continuous(trans='log10')+ 
        xlab("Subgroup")+
        ylab("Microbe/cyst/oocyst concentration per gram of wet feces")+
        theme_classic()+
        theme(legend.position = "bottom",
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              text = element_text(size = 60),
              legend.title = element_blank(),
              legend.text = element_text(margin = margin(r = 4, unit = "cm")),
              legend.key.size = unit(40, unit = "mm"),
              plot.title = element_text(size = 80))+
        scale_fill_grey(start = 0, end = 1,guide = guide_legend(label.hjust=5))
)

dev.off()


