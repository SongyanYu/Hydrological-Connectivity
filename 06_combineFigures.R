#---
# R tutorial for advanced image processing with package magick.
# tutoiral website: https://cran.r-project.org/web/packages/magick/vignettes/intro.html
# date: 30/01/2022
#---

#install.packages("magick")
library(magick)

#install.packages("magrittr")
library(magrittr)  # to activate the pipe syntax

bc <- image_read(path = "../../Figures/Fig02_BC.jpg") %>%
  image_crop("3300x2900+1200+200") %>%
  image_annotate("(a)", size = 200, location = "+2900+100")

#image_browse(bc)

sprep <- image_read(path = "../../Figures/Fig02_SpRichness.jpg")%>%
  image_crop("3300x2900+1200+200") %>%
  image_annotate("(b)", size = 200, location = "+2900+100")

years <- image_read(path = "../../Figures/Fig02_#years.jpg") %>%
  image_crop("3300x2900+1200+200") %>%
  image_annotate("(c)", size = 200, location = "+2900+100")

rdi <- image_read(path = "../../Figures/Fig02_RDI.jpg") %>%
  image_crop("3300x2900+1200+200") %>%
  image_annotate("(d)", size = 200, location = "+2900+100")


gc()

top.row <- image_append(c(bc, sprep))
bottome.row <- image_append(c(years, rdi))

image_all <- image_append(c(top.row, bottome.row), stack = TRUE)
image_write(image_all, path = '../../Figures/Fig02_ALL_r.jpg', format = 'jpeg',
            quality = 75)
gc()

# priority network + selection frequency

priority.015 <- 
  image_read(path = "../../Figures/Fig03_PriorityNetwork_Target15.jpg") %>%
  image_crop("2500x2200+800+200") %>%
  image_annotate("(a)", size = 200, location = "+2200+100")

#image_browse(priority.015)

priority.025 <- 
  image_read(path = "../../Figures/Fig03_PriorityNetwork_Target25.jpg") %>%
  image_crop("2500x2200+800+200") %>%
  image_annotate("(c)", size = 200, location = "+2200+100")

priority.035 <- 
  image_read(path = "../../Figures/Fig03_PriorityNetwork_Target35.jpg") %>%
  image_crop("2500x2200+800+200") %>%
  image_annotate("(e)", size = 200, location = "+2200+100")

gc()

selecFreq.015 <- 
  image_read(path = "../../Figures/Fig03_SelecFreq_Target15.jpg") %>%
  image_crop("2500x2200+800+200") %>%
  image_annotate("(b)", size = 200, location = "+2200+100")

selecFreq.025 <- 
  image_read(path = "../../Figures/Fig03_SelecFreq_Target25.jpg") %>%
  image_crop("2500x2200+800+200") %>%
  image_annotate("(d)", size = 200, location = "+2200+100")

selecFreq.035 <- 
  image_read(path = "../../Figures/Fig03_SelecFreq_Target35.jpg") %>%
  image_crop("2500x2200+800+200") %>%
  image_annotate("(f)", size = 200, location = "+2200+100")

target_015 <- image_append(c(priority.015, selecFreq.015))
target_025 <- image_append(c(priority.025, selecFreq.025))
target_035 <- image_append(c(priority.035, selecFreq.035))

image_priority_all <- image_append(c(target_015, target_025, target_035), stack = TRUE)

image_write(image_priority_all, path = '../../Figures/Fig03_ALL_r.jpg', format = 'jpeg',
            quality = 50)

# priority network + selection frequency (no locked in areas)

priority.015 <- 
  image_read(path = "../../Figures/Fig03_PriorityNetwork_Target15_NoLockedIn.jpg") %>%
  image_crop("2500x2200+800+200") %>%
  image_annotate("(a)", size = 200, location = "+2200+100")

priority.025 <- 
  image_read(path = "../../Figures/Fig03_PriorityNetwork_Target25_NoLockedIn.jpg") %>%
  image_crop("2500x2200+800+200") %>%
  image_annotate("(c)", size = 200, location = "+2200+100")

priority.035 <- 
  image_read(path = "../../Figures/Fig03_PriorityNetwork_Target35_NoLockedIn.jpg") %>%
  image_crop("2500x2200+800+200") %>%
  image_annotate("(e)", size = 200, location = "+2200+100")

gc()

selecFreq.015 <- 
  image_read(path = "../../Figures/Fig03_SelecFreq_Target15_NoLockedIn.jpg") %>%
  image_crop("2500x2200+800+200") %>%
  image_annotate("(b)", size = 200, location = "+2200+100")

selecFreq.025 <- 
  image_read(path = "../../Figures/Fig03_SelecFreq_Target25_NoLockedIn.jpg") %>%
  image_crop("2500x2200+800+200") %>%
  image_annotate("(d)", size = 200, location = "+2200+100")

selecFreq.035 <- 
  image_read(path = "../../Figures/Fig03_SelecFreq_Target35_NoLockedIn.jpg") %>%
  image_crop("2500x2200+800+200") %>%
  image_annotate("(f)", size = 200, location = "+2200+100")

target_015 <- image_append(c(priority.015, selecFreq.015))
target_025 <- image_append(c(priority.025, selecFreq.025))
target_035 <- image_append(c(priority.035, selecFreq.035))

image_priority_all <- image_append(c(target_015, target_025, target_035), stack = TRUE)

image_write(image_priority_all, path = '../../Figures/Fig03_ALL_NoAreasLockedIn_r.jpg', format = 'jpeg',
            quality = 50)



