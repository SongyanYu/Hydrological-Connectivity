#---
# R tutorial for advanced image processing with package magick.
# tutoiral website: https://cran.r-project.org/web/packages/magick/vignettes/intro.html
# date: 30/01/2022
#---

#install.packages("magick")
library(magick)

#install.packages("magrittr")
library(magrittr)  # to activate the pipe syntax

bc <- image_read(path = "../../Figures/Fig03") %>%
  image_crop("7000x6000+2200+500") %>%
  image_annotate("(a)", size = 300, location = "+6000+100")

sprep <- image_read(path = "fig/Fig02_SpRichness.jpg")%>%
  image_crop("7000x6000+2200+500") %>%
  image_annotate("(b)", size = 300, location = "+6000+100")

years <- image_read(path = 'fig/Fig02_#years.jpg') %>%
  image_crop("7000x6000+2200+500") %>%
  image_annotate("(c)", size = 300, location = "+6000+100")

rdi <- image_read(path = "fig/Fig02_RDI.jpg") %>%
  image_crop("7000x6000+2200+500") %>%
  image_annotate("(d)", size = 300, location = "+6000+100")

gc()

top.row <- image_append(c(bc, sprep))
bottome.row <- image_append(c(years, rdi))

image_all <- image_append(c(top.row, bottome.row), stack = TRUE)
image_write(image_all, path = 'fig/Fig02_ALL_r.jpg', format = 'jpeg',
            quality = 100)
gc()
