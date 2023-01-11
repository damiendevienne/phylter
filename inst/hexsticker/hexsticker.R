#' 
#' Create an Hexagonal Sticker for the Package
#' 
library(hexSticker) # hexSticker generator
library(magick)     # Advanced image processing

sticker(
  subplot = magick::image_read(here::here("inst", "hexsticker", "icon.png")),
  package = "phylter",
  s_width = 1.7,              # subplot width
  s_height = 1.7,             # subplot height
  s_x = 1,                    # subplot x-position (left/right position)
  s_y = 1,                    # subplot y-position (up/down position)
  p_size = 0,                 # package name font size
  h_fill = 'white',           # hexSticker background color
  h_color = 'orange',         # hexsticker border color
  h_size = 1.5,               # hexSticker size
  url = "https://damiendevienne.github.io/phylter/",
  u_size = 3,                 # url font size
  u_color = 'black',          # url font color
  filename = here::here("man", "figures", "hexsticker.png")
)

