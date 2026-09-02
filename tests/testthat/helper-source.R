

for (f in list.files("R/R", pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}