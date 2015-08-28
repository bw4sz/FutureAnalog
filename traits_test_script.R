require(dplyr)

# load amnat traits list
amnat_traits <- read.csv("InputData/hummingbird_traits.csv")
amnat_traits$species <- gsub(" ", ".", amnat_traits$species)

# load original traits list
morph <- read.csv("InputData/MorphologyShort.csv", na.strings="9999")

#just get males & the 3 traits of interest
mon <- filter(morph, Sex == "Macho") %>%
  select(SpID, ExpC, Peso, AlCdo) %>%
  group_by(SpID) %>%
  summarise_each(funs(mean(., na.rm = TRUE))) %>%
  filter(complete.cases(.))

mon <- data.frame(mon)
colnames(mon) <- c("species","Bill","Mass","WingChord")
mon$species <- gsub(" ", ".", mon$species)

# load modelled species
load("sppXsite/current.rda")


w <- apply(sppXsite, 2, function(x)all(!is.na(x)))
if (any(w)) {
  modelled.species <- names(which(w))
}

modelled.species <- data.frame(species=modelled.species[2:92]) # get rid of additional columns

test <- merge(modelled.species, mon, all.x=TRUE)
test <- merge(test, amnat_traits, all.x=TRUE)

