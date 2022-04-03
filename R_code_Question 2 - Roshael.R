library(tidyverse)
library(ggplot2)

getwd()
MyData <- read.csv("./Data/GreenhouseData.csv")

# subset the dataset to make smaller dataset with PlantSpecies, SoilHostSpecies, Age, and Total_Dry_mg
SubsetMyData <- MyData %>% select(c("PlantSpecies", "SoilHostSpecies", "Age", "Total_Dry_mg"))

# subset further for each plant species
AaMDSubset <- SubsetMyData %>% filter(PlantSpecies == "Aa")
BpMDSubset <- SubsetMyData %>% filter(PlantSpecies == "Bp")
CeMDSubset <- SubsetMyData %>% filter(PlantSpecies == "Ce")
LaMDSubset <- SubsetMyData %>% filter(PlantSpecies == "La")

# create stacked bar graphs to show the variation in soil host species composition as age increases per plant species
ggplot(data = AaMDSubset, aes(x = Age, fill = SoilHostSpecies)) +
  geom_bar() +
  scale_x_continuous(breaks = seq(0, 28, by = 2)) + 
  scale_y_continuous(breaks = seq(0, 20, by = 2)) +
  ggtitle("Aa Soil Host Species Composition") +
  xlab("Age (days)") + ylab("Number of Samples") +
  guides(fill=guide_legend(title="Soil Host Species")) +
  theme_bw()

ggplot(data = BpMDSubset, aes(x = Age, fill = SoilHostSpecies)) +
  geom_bar() +
  scale_x_continuous(breaks = seq(0, 28, by = 2)) +  
  scale_y_continuous(breaks = seq(0, 20, by = 2)) +
  ggtitle("Bp Soil Host Species Composition") +
  xlab("Age (days)") + ylab("Number of Samples") +
  guides(fill=guide_legend(title="Soil Host Species")) +
  theme_bw()

ggplot(data = CeMDSubset, aes(x = Age, fill = SoilHostSpecies)) +
  geom_bar() +
  scale_x_continuous(breaks = seq(0, 28, by = 2)) +  
  scale_y_continuous(breaks = seq(0, 20, by = 2)) +
  ggtitle("Ce Soil Host Species Composition") +
  xlab("Age (days)") + ylab("Number of Samples") +
  guides(fill=guide_legend(title="Soil Host Species")) +
  theme_bw()

ggplot(data = LaMDSubset, aes(x = Age, fill = SoilHostSpecies)) +
  geom_bar() +
  scale_x_continuous(breaks = seq(0, 28, by = 2)) +  
  scale_y_continuous(breaks = seq(0, 20, by = 2)) +
  ggtitle("La Soil Host Species Composition") +
  xlab("Age (days)") + ylab("Number of Samples") +
  guides(fill=guide_legend(title="Soil Host Species")) +
  theme_bw()
