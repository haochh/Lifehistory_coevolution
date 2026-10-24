#----------------------------------------------------------------------------#
#                                                                            #
#        Zhang et al. Code for generating figures in the main text           #
#                                                                            #
#----------------------------------------------------------------------------#



##### Figure 2: Data overview -----               


#----------------------------------------------------------#
#  Description: This figure is to display the overview of  #
#  the life history trait value across species for mammals,#
#  birds and reptiles, respectively.                       #
#----------------------------------------------------------#


# Load packages

library(tidyverse)
library(ape)
library(ggtree)
library(treeio)
library(ggtreeExtra)
library(paletteer)
library(rphylopic)
library(ggimage)
library(ggnewscale)
library(deeptime) 
library(patchwork)

#---- ---- ---- ---- ---- ---- ----# 
#              Mammal              #
#---- ---- ---- ---- ---- ---- ----# 
# Data is not scaled
mammal_data <- read.csv("all_mammal_data.csv")
# Tree is not scaled
mammal_tree <- read.tree("last_mammal_tree.tre")

##### Layer 1: Tree + Geologic timescale -----
# Calculate maximum tree height (root age)
max(node.depth.edgelength(mammal_tree)) # = 179.2078

p <- ggtree(mammal_tree, linewidth = 0.3, color = "gray30") +
  coord_geo(xlim = c(-180, 0), 
            ylim = c(-5, Ntip(mammal_tree)), 
            neg = TRUE, 
            abbrv = F,  
            pos = "bottom") +
  scale_x_continuous(breaks = seq(-180, 0, 40), 
                     labels = abs(seq(-180, 0, 40)),
                     name = "Time (Mya)") +
  theme_tree2()

p2 <- revts(p)


##### Layer 2: Major Order visualization -----
head(mammal_data)
### Step 1: Get the major Orders that we would like to visualize
# The major Orders are n >= 50, n-Order = 7
order.m <- mammal_data %>%
  group_by(Order) %>%
  summarise(n = n()) %>%
  filter(n >= 50)

### Step 2: Get the node label of them
clade_info.m <- mammal_data %>%
  filter(Order %in% order.m$Order) %>%
  group_by(Order) %>%
  summarize(Species_list = list(Species)) %>%
  mutate(
    node = sapply(Species_list, function(sp) {
      if (length(sp) == 1) {
        # Only one species: use its tip index
        which(mammal_tree$tip.label == sp)
      } else {
        # Multiple species: get the MRCA node
        getMRCA(mammal_tree, tip = sp)
      }
    }),
    group_type = sapply(Species_list, function(sp) {
      if (length(sp) == 1) "tip" else "node"
    })
  ) %>%
  ungroup() 

head(clade_info.m)

p1.2 <- ggtree(mammal_tree, linewidth = 0.3, color = "gray30")
### Step 3: Add color to each Order
p3 <- p1.2 +
  geom_hilight(node=785, fill="#8DD3C7FF", align="both", alpha=.3)+
  geom_hilight(node=877, fill="#FFFFB3FF", align="both", alpha=.3)+
  geom_hilight(node=937, fill="#BEBADAFF", align="both", alpha=.3)+
  geom_hilight(node=1386, fill="#EAC862", align="both", alpha=.3)+
  geom_hilight(node=724, fill="#80B1D3FF", align="both", alpha=.3)+
  geom_hilight(node=998, fill="#BC80BDFF", align="both", alpha=.3)+
  geom_hilight(node=1116, fill="#FF7F00", align="both", alpha=.3)


clade_info.m <- clade_info.m %>%
  mutate(color = case_when(
    node == 785 ~ "#8DD3C7FF",
    node == 877 ~ "#FFFFB3FF",
    node == 937 ~ "#BEBADAFF",
    node == 1386 ~ "#EAC862",
    node == 724 ~ "#80B1D3FF",
    node == 998 ~ "#BC80BDFF",
    node == 1116 ~ "#FF7F00"
  ))

### Step 4: Add Order images
# Retrieve the image id:Phylopic UID (https://www.phylopic.org/)
uid.m <- data.frame()
for (i in seq_along(order.m$Order)) {
  order <- order.m$Order[i]
  id <- rphylopic::get_uuid(order)
  df <- data.frame(Order = order,
                   UID = id)
  uid.m <- rbind(uid.m, df)
}

clade_info.m <- left_join(clade_info.m, uid.m)

p4 <- p3 +
  geom_cladelab(
    data = filter(clade_info.m, group_type == "node"),
    mapping = aes(
      node = node, 
      label = Order, 
      image = UID
    ),
    geom = "phylopic",
    imagecolor = "grey20",
    barsize = 0,
    imagesize = 0.1,
    barcolor = "transparent",
    align = TRUE,
    show.legend = FALSE,
    offset = -100
  ) + ggtree::vexpand(0.05, -1)

p4

##### Layer 3: Life History data -----
# Load scaled data
load("mammal_data_ordered.RData")
head(mammal_data_ordered)

mammal_data_ordered <- mammal_data_ordered %>%
  rownames_to_column(var = "Species")

# DT
pDT <- p4 +
  geom_fruit(
    data = subset(mammal_data, Development.time >= 0 & Development.time <= 3000), # Filter data here: (0, 3000)
    geom = geom_bar,
    mapping = aes(x = Development.time, y = Species, fill = Development.time),
    offset = 0.1,
    pwidth = 0.3,
    stat = "identity",
    grid.params = list(
      size = 0.1,
      alpha = 0.7
    ),
    axis.params = list(
      axis = "x",
      text.size = 4.5,
      hjust = 0.8, 
      vjust = 1.5
    )
  ) +
  scale_fill_paletteer_c(
    palette = "grDevices::Purple-Yellow",
    direction = -1,
    limits = c(0, 3000),
    name = "Development time",
    guide = guide_colorbar( 
      barheight = unit(4, "cm"), 
      barwidth = unit(1, "cm"), 
      title.position = "top",
      title.hjust = 0.5
    )
  ) + ggtree::vexpand(0.05, -1) 


# OS
pOS <- pDT +
  new_scale_fill() +
  geom_fruit(
    data = subset(mammal_data, Offspring.size  >= 0 & Offspring.size  <= 5000), # Filter data here: (0, 5000)
    geom = geom_bar,
    mapping = aes(x = Offspring.size, y = Species, fill = Offspring.size),
    offset = 0.05,
    pwidth = 0.3,
    stat = "identity",
    grid.params = list(
      size = 0.1,
      alpha = 0.7
    ),
    axis.params = list(
      axis = "x",
      text.size = 4.5,
      hjust = 0.8, 
      vjust = 1.5
    )
  ) +
  scale_fill_paletteer_c(
    palette = "grDevices::SunsetDark",
    direction = -1,
    limits = c(0, 5000),
    name = "Offspring size",
    guide = guide_colorbar( 
      barheight = unit(4, "cm"), 
      barwidth = unit(1, "cm"), 
      title.position = "top",
      title.hjust = 0.5
    )
  ) 

# AL
head(mammal_data)
hist(mammal_data$Adult.lifespan, breaks = 50)

pAL <- pOS +
  new_scale_fill() +
  geom_fruit(
    data = subset(mammal_data, Adult.lifespan  >= 0 & Adult.lifespan  <= 20000), # Filter data here: (0, 20000)
    geom = geom_bar,
    mapping = aes(x = Adult.lifespan, y = Species, fill = Adult.lifespan),
    offset = 0.05,
    pwidth = 0.3,
    stat = "identity",
    grid.params = list(
      size = 0.1,
      alpha = 0.7
    ),
    axis.params = list(
      axis = "x",
      text.size = 4.5,
      hjust = 0.8, 
      vjust = 1.5
    )
  ) +
  scale_fill_paletteer_c(
    palette = "grDevices::Blues 3",
    direction = -1,
    limits = c(0, 20000),
    name = "Adult lifespan",
    guide = guide_colorbar( 
      barheight = unit(4, "cm"), 
      barwidth = unit(1, "cm"), 
      title.position = "top",
      title.hjust = 0.5
    )
  ) 

# CS
head(mammal_data)
hist(mammal_data$Clutch.size, breaks = 50)

pCS <- pAL +
  new_scale_fill() +
  geom_fruit(
    data = mammal_data,
    #data = subset(mammal_data, Clutch.size  >= 0 & Clutch.size  <= 20000), # Filter data here: (0, 20000)
    geom = geom_bar,
    mapping = aes(x = Clutch.size, y = Species, fill = Clutch.size),
    offset = 0.05,
    pwidth = 0.3,
    stat = "identity",
    grid.params = list(
      size = 0.1,
      alpha = 0.7
    ),
    axis.params = list(
      axis = "x",
      text.size = 4.5,
      hjust = 0.8, 
      vjust = 1.5
    )
  ) +
  scale_fill_paletteer_c(
    palette = "grDevices::Red-Purple",
    direction = -1,
    #limits = c(0, 20000),
    name = "Clutch size",
    guide = guide_colorbar( 
      barheight = unit(4, "cm"), 
      barwidth = unit(1, "cm"), 
      title.position = "top",
      title.hjust = 0.5
    )
  ) 


# CS
head(mammal_data)
hist(mammal_data$Clutch.frequency, breaks = 50)

pCF <- pCS +
  new_scale_fill() +
  geom_fruit(
    data = mammal_data,
    #data = subset(mammal_data, Clutch.size  >= 0 & Clutch.size  <= 20000), # Filter data here: (0, 20000)
    geom = geom_bar,
    mapping = aes(x = Clutch.frequency, y = Species, fill = Clutch.frequency),
    offset = 0.05,
    pwidth = 0.3,
    stat = "identity",
    grid.params = list(
      size = 0.1,
      alpha = 0.7
    ),
    axis.params = list(
      axis = "x",
      text.size = 4.5,
      hjust = 0.8, 
      vjust = 1.5
    )
  ) +
  scale_fill_paletteer_c(
    palette = "grDevices::Oranges",
    direction = -1,
    #limits = c(0, 20000),
    name = "Clutch frequency",
    guide = guide_colorbar( 
      barheight = unit(4, "cm"), 
      barwidth = unit(1, "cm"), 
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  theme(legend.position = "none")

pMammal <- pCF
pdf(file = "/Users/user/Chunhui/Zhang_mammal/revision/Visualization/fig1_test.pdf", height = 14.36, width = 22)
pMammal
dev.off()

#---- ---- ---- ---- ---- ---- ----# 
#               Bird               #
#---- ---- ---- ---- ---- ---- ----# 
# Data is not scaled
bird_data <- read.csv("all_bird_data.csv")
# Tree is not scaled
bird_tree <- read.tree("last_bird_tree.tre")

##### Layer 1: Tree + Geologic timescale -----
# Calculate maximum tree height (root age)
max(node.depth.edgelength(bird_tree)) # = 101.1764

p <- ggtree(mammal_tree, linewidth = 0.3, color = "gray30") +
  coord_geo(xlim = c(-180, 0), 
            ylim = c(-5, Ntip(mammal_tree)), 
            neg = TRUE, 
            abbrv = F,  
            pos = "bottom") +
  scale_x_continuous(breaks = seq(-180, 0, 40), 
                     labels = abs(seq(-180, 0, 40)),
                     name = "Time (Mya)") +
  theme_tree2()

p2 <- revts(p)

##### Layer 2: Major Order visualization -----
head(bird_data)
### Step 1: Get the major Orders that we would like to visualize
# The major Orders are n >= 20, n-Order = 2
order.b <- bird_data %>%
  group_by(Order) %>%
  summarise(n = n()) %>%
  filter(n >= 20)

### Step 2: Get the node label of them
clade_info.b <- bird_data %>%
  filter(Order %in% order.b$Order) %>%
  group_by(Order) %>%
  summarize(Species_list = list(Species)) %>%
  mutate(
    node = sapply(Species_list, function(sp) {
      if (length(sp) == 1) {
        # Only one species: use its tip index
        which(bird_tree$tip.label == sp)
      } else {
        # Multiple species: get the MRCA node
        getMRCA(bird_tree, tip = sp)
      }
    }),
    group_type = sapply(Species_list, function(sp) {
      if (length(sp) == 1) "tip" else "node"
    })
  ) %>%
  ungroup() 

head(clade_info.b)

p1.2.b <- ggtree(bird_tree, linewidth = 0.3, color = "gray30")
### Step 3: Add color to each Order
p3.b <- p1.2.b +
  geom_hilight(node=151, fill="#FFED6FFF", align="both", alpha=.3)+
  geom_hilight(node=171, fill="#BC80BDFF", align="both", alpha=.3)


### Step 4: Add Order images
# Retrieve the image id:Phylopic UID (https://www.phylopic.org/)
uid.b <- data.frame()
for (i in seq_along(order.b$Order)) {
  order <- order.b$Order[i]
  id <- rphylopic::get_uuid(order)
  df <- data.frame(Order = order,
                   UID = id)
  uid.b <- rbind(uid.b, df)
}

clade_info.b <- left_join(clade_info.b, uid.b)

p4.b <- p3.b +
  geom_cladelab(
    data = filter(clade_info.b, group_type == "node"),
    mapping = aes(
      node = node, 
      label = Order, 
      image = UID
    ),
    geom = "phylopic",
    imagecolor = "grey20",
    barsize = 0,
    imagesize = 0.5,
    barcolor = "transparent",
    align = TRUE,
    show.legend = FALSE,
    offset = -100
  ) + ggtree::vexpand(0.05, -1)

p4.b



##### Layer 3: Life History data -----
# DT
head(bird_data)
hist(bird_data$DT, breaks = 50)

pDT.b <- p4.b +
  geom_fruit(
    data = subset(bird_data, DT >= 0 & DT <= 3000), # Filter data here: (0, 3000)
    geom = geom_bar,
    mapping = aes(x = DT, y = Species, fill = DT),
    offset = 0.1,
    pwidth = 0.3,
    stat = "identity",
    grid.params = list(
      size = 0.1,
      alpha = 0.7
    ),
    axis.params = list(
      axis = "x",
      text.size = 4.5,
      hjust = 0.8, 
      vjust = 1.5
    )
  ) +
  scale_fill_paletteer_c(
    palette = "grDevices::Purple-Yellow",
    direction = -1,
    limits = c(0, 3000),
    name = "Development time",
    guide = guide_colorbar( 
      barheight = unit(4, "cm"), 
      barwidth = unit(1, "cm"), 
      title.position = "top",
      title.hjust = 0.5
    )
  ) + ggtree::vexpand(0.05, -1) 

pDT.b

# OS
head(bird_data)
hist(bird_data$OS, breaks = 50)

pOS.b <- pDT.b +
  new_scale_fill() +
  geom_fruit(
    data = subset(bird_data, OS  >= 0 & OS  <= 1000), # Filter data here: (0, 1000)
    geom = geom_bar,
    mapping = aes(x = OS, y = Species, fill = OS),
    offset = 0.05,
    pwidth = 0.3,
    stat = "identity",
    grid.params = list(
      size = 0.1,
      alpha = 0.7
    ),
    axis.params = list(
      axis = "x",
      text.size = 4.5,
      hjust = 0.8, 
      vjust = 1.5
    )
  ) +
  scale_fill_paletteer_c(
    palette = "grDevices::SunsetDark",
    direction = -1,
    limits = c(0, 1000),
    name = "Offspring size",
    guide = guide_colorbar( 
      barheight = unit(4, "cm"), 
      barwidth = unit(1, "cm"), 
      title.position = "top",
      title.hjust = 0.5
    )
  ) 

pOS.b

# AL
head(bird_data)
hist(bird_data$AL, breaks = 50)

pAL.b <- pOS.b +
  new_scale_fill() +
  geom_fruit(
    data = bird_data,
    geom = geom_bar,
    mapping = aes(x = AL, y = Species, fill = AL),
    offset = 0.05,
    pwidth = 0.3,
    stat = "identity",
    grid.params = list(
      size = 0.1,
      alpha = 0.7
    ),
    axis.params = list(
      axis = "x",
      text.size = 4.5,
      hjust = 0.8, 
      vjust = 1.5
    )
  ) +
  scale_fill_paletteer_c(
    palette = "grDevices::Blues 3",
    direction = -1,
    #limits = c(0, 20000),
    name = "Adult lifespan",
    guide = guide_colorbar( 
      barheight = unit(4, "cm"), 
      barwidth = unit(1, "cm"), 
      title.position = "top",
      title.hjust = 0.5
    )
  ) 

# CS
head(bird_data)
hist(bird_data$CS, breaks = 50)

pCS.b <- pAL.b +
  new_scale_fill() +
  geom_fruit(
    data = bird_data,
    #data = subset(mammal_data, Clutch.size  >= 0 & Clutch.size  <= 20000), # Filter data here: (0, 20000)
    geom = geom_bar,
    mapping = aes(x = CS, y = Species, fill = CS),
    offset = 0.05,
    pwidth = 0.3,
    stat = "identity",
    grid.params = list(
      size = 0.1,
      alpha = 0.7
    ),
    axis.params = list(
      axis = "x",
      text.size = 4.5,
      hjust = 0.8, 
      vjust = 1.5
    )
  ) +
  scale_fill_paletteer_c(
    palette = "grDevices::Red-Purple",
    direction = -1,
    #limits = c(0, 20000),
    name = "Clutch size",
    guide = guide_colorbar( 
      barheight = unit(4, "cm"), 
      barwidth = unit(1, "cm"), 
      title.position = "top",
      title.hjust = 0.5
    )
  ) 


# CF
head(bird_data)
hist(bird_data$CF, breaks = 50)

pCF.b <- pCS.b +
  new_scale_fill() +
  geom_fruit(
    data = bird_data,
    #data = subset(mammal_data, Clutch.size  >= 0 & Clutch.size  <= 20000), # Filter data here: (0, 20000)
    geom = geom_bar,
    mapping = aes(x = CF, y = Species, fill = CF),
    offset = 0.05,
    pwidth = 0.3,
    stat = "identity",
    grid.params = list(
      size = 0.1,
      alpha = 0.7
    ),
    axis.params = list(
      axis = "x",
      text.size = 4.5,
      hjust = 0.8, 
      vjust = 1.5
    )
  ) +
  scale_fill_paletteer_c(
    palette = "grDevices::Oranges",
    direction = -1,
    #limits = c(0, 20000),
    name = "Clutch frequency",
    guide = guide_colorbar( 
      barheight = unit(4, "cm"), 
      barwidth = unit(1, "cm"), 
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  theme(legend.position = "none")

pBird <- pCF.b

p.final <- (pBird/pMammal)+
  plot_layout(heights = c(1.44, 7.18))

pdf(file = "/Users/user/Chunhui/Zhang_mammal/revision/Visualization/fig1_test.pdf", height = 17.24, width = 23)
p.final
dev.off()

#---- ---- ---- ---- ---- ---- ----# 
#              Reptile             #
#---- ---- ---- ---- ---- ---- ----# 
# Data is not scaled
reptile_data <- read.csv("all_reptile_data.csv")
# Tree is not scaled
reptile_tree <- read.tree("last_reptile_tree.tre")

##### Layer 1: Tree + Geologic timescale -----
# Calculate maximum tree height (root age)
max(node.depth.edgelength(reptile_tree)) # = 178.1413

p <- ggtree(mammal_tree, linewidth = 0.3, color = "gray30") +
  coord_geo(xlim = c(-180, 0), 
            ylim = c(-5, Ntip(mammal_tree)), 
            neg = TRUE, 
            abbrv = F, 
            pos = "bottom") +
  scale_x_continuous(breaks = seq(-180, 0, 40), 
                     labels = abs(seq(-180, 0, 40)),
                     name = "Time (Mya)") +
  theme_tree2()

p2 <- revts(p)

##### Layer 2: Major Families visualization -----
head(reptile_data)
### Step 1: Get the major Families that we would like to visualize
# The major Orders are n >= 25, n-Family = 7
family.r <- reptile_data %>%
  group_by(Family) %>%
  summarise(n = n()) %>%
  filter(n >= 25)

### Step 2: Get the node label of them
clade_info.r <- reptile_data %>%
  filter(Family %in% family.r$Family) %>%
  group_by(Family) %>%
  summarize(Species_list = list(Species)) %>%
  mutate(
    node = sapply(Species_list, function(sp) {
      if (length(sp) == 1) {
        # Only one species: use its tip index
        which(reptile_tree$tip.label == sp)
      } else {
        # Multiple species: get the MRCA node
        getMRCA(reptile_tree, tip = sp)
      }
    }),
    group_type = sapply(Species_list, function(sp) {
      if (length(sp) == 1) "tip" else "node"
    })
  ) %>%
  ungroup() 

head(clade_info.r)

p1.2.r <- ggtree(reptile_tree, linewidth = 0.3, color = "gray30")
### Step 3: Add color to each Order
p3.r <- p1.2.r +
  geom_hilight(node=473, fill="#8DD3C7FF", align="both", alpha=.5)+
  geom_hilight(node=642, fill="#FFFFB3FF", align="both", alpha=.5)+
  geom_hilight(node=772, fill="#BEBADAFF", align="both", alpha=.5)+
  geom_hilight(node=702, fill="#EAC862", align="both", alpha=.5)+
  geom_hilight(node=558, fill="#80B1D3FF", align="both", alpha=.5)+
  geom_hilight(node=415, fill="#E31A1C", align="both", alpha=.5)+
  geom_hilight(node=607, fill="#FF7F00", align="both", alpha=.5)

### Step 4: Add Order images
# Retrieve the image id:Phylopic UID (https://www.phylopic.org/)
uid.r <- data.frame()
for (i in seq_along(family.r$Family)) {
  family <- family.r$Family[i]
  id <- rphylopic::get_uuid(family)
  df <- data.frame(Family = family,
                   UID = id)
  uid.r <- rbind(uid.r, df)
}

clade_info.r <- left_join(clade_info.r, uid.r)

p4.r <- p3.r +
  geom_cladelab(
    data = filter(clade_info.r, group_type == "node"),
    mapping = aes(
      node = node, 
      label = Family, 
      image = UID
    ),
    geom = "phylopic",
    imagecolor = "grey20",
    barsize = 0,
    imagesize = 0.12,
    barcolor = "transparent",
    align = TRUE,
    show.legend = FALSE,
    offset = -100
  ) + ggtree::vexpand(0.05, -1)

##### Layer 3: Life History data -----

# DT
head(reptile_data)
hist(reptile_data$DT, breaks = 50)

pDT.r <- p4.r +
  geom_fruit(
    data = subset(reptile_data, DT >= 0 & DT <= 2000), # Filter data here: (0, 2000)
    geom = geom_bar,
    mapping = aes(x = DT, y = Species, fill = DT),
    offset = 0.1,
    pwidth = 0.3,
    stat = "identity",
    grid.params = list(
      size = 0.1,
      alpha = 0.7
    ),
    axis.params = list(
      axis = "x",
      text.size = 4.5,
      hjust = 0.8, 
      vjust = 1.5
    )
  ) +
  scale_fill_paletteer_c(
    palette = "grDevices::Purple-Yellow",
    direction = -1,
    limits = c(0, 2000),
    name = "Development time",
    guide = guide_colorbar( 
      barheight = unit(4, "cm"), 
      barwidth = unit(1, "cm"), 
      title.position = "top",
      title.hjust = 0.5
    )
  ) + ggtree::vexpand(0.05, -1) 

pDT.r

# OS
head(reptile_data)
hist(reptile_data$OS, breaks = 50)

pOS.r <- pDT.r +
  new_scale_fill() +
  geom_fruit(
    data = subset(reptile_data, OS  >= 0 & OS  <= 50), # Filter data here: (0, 50)
    geom = geom_bar,
    mapping = aes(x = OS, y = Species, fill = OS),
    offset = 0.05,
    pwidth = 0.3,
    stat = "identity",
    grid.params = list(
      size = 0.1,
      alpha = 0.7
    ),
    axis.params = list(
      axis = "x",
      text.size = 4.5,
      hjust = 0.8, 
      vjust = 1.5
    )
  ) +
  scale_fill_paletteer_c(
    palette = "grDevices::SunsetDark",
    direction = -1,
    limits = c(0, 50),
    name = "Offspring size",
    guide = guide_colorbar( 
      barheight = unit(4, "cm"), 
      barwidth = unit(1, "cm"), 
      title.position = "top",
      title.hjust = 0.5
    )
  ) 

pOS.r

# AL
head(reptile_data)
hist(reptile_data$AL, breaks = 50)

pAL.r <- pOS.r +
  new_scale_fill() +
  geom_fruit(
    data = reptile_data,
    geom = geom_bar,
    mapping = aes(x = AL, y = Species, fill = AL),
    offset = 0.05,
    pwidth = 0.3,
    stat = "identity",
    grid.params = list(
      size = 0.1,
      alpha = 0.7
    ),
    axis.params = list(
      axis = "x",
      text.size = 4.5,
      hjust = 0.8, 
      vjust = 1.5
    )
  ) +
  scale_fill_paletteer_c(
    palette = "grDevices::Blues 3",
    direction = -1,
    #limits = c(0, 20000),
    name = "Adult lifespan",
    guide = guide_colorbar( 
      barheight = unit(4, "cm"), 
      barwidth = unit(1, "cm"), 
      title.position = "top",
      title.hjust = 0.5
    )
  ) 

pAL.r

# CS
head(reptile_data)
hist(reptile_data$CS, breaks = 50)

pCS.r <- pAL.r +
  new_scale_fill() +
  geom_fruit(
    data = reptile_data,
    #data = subset(mammal_data, Clutch.size  >= 0 & Clutch.size  <= 20000), # Filter data here: (0, 20000)
    geom = geom_bar,
    mapping = aes(x = CS, y = Species, fill = CS),
    offset = 0.05,
    pwidth = 0.3,
    stat = "identity",
    grid.params = list(
      size = 0.1,
      alpha = 0.7
    ),
    axis.params = list(
      axis = "x",
      text.size = 4.5,
      hjust = 0.8, 
      vjust = 1.5
    )
  ) +
  scale_fill_paletteer_c(
    palette = "grDevices::Red-Purple",
    direction = -1,
    #limits = c(0, 20000),
    name = "Clutch size",
    guide = guide_colorbar( 
      barheight = unit(4, "cm"), 
      barwidth = unit(1, "cm"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) 

pCS.r

# CF
head(reptile_data)
hist(reptile_data$CF, breaks = 50)

pCF.r <- pCS.r +
  new_scale_fill() +
  geom_fruit(
    data = reptile_data,
    #data = subset(mammal_data, Clutch.size  >= 0 & Clutch.size  <= 20000), # Filter data here: (0, 20000)
    geom = geom_bar,
    mapping = aes(x = CF, y = Species, fill = CF),
    offset = 0.05,
    pwidth = 0.3,
    stat = "identity",
    grid.params = list(
      size = 0.1,
      alpha = 0.7
    ),
    axis.params = list(
      axis = "x",
      text.size = 4.5,
      hjust = 0.8, 
      vjust = 1.5
    )
  ) +
  scale_fill_paletteer_c(
    palette = "grDevices::Oranges",
    direction = -1,
    #limits = c(0, 20000),
    name = "Clutch frequency",
    guide = guide_colorbar( 
      barheight = unit(4, "cm"), 
      barwidth = unit(1, "cm"), 
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  theme(legend.position = "none")

pCF.r

pReptile <- pCF.r

p.final <- (pMammal/pBird/pReptile)+
  plot_layout(heights = c(7.18,1.44,4.11))

p.final <- (pBird/pReptile/pMammal)+
  plot_layout(heights = c(1.44,4.11,7.18))


##### Get the final version of figure 1
pdf(file = "Figure1.pdf", height = 25.46, width = 27)
p.final
dev.off()


##### Figure 4:  PPA output -----                     


## Load packages

library(ggplot2)
library(dplyr)
library(patchwork)
library(grid)

load("mammal_PPA_results.RData")
load("bird_PPA_results.RData")
load("reptile_PPA_results.RData")
ls()

average_model_mammals <- average(PPA_mammals, avg_method = "full")
average_model_birds <- average(PPA_birds, avg_method = "full")
average_model_reptiles <- average(PPA_reptiles, avg_method = "full")

manual_layout_data <- data.frame(
  name = c("OS", "CS", "DT", "CF", "AL"),  
  x = c(8, 20, 8, 20, 14),               
  y = c(9, 9, 7, 7, 5.5)                
)


# mammals

PPA_mammal <- plot(
  average_model_mammals, 
  type = "width", 
  labels = NULL, 
  algorithm = "sugiyama",
  manual_layout = manual_layout_data,
  text_size = 10,
  box_x = 18,
  box_y = 14,
  edge_width = 20,
  curvature = 0,
  rotation = 0,
  flip_x = FALSE,
  flip_y = FALSE,
  arrow = grid::arrow(type = "closed", 18, grid::unit(15, "points")),
  colors = c("#2E5A87FF", "#A90C38FF"),
  show.legend = TRUE
) + theme(legend.position = "none") + 
  theme(aspect.ratio = 0.7)   

print(PPA_mammal)

# birds

PPA_bird <- plot(
  average_model_birds, 
  type = "width", 
  labels = NULL, 
  algorithm = "sugiyama",
  manual_layout = manual_layout_data,
  text_size = 10,
  box_x = 18,
  box_y = 14,
  edge_width = 20,
  curvature = 0,
  rotation = 0,
  flip_x = FALSE,
  flip_y = FALSE,
  arrow = grid::arrow(type = "closed", 18, grid::unit(15, "points")),
  colors = c("#2E5A87FF", "#A90C38FF"),
  show.legend = TRUE
) + theme(legend.position = "none") + 
  theme(aspect.ratio = 0.7)   


# reptiles

PPA_reptile <- plot(
  average_model_reptiles, 
  type = "width", 
  labels = NULL, 
  algorithm = "sugiyama",
  manual_layout = manual_layout_data,
  text_size = 10,
  box_x = 18,
  box_y = 14,
  edge_width = 20,
  curvature = 0,
  rotation = 0,
  flip_x = FALSE,
  flip_y = FALSE,
  arrow = grid::arrow(type = "closed", 18, grid::unit(15, "points")),
  colors = c("#2E5A87FF", "#A90C38FF"),
  show.legend = TRUE
) + theme(legend.position = "none") + 
  theme(aspect.ratio = 0.7)   


#----------------------------------------------------------#
#  Construct the combined dataset using data from Table S6 #
#----------------------------------------------------------#


combined_data$Class <- factor(combined_data$Class,
                              levels = c("Mammals", "Birds", "Reptiles"))

combined_data <- combined_data %>%
  mutate(
    Predictor = factor(Predictor, levels = rev(unique(Predictor))),  
    id = as.numeric(Predictor),
    bg_fill = ifelse(id %% 2 == 0, "gray95", "white"),
    color_by_CI = case_when(
      lower > 0 & upper > 0 ~ "#A90C38FF",    
      lower < 0 & upper < 0 ~ "#2E5A87FF",    
      lower < 0 & upper > 0 ~ "grey70",      
      TRUE ~ "grey70"
    )
  )


forest_all <- ggplot(combined_data, aes(x = Predictor, y = Estimate, 
                                        shape = Model, color = color_by_CI)) +
  geom_rect(data = distinct(combined_data, Predictor, id, bg_fill),
            aes(xmin = id - 0.5, xmax = id + 0.5,
                ymin = -Inf, ymax = Inf,
                fill = bg_fill),
            color = NA,
            alpha = 0.4,
            inherit.aes = FALSE) +
  scale_fill_identity() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.6),
                  size = 1.1, linewidth = 1) +
  coord_flip() +
  scale_color_identity() +
  scale_shape_manual(values = c("Multivariate" = 16, "Univariate" = 17)) +
  facet_wrap(~ Class, nrow = 1, strip.position = "top") +
  labs(x = "Correlation",
       y = "Coefficient (β) with 95% CI",
       shape = "Model") +
  theme_bw(base_size = 15) +
  theme(
    strip.text = element_text(size = 30, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),  
    axis.text.y = element_text(size = 22),
    axis.text.x = element_text(size = 22),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  ) 

PPA_all <- PPA_mammal | PPA_bird | PPA_reptile 

fig3 <- PPA_all / forest_all


##### Get the final version of figure 3
pdf(file = "Figure3.pdf", height = 25.46, width = 27)
fig3
dev.off()




##### Figure 5: mvOU model output -----


#----------------------------------------------------------#
#  Description: This figure is to present the output of    #
#  the best fitted mvOU model for mammals, birds and       #
#  reptiles, respectively.                                 #
#                                                          #
# Figure 5a: Model selection for Syy matrix.               #
# Figure 5b: The best model structure for each group.      #
# Figure 5c: The final best model for each group.          #
#----------------------------------------------------------#

## Load packages

library(tidyverse)
library(ape)
library(mvSLOUCH)
library(phylopath)
library(PCMBaseCpp) 
library(ggpubr)
library(patchwork)


##### Figure 5a: Syy Matrix, Diagonal vs. UpperTri -----

### AICc: Diagonal
load("Mammal_SyyDiagonal_OUOU_fit_niter500.RData")
mini_mammal_Diagonal <- Mammal_SyyDiagonal_OUOU_fit_niter500 %>%
  group_by(Model) %>%
  summarise(min_AICc = min(AICc)) %>%
  arrange(min_AICc) %>%
  mutate(delta_AICc = min_AICc - min(min_AICc)) %>%  # Calculate ΔAICc
  mutate(Best = ifelse(delta_AICc <= 2, TRUE, FALSE)) %>%
  mutate(Syy = "Diagonal")

### AICc: UpperTri
load("Mammal_SyyUpperTri_OUOU_fit_niter500.RData")
Mammal_SyyUpperTri_OUOU_fit_niter500 <- Mammal_SyyUpperTri_OUOU_fit_niter500 %>%
  filter(Model != "model_null" & Model != "model_full")

mini_mammal_UpperTri <- Mammal_SyyUpperTri_OUOU_fit_niter500 %>%
  group_by(Model) %>%
  summarise(min_AICc = min(AICc)) %>%
  arrange(min_AICc) %>%
  mutate(delta_AICc = min_AICc - min(min_AICc)) %>%  # Calculate ΔAICc
  mutate(Best = ifelse(delta_AICc <= 2, TRUE, FALSE)) %>%
  mutate(Syy = "UpperTri")

### Visualization
Syy_comparision_mammal <- rbind(mini_mammal_Diagonal, mini_mammal_UpperTri)

Fig.syy.m <- ggpaired(Syy_comparision_mammal, x = "Syy", y = "min_AICc", id = "Model",
                      color = "Syy", line.color = "gray", line.size = 0.4, point.size = 1.5,
                      palette = "jco",
                      xlab = NULL,
                      ylab = "AICc",
                      title = "Mammal") +
  xlab(NULL)+
  ylab("AICc")+
  ylim(10200, 11600) +
  stat_compare_means(
    paired = TRUE,
    vjust = 2,
    size = 6
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold")
  )


##### Figure 5b: The best model structure for each group -----
### 1. Mammal    
load("Model6_models.RData")
# Name the model list
Model_6_list <- all_models
# Get a random matrix
mat.m <- Model_6_list[[111]]
# Modify the matrix
mat.m[mat.m=="0"] <- 0
mat.m[mat.m=="+"] <- 0
mat.m[mat.m=="NA"] <- 5
# Create the igraph object
g.m <- graph_from_adjacency_matrix(mat.m, 
                                   mode = "undirected",
                                   weighted = T,
                                   diag = FALSE)

E(g.m)$width <- abs(E(g.m)$weight) * 2 # Scale width for visibility
# Set vertex shape to none, i.e., only display the label
V(g.m)$shape <- "none"
# The vertex label size
V(g.m)$size <- 40  # Vertex size

## Visualization
plot(g.m,
     layout = layout_in_circle(g.m),
     vertex.label.color = "black",
     vertex.label.cex = 3.5,
     edge.arrow.size = 1,
     edge.curved = 0)

### 2. Bird    
load("Model3_models.RData")
# Name the model list
Model_3_list <- all_models
# Get a random matrix
mat.b <- Model_3_list[[11]]
# Modify the matrix
mat.b[mat.b=="0"] <- 0
mat.b[mat.b=="+"] <- 0
mat.b[mat.b=="NA"] <- 5

# Create the igraph object
g.b <- graph_from_adjacency_matrix(mat.b, 
                                   mode = "undirected",
                                   weighted = T,
                                   diag = FALSE)

E(g.b)$width <- abs(E(g.b)$weight) * 2 # Scale width for visibility

# Set vertex shape to none, i.e., only display the label
V(g.b)$shape <- "none"
# The vertex label size
V(g.b)$size <- 40  # Vertex size
## Visualization

plot(g.b,
     layout = layout_in_circle(g.b),
     vertex.label.color = "black",
     vertex.label.cex = 3.5,
     edge.arrow.size = 1,
     edge.curved = 0)


### 3. Reptile    
load("Model10_models.RData")
# Name the model list
Model_10_list <- all_models
# Get a random matrix
mat.r <- Model_10_list[[121]]
# Modify the matrix
mat.r[mat.r=="0"] <- 0
mat.r[mat.r=="+"] <- 0
mat.r[mat.r=="NA"] <- 5

# Create the igraph object
g.r <- graph_from_adjacency_matrix(mat.r, 
                                   mode = "undirected",
                                   weighted = T,
                                   diag = FALSE)

E(g.r)$width <- abs(E(g.r)$weight) * 2 # Scale width for visibility

# Set vertex shape to none, i.e., only display the label
V(g.r)$shape <- "none"
# The vertex label size
V(g.r)$size <- 40  # Vertex size

## Visualization
plot(g.r,
     layout = layout_in_circle(g.r),
     vertex.label.color = "black",
     vertex.label.cex = 3.5,
     edge.arrow.size = 1,
     edge.curved = 0)

##### Figure 5c: The final best model for each group -----

### 1. Mammal    
load("result_constrained_Model_6_mammal_OU_niter100.RData")
# Get the summary of AICc for each run
causal_constrained_mammal_OU_dignosis <- data.frame()
for (n in 1:100) {
  for (i in 1:128) {
    df <- data.frame(
      AICc = result_constrained_Model_6_mammal_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["aic.c"]],
      AIC = result_constrained_Model_6_mammal_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["aic"]],
      BIC = result_constrained_Model_6_mammal_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["bic"]],
      LogLik = result_constrained_Model_6_mammal_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["LogLik"]],
      R2 = result_constrained_Model_6_mammal_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["RSS"]][["R2"]],
      DOF = result_constrained_Model_6_mammal_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["dof"]],
      Iteration = n,
      Model = i,
      Group = "Mammal"
    )
    causal_constrained_mammal_OU_dignosis <- rbind(causal_constrained_mammal_OU_dignosis, df)
  }
}
### The best model is when Iteration = 5, Model = 55
causal_constrained_mammal_OU_best <- result_constrained_Model_6_mammal_OU_niter100[[5]][[55]]

### Get the best A matrix
A.causal_constrained_mammal_OU_best <- causal_constrained_mammal_OU_best[["FinalFound"]][["ParamsInModel"]][["A"]]
# Visualization
A.causal_constrained_mammal_OU_best <- t(A.causal_constrained_mammal_OU_best)
# Create directed graph
g.m <- graph_from_adjacency_matrix(A.causal_constrained_mammal_OU_best, 
                                   mode = "directed",
                                   weighted = TRUE,
                                   diag = FALSE)

# Set edge attributes
E(g.m)$color <- ifelse(E(g.m)$weight < 0, "#2E5A87FF", "#A90C38FF")  # Blue/red
w <- abs(E(g.m)$weight)
E(g.m)$width <- scales::rescale(w, to = c(1, 8))* 2  # min=1, max=8

# Set vertex attributes for pentagrams (5-pointed stars)
V(g.m)$shape <- "none"
V(g.m)$size <- 40  # Vertex size

# Plot with circular layout
plot(g.m,
     layout = layout_in_circle(g.m),
     vertex.label.color = "black",
     vertex.label.cex = 3.5,
     edge.arrow.size = 2.5,
     edge.curved = 0)

### Get the eigen-decomposition results of the A matrix
half.t.mammal_OU <- rbind(causal_constrained_mammal_OU_best[["FinalFound"]][["ParamSummary"]][["phyl.halflife"]][["directions"]],
                          causal_constrained_mammal_OU_best[["FinalFound"]][["ParamSummary"]][["phyl.halflife"]][["halflives"]])

### 2. Bird    
load("result_constrained_Model_3_bird_OU_niter100.RData")
# Get the summary of AICc for each run
causal_constrained_bird_OU_dignosis <- data.frame()
for (n in 1:100) {
  for (i in 1:64) {
    df <- data.frame(
      AICc = result_constrained_Model_3_bird_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["aic.c"]],
      AIC = result_constrained_Model_3_bird_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["aic"]],
      BIC = result_constrained_Model_3_bird_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["bic"]],
      LogLik = result_constrained_Model_3_bird_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["LogLik"]],
      R2 = result_constrained_Model_3_bird_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["RSS"]][["R2"]],
      DOF = result_constrained_Model_3_bird_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["dof"]],
      Iteration = n,
      Model = i,
      Group = "Bird"
    )
    causal_constrained_bird_OU_dignosis <- rbind(causal_constrained_bird_OU_dignosis, df)
  }
}
### The best model is when Iteration = 64, Model = 58
causal_constrained_bird_OU_dignosis %>%
  slice_min(order_by = AICc, n = 1, with_ties = TRUE)

causal_constrained_bird_OU_best <- result_constrained_Model_3_bird_OU_niter100[[64]][[58]]

### Get the best A matrix
A.causal_constrained_bird_OU_best <- causal_constrained_bird_OU_best[["FinalFound"]][["ParamsInModel"]][["A"]]
# Visualization
A.causal_constrained_bird_OU_best <- t(A.causal_constrained_bird_OU_best)
# Create directed graph
g.b <- graph_from_adjacency_matrix(A.causal_constrained_bird_OU_best, 
                                   mode = "directed",
                                   weighted = TRUE,
                                   diag = FALSE)

# Set edge attributes
E(g.b)$color <- ifelse(E(g.b)$weight < 0, "#2E5A87FF", "#A90C38FF")  # Blue/red
w <- abs(E(g.b)$weight)
E(g.b)$width <- scales::rescale(w, to = c(1, 6))* 3  # min=1, max=5

# Set vertex attributes for pentagrams (5-pointed stars)
V(g.b)$shape <- "none"
V(g.b)$size <- 40  # Vertex size

# Plot with circular layout
plot(g.b,
     layout = layout_in_circle(g.b),
     vertex.label.color = "black",
     vertex.label.cex = 3.5,
     edge.arrow.size = 2.5,
     edge.curved = 0)
### Get the eigen-decomposition results of the A matrix
half.t.bird_OU <- rbind(causal_constrained_bird_OU_best[["FinalFound"]][["ParamSummary"]][["phyl.halflife"]][["directions"]],
                        causal_constrained_bird_OU_best[["FinalFound"]][["ParamSummary"]][["phyl.halflife"]][["halflives"]])

### 3. Reptile
load("result_constrained_Model_10_reptile_OU_niter100.RData")
# Get the summary of AICc for each run
causal_constrained_reptile_OU_dignosis <- data.frame()
for (n in 1:100) {
  for (i in 1:128) {
    df <- data.frame(
      AICc = result_constrained_Model_10_reptile_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["aic.c"]],
      AIC = result_constrained_Model_10_reptile_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["aic"]],
      BIC = result_constrained_Model_10_reptile_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["bic"]],
      LogLik = result_constrained_Model_10_reptile_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["LogLik"]],
      R2 = result_constrained_Model_10_reptile_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["RSS"]][["R2"]],
      DOF = result_constrained_Model_10_reptile_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["dof"]],
      Iteration = n,
      Model = i,
      Group = "Reptile"
    )
    causal_constrained_reptile_OU_dignosis <- rbind(causal_constrained_reptile_OU_dignosis, df)
  }
}
### The best model is when Iteration = 52, Model = 83
causal_constrained_reptile_OU_dignosis %>%
  slice_min(order_by = AICc, n = 1, with_ties = TRUE)

causal_constrained_reptile_OU_best <- result_constrained_Model_10_reptile_OU_niter100[[52]][[83]]

### Get the best A matrix
A.causal_constrained_reptile_OU_best <- causal_constrained_reptile_OU_best[["FinalFound"]][["ParamsInModel"]][["A"]]
# Visualization
A.causal_constrained_reptile_OU_best <- t(A.causal_constrained_reptile_OU_best)
# Create directed graph
g.r <- graph_from_adjacency_matrix(A.causal_constrained_reptile_OU_best, 
                                   mode = "directed",
                                   weighted = TRUE,
                                   diag = FALSE)

# Set edge attributes
E(g.r)$color <- ifelse(E(g.r)$weight < 0, "#2E5A87FF", "#A90C38FF")  # Blue/red
w <- abs(E(g.r)$weight)
E(g.r)$width <- scales::rescale(w, to = c(1, 6))* 3  # min=1, max=5

# Set vertex attributes for pentagrams (5-pointed stars)
V(g.r)$shape <- "none"
V(g.r)$size <- 40  # Vertex size

# Plot with circular layout
plot(g.r,
     layout = layout_in_circle(g.r),
     vertex.label.color = "black",
     vertex.label.cex = 3.5,
     edge.arrow.size = 2.5,
     edge.curved = 0)
### Get the eigen-decomposition results of the A matrix
half.t.reptile_OU <- rbind(causal_constrained_reptile_OU_best[["FinalFound"]][["ParamSummary"]][["phyl.halflife"]][["directions"]],
                           causal_constrained_reptile_OU_best[["FinalFound"]][["ParamSummary"]][["phyl.halflife"]][["halflives"]])





