# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, FactoMineR, factoextra, labelled, corrplot, plotly, ggpmisc, vroom, rstatix, modeest, multimode))

root <- '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO'

isos <- c('KEN','MLI','NGA','SDN','SEN','UGA','ZWE')

iso <- 'KEN'

# Climate clusters using descriptive statistics

# Load descriptive statistics per country and put them all in a single table
stts <- isos %>%
  purrr::map(.f = function(iso){
    stts <- read.csv(paste0(root,'/data/',iso,'/_results/cluster_results/climate/climate_reg_cluster_statistics.csv'))
    stts <- stts %>% dplyr::select(variable, clust, mean, median, sd, mad, min, max)
    stts$iso <- iso
    return(stts)
  }) %>%
  dplyr::bind_rows()

vrbl <- list()
for(i in 1:length(isos)){
  vrbl[[i]] <- stts[stts$iso == isos[i],'variable']
}

unq_vars <- Reduce(intersect, vrbl); rm(vrbl, i)

stts <- stts[stts$variable %in% unq_vars,]

stts <- stts %>%
  dplyr::select(variable, clust, iso, median) %>%
  tidyr::pivot_wider(names_from = 'variable', values_from = median)
stts$key <- paste0(stts$iso,'-',stts$clust)
stts$iso <- stts$clust <- NULL
stts <- base::as.data.frame(stts)
rownames(stts) <- stts$key; stts$key <- NULL

# stts$key <- paste0(stts$variable,'-',stts$clust,'-',stts$iso)
# stts$iso <- stts$variable <- stts$clust <- NULL
# stts <- base::as.data.frame(stts)
# rownames(stts) <- stts$key; stts$key <- NULL

# Perform a PCA over this table
pca <- stts %>%
  FactoMineR::PCA(X = ., scale.unit = T, ncp = ncol(.), graph = T)

corrplot::corrplot(pca$var$contrib, is.corr = F) # Variable contributions
corrplot::corrplot(pca$var$cos2, is.corr = F) # Quality of representation

fviz_pca_var(pca,
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

fviz_pca_biplot(pca,
                # Individuals
                pointsize = 2,
                palette = "jco",
                addEllipses = F,
                # Variables
                col.var = "contrib",
                gradient.cols = "RdYlBu",
                legend.title = list(color = "Contrib")
)

# Perform a hierarchical cluster over the PC axes
hcpc <- FactoMineR::HCPC(res = pca, nb.clust = -1, consol = T, graph = T)
fviz_dend(hcpc, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)
fviz_cluster(hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)
hcpc$desc.var$quanti

stts %>%
  dplyr::mutate(ind = rownames(.)) %>%
  dplyr::filter(ind %in% paste0('ZWE-',1:4)) %>%
  ggplot2::ggplot(aes(x = ind, y = medn_prec)) +
  ggplot2::geom_point() +
  ggplot2::coord_flip()

# Climate clusters using relative changes

# Load relative changes per country and put them all in a single table
rchg <- isos %>%
  purrr::map(.f = function(iso){
    rchg <- read.csv(paste0(root,'/data/',iso,'/_results/cluster_results/climate/climate_cluster_rel_change.csv'))
    rchg <- rchg[,c(1,grep(pattern = 'rel_change', x = names(rchg)))]
    # rchg <- rchg[rchg$Variables %in% unq_vars,]
    rchg$iso <- iso
    rchg <- rchg %>%
      tidyr::pivot_longer(cols      = grep(pattern = 'rel_change', x = names(rchg)),
                          names_to  = 'clust',
                          values_to = 'change') %>%
      tidyr::pivot_wider(names_from = 'Variables', values_from = 'change')
    return(rchg)
  }) %>%
  dplyr::bind_rows()

sds <- stts %>%
  dplyr::filter(iso == 'KEN') %>%
  dplyr::select(variable, clust, median) %>%
  tidyr::pivot_wider(names_from = 'variable', values_from = 'median') %>%
  dplyr::mutate(clust = factor(clust)) %>%
  apply(X = ., MARGIN = 2, FUN = sd, na.rm = T)
avg <- stts %>%
  dplyr::filter(iso == 'KEN') %>%
  dplyr::select(variable, clust, median) %>%
  tidyr::pivot_wider(names_from = 'variable', values_from = 'median') %>%
  colMeans(x = ., na.rm = T)

sort(sds/avg * 100, decreasing = T)

rchg %>%
  dplyr::filter(iso == 'KEN') %>%
  .[,3:ncol(.)] %>%
  base::colMeans(x = ., na.rm = T) %>%
  sort(decreasing = T) %>%
  barplot()

rchg$key <- paste0(rchg$iso,'-',rchg$clust)
rchg$iso <- rchg$clust <- NULL
rchg <- base::as.data.frame(rchg)
rownames(rchg) <- rchg$key; rchg$key <- NULL

# Perform a PCA over this table
pca <- rchg %>%
  FactoMineR::PCA(X = ., scale.unit = T, ncp = ncol(.), graph = T)

corrplot::corrplot(pca$var$contrib, is.corr = F) # Variable contributions
corrplot::corrplot(pca$var$cos2, is.corr = F) # Quality of representation

fviz_pca_var(pca,
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

fviz_pca_biplot(pca,
                # Individuals
                pointsize = 2,
                palette = "jco",
                addEllipses = F,
                # Variables
                col.var = "contrib",
                gradient.cols = "RdYlBu",
                legend.title = list(color = "Contrib")
)

# Perform a hierarchical cluster over the PC axes
hcpc <- FactoMineR::HCPC(res = pca, nb.clust = -1, consol = T, graph = T)
fviz_dend(hcpc, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)
fviz_cluster(hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)
hcpc$desc.var$quanti

# Raster extracted values

iso <- 'KEN'
rawv <- vroom::vroom(paste0(root,'/data/',iso,'/_results/cluster_results/climate/climate_reg_cluster_values_extracted.csv'), delim = ',')
prob <- table(rawv$clust)/nrow(rawv)
n <- 500
ns <- round(prob * n, 0)
set.seed(1235)
strat_sample <- DescTools::Strata(x = rawv, stratanames = 'clust', size = ns, method = 'srswor')
strat_sample$stratum <- strat_sample$size <- strat_sample$id <- NULL

names(strat_sample)

strat_sample$clust <- factor(strat_sample$clust)

statip::mfv(x = strat_sample$trnd_cwdf[strat_sample$clust == 2])
statip::mfv1(x = strat_sample$medn_prec[strat_sample$clust == 2])
modeest::mlv(x = strat_sample$trnd_cwdf[strat_sample$clust == 2], method = 'mfv')
modeest::mlv(x = strat_sample$medn_prec[strat_sample$clust == 2], method = 'meanshift')
median(x = strat_sample$medn_prec[strat_sample$clust == 2])
mean(x = strat_sample$medn_prec[strat_sample$clust == 2])
hist(strat_sample$medn_prec[strat_sample$clust == 2], probability = T)
lines(density(strat_sample$medn_prec[strat_sample$clust == 2]), col = 4)
abline(v = modeest::mlv(x = strat_sample$medn_prec[strat_sample$clust == 2], method = 'mfv'), col = 'red')

multimode::locmodes(data = strat_sample$medn_prec[strat_sample$clust == 2], mod0 = 1)

pca <- strat_sample %>%
  dplyr::select(-clust) %>%
  FactoMineR::PCA(X = ., scale.unit = T, ncp = ncol(.), graph = F)

# Eigen values
factoextra::fviz_eig(pca, addlabels = T)

# Correlation circle
factoextra::fviz_pca_var(pca,
                         col.var = "cos2",
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                         repel = T # Avoid text overlapping
)

# Individuals map
factoextra::fviz_pca_ind(pca,
                         col.ind = "cos2",
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = TRUE # Avoid text overlapping (slow if many points)
)

factoextra::fviz_pca_ind(pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = strat_sample$clust, # color by groups
                         palette = MetBrewer::met.brewer(name = 'Egypt', n = length(levels(strat_sample$clust))),
                         addEllipses = T, # Concentration ellipses
                         legend.title = "Groups"
)

car::scatter3d(x = pca$ind$coord[,1],
               y = pca$ind$coord[,2],
               z = pca$ind$coord[,3],
               groups = strat_sample$clust,
               surface = F,
               grid = F,
               ellipsoid = T,
               surface.col = MetBrewer::met.brewer(name = 'Egypt', n = length(levels(strat_sample$clust))),
               xlab = "PC1",
               ylab = "PC2",
               zlab = "PC3")

testRes <- corrplot::cor.mtest(mat = strat_sample %>%
                                 dplyr::select(-clust),
                               conf.level = 0.95,
                               method = 'pearson',
                               use = 'pairwise.complete.obs')
M <- strat_sample %>%
  dplyr::select(-clust) %>%
  cor(method = 'pearson', use = 'pairwise.complete.obs')

corrplot::corrplot(M,
                   p.mat = testRes$p,
                   lowCI = testRes$lowCI,
                   uppCI = testRes$uppCI,
                   addrect = 3,
                   rect.col = 'navy',
                   plotC = 'rect',
                   cl.pos = 'n')

strat_sample %>%
  dplyr::select(clust,medn_prec:NWLD_trnd) %>%
  tidyr::pivot_longer(cols = medn_prec:NWLD_trnd, names_to = 'Indicator', values_to = 'Value') %>%
  ggplot2::ggplot(aes(x = clust, y = Value)) +
  ggplot2::geom_boxplot() +
  ggplot2::facet_wrap(~Indicator, scales = 'free')

strat_sample %>%
  dplyr::filter(clust == 4) %>%
  dplyr::select(medn_prec:NWLD_trnd) %>%
  cor(method = 'kendall', use = 'pairwise.complete.obs') %>%
  corrplot::corrplot()

ade4::mantel.rtest(m1 = 1 - strat_sample %>%
                     dplyr::filter(clust == 1) %>%
                     dplyr::select(medn_prec:NWLD_trnd) %>%
                     cor(method = 'kendall', use = 'pairwise.complete.obs') %>%
                     as.dist(),
                   m2 = 1 - strat_sample %>%
                     dplyr::filter(clust == 4) %>%
                     dplyr::select(medn_prec:NWLD_trnd) %>%
                     cor(method = 'kendall', use = 'pairwise.complete.obs') %>%
                     as.dist(),
                   nrepet = 1000)

rawv_medn_prec %>%
  ggplot2::ggplot(aes(x = factor(clust), y = trnd_cwdf)) +
  ggplot2::geom_boxplot()
strat_sample %>%
  ggplot2::ggplot(aes(x = factor(clust), y = trnd_cwdf)) +
  ggplot2::geom_boxplot()

agricolae::kruskal(y     = rawv_medn_prec$trnd_cwdf,
                   trt   = factor(rawv_medn_prec$clust),
                   alpha = 0.05,
                   p.adj = 'bonferroni', console = T)
agricolae::kruskal(y     = strat_sample$trnd_cwdf,
                   trt   = factor(strat_sample$clust),
                   alpha = 0.05,
                   p.adj = 'bonferroni', console = T)

rstatix::kruskal_effsize()

strat_sample$clust <- factor(strat_sample$clust)

gg <- strat_sample %>%
  ggplot2::ggplot(aes(x = trnd_cwdf, color = clust, fill = clust)) +
  ggplot2::geom_density(alpha = 0.6) +
  viridis::scale_fill_viridis(discrete = TRUE) +
  viridis::scale_color_viridis(discrete = TRUE) +
  hrbrthemes::theme_ipsum() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::ylab("Density") +
  ggplot2::xlab("Slope")

strat_sample %>%
  ggplot(aes(x = medn_aet)) +
  geom_histogram(aes(y = ..density.., color = clust), 
                   fill = "white",
                   position = "identity")+
  geom_density(aes(color = clust), size = 1)
