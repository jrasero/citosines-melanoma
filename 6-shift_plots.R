require(openxlsx)
require(ggplot2)
require(cowplot)
require(rogme)

my_plot_sf<-function(sf){
  xplot<- seq(0.1, 0.9, 0.1)
  yplot<- sf[[1]]$difference
  
  ylim <- max(max(abs(sf[[1]]$ci_upper)), max(abs(sf[[1]]$ci_lower)))
  ylim <- c(-ylim, ylim)
  midpt <- (nrow(sf[[1]]) - 1)/2 + 1
  xintercept<-0.5
  yintercept<-0
  sf.df <- data.frame("x" = xplot, "difference"= yplot)
  sf.df$sign <- sign(sf.df$difference)
  sf.df$ci_lower<- sf[[1]]$ci_lower
  sf.df$ci_upper<- sf[[1]]$ci_upper
  sf.df$deco <- c(seq(1, midpt), seq(midpt - 1, 1))
  q_line_alpha <- 0.5
  q_line_size <- 1.5
  symb_fill = c("darkred", "darkblue")
  symb_col <- c("black")
  q_line_col <- "grey50"
  theme2_alpha <- c(0.4, 1)
  symb_size<-5
  
  psf <- ggplot(sf.df, aes_string(x = xplot, y = "difference")) + 
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) + 
    geom_vline(xintercept = xintercept, linetype = 2, 
               alpha = 0.5) + theme_bw() + 
    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), 
          axis.title.x = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 16, face = "bold")) + 
    scale_y_continuous(limits = ylim)
  
  psf <- psf + geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), 
                              colour = "white", size = 0.5) + 
    geom_linerange(aes(ymin = ci_lower, ymax = ci_upper, 
                       colour = factor(sign), alpha = factor(deco)), 
                   size = 0.5) + 
    scale_color_manual(values = symb_fill, guide = FALSE) + 
    scale_alpha_discrete(range = theme2_alpha, guide = FALSE) + 
    geom_line(colour = q_line_col, alpha = q_line_alpha, 
              linetype = "solid", size = q_line_size) + 
    geom_point(colour = "black", fill = "white", 
               size = symb_size, shape = 21) + 
    geom_point(aes(fill = factor(sign), 
                   alpha = factor(deco)), 
               colour = symb_col, 
               size = symb_size, 
               shape = 21) + 
    scale_fill_manual(values = symb_fill, guide = FALSE) + 
    scale_alpha_discrete(range = theme2_alpha, guide = FALSE)
  
  psf
}

sub<-read.csv("data/ml_data.csv")

vars.interest<-names(sub)[c(1, 3:10)]
names.vars<- list("BRESLOW", "GM-CSF", "IL-4", "IL-6", 
               "IL-10", "IL-17A", "IFN-\u03b3",
               "TGF-\u03b2","DCD")
names(names.vars)<-vars.interest

unit.vars<-list("mm", "pg/mL", 
                "pg/mL", "pg/mL", 
                "pg/mL", "pg/mL",
                "pg/mL",
                "ng/mL", "ng/mL")

names(unit.vars)<-vars.interest

for(newVar in vars.interest){
  g1<-na.omit(sub[sub$EVOL_METASTASIS==0, newVar])
  g2<-na.omit(sub[sub$EVOL_METASTASIS==1, newVar])
  df <- mkt2(g1, g2)
  
  # Compute shift function
  sf <- shifthd_pbci(data = df, formula = obs ~ gr)
  
  #PLOT DISTRIBUTION AND DECILES
  p <- plot_scat2(df,
                  xlabel = "",
                  ylabel = unit.vars[newVar],
                  alpha = .3,
                  shape = 21,
                  colour = "grey10",
                  fill = "grey90")# scatterplots
  
  
  p <- plot_hd_bars(p, col = "grey21", width = 0.5, q_size = 0.5, 
                     md_size = 1.5)
  p<-plot_hd_links(p, sf[[1]], link_col = c("darkred", "darkblue"))
  
  p <- p + 
    coord_flip() + # flip axes 
    scale_x_discrete(labels=c("Group2" = "MET", "Group1" = "FREE")) + 
    theme(axis.text.x = element_text(size=15, face="bold"),
          axis.text.y = element_text(size=15, face="bold"),
          axis.title.y = element_text(size=25, face="bold")) +
    ggtitle(names.vars[newVar]) + 
    theme(plot.title = element_text(size = 25, face="bold", hjust = 0.5))
    
    save_plot(paste0("plots/deciles_", paste0(names.vars[newVar], ".svg")),
              p, base_height = 7, base_width = 7)

  ## PLOT SHIFT FUNCTION
  psf<-my_plot_sf(sf)
  
  psf <- psf + 
    xlab("Group Deciles") + 
    ylab(paste0("FREE - MET \n decile difference \n", 
               paste0(paste0("(", unit.vars[newVar]), ")"))) +
    scale_x_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), 
                       labels = c("1st", "2nd", "3rd", "4th", "5th \n (median)", 
                                  "6th", "7th", "8th", "9th"))
  save_plot(paste0("plots/shift_", paste0(names.vars[newVar], ".svg")),
            psf, base_height = 7, base_width = 7)
  
  # Arrange both plots
  both.p<-plot_grid(p, psf, ncol = 1, nrow = 2)
  save_plot(paste0("plots/decile_shift_", paste0(names.vars[newVar], ".svg")),
            both.p, base_height = 7, base_width = 7)
  
}

