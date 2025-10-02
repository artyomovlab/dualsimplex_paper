library(plotly)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)



colors_for_genes <-  brewer.pal(9, "Blues")
colors_for_samples <-  brewer.pal(9, "Oranges")

#util method to add triangles
add_points_to_plot <- function(plt, points, with_lines=T,  with_points =T,
                               color_points="purple", color_lines="blue", line_type="dashed",
                               points_label = "True solution") {
  points <- as.data.frame(points)
  points$color <- points_label
  
  if (with_lines) {
    plt <-  plt + 
      geom_polygon(
        data = points,
        size = 1,
        fill = NA,
        color = color_lines,
        linetype = line_type,
        aes(fill = points_label)
      ) 
  }
  if (with_points) {
    plt <- plt +
      geom_point(
        data = points,
        color = color_points,
        size = 3,
        aes(fill = points_label)
      )  
  }
  return(plt)
}



get_metric_plot_title <- function(metric) {
  per_metric_titles <- c()
  per_metric_titles["rmse_loss"] <-  expression("RMSE")
  per_metric_titles["pearson_loss"] <-  expression(1 - pearson^2 )
  per_metric_titles["spearman_loss"] <-   expression(1 - spearman^2 )
  return(per_metric_titles[[metric]])
}

plot_single_boxplot <-  function(metric_results_for_datasets_and_methods, 
                    comparisons, 
                    title,
                    reference_group ="DualSimplex",
                    metric_column="metric_value", metric_title = "metric_value", 
                    custom_colors = brewer.pal(length(unique(metric_results_for_datasets_and_methods$method)),"Set1"),
                    log_scale = F, test="wilcox.test") {
  
  result_plot <- ggplot(total_result, aes(x=method, y=.data[[metric_column]], fill=method)) +
    geom_boxplot(position = position_dodge2(width = 0.1)) +
    stat_compare_means(comparisons=comparisons, label = "p.signif", tip.length = 0,  ref.group = "DualSimplex", size=6, method=test) +
    geom_point(position =position_jitterdodge(jitter.width = 0.1)) +
    scale_fill_manual(values=custom_colors)+
    theme_classic(base_size=25,base_family = 'sans') +
    theme(plot.title = element_text(face='bold', size=21)) + 
    theme(axis.title.x  = element_text(size=20)) +  
    theme(axis.title.y  = element_text(size=20))+
    theme(axis.text.x  = element_text(size=18)) +  
    theme(axis.text.y  = element_text(size=14)) + 
    theme(legend.position = "none", legend.text=element_text(size=16))+
    theme(legend.title=element_text(size=16)) +
    theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1)) + 
    labs(title = title,x=NULL, y=metric_title
    ) + 
    theme(plot.margin = margin(t = 3, r = 0, b = 3, l = 5, unit = "pt")) 
    if (log_scale) {
      result_plot <-  result_plot + scale_y_continuous(trans='log10', limits = c(NA, NA))
    }
    else {
    result_plot <-  result_plot + ylim(0, NA)
    } 
  return(result_plot)
  
}








plot_3d <- function(projection_points, axes="V", size =5, colors=brewer.pal(9,"Greys"), shift=0) {
  t <- list(
    family = "Helvetica",
    size = 30)
  dims <-  dim(projection_points)[[2]]
  to_plot_points <- as.data.frame(projection_points[,1:dims])
  names_for_col <- paste0(axes, c(1:dims))
  colnames(to_plot_points) <- names_for_col
  fig_Omega <- plot_ly(to_plot_points,  
                      width = 500, height = 500,   
                       mode   = 'markers',  type   = 'scatter3d') 
  
  fig_Omega <-fig_Omega %>% add_trace( x = to_plot_points[[names_for_col[[shift+1]]]],
                                       y = to_plot_points[[names_for_col[[shift+2]]]],
                                       z = to_plot_points[[names_for_col[[shift+3]]]]
                                       #, 
                                      # size=size 
                                       ,marker = list(symbol = 'circle', 
                                                      color = colors[4],
                                                      size=size,
                                                      opacity = 0.9, 
                                                      line = list(
                                                        color = colors[7],
                                                        opacity=0.9,
                                                        width = 2
                                       )) 
                                       )
  
  fig_Omega <- fig_Omega %>% layout(showlegend = F,
                                    
                                    scene = list(
                                      aspectmode='cube',
                                      xaxis = list(
                                        title = names_for_col[[shift+1]],
                                      titlefont = t
                                      
                                      ),
                                     yaxis = list(
                                        title =names_for_col[[shift+2]],
                                       titlefont = t
                                      
                                       ),
                                     zaxis = list(
                                        title = names_for_col[[shift+3]],
                                       titlefont = t
                                     
                                       ))
  ) %>% config(toImageButtonOptions = list(format = "svg", width = 600,
                                           height = 600))
  
  return(fig_Omega)
}
