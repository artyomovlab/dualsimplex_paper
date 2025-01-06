plot_3d <- function(projection_points, axes="V", size =5, colors=brewer.pal(9,"Greys"), shift=0) {
  #projection_points <-  round(projection_points,15)
  
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
