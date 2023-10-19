
test_that("vr_graphs can detect incorrect parameters properly",{
  
  expect_error(vr_graphs(X = 2,eps = c(-1,0)),"eps")
  expect_error(vr_graphs(X = 2,eps = c(1),distance_mat = NULL),"NULL")
  expect_error(vr_graphs(X = 2,eps = c(1),distance_mat = c(T,F)),"single")
  expect_error(vr_graphs(X = 2,eps = c(1),distance_mat = NA),"NA")
  
  expect_error(vr_graphs(X = 2,eps = c(1)),"matrix")
  expect_error(vr_graphs(X = data.frame(),eps = c(1)),"one")
  expect_error(vr_graphs(X = data.frame(x = c(1,2,3),y = c(1,NA,3)),eps = c(1)),"missing")
  expect_error(vr_graphs(X = data.frame(x = 1:3,y = 1:3),distance_mat = T,eps = c(1)),"distance")
  expect_error(vr_graphs(X = data.frame(x = 1:3,y = as.character(1:3)),eps = c(1)),"numeric")
  
  expect_error(vr_graphs(X = data.frame(x = 1:10,y = 1:10),eps = c(1),return_clusters = NULL),"NULL")
  expect_error(vr_graphs(X = data.frame(x = 1:10,y = 1:10),eps = c(1),return_clusters = c(T,F)),"single")
  expect_error(vr_graphs(X = data.frame(x = 1:10,y = 1:10),eps = c(1),return_clusters = NA),"NA")
  
})

test_that("vr_graphs is working properly",{
  
  # simulate data from the unit circle and calculate
  # its diagram
  df <- TDA::circleUnif(n = 25)
  diag <- TDA::ripsDiag(df,maxdimension = 1,maxscale = 2)

  # get minimum death radius of any data cluster
  min_death_H0 <- min(diag$diagram[which(diag$diagram[,1] == 0),3L])

  # get birth and death radius of the loop
  loop_birth <- as.numeric(diag$diagram[nrow(diag$diagram),2L])
  loop_death <- as.numeric(diag$diagram[nrow(diag$diagram),3L])

  # compute Rips-Vietoris complexes at radii half of
  # min_death_H0 and the mean of loop_birth and
  # loop_death, returning clusters
  comp <- vr_graphs(X = df,eps = c(0.5*min_death_H0,(loop_birth + loop_death)/2))

  # verify that there are 25 clusters for the smaller radius
  # and 1 cluster for the larger radius
  expect_identical(length(comp$graphs[[1]]$clusters),25L)
  expect_identical(length(comp$graphs[[2]]$clusters),1L)
  
  # get value between min_death_H0 and next value
  d <- as.matrix(dist(df))
  inds <- which(d > min_death_H0,arr.ind = T)
  vals <- d[inds]
  vals <- vals[order(vals)]
  vals <- unique(vals)
  next_val <- vals[[1]]
  next_next_val <- vals[[2]]
  
  # check that the rips complex at a next radius has at most 24 clusters
  comp <- vr_graphs(X = df,eps = c(0.5*(next_val + next_next_val)))
  expect_lte(length(comp$complexes[[1]]$clusters),24L)
  
  # make sure that when data has rownames that they get saved properly
  rownames(df) <- paste0("V",1:25)
  comp <- vr_graphs(X = df,eps = c(0.5*min_death_H0,(loop_birth + loop_death)/2))
  
})

test_that("plot_vr_graph can detect incorrect parameters properly",{
  
  skip_if_not_installed("igraph")
  
  expect_error(plot_vr_graph(graphs = c()),"vr_graphs")
  expect_error(plot_vr_graph(graphs = list(1,2,3)),"vr_graphs")
  expect_error(plot_vr_graph(graphs = list(data = 1,graphs = 2)),"vr_graphs")
  expect_error(plot_vr_graph(graphs = list(data = c(1,2,3),graphs = c())),"vr_graphs")
  
  # simulate data from the unit circle and calculate
  # its diagram
  df <- TDA::circleUnif(n = 25)
  diag <- TDA::ripsDiag(df,maxdimension = 1,maxscale = 2)
  
  # get minimum death radius of any data cluster
  min_death_H0 <- min(diag$diagram[which(diag$diagram[,1] == 0),3L])
  
  # get birth and death radius of the loop
  loop_birth <- as.numeric(diag$diagram[nrow(diag$diagram),2L])
  loop_death <- as.numeric(diag$diagram[nrow(diag$diagram),3L])
  
  # compute Rips-Vietoris complexes at radii half of
  # min_death_H0 and the mean of loop_birth and
  # loop_death, returning clusters
  comp <- vr_graphs(X = df,eps = c(0.5*min_death_H0,(loop_birth + loop_death)/2))
  
  expect_error(plot_vr_graph(graphs = comp,eps = 0),"positive")
  expect_error(plot_vr_graph(graphs = comp,eps = 3),"scale")
  
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,plot_isolated_vertices = NULL),"NULL")
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,plot_isolated_vertices = c(T,F)),"single")
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,plot_isolated_vertices = 2),"boolean")
  
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,return_layout = NULL),"NULL")
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,return_layout = c(T,F)),"single")
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,return_layout = 2),"boolean")
  
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,vertex_labels = NULL),"NULL")
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,vertex_labels = c(T,F)),"single")
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,vertex_labels = "2"),"boolean")
  
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,cols = NA),"NA")
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,cols = c(1,2)),"character")
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,cols = c("1",NA)),"missing")
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,cols = c("1","2")),"vertices")
  
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,component_of = c(1,2)),"single")
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,component_of = NaN),"NaN")
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,component_of = 26),"vertices")
  
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,title = NA),"character")
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,title = c("1","2")),"single")
  
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,layout = NA),"matrix")
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,layout = matrix(data = 1)),"columns")
  expect_error(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0,layout = matrix(data = c(1,2),ncol = 2)),"rows")
  
  expect_warning(plot_vr_graph(graphs = comp,eps = 0.5*min_death_H0, component_of = 1),"empty")

})

test_that("plot_vr_graph is working properly",{
  
  # simulate data from the unit circle and calculate
  # its diagram
  df <- TDA::circleUnif(n = 25)
  diag <- TDA::ripsDiag(df,maxdimension = 1,maxscale = 2)
  
  # get minimum death radius of any data cluster
  min_death_H0 <- min(diag$diagram[which(diag$diagram[,1] == 0),3L])
  
  # get birth and death radius of the loop
  loop_birth <- as.numeric(diag$diagram[nrow(diag$diagram),2L])
  loop_death <- as.numeric(diag$diagram[nrow(diag$diagram),3L])
  
  # compute Rips-Vietoris complexes at radii half of
  # min_death_H0 and the mean of loop_birth and
  # loop_death, returning clusters
  comp <- vr_graphs(X = df,eps = c(0.5*min_death_H0,(loop_birth + loop_death)/2))
  
  # check layout functionality
  suppressWarnings({
    
    layout <- plot_vr_graph(comp,eps = 0.5*min_death_H0,return_layout = T)
    
  })
  expect_equal(dim(layout),c(0,2))
  layout <- plot_vr_graph(comp,eps = 0.5*min_death_H0,return_layout = T,plot_isolated_vertices = T)
  expect_equal(nrow(layout),25)
  layout <- plot_vr_graph(comp,eps = 0.5*min_death_H0,return_layout = T,plot_isolated_vertices = T,component_of = 1)
  expect_equal(nrow(layout),1)
  layout <- plot_vr_graph(comp,eps = (loop_birth + loop_death)/2,return_layout = T,component_of = 1)
  expect_equal(nrow(layout),25)
  rownames(df) <- paste("V",1:25,sep = "")
  comp <- vr_graphs(X = df,eps = c(0.5*min_death_H0,(loop_birth + loop_death)/2))
  layout <- plot_vr_graph(comp,eps = 0.5*min_death_H0,return_layout = T,plot_isolated_vertices = T)
  expect_equal(nrow(layout),25)
  layout <- plot_vr_graph(comp,eps = 0.5*min_death_H0,return_layout = T,plot_isolated_vertices = T,component_of = "V1")
  expect_equal(nrow(layout),1)
  layout <- plot_vr_graph(comp,eps = 0.5*min_death_H0,return_layout = T,plot_isolated_vertices = T,layout = layout)
  expect_equal(nrow(layout),1)
  
})

