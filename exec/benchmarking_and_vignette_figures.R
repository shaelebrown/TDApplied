
# benchmarking script for generating vignette figures
# uses library bench

#### diagram_distance vs. TDA's wasserstein ####

# generate persistence diagrams from Tori and spheres with  100,200,...,1000 data points.
runtimes_distance <- data.frame(n_row = numeric(),package = character(),time_in_sec = numeric())
for(n_row in seq(100,1000,100)){
  for(iteration in 1:10)
  {
    # simulate pair of diagrams from the desired shapes
    diagram_torus = ripsDiag(X = TDA::torusUnif(n = n_row,a = 1,c = 2),
                             maxdimension = 2,maxscale = 2)
    diagram_sphere = ripsDiag(X = TDA::sphereUnif(n = n_row,d = 2,r = 1),
                              maxdimension = 2,maxscale = 1)
    # compute their wasserstein distances in all dimensions and benchmark
    start_time_TDApplied = Sys.time()
    diagram_distance(D1 = diagram_torus,D2 = diagram_sphere,dim = 0,
                     p = 2,distance = "wasserstein")
    diagram_distance(D1 = diagram_torus,D2 = diagram_sphere,dim = 1,
                     p = 2,distance = "wasserstein")
    diagram_distance(D1 = diagram_torus,D2 = diagram_sphere,dim = 2,
                     p = 2,distance = "wasserstein")
    end_time_TDApplied = Sys.time()
    time_diff_TDApplied = as.numeric(end_time_TDApplied - start_time_TDApplied,units = "secs")
    start_time_TDA = Sys.time()
    TDA::wasserstein(Diag1 = diagram_torus$diagram,Diag2 = diagram_sphere$diagram,
                     dimension = 0,p = 2)
    TDA::wasserstein(Diag1 = diagram_torus$diagram,Diag2 = diagram_sphere$diagram,
                     dimension = 1,p = 2)
    TDA::wasserstein(Diag1 = diagram_torus$diagram,Diag2 = diagram_sphere$diagram,
                     dimension = 2,p = 2)
    end_time_TDA = Sys.time()
    time_diff_TDA = as.numeric(end_time_TDA - start_time_TDA,units = "secs")
    runtimes_distance = rbind(runtimes_distance,data.frame(n_row = n_row,
                                                     package = "TDApplied",
                                                     time_in_sec = time_diff_TDApplied))
    runtimes_distance = rbind(runtimes_distance,data.frame(n_row = n_row,
                                                     package = "TDA",
                                                     time_in_sec = time_diff_TDA))
  }
  print(paste0("Done ",n_row," rows"))
}
# compute means and sd's at each value of rows for both packages
summary_table_distance = data.frame(n_row = numeric(),mean = numeric(),sd = numeric(),
                           package = character())
for(n_row in seq(100,1000,100))
{
  for(p in c("TDApplied","TDA"))
  {
    result = data.frame(n_row = n_row,
                        mean = mean(runtimes_distance[which(runtimes_distance$n_row == n_row
                                                         & runtimes_distance$package == p),
                                                   3]),
                        sd = sd(runtimes_distance[which(runtimes_distance$n_row == n_row 
                                                     & runtimes_distance$package == p),
                                               3]),
                        package = p)
    summary_table_distance = rbind(summary_table_distance,result)
  }
}
# plot table
plot(summary_table_distance$n_row[summary_table_distance$package=="TDA"], 
     summary_table_distance$mean[summary_table_distance$package=="TDA"], 
     type="b",
     xlim=range(summary_table_distance$n_row),
     ylim=range(0,summary_table_distance$mean+1.96*summary_table_distance$sd/sqrt(10)),
     xlab = "Points in shape",ylab = "Mean execution time (sec)")
lines(summary_table_distance$n_row[summary_table_distance$package=="TDApplied"],
      summary_table_distance$mean[summary_table_distance$package=="TDApplied"], 
      col=2, type="b")
legend(x = 200,y = 2000,legend = c("TDApplied","TDA"),
       col = c("red","black"),lty = c(1,1),cex = 0.8)
arrows(summary_table_distance$n_row[summary_table_distance$package == "TDApplied"], 
       summary_table_distance$mean[summary_table_distance$package == "TDApplied"]
       -1.96*summary_table_distance$sd[summary_table_distance$package == "TDApplied"]/sqrt(10),
       summary_table_distance$n_row[summary_table_distance$package == "TDApplied"], 
       summary_table_distance$mean[summary_table_distance$package == "TDApplied"]
       +1.96*summary_table_distance$sd[summary_table_distance$package == "TDApplied"]/sqrt(10), 
       length=0.05, angle=90, code=3,col = "red")
arrows(summary_table_distance$n_row[summary_table_distance$package == "TDA"], 
       summary_table_distance$mean[summary_table_distance$package == "TDA"]
       -1.96*summary_table_distance$sd[summary_table_distance$package == "TDA"]/sqrt(10), 
       summary_table_distance$n_row[summary_table_distance$package == "TDA"], 
       summary_table_distance$mean[summary_table_distance$package == "TDA"]
       +1.96*summary_table_distance$sd[summary_table_distance$package == "TDA"]/sqrt(10), 
       length=0.05, angle=90, code=3,col = "black")

#### diagram_distance vs persim's wasserstein ####

# import two python modules
# reticulate::py_install("persim") # if running for the first time
persim <- reticulate::import("persim")
ripser <- import_ripser()

# generate persistence diagrams from Tori and spheres with  100,200,...,1000 data points.
runtimes_language <- data.frame(n_row = numeric(),package = character(),
                                time_in_sec = numeric())
for(n_row in seq(100,1000,100)){
  
  for(iteration in 1:10)
  {
    # simulate pair of diagrams from the desired shapes
    torus = TDA::torusUnif(n = n_row,a = 1,c = 2)
    sphere = TDA::sphereUnif(n = n_row,d = 2,r = 1)
    diagram_torus = ripsDiag(X = torus,
                             maxdimension = 2,maxscale = 2)
    diagram_sphere = ripsDiag(X = sphere,
                              maxdimension = 2,maxscale = 1)
    diagram_torus_py = ripser$ripser(torus,maxdim = 2,thresh = 2)$dgms
    diagram_torus_py[[1]][which(diagram_torus_py[[1]][,2] == Inf),2] = 2
    diagram_torus_py[[2]][which(diagram_torus_py[[2]][,2] == Inf),2] = 2
    diagram_torus_py[[3]][which(diagram_torus_py[[3]][,2] == Inf),2] = 2
    diagram_sphere_py = ripser$ripser(sphere,maxdim = 2,thresh = 1)$dgms
    diagram_sphere_py[[1]][which(diagram_sphere_py[[1]][,2] == Inf),2] = 2
    diagram_sphere_py[[2]][which(diagram_sphere_py[[2]][,2] == Inf),2] = 2
    diagram_sphere_py[[3]][which(diagram_sphere_py[[3]][,2] == Inf),2] = 2
    
    # compute their wasserstein distances in all dimensions and benchmark
    start_time_TDApplied = Sys.time()
    diagram_distance(D1 = diagram_torus,D2 = diagram_sphere,dim = 0,
                     p = 2,distance = "wasserstein")
    diagram_distance(D1 = diagram_torus,D2 = diagram_sphere,dim = 1,
                     p = 2,distance = "wasserstein")
    diagram_distance(D1 = diagram_torus,D2 = diagram_sphere,dim = 2,
                     p = 2,distance = "wasserstein")
    end_time_TDApplied = Sys.time()
    time_diff_TDApplied = as.numeric(end_time_TDApplied - start_time_TDApplied,units = "secs")
    
    start_time_persim = Sys.time()
    persim$wasserstein(diagram_torus_py[[1]],diagram_sphere_py[[1]])
    persim$wasserstein(diagram_torus_py[[2]],diagram_sphere_py[[2]])
    persim$wasserstein(diagram_torus_py[[3]],diagram_sphere_py[[3]])
    end_time_persim = Sys.time()
    time_diff_persim = as.numeric(end_time_persim - start_time_persim,units = "secs")
    
    runtimes_language = rbind(runtimes_language,data.frame(n_row = n_row,
                                                           package = "TDApplied",
                                                           time_in_sec = time_diff_TDApplied))
    runtimes_language = rbind(runtimes_language,data.frame(n_row = n_row,
                                                           package = "persim",
                                                           time_in_sec = time_diff_persim))
    
  }
  print(paste0("Done ",n_row," rows"))
  
}
# compute means and sd's at each value of rows for both packages
summary_table_language = data.frame(n_row = numeric(),mean = numeric(),sd = numeric(),
                           package = character())
for(n_row in seq(100,1000,100))
{
  for(p in c("TDApplied","persim"))
  {
    result = data.frame(n_row = n_row,
                        mean = mean(runtimes_language[which(runtimes_language$n_row == n_row
                                                            & runtimes_language$package == p),
                                                      3]),
                        sd = sd(runtimes_language[which(runtimes_language$n_row == n_row 
                                                        & runtimes_language$package == p),
                                                  3]),
                        package = p)
    summary_table_language = rbind(summary_table_language,result)
  }
}
# plot table
plot(summary_table_language$n_row[summary_table_language$package=="TDApplied"],
     summary_table_language$mean[summary_table_language$package=="persim"], type="b",
     xlim=range(summary_table_language$n_row),
     ylim=range(0,summary_table_language$mean+1.96*summary_table_language$sd/sqrt(10)),
     xlab = "Points in shape",ylab = "Mean execution time (sec)")
lines(summary_table_language$n_row[summary_table_language$package=="TDApplied"],
      summary_table_language$mean[summary_table_language$package=="TDApplied"], 
      col="red", type="b")
lines(summary_table_language$n_row[summary_table_language$package=="persim"],
      summary_table_language$mean[summary_table_language$package=="persim"], 
      col="black", type="b")
legend(x = 200,y = 20,legend = c("TDApplied","persim"),
       col = c("red","black"),lty = c(1,1),cex = 0.8)
arrows(summary_table_language$n_row[summary_table_language$package == "TDApplied"], 
       summary_table_language$mean[summary_table_language$package == "TDApplied"]
       -1.96*summary_table_language$sd[summary_table_language$package == "TDApplied"]/sqrt(10),
       summary_table_language$n_row[summary_table_language$package == "TDApplied"], 
       summary_table_language$mean[summary_table_language$package == "TDApplied"]
       +1.96*summary_table_language$sd[summary_table_language$package == "TDApplied"]/sqrt(10), 
       length=0.05, angle=90, code=3,col = "red")
arrows(summary_table_language$n_row[summary_table_language$package == "persim"], 
       summary_table_language$mean[summary_table_language$package == "persim"]
       -1.96*summary_table_language$sd[summary_table_language$package == "persim"]/sqrt(10),
       summary_table_language$n_row[summary_table_language$package == "persim"], 
       summary_table_language$mean[summary_table_language$package == "persim"]
       +1.96*summary_table_language$sd[summary_table_language$package == "persim"]/sqrt(10), 
       length=0.05, angle=90, code=3,col = "black")

#### PyH vs. calculate_homology ####

# generate persistence diagrams from circles with  100,200,...,1000 data points.
runtimes_circle <- data.frame(n_row = numeric(),package = character(),time_in_sec = numeric())
for(n_row in seq(100,1000,100)){
  
  for(iteration in 1:10)
  {
    # simulate a circle
    circ <- TDA::circleUnif(n = n_row)
    
    # compute their wasserstein distances in all dimensions and benchmark
    start_time_TDApplied = Sys.time()
    phom_TDApplied <- PyH(circ,maxdim = 1,thresh = 1,ripser = ripser)
    end_time_TDApplied = Sys.time()
    time_diff_TDApplied = as.numeric(end_time_TDApplied - start_time_TDApplied,units = "secs")
    
    start_time_TDAstats = Sys.time()
    start_time_TDAstats = Sys.time()
    phom_TDAstats <- TDAstats::calculate_homology(circ,threshold = 1)
    end_time_TDAstats = Sys.time()
    time_diff_TDAstats = as.numeric(end_time_TDAstats - start_time_TDAstats,units = "secs")
    
    runtimes_circle = rbind(runtimes_circle,data.frame(n_row = n_row,
                                                       package = "TDApplied",
                                                       time_in_sec = time_diff_TDApplied))
    runtimes_circle = rbind(runtimes_circle,data.frame(n_row = n_row,
                                                       package = "TDAstats",
                                                       time_in_sec = time_diff_TDAstats))
    
  }
  print(paste0("Done ",n_row," rows"))
  
}
# compute means and sd's at each value of rows for both packages
summary_table_circle = data.frame(n_row = numeric(),mean = numeric(),sd = numeric(),
                                  package = character())
for(n_row in seq(100,1000,100))
{
  for(p in c("TDApplied","TDAstats"))
  {
    result = data.frame(n_row = n_row,
                        mean = mean(runtimes_circle[which(runtimes_circle$n_row == n_row
                                                          & runtimes_circle$package == p),
                                                    3]),
                        sd = sd(runtimes_circle[which(runtimes_circle$n_row == n_row 
                                                      & runtimes_circle$package == p),
                                                3]),
                        package = p)
    summary_table_circle = rbind(summary_table_circle,result)
  }
}
# generate persistence diagrams from Tori with  100,200,...,1000 data points.
runtimes_torus <- data.frame(n_row = numeric(),package = character(),time_in_sec = numeric())
for(n_row in seq(100,1000,100)){
  
  for(iteration in 1:10)
  {
    # simulate a torus
    torus <- TDA::torusUnif(n = n_row,a = 0.25,c = 0.75)
    
    # compute their wasserstein distances in all dimensions and benchmark
    start_time_TDApplied = Sys.time()
    phom_TDApplied <- PyH(torus,maxdim = 2,thresh = 1,ripser = ripser)
    end_time_TDApplied = Sys.time()
    time_diff_TDApplied = as.numeric(end_time_TDApplied - start_time_TDApplied,units = "secs")
    
    start_time_TDAstats = Sys.time()
    phom_TDAstats <- TDAstats::calculate_homology(torus,threshold = 1,dim = 2)
    end_time_TDAstats = Sys.time()
    time_diff_TDAstats = as.numeric(end_time_TDAstats - start_time_TDAstats,units = "secs")
    
    runtimes_torus = rbind(runtimes_torus,data.frame(n_row = n_row,
                                                     package = "TDApplied",
                                                     time_in_sec = time_diff_TDApplied))
    runtimes_torus = rbind(runtimes_torus,data.frame(n_row = n_row,
                                                     package = "TDAstats",
                                                     time_in_sec = time_diff_TDAstats))
    
  }
  print(paste0("Done ",n_row," rows"))
  
}
# compute means and sd's at each value of rows for both packages
summary_table_torus = data.frame(n_row = numeric(),mean = numeric(),sd = numeric(),
                                 package = character())
for(n_row in seq(100,1000,100))
{
  for(p in c("TDApplied","TDAstats"))
  {
    result = data.frame(n_row = n_row,
                        mean = mean(runtimes_torus[which(runtimes_torus$n_row == n_row
                                                         & runtimes_torus$package == p),
                                                   3]),
                        sd = sd(runtimes_torus[which(runtimes_torus$n_row == n_row 
                                                     & runtimes_torus$package == p),
                                               3]),
                        package = p)
    summary_table_torus = rbind(summary_table_torus,result)
  }
}
# generate persistence diagrams from spheres with  100,200,...,1000 data points.
runtimes_sphere <- data.frame(n_row = numeric(),package = character(),time_in_sec = numeric())
for(n_row in seq(100,1000,100)){
  
  for(iteration in 1:10)
  {
    # simulate a sphere
    circ <- TDA::sphereUnif(n = n_row,d = 2)
    
    # compute their wasserstein distances in all dimensions and benchmark
    start_time_TDApplied = Sys.time()
    phom_TDApplied <- PyH(circ,maxdim = 2,thresh = 1,ripser = ripser)
    end_time_TDApplied = Sys.time()
    time_diff_TDApplied = as.numeric(end_time_TDApplied - start_time_TDApplied,units = "secs")
    
    start_time_TDAstats = Sys.time()
    start_time_TDAstats = Sys.time()
    phom_TDAstats <- TDAstats::calculate_homology(circ,threshold = 1,dim = 2)
    end_time_TDAstats = Sys.time()
    time_diff_TDAstats = as.numeric(end_time_TDAstats - start_time_TDAstats,units = "secs")
    
    runtimes_sphere = rbind(runtimes_sphere,data.frame(n_row = n_row,
                                                       package = "TDApplied",
                                                       time_in_sec = time_diff_TDApplied))
    runtimes_sphere = rbind(runtimes_sphere,data.frame(n_row = n_row,
                                                       package = "TDAstats",
                                                       time_in_sec = time_diff_TDAstats))
    
  }
  print(paste0("Done ",n_row," rows"))
  
}
# compute means and sd's at each value of rows for both packages
summary_table_sphere = data.frame(n_row = numeric(),mean = numeric(),sd = numeric(),
                                  package = character())
for(n_row in seq(100,1000,100))
{
  for(p in c("TDApplied","TDAstats"))
  {
    result = data.frame(n_row = n_row,
                        mean = mean(runtimes_sphere[which(runtimes_sphere$n_row == n_row
                                                          & runtimes_sphere$package == p),
                                                    3]),
                        sd = sd(runtimes_sphere[which(runtimes_sphere$n_row == n_row 
                                                      & runtimes_sphere$package == p),
                                                3]),
                        package = p)
    summary_table_sphere = rbind(summary_table_sphere,result)
  }
}

plot(summary_table_circle$n_row[summary_table_circle$package=="TDAstats"], 
     summary_table_circle$mean[summary_table_circle$package=="TDAstats"], 
     type="b",
     xlim=range(summary_table_circle$n_row),
     ylim=range(0,summary_table_circle$mean+1.96*summary_table_circle$sd/sqrt(10)),
     xlab = "Points in shape",ylab = "Mean execution time (sec)",
     main = "Circles")
lines(summary_table_circle$n_row[summary_table_circle$package=="TDApplied"],
      summary_table_circle$mean[summary_table_circle$package=="TDApplied"], 
      col=2, type="b")
legend(x = 200,y = 1.5,legend = c("TDApplied","TDAstats"),
       col = c("red","black"),lty = c(1,1),cex = 0.8)
arrows(summary_table_circle$n_row[summary_table_circle$package == "TDApplied"], 
       summary_table_circle$mean[summary_table_circle$package == "TDApplied"]
       -1.96*summary_table_circle$sd[summary_table_circle$package == "TDApplied"]/sqrt(10),
       summary_table_circle$n_row[summary_table_circle$package == "TDApplied"], 
       summary_table_circle$mean[summary_table_circle$package == "TDApplied"]
       +1.96*summary_table_circle$sd[summary_table_circle$package == "TDApplied"]/sqrt(10), 
       length=0.05, angle=90, code=3,col = "red")
arrows(summary_table_circle$n_row[summary_table_circle$package == "TDAstats"], 
       summary_table_circle$mean[summary_table_circle$package == "TDAstats"]
       -1.96*summary_table_circle$sd[summary_table_circle$package == "TDAstats"]/sqrt(10), 
       summary_table_circle$n_row[summary_table_circle$package == "TDAstats"], 
       summary_table_circle$mean[summary_table_circle$package == "TDAstats"]
       +1.96*summary_table_circle$sd[summary_table_circle$package == "TDAstats"]/sqrt(10), 
       length=0.05, angle=90, code=3,col = "black")
plot(summary_table_torus$n_row[summary_table_torus$package=="TDAstats"], 
     summary_table_torus$mean[summary_table_torus$package=="TDAstats"], 
     type="b",
     xlim=range(summary_table_torus$n_row),
     ylim=range(0,summary_table_torus$mean+1.96*summary_table_torus$sd/sqrt(10)),
     xlab = "Points in shape",ylab = "Mean execution time (sec)",
     main = "Tori")
lines(summary_table_torus$n_row[summary_table_torus$package=="TDApplied"],
      summary_table_torus$mean[summary_table_torus$package=="TDApplied"], 
      col=2, type="b")
legend(x = 200,y = 120,legend = c("TDApplied","TDAstats"),
       col = c("red","black"),lty = c(1,1),cex = 0.8)
arrows(summary_table_torus$n_row[summary_table_torus$package == "TDApplied"], 
       summary_table_torus$mean[summary_table_torus$package == "TDApplied"]
       -1.96*summary_table_torus$sd[summary_table_torus$package == "TDApplied"]/sqrt(10),
       summary_table_torus$n_row[summary_table_torus$package == "TDApplied"], 
       summary_table_torus$mean[summary_table_torus$package == "TDApplied"]
       +1.96*summary_table_torus$sd[summary_table_torus$package == "TDApplied"]/sqrt(10), 
       length=0.05, angle=90, code=3,col = "red")
arrows(summary_table_torus$n_row[summary_table_torus$package == "TDAstats"], 
       summary_table_torus$mean[summary_table_torus$package == "TDAstats"]
       -1.96*summary_table_torus$sd[summary_table_torus$package == "TDAstats"]/sqrt(10), 
       summary_table_torus$n_row[summary_table_torus$package == "TDAstats"], 
       summary_table_torus$mean[summary_table_torus$package == "TDAstats"]
       +1.96*summary_table_torus$sd[summary_table_torus$package == "TDAstats"]/sqrt(10), 
       length=0.05, angle=90, code=3,col = "black")
plot(summary_table_sphere$n_row[summary_table_sphere$package=="TDAstats"], 
     summary_table_sphere$mean[summary_table_sphere$package=="TDAstats"], 
     type="b",
     xlim=range(summary_table_sphere$n_row),
     ylim=range(0,summary_table_sphere$mean+1.96*summary_table_sphere$sd/sqrt(10)),
     xlab = "Points in shape",ylab = "Mean execution time (sec)",
     main = "Spheres")
lines(summary_table_sphere$n_row[summary_table_sphere$package=="TDApplied"],
      summary_table_sphere$mean[summary_table_sphere$package=="TDApplied"], 
      col=2, type="b")
legend(x = 200,y = 45,legend = c("TDApplied","TDAstats"),
       col = c("red","black"),lty = c(1,1),cex = 0.8)
arrows(summary_table_sphere$n_row[summary_table_sphere$package == "TDApplied"], 
       summary_table_sphere$mean[summary_table_sphere$package == "TDApplied"]
       -1.96*summary_table_sphere$sd[summary_table_sphere$package == "TDApplied"]/sqrt(10),
       summary_table_sphere$n_row[summary_table_sphere$package == "TDApplied"], 
       summary_table_sphere$mean[summary_table_sphere$package == "TDApplied"]
       +1.96*summary_table_sphere$sd[summary_table_sphere$package == "TDApplied"]/sqrt(10), 
       length=0.05, angle=90, code=3,col = "red")
arrows(summary_table_sphere$n_row[summary_table_sphere$package == "TDAstats"], 
       summary_table_sphere$mean[summary_table_sphere$package == "TDAstats"]
       -1.96*summary_table_sphere$sd[summary_table_sphere$package == "TDAstats"]/sqrt(10), 
       summary_table_sphere$n_row[summary_table_sphere$package == "TDAstats"], 
       summary_table_sphere$mean[summary_table_sphere$package == "TDAstats"]
       +1.96*summary_table_sphere$sd[summary_table_sphere$package == "TDAstats"]/sqrt(10), 
       length=0.05, angle=90, code=3,col = "black")

#### bootstrap_persistence_thresholds vs. bootstrapDiagram ####

# currently unfinished, to be done for the next version of TDApplied

# # generate persistence diagrams from circles with  100,200,...,1000 data points.
# benchmark_bootstrap_circle <- data.frame(n_row = numeric(),package = character(),time_in_sec = numeric(),Mb = numeric())
# for(n_row in seq(100,1000,100)){
#   
#   for(iteration in 1:10)
#   {
#     # simulate a circle
#     circ <- TDA::circleUnif(n = n_row)
#     
#     # bootstrap and benchmark
#     
#     # TDApplied
#     bootstrap_circ_TDApplied = as.data.frame(bench::mark(bootstrap_persistence_thresholds(X = circ,FUN = "calculate_homology",maxdim = 1,thresh = 1,global_threshold = F)))
#     
#     # TDAstats
#     # following the example on the help page of bootstrapDiagram
#     bootstrap_circ_TDA_dim_0 = as.data.frame(bench::mark(TDA::bootstrapDiagram(X = circ,FUN = TDA::distFct,lim = cbind(c(-1,1),c(-1,1)),by = 0.05,sublevel = T,B = 30,alpha = 0.05,dimension = 0,printProgress = F)))
#     
#     bootstrap_circ_TDA_dim_1 = as.data.frame(bench::mark(TDA::bootstrapDiagram(X = circ,FUN = TDA::distFct,lim = cbind(c(-1,1),c(-1,1)),by = 0.05,sublevel = T,B = 30,alpha = 0.05,dimension = 1,printProgress = F)))
#     
#     benchmark_bootstrap_circle = rbind(benchmark_bootstrap_circle,data.frame(n_row = n_row,
#                                                                              package = "TDApplied",
#                                                                              time_in_sec = as.numeric(bootstrap_circ_TDApplied$min,units = "sec"),
#                                                                              Mb = as.numeric(bootstrap_circ_TDApplied$mem_alloc)/1000000))
#     benchmark_bootstrap_circle = rbind(benchmark_bootstrap_circle,data.frame(n_row = n_row,
#                                                                              package = "TDA",
#                                                                              time_in_sec = as.numeric(bootstrap_circ_TDA_dim_0$min + bootstrap_circ_TDA_dim_1$min,units = "sec"),
#                                                                              Mb = as.numeric(bootstrap_circ_TDA_dim_0$mem_alloc + bootstrap_circ_TDA_dim_1$mem_alloc)/1000000))
#     
#   }
#   
#   print(paste0("Done ",n_row," rows"))
#   
# }
# 
# # compute means and sd's at each value of rows for both packages
# summary_table_bootstrap_circle = data.frame(n_row = numeric(),mean_time = numeric(),sd_time = numeric(),
#                                             mean_memory = numeric(),sd_memory = numeric(),package = character())
# for(n_row in seq(100,1000,100))
# {
#   for(p in c("TDApplied","TDA"))
#   {
#     df_subset = benchmark_circle[which(benchmark_circle$n_row == n_row & benchmark_circle$package == p),]
#     summary_table_bootstrap_circle = rbind(summary_table_bootstrap_circle,data.frame(n_row = n_row,
#                                                                                      mean_time = mean(df_subset$time_in_sec),
#                                                                                      sd_time = sd(df_subset$time_in_sec),
#                                                                                      mean_memory = mean(df_subset$Mb),
#                                                                                      sd_memory = sd(df_subset$Mb),
#                                                                                      package = p))
#   }
# }
# 
# # generate persistence diagrams from tori with  100,200,...,1000 data points.
# benchmark_bootstrap_torus <- data.frame(n_row = numeric(),package = character(),time_in_sec = numeric(),Mb = numeric())
# for(n_row in seq(100,1000,100)){
#   
#   for(iteration in 1:10)
#   {
#     # simulate a torus
#     torus <- TDA::torusUnif(n = n_row,a = 0.25,c = 0.75)
#     
#     # bootstrap and benchmark
#     
#     # TDApplied
#     bootstrap_torus_TDApplied = as.data.frame(bench::mark(bootstrap_persistence_thresholds(X = torus,FUN = "calculate_homology",maxdim = 2,thresh = 1,global_threshold = F)))
#     
#     # TDAstats
#     # following the example on the help page of bootstrapDiagram
#     bootstrap_torus_TDA_dim_0 = as.data.frame(bench::mark(TDA::bootstrapDiagram(X = torus,FUN = TDA::distFct,lim = cbind(c(-1,1),c(-1,1),c(-1,1)),by = 0.05,sublevel = T,B = 30,alpha = 0.05,dimension = 0,printProgress = F)))
#     
#     bootstrap_torus_TDA_dim_1 = as.data.frame(bench::mark(TDA::bootstrapDiagram(X = torus,FUN = TDA::distFct,lim = cbind(c(-1,1),c(-1,1),c(-1,1)),by = 0.05,sublevel = T,B = 30,alpha = 0.05,dimension = 1,printProgress = F)))
#     
#     bootstrap_torus_TDA_dim_2 = as.data.frame(bench::mark(TDA::bootstrapDiagram(X = torus,FUN = TDA::distFct,lim = cbind(c(-1,1),c(-1,1),c(-1,1)),by = 0.05,sublevel = T,B = 30,alpha = 0.05,dimension = 2,printProgress = F)))
#     
#     benchmark_bootstrap_torus = rbind(benchmark_bootstrap_torus,data.frame(n_row = n_row,
#                                                                              package = "TDApplied",
#                                                                              time_in_sec = as.numeric(bootstrap_torus_TDApplied$min,units = "sec"),
#                                                                              Mb = as.numeric(bootstrap_torus_TDApplied$mem_alloc)/1000000))
#     benchmark_bootstrap_torus = rbind(benchmark_bootstrap_torus,data.frame(n_row = n_row,
#                                                                              package = "TDA",
#                                                                              time_in_sec = as.numeric(bootstrap_torus_TDA_dim_0$min + bootstrap_torus_TDA_dim_1$min + bootstrap_torus_TDA_dim_2$min,units = "sec"),
#                                                                              Mb = as.numeric(bootstrap_torus_TDA_dim_0$mem_alloc + bootstrap_torus_TDA_dim_1$mem_alloc + bootstrap_torus_TDA_dim_2$mem_alloc)/1000000))
#     
#   }
#   
#   print(paste0("Done ",n_row," rows"))
#   
# }
# 
# # compute means and sd's at each value of rows for both packages
# summary_table_bootstrap_torus = data.frame(n_row = numeric(),mean_time = numeric(),sd_time = numeric(),
#                                             mean_memory = numeric(),sd_memory = numeric(),package = character())
# for(n_row in seq(100,1000,100))
# {
#   for(p in c("TDApplied","TDA"))
#   {
#     df_subset = benchmark_torus[which(benchmark_torus$n_row == n_row & benchmark_torus$package == p),]
#     summary_table_bootstrap_torus = rbind(summary_table_bootstrap_torus,data.frame(n_row = n_row,
#                                                                                      mean_time = mean(df_subset$time_in_sec),
#                                                                                      sd_time = sd(df_subset$time_in_sec),
#                                                                                      mean_memory = mean(df_subset$Mb),
#                                                                                      sd_memory = sd(df_subset$Mb),
#                                                                                      package = p))
#   }
# }
# 
# # generate persistence diagrams from tori with  100,200,...,1000 data points.
# benchmark_bootstrap_sphere <- data.frame(n_row = numeric(),package = character(),time_in_sec = numeric(),Mb = numeric())
# for(n_row in seq(100,1000,100)){
#   
#   for(iteration in 1:10)
#   {
#     # simulate a sphere
#     sphere <- TDA::sphereUnif(n = n_row,d = 2)
#     
#     # bootstrap and benchmark
#     
#     # TDApplied
#     bootstrap_sphere_TDApplied = as.data.frame(bench::mark(bootstrap_persistence_thresholds(X = sphere,FUN = "calculate_homology",maxdim = 2,thresh = 1,global_threshold = F)))
#     
#     # TDAstats
#     # following the example on the help page of bootstrapDiagram
#     bootstrap_sphere_TDA_dim_0 = as.data.frame(bench::mark(TDA::bootstrapDiagram(X = sphere,FUN = TDA::distFct,lim = cbind(c(-1,1),c(-1,1),c(-1,1)),by = 0.05,sublevel = T,B = 30,alpha = 0.05,dimension = 0,printProgress = F)))
#     
#     bootstrap_sphere_TDA_dim_1 = as.data.frame(bench::mark(TDA::bootstrapDiagram(X = sphere,FUN = TDA::distFct,lim = cbind(c(-1,1),c(-1,1)),by = 0.05,sublevel = T,B = 30,alpha = 0.05,dimension = 1,printProgress = F)))
#     
#     bootstrap_sphere_TDA_dim_2 = as.data.frame(bench::mark(TDA::bootstrapDiagram(X = sphere,FUN = TDA::distFct,lim = cbind(c(-1,1),c(-1,1)),by = 0.05,sublevel = T,B = 30,alpha = 0.05,dimension = 2,printProgress = F)))
#     
#     benchmark_bootstrap_sphere = rbind(benchmark_bootstrap_sphere,data.frame(n_row = n_row,
#                                                                            package = "TDApplied",
#                                                                            time_in_sec = as.numeric(bootstrap_sphere_TDApplied$min,units = "sec"),
#                                                                            Mb = as.numeric(bootstrap_sphere_TDApplied$mem_alloc)/1000000))
#     benchmark_bootstrap_sphere = rbind(benchmark_bootstrap_sphere,data.frame(n_row = n_row,
#                                                                            package = "TDA",
#                                                                            time_in_sec = as.numeric(bootstrap_sphere_TDA_dim_0$min + bootstrap_sphere_TDA_dim_1$min + bootstrap_sphere_TDA_dim_2$min,units = "sec"),
#                                                                            Mb = as.numeric(bootstrap_sphere_TDA_dim_0$mem_alloc + bootstrap_sphere_TDA_dim_1$mem_alloc + bootstrap_sphere_TDA_dim_2$mem_alloc)/1000000))
#     
#   }
#   
#   print(paste0("Done ",n_row," rows"))
#   
# }
# 
# # compute means and sd's at each value of rows for both packages
# summary_table_bootstrap_sphere = data.frame(n_row = numeric(),mean_time = numeric(),sd_time = numeric(),
#                                             mean_memory = numeric(),sd_memory = numeric(),package = character())
# for(n_row in seq(100,1000,100))
# {
#   for(p in c("TDApplied","TDA"))
#   {
#     df_subset = benchmark_sphere[which(benchmark_sphere$n_row == n_row & benchmark_sphere$package == p),]
#     summary_table_bootstrap_sphere = rbind(summary_table_bootstrap_sphere,data.frame(n_row = n_row,
#                                                                                    mean_time = mean(df_subset$time_in_sec),
#                                                                                    sd_time = sd(df_subset$time_in_sec),
#                                                                                    mean_memory = mean(df_subset$Mb),
#                                                                                    sd_memory = sd(df_subset$Mb),
#                                                                                    package = p))
#   }
# }
# 
# # plot time benchmarking for boostrapping circles
# plot(summary_table_bootstrap_circle$n_row[summary_table_bootstrap_circle$package=="TDAstats"], 
#      summary_table_bootstrap_circle$mean_time[summary_table_bootstrap_circle$package=="TDAstats"], 
#      type="b",
#      xlim=range(summary_table_bootstrap_circle$n_row),
#      ylim=range(0,summary_table_bootstrap_circle$mean_time+1.96*summary_table_bootstrap_circle$sd_time/sqrt(10)),
#      xlab = "Points in shape",ylab = "Mean execution time (sec)",
#      main = "Circles")
# lines(summary_table_bootstrap_circle$n_row[summary_table_bootstrap_circle$package=="TDApplied"],
#       summary_table_bootstrap_circle$mean_time[summary_table_bootstrap_circle$package=="TDApplied"], 
#       col=2, type="b")
# legend(x = 200,y = 1.5,legend = c("TDApplied","TDAstats"),
#        col = c("red","black"),lty = c(1,1),cex = 0.8)
# arrows(summary_table_bootstrap_circle$n_row[summary_table_bootstrap_circle$package == "TDApplied"], 
#        summary_table_bootstrap_circle$mean_time[summary_table_bootstrap_circle$package == "TDApplied"]
#        -1.96*summary_table_bootstrap_circle$sd_time[summary_table_bootstrap_circle$package == "TDApplied"]/sqrt(10),
#        summary_table_bootstrap_circle$n_row[summary_table_bootstrap_circle$package == "TDApplied"], 
#        summary_table_bootstrap_circle$mean_time[summary_table_bootstrap_circle$package == "TDApplied"]
#        +1.96*summary_table_bootstrap_circle$sd_time[summary_table_bootstrap_circle$package == "TDApplied"]/sqrt(10), 
#        length=0.05, angle=90, code=3,col = "red")
# arrows(summary_table_bootstrap_circle$n_row[summary_table_bootstrap_circle$package == "TDAstats"], 
#        summary_table_bootstrap_circle$mean_time[summary_table_bootstrap_circle$package == "TDAstats"]
#        -1.96*summary_table_bootstrap_circle$sd_time[summary_table_bootstrap_circle$package == "TDAstats"]/sqrt(10), 
#        summary_table_bootstrap_circle$n_row[summary_table_bootstrap_circle$package == "TDAstats"], 
#        summary_table_bootstrap_circle$mean_time[summary_table_bootstrap_circle$package == "TDAstats"]
#        +1.96*summary_table_bootstrap_circle$sd_time[summary_table_bootstrap_circle$package == "TDAstats"]/sqrt(10), 
#        length=0.05, angle=90, code=3,col = "black")
# 
# # plot memory benchmarking for bootstrapping circles
# plot(summary_table_bootstrap_circle$n_row[summary_table_bootstrap_circle$package=="TDAstats"], 
#      summary_table_bootstrap_circle$mean_memory[summary_table_bootstrap_circle$package=="TDAstats"], 
#      type="b",
#      xlim=range(summary_table_bootstrap_circle$n_row),
#      ylim=range(0,summary_table_bootstrap_circle$mean_memory+1.96*summary_table_bootstrap_circle$sd_memory/sqrt(10)),
#      xlab = "Points in shape",ylab = "Mean memory allocated (Mb)",
#      main = "Circles")
# lines(summary_table_bootstrap_circle$n_row[summary_table_bootstrap_circle$package=="TDApplied"],
#       summary_table_bootstrap_circle$mean_memory[summary_table_bootstrap_circle$package=="TDApplied"], 
#       col=2, type="b")
# legend(x = 200,y = 1.5,legend = c("TDApplied","TDAstats"),
#        col = c("red","black"),lty = c(1,1),cex = 0.8)
# arrows(summary_table_bootstrap_circle$n_row[summary_table_bootstrap_circle$package == "TDApplied"], 
#        summary_table_bootstrap_circle$mean_memory[summary_table_bootstrap_circle$package == "TDApplied"]
#        -1.96*summary_table_bootstrap_circle$sd_memory[summary_table_bootstrap_circle$package == "TDApplied"]/sqrt(10),
#        summary_table_bootstrap_circle$n_row[summary_table_bootstrap_circle$package == "TDApplied"], 
#        summary_table_bootstrap_circle$mean_memory[summary_table_bootstrap_circle$package == "TDApplied"]
#        +1.96*summary_table_bootstrap_circle$sd_memory[summary_table_bootstrap_circle$package == "TDApplied"]/sqrt(10), 
#        length=0.05, angle=90, code=3,col = "red")
# arrows(summary_table_bootstrap_circle$n_row[summary_table_bootstrap_circle$package == "TDAstats"], 
#        summary_table_bootstrap_circle$mean_memory[summary_table_bootstrap_circle$package == "TDAstats"]
#        -1.96*summary_table_bootstrap_circle$sd_memory[summary_table_bootstrap_circle$package == "TDAstats"]/sqrt(10), 
#        summary_table_bootstrap_circle$n_row[summary_table_bootstrap_circle$package == "TDAstats"], 
#        summary_table_bootstrap_circle$mean_memory[summary_table_bootstrap_circle$package == "TDAstats"]
#        +1.96*summary_table_bootstrap_circle$sd_memory[summary_table_bootstrap_circle$package == "TDAstats"]/sqrt(10), 
#        length=0.05, angle=90, code=3,col = "black")
# 
# # plot time benchmarking for bootstrapping torus
# plot(summary_table_bootstrap_torus$n_row[summary_table_bootstrap_torus$package=="TDAstats"], 
#      summary_table_bootstrap_torus$mean_time[summary_table_bootstrap_torus$package=="TDAstats"], 
#      type="b",
#      xlim=range(summary_table_bootstrap_torus$n_row),
#      ylim=range(0,summary_table_bootstrap_torus$mean_time+1.96*summary_table_bootstrap_torus$sd_time/sqrt(10)),
#      xlab = "Points in shape",ylab = "Mean execution time (sec)",
#      main = "Tori")
# lines(summary_table_bootstrap_torus$n_row[summary_table_bootstrap_torus$package=="TDApplied"],
#       summary_table_bootstrap_torus$mean_time[summary_table_bootstrap_torus$package=="TDApplied"], 
#       col=2, type="b")
# legend(x = 200,y = 120,legend = c("TDApplied","TDAstats"),
#        col = c("red","black"),lty = c(1,1),cex = 0.8)
# arrows(summary_table_bootstrap_torus$n_row[summary_table_bootstrap_torus$package == "TDApplied"], 
#        summary_table_bootstrap_torus$mean_time[summary_table_bootstrap_torus$package == "TDApplied"]
#        -1.96*summary_table_bootstrap_torus$sd_time[summary_table_bootstrap_torus$package == "TDApplied"]/sqrt(10),
#        summary_table_bootstrap_torus$n_row[summary_table_bootstrap_torus$package == "TDApplied"], 
#        summary_table_bootstrap_torus$mean_time[summary_table_bootstrap_torus$package == "TDApplied"]
#        +1.96*summary_table_bootstrap_torus$sd_time[summary_table_bootstrap_torus$package == "TDApplied"]/sqrt(10), 
#        length=0.05, angle=90, code=3,col = "red")
# arrows(summary_table_bootstrap_torus$n_row[summary_table_bootstrap_torus$package == "TDAstats"], 
#        summary_table_bootstrap_torus$mean_time[summary_table_bootstrap_torus$package == "TDAstats"]
#        -1.96*summary_table_bootstrap_torus$sd_time[summary_table_bootstrap_torus$package == "TDAstats"]/sqrt(10), 
#        summary_table_bootstrap_torus$n_row[summary_table_bootstrap_torus$package == "TDAstats"], 
#        summary_table_bootstrap_torus$mean_time[summary_table_bootstrap_torus$package == "TDAstats"]
#        +1.96*summary_table_bootstrap_torus$sd_time[summary_table_bootstrap_torus$package == "TDAstats"]/sqrt(10), 
#        length=0.05, angle=90, code=3,col = "black")
# 
# # plot memory benchmarking for bootstrapping torus
# plot(summary_table_bootstrap_torus$n_row[summary_table_bootstrap_torus$package=="TDAstats"], 
#      summary_table_bootstrap_torus$mean_memory[summary_table_bootstrap_torus$package=="TDAstats"], 
#      type="b",
#      xlim=range(summary_table_bootstrap_torus$n_row),
#      ylim=range(0,summary_table_bootstrap_torus$mean_memory+1.96*summary_table_bootstrap_torus$sd_memory/sqrt(10)),
#      xlab = "Points in shape",ylab = "Mean memory allocated (Mb)",
#      main = "Tori")
# lines(summary_table_bootstrap_torus$n_row[summary_table_bootstrap_torus$package=="TDApplied"],
#       summary_table_bootstrap_torus$mean_memory[summary_table_bootstrap_torus$package=="TDApplied"], 
#       col=2, type="b")
# legend(x = 200,y = 120,legend = c("TDApplied","TDAstats"),
#        col = c("red","black"),lty = c(1,1),cex = 0.8)
# arrows(summary_table_bootstrap_torus$n_row[summary_table_bootstrap_torus$package == "TDApplied"], 
#        summary_table_bootstrap_torus$mean_memory[summary_table_bootstrap_torus$package == "TDApplied"]
#        -1.96*summary_table_bootstrap_torus$sd_memory[summary_table_bootstrap_torus$package == "TDApplied"]/sqrt(10),
#        summary_table_bootstrap_torus$n_row[summary_table_bootstrap_torus$package == "TDApplied"], 
#        summary_table_bootstrap_torus$mean_memory[summary_table_bootstrap_torus$package == "TDApplied"]
#        +1.96*summary_table_bootstrap_torus$sd_memory[summary_table_bootstrap_torus$package == "TDApplied"]/sqrt(10), 
#        length=0.05, angle=90, code=3,col = "red")
# arrows(summary_table_bootstrap_torus$n_row[summary_table_bootstrap_torus$package == "TDAstats"], 
#        summary_table_bootstrap_torus$mean_memory[summary_table_bootstrap_torus$package == "TDAstats"]
#        -1.96*summary_table_bootstrap_torus$sd_memory[summary_table_bootstrap_torus$package == "TDAstats"]/sqrt(10), 
#        summary_table_bootstrap_torus$n_row[summary_table_bootstrap_torus$package == "TDAstats"], 
#        summary_table_bootstrap_torus$mean_memory[summary_table_bootstrap_torus$package == "TDAstats"]
#        +1.96*summary_table_bootstrap_torus$sd_memory[summary_table_bootstrap_torus$package == "TDAstats"]/sqrt(10), 
#        length=0.05, angle=90, code=3,col = "black")
# 
# # plot time benchmarking for bootstrapping sphere
# plot(summary_table_bootstrap_sphere$n_row[summary_table_bootstrap_sphere$package=="TDAstats"], 
#      summary_table_bootstrap_sphere$mean_time[summary_table_bootstrap_sphere$package=="TDAstats"], 
#      type="b",
#      xlim=range(summary_table_bootstrap_sphere$n_row),
#      ylim=range(0,summary_table_bootstrap_sphere$mean_time+1.96*summary_table_bootstrap_sphere$sd_time/sqrt(10)),
#      xlab = "Points in shape",ylab = "Mean execution time (sec)",
#      main = "Spheres")
# lines(summary_table_bootstrap_sphere$n_row[summary_table_bootstrap_sphere$package=="TDApplied"],
#       summary_table_bootstrap_sphere$mean_time[summary_table_bootstrap_sphere$package=="TDApplied"], 
#       col=2, type="b")
# legend(x = 200,y = 45,legend = c("TDApplied","TDAstats"),
#        col = c("red","black"),lty = c(1,1),cex = 0.8)
# arrows(summary_table_bootstrap_sphere$n_row[summary_table_bootstrap_sphere$package == "TDApplied"], 
#        summary_table_bootstrap_sphere$mean_time[summary_table_bootstrap_sphere$package == "TDApplied"]
#        -1.96*summary_table_bootstrap_sphere$sd_time[summary_table_bootstrap_sphere$package == "TDApplied"]/sqrt(10),
#        summary_table_bootstrap_sphere$n_row[summary_table_bootstrap_sphere$package == "TDApplied"], 
#        summary_table_bootstrap_sphere$mean_time[summary_table_bootstrap_sphere$package == "TDApplied"]
#        +1.96*summary_table_bootstrap_sphere$sd_time[summary_table_bootstrap_sphere$package == "TDApplied"]/sqrt(10), 
#        length=0.05, angle=90, code=3,col = "red")
# arrows(summary_table_bootstrap_sphere$n_row[summary_table_bootstrap_sphere$package == "TDAstats"], 
#        summary_table_bootstrap_sphere$mean_time[summary_table_bootstrap_sphere$package == "TDAstats"]
#        -1.96*summary_table_bootstrap_sphere$sd_time[summary_table_bootstrap_sphere$package == "TDAstats"]/sqrt(10), 
#        summary_table_bootstrap_sphere$n_row[summary_table_bootstrap_sphere$package == "TDAstats"], 
#        summary_table_bootstrap_sphere$mean_time[summary_table_bootstrap_sphere$package == "TDAstats"]
#        +1.96*summary_table_bootstrap_sphere$sd_time[summary_table_bootstrap_sphere$package == "TDAstats"]/sqrt(10), 
#        length=0.05, angle=90, code=3,col = "black")
# 
# # plot memory benchmarking for bootstrapping sphere
# plot(summary_table_bootstrap_sphere$n_row[summary_table_bootstrap_sphere$package=="TDAstats"], 
#      summary_table_bootstrap_sphere$mean_memory[summary_table_bootstrap_sphere$package=="TDAstats"], 
#      type="b",
#      xlim=range(summary_table_bootstrap_sphere$n_row),
#      ylim=range(0,summary_table_bootstrap_sphere$mean_memory+1.96*summary_table_bootstrap_sphere$sd_memory/sqrt(10)),
#      xlab = "Points in shape",ylab = "Mean execution time (sec)",
#      main = "Spheres")
# lines(summary_table_bootstrap_sphere$n_row[summary_table_bootstrap_sphere$package=="TDApplied"],
#       summary_table_bootstrap_sphere$mean_memory[summary_table_bootstrap_sphere$package=="TDApplied"], 
#       col=2, type="b")
# legend(x = 200,y = 45,legend = c("TDApplied","TDAstats"),
#        col = c("red","black"),lty = c(1,1),cex = 0.8)
# arrows(summary_table_bootstrap_sphere$n_row[summary_table_bootstrap_sphere$package == "TDApplied"], 
#        summary_table_bootstrap_sphere$mean_memory[summary_table_bootstrap_sphere$package == "TDApplied"]
#        -1.96*summary_table_bootstrap_sphere$sd_memory[summary_table_bootstrap_sphere$package == "TDApplied"]/sqrt(10),
#        summary_table_bootstrap_sphere$n_row[summary_table_bootstrap_sphere$package == "TDApplied"], 
#        summary_table_bootstrap_sphere$mean_memory[summary_table_bootstrap_sphere$package == "TDApplied"]
#        +1.96*summary_table_bootstrap_sphere$sd_memory[summary_table_bootstrap_sphere$package == "TDApplied"]/sqrt(10), 
#        length=0.05, angle=90, code=3,col = "red")
# arrows(summary_table_bootstrap_sphere$n_row[summary_table_bootstrap_sphere$package == "TDAstats"], 
#        summary_table_bootstrap_sphere$mean_memory[summary_table_bootstrap_sphere$package == "TDAstats"]
#        -1.96*summary_table_bootstrap_sphere$sd_memory[summary_table_bootstrap_sphere$package == "TDAstats"]/sqrt(10), 
#        summary_table_bootstrap_sphere$n_row[summary_table_bootstrap_sphere$package == "TDAstats"], 
#        summary_table_bootstrap_sphere$mean_memory[summary_table_bootstrap_sphere$package == "TDAstats"]
#        +1.96*summary_table_bootstrap_sphere$sd_memory[summary_table_bootstrap_sphere$package == "TDAstats"]/sqrt(10), 
#        length=0.05, angle=90, code=3,col = "black")