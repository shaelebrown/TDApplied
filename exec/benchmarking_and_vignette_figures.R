
# benchmarking

#### diagram_distance vs. fast approximation ####
ripser <- import_ripser()

# generate persistence diagrams from spheres and tori with  100,200,...,1000 data points.
runtimes_approx <- data.frame(n_row = numeric(),method = character(),time_in_sec = numeric())
for(n_row in seq(100,1000,100)){
  
  for(iteration in 1:10)
  {
    # simulate a sphere and a torus
    sphere <- TDA::sphereUnif(n = n_row,d = 2)
    torus <- TDA::torusUnif(n = n_row,a = 0.25,c = 0.75)
    
    # compute their persistence diagrams
    sphere_diag <- PyH(sphere,maxdim = 2,thresh = 1,ripser = ripser)
    torus_diag <- PyH(torus,maxdim = 2,thresh = 1,ripser = ripser)
    
    # compute their Fisher information distances in all dimensions and benchmark
    start_time_exact = Sys.time()
    diagram_distance(sphere_diag,torus_diag,dim = 0,distance = "fisher",sigma = 1)
    diagram_distance(sphere_diag,torus_diag,dim = 1,distance = "fisher",sigma = 1)
    diagram_distance(sphere_diag,torus_diag,dim = 2,distance = "fisher",sigma = 1)
    end_time_exact = Sys.time()
    time_diff_exact = as.numeric(end_time_exact - start_time_exact,units = "secs")
    
    start_time_approximation = Sys.time()
    diagram_distance(sphere_diag,torus_diag,dim = 0,distance = "fisher",sigma = 1,rho = 0.001)
    diagram_distance(sphere_diag,torus_diag,dim = 1,distance = "fisher",sigma = 1,rho = 0.001)
    diagram_distance(sphere_diag,torus_diag,dim = 2,distance = "fisher",sigma = 1,rho = 0.001)
    end_time_approximation = Sys.time()
    time_diff_approximation = as.numeric(end_time_approximation - start_time_approximation,units = "secs")
    
    runtimes_approx = rbind(runtimes_approx,data.frame(n_row = n_row,
                                                       method = "Exact",
                                                       time_in_sec = time_diff_exact))
    runtimes_approx = rbind(runtimes_approx,data.frame(n_row = n_row,
                                                       method = "Approximation",
                                                       time_in_sec = time_diff_approximation))
    
  }
  print(paste0("Done ",n_row," rows"))
  
}
# compute means and sd's at each value of rows for both methods
summary_table_approx = data.frame(n_row = numeric(),mean = numeric(),sd = numeric(),
                                  method = character())
for(n_row in seq(100,1000,100))
{
  for(m in c("Exact","Approximation"))
  {
    result = data.frame(n_row = n_row,
                        mean = mean(runtimes_approx[which(runtimes_approx$n_row == n_row
                                                          & runtimes_approx$method == m),
                                                    3]),
                        sd = sd(runtimes_approx[which(runtimes_approx$n_row == n_row
                                                      & runtimes_approx$method == m),
                                                3]),
                        method = m)
    summary_table_approx = rbind(summary_table_approx,result)
  }
}

plot(summary_table_approx$n_row[summary_table_approx$method=="Exact"],
     summary_table_approx$mean[summary_table_approx$method=="Exact"],
     type="b",
     xlim=range(summary_table_approx$n_row),
     ylim=range(0,summary_table_approx$mean+1.96*summary_table_approx$sd/sqrt(10)),
     xlab = "Points in shape",ylab = "Mean execution time (sec)",
     main = "Approximation Benchmarking")
lines(summary_table_approx$n_row[summary_table_approx$method=="Approximation"],
      summary_table_approx$mean[summary_table_approx$method=="Approximation"],
      col=2, type="b")
legend(x = 200,y = 350,legend = c("Approximation","Exact"),
       col = c("red","black"),lty = c(1,1),cex = 0.8)
arrows(summary_table_approx$n_row[summary_table_approx$method == "Approximation"],
       summary_table_approx$mean[summary_table_approx$method == "Approximation"]
       -1.96*summary_table_approx$sd[summary_table_approx$method == "Approximation"]/sqrt(10),
       summary_table_approx$n_row[summary_table_approx$method == "Approximation"],
       summary_table_approx$mean[summary_table_approx$method == "Approximation"]
       +1.96*summary_table_approx$sd[summary_table_approx$method == "Approximation"]/sqrt(10),
       length=0.05, angle=90, code=3,col = "red")
arrows(summary_table_approx$n_row[summary_table_approx$method == "Exact"],
       summary_table_approx$mean[summary_table_approx$method == "Exact"]
       -1.96*summary_table_approx$sd[summary_table_approx$method == "Exact"]/sqrt(10),
       summary_table_approx$n_row[summary_table_approx$method == "Exact"],
       summary_table_approx$mean[summary_table_approx$method == "Exact"]
       +1.96*summary_table_approx$sd[summary_table_approx$method == "Exact"]/sqrt(10),
       length=0.05, angle=90, code=3,col = "black")

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
     xlab = "Points in shape",ylab = "Mean execution time (sec)",col = "darkgreen")
lines(summary_table_distance$n_row[summary_table_distance$package=="TDApplied"],
      summary_table_distance$mean[summary_table_distance$package=="TDApplied"], 
      col=2, type="b")
legend(x = 200,y = 2000,legend = c("TDApplied","TDA"),
       col = c("red","darkgreen"),lty = c(1,1),cex = 0.8)
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
       length=0.05, angle=90, code=3,col = "darkgreen")

#### diagram_distance vs persim's wasserstein ####

if(requireNamespace("reticulate",quietly = T) == T)
{
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
        col="darkorange", type="b")
  legend(x = 200,y = 20,legend = c("TDApplied","persim"),
         col = c("red","darkorange"),lty = c(1,1),cex = 0.8)
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
         length=0.05, angle=90, code=3,col = "darkorange") 
}

#### diagram_distance vs rgudhi's PersistenceFisherDistance ####

# generate persistence diagrams from Tori and spheres with  100,200,...,1000 data points.
runtimes_rgudhi <- data.frame(n_row = numeric(),package = character(),approx = logical(),time_in_sec = numeric())
dis_fish <- rgudhi::PersistenceFisherDistance$new() # default sigma is 1
# approx = reticulate::import("sklearn.kernel_approximation") # these lines are not working
# dis_fish_approx <- rgudhi::PersistenceFisherDistance$new(kernel_approx = approx$RBFSampler)
for(n_row in seq(100,1000,100)){
  for(iteration in 1:10)
  {
    # simulate pair of diagrams from the desired shapes
    diagram_torus = PyH(X = TDA::torusUnif(n = n_row,a = 1,c = 2),
                        maxdim = 2,thresh = 2,ripser = ripser)
    diagram_sphere = PyH(X = TDA::sphereUnif(n = n_row,d = 2,r = 1),
                         maxdim = 2,thresh = 2,ripser = ripser)
    # compute their wasserstein distances in all dimensions and benchmark
    start_time_TDApplied = Sys.time()
    diagram_distance(D1 = diagram_torus,D2 = diagram_sphere,dim = 0,
                     distance = "fisher",sigma = 1)
    diagram_distance(D1 = diagram_torus,D2 = diagram_sphere,dim = 1,
                     distance = "fisher",sigma = 1)
    diagram_distance(D1 = diagram_torus,D2 = diagram_sphere,dim = 2,
                     distance = "fisher",sigma = 1)
    end_time_TDApplied = Sys.time()
    time_diff_TDApplied = as.numeric(end_time_TDApplied - start_time_TDApplied,units = "secs")
    
    start_time_TDApplied_approx = Sys.time()
    diagram_distance(D1 = diagram_torus,D2 = diagram_sphere,dim = 0,
                     distance = "fisher",sigma = 1,rho = 0.001)
    diagram_distance(D1 = diagram_torus,D2 = diagram_sphere,dim = 1,
                     distance = "fisher",sigma = 1,rho = 0.001)
    diagram_distance(D1 = diagram_torus,D2 = diagram_sphere,dim = 2,
                     distance = "fisher",sigma = 1,rho = 0.001)
    end_time_TDApplied_approx = Sys.time()
    time_diff_TDApplied_approx = as.numeric(end_time_TDApplied_approx - start_time_TDApplied_approx,units = "secs")
    
    start_time_rgudhi = Sys.time()
    dis_fish$apply(diagram_torus[which(diagram_torus[,1] == 0),2:3],diagram_sphere[which(diagram_sphere[,1] == 0),2:3])
    dis_fish$apply(diagram_torus[which(diagram_torus[,1] == 1),2:3],diagram_sphere[which(diagram_sphere[,1] == 1),2:3])
    dis_fish$apply(diagram_torus[which(diagram_torus[,1] == 2),2:3],diagram_sphere[which(diagram_sphere[,1] == 2),2:3])
    end_time_rgudhi = Sys.time()
    time_diff_rgudhi = as.numeric(end_time_rgudhi - start_time_rgudhi,units = "secs")
    
    # start_time_rgudhi_approx = Sys.time()
    # dis_fish_approx$apply(diagram_torus[which(diagram_torus[,1] == 0),2:3],diagram_sphere[which(diagram_sphere[,1] == 0),2:3])
    # dis_fish_approx$apply(diagram_torus[which(diagram_torus[,1] == 1),2:3],diagram_sphere[which(diagram_sphere[,1] == 1),2:3])
    # dis_fish_approx$apply(diagram_torus[which(diagram_torus[,1] == 2),2:3],diagram_sphere[which(diagram_sphere[,1] == 2),2:3])
    # end_time_rgudhi_approx = Sys.time()
    # time_diff_rgudhi_approx = as.numeric(end_time_rgudhi_approx - start_time_rgudhi_approx,units = "secs")
    
    runtimes_rgudhi = rbind(runtimes_rgudhi,data.frame(n_row = n_row,
                                                       package = "TDApplied",
                                                       approx = F,
                                                       time_in_sec = time_diff_TDApplied))
    runtimes_rgudhi = rbind(runtimes_rgudhi,data.frame(n_row = n_row,
                                                       package = "TDApplied",
                                                       approx = T,
                                                       time_in_sec = time_diff_TDApplied_approx))
    runtimes_rgudhi = rbind(runtimes_rgudhi,data.frame(n_row = n_row,
                                                       package = "rgudhi",
                                                       approx = F,
                                                       time_in_sec = time_diff_rgudhi))
    # runtimes_rgudhi = rbind(runtimes_rgudhi,data.frame(n_row = n_row,
    #                                                    package = "rgudhi",
    #                                                    approx = T,
    #                                                    time_in_sec = time_diff_rgudhi_approx))
  }
  print(paste0("Done ",n_row," rows"))
}
# compute means and sd's at each value of rows for both packages
summary_table_rgudhi = data.frame(n_row = numeric(),mean = numeric(),sd = numeric(),
                                    package = character(),approx = logical())
for(n_row in seq(100,1000,100))
{
  for(p in c("TDApplied","rgudhi"))
  {
    approx_vec <- list(FALSE)
    if(p == "TDApplied")
    {
      approx_vec <- list(TRUE,FALSE)
    }
    for(a in approx_vec)
    {
      result = data.frame(n_row = n_row,
                          mean = mean(runtimes_rgudhi[which(runtimes_rgudhi$n_row == n_row
                                                            & runtimes_rgudhi$package == p
                                                            & runtimes_rgudhi$approx == a),
                                                      4]),
                          sd = sd(runtimes_rgudhi[which(runtimes_rgudhi$n_row == n_row 
                                                        & runtimes_rgudhi$package == p
                                                        & runtimes_rgudhi$approx == a),
                                                  4]),
                          package = p,approx = a)
      summary_table_rgudhi = rbind(summary_table_rgudhi,result) 
    }
  }
}
# plot table just with TDApplied approximation and rgudhi
summary_table_rgudhi <- summary_table_rgudhi[which(summary_table_rgudhi$package == "rgudhi" | summary_table_rgudhi$approx == T),1:4]
plot(summary_table_rgudhi$n_row[summary_table_rgudhi$package=="rgudhi"], 
     summary_table_rgudhi$mean[summary_table_rgudhi$package=="rgudhi"], 
     type="b",
     xlim=range(summary_table_rgudhi$n_row),
     ylim=range(0,summary_table_rgudhi$mean+1.96*summary_table_rgudhi$sd/sqrt(10)),
     xlab = "Points in shape",ylab = "Mean execution time (sec)",col = "blue")
lines(summary_table_rgudhi$n_row[summary_table_rgudhi$package=="TDApplied"],
      summary_table_rgudhi$mean[summary_table_rgudhi$package=="TDApplied"], 
      col=2, type="b")
legend(x = 200,y = 0.7,legend = c("TDApplied","rgudhi"),
       col = c("red","blue"),lty = c(1,1),cex = 0.8)
arrows(summary_table_rgudhi$n_row[summary_table_rgudhi$package == "TDApplied"], 
       summary_table_rgudhi$mean[summary_table_rgudhi$package == "TDApplied"]
       -1.96*summary_table_rgudhi$sd[summary_table_rgudhi$package == "TDApplied"]/sqrt(10),
       summary_table_rgudhi$n_row[summary_table_rgudhi$package == "TDApplied"], 
       summary_table_rgudhi$mean[summary_table_rgudhi$package == "TDApplied"]
       +1.96*summary_table_rgudhi$sd[summary_table_rgudhi$package == "TDApplied"]/sqrt(10), 
       length=0.05, angle=90, code=3,col = "red")
arrows(summary_table_rgudhi$n_row[summary_table_rgudhi$package == "rgudhi"], 
       summary_table_rgudhi$mean[summary_table_rgudhi$package == "rgudhi"]
       -1.96*summary_table_rgudhi$sd[summary_table_rgudhi$package == "rgudhi"]/sqrt(10), 
       summary_table_rgudhi$n_row[summary_table_rgudhi$package == "rgudhi"], 
       summary_table_rgudhi$mean[summary_table_rgudhi$package == "rgudhi"]
       +1.96*summary_table_rgudhi$sd[summary_table_rgudhi$package == "rgudhi"]/sqrt(10), 
       length=0.05, angle=90, code=3,col = "blue")

#### PyH vs. calculate_homology vs. compute_homology ####

if(requireNamespace("reticulate",quietly = T) == T)
{
  # generate persistence diagrams from circles with  100,200,...,1000 data points.
  runtimes_circle <- data.frame(n_row = numeric(),package = character(),time_in_sec = numeric())
  for(n_row in seq(100,1000,100)){
    
    for(iteration in 1:10)
    {
      # simulate a circle
      circ <- TDA::circleUnif(n = n_row)
      
      # compute their diagrams in all dimensions and benchmark
      start_time_TDApplied = Sys.time()
      phom_TDApplied <- PyH(circ,maxdim = 1,thresh = 1,ripser = ripser)
      end_time_TDApplied = Sys.time()
      time_diff_TDApplied = as.numeric(end_time_TDApplied - start_time_TDApplied,units = "secs")
      
      start_time_TDAstats = Sys.time()
      phom_TDAstats <- TDAstats::calculate_homology(circ,threshold = 1)
      end_time_TDAstats = Sys.time()
      time_diff_TDAstats = as.numeric(end_time_TDAstats - start_time_TDAstats,units = "secs")
      
      start_time_rgudhi = Sys.time()
      rc <- rgudhi::RipsComplex$new(data = circ, max_edge_length = 1)
      st <- rc$create_simplex_tree(1)
      st$compute_persistence()
      dim_0 <- as.data.frame(st$persistence_intervals_in_dimension(0))
      dim_0$dimension = 0
      dim_1 <- as.data.frame(st$persistence_intervals_in_dimension(1))
      if(nrow(dim_1) > 0)
      {
        dim_1$dimension = 1
      }else
      {
        dim_1$dimension = numeric() 
      }
      phom_rgudhi <- rbind(dim_0,dim_1)
      phom_rgudhi <- phom_rgudhi[,c("dimension","birth","death")]
      end_time_rgudhi = Sys.time()
      time_diff_rgudhi = as.numeric(end_time_rgudhi - start_time_rgudhi,units = "secs")
      
      runtimes_circle = rbind(runtimes_circle,data.frame(n_row = n_row,
                                                         package = "TDApplied",
                                                         time_in_sec = time_diff_TDApplied))
      runtimes_circle = rbind(runtimes_circle,data.frame(n_row = n_row,
                                                         package = "TDAstats",
                                                         time_in_sec = time_diff_TDAstats))
      runtimes_circle = rbind(runtimes_circle,data.frame(n_row = n_row,
                                                         package = "rgudhi",
                                                         time_in_sec = time_diff_rgudhi))
      
    }
    print(paste0("Done ",n_row," rows"))
    
  }
  # compute means and sd's at each value of rows for both packages
  summary_table_circle = data.frame(n_row = numeric(),mean = numeric(),sd = numeric(),
                                    package = character())
  for(n_row in seq(100,1000,100))
  {
    for(p in c("TDApplied","TDAstats","rgudhi"))
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
      
      # compute their diagrams in all dimensions and benchmark
      start_time_TDApplied = Sys.time()
      phom_TDApplied <- PyH(torus,maxdim = 2,thresh = 1,ripser = ripser)
      end_time_TDApplied = Sys.time()
      time_diff_TDApplied = as.numeric(end_time_TDApplied - start_time_TDApplied,units = "secs")
      
      start_time_TDAstats = Sys.time()
      phom_TDAstats <- TDAstats::calculate_homology(torus,threshold = 1,dim = 2)
      end_time_TDAstats = Sys.time()
      time_diff_TDAstats = as.numeric(end_time_TDAstats - start_time_TDAstats,units = "secs")
      
      start_time_rgudhi = Sys.time()
      rc <- rgudhi::RipsComplex$new(data = torus, max_edge_length = 1)
      st <- rc$create_simplex_tree(2)
      st$compute_persistence()
      dim_0 <- as.data.frame(st$persistence_intervals_in_dimension(0))
      dim_0$dimension = 0
      dim_1 <- as.data.frame(st$persistence_intervals_in_dimension(1))
      if(nrow(dim_1) > 0)
      {
        dim_1$dimension = 1
      }else
      {
        dim_1$dimension = numeric() 
      }
      dim_2 <- as.data.frame(st$persistence_intervals_in_dimension(2))
      if(nrow(dim_2) > 0)
      {
        dim_2$dimension = 2
      }else
      {
        dim_2$dimension = numeric() 
      }
      phom_rgudhi <- rbind(dim_0,dim_1,dim_2)
      phom_rgudhi <- phom_rgudhi[,c("dimension","birth","death")]
      end_time_rgudhi = Sys.time()
      time_diff_rgudhi = as.numeric(end_time_rgudhi - start_time_rgudhi,units = "secs")
      
      runtimes_torus = rbind(runtimes_torus,data.frame(n_row = n_row,
                                                       package = "TDApplied",
                                                       time_in_sec = time_diff_TDApplied))
      runtimes_torus = rbind(runtimes_torus,data.frame(n_row = n_row,
                                                       package = "TDAstats",
                                                       time_in_sec = time_diff_TDAstats))
      runtimes_torus = rbind(runtimes_torus,data.frame(n_row = n_row,
                                                       package = "rgudhi",
                                                       time_in_sec = time_diff_rgudhi))
      
    }
    print(paste0("Done ",n_row," rows"))
    
  }
  # compute means and sd's at each value of rows for both packages
  summary_table_torus = data.frame(n_row = numeric(),mean = numeric(),sd = numeric(),
                                   package = character())
  for(n_row in seq(100,1000,100))
  {
    for(p in c("TDApplied","TDAstats","rgudhi"))
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
      sphere <- TDA::sphereUnif(n = n_row,d = 2)
      
      # compute their diagrams in all dimensions and benchmark
      start_time_TDApplied = Sys.time()
      phom_TDApplied <- PyH(sphere,maxdim = 2,thresh = 1,ripser = ripser)
      end_time_TDApplied = Sys.time()
      time_diff_TDApplied = as.numeric(end_time_TDApplied - start_time_TDApplied,units = "secs")
      
      start_time_TDAstats = Sys.time()
      phom_TDAstats <- TDAstats::calculate_homology(sphere,threshold = 1,dim = 2)
      end_time_TDAstats = Sys.time()
      time_diff_TDAstats = as.numeric(end_time_TDAstats - start_time_TDAstats,units = "secs")
      
      start_time_rgudhi = Sys.time()
      rc <- rgudhi::RipsComplex$new(data = sphere, max_edge_length = 1)
      st <- rc$create_simplex_tree(2)
      st$compute_persistence()
      dim_0 <- as.data.frame(st$persistence_intervals_in_dimension(0))
      dim_0$dimension = 0
      dim_1 <- as.data.frame(st$persistence_intervals_in_dimension(1))
      if(nrow(dim_1) > 0)
      {
        dim_1$dimension = 1
      }else
      {
        dim_1$dimension = numeric() 
      }
      dim_2 <- as.data.frame(st$persistence_intervals_in_dimension(2))
      if(nrow(dim_2) > 0)
      {
        dim_2$dimension = 2
      }else
      {
        dim_2$dimension = numeric() 
      }
      phom_rgudhi <- rbind(dim_0,dim_1,dim_2)
      phom_rgudhi <- phom_rgudhi[,c("dimension","birth","death")]
      end_time_rgudhi = Sys.time()
      time_diff_rgudhi = as.numeric(end_time_rgudhi - start_time_rgudhi,units = "secs")
      
      runtimes_sphere = rbind(runtimes_sphere,data.frame(n_row = n_row,
                                                         package = "TDApplied",
                                                         time_in_sec = time_diff_TDApplied))
      runtimes_sphere = rbind(runtimes_sphere,data.frame(n_row = n_row,
                                                         package = "TDAstats",
                                                         time_in_sec = time_diff_TDAstats))
      runtimes_sphere = rbind(runtimes_sphere,data.frame(n_row = n_row,
                                                         package = "rgudhi",
                                                         time_in_sec = time_diff_rgudhi))
      
    }
    print(paste0("Done ",n_row," rows"))
    
  }
  # compute means and sd's at each value of rows for both packages
  summary_table_sphere = data.frame(n_row = numeric(),mean = numeric(),sd = numeric(),
                                    package = character())
  for(n_row in seq(100,1000,100))
  {
    for(p in c("TDApplied","TDAstats","rgudhi"))
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
  lines(summary_table_circle$n_row[summary_table_circle$package=="rgudhi"],
        summary_table_circle$mean[summary_table_circle$package=="rgudhi"], 
        col=4, type="b")
  legend(x = 200,y = 1.5,legend = c("TDApplied","TDAstats","rgudhi"),
         col = c("red","black","blue"),lty = c(1,1),cex = 0.8)
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
  arrows(summary_table_circle$n_row[summary_table_circle$package == "rgudhi"], 
         summary_table_circle$mean[summary_table_circle$package == "rgudhi"]
         -1.96*summary_table_circle$sd[summary_table_circle$package == "rgudhi"]/sqrt(10), 
         summary_table_circle$n_row[summary_table_circle$package == "rgudhi"], 
         summary_table_circle$mean[summary_table_circle$package == "rgudhi"]
         +1.96*summary_table_circle$sd[summary_table_circle$package == "rgudhi"]/sqrt(10), 
         length=0.05, angle=90, code=3,col = "blue")
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
  lines(summary_table_torus$n_row[summary_table_torus$package=="rgudhi"],
        summary_table_torus$mean[summary_table_torus$package=="rgudhi"], 
        col=4, type="b")
  legend(x = 200,y = 120,legend = c("TDApplied","TDAstats","rgudhi"),
         col = c("red","black","blue"),lty = c(1,1),cex = 0.8)
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
  arrows(summary_table_torus$n_row[summary_table_torus$package == "rgudhi"], 
         summary_table_torus$mean[summary_table_torus$package == "rgudhi"]
         -1.96*summary_table_torus$sd[summary_table_torus$package == "rgudhi"]/sqrt(10), 
         summary_table_torus$n_row[summary_table_torus$package == "rgudhi"], 
         summary_table_torus$mean[summary_table_torus$package == "rgudhi"]
         +1.96*summary_table_torus$sd[summary_table_torus$package == "rgudhi"]/sqrt(10), 
         length=0.05, angle=90, code=3,col = "blue")
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
  lines(summary_table_sphere$n_row[summary_table_sphere$package=="rgudhi"],
        summary_table_sphere$mean[summary_table_sphere$package=="rgudhi"], 
        col=4, type="b")
  legend(x = 200,y = 45,legend = c("TDApplied","TDAstats","rgudhi"),
         col = c("red","black","blue"),lty = c(1,1),cex = 0.8)
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
  arrows(summary_table_sphere$n_row[summary_table_sphere$package == "rgudhi"], 
         summary_table_sphere$mean[summary_table_sphere$package == "rgudhi"]
         -1.96*summary_table_sphere$sd[summary_table_sphere$package == "rgudhi"]/sqrt(10), 
         summary_table_sphere$n_row[summary_table_sphere$package == "rgudhi"], 
         summary_table_sphere$mean[summary_table_sphere$package == "rgudhi"]
         +1.96*summary_table_sphere$sd[summary_table_sphere$package == "rgudhi"]/sqrt(10), 
         length=0.05, angle=90, code=3,col = "blue")
  
}


#### comparing distance calcs ####
# requires additional package rgudhi with python!

# create diagrams
D1 = data.frame(dimension = c(0),birth = c(2),death = c(3))
D2 = data.frame(dimension = c(0),birth = c(2,0),death = c(3.3,0.5))
D3 = data.frame(dimension = c(0),birth = c(0),death = c(0.5))

# format
D1_TDA <- as.matrix(D1)
colnames(D1_TDA) <- NULL
D2_TDA <- as.matrix(D2)
colnames(D2_TDA) <- NULL
D3_TDA <- as.matrix(D3)
colnames(D3_TDA) <- NULL

# create rgudhi distance objects
dis_bot <- rgudhi::BottleneckDistance$new()
dis_wass <- rgudhi::WassersteinDistance$new()
dis_fish <- rgudhi::PersistenceFisherDistance$new() # default sigma is 1

# calculate tables
bottleneck_comparison <- data.frame(pair = c("D1 and D2","D1 and D3","D2 and D3"),ground_truth = c(0.3,0.5,0.65),TDApplied = c(diagram_distance(D1,D2,p = Inf),diagram_distance(D1,D3,p = Inf),diagram_distance(D2,D3,p = Inf)),TDA = c(TDA::bottleneck(D1_TDA,D2_TDA,dimension = 0),TDA::bottleneck(D1_TDA,D3_TDA,dimension = 0),TDA::bottleneck(D2_TDA,D3_TDA,dimension = 0)),TDAstats = c(TDAstats::phom.dist(D1_TDA,D2_TDA,limit.num = 0)[[1]],TDAstats::phom.dist(D1_TDA,D3_TDA,limit.num = 0)[[1]],TDAstats::phom.dist(D2_TDA,D3_TDA,limit.num = 0)[[1]]),rgudhi = c(dis_bot$apply(D1[,2:3],D2[,2:3]),dis_bot$apply(D1[,2:3],D3[,2:3]),dis_bot$apply(D2[,2:3],D3[,2:3])))

wasserstein_comparison <- data.frame(pair = c("D1 and D2","D1 and D3","D2 and D3"),ground_truth = c(0.3,0.5,0.65),TDApplied = c(diagram_distance(D1,D2),diagram_distance(D1,D3),diagram_distance(D2,D3)),TDA = c(TDA::wasserstein(D1_TDA,D2_TDA,dimension = 0),TDA::wasserstein(D1_TDA,D3_TDA,dimension = 0),TDA::wasserstein(D2_TDA,D3_TDA,dimension = 0)),TDAstats = c(TDAstats::phom.dist(D1_TDA,D2_TDA,limit.num = 0)[[1]],TDAstats::phom.dist(D1_TDA,D3_TDA,limit.num = 0)[[1]],TDAstats::phom.dist(D2_TDA,D3_TDA,limit.num = 0)[[1]]),rgudhi = c(dis_wass$apply(D1[,2:3],D2[,2:3]),dis_wass$apply(D1[,2:3],D3[,2:3]),dis_wass$apply(D2[,2:3],D3[,2:3])))

fisher_comparison <- data.frame(pair = c("D1 and D2","D1 and D3","D2 and D3"),ground_truth = c(0.02354624,0.08821907,0.1139891),TDApplied = c(diagram_distance(D1,D2,distance = "fisher",sigma = 1),diagram_distance(D1,D3,distance = "fisher",sigma = 1),diagram_distance(D2,D3,distance = "fisher",sigma = 1)),rgudhi = c(dis_fish$apply(D1[,2:3],D2[,2:3]),dis_fish$apply(D1[,2:3],D3[,2:3]),dis_fish$apply(D2[,2:3],D3[,2:3])))
