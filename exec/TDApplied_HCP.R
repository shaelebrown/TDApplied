
# This script automatically downloads HCP fMRI data from 100 subjects
# and analyzes the files with TDApplied as described in the package
# vignette. It is currently configured to run on Windows or Mac OS machines,
# but can run on other operating systems by changing the commands
# which create directories to your desired system commands (lines 105, 372-383).
# If run from top to bottom (once updating the connectomedb login and desired path
# on lines 41 and 958 respectively) this script will perform the same
# analyses as in the vignette for the same subjects, although
# there is the option to run the analysis for 100 randomly selected
# subjects from HCP (see the notes in the analyze_HCP function). 
# To analyze a different number of subjects other changes must be made to the code.
# Note that the desired path should be a full path
# without the ~ symbol (which can be used on Mac OS), otherwise lines
# 491 and 494 may throw an error.

# The plotting of the HCP resting state 1 loop requires python to be configured
# as well as the additional modules nibabel, nilearn, matplotlib and hcp_utils
# to be downloaded (with reticulate::py_install("module_name"), except for hcp_utils which must be downloaded
# with pip in the terminal: pip install hcp_utils).
# If this is not possible or desired, simply comment lines 491 and 494 to avoid generating an error.
# Also lines 491 and 494 seem to throw an error on
# Windows machines, but run fine on Mac OS (so comment them if need be).

# Once a 'directory_for_subjects' path is defined, which should have no files in it, each subject's persistence
# diagrams are saved in 'directory_for_subjects/subject_ID', which is created by this script.
# The results of applied analyses are saved in 'directory_for_subjects/analysis_results'. 
# The plots for the resting state 1 loop are saved in directory_for_subjects directly.

# initialize script for memory and reproducibility
set.seed(123)
memory.limit(850000) # need large memory limit to work with large files

# REQUIRED LIBS
library(neurohcp)
library(RNiftyReg)
library(TDApplied)

# sign on to aws using personal access key and secret key
# need connectomedb account @ https://db.humanconnectome.org
set_aws_api_key(access_key = "XXXXXXXXXXXXXXXXXXXX", secret_key = "YYYYYYYYYYYYYYYYYYYYYYYYYYYY")

# FUNCTIONS

# try to download one HCP fMRI file 'f' via AWS
try_to_download_one_file <- function(f,directory_for_subjects){
  
  # in case certain files don't exist for some subjects
  # wrap in tryCatch
  fl <- strsplit(f,split = "/")[[1]][[6]]
  s <- strsplit(f,split = "/")[[1]][[2]]
  m <- tryCatch(expr = {neurohcp::download_hcp_file(path_to_file = f,destfile = paste0(directory_for_subjects,"/",s,"/",fl))
    return()},error = function(e){
      
      unlink(paste0(directory_for_subjects,"/",s,"/",fl))
      
      return(e)})
  
  return(m)
  
}

# download a subject's stats file via AWS
download_RL_emotion_stats_file <- function(s,directory_for_subjects){
  
  f <- paste0("HCP_1200/",s,"/MNINonLinear/Results/tfMRI_EMOTION_RL/EVs/EMOTION_Stats.csv")
  m <- tryCatch(expr = {neurohcp::download_hcp_file(path_to_file = f,destfile = paste0(directory_for_subjects,"/",s,"/stats.csv"))
    return()},error = function(e){
      
      unlink(paste0(directory_for_subjects,"/",s,"/",fl))
      
      return(e)})
  
  return(m)
  
}

# scale a vector to colors, between blue (lowest), white (middle) and red (highest)
scale_between_r_w_and_b <- function(vec){
  
  min_val <- min(vec,na.rm = T)
  max_val <- max(vec[which(is.finite(vec))],na.rm = T)
  return(unlist(lapply(X = vec,FUN = function(X){
    
    if(is.infinite(X))
    {
      return(rgb(0,0,0,1))
    }
    X <- (X - min_val)/(max_val - min_val)
    if(X < 0.5)
    {
      return(rgb(red = 2*X, green = 2*X, blue = 1, alpha = 1))
    }
    return(rgb(red = 1, green = 2*(1 - X), blue = 2*(1 - X), alpha = 1))
    
  })))
  
}

# normalize a vector between 0 and 1
normalize <- function(v){
  
  return((v - min(v))/(max(v) - min(v)))
  
}

# temporary
directory_for_subjects = "/Users/jibaccount/Downloads/HCP_test"
devtools::load_all()

# carry out analysis
analyze_HCP <- function(directory_for_subjects){
  
  # download the specific file needed
  try_to_download_one_file("HCP_1200/103111/MNINonLinear/Results/tfMRI_EMOTION_RL/tfMRI_EMOTION_RL_Atlas_MSMAll.dtseries.nii",directory_for_subjects = directory_for_subjects)
  
  # set file names
  f = paste0(directory_for_subjects,"/103111/tfMRI_EMOTION_RL_Atlas_MSMAll.dtseries.nii")
  
  # read in data
  emotion = readNifti(f)
  emotion <- emotion[1,1,1,1,1:dim(emotion)[[5]],1:91282]
  
  # create representational dissimilarity matrix
  RDM <- sqrt(2*(matrix(data = 1,nrow = nrow(emotion),ncol = nrow(emotion)) - cor(t(emotion))))
  
  # compute diagram up to dimension 1 with representatives
  boot <- bootstrap_persistence_thresholds(X = RDM,FUN_diag = 'ripsDiag',maxdim = 1,thresh = 2,calculate_representatives = T,distance_mat = T,return_subsetted = T,return_diag = T)
  PH <- TDA::ripsDiag(X = RDM,maxdimension = 1,maxscale = max(RDM),dist = "arbitrary",library = "dionysus",location = T)
  
  # get most persistent loop
  diag <- diagram_to_df(PH)
  repr <- 296
  loop <- as.numeric(diag[repr,2:3])
  cycle <- unique(as.numeric(boot$representatives[[repr]]))
  secondary_cycle <- unique(as.numeric(PH$cycleLocation[[207]]))
  
  # plot diagram
  plot_diagram(diag,title = "Subject 103111 emotion task RL",thresholds = boot$thresholds,max_radius = 0.26)
  
  # compute rips graphs and plot
  RDM_cycle = RDM[cycle,cycle]
  colnames(RDM_cycle) <- cycle
  rownames(RDM_cycle) <- cycle
  cols <- rep("lightblue",ncol(RDM))
  cols[cycle] <- "red"
  g <- vr_graphs(X = RDM_cycle,distance_mat = T,eps = loop[[1]],return_clusters = F)
  plot_vr_graph(g,eps = loop[[1]],vertex_labels = T,title = "Rips graph of representative cycle") # clear loop
  g <- vr_graphs(X = RDM,distance_mat = T,eps = loop[[1]],return_clusters = F)
  l <- plot_vr_graph(g,eps = loop[[1]],component_of = cycle[[1]],vertex_labels = F,cols = cols,return_layout = T,title = "Rips graph of all data") # two clear loops! One big, one dense
  
  # gather stimulus onset and physiological data
  # and summarize by TR
  stimulus_timing <- data.frame(condition = c("fear","fear","fear","shape","shape","shape"),onset = c(32.066, 74.209, 116.352, 10.995, 53.138, 95.28))
  stimulus_timing <- stimulus_timing[order(stimulus_timing$onset),]
  stimulus_timing$end <- stimulus_timing$onset + 18
  neurohcp::download_hcp_file(path_to_file = "/HCP_1200/103111/MNINonLinear/Results/tfMRI_EMOTION_RL/tfMRI_EMOTION_RL_Physio_log.txt",destfile = paste0(directory_for_subjects,"/tfMRI_EMOTION_RL_Physio_log.txt"))
  physio_data <- read.table(paste0(directory_for_subjects,"/tfMRI_EMOTION_RL_Physio_log.txt"))
  colnames(physio_data) <- c("trigger_pulse","respiratory","pulse_oximeter")
  physio_data$TR <- floor(c(0:(nrow(physio_data) - 1))/400/0.72) + 1 # fMRI data was collected with 0.72s TR
  physio_data[nrow(physio_data),4L] <- 176
  physio_data$timing <- c(0:(nrow(physio_data) - 1))/400
  suppressWarnings({
    
    tr_data <- do.call(rbind,lapply(X = unique(physio_data$TR),FUN = function(X){
      
      df <- physio_data[which(physio_data$TR == X),]
      onset <- 0.72*(X - 1)
      return(data.frame(tr = X,onset = onset,trigger_pulse = mean(df$trigger_pulse),
                        respiratory = mean(df$respiratory),pulse_oximeter = mean(df$pulse_oximeter),
                        time_since_last_block = onset - max(stimulus_timing[which(stimulus_timing$onset <= onset),2L]),
                        time_since_last_fear_block = onset - max(stimulus_timing[which(stimulus_timing$onset <= onset & stimulus_timing$condition == "fear"),2L]),
                        time_since_last_shape_block = onset - max(stimulus_timing[which(stimulus_timing$onset <= onset & stimulus_timing$condition == "shape"),2L])))
      
    }))
    
  })
  
  # create colors based on TR data
  cols_trigger_pulse <- scale_between_r_w_and_b(tr_data$trigger_pulse) 
  cols_respiratory <- scale_between_r_w_and_b(tr_data$respiratory) 
  cols_pulse_oximeter <- scale_between_r_w_and_b(tr_data$pulse_oximeter)
  cols_time_since_last_block <- scale_between_r_w_and_b(tr_data$time_since_last_block)
  cols_time_since_last_fear_block <- scale_between_r_w_and_b(tr_data$time_since_last_fear_block)
  cols_time_since_last_shape_block <- scale_between_r_w_and_b(tr_data$time_since_last_shape_block)
  
  # plot rips graphs with each color scheme
  plot_vr_graph(g,eps = loop[[1]],cols = cols_trigger_pulse,component_of = cycle[[1]],vertex_labels = F,layout = l,title = "Trigger pulse") # not super interesting
  plot_vr_graph(g,eps = loop[[1]],cols = cols_respiratory,component_of = cycle[[1]],vertex_labels = F,layout = l,title = "Respiratory") # not super interesting
  plot_vr_graph(g,eps = loop[[1]],cols = cols_pulse_oximeter,component_of = cycle[[1]],vertex_labels = F,layout = l,title = "Pulse Oximeter") # not super interesting
  plot_vr_graph(g,eps = loop[[1]],cols = cols_time_since_last_block,component_of = cycle[[1]],vertex_labels = F,layout = l,title = "Time since last block") # interesting!!! same as time since last shape block
  plot_vr_graph(g,eps = loop[[1]],cols = cols_time_since_last_fear_block,component_of = cycle[[1]],vertex_labels = F,layout = l,title = "Time since last fear block") # not interesting
  plot_vr_graph(g,eps = loop[[1]],cols = cols_time_since_last_shape_block,component_of = cycle[[1]],vertex_labels = F,layout = l,title = "Time since last shape black") # interesting!! same as time since last block
  
  # main loop is the first stimulus block (shape)
  # is the secondary loop also time since shape block?
  secondary_loop <- g
  secondary_loop$vertices <- setdiff(g$vertices,as.character(cycle))
  secondary_loop$graphs[[1]][[1]] <- secondary_loop$graphs[[1]][[1]][which(secondary_loop$graphs[[1]][[1]][,1L] %in% cycle == F & secondary_loop$graphs[[1]][[1]][,2L] %in% cycle == F),]
  info <- plot_vr_graph(secondary_loop,eps = loop[[1]],component_of = 99,vertex_labels = F,cols = cols_time_since_last_block[as.numeric(secondary_loop$vertices)],return_layout = T,title = "Rips graph of secondary loop")
  
  # compute position around loop, theta, and
  # distance to loop center, r
  theta <- atan2(y = info[,2],x = info[,1])
  r <- sqrt(info[,1]^2 + info[,2]^2)
  
  # was this secondary loop the second most persistent loop?
  PH_secondary <- TDA::ripsDiag(X = RDM[as.numeric(rownames(info)),as.numeric(rownames(info))],dist = "arbitrary",library = "dionysus",maxdimension = 1,maxscale = max(RDM))
  diag_secondary <- diagram_to_df(PH_secondary)
  diag_secondary$pers <- diag_secondary$death - diag_secondary$birth
  diag_secondary <- diag_secondary[order(diag_secondary$pers,decreasing = T),]
  diag_secondary[min(which(diag_secondary$dimension == 1)),1:3]
  diag[207,]
  diag$pers <- diag$death - diag$birth
  diag <- diag[order(diag$pers,decreasing = T),]
  diag[which(diag$dimension == 1)[1:2],1:3] # it is the second biggest loop!!
  
  # modelling loop activity with theta and r
  stimulus_timing_TR <- stimulus_timing
  stimulus_timing_TR$onset <- 1 + stimulus_timing_TR$onset/0.72
  stimulus_timing_TR$end <- 1 + stimulus_timing_TR$end/0.72
  p_val_thresh <- 0.05/ncol(emotion)/2
  theta_coefs <- unlist(lapply(X = 1:ncol(emotion),FUN = function(X){
    
    mod <- lm(data = data.frame(x = theta,y = emotion[as.numeric(rownames(info)),X]),formula = y ~ cos(x) + sin(x))
    tab <- summary(mod)$coefficients
    coefs <- as.numeric(tab[2:3,1L])*ifelse(test = as.numeric(tab[2:3,4L]) < p_val_thresh,yes = 1,no = 0)
    names(coefs) <- c(paste0("cos_",X),paste0("sin_",X))
    return(coefs)
    
  }))
  r_coefs <- unlist(lapply(X = 1:ncol(emotion),FUN = function(X){
    
    mod <- lm(data = data.frame(x = r,y = emotion[as.numeric(rownames(info)),X]),formula = y ~ x)
    tab <- summary(mod)$coefficients
    return(as.numeric(tab[2,1])*ifelse(test = as.numeric(tab[2,4L]) < p_val_thresh,yes = 1,no = 0))
    
  }))
  length(which(theta_coefs != 0))/length(theta_coefs) # only about 0.8% of model coefficients were non-zero
  length(which(theta_coefs != 0)) # 1475
  length(which(r_coefs != 0))/length(r_coefs) # only about 0.06% of model coefficients were non-zero
  length(which(r_coefs != 0)) # 53
  non_zero_theta_nodes <- unlist(lapply(X = 1:ncol(emotion),FUN = function(X){
    
    return(length(which(theta_coefs[(2*X - 1):(2*X)] != 0)) > 0)
    
  }))
  non_zero_r_nodes <- ifelse(r_coefs != 0,yes = T,no = F)
  length(which(non_zero_theta_nodes == T))/length(non_zero_theta_nodes) # only about 1.6% of nodes had non-zero effects across the loop
  length(which(non_zero_theta_nodes == T)) # 1442
  
  length(setdiff(which(non_zero_r_nodes),which(non_zero_theta_nodes))) # 52! so pretty much all nodes orthogonal to the loop dimension
  length(setdiff(which(non_zero_theta_nodes),which(non_zero_r_nodes))) # 1441
  
  # let's plot the nodes in the loop
  hcp = reticulate::import("hcp_utils")
  atlas <- hcp$hcp_utils$mmp$map_all
  LUT <- unlist(hcp$hcp_utils$mmp$labels)
  nib = reticulate::import("nibabel")
  plotting = reticulate::import("nilearn.plotting")
  np = reticulate::import("numpy")
  
  # load surface mesh
  theta_nodes <- ifelse(non_zero_theta_nodes == T,yes = 1,no = 0)[1:64984]
  r_nodes <- ifelse(non_zero_r_nodes == T,yes = 1,no = 0)[1:64984]
  combined_nodes <- theta_nodes + r_nodes
  combined_nodes[which(combined_nodes > 1)] <- 1
  
  # plot loop nodes on left hemisphere
  plotting$plot_surf(hcp$mesh$inflated,reticulate::np_array(r_nodes),bg_map = hcp$mesh$sulc,hemi = "left",output_file = paste0(directory_for_subjects,"/r_nodes.png"),threshold = 0.001,colorbar = F)
  plotting$plot_surf(hcp$mesh$inflated,reticulate::np_array(theta_nodes),bg_map = hcp$mesh$sulc,hemi = "left",output_file = paste0(directory_for_subjects,"/theta_nodes.png"),threshold = 0.001,colorbar = F)
  plotting$plot_surf(hcp$mesh$inflated,reticulate::np_array(combined_nodes),bg_map = hcp$mesh$sulc,hemi = "left",output_file = paste0(directory_for_subjects,"/combined_nodes.png"),threshold = 0.001,colorbar = F)
  
  # now let's plot the loop activity around the loop for left hemi
  sampled_theta_inds <- unlist(lapply(X = c(0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4),FUN = function(X){
    
    theta2 <- ifelse(theta < 0,yes = theta + 2*pi,no = theta)
    dists <- abs(theta2 - X)
    return(which(dists == min(dists)))
    
  }))
  scaled_emotion <- apply(emotion[as.numeric(rownames(info)),1:64984],MARGIN = 2,FUN = function(X){
    
    return((X - min(X))/(max(X) - min(X)))
    
  })
  plotting$plot_surf(hcp$mesh$inflated,reticulate::np_array(theta_nodes*scaled_emotion[sampled_theta_inds[[1]],1:64984]),bg_map = hcp$mesh$sulc,hemi = "left",output_file = paste0(directory_for_subjects,"/loop_0.png"),threshold = 0.000002)
  plotting$plot_surf(hcp$mesh$inflated,reticulate::np_array(theta_nodes*scaled_emotion[sampled_theta_inds[[2]],1:64984]),bg_map = hcp$mesh$sulc,hemi = "left",output_file = paste0(directory_for_subjects,"/loop_pi_over_4.png"),threshold = 0.000002)
  plotting$plot_surf(hcp$mesh$inflated,reticulate::np_array(theta_nodes*scaled_emotion[sampled_theta_inds[[3]],1:64984]),bg_map = hcp$mesh$sulc,hemi = "left",output_file = paste0(directory_for_subjects,"/loop_pi_over_2.png"),threshold = 0.000002)
  plotting$plot_surf(hcp$mesh$inflated,reticulate::np_array(theta_nodes*scaled_emotion[sampled_theta_inds[[4]],1:64984]),bg_map = hcp$mesh$sulc,hemi = "left",output_file = paste0(directory_for_subjects,"/loop_3_pi_over_4.png"),threshold = 0.000002)
  plotting$plot_surf(hcp$mesh$inflated,reticulate::np_array(theta_nodes*scaled_emotion[sampled_theta_inds[[5]],1:64984]),bg_map = hcp$mesh$sulc,hemi = "left",output_file = paste0(directory_for_subjects,"/loop_pi.png"),threshold = 0.000002)
  plotting$plot_surf(hcp$mesh$inflated,reticulate::np_array(theta_nodes*scaled_emotion[sampled_theta_inds[[6]],1:64984]),bg_map = hcp$mesh$sulc,hemi = "left",output_file = paste0(directory_for_subjects,"/loop_5_pi_over_4.png"),threshold = 0.000002)
  plotting$plot_surf(hcp$mesh$inflated,reticulate::np_array(theta_nodes*scaled_emotion[sampled_theta_inds[[7]],1:64984]),bg_map = hcp$mesh$sulc,hemi = "left",output_file = paste0(directory_for_subjects,"/loop_3_pi_over_2.png"),threshold = 0.000002)
  plotting$plot_surf(hcp$mesh$inflated,reticulate::np_array(theta_nodes*scaled_emotion[sampled_theta_inds[[8]],1:64984]),bg_map = hcp$mesh$sulc,hemi = "left",output_file = paste0(directory_for_subjects,"/loop_7_pi_over_4.png"),threshold = 0.000002)

  plot(c(-1, 1), c(-1, 1), type = "n",xlab = "",ylab = "",xaxt = "n",yaxt = "n",bty = "n",xlim = c(-2.6,2.6),ylim = c(-2.6,2.6),main = "Activity with respect to theta")
  radius = 1
  center_x = 0
  center_y = 0
  theta = seq(0, 2 * pi, length = 200)
  lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y)
  
  angles <- c(0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4)
  angle_names <- c("0","pi_over_4","pi_over_2","3_pi_over_4","pi","5_pi_over_4","3_pi_over_2","7_pi_over_4")
  # this plot will only generate if the png library is installed
  if(requireNamespace("png"))
  {
    graphics::rasterImage(png::readPNG(paste0(directory_for_subjects,"/loop_",angle_names[[1]],".png")),xleft = 1.1,xright = 2.6,ybottom = -0.75,ytop = 0.75)
    graphics::rasterImage(png::readPNG(paste0(directory_for_subjects,"/loop_",angle_names[[2]],".png")),xleft = 0.8,xright = 2.3,ybottom = 0.6,ytop = 2.1)
    graphics::rasterImage(png::readPNG(paste0(directory_for_subjects,"/loop_",angle_names[[3]],".png")),xleft = -0.75,xright = 0.75,ybottom = 1.1,ytop = 2.6)
    graphics::rasterImage(png::readPNG(paste0(directory_for_subjects,"/loop_",angle_names[[4]],".png")),xleft = -2.3,xright = -0.8,ybottom = 0.6,ytop = 2.1)
    graphics::rasterImage(png::readPNG(paste0(directory_for_subjects,"/loop_",angle_names[[5]],".png")),xleft = -2.6,xright = -1.1,ybottom = -0.75,ytop = 0.75)
    graphics::rasterImage(png::readPNG(paste0(directory_for_subjects,"/loop_",angle_names[[6]],".png")),xleft = -2.3,xright = -0.8,ybottom = -2.1,ytop = -0.6)
    graphics::rasterImage(png::readPNG(paste0(directory_for_subjects,"/loop_",angle_names[[7]],".png")),xleft = -0.75,xright = 0.75,ybottom = -2.6,ytop = -1.1)
    graphics::rasterImage(png::readPNG(paste0(directory_for_subjects,"/loop_",angle_names[[8]],".png")),xleft = 0.8,xright = 2.3,ybottom = -2.1,ytop = -0.6)
  }
  
  # now let's plot activity over different values of r
  min_ind <- which(r == min(r))
  max_ind <- which(r == max(r))
  med_ind <- which(abs(r - mean(r)) == min(abs(r - mean(r))))[[1]]
  plotting$plot_surf(hcp$mesh$inflated,reticulate::np_array(r_nodes*scaled_emotion[min_ind,1:64984]),bg_map = hcp$mesh$sulc,hemi = "left",output_file = paste0(directory_for_subjects,"/r_min.png"),threshold = 0.000002)
  plotting$plot_surf(hcp$mesh$inflated,reticulate::np_array(r_nodes*scaled_emotion[med_ind,1:64984]),bg_map = hcp$mesh$sulc,hemi = "left",output_file = paste0(directory_for_subjects,"/r_mean.png"),threshold = 0.000002)
  plotting$plot_surf(hcp$mesh$inflated,reticulate::np_array(r_nodes*scaled_emotion[max_ind,1:64984]),bg_map = hcp$mesh$sulc,hemi = "left",output_file = paste0(directory_for_subjects,"/r_max.png"),threshold = 0.000002)
  
  plot(c(-1, 1), c(-1, 1), type = "n",xlab = "",ylab = "",xaxt = "n",yaxt = "n",bty = "n",main = "Activity with respect to r",xlim = c(-4.5,4.5),ylim = c(-1.5,1.5))
  axis(1, at=c(-3,0,3), las=2,labels = c("Min","Mean","Max"))
  # this plot will only generate if the png library is installed
  if(requireNamespace("png"))
  {
    graphics::rasterImage(png::readPNG(paste0(directory_for_subjects,"/r_min.png")),xleft = -4.5,xright = -1.5,ybottom = -1.5,ytop = 1.5)
    graphics::rasterImage(png::readPNG(paste0(directory_for_subjects,"/r_mean.png")),xleft = -1.5,xright = 1.5,ybottom = -1.5,ytop = 1.5)
    graphics::rasterImage(png::readPNG(paste0(directory_for_subjects,"/r_max.png")),xleft = 1.5,xright = 4.5,ybottom = -1.5,ytop = 1.5)
  }
  
  # comparing values of r across shapes and faces
  trs <- as.numeric(rownames(info))
  tr_onset <- stimulus_timing$onset/0.72 + 1
  tr_end <- stimulus_timing$end/0.72 + 1
  shape_inds <- trs[c(which(trs >= tr_onset[[1]] & trs <= tr_end[[1]]),which(trs >= tr_onset[[3]] & trs <= tr_end[[3]]),which(trs >= tr_onset[[5]] & trs <= tr_end[[5]]))]
  face_inds <- trs[c(which(trs >= tr_onset[[2]] & trs <= tr_end[[2]]),which(trs >= tr_onset[[4]] & trs <= tr_end[[4]]),which(trs >= tr_onset[[6]] & trs <= tr_end[[6]]))]
  shape_inds <- intersect(shape_inds,as.numeric(rownames(info)))
  face_inds <- intersect(face_inds,as.numeric(rownames(info)))
  shape_inds_2 <- c()
  face_inds_2 <- c()
  for(i in shape_inds)
  {
    shape_inds_2 <- c(shape_inds_2,which(trs == i))
  }
  for(i in face_inds)
  {
    face_inds_2 <- c(face_inds_2,which(trs == i))
  }
  t.test(r[shape_inds_2],r[face_inds_2],alternative = "two.sided")
  boxplot(r[shape_inds_2],r[face_inds_2],main = "r across conditions, p < 0.01",xlab = "Condition",ylab = "r",names = c("Shape","Face"))
  
  # now do embedding with 100 emotion RL diagrams
  source("./exec/parallel_with_approximation.R")
  ripser <- import_ripser()
  subjects = c(103111,103212,105620,106521,108222,110007,110613,113316,117122,118023,118528,118932,119126,122317,123420,123521,123723,129028,129634,130417,130720,131419,133827,135528,136631,136833,138332,138837,140319,140824,143224,143325,147030,147636,151324,153631,153934,154229,154835,156536,156637,158843,160729,162228,162329,173839,173940,176037,177241,178748,179245,180735,185947,186040,186141,187850,192237,198350,200008,202113,205826,206222,213522,219231,237334,239944,255639,299760,305830,329440,334635,341834,353740,395958,406432,424939,433839,456346,479762,510225,545345,555954,561444,562446,567052,571144,579665,579867,580751,586460,590047,599065,599671,657659,663755,687163,786569,788674,800941,818455)
  diags <- lapply(X = subjects,FUN = function(X){
    
    f <- paste0("HCP_1200/",X,"/MNINonLinear/Results/tfMRI_EMOTION_RL/tfMRI_EMOTION_RL_Atlas_MSMAll.dtseries.nii")
    try_to_download_one_file(f = f,directory_for_subjects = directory_for_subjects)
    f = paste0(directory_for_subjects,"/",X,"/tfMRI_EMOTION_RL_Atlas_MSMAll.dtseries.nii")
    dt = readNifti(f)
    dt <- dt[1,1,1,1,1:dim(dt)[[5]],1:91282]
    RDM <- sqrt(2*(matrix(data = 1,nrow = nrow(dt),ncol = nrow(dt)) - cor(t(dt))))
    return(PyH(X = RDM,thresh = max(RDM),distance_mat = T,ripser = ripser))
    
  })
  K <- parallel_approx_gram_matrix(diagrams = diags,dim = 1,sigma = 0.05,rho = 1e-4) # occasionally needs to be rerun
  emb <- diagram_kpca(diagrams = diags,K = K,dim = 1,t = 1,sigma = 0.05,rho = 1e-4,features = 2,th = 1e-6) # only the first dimension was used and plotted in the HCP analysis vignette
  
  # download subject emotion stats for plotting
  for(s in subjects)
  {
    download_RL_emotion_stats_file(s = s,directory_for_subjects = directory_for_subjects)
  }
  stats <- do.call(rbind,lapply(X = subjects,FUN = function(X){
    
    df <- read.csv(paste0(directory_for_subjects,"/",X,"/stats.csv"))
    df$X <- NULL
    return(data.frame(subject = X,face_acc = df[1,1L],face_rt = df[2,1L],shape_acc = df[3,1L],shape_rt = df[4,1L]))
    
  }))
  plot(emb$pca@rotated[,1],emb$pca@rotated[,2],xlab = "Embedding dim 1",ylab = "Embedding dim 2",main = "100 fMRI Persistence Diagrams",col = rgb(red = 1,green = 0,blue = 0,alpha = normalize(stats$shape_rt)))
  plot(emb$pca@rotated[,1],stats$shape_rt,xlab = "Embedding dim 1",ylab = "Mean Shape Block Reaction Time (ms)",main = "Topology-Behavior Relationship")
  l <- lm(formula = stats$shape_rt ~ emb$pca@rotated[,1])
  coefficients <- as.numeric(coef(l))
  graphics::lines(x = c(min(emb$pca@rotated[,1]),max(emb$pca@rotated[,1])),y = coefficients[[1]] + coefficients[[2]]*c(min(emb$pca@rotated[,1]),max(emb$pca@rotated[,1])))
  cor_val <- round(cor(emb$pca@rotated[,1],stats$shape_rt),digits = 2)
  graphics::text(x = 0,y = 1400,paste0("Correlation = ",cor_val))
  
}

# EXECUTION
analyze_HCP(directory_for_subjects = "desired/path")
