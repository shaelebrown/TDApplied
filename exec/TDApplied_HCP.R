
# This script automatically downloads HCP fMRI data from 100 subjects
# and analyzes the files with TDApplied as described in the package
# vignette. It is currently configured to run on Windows machines,
# but can run on other operating systems by changing the commands
# which create directories to your desired system commands (lines 107, 410-421).
# If run from top to bottom (once updating the connectomedb login and desired path
# on lines 41 and 909 respectively) this script will perform the same
# analyses as in the vignette for the same subjects, although
# there is the option to run the analysis for 100 randomly selected
# subjects from HCP (see the notes in the analyze_HCP function). 
# To analyze a different number of subjects other changes must be made to the code.
# Note that the desired path should be a full path
# without the ~ symbol (which can be used on Mac OS), otherwise lines
# 528 and 531 may throw an error.

# The plotting of HCP resting state 1 loops requires python to be configured
# as well as the additional modules nibabel, nilearn, matplotlib and hcp_utils
# to be downloaded (with reticulate::py_install("module_name"), except for hcp_utils which must be downloaded
# with pip in the terminal: pip install hcp_utils).
# If this is not possible or desired, simply comment lines 528 and 531 to avoid generating an error.
# Also lines 528 and 531 seem to throw an error on
# Windows machines, but run fine on Mac OS (so comment them if need be).

# Once a 'directory_for_subjects' path is defined, which should have no files in it, each subject's persistence
# diagrams are saved in 'directory_for_subjects/subject_ID', which is created by this script.
# The results of applied analyses are saved in 'directory_for_subjects/analysis_results'. 
# The plots for resting state 1 loops are saved in directory_for_subjects directly.

# initialize script for memory and reproducibility
set.seed(123)
memory.limit(850000) # need large memory limit to work with large files

# REQUIRED LIBS
library(neurohcp)
library(RNiftyReg)
library(TDApplied)

# sign on to aws using personal access key and secret key
# need connectomedb account @ https://db.humanconnectome.org
set_aws_api_key(access_key = "XXXXXXXXXXXXXXXXXXXX", secret_key = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

# FUNCTIONS

# compute persistence diagram of HCP fMRI file 'f'
PH_from_file <- function(f){
  
  # read in fMRI data, each column is a time point and each row is a cortical surface node
  fMRI <- readNifti(f)
  fMRI <- t(fMRI[1,1,1,1,1:dim(fMRI)[[5]],1:91282])
  
  # create representational dissimilarity matrix, which is 1 - correlation between
  # each pair of time points
  RSM <- cor(fMRI)
  RDM <- matrix(data = 1,nrow = nrow(RSM),ncol = ncol(RSM)) - RSM
  
  # compute diagram up to dimension 1 via fast bootstrapping
  # with 30 bootstrap iterations
  diag <- bootstrap_persistence_thresholds(X = RDM,FUN = "calculate_homology",maxdim = 1,thresh = max(RDM),num_samples = 30,return_subsetted = T,return_diag = T,distance_mat = T)$subsetted_diag
  
  # subset for dimenion 1
  diag <- diag[which(diag$dimension == 1),]
  
  # save results
  prefix <- strsplit(f,split = "_Atlas_MSMAll.dtseries.nii")[[1]][[1]]
  write.csv(diag,paste0(prefix,"_diagram.csv"))
  
  # remove raw file for clean up
  file.remove(f)
  
}

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

# try to download one HCP fMRI file 'f' via AWS and compute its persistence diagram
analyze_one_file <- function(f,directory_for_subjects){
  
  try_to_download_one_file(f,directory_for_subjects)
  s <- strsplit(f,"/")[[1]][[2]]
  prefix <- strsplit(f,paste0("HCP_1200/",s,"/MNINonLinear/Results/"))[[1]][[2]]
  prefix <- strsplit(prefix,"/")[[1]][[2]]
  PH_from_file(f = paste0(directory_for_subjects,"/",s,"/",prefix))
  
}

# try to download all HCP fMRI files via AWS for one subject 's'
# and compute their persistence diagrams
PH_subject_fMRI_data <- function(s,directory_for_subjects){

  # make directory for the subject
  system(command = paste0("powershell -command mkdir ",directory_for_subjects,"/",s))
  
  # download the subject's data and compute persistence diagrams
  
  # resting-state fMRI data
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/rfMRI_REST2_LR/rfMRI_REST2_LR_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/rfMRI_REST2_RL/rfMRI_REST2_RL_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)

  # task-based fMRI data
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/tfMRI_EMOTION_LR/tfMRI_EMOTION_LR_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/tfMRI_EMOTION_RL/tfMRI_EMOTION_RL_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/tfMRI_GAMBLING_LR/tfMRI_GAMBLING_LR_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/tfMRI_GAMBLING_RL/tfMRI_GAMBLING_RL_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/tfMRI_LANGUAGE_LR/tfMRI_LANGUAGE_LR_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/tfMRI_LANGUAGE_RL/tfMRI_LANGUAGE_RL_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/tfMRI_MOTOR_LR/tfMRI_MOTOR_LR_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/tfMRI_MOTOR_RL/tfMRI_MOTOR_RL_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/tfMRI_RELATIONAL_LR/tfMRI_RELATIONAL_LR_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/tfMRI_RELATIONAL_RL/tfMRI_RELATIONAL_RL_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/tfMRI_SOCIAL_LR/tfMRI_SOCIAL_LR_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/tfMRI_SOCIAL_RL/tfMRI_SOCIAL_RL_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/tfMRI_WM_LR/tfMRI_WM_LR_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)
  analyze_one_file(paste0("HCP_1200/",s,"/MNINonLinear/Results/tfMRI_WM_RL/tfMRI_WM_RL_Atlas_MSMAll.dtseries.nii"),directory_for_subjects)

}

# helper function for determining if a subject has all necessary fMRI files
# for analysis
subject_has_all_files <- function(subj){
  
  dirs <- unname(unlist(parse_list_files(hcp_list_dirs(paste0("HCP_1200/",subj,"/MNINonLinear/Results/")))$prefixes))
  if(is.null(dirs))
  {
    return(F)
  }
  dirs <- unlist(strsplit(dirs,split = paste0("HCP_1200/",subj,"/MNINonLinear/Results/")))
  
  if(length(dirs) == 50)
  {
    return(T)
  }
  return(F)
  
}

# calculate persistence diagrams for all fMRI files for 100 subjects
calculate_diags <- function(directory_for_subjects,subjects){
  
  # including timing to track progress over subjects
  S <- Sys.time()

  for(subj in subjects)
  {
    print(paste0("Time before starting subject ",subj,": ",Sys.time()))
    PH_subject_fMRI_data(subj,directory_for_subjects)
    print(paste0("Time after finishing subject ",subj,": ",Sys.time()))
  }
  
  print(paste0("Time taken for 100 subjects was ",Sys.time()-S))
  
}

# visualize a list of persistence diagrams via a heat kernel of parameter 'sigma'
visualize_group_of_diagrams <- function(diagrams,sigma,plot_title){
  
  diag <- do.call(rbind,diagrams)
  if(nrow(diag) == 0 | length(which(complete.cases(diag) == F)) == nrow(diag))
  {
    plot(1, type="n", xlab="Birth", ylab="Death", xlim=c(0, 1), ylim=c(0, 1),main = plot_title)
    abline(a = 0,b = 1)
    return()
  }
  min_x <- min(diag$birth) - 3*sigma
  max_x <- max(diag$birth) + 3*sigma
  min_y <- min(diag$death) - 3*sigma
  max_y <- max(diag$death) + 3*sigma
  x <- seq(min_x,max_x,(max_x-min_x)/100)
  y <- seq(min_y,max_y,(max_y-min_y)/100)
  z <- outer(x,y,FUN = function(x,y){
    
    sum <- 0
    for(i in 1:nrow(diag))
    {
      sum <- sum + (exp(-((x-diag[i,2])^2+(y-diag[i,3])^2)/(2*sigma^2)))/sqrt(2*pi*sigma^2)
    }
    return(sum)
    
  })
  z <- z/sum(z)
  image(x = x,y = y,z,main = plot_title,xlab = "Birth",xlim = c(min_x,max_x),ylim = c(min_y,max_y),ylab = "Death")
  abline(a = 0,b = 1)
  
}

# plot time point, smoothed spatially in neighborhoods of
# size sigma and thresholded by absolute value th
# nib, plotting, np and hcp are required python modules
plot_smooth_timepoint <- function(f,tp,output_file,time_point,sigma,th,nib,plotting,np,hcp){
  
  # load surface mesh
  coords = hcp$mesh$pial[[1]]
  
  # load image and normalize
  img = nib$load(f)
  X = img$get_fdata()
  Xn = hcp$normalize(X)
  
  # subset for specific time point
  spatial_pattern = Xn[tp,]
  
  # smooth a time point in parallel
  cl = parallel::makeCluster(parallel::detectCores() - 1)
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl,library("rdist"))
  parallel::clusterExport(cl,c("spatial_pattern","coords","sigma"),envir = environment())
  smoothed = foreach::`%dopar%`(foreach::foreach(n = 1:64984,.combine = c),{
    
    distances = rdist::cdist(coords,t(as.matrix(coords[n,])))
    return(mean(spatial_pattern[which(distances > 0 & distances <= sigma)]))
    
  })
  
  # clean up
  parallel::stopCluster(cl)
  
  # plot 
  plotting$plot_surf(hcp$mesh$inflated,reticulate::np_array(smoothed),threshold = th,bg_map = hcp$mesh$sulc,hemi = "right",output_file = paste0(output_file,"_hemi_right.png"))
  plotting$plot_surf(hcp$mesh$inflated,reticulate::np_array(smoothed),threshold = th,bg_map = hcp$mesh$sulc,hemi = "left",output_file = paste0(output_file,"_hemi_left.png"))
  
}

# analyze one loop representative using PyH
# mainly internal function
loop_scatterplot <- function(f,directory_for_subjects,ripser,loop_no,time_points,plot_title,output_file){
  
  # analyze file
  fMRI <- readNifti(f)
  fMRI <- t(fMRI[1,1,1,1,1:dim(fMRI)[[5]],1:91282])
  RSM <- cor(fMRI)
  RDM <- matrix(data = 1,nrow = nrow(RSM),ncol = ncol(RSM)) - RSM
  diag <- PyH(X = RDM,maxdim = 1,distance_mat = T,ripser = ripser,thresh = max(RDM),calculate_representatives = T)
  repr <- diag$representatives[[2]][[loop_no]]
  
  # get birth value of the representative and all indices
  # in the RDM which have distance at most the birth value
  birth_val <- min(RDM[repr])
  indices <- which(RDM <= birth_val,arr.ind = T)
  indices <- indices[which(indices[,1] < indices[,2]),]
  
  # find all time points which are within distance birth value
  # of each other and include the birth edge time points
  tp <- c(repr[1,1],repr[1,2])
  visited <- c()
  edges <- data.frame(from = numeric(),to = numeric())
  while(length(tp) > 0)
  {
    pt <- tp[[length(tp)]]
    visited <- c(visited,pt)
    adj <- setdiff(unique(c(indices[which(indices[,1] == pt),2],indices[which(indices[,2] == pt),1])),visited)
    tp <- c(tp,adj)
    tp <- tp[which(tp != pt)]
    if(length(adj) > 0)
    {
      edges <- rbind(edges,data.frame(from = rep(pt,length(adj)),to = adj)) 
    }
  }
  
  grDevices::png(output_file)
  
  # compute embedding and plot
  if(plot_title == "Loop 1")
  {
    emb <- cmdscale(RDM[visited,visited],k = 2)
    plot(x = emb[,1],y = emb[,2],xlab = "Embedding coordinate 1",ylab = "Embedding coordinate 2",pch = 19,main = plot_title)
    i <- 1
    j <- 2
  }else
  {
    # need 4 dimensions to clearly see the second loop
    emb <- cmdscale(RDM[visited,visited],k = 4)
    plot(x = emb[,3],y = emb[,4],xlab = "Embedding coordinate 3",ylab = "Embedding coordinate 4",pch = 19,main = plot_title)
    i <- 3
    j <- 4
  }
  # add red points for specific time points along the loop
  points(x = emb[time_points,i],y = emb[time_points,j],col = "red",pch = 19)
  
  dev.off()
  
}

# analyze two specific loops
analyze_two_loops <- function(directory_for_subjects,output_file_loop_1,output_file_loop_2,th,sigma){
  
  # download the two specific files needed
  try_to_download_one_file("HCP_1200/162228/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_Atlas_MSMAll.dtseries.nii",directory_for_subjects = directory_for_subjects)
  try_to_download_one_file("HCP_1200/341834/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_MSMAll.dtseries.nii",directory_for_subjects = directory_for_subjects)
  
  # set file names
  f1 = paste0(directory_for_subjects,"/162228/rfMRI_REST1_RL_Atlas_MSMAll.dtseries.nii")
  f2 = paste0(directory_for_subjects,"/341834/rfMRI_REST1_LR_Atlas_MSMAll.dtseries.nii")
  
  # plot scatter plots for two loops
  ripser = import_ripser()
  loop_scatterplot(f = f1,directory_for_subjects = directory_for_subjects,ripser = ripser,loop_no = 1233,time_points = c(62,46,10,101),output_file = paste0(output_file_loop_1,"_scatter.png"),plot_title = "Loop 1")
  loop_scatterplot(f = f2,directory_for_subjects = directory_for_subjects,ripser = ripser,loop_no = 1014,time_points = c(50,157,296,20),output_file = paste0(output_file_loop_2,"_scatter.png"),plot_title = "Loop 2")
  
  # import python modules
  nib = reticulate::import("nibabel")
  plotting = reticulate::import("nilearn.plotting")
  np = reticulate::import("numpy")
  hcp = reticulate::import("hcp_utils")
  
  # first loop
  # time points were 62 for the top, 46 for the right, 10 for the bottom and 101 for the left
  plot_smooth_timepoint(f = f1,tp = 62,output_file = paste0(output_file_loop_1,"_top"),th = th,sigma = sigma,nib = nib,plotting = plotting,np = np,hcp = hcp)
  plot_smooth_timepoint(f = f1,tp = 46,output_file = paste0(output_file_loop_1,"_right"),th = th,sigma = sigma,nib = nib,plotting = plotting,np = np,hcp = hcp)
  plot_smooth_timepoint(f = f1,tp = 10,output_file = paste0(output_file_loop_1,"_bottom"),th = th,sigma = sigma,nib = nib,plotting = plotting,np = np,hcp = hcp)
  plot_smooth_timepoint(f = f1,tp = 101,output_file = paste0(output_file_loop_1,"_left"),th = th,sigma = sigma,nib = nib,plotting = plotting,np = np,hcp = hcp)
  
  # second loop
  # time points were 50 for the top, 157 for the right, 296 for the bottom and 20 for the left
  plot_smooth_timepoint(f = f2,tp = 50,output_file = paste0(output_file_loop_2,"_top"),th = th,sigma = sigma,nib = nib,plotting = plotting,np = np,hcp = hcp)
  plot_smooth_timepoint(f = f2,tp = 157,output_file = paste0(output_file_loop_2,"_right"),th = th,sigma = sigma,nib = nib,plotting = plotting,np = np,hcp = hcp)
  plot_smooth_timepoint(f = f2,tp = 296,output_file = paste0(output_file_loop_2,"_bottom"),th = th,sigma = sigma,nib = nib,plotting = plotting,np = np,hcp = hcp)
  plot_smooth_timepoint(f = f2,tp = 20,output_file = paste0(output_file_loop_2,"_left"),th = th,sigma = sigma,nib = nib,plotting = plotting,np = np,hcp = hcp)
  
  # clean up files
  unlink(f1)
  unlink(f2)
  
}

# plot parcellation
view_parcellation <- function(output_file){
  
  # based on the hcp_utils function view_parcellation
  # initialize variables for yeo7 parcellation
  # can be adapted for other parcellations easily, though
  hcp <- reticulate::import("hcp_utils")
  np <- reticulate::import("numpy")
  matplotlib <- reticulate::import("matplotlib")
  plotting <- reticulate::import("nilearn.plotting")
  meshLR <- hcp$mesh$inflated
  parcellation <- hcp$yeo7
  
  cortex_map <- hcp$cortex_data(parcellation$map_all)
  ids <- np$unique(cortex_map)
  normalized_cortex_map <- np$zeros_like(cortex_map)
  rgba <- np$zeros(reticulate::tuple(length(ids),4L))
  for(i in 0:(length(ids)-1))
  {
    ind <- which(cortex_map == ids[[i + 1]])
    normalized_cortex_map[ind] = i
    rgba[i + 1,] <- get(as.character(i),parcellation$rgba)
  }
  
  cmap = matplotlib$colors$ListedColormap(rgba)
  
  # view parcellation on surface
  plotting$plot_surf(meshLR,normalized_cortex_map,threshold = 0,bg_map = hcp$mesh$sulc,hemi = "left",output_file = paste0(output_file,"_hemi_left.png"),cmap=cmap)
  plotting$plot_surf(meshLR,normalized_cortex_map,threshold = 0,bg_map = hcp$mesh$sulc,hemi = "right",output_file = paste0(output_file,"_hemi_right.png"),cmap=cmap)
  
  # view parcellation labels
  hcp$parcellation_labels(parcellation)
  matplotlib$pyplot$savefig(paste0(output_file,"_labels.png"))
  
}

# calculate persistence diagrams for all fMRI files for 100 subjects
# and analyze them
# if desired subjects can be set to NULL to randomly select 100 new subjects
analyze_HCP <- function(directory_for_subjects){
  
  subjects = c(103111,103212,105620,106521,108222,110007,110613,113316,117122,118023,118528,118932,119126,122317,123420,123521,123723,129028,129634,130417,130720,131419,133827,135528,136631,136833,138332,138837,140319,140824,143224,143325,147030,147636,151324,153631,153934,154229,154835,156536,156637,158843,160729,162228,162329,173839,173940,176037,177241,178748,179245,180735,185947,186040,186141,187850,192237,198350,200008,202113,205826,206222,213522,219231,237334,239944,255639,299760,305830,329440,334635,341834,353740,395958,406432,424939,433839,456346,479762,510225,545345,555954,561444,562446,567052,571144,579665,579867,580751,586460,590047,599065,599671,657659,663755,687163,786569,788674,800941,818455)
  
  # generate 100 new subjects for analysis if desired
  # (i.e. not the same 100 as were analyzed in the vignette)
  # note that subjects 162228 and 341834 must be analyzed in order
  # to not cause an error on lines X and Y.
  # if(is.null(subjects))
  # {
  #   possible_subjects <- unlist(strsplit(unname(unlist(parse_list_files(hcp_list_dirs("HCP_1200/"))$prefixes)),split = "HCP_1200/"))
  #   possible_subjects <- possible_subjects[seq(2,length(possible_subjects),2)]
  #   possible_subjects <- unlist(strsplit(possible_subjects,split = "/"))
  #   
  #   subjects_with_all_files <- unlist(lapply(X = possible_subjects,FUN = function(X){
  #     
  #     return(subject_has_all_files(X))
  #     
  #   }))
  #   
  #   possible_subjects <- possible_subjects[subjects_with_all_files]
  #   
  #   subjects <- sample(possible_subjects,size = 100,replace = F)
  # }
  
  # start by reading in data, then analyze:
  calculate_diags(directory_for_subjects = directory_for_subjects,subjects = subjects)
  
  # create directories to store results
  system(paste0("powershell -command mkdir ",directory_for_subjects,"/analysis_results"))
  system(paste0("powershell -command mkdir ",directory_for_subjects,"/analysis_results/permutation_tests"))
  system(paste0("powershell -command mkdir ",directory_for_subjects,"/analysis_results/permutation_tests/tasks"))
  system(paste0("powershell -command mkdir ",directory_for_subjects,"/analysis_results/permutation_tests/subjects"))
  system(paste0("powershell -command mkdir ",directory_for_subjects,"/analysis_results/independence_tests"))
  system(paste0("powershell -command mkdir ",directory_for_subjects,"/analysis_results/independence_tests/tasks"))
  system(paste0("powershell -command mkdir ",directory_for_subjects,"/analysis_results/independence_tests/subjects"))
  system(paste0("powershell -command mkdir ",directory_for_subjects,"/analysis_results/mds"))
  system(paste0("powershell -command mkdir ",directory_for_subjects,"/analysis_results/kmeans"))
  system(paste0("powershell -command mkdir ",directory_for_subjects,"/analysis_results/visualizations"))
  system(paste0("powershell -command mkdir ",directory_for_subjects,"/analysis_results/visualizations/tasks"))
  system(paste0("powershell -command mkdir ",directory_for_subjects,"/analysis_results/visualizations/subjects"))
  
  # read in subject meta data and format
  meta_data <- data.frame(Subject = subjects)

  # read in diagrams, order is rest1, rest2, emotion, gambling, language, motor, relational, social, wm
  diagrams <- lapply(meta_data$Subject,FUN = function(X){
    
    ret_list <- list()
    ret_list[[1]] <- read.csv(paste0(directory_for_subjects,"/",X,"/rfMRI_REST1_LR_diagram.csv"))
    ret_list[[1]]$X <- NULL
    ret_list[[2]] <- read.csv(paste0(directory_for_subjects,"/",X,"/rfMRI_REST1_RL_diagram.csv"))
    ret_list[[2]]$X <- NULL
    ret_list[[3]] <- read.csv(paste0(directory_for_subjects,"/",X,"/rfMRI_REST2_LR_diagram.csv"))
    ret_list[[3]]$X <- NULL
    ret_list[[4]] <- read.csv(paste0(directory_for_subjects,"/",X,"/rfMRI_REST2_RL_diagram.csv"))
    ret_list[[4]]$X <- NULL
    ret_list[[5]] <- read.csv(paste0(directory_for_subjects,"/",X,"/tfMRI_EMOTION_LR_diagram.csv"))
    ret_list[[5]]$X <- NULL
    ret_list[[6]] <- read.csv(paste0(directory_for_subjects,"/",X,"/tfMRI_EMOTION_RL_diagram.csv"))
    ret_list[[6]]$X <- NULL
    ret_list[[7]] <- read.csv(paste0(directory_for_subjects,"/",X,"/tfMRI_GAMBLING_LR_diagram.csv"))
    ret_list[[7]]$X <- NULL
    ret_list[[8]] <- read.csv(paste0(directory_for_subjects,"/",X,"/tfMRI_GAMBLING_RL_diagram.csv"))
    ret_list[[8]]$X <- NULL
    ret_list[[9]] <- read.csv(paste0(directory_for_subjects,"/",X,"/tfMRI_LANGUAGE_LR_diagram.csv"))
    ret_list[[9]]$X <- NULL
    ret_list[[10]] <- read.csv(paste0(directory_for_subjects,"/",X,"/tfMRI_LANGUAGE_RL_diagram.csv"))
    ret_list[[10]]$X <- NULL
    ret_list[[11]] <- read.csv(paste0(directory_for_subjects,"/",X,"/tfMRI_MOTOR_LR_diagram.csv"))
    ret_list[[11]]$X <- NULL
    ret_list[[12]] <- read.csv(paste0(directory_for_subjects,"/",X,"/tfMRI_MOTOR_RL_diagram.csv"))
    ret_list[[12]]$X <- NULL
    ret_list[[13]] <- read.csv(paste0(directory_for_subjects,"/",X,"/tfMRI_RELATIONAL_LR_diagram.csv"))
    ret_list[[13]]$X <- NULL
    ret_list[[14]] <- read.csv(paste0(directory_for_subjects,"/",X,"/tfMRI_RELATIONAL_RL_diagram.csv"))
    ret_list[[14]]$X <- NULL
    ret_list[[15]] <- read.csv(paste0(directory_for_subjects,"/",X,"/tfMRI_SOCIAL_LR_diagram.csv"))
    ret_list[[15]]$X <- NULL
    ret_list[[16]] <- read.csv(paste0(directory_for_subjects,"/",X,"/tfMRI_SOCIAL_RL_diagram.csv"))
    ret_list[[16]]$X <- NULL
    ret_list[[17]] <- read.csv(paste0(directory_for_subjects,"/",X,"/tfMRI_WM_LR_diagram.csv"))
    ret_list[[17]]$X <- NULL
    ret_list[[18]] <- read.csv(paste0(directory_for_subjects,"/",X,"/tfMRI_WM_RL_diagram.csv"))
    ret_list[[18]]$X <- NULL
    return(ret_list)
    
  })
  tasks <- c("rest1","rest2","emotion","gambling","language","motor","relational","social","wm")
  unlisted_diagrams <- list()
  for(i in 1:length(diagrams))
  {
    for(j in 1:18)
    {
      if(nrow(diagrams[[i]][[j]]) == 0)
      {
        unlisted_diagrams[[length(unlisted_diagrams) + 1]] <- data.frame(dimension = numeric(),birth = numeric(),death = numeric())
      }else
      {
        unlisted_diagrams[[length(unlisted_diagrams)+1]] <- diagrams[[i]][[j]] 
      }
    }
  }
  rm(diagrams)
  gc()
  
  # determine kernel parameter values
  print(paste0("Determining kernel parameters at ",Sys.time()))
  sigma_vals <- c(0.001,0.01,0.1)
  
  # find percentiles of Fisher information metrics between all diagrams
  D_fisher_0.001 <- distance_matrix(diagrams = unlisted_diagrams,dim = 1,distance = "fisher",sigma = 0.001)
  D_fisher_0.01 <- distance_matrix(diagrams = unlisted_diagrams,dim = 1,distance = "fisher",sigma = 0.01)
  D_fisher_0.1 <- distance_matrix(diagrams = unlisted_diagrams,dim = 1,distance = "fisher",sigma = 0.1)
  vals_0.001 <- D_fisher_0.001[which(upper.tri(D_fisher_0.001),arr.ind = T)]
  vals_0.001 <- vals_0.001[which(vals_0.001!=0)]
  vals_0.01 <- D_fisher_0.01[which(upper.tri(D_fisher_0.01),arr.ind = T)]
  vals_0.01 <- vals_0.01[which(vals_0.01!=0)]
  vals_0.1 <- D_fisher_0.1[which(upper.tri(D_fisher_0.1),arr.ind = T)]
  vals_0.1 <- vals_0.1[which(vals_0.1!=0)]
  t_vals_0.001 <- quantile(vals_0.001,probs = c(0.01,0.02,0.05,0.1,0.25,0.5))
  t_vals_0.01 <- quantile(vals_0.01,probs = c(0.01,0.02,0.05,0.1,0.25,0.5))
  t_vals_0.1 <- quantile(vals_0.1,probs = c(0.01,0.02,0.05,0.1,0.25,0.5))
  
  # combine
  kernel_parameters <- rbind(expand.grid(t = t_vals_0.001,sigma = sigma_vals),rbind(expand.grid(t = t_vals_0.01,sigma = sigma_vals),expand.grid(t = t_vals_0.1,sigma = sigma_vals)))
  
  # visualize all subject and task diagram distributions:
  print(paste0("Working on visualizing diagrams at ",Sys.time()))
  for(sigma in sigma_vals)
  {
    for(s in 1:nrow(meta_data))
    {
      png(filename = paste0(directory_for_subjects,"/analysis_results/visualizations/subjects/",meta_data$Subject[[s]],"_",sigma,".png"),width = 400,height = 450)
      visualize_group_of_diagrams(diagrams = unlisted_diagrams[(18*(s-1)+1):(18*s)],sigma = sigma,plot_title = paste0("Subject ",meta_data$Subject[[s]]))
      dev.off()
    }
    
    for(task in 1:length(tasks))
    {
      png(filename = paste0(directory_for_subjects,"/analysis_results/visualizations/tasks/",tasks[[task]],"_",sigma,".png"),width = 400,height = 450)
      visualize_group_of_diagrams(diagrams = unlisted_diagrams[c(seq(2*task-1,length(unlisted_diagrams),18),seq(2*task,length(unlisted_diagrams),18))],sigma = sigma,plot_title = paste0("Task ",tasks[[task]]))
      dev.off()
    }
  }
  
  # plot specifically the two representative loops of resting state 1 task
  analyze_two_loops(directory_for_subjects = directory_for_subjects,output_file_loop_1 = paste0(directory_for_subjects,"/loop1"),output_file_loop_2 = paste0(directory_for_subjects,"/loop2"),th = 0.75,sigma = 5)
  
  # plot surface parcellation
  view_parcellation(output_file = paste0(directory_for_subjects,"/parcellation"))
  
  # look for statistical differences in the number of loops in the diagrams between subjects and tasks
  num_loops <- unlist(lapply(unlisted_diagrams,FUN = function(X){return(nrow(X))}))
  loop_test_subjs <- matrix(data = 0,nrow = 100,ncol = 100)
  for(i in 1:(nrow(meta_data) - 1))
  {
    for(j in (i+1):nrow(meta_data))
    {
      v <- 0
      t <- t.test(x = num_loops[seq(18*(i-1)+1,18*i,1)],y = num_loops[seq(18*(j-1)+1,18*j,1)],paired = T,alternative = "two.sided")
      if(!is.na(t$p.value) & t$p.value < 0.05)
      {
        if(t$conf.int[[1]] + t$conf.int[[2]] < 0)
        {
          v <- -1 # lower mean in row than col
        }else
        {
          v <- 1 # higher mean in row than col
        }
      }
      loop_test_subjs[i,j] <- v
      loop_test_subjs[j,i] <- v
      
    }
  }
  
  loop_test_tasks <- matrix(data = 0,nrow = 9,ncol = 9)
  for(i in 1:8)
  {
    for(j in (i+1):9)
    {
      v <- 0
      t <- t.test(x = num_loops[c(seq(2*i-1,length(unlisted_diagrams),18),seq(2*i,length(unlisted_diagrams),18))],y = num_loops[c(seq(2*j-1,length(unlisted_diagrams),18),seq(2*j,length(unlisted_diagrams),18))],paired = T,alternative = "two.sided")
      if(!is.na(t$p.value) & t$p.value < 0.05)
      {
        if(t$conf.int[[1]] + t$conf.int[[2]] < 0)
        {
          v <- -1 
        }else
        {
          v <- 1
        }
      }
      loop_test_tasks[i,j] <- v
      loop_test_tasks[j,i] <- v
      
    }
  }
  
  # look for statistical differences in mean persistence between subjects and tasks
  mean_persistence <- unlist(lapply(unlisted_diagrams,FUN = function(X){
    
    if(nrow(X) == 0)
    {
      return(0)
    }
    return(mean(X[,3L] - X[,2L]))
    
  }))
  persistence_test_subjs <- matrix(data = 0,nrow = 100,ncol = 100)
  for(i in 1:(nrow(meta_data) - 1))
  {
    for(j in (i+1):nrow(meta_data))
    {
      v <- 0
      t <- t.test(x = mean_persistence[seq(18*(i-1)+1,18*i,1)],y = mean_persistence[seq(18*(j-1)+1,18*j,1)],paired = T,alternative = "two.sided")
      if(!is.na(t$p.value) & t$p.value < 0.05)
      {
        if(t$conf.int[[1]] + t$conf.int[[2]] < 0)
        {
          v <- -1
        }else
        {
          v <- 1
        }
      }
      persistence_test_subjs[i,j] <- v
      persistence_test_subjs[j,i] <- v
      
    }
  }
  
  persistence_test_tasks <- matrix(data = 0,nrow = 9,ncol = 9)
  for(i in 1:8)
  {
    for(j in (i+1):9)
    {
      v <- 0
      t <- t.test(x = mean_persistence[c(seq(2*i-1,length(unlisted_diagrams),18),seq(2*i,length(unlisted_diagrams),18))],y = mean_persistence[c(seq(2*j-1,length(unlisted_diagrams),18),seq(2*j,length(unlisted_diagrams),18))],paired = T,alternative = "two.sided")
      if(!is.na(t$p.value) & t$p.value < 0.05)
      {
        if(t$conf.int[[1]] + t$conf.int[[2]] < 0)
        {
          v <- -1
        }else
        {
          v <- 1
        }
      }
      persistence_test_tasks[i,j] <- v
      persistence_test_tasks[j,i] <- v
      
    }
  }

  # test if each task is different using the permutation test
  print(paste0("Working on permutation tests for tasks at ",Sys.time()))
  perm_test_mat_wass <- matrix(data = 1,nrow = 9,ncol = 9)
  perm_test_mat_bottleneck <- matrix(data = 1,nrow = 9,ncol = 9)
  for(i in 1:8)
  {
    for(j in (i+1):9)
    {
      perm_test_mat_wass[i,j] <- permutation_test(unlisted_diagrams[c(seq(2*i-1,length(unlisted_diagrams),18),seq(2*i,length(unlisted_diagrams),18))],unlisted_diagrams[c(seq(2*j-1,length(unlisted_diagrams),18),seq(2*j,length(unlisted_diagrams),18))],p = 2,iterations = 1000,dims = c(1))$p_values[[1]]
      perm_test_mat_wass[j,i] <- perm_test_mat_wass[i,j]
      perm_test_mat_bottleneck[i,j] <- permutation_test(unlisted_diagrams[c(seq(2*i-1,length(unlisted_diagrams),18),seq(2*i,length(unlisted_diagrams),18))],unlisted_diagrams[c(seq(2*j-1,length(unlisted_diagrams),18),seq(2*j,length(unlisted_diagrams),18))],p = Inf,iterations = 1000,dims = c(1))$p_values[[1]]
      perm_test_mat_bottleneck[j,i] <- perm_test_mat_bottleneck[i,j]
    }
  }
  write.csv(perm_test_mat_wass,paste0(directory_for_subjects,"/analysis_results/permutation_tests/tasks/perm_test_tasks_wass.csv"))
  write.csv(perm_test_mat_bottleneck,paste0(directory_for_subjects,"/analysis_results/permutation_tests/tasks/perm_test_tasks_bottleneck.csv"))
  
  # what about subject differences?
  print(paste0("Working on permutation tests for subjects at ",Sys.time()))
  perm_test_subj_wass <- matrix(data = 1,nrow = nrow(meta_data),ncol = nrow(meta_data))
  perm_test_subj_bottleneck <- matrix(data = 1,nrow = nrow(meta_data),ncol = nrow(meta_data))
  for(i in 1:(nrow(meta_data) - 1))
  {
    for(j in (i+1):nrow(meta_data))
    {
      perm_test_subj_wass[i,j] <- permutation_test(unlisted_diagrams[seq(18*(i-1)+1,18*i,1)],unlisted_diagrams[seq(18*(j-1)+1,18*j,1)],p = 2,iterations = 1000,dims = c(1))$p_values[[1]]
      perm_test_subj_bottleneck[i,j] <- permutation_test(unlisted_diagrams[seq(18*(i-1)+1,18*i,1)],unlisted_diagrams[seq(18*(j-1)+1,18*j,1)],p = Inf,iterations = 1000,dims = c(1))$p_values[[1]]
      perm_test_subj_wass[j,i] <- perm_test_subj_wass[i,j]
      perm_test_subj_bottleneck[j,i] <- perm_test_subj_bottleneck[i,j]
    }
  }
  write.csv(perm_test_subj_wass,paste0(directory_for_subjects,"/analysis_results/permutation_tests/subjects/perm_test_subjs_wass.csv"))
  write.csv(perm_test_subj_bottleneck,paste0(directory_for_subjects,"/analysis_results/permutation_tests/subjects/perm_test_subjs_bottleneck.csv"))

  # do independence tests between tasks and subjects
  print(paste0("Working on independence tests for tasks at ",Sys.time()))
  for(k in 1:nrow(kernel_parameters))
  {
    indep_mat <- matrix(data = 1,nrow = 9,ncol = 9)
    for(i in 1:8)
    {
      for(j in (i+1):9)
      {
        indep_mat[i,j] <- independence_test(g1 = unlisted_diagrams[c(seq(2*i-1,length(unlisted_diagrams),18),seq(2*i,length(unlisted_diagrams),18))],g2 = unlisted_diagrams[c(seq(2*j-1,length(unlisted_diagrams),18),seq(2*j,length(unlisted_diagrams),18))],dims = c(1),sigma = kernel_parameters[k,2],t = kernel_parameters[k,1])$p_values[[1]]
        indep_mat[j,i] <- indep_mat[i,j]
      }
    }
    write.csv(indep_mat,paste0(directory_for_subjects,"/analysis_results/independence_tests/tasks/indep_mat_task_t_",kernel_parameters[k,1],"_sigma_",kernel_parameters[k,2],".csv"))
  }

  # now check for dependence within subjects
  print(paste0("Working on independence tests for subjects at ",Sys.time()))
  for(k in 1:nrow(kernel_parameters))
  {
    indep_mat <- matrix(data = 1,nrow = nrow(meta_data),ncol = nrow(meta_data))
    for(i in 1:(nrow(meta_data) - 1))
    {
      for(j in (i+1):nrow(meta_data))
      {
        indep_mat[i,j] <- independence_test(g1 = unlisted_diagrams[seq(18*(i-1)+1,18*i,1)],g2 = unlisted_diagrams[seq(18*(j-1)+1,18*j,1)],dims = c(1),sigma = kernel_parameters[k,2],t = kernel_parameters[k,1])$p_values[[1]]
        indep_mat[j,i] <- indep_mat[i,j]
      }
    }
    write.csv(indep_mat,paste0(directory_for_subjects,"/analysis_results/independence_tests/subjects/indep_mat_subj_t_",kernel_parameters[k,1],"_sigma_",kernel_parameters[k,2],".csv"))
  }
  
  print(paste0("Working on comparing independence of tasks and subjects at ",Sys.time()))
  
  # now verify that some tasks/subjects were consistently dependent:
  indep_task_mats <- lapply(X = list.files(paste0(directory_for_subjects,"/analysis_results/independence_tests/tasks"),pattern = "indep_mat_task"),FUN = function(X){
    
    df <- read.csv(paste0(directory_for_subjects,"/analysis_results/independence_tests/tasks/",X))
    df$X <- NULL
    return(df)
    
  })
  indep_subj_mats <- lapply(X = list.files(paste0(directory_for_subjects,"/analysis_results/independence_tests/subjects"),pattern = "indep_mat_subj"),FUN = function(X){
    
    df <- read.csv(paste0(directory_for_subjects,"/analysis_results/independence_tests/subjects/",X))
    df$X <- NULL
    return(df)
    
  })
  
  # subset
  indep_task_mats_good <- list()
  indep_subj_mats_good <- list()
  indices_task_good <- c()
  indices_subj_good <- c()
  for(i in 1:nrow(kernel_parameters))
  {
    if(length(which(as.vector(indep_task_mats[[i]]) < 0.05)) > 0)
    {
      indep_task_mats_good[[length(indep_task_mats_good) + 1]] <- indep_task_mats[[i]]
      indices_task_good <- c(indices_task_good,i)
    }
    
    if(length(which(as.vector(indep_subj_mats[[i]]) < 0.05)) > 0)
    {
      indep_subj_mats_good[[length(indep_subj_mats_good) + 1]] <- indep_subj_mats[[i]]
      indices_subj_good <- c(indices_subj_good,i)
    }
  }
  
  # get percent significant across all kernel parameters for each p value calculation for tasks
  indep_task_mats_good <- lapply(indep_task_mats_good,FUN = function(X){
    
    # switch matrices to binary - 1 if p value less than 0.05 and 0 otherwise
    return(ifelse(test = X < 0.05,yes = 1,no = 0))
    
  })
  
  # take mean of all matrices
  indep_task_mat <- matrix(data = 0,nrow = 9,ncol = 9)
  for(i in 1:8)
  {
    for(j in (i+1):9)
    {
      indep_task_mat[i,j] <- mean(unlist(lapply(X = 1:length(indep_task_mats_good),FUN = function(X){return(indep_task_mats_good[[X]][i,j])})),na.rm = T)
      indep_task_mat[j,i] <- indep_task_mat[i,j]
    }
  }
  write.csv(indep_task_mat,paste0(directory_for_subjects,"/analysis_results/independence_tests/tasks/indep_task_mat_percent.csv"))
  
  # same but for subjects
  indep_subj_mats_good <- lapply(indep_subj_mats_good,FUN = function(X){
    
    return(ifelse(test = X < 0.05,yes = 1,no = 0))
    
  })
  indep_subj_mat <- matrix(data = 0,nrow = nrow(meta_data),ncol = nrow(meta_data))
  for(i in 1:(nrow(meta_data)-1))
  {
    for(j in (i+1):nrow(meta_data))
    {
      indep_subj_mat[i,j] <- mean(unlist(lapply(X = 1:length(indep_subj_mats_good),FUN = function(X){return(indep_subj_mats_good[[X]][i,j])})),na.rm = T)
      indep_subj_mat[j,i] <- indep_subj_mat[i,j]
    }
  }
  write.csv(indep_subj_mat,paste0(directory_for_subjects,"/analysis_results/independence_tests/subjects/indep_subj_mat_percent.csv"))
  
  # t tests for comparing percent kernel parameter non-independence for topologically different diagrams and not
  t.test(x = indep_subj_mat[which(perm_test_subj_bottleneck < 0.05,arr.ind = T)],y = indep_subj_mat[which(perm_test_subj_bottleneck >= 0.05,arr.ind = T)],alternative = "two.sided",paired = F)
  t.test(x = indep_subj_mat[which(perm_test_subj_wass < 0.05,arr.ind = T)],y = indep_subj_mat[which(perm_test_subj_wass >= 0.05,arr.ind = T)],alternative = "two.sided",paired = F)
  
  # run clustering to see if tasks or subjects can be discerned
  print(paste0("Working on clustering at ",Sys.time()))
  # parameter grid for k in kkmeans
  kvals <- c(2:18,25,50,100)
  for(i in 1:nrow(kernel_parameters))
  {
    kmeans_results <- data.frame(k = kvals,withinss = unlist(lapply(X = kvals,FUN = function(X){
      
      s <- tryCatch(expr = {
      
        R.utils::withTimeout(sum(diagram_kkmeans(diagrams = unlisted_diagrams,centers = X,dim = 1,t = kernel_parameters[i,1L],sigma = kernel_parameters[i,2L])$clustering@withinss),timeout = 60)
        },
               TimeoutException = function(ex){return(NA)},
               error = function(e){return(NA)})
      return(s)
      
    })))
    kmeans_results <- kmeans_results[which(is.na(kmeans_results$withinss) == F),]
    write.csv(kmeans_results,paste0(directory_for_subjects,"/analysis_results/kmeans/kmeans_results_t_",kernel_parameters[i,1],"_sigma_",kernel_parameters[i,2],".csv"))
  }
  
  # get optimal number of clusters according to min within cluster sum of squares
  kmeans_best_num_clusters <- unlist(lapply(X = list.files(paste0(directory_for_subjects,"/analysis_results/kmeans")),FUN = function(X){
    
    df <- read.csv(paste0(directory_for_subjects,"/analysis_results/kmeans/",X))
    df$X <- NULL
    return(df[which(df$withinss == min(df$withinss)),1L])
    
  }))
  
  # barplot distribution of best number of clusters
  barplot(table(kmeans_best_num_clusters))
  
  # compute cluster membership adjacency matrices for each kernel parameters
  # for three clusters then four clusters
  three_clusters <- list()
  for(X in 1:nrow(kernel_parameters))
  {
    ret <- tryCatch(expr = {
      
      R.utils::withTimeout(diagram_kkmeans(diagrams = unlisted_diagrams,dim = 1,t = kernel_parameters[X,1L],sigma = kernel_parameters[X,2L],centers = 3)$clustering@.Data,timeout = 60)
    },
    TimeoutException = function(ex){return(matrix(data = 1,nrow = 0,ncol = 1800))},
    error = function(e){return(matrix(data = 1,nrow = 0,ncol = 1800))})
    ret_mat <- matrix(data = 0,nrow = 1800,ncol = 1800)
    if(!is.matrix(ret))
    {
      for(i in 1:3)
      {
        inds <- which(ret == i)
        ret_mat[as.matrix(expand.grid(inds,inds))] <- 1
      }
    }
    three_clusters[[length(three_clusters) + 1]] <- ret_mat
  }
  
  four_clusters <- list()
  for(X in 1:nrow(kernel_parameters))
  {
    ret <- tryCatch(expr = {
      
      R.utils::withTimeout(diagram_kkmeans(diagrams = unlisted_diagrams,dim = 1,t = kernel_parameters[X,1L],sigma = kernel_parameters[X,2L],centers = 4)$clustering@.Data,timeout = 60)
    },
    TimeoutException = function(ex){return(matrix(data = 1,nrow = 0,ncol = 1800))},
    error = function(e){return(matrix(data = 1,nrow = 0,ncol = 1800))})
    ret_mat <- matrix(data = 0,nrow = 1800,ncol = 1800)
    if(!is.matrix(ret))
    {
      for(i in 1:4)
      {
        inds <- which(ret == i)
        ret_mat[as.matrix(expand.grid(inds,inds))] <- 1
      }
    }
    four_clusters[[length(four_clusters) + 1]] <- ret_mat
  }
  
  # compute mean cluster membership matrix for 3 and 4 clusters
  three_clusters <- Reduce("+",three_clusters)/nrow(kernel_parameters)
  four_clusters <- Reduce("+",four_clusters)/nrow(kernel_parameters)
  
  # plot histograms of percent same-cluster membership
  hist(three_clusters[which(upper.tri(three_clusters),arr.ind = T)],xlab = "Percent same cluster membership across all pairs of diagrams",main = "Cluster stability with k = 3")
  hist(four_clusters[which(upper.tri(four_clusters),arr.ind = T)],xlab = "Percent same cluster membership across all pairs of diagrams",main = "Cluster stability with k = 4")
  
  # cluster three_clusters and four_clusters to get static labels
  class(three_clusters) <- "kernelMatrix"
  class(four_clusters) <- "kernelMatrix"
  three_cluster_labels <- kernlab::kkmeans(x = three_clusters,centers = 3)@.Data
  # four_cluster_labels <- kernlab::kkmeans(x = four_clusters,centers = 4) # did not converge
  
  # analyze the number of rows in the diagrams in each cluster
  table(num_loops[which(three_cluster_labels == 1)])
  table(num_loops[which(three_cluster_labels == 2)])
  table(num_loops[which(three_cluster_labels == 3)])
  
  # analyze the task composition of the three clusters
  tasks_long <- rep(rep(tasks,each = 2),100)
  table(tasks_long[which(three_cluster_labels == 1)])/length(which(three_cluster_labels == 1))
  table(tasks_long[which(three_cluster_labels == 2)])/length(which(three_cluster_labels == 2))
  table(tasks_long[which(three_cluster_labels == 3)])/length(which(three_cluster_labels == 3))
  
  # visualize three found clusters
  visualize_group_of_diagrams(diagrams = unlisted_diagrams[which(three_cluster_labels == 1)],sigma = 0.001,plot_title = "Cluster 1")
  visualize_group_of_diagrams(diagrams = unlisted_diagrams[which(three_cluster_labels == 2)],sigma = 0.001,plot_title = "Cluster 2")
  visualize_group_of_diagrams(diagrams = unlisted_diagrams[which(three_cluster_labels == 3)],sigma = 0.001,plot_title = "Cluster 3")
  
  # embed all diagrams into 2D using mds
  print(paste0("Working on MDS embeddings at ",Sys.time()))
  mds_embedding_wass <- diagram_mds(unlisted_diagrams,k = 2,dim = 1)
  mds_embedding_bottleneck <- diagram_mds(unlisted_diagrams,k = 2,dim = 1)
  
  # add subject and task identifiers to each row
  mds_embedding_wass <- as.data.frame(mds_embedding_wass)
  mds_embedding_bottleneck <- as.data.frame(mds_embedding_bottleneck)
  mds_embedding_wass$Subject <- rep(meta_data$Subject,each = 18)
  mds_embedding_wass$task <- rep(rep(tasks,each = 2),nrow(meta_data))
  mds_embedding_bottleneck$Subject <- rep(meta_data$Subject,each = 18)
  mds_embedding_bottleneck$task <- rep(rep(tasks,each = 2),nrow(meta_data))
  
  write.csv(mds_embedding_wass,paste0(directory_for_subjects,"/analysis_results/mds/mds_wass.csv"))
  write.csv(mds_embedding_bottleneck,paste0(directory_for_subjects,"/analysis_results/mds/mds_bottleneck.csv"))
  
}

# EXECUTION
# with default 100 subjects
analyze_HCP(directory_for_subjects = "desired/path")
