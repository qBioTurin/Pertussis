library(parallel)
library(ggplot2)

compute <- function(x, par_1, par_2, distance_fname, out_fname, data_dir, fun, config, reference)
{
  # Load the distance function
  source(distance_fname)
  trace <- read.csv(file = paste0(data_dir, out_fname, "-", fun, "-", x, ".trace"), sep = "", header = TRUE)
  return(data.frame(score = do.call(tools::file_path_sans_ext(basename(distance_fname)),
                                    list(reference, trace)),
                    par_1 = config[[par_1$cfg.idx]][[x]][[3]][[par_1$par.idx]],
                    par_2 = config[[par_2$cfg.idx]][[x]][[3]][[par_2$par.idx]]))
}


# Script configuration
distance_fname <- "~/Rfunction/msqd.R"
data_dir <- "~/results_sensitivity_analysis_reduced_probability_range/"
out_fname <- "Pertussis"
reference <- as.data.frame(t(read.csv("~/input/reference_data.csv", header = FALSE, sep = "")))
fun <- "sensitivity"
par_1 <- data.frame(cfg.idx=5, par.idx=2)
par_2 <- data.frame(cfg.idx=5, par.idx=3)
par_processors <- 16
# Load configuration
load(paste0(data_dir, "prcc_", out_fname, "-", fun, ".RData"))
load(paste0(data_dir,paste0(out_fname, "-", fun, ".RData")))
# Get the number of run
n_run <- length(config[[1]])
# Make cluster
cl <- makeForkCluster(nnodes = par_processors)

# apply the 
rnk<- parLapply(X = c(1:n_run),
                fun = compute,
                par_1,
                par_2,
                distance_fname,
                out_fname,
                data_dir,
                fun,
                config,
                reference,
                cl=cl)
stopCluster(cl = cl)

rnk <- do.call("rbind", rnk)

names(rnk) <- c("score","prob_infection_S", "prob_infection_R_l1")

rnk<- rnk[order(rnk$score, decreasing = TRUE),]
plt <- ggplot(data = rnk,
              aes(x=prob_infection_S,y=prob_infection_R_l1,col= (min(score)-score)/ (min(score)-max(score)) ) )+
  geom_point(size=2.5)+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=18),
        legend.title=element_text(size=20,face="bold")) +
  labs(title="")+
  scale_colour_gradientn("Error",colours = c("black","deepskyblue2","cyan"),
                         values= c(0,.000001875,1),
                         breaks=c(0,.999),
                         labels=c("min","max"))

ggsave(plot = plt, filename = "~/Plots/scatterplt.pdf", device = "pdf", dpi="retina", width=10, height=7.5, units="in")