coarse.cell.types <-
  c("B.cells", "CD4.T.cells", "CD8.T.cells", "NK.cells", "neutrophils", "monocytic.lineage",
    "fibroblasts", "endothelial.cells")

output.plot <- function(g, file.suffix, plot.types = c("pdf", "png", "svg"),
                        pdf.delta.width = 1, pdf.delta.height = 1,
                        png.delta.width = 1, png.delta.height = 1,
                        svg.delta.width = 1, svg.delta.height = 1) {
  pdf.width <- 7
  pdf.height <- 7
  png.width <- 480
  png.height <- 480
  svg.width <- 7
  svg.height <- 7
  
  if("pdf" %in% plot.types) {
    pdf(file = paste0(file.suffix, ".pdf"), width = pdf.width * pdf.delta.width, height = pdf.height * pdf.delta.height)
    print(g)
    d <- dev.off()
  }
  if("png" %in% plot.types) {
    png(file = paste0(file.suffix, ".png"), width = png.width * png.delta.width, height = png.height * png.delta.height)
    print(g)
    d <- dev.off()
  }
  if("svg" %in% plot.types) {
    svg(file = paste0(file.suffix, ".svg"), width = svg.width * svg.delta.width, height = svg.height * png.delta.height)
    print(g)
    d <- dev.off()
  }
}


## See https://stackoverflow.com/questions/27803710/ggplot2-divide-legend-into-two-columns-each-with-its-own-title
plot.anno.heatmap.with.multiple.legends <-
  function(df, id.col, anno.columns, anno.pals) {
    
    suppressPackageStartupMessages(p_load("RColorBrewer"))
    df <- df[, c(id.col, anno.columns)]
    
    ## Assume annotations are characters
    ## NB: id.col is a factor
    for(col in c(anno.columns)) {
      df[, col] <- as.character(df[, col])
    }
    
    columns <- 1:length(anno.columns)
    names(columns) <- anno.columns
    
    color.vecs <-
      llply(columns,
            .fun = function(idx) {
              anno.col <- anno.columns[idx]
              vec <- unique(df[, anno.col])
              len <- length(vec)
              colors <- brewer.pal(len, anno.pals[idx])
              names(colors) <- vec
              colors
            })
    
    all.colors <- Reduce("c", color.vecs)
    names(all.colors) <- Reduce("c", unlist(lapply(color.vecs, names)))
    
    names(anno.columns) <- anno.columns
    anno.df <- ldply(anno.columns,
                     .fun = function(anno.col) {
                       data.frame(val = df[, anno.col], id = df[, id.col])
                     })
    colnames(anno.df)[1] <- "type"
    
    full.plot <-
      ggplot(anno.df, aes(y = id, x = type, fill = val)) + geom_tile() +
      scale_fill_manual(values = all.colors) +
      theme(legend.position="none")
    
    full.plot <- full.plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                                   axis.ticks.y = element_blank(), text = element_text(size = 18),
                                   axis.text.x = element_text(angle = 45, hjust = 1),
                                   axis.title.x = element_blank())
    
    
    legends <-
      llply(anno.columns, .parallel = FALSE,
            .fun = function(anno.col) {
              flag <- anno.df$type == anno.col
              leg.colors <- all.colors[names(all.colors) %in% anno.df[flag,"val"]]
              g <- ggplot(anno.df[flag, ], aes_string(x = "id", y = "type", fill = "val"))
              g <- g + geom_tile()
              g <- g + scale_fill_manual(values = leg.colors, name = anno.col)
            })
    return(list("full.plot" = full.plot, "legends" = legends))
  }

get.method.annotations <- function() {
  data_path_to_read <-  paste0(validation_path_base, "/deconvolution-method-annotation.csv")  #"syn23395242"
  method.anno <- read.csv(data_path_to_read, sep=",")
  round.col <- "submission"
  subchallenge.col <- "subchallenge"
  
  method.anno$method.type <- as.character(method.anno$method.type)
  method.anno$output.type <- as.character(method.anno$output.type)
  method.anno[, round.col] <- as.character(method.anno[, round.col])
  flag <- is.na(method.anno[, round.col])
  method.anno[flag, round.col] <- NA
  method.anno[, subchallenge.col] <- as.character(method.anno[, subchallenge.col])
  flag <- is.na(method.anno[, subchallenge.col])
  method.anno[flag, subchallenge.col] <- NA
  
  method.rename.list <-
    list("NNLS" = "NNLS",
         "summary" = "SUM",
         "other" = "OTH",
         "other regression" = "REG",
         "unknown" = "UNK",
         "SVR" = "SVR",
         "DNN" = "DNN",
         "ensemble" = "ENS",
         "NMF" = "NMF",
         "probabilistic inference" = "PI")
  method.rename.df <- data.frame(method.type = names(method.rename.list), Method = as.character(method.rename.list))
  
  output.rename.list <-
    list("fraction" = "Frac",
         "proportion" = "Prop",
         "normalized.score" = "Norm",
         "score" = "Score")
  output.rename.df <- data.frame(output.type = names(output.rename.list), Output = as.character(output.rename.list))
  
  method.anno <- merge(method.anno, method.rename.df)
  method.anno <- merge(method.anno, output.rename.df)
  method.anno$Output <- as.character(method.anno$Output)
  method.anno$Method <- as.character(method.anno$Method)
  method.anno
}


get.round.specific.annotations <- function(method.anno, round) {
  method.anno.round <-
    ddply(method.anno,
          .variables = c(method.name.col, subchallenge.col),
          .fun = function(df) {
            flag <- is.na(df[, round.col]) | (df[, round.col] == round)
            if(any(flag)) {
              ret <- df[flag,,drop=F]
              if(nrow(ret) != 1) { print(ret); print(df); stop(paste0("Wrong number of rows for ",
                                                                      df[1,method.name.col], " in ", df[1, subchallenge.col], "!\n")); }
              return(ret)
            }
            flag <- df[, round.col] == "latest"
            ret <- df[flag,,drop=F]
            if(nrow(ret) != 1) { print(ret); print(df); stop(paste0("Wrong number of rows2 for ",
                                                                    df[1,method.name.col], " in ", df[1, subchallenge.col], "!\n")); }
            return(ret)
          })
  method.anno.round
}

get.cell.type.translation <- function() {
  list("naive.B.cells" = "naive B",
       "memory.B.cells" = "memory B",
       "fibroblasts" = "fibroblasts",
       "neutrophils" = "neutrophils",
       "endothelial.cells" = "endothelial",
       "monocytes" = "monocytes",
       "NK.cells" = "NK",
       "macrophages" = "macrophages",
       "memory.CD8.T.cells" = "memory CD8 T",
       "regulatory.T.cells" = "Tregs",
       "B.cells" = "B",
       "naive.CD8.T.cells" = "naive CD8 T",
       "monocytic.lineage" = "monocytic lineage",
       "naive.CD4.T.cells" = "naive CD4 T",
       "myeloid.dendritic.cells" = "myeloid DCs",
       "CD8.T.cells" = "CD8 T",
       "memory.CD4.T.cells" = "memory CD4 T",
       "CD4.T.cells" = "CD4 T")         
}

rename.cell.types <- function(df, from.col = "cell.type", to.col = "cell.type") {
  cell.type.trans <- get.cell.type.translation()
  for(nm in names(cell.type.trans)) {
    flag <- (grepl(df[, from.col], pattern = nm))
    df[flag, to.col] <- cell.type.trans[[nm]]
  }
  df
}

safe.merge <- function(x, y, ...) {
  
  orig.nrow <- nrow(x)
  x <- merge(x, y, ...)
  new.nrow <- nrow(x)
  if(orig.nrow != new.nrow) {
    stop("Changed rows\n")
  }
  x
}

## Assign team names and rounds to results
assign.result.team.names.and.rounds <- function(res, error.fun = stop) {
  ## Read in the final predictions from the competitive phase (i.e., when coarse- and fine-grained
  ## datasets differed)
  synId <- "syn22149603"
  obj <- synGet(synId, downloadFile=TRUE)
  res.comp <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)
  
  ## Use the competitive phase results to translate objectIds to teams / submitters
  ## The final name in the "path" of the repo_name for the post-competitive phase is the
  ## objectId in the competitive phase (except for the baselines)
  cat("Assigning objectId\n")
  map <- data.frame(repo_name = as.character(unique(res$repo_name)), stringsAsFactors = FALSE)
  
  map$objectId <-
    unlist(lapply(map$repo_name,
                  function(str) {
                    strs <- unlist(strsplit(str, "/"))
                    strs[length(strs)]
                  }))
  
  cat("Assigning comparator\n")
  map$comparator <-
    unlist(lapply(as.character(map$objectId),
                  function(str) {
                    comps <- get.comparators()
                    for(comp in comps) {
                      if(grepl(str, pattern=comp, ignore.case=TRUE)) { return(TRUE) }
                    }
                    return(FALSE)
                  }))
  
  res <- res[, !(colnames(res) %in% c("objectId", "comparator"))]
  res <- safe.merge(res, map, by = c("repo_name"))
  
  cat("Ensuring\n")
  
  ## Ensure that all objectIds match between competitive and post-competitive
  flag <- grepl(res.comp$repo_name, pattern="baseline") | (res.comp$objectId %in% res$objectId)
  if(!all(flag)) {
    print(table(flag))
    print(unique(res.comp[!flag, c("repo_name", "objectId", "submitterId")]))
    l_ply(unique(res.comp[!flag, "submitterId"]),
          .fun = function(id) print(translate.submitterId(id)))
    error.fun("Some competitive objectIds are not in post-competitive results\n")
  }
  
  flag <- (res$comparator == TRUE) | (res$objectId %in% res.comp$objectId)
  if(!all(flag)) {
    error.fun("Some post-competitive objectIds are not in competitive results\n")
  }
  
  cat("Merging with res.comp\n")
  res <- safe.merge(res, unique(res.comp[, c("objectId", "submitterId")]), all.x=TRUE)
  
  ## Assign the team name
  
  cat("Defining team name tbl\n")
  team.name.tbl <- unique(res[, c("objectId", "subchallenge", "submitterId", "comparator")])
  flag <- !(team.name.tbl$comparator == TRUE)
  
  team.name.tbl$method.name <- NA
  team.name.tbl[flag, "method.name"] <-
    unlist(lapply(as.character(team.name.tbl[flag, "submitterId"]),
                  function(str) translate.submitterId(str)))
  
  team.name.tbl <-
    assign.baseline.names(team.name.tbl, from.col = "objectId", to.col = "method.name")
  
  
  team.name.tbl <- simplify.submitter.names(team.name.tbl, col = "method.name")
  
  ## Assign the round
  team.name.tbl <-
    assign.submission.rounds(team.name.tbl, object.id.col = "objectId",
                             context.cols = c("subchallenge", "submitterId"),
                             method.name.col = "method.name")
  
  cat("Merging team name tbl\n")
  print(colnames(team.name.tbl))
  res <- merge(res, team.name.tbl, all.x = TRUE)
  cat("Done merging team name tbl\n")
  
  res
}

get.comparators <- function() {
  c("cibersort", "mcp", "quantiseq", "xcell", "epic", "timer", "cibersortx")
}

get.comparators.cap <- function() {
  c("CIBERSORT", "MCP-counter", "quanTIseq", "xCell", "EPIC", "TIMER", "CIBERSORTx")
}

get.baseline.translation <- function() {
  baseline.method.trans <-
    list("baseline_method1" = "CIBERSORT",
         "baseline_method2" = "MCP-counter",
         "baseline_method3" = "quanTIseq",
         "baseline_method4" = "xCell",
         "baseline_method5" = "EPIC",
         "baseline_method6" = "TIMER",
         "baseline_method7" = "CIBERSORTx",
         "cibersortx" = "CIBERSORTx",
         "CIBERSORTx" = "CIBERSORTx",             
         "cibersort" = "CIBERSORT",
         "mcp" = "MCP-counter",
         "quantiseq" = "quanTIseq",
         "xcell" = "xCell",
         "epic" = "EPIC",
         "timer" = "TIMER")
  baseline.method.trans
}

assign.baseline.names <- function(df, from.col = "repo_name", to.col = "repo_name") {
  baseline.method.trans <- get.baseline.translation()
  for(nm in names(baseline.method.trans)) {
    flag <- (grepl(df[, from.col], pattern = nm))
    df[flag, to.col] <- baseline.method.trans[[nm]]
  }
  df
}

define.baseline.method.flag <- function(res, method.name.col) {
  baseline.method.trans <- get.baseline.translation()    
  baseline.methods <- unname(unlist(baseline.method.trans))
  baseline.method.flag <- grepl(res[, method.name.col], pattern="baseline")
  for(baseline.method in baseline.methods) {
    baseline.method.flag <- baseline.method.flag |
      grepl(res[, method.name.col], pattern=baseline.method, ignore.case=TRUE)
  }
  baseline.method.flag
}

## Assign rounds to assignments based on objectId < object Id Y then X was
## submitted before Y

## context.cols: the combination of columns that select all versions/submissions of a method
## i.e., for the raw prediction results, context.cols = c(subchallenge.col, submitter.id.col)
## For the leaderboard results (already separated by leaderboard), it is
## context.cols = c(submitter.id.col)
## method.name.col: probably repo_name
assign.submission.rounds <- function(res, object.id.col, context.cols, method.name.col,
                                     assign.latest = TRUE) {
  baseline.method.flag <-
    define.baseline.method.flag(res, method.name.col)
  tmp <- unique(res[!baseline.method.flag, c(context.cols, object.id.col)])
  ret <-
    ddply(tmp,
          .variables = context.cols,
          .fun = function(df) {
            if(length(unique(df[, object.id.col])) != nrow(df)) {
              stop("Was only expecting unique object IDs\n")
            }
            o <- order(df[, object.id.col], decreasing = FALSE)
            df <- df[o, ]
            ret.df <- data.frame(df, submission = 1:nrow(df))
            if(assign.latest) {
              ret.df <- rbind(ret.df, ret.df[nrow(ret.df),,drop=F])
              ret.df[nrow(ret.df), "submission"] <- "latest"
            }
            ret.df
          })
  ret <- merge(res[!baseline.method.flag, ], ret, all.x = TRUE)
  
  ## Handle baselines, where all submissions of the same method
  ## should have the same method name
  tmp <- unique(res[baseline.method.flag, c(context.cols, object.id.col, method.name.col)])
  baseline.ret <-
    ddply(tmp,
          .variables = c(context.cols, method.name.col),
          .fun = function(df) {
            if(length(unique(df[, object.id.col])) != nrow(df)) {
              stop("Was only expecting unique object IDs\n")
            }
            o <- order(df[, object.id.col], decreasing = FALSE)
            df <- df[o, ]
            ret.df <- data.frame(df, submission = 1:nrow(df))
            if(assign.latest) {
              ret.df <- rbind(ret.df, ret.df[nrow(ret.df),,drop=F])
              ret.df[nrow(ret.df), "submission"] <- "latest"
            }
            ret.df
          })
  baseline.ret <- merge(res[baseline.method.flag, ], baseline.ret, all.x = TRUE)    
  
  rbind(ret, baseline.ret)
}

simplify.submitter.names <- function(df, col = "submitter") {
  df[, col] <- as.character(df[, col])
  trans <-
    list("Northwestern Polytechnical University" = "NPU",
         "TJU and the renegade mouse" = "TJU")
  for(nm in names(trans)) {
    flag <- grepl(df[, col], pattern = nm)
    df[flag, col] <- trans[[nm]]
  }
  df
  
}

# ictd and cancer deconv were both submitted by the same group and
# are essentially the same method. 
# drop ICTD and rename Cancer_Deconv to ICTD
correct.ictd.and.cancer.deconv <- function(df, col = "method.name") {
  df[, col] <- as.character(df[, col])
  flag <- !is.na(df[,col]) & (df[,col] == "ICTD")
  df <- df[!flag,]
  flag <- !is.na(df[,col]) & (df[,col] == "Cancer_Decon")
  df[flag,col] <- "ICTD"
  df
}

translate.submitterId <- function(submitterId) {
  tryCatch(synGetTeam(submitterId)$name,
           error = function(e) synGetUserProfile(submitterId)$userName)
}

ggplot_smooth_scatter <- function(data, mapping, pwr = 0.25, n = 200){
  p <- ggplot(data = data, mapping = mapping) + 
    stat_density2d(aes(fill=..density..^pwr), geom="tile", contour = FALSE, n = n) +
    scale_fill_continuous("", low = "white", high = "dodgerblue4")
  ##        scale_fill_gradientn(colours=rainbow(100))
  p
}

## Calculate the mean / variance loess fit for the genes (rows) in mat
## Expression values should be in log space
calc.mean.variance <- function(mat) {
  means <- unlist(apply(mat, 1, mean))
  sds <- unlist(apply(mat, 1, sd))
  x <- means
  y <- sqrt(sds)
  loess.fit <- loess(y~x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct")
  se <- FALSE
  loess.pred <- predict(loess.fit, data.frame(x=as.numeric(x)), se=se)
  
  loess.df <- data.frame(mean = x, sqrt.std = y, 
                         residual = as.numeric(residuals(loess.fit)))
  if(se) {
    loess.df$sqrt.std.pred = as.numeric(loess.pred$fit)
  } else {
    loess.df$sqrt.std.pred = as.numeric(loess.pred)
  }
  
  ## See from stats:::predLoess that se.fit is s(x) = s * \sqrt(\sum_{i=1}^n l_i^2(x)) on page 4 of Cleveland and Grosse (1991) 
  ## Statistics and Computing.  There the authors state that [ g^(x) - g(x) ] / s(x) has a t-distribution.
  ## g(x) is the original function and g^(x) is the loess prediction, so that g^(x) - g(x) is the residual.
  ## We considered selecting highly-variable genes based on the standardized residuals, but this penalizes areas of the
  ## curve in which the prediction is poor because of few points.  So, just defined highly variable points based
  ## on residual.
  ##
  ## Define standardized residuals
  if(se) { 
    loess.df$se.fit <- as.numeric(less.pred$se.fit)
    loess.df$residual.std <- loess.df$residual / loess.df$se.fit
  }
  loess.df <- loess.df[order(loess.df$mean), ]
  
  loess.df
}


## Assumes that expressed genes are represented by the maximal-density
## peak to the right of zero.
z_fpkm<-function(log.expr){
  if(all(log.expr>0)) stop('Input not log2 transformed.')
  pos.log.expr <- log.expr[log.expr > 0]
  density.peak <- max(density(pos.log.expr,na.rm=T)$y)
  my<-density(pos.log.expr,na.rm=T)$x[which.max(density(pos.log.expr,na.rm=T)$y)]
  U<-mean(log.expr[log.expr>my],na.rm=T)
  sigma<-(U-my)*(.5*pi)^.5
  z<-(log.expr-my)/sigma
  z[z< -3]<-NA
  return(list("mean" = my, "sd" = sigma, "z" = z, "density.peak" = density.peak))
}

z_fpkm_rightmost_peak <-function(log.expr, min.boundary = 0, max.boundary = 6){
  if(all(log.expr>0)) stop('Input not log2 transformed.')
  use.derivative.method <- TRUE
  ## I use this first derivative approach, because otherwise I need to closely
  ## bracket the peak (between min/max boundaries).  e.g., the max density
  ## in TCGA sample TCGA−AB−2929−03A−01T−0735−13 occurs at 0, whereas the
  ## desired peak occurs at 5. The latter is more easily found using the
  ## derivative approach. The max density approach would require setting the
  ## min boundary in the trough between peaks or up the left hand side of the
  ## max density peak.
  my <- NA
  density.peak <- NA
  den <- density(log.expr, na.rm=TRUE)
  df <- data.frame(x = den$x, y = den$y)
  if(use.derivative.method) {
    sp <- splinefun(den$x, den$y)
    df$deriv = sp(den$x, deriv=1)
    df <- subset(df, (x > min.boundary) & (x < max.boundary))
    for(i in (nrow(df)-1):1) {
      if(sign(df$deriv[i+1]) != sign(df$deriv[i])) {
        my <- df$x[i]
        density.peak <- df$y[i]
        break
      }
    }
  } else {
    df <- subset(df, (x > min.boundary) & (x < max.boundary))
    flag <- which.max(df$y)
    density.peak <- df$y[flag]
    my <- df$x[flag]
  }
  U<-mean(log.expr[log.expr>my],na.rm=T)
  sigma<-(U-my)*(.5*pi)^.5
  z<-(log.expr-my)/sigma
  z[z< -3]<-NA
  return(list("mean" = my, "sd" = sigma, "z" = z, "density.peak" = density.peak))
}

calc.density.of.expressed.genes <- function(log.expr, num.sds = NULL, z = NULL) {
  if(is.null(z)) {
    z <- z_fpkm(log.expr)
  }
  den <- density(log.expr)
  data.density <- data.frame(x = den$x, y = den$y)
  dn <- dnorm(den$x, mean = z$mean, sd = z$sd)
  density.peak <- den$y[which.min((den$x - z$mean)^2)]
  fit <- data.frame(x = den$x, y = density.peak * dn / max(dn)) 
  lst <- list("density" = data.density, "fit" = fit, "z" = z)
  lst
}

plot.density.of.expressed.genes <- function(log.expr, num.sds = NULL, z = NULL) {
  suppressPackageStartupMessages(p_load(latex2exp))
  lst <- calc.density.of.expressed.genes(log.expr, num.sds, z)
  data.density <- lst[["density"]]
  fit <- lst[["fit"]]
  z <- lst[["z"]]
  g <- ggplot()
  g <- g + geom_line(data = data.density, aes(x = x, y = y))
  g <- g + geom_line(data = fit, aes(x = x, y = y), linetype = "dashed", col = "blue")
  if(!is.null(num.sds)) {
    g <- g + geom_vline(xintercept = z$mean - num.sds * z$sd, linetype = "dashed")
  }
  ##  g <- g + xlab("Log CPM")
  g <- g + xlab(TeX("$Log_2(CPM)$"))
  g <- g + ylab("Density")
  g
}


get.deconvolution.genes <- function() {
  suppressPackageStartupMessages(p_load("immunedeconv"))
  suppressPackageStartupMessages(p_load("MCPcounter"))
  
  xcell.deconv.genes <-
    sort(unique(unlist(lapply(xCell.data$signatures, function(x) x@geneIds))))
  
  lm22 <- NULL
  files <- c("/home/bwhite/LM22.txt", "/Users/Brian/Downloads/new-cibersort-code/LM22.txt", "/home/whitebr/LM22.txt")
  for(file in files) {
    if(file.exists(file)) {
      lm22 <- read.table(file, sep="\t", header=TRUE)
    }
  }
  cibersort.deconv.genes <- as.character(lm22$Gene.symbol)
  mcp.deconv.genes <-
    sort(unique(
      read.table(curl:::curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"), 
                 sep = "\t", stringsAsFactors = FALSE, header = TRUE, 
                 colClasses = "character", check.names = FALSE)$`HUGO symbols`))
  
  epic.deconv.genes <- sort(unique(c(TRef$sigGenes, BRef$sigGenes)))
  
  signame <- "TIL10"
  sig.mat.file <- system.file("extdata", "quantiseq",
                              paste0(signame, 
                                     "_signature.txt"), package = "immunedeconv", mustWork = TRUE)
  df <- read.table(sig.mat.file, sep="\t", header=TRUE)
  quantiseq.deconv.genes <- as.character(unique(df$ID))
  
  deconv.genes <-
    sort(unique(c(cibersort.deconv.genes, mcp.deconv.genes, epic.deconv.genes, quantiseq.deconv.genes)))
  
  deconv.genes
}

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

fix.string <- function(str) {
  str <- gsub(str, pattern="\\.", replacement=" ")
  str <- unlist(strsplit(str, split=" "))
  str <- lapply(str, function(x) firstup(x))
  str <- paste0(str, collapse=" ")
  str
}

get.in.silico.metadata <- function() {
  
  synId <- "syn21773017"
  metadata.file <- synGet(synId, downloadFile = TRUE)$path
  df <- read.table(metadata.file, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  return(df)
  
  l <-
    list(
      "A" = list("mixture.type" = "Random", "tumor.type" = NA, "subchallenge" = "fine"),
      "B" = list("mixture.type" = "Biological", "tumor.type" = NA, "subchallenge" = "fine"),            
      "G" = list("mixture.type" = "Random", "tumor.type" = NA, "subchallenge" = "coarse"),
      "H" = list("mixture.type" = "Biological", "tumor.type" = NA, "subchallenge" = "coarse"),
      "C" = list("mixture.type" = "Random", "tumor.type" = "CRC", "subchallenge" = "fine"),
      "E" = list("mixture.type" = "Biological", "tumor.type" = "CRC", "subchallenge" = "fine"),            
      "D" = list("mixture.type" = "Random", "tumor.type" = "BRCA", "subchallenge" = "fine"),
      "F" = list("mixture.type" = "Biological", "tumor.type" = "BRCA", "subchallenge" = "fine"),            
      "I" = list("mixture.type" = "Random", "tumor.type" = "CRC", "subchallenge" = "coarse"),
      "K" = list("mixture.type" = "Biological", "tumor.type" = "CRC", "subchallenge" = "coarse"),            
      "J" = list("mixture.type" = "Random", "tumor.type" = "BRCA", "subchallenge" = "coarse"),
      "L" = list("mixture.type" = "Biological", "tumor.type" = "BRCA", "subchallenge" = "coarse")
    )
  df <- ldply(l, .fun = function(entry) as.data.frame(entry, stringsAsFactors = FALSE))
  colnames(df)[1] <- "dataset.name"
  df
}

get.validation.metadata <- function() {
  ## Our original admixture specification includes the vendor for each sample
  synId <- "syn21577258"
  obj <- synGet(synId, downloadFile = TRUE)
  bm1 <- read.xlsx(obj$path, sheet = 1)
  bm1 <- bm1[-c(1:3),]
  bm2 <- read.xlsx(obj$path, sheet = 2)
  bm2 <- bm2[-c(1:3),]
  rm1 <- read.xlsx(obj$path, sheet = 3)
  rm1 <- rm1[-c(1:3),]
  rm2 <- read.xlsx(obj$path, sheet = 4)
  rm2 <- rm2[-c(1:3),]
  
  lst <- list(bm1, bm2, rm1, rm2)
  
  df <-
    do.call(rbind, lapply(lst,
                          function(entry) {
                            label = ifelse(entry$CRC > 0, "CRC",
                                           ifelse(entry$breast > 0, "BRCA", ""))                              
                            data.frame(id = entry[,1], tumor.type = label)
                          }))
  df$batch <- "empty"
  flag <- df$id %in% bm1[,1]
  df$batch[flag] <- "BM1"
  flag <- df$id %in% bm2[,1]
  df$batch[flag] <- "BM2"
  flag <- df$id %in% rm1[,1]
  df$batch[flag] <- "RM1"
  flag <- df$id %in% rm2[,1]
  df$batch[flag] <- "RM2"
  
  df$mixture.type <- "empty"
  flag <- ( df$id %in% bm1[,1] ) | ( df$id %in% bm2[,1] )
  ##    df$mixture.type[flag] <- "BM"
  df$mixture.type[flag] <- "Biological"
  flag <- ( df$id %in% rm1[,1] ) | ( df$id %in% rm2[,1] )
  ## df$mixture.type[flag] <- "RM"
  df$mixture.type[flag] <- "Random"
  
  df$dataset <- "empty"
  flag <- ( df$tumor.type == "BRCA" ) & ( df$mixture.type == "Biological" )
  df$dataset[flag] <- "DS1"
  flag <- ( df$tumor.type == "CRC" ) & ( df$mixture.type == "Biological" )
  df$dataset[flag] <- "DS2"
  flag <- ( df$tumor.type == "BRCA" ) & ( df$mixture.type == "Random" )
  df$dataset[flag] <- "DS3"
  flag <- ( df$tumor.type == "CRC" ) & ( df$mixture.type == "Random" )
  df$dataset[flag] <- "DS4"
  
  df
}

## This is the ground truth we _specified_ to the Stanford core, not that we actually had created
get.specified.validation.ground.truth <- function() {
  ## Our original admixture specification includes the vendor for each sample
  synId <- "syn21577258"
  obj <- synGet(synId, downloadFile = TRUE)
  bm1 <- read.xlsx(obj$path, sheet = 1)
  bm1 <- bm1[-c(1:3),]
  bm2 <- read.xlsx(obj$path, sheet = 2)
  bm2 <- bm2[-c(1:3),]
  rm1 <- read.xlsx(obj$path, sheet = 3)
  rm1 <- rm1[-c(1:3),]
  rm2 <- read.xlsx(obj$path, sheet = 4)
  rm2 <- rm2[-c(1:3),]
  
  lst <- list(bm1, bm2, rm1, rm2)
  
  df <-
    do.call(rbind, lapply(lst,
                          function(entry) {
                            rownames(entry) <- entry[,1]
                            entry <- entry[,-1]
                            ret <- melt(as.matrix(entry))
                            colnames(ret) <- c("sample.id", "cell.type", "measured")
                            ret$sample.id <- as.character(ret$sample.id)
                            ret$cell.type <- as.character(ret$cell.type)
                            ret$measured <- as.numeric(as.character(ret$measured))
                            ret
                          }))
  
  md <- get.validation.metadata()
  colnames(md) <- c("sample.id", "tumor.type", "batch", "mixture.type", "dataset.name")
  df <- merge(df, md)
  
  df
}

get.purified.sample.translation.table <- function() {
  fine.grained.translation <- list(
    "Naive_B_cells_1" = "naive.B.cells",
    "Macrophages_2" = "macrophages",
    "Dendritic_cells_1" = "myeloid.dendritic.cells",
    "Macrophages_1" = "macrophages",
    "Fibroblasts" = "fibroblasts",
    "Endothelial_cells" = "endothelial.cells",
    "Memory_CD8_T_cells_1" = "memory.CD8.T.cells",
    "Monocytes_1" = "monocytes",
    "Neutrophils_2" = "neutrophils",
    "Dendritic_cells_2" = "myeloid.dendritic.cells",
    "NK_cells_2" = "NK.cells",
    "Monocytes_2" = "monocytes",
    "NK_cells_1" = "NK.cells",
    "Memory_CD8_T_cells_2" = "memory.CD8.T.cells",
    "Memory_CD4_T_cells_2" = "memory.CD4.T.cells",
    "Naive_CD4_T_cells_1" = "naive.CD4.T.cells",
    "Tregs" = "regulatory.T.cells",
    "Memory_CD4_T_cells_1" = "memory.CD4.T.cells",
    "Naive_CD8_T_cells_2" = "naive.CD8.T.cells",
    "Naive_CD4_T_cells_2" = "naive.CD4.T.cells")
  
  coarse.grained.translation <- list(
    "Naive_B_cells_1" = "B.cells",
    "Macrophages_2" = "monocytic.lineage",
    "Dendritic_cells_1" = "monocytic.lineage",
    "Macrophages_1" = "monocytic.lineage",
    "Fibroblasts" = "fibroblasts",
    "Endothelial_cells" = "endothelial.cells",
    "Memory_CD8_T_cells_1" = "CD8.T.cells",
    "Monocytes_1" = "monocytic.lineage",
    "Neutrophils_2" = "neutrophils",
    "Dendritic_cells_2" = "monocytic.lineage",
    "NK_cells_2" = "NK.cells",
    "Monocytes_2" = "monocytic.lineage",
    "NK_cells_1" = "NK.cells",
    "Memory_CD8_T_cells_2" = "CD8.T.cells",
    "Memory_CD4_T_cells_2" = "CD4.T.cells",
    "Naive_CD4_T_cells_1" = "CD4.T.cells",
    "Tregs" = "CD4.T.cells",
    "Memory_CD4_T_cells_1" = "CD4.T.cells",
    "Naive_CD8_T_cells_2" = "CD8.T.cells",
    "Naive_CD4_T_cells_2" = "CD4.T.cells")
  
  fine.grained.tbl <- data.frame(sample = names(fine.grained.translation),
                                 fine.grained.cell.type = unname(unlist(fine.grained.translation)))
  
  coarse.grained.tbl <- data.frame(sample = names(coarse.grained.translation),
                                   coarse.grained.cell.type = unname(unlist(coarse.grained.translation)))
  tbl <- merge(fine.grained.tbl, coarse.grained.tbl, by = "sample")
  tbl
}





## This is the ground truth that the Stanford core actually created
## NB: this function returns ground truth specified at the level of cell types
## (e.g., naive.B.cells) not sample names (e.g., Naive_B_cells_2)
get.actual.validation.ground.truth <- function() {
  
  ## Get the final admixture ratios
  ## Note these rows to not sum to 1, so need to renormalize
  ## This is the old/first "ratio" file that John Coller sent us "Final_Mixture_Ratios.xlsx"
  ## Instead use the subsequent file he sent "Final_Mixture_Fractions.xlsx" in which
  ## he appears to have just decided the rows by the row sum (which is what I do for the ratio file below).
  if(FALSE) {
    synId <- "syn21577248"
    obj <- synGet(synId, downloadFile = TRUE)
    stop("Need to read _both_ sheets -- this code is not doing that. See below for fractions file")
    ratios <- read.xlsx(obj$path, sheet = 1, startRow = 2)
    ratios <- ratios[-1,]
    nms <- ratios$X1
    ratios <- ratios[,-1]
    ratios <-
      ldply(1:nrow(ratios), .fun = function(i) ratios[i,] / sum(ratios[i,]))
    rownames(ratios) <- nms
  }
  
  synId <- "syn21598638"
  obj <- synGet(synId, downloadFile = TRUE)
  
  eps <- 10^-4
  
  old.col.names <- c("breast", "CRC", "Fibroblasts", "Endothelial_cells", "Dendritic_cells",
                     "Monocytes", "Macrophages", "NK_cells", "Tregs", "Naive_CD4_T_cells",
                     "Memory_CD4_T_cells", "Memory_CD8_T_cells", "Naive_B_cells")
  new.col.names <- c("BRCA", "CRC", "fibroblasts", "endothelial.cells", "DC", "monocytes",
                     "macrophages", "NK.cells", "regulatory.T.cells", "naive.CD4.T.cells",
                     "memory.CD4.T.cells", "memory.CD8.T.cells", "naive.B.cells")
  
  
  
  ratios <-
    ldply(1:2,
          .fun = function(sheetIndex) {
            df <- read.xlsx(obj$path, sheet = sheetIndex, startRow = 2)
            df <- df[-1,]
            nms <- df$X1
            df <- df[,-1]
            df <- df[, old.col.names]
            colnames(df) <- new.col.names
            df$sample <- nms
            df
          })
  
  rownames(ratios) <- ratios$sample
  ratios <- ratios[, !(colnames(ratios) == "sample")]
  for(col in 1:ncol(ratios)) { ratios[,col] <- as.numeric(ratios[,col]) }
  
  l_ply(1:nrow(ratios),
        .fun = function(i) if(abs(sum(as.numeric(ratios[i,])) - 1) > eps) { stop("Rows don't sum to one") })
  
  melted.ratios <- melt(as.matrix(ratios))
  colnames(melted.ratios) <- c("sample.id", "cell.type", "measured")
  melted.ratios$cell.type <- as.character(melted.ratios$cell.type)
  
  md <- get.validation.metadata()
  colnames(md) <- c("sample.id", "tumor.type", "batch", "mixture.type", "dataset.name")
  melted.ratios <- merge(melted.ratios, md)
  
  return(melted.ratios)
}

## This is the ground truth that the Stanford core actually created
## NB: this function returns ground truth specified at the level of
## sample names (e.g., Naive_B_cells_2) not cell types (e.g., naive.B.cells) 
get.actual.validation.ground.truth.sample.names <- function() {
  ## Get the final admixture ratios
  ## Note these rows to not sum to 1, so need to renormalize
  ## This is the old/first "ratio" file that John Coller sent us "Final_Mixture_Ratios.xlsx"
  ## Instead use the subsequent file he sent "Final_Mixture_Fractions.xlsx" in which
  ## he appears to have just decided the rows by the row sum (which is what I do for the ratio file below).
  if(FALSE) {
    synId <- "syn21577248"
    obj <- synGet(synId, downloadFile = TRUE)
    stop("Need to read _both_ sheets -- this code is not doing that. See below for fractions file")
    ratios <- read.xlsx(obj$path, sheet = 1, startRow = 2)
    ratios <- ratios[-1,]
    nms <- ratios$X1
    ratios <- ratios[,-1]
    ratios <-
      ldply(1:nrow(ratios), .fun = function(i) ratios[i,] / sum(ratios[i,]))
    rownames(ratios) <- nms
  }
  
  synId <- "syn21598638"
  obj <- synGet(synId, downloadFile = TRUE)
  old.col.names <- c("breast", "CRC", "Fibroblasts", "Endothelial_cells", "Dendritic_cells",
                     "Monocytes", "Macrophages", "NK_cells", "Tregs", "Naive_CD4_T_cells",
                     "Memory_CD4_T_cells", "Memory_CD8_T_cells", "Naive_B_cells")
  
  exclude.cols <- c("Tregs", "Endothelial_cells", "breast", "CRC", "Fibroblasts")
  
  ## Leave off the second sheet (note 1:1 below) since Naive_B_cells_2
  ## appear to have been used in admixtures, but not sequenced.
  sample.ratios <-
    ldply(1:1,
          .fun = function(sheetIndex) {
            df <- read.xlsx(obj$path, sheet = sheetIndex, startRow = 2)
            df <- df[-1,]
            nms <- df$X1
            df <- df[,-1]
            df <- df[, old.col.names]
            colnames(df) <-
              unlist(lapply(colnames(df),
                            function(str) ifelse(str %in% exclude.cols,
                                                 str,
                                                 paste0(str, "_", sheetIndex))))
            rownames(df) <- nms
            for(col in 1:ncol(df)) { df[,col] <- as.numeric(df[,col]) }		 
            m <- melt(as.matrix(df))
            colnames(m) <- c("sample", "cell.type", "actual")
            m
            
          })
  sample.ratios$cell.type <- as.character(sample.ratios$cell.type)
  sample.ratios$sample <- as.character(sample.ratios$sample)
  
  sample.ratios[sample.ratios$cell.type == "breast", "cell.type"] <- "Breast"
  eps <- 10^-4
  chk <- ddply(sample.ratios, .variables = c("sample"), .fun = function(df) sum(df$actual))
  if(any(abs(chk$V1 - 1) > eps)) { stop("Something doesn't sum to one") }
  sample.ratios
}

## Our original admixture specification includes the vendor for each sample
## Return a data frame with columns sample and vendor, given the vendor associated with each sample
get.vendor.sample.assignments <- function() {
  synId <- "syn21577258"
  obj <- synGet(synId, downloadFile = TRUE)
  vendors1 <- read.xlsx(obj$path, sheet = 1)[3,,drop=F]
  vendors2 <- read.xlsx(obj$path, sheet = 2)[3,,drop=F]
  vendors1 <- vendors1[,-1,drop=F]
  vendors2 <- vendors2[,-1,drop=F]
  exclude.cols <- c("Tregs", "Endothelial_cells", "breast", "CRC", "Fibroblasts")
  colnames(vendors1) <-
    unlist(lapply(colnames(vendors1),
                  function(str) ifelse(str %in% exclude.cols,
                                       str,
                                       paste0(str, "_1"))))
  colnames(vendors2) <-
    unlist(lapply(colnames(vendors2),
                  function(str) ifelse(str %in% exclude.cols,
                                       str,
                                       paste0(str, "_2"))))
  vendors <- cbind(vendors1, vendors2[, !(colnames(vendors2) %in% exclude.cols), drop=F])
  vendors <- t(rename.samples(vendors))
  rownames(vendors)[grepl(rownames(vendors), pattern="breast")] <- "Breast"
  colnames(vendors) <- c("vendor")
  vendors <- data.frame(sample = rownames(vendors), vendors)
  vendors
}

plot.admixtures <- function(mat) {
  
  hc <- hclust(dist(mat))
  labels <- hc$labels[hc$order]
  if("CRC" %in% labels) { labels <- c("CRC", labels[labels != "CRC"]) }
  if("breast" %in% labels) { labels <- c("breast", labels[labels != "breast"]) }
  if("Breast" %in% labels) { labels <- c("Breast", labels[labels != "Breast"]) }
  if("BRCA" %in% labels) { labels <- c("BRCA", labels[labels != "BRCA"]) }        
  
  mat <- mat[labels,]
  cols <- c("CRC", "breast", "Breast", "BRCA")
  cols <- cols[cols %in% labels]
  if(length(cols) == 2) {
    o <- order(mat[cols[1],,drop=T], mat[cols[2],,drop=T])
    mat <- mat[, o]
  } else {
    hc <- hclust(dist(t(mat)))
    mat <- mat[,hc$labels[hc$order]]
  }
  
  row.levels <- rownames(mat)
  col.levels <- colnames(mat)
  
  df <- melt(mat)
  df$Var1 <- factor(df$Var1, levels = row.levels)
  df$Var2 <- factor(df$Var2, levels = col.levels)
  
  g <- ggplot(data = df, aes(x = Var2, y = Var1, fill = value))
  g <- g + geom_tile()
  g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 text = element_text(size=15))
  g <- g + scale_fill_gradient2("Proportion", limits = c(0,1), low = "red", high = "blue", mid = "white", na.value = "black")
  g <- g + xlab("Admixture") + ylab("Cell Type")
  g <- g + theme(axis.text.x = element_blank())
  g
}

# whether or not to order methods in the heatmap using both a summary statistic (e.g., mean)
# that excludes nulls and one that does not or only to use that excludes or includes ...
use.include.and.exclude.null.in.summary <- FALSE
# ... and if only one, which to use
na.rm.summary <- TRUE


calculate.method.levels <- function(df, id.var = "modelId", cell.type.var = "cell.type", cor.var = "cor.p",
                                    row.summary.fun = "mean", col.summary.fun = "max",
                                    order.decreasing = FALSE) {
  
  orig.df <- df
  df <- df[, c(id.var, cell.type.var, cor.var)]
  df[, id.var] <- as.character(df[, id.var])
  df[, cell.type.var] <- as.character(df[, cell.type.var])
  
  # Calculate the method summary over cell types (e.g., mean) with and without excluding NAs (if use.include.and.exclude.null.in.summary == TRUE)
  # e.g., "mean" = mean without excluding NAs
  #       "nnmean" = non-null mean (mean excluding NAs)
  # i.e., always prepend "nn" to summary
  # Order NAs from highest to lowest based on nnmean and then non-NAs from highest to lowest based on mean.
  # Admittedly, this is a strange rule, but I like the end result in the heatmaps -- it makes a block of
  # the comparators at the top that have NAs, with the next score the top performing method (comparator or
  # or otherwise). This makes that top performer stand out more than it would if the blanks/NAs were 
  # scattered throughout.
  nnrow.summary.fun <- paste0("nn", row.summary.fun)
  method.summaries <- ddply(df, .variables = c(id.var),
                            .fun = function(tmp) {
                              if(use.include.and.exclude.null.in.summary) {
                                nret <- data.frame(cell.type = row.summary.fun, cor = do.call(row.summary.fun, list(tmp[, cor.var], na.rm=FALSE)))
                                nnret <- data.frame(cell.type = nnrow.summary.fun, cor = do.call(row.summary.fun, list(tmp[, cor.var], na.rm=TRUE)))
                                ret <- rbind(nret, nnret)
                              } else {
                                ret <- data.frame(cell.type = row.summary.fun, cor = do.call(row.summary.fun, list(tmp[, cor.var], na.rm=na.rm.summary)))
                              }
                              colnames(ret) <- c(cell.type.var, cor.var)
                              ret
                            })
  
  if(use.include.and.exclude.null.in.summary) {
    na.methods <- method.summaries[is.na(method.summaries[, cor.var]),id.var]
    na.summaries <- method.summaries[(method.summaries[, cell.type.var] == nnrow.summary.fun) & (method.summaries[, id.var] %in% na.methods),]
    na.summaries <- na.summaries[order(na.summaries[, cor.var], decreasing = order.decreasing),]
    flag1 <- !is.na(method.summaries[, cor.var])
    flag2 <- !(method.summaries[, id.var] %in% na.summaries[, id.var])
    flag3 <- (method.summaries[, cell.type.var] == row.summary.fun)
    non.na.summaries <- method.summaries[flag1 & flag2 & flag3, ]
    non.na.summaries <- non.na.summaries[order(non.na.summaries[, cor.var], decreasing = order.decreasing),]
    method.summaries <- rbind(non.na.summaries, na.summaries)
  } else {
    method.summaries <- method.summaries[order(method.summaries[, cor.var], decreasing = order.decreasing),]
  }
  method.levels <- method.summaries[method.summaries[, id.var] != col.summary.fun, id.var]
  
  method.levels
}

calculate.cell.type.levels <- function(df, id.var = "modelId", cell.type.var = "cell.type", cor.var = "cor.p",
                                       row.summary.fun = "mean", col.summary.fun = "max",
                                       order.decreasing = FALSE) {
  
  orig.df <- df
  df <- df[, c(id.var, cell.type.var, cor.var)]
  df[, id.var] <- as.character(df[, id.var])
  df[, cell.type.var] <- as.character(df[, cell.type.var])
  
  na.rm <- FALSE
  
  cell.type.summaries <- ddply(df, .variables = cell.type.var,
                               .fun = function(tmp) {
                                 ret <- data.frame(id = col.summary.fun, cor = do.call(col.summary.fun, list(tmp[, cor.var], na.rm=TRUE)))
                                 colnames(ret)[1] <- id.var
                                 colnames(ret)[2] <- cor.var                                 
                                 ret
                               })
  cell.type.summaries <- cell.type.summaries[order(cell.type.summaries[, cor.var], decreasing = order.decreasing),]
  
  cell.type.levels <- cell.type.summaries[cell.type.summaries[, cell.type.var] != row.summary.fun, cell.type.var]
  
  cell.type.levels
}


# See
# https://stackoverflow.com/questions/61733297/apply-bold-font-on-specific-axis-ticks
#library(ggtext)
#library(glue)
#bold.highlight <- function(x, strings.to.bolden) {
#  # ifelse(grepl(pat, x), glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
#  ifelse(x %in% strings.to.bolden, glue("<b>{x}</b>"), x)
#}

plot.cell.type.correlation.heatmap <- function(df, show.corr.text = FALSE, id.var = "modelId", cell.type.var = "cell.type", cor.var = "cor.p",
                                               cor.type.label = "Pearson\nCorrelation", limits = c(-1, 1),
                                               pval.var = NULL, row.summary.fun = "mean", col.summary.fun = "max",
                                               second.col.summary.fun = NULL,
                                               order.decreasing = FALSE,
                                               method.levels = NULL,
                                               cell.type.levels = NULL,
                                               highlight.fun = "max",
                                               ids.to.bold = NULL,
                                               formatter = function(x) formatC(x, format="f", digits=2)) {
  orig.df <- df
  df <- df[, c(id.var, cell.type.var, cor.var)]
  df[, id.var] <- as.character(df[, id.var])
  df[, cell.type.var] <- as.character(df[, cell.type.var])
  
  ## Round the correlation values (not just the text below), so that when we highlight
  ## the max values we will highlight any ties (up to rounding)
  df[, cor.var] <- unlist(lapply(df[, cor.var], function(x) as.numeric(formatter(x))))
  
  ## Add NAs for any missing entries
  df <- acast(df, as.formula(paste0(id.var, " ~ ", cell.type.var)), value.var = cor.var, fill = NA)
  df <- reshape2::melt(as.matrix(df))
  colnames(df) <- c(id.var, cell.type.var, cor.var)
  df[, id.var] <- as.character(df[, id.var])
  df[, cell.type.var] <- as.character(df[, cell.type.var])
  # highlight the max performer
  df <- ddply(df, .variables=c(cell.type.var), 
              .fun = function(sub.df) { 
                sub.df$highlight <- !is.na(sub.df[,cor.var]) & ( (sub.df[,cor.var] == do.call(highlight.fun, list(sub.df[,cor.var], na.rm=TRUE))))
                sub.df
              })
  
  # na.rm <- FALSE
  cell.type.summaries <- ddply(df, .variables = cell.type.var,
                               .fun = function(tmp) {
                                 ret <- data.frame(id = col.summary.fun, cor = do.call(col.summary.fun, list(tmp[, cor.var], na.rm=TRUE)), highlight = FALSE)
                                 if(!is.null(second.col.summary.fun)) {
                                   tmp <- data.frame(id = second.col.summary.fun, cor = do.call(second.col.summary.fun, list(tmp[, cor.var], na.rm=TRUE)), highlight = FALSE)
                                   ret <- rbind(ret, tmp)
                                 }
                                 colnames(ret)[1] <- id.var
                                 colnames(ret)[2] <- cor.var                                 
                                 ret
                               })
  
  cell.type.summaries <- cell.type.summaries[order(cell.type.summaries[, cor.var], decreasing = order.decreasing),]
  
  # Calculate the method summary over cell types (e.g., mean) with and without excluding NAs (if use.include.and.exclude.null.in.summary == TRUE)
  # e.g., "mean" = mean without excluding NAs
  #       "nnmean" = non-null mean (mean excluding NAs)
  # i.e., always prepend "nn" to summary
  nnrow.summary.fun <- paste0("nn", row.summary.fun)
  method.summaries <- ddply(df, .variables = c(id.var),
                            .fun = function(tmp) {
                              if(use.include.and.exclude.null.in.summary) {
                                ret <- data.frame(cell.type = row.summary.fun, cor = do.call(row.summary.fun, list(tmp[, cor.var], na.rm=FALSE)), highlight = FALSE)
                                nnret <- data.frame(cell.type = nnrow.summary.fun, cor = do.call(row.summary.fun, list(tmp[, cor.var], na.rm=TRUE)), highlight = FALSE)
                                ret <- rbind(ret, nnret)
                              } else {
                                ret <- data.frame(cell.type = row.summary.fun, cor = do.call(row.summary.fun, list(tmp[, cor.var], na.rm=na.rm.summary)), highlight = FALSE)
                              }
                              colnames(ret)[1:2] <- c(cell.type.var, cor.var)
                              ret
                            })
  
  method.summaries <- method.summaries[order(method.summaries[, cor.var], decreasing = order.decreasing),]
  
  
  if(is.null(cell.type.levels)) {
    # cell.type.levels <- c(grep(x = cell.type.summaries[, cell.type.var], pattern = row.summary.fun, values = TRUE), nnrow.summary.fun, row.summary.fun)
    cell.type.levels <- c(cell.type.summaries[!grepl(x = cell.type.summaries[, cell.type.var], pattern = row.summary.fun), cell.type.var], nnrow.summary.fun, row.summary.fun)
  } else {
    cell.type.levels <- c(cell.type.levels, nnrow.summary.fun, row.summary.fun)
  }
  
  col.summary.funs <- col.summary.fun
  if(!is.null(second.col.summary.fun)) { col.summary.funs <- c(col.summary.fun, second.col.summary.fun) }
  ## method.levels <- c(col.summary.funs, method.summaries[!(method.summaries[, id.var] %in% col.summary.funs), id.var])
  if(is.null(method.levels)) {
    method.levels <- c(col.summary.funs, method.summaries[!(method.summaries[, id.var] %in% col.summary.funs), id.var])
  } else {
    method.levels <- c(col.summary.funs, method.levels)
  }
  
  df <- rbind(df, cell.type.summaries[, c(id.var, cell.type.var, cor.var, "highlight")])
  df <- rbind(df, method.summaries[, c(id.var, cell.type.var, cor.var, "highlight")])    
  
  ## df$cor.label <- formatC(df[, cor.var], format="f", digits=digits)
  ## formatter <- function(x) formatC(x, format="f", digits=2)
  df$cor.label <- formatter(df[, cor.var])    
  if(!is.null(pval.var)) {
    p_load(gtools)
    df <- merge(df, orig.df[, c(id.var, cell.type.var, pval.var)], all.x = TRUE)
    flag <- is.na(df[, pval.var])
    df[flag, pval.var] <- 1
    df$cor.label <- paste0(stars.pval(df[, pval.var]), "\n", df$cor.label)
  }
  df[, cell.type.var] <- factor(df[, cell.type.var], levels = cell.type.levels)
  df[, id.var] <- factor(df[, id.var], levels = method.levels)
  
  g <- ggplot(data = df, aes_string(y = id.var, x = cell.type.var, fill = cor.var))
  g <- g + geom_tile()
  
  if(show.corr.text) {
    #        g <- g + geom_text(aes(label = cor.label))
    g <- g + geom_text(data = subset(df, highlight == TRUE), aes(label = cor.label), fontface = "bold.italic")
    g <- g + geom_text(data = subset(df, highlight != TRUE), aes(label = cor.label))
  }
  g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 text = element_text(size=15))
  if(!is.null(ids.to.bold)) {
    # g <- g + scale_y_discrete(labels = function(x) bold.highlight(x, ids.to.bold)) + theme(axis.text.y=element_markdown())
    bold.labels <- ifelse(levels(df[,id.var]) %in% ids.to.bold, yes = "bold", no = "plain")
    g <- g + theme(axis.text.y = element_text(face = bold.labels))
  }
  ## g <- g + ylab("Method") + xlab("")
  g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  ## g <- g + scale_fill_continuous("Pearson\ncorrelation", limits = c(-1,1))
  ## g <- g + scale_fill_gradient2("Pearson\ncorrelation", limits = c(-1,1),
  ##                               low = "red", high = "blue", mid = "white", na.value = "black")
  ##    g <- g + scale_fill_gradient2(paste0(cor.type.label, "\nCorrelation"),
  ##                                  limits = c(-1,1), low = "red", high = "blue", mid = "white", na.value = "black")
  g <- g + scale_fill_gradient2(cor.type.label, limits = limits,
                                low = "blue", high = "red", mid = "white", na.value = "black")
  ## g <- g + theme(text = element_text(size=20))
  g
}

proportion.labels <- function(x) ifelse(x == 0, "0", ifelse(x == 1, "1", x))

plot.strip.plots <- function(df, id.var = "modelId", cell.type.var = "cell.type", var = "cor.p",
                             label = "Pearson\nCorrelation", digits = 2, limits = c(-1, 1),
                             pval.var = NULL, col.summary.fun = "max", row.summary.fun = "mean", order.decreasing = FALSE,
                             method.levels = NULL,
                             cell.type.levels = NULL) {
  orig.df <- df
  df <- df[, c(id.var, cell.type.var, var)]
  df[, id.var] <- as.character(df[, id.var])
  df[, cell.type.var] <- as.character(df[, cell.type.var])
  
  na.rm <- TRUE
  if(is.null(cell.type.levels)) {
    cell.type.summaries <- ddply(df, .variables = cell.type.var,
                                 .fun = function(tmp) {
                                   ret <- data.frame(id = col.summary.fun, cor = do.call(col.summary.fun, list(tmp[, var], na.rm=TRUE)))
                                   colnames(ret)[1] <- id.var
                                   colnames(ret)[2] <- var                                 
                                   ret
                                 })
    cell.type.summaries <- cell.type.summaries[order(cell.type.summaries[, var], decreasing = order.decreasing),]
    
    ##    cell.type.levels <- c(cell.type.summaries[cell.type.summaries[, cell.type.var] != row.summary.fun, cell.type.var], row.summary.fun)
    cell.type.levels <- cell.type.summaries[cell.type.summaries[, cell.type.var] != row.summary.fun, cell.type.var]
  } 
  
  if(is.null(method.levels)) {
    method.summaries <- ddply(df, .variables = c(id.var),
                              .fun = function(tmp) {
                                ret <- data.frame(cell.type = row.summary.fun, cor = do.call(row.summary.fun, list(tmp[, var], na.rm=na.rm)))
                                colnames(ret) <- c(cell.type.var, var)
                                ret
                              })
    
    method.summaries <- method.summaries[order(method.summaries[, var], decreasing = order.decreasing),]
    
    method.levels <- method.summaries[method.summaries[, id.var] != col.summary.fun, id.var]
  }
  
  
  df$cor.label <- formatC(df[, var], format="f", digits=digits)
  if(!is.null(pval.var)) {
    p_load(gtools)
    df <- merge(df, orig.df[, c(id.var, cell.type.var, pval.var)], all.x = TRUE)
    flag <- is.na(df[, pval.var])
    df[flag, pval.var] <- 1
    df$cor.label <- paste0(stars.pval(df[, pval.var]), "\n", df$cor.label)
  }
  cat("ggplot'ing\n")
  ## This re-ordering o works with as.table = FALSE to respect the
  ## the ordering we want from cell.type.levels
  ntot <- length(cell.type.levels)
  nrow <- 3
  ncol <- ceiling(ntot / nrow)
  o <- unlist(llply(1:(nrow-1), .fun = function(i) (i*ncol):(1+((i-1)*ncol))))
  o <- c(o, ntot:(max(o)+1))
  df[, cell.type.var] <- factor(df[, cell.type.var], levels = cell.type.levels[o])
  df[, id.var] <- factor(df[, id.var], levels = method.levels)
  g <- ggplot(data = df, aes_string(y = id.var, x = var))
  g <- g + geom_boxplot(outlier.shape = NA)
  g <- g + facet_wrap(cell.type.var, as.table = FALSE, nrow = nrow)
  ##    g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), text = element_text(size=15))
  strip.text.sz <- 15
  if(ntot > 10) { strip.text.sz <- 8 }
  g <- g + theme(axis.text.y = element_text(size = 6), text = element_text(size=15), strip.text = element_text(size = strip.text.sz))
  
  g <- g + scale_x_continuous(labels = proportion.labels)
  g <- g + xlab(label) + ylab("")
  return(list("g" = g, "df" = df))
}

plot.correlation <- function(x, y, labels = NULL, colors = NULL, display.r2 = FALSE, method = "pearson", display.pval = FALSE, xoffset = 0.5, ...) {
  df <- data.frame(x = x, y = y)
  if(!is.null(labels)) {
    df$labels <- labels
  }
  g <- NULL
  if(is.null(labels)) {
    g <- ggplot(df, aes(x = x, y = y))
  } else {
    g <- ggplot(df, aes(x = x, y = y, label = labels))
  }
  if(!is.null(colors)) {
    g <- g + geom_point(aes(colour = colors))
  } else {
    g <- g + geom_point()
  }
  if(!is.null(labels)) {
    g <- g + geom_text(vjust = "inward", hjust = "inward")
    ##    suppressPackageStartupMessages(p_load(ggrepel))
    ##    g <- g + geom_text_repel(point.padding = NA, box.padding = 1)
  }
  ##  g <- g + theme(legend.position="none")
  g <- g + geom_smooth(data = df, aes(x = x, y = y), method='lm')
  x.min <- min(df$x, na.rm=TRUE)
  x.max <- max(df$x, na.rm=TRUE)
  y.min <- min(df$y, na.rm=TRUE)
  y.max <- max(df$y, na.rm=TRUE)
  
  ylimits <- NULL
  if(FALSE) {
    use.ggplot.2.2.1.limit.code <- TRUE
    if(use.ggplot.2.2.1.limit.code) {
      ylimits <- ggplot_build(g)$layout$panel_ranges[[1]]$y.range
      xlimits <- ggplot_build(g)$layout$panel_ranges[[1]]$x.range
    } else {
      ylimits <- ggplot_build(g)$layout$panel_params[[1]]$y.range
      xlimits <- ggplot_build(g)$layout$panel_params[[1]]$x.range
    }
  }
  xlimits <- ggplot_build(g)$layout$panel_params[[1]]$x.range
  ylimits <- ggplot_build(g)$layout$panel_params[[1]]$y.range
  
  ## to see why geom_text(size = sz) sz is different than in theme see: ratio of 14/5
  ## https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control/25062509 
  ##  g <- g + geom_text(x = x.min + 0.5 * (x.max - x.min), y = y.min + 1 * (y.max - y.min), label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
  ##  g <- g + geom_text(x = x.min + 0.5 * (x.max - x.min), y = y.min + 0.8 * (y.max - y.min), label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
  ##  g <- g + geom_text(x = x.min + 0.5 * (x.max - x.min), y = 0.8 * ylimits[2], label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
  sz <- 25
  g <- g + geom_text(x = xlimits[1] + xoffset * (xlimits[2] - xlimits[1]), y = ylimits[1] + 0.8 * (ylimits[2] - ylimits[1]), label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
  g <- g +theme(text = element_text(size = sz),
                axis.text.x = element_text(size=sz),
                axis.text.y = element_text(size=sz),
                axis.title.x = element_text(size=sz),
                axis.title.y = element_text(size=sz),
                title = element_text(size=sz),
                plot.title = element_text(hjust = 0.5, size=sz))
  g
}

limit.matrix.to.protein.coding <- function(mat, use.symbols = TRUE) {
  
  use.biomaRt <- TRUE
  exclude.mt <- FALSE
  # Need to do follow this
  # https://support.bioconductor.org/p/p132709/#p133562
  # to resolve an error in _filter
  #Install previous version of dplyr (0.8) 
  # devtools::install_github("hadley/dplyr@v0.8.0")
  
  #Install previous version of dbplyr (1.3)
  #devtools::install_url("https://cran.r-project.org/src/contrib/Archive/dbplyr/dbplyr_1.3.0.tar.gz")
  
  #Install AnnotationHub for Bioconductor 3.10 
  #BiocManager::install("AnnotationHub", version = "3.10", dependencies = TRUE)
  if(use.biomaRt) {
    suppressPackageStartupMessages(p_load(biomaRt))
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    pc.tbl <- getBM(attributes = c("ensembl_gene_id","external_gene_name","description","chromosome_name"),
                    filters = 'biotype', values = c('protein_coding'), mart = ensembl)
    if(exclude.mt) {
      pc.tbl <- subset(pc.tbl, chrosome_name != "MT")
    }
    keys <- as.character(pc.tbl$external_gene_name)
    if(!use.symbols) { keys <- as.character(pc.tbl$ensembl_gene_id) }
    mat <- mat[rownames(mat) %in% keys,]
    mat <- sweep(mat, 2, colSums(mat), `/`) * 10^6
    return(mat)
  }
  if(exclude.mt) { stop("Have not implemented exclude.mt for annotationHub -- but is easy\n") }
  suppressPackageStartupMessages(p_load(AnnotationHub))
  suppressPackageStartupMessages(p_load(ensembldb))
  ah <- AnnotationHub()
  flag <- (ah$species == "Homo sapiens") & (ah$genome == "GRCh38") & (ah$dataprovider == "Ensembl") & (ah$rdataclass == "EnsDb")
  ah2 <- ah[flag, ]
  ## as.data.frame(mcols(ah2))[1:10,c("title"),drop=FALSE]
  edb <- ah2[["AH73881"]]
  
  ## keytypes(edb)
  ## columns(edb)
  keys <- keys(edb, "GENENAME")
  columns <- c("GENEID", "ENTREZID", "GENEBIOTYPE")
  tbl <- ensembldb::select(edb, keys, columns, keytype = "GENENAME")
  pc.tbl <- subset(tbl, GENEBIOTYPE == "protein_coding")
  
  keys <- as.character(pc.tbl$GENENAME)
  if(!use.symbols) { keys <- as.character(pc.tbl$GENEID) }
  mat <- mat[rownames(mat) %in% keys,]
  mat <- sweep(mat, 2, colSums(mat), `/`) * 10^6
  mat
}


## Create in silico admixtures from the purified profiles
## sample.ratios is a data frame with columns sample.col (indicating the name of a sample
## admixture), cell.type.col (indicating a cell type within the sample), and fraction.col
## (indicating the fraction of the corresponding cell type within the sample). Note that
## the cell types indicated in cell.type.col should exist in mat, which provides the
## profiles that are mixed together.
create.in.silico.admixtures <- function(mat, sample.ratios,
                                        sample.col = "sample", cell.type.col = "cell.type",
                                        fraction.col = "actual") {
  insilico.admixtures <-
    dlply(sample.ratios,
          .variables = sample.col,
          .fun = function(df) {
            cols <- as.character(df[, cell.type.col])
            fracs <- as.numeric(df[, fraction.col])
            tmp.mat <- mat[, cols]
            ret <- as.matrix(tmp.mat) %*% fracs
            colnames(ret) <- df[1, sample.col]
            ret
          })
  insilico.admixtures <- do.call(cbind, insilico.admixtures)
  insilico.admixtures
}

# Compute predictions for a single method and subchallenge!
# i.e., preds are assumed to be for a specific method and subchallenge
score.preds <- function(preds) {
  # As described in Challenge here
  # Overall score is computed by:
  # 1. Computing correlation/RMSE independently for each dataset/cell type
  # 2. Averaging over cell type for each dataset
  # 3. Averaging over dataset
  res <- 
    ddply(preds, .variables = c("cell.type", "dataset.name"),
          .fun = function(df) {
            cor.p <- cor(df$prediction, df$measured, method="pearson")
            cor.s <- cor(df$prediction, df$measured, method="spearman")
            rmse <- sqrt(mean((df$prediction - df$measured)^2))
            data.frame(cor.p = cor.p, cor.s = cor.s, rmse = rmse)
          })
  
  means.across.cell.types <-
    ddply(res, .variables = c("dataset.name"),
          .fun = function(df) {
            data.frame(cor.p = mean(df$cor.p), cor.s = mean(df$cor.s), rmse = mean(df$rmse))
          })
  
  means.across.datasets <-
    ddply(res, .variables = c("cell.type"),
          .fun = function(df) {
            data.frame(cor.p = mean(df$cor.p), cor.s = mean(df$cor.s), rmse = mean(df$rmse))
          })
  
  overall <- 
    data.frame(cor.p = mean(means.across.cell.types$cor.p),
               cor.s = mean(means.across.cell.types$cor.s),
               rmse = mean(means.across.cell.types$rmse))
  
  list(all = res, means.across.datasets = means.across.datasets, overall = overall)
  
}

load.challenge.validation.data <- function(validation_path_base) { 
  # Location of Challenge validation admixtures on Synapse
  total_data_path <- paste0(validation_path_base,"/total_validation_data" )
  # List all files in the Challenge validation folder
  children <- list.files(total_data_path, full.names = F)
  l <- as.list(children)
  df <- do.call(rbind.data.frame, l)
  colnames(df) <-  c("name")
  # Load the "input.csv" file that lists the expression matrices in each
  # of the Challenge validation datasets
  input.tbl <- read.table(paste0(total_data_path, "/input.csv"), sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)
  # Subset to the expression data with Hugo symbols
  input.files <- input.tbl[, c("dataset.name", "hugo.expr.file")]
  df <- merge(df, input.files, by.x = "name", by.y = "hugo.expr.file")
  
  # Download each dataset from Synapse and read it in
  datasets <- 
    llply(df$name, 
          .fun = function(file_name) {
            obj <-read.csv(paste0(total_data_path, "/", file_name), row.names = 1)
            mat <- as.data.frame(obj)
          })
  names(datasets) <- df$dataset.name
  datasets
}

# this defines methods.to.exclude
methods.to.exclude <- c("TIMER", "timer")

# Return a list with named entries:
# predictions: the predictions from the comparator and participant methods, with columns
#              dataset.name, subchallenge, sample.id, cell.type, prediction, method.name, and submission (1, 2, 3, or latest)
# dataset.anno: the description of each dataset, with columns
#              dataset.name, tumor.type, distribution.type, mixture.type
# ground.truth: the ground truth, with columns:
#              dataset.name, subchallenge, sample.id, cell.type, measured
# In all cases, limit results to target.submissions and target.subchallenges
load.challenge.validation.results <- function(validation_path_base, target.submissions = c("1", "2", "3", "latest"), 
                                              target.subchallenges = c("fine", "coarse"),
                                              method.name.col = "method.name") {
  data_path_to_read <-  paste0(validation_path_base, "/rerun-validation-sanitized-predictions.csv")  #"syn22320329"
  obj <- read.csv(data_path_to_read)
  res.all <- as.data.frame(obj)
  
  flag <- res.all[,method.name.col] %in% methods.to.exclude
  cat(paste0("Excluding methods: ", paste(unique(res.all[flag, method.name.col]), collapse = ", "), "\n"))
  res.all <- res.all[!flag,]
  
  res.all <- subset(res.all, submission %in% target.submissions)
  res.all <- subset(res.all, subchallenge %in% target.subchallenges)
  
  preds <- res.all[, c("dataset.name", "subchallenge", "sample.id", "cell.type", "prediction", "method.name", "submission")]
  anno <- unique(res.all[, c("dataset.name", "tumor.type", "distribution.type", "mixture.type")])
  rownames(anno) <- NULL
  gt <- unique(res.all[, c("dataset.name", "subchallenge", "sample.id", "cell.type", "measured")])
  rownames(gt) <- NULL
  
  ret <- list("predictions" = preds, "dataset.anno" = anno, "ground.truth" = gt)
  ret
}

## Ensure we have a prediction (even if it is NA) for all cell types in all datasets by all methods
make.missing.predictions.na <- function(res.all, cell.type.col = "cell.type", prediction.col = "prediction",
                                        measured.col = "measured", subchallenge.col = "subchallenge",
                                        sample.id.col = "sample.id", dataset.name.col = "dataset.name",
                                        round.col = "submission") {
  tmp <- unique(res.all[, !(colnames(res.all) %in% c(cell.type.col, prediction.col, measured.col))])
  cell.types.by.sub.challenge <- unique(res.all[, c(cell.type.col, subchallenge.col, sample.id.col, dataset.name.col)])
  tmp <- merge(tmp, cell.types.by.sub.challenge, all = TRUE)
  measured.tbl <- unique(res.all[, c(cell.type.col, subchallenge.col, sample.id.col, dataset.name.col, measured.col)])
  prediction.tbl <- res.all[, c(cell.type.col, subchallenge.col, sample.id.col, dataset.name.col, method.name.col, round.col, prediction.col)]
  tmp <- merge(tmp, measured.tbl, all = TRUE)
  tmp <- merge(tmp, prediction.tbl, all = TRUE)
  res.all <- tmp
  flag <- !is.na(res.all[, measured.col])
  res.all <- res.all[flag, ]
  
  res.all[, cell.type.col] <- as.character(res.all[, cell.type.col])
  res.all
}

# Aggregate predictions from raw cell types as output by a deconvolution method
# into those cell types expected of the Challenge.
# Aggregation is performed via summing.
# The translation between raw and Challenge subtypes is encoded in the translation table,
# where the raw.cell.type column gives the raw cell type name and cell.type gives the Challenge cell type.
# Results in res_df are described by columns method.name, dataset.name, sample.id, cell.type, and prediction
aggregate.cell.type.predictions <- function(res_df, translation_df) {
  nm <- unique(res_df$method.name)
  res_df <- res_df %>%
    dplyr::inner_join(translation_df) %>%
    dplyr::select(dataset.name, sample.id, cell.type, prediction) %>%
    dplyr::group_by(dataset.name, sample.id, cell.type) %>%
    dplyr::summarise(prediction = sum(prediction)) %>%
    dplyr::ungroup()
  res_df <- as.data.frame(res_df)
  res_df <- cbind(method.name = nm, res_df)
}

sub.challenges <- list("coarse" = "coarse", "fine" = "fine")

calculate.empirical.bayes <-
  function(df, col.id, numerator.id, denominator.id, sample.id.cols, score.col) {
    flag <- df[, col.id] %in% c(numerator.id, denominator.id)
    tmp <- df[flag, ]
    n.tot <- nrow(tmp) / 2
    n.num <-
      sum(unlist(
        dlply(tmp, .variables = sample.id.cols,
              .fun = function(df.2) {
                if(nrow(df.2) != 2) {
                  print(df.2)
                  print(c(numerator.id, denominator.id))
                  stop("Was expecting 2 rows\n")
                }
                if(any(is.na(df.2[, score.col]))) { stop("Got NA scores\n") }
                rownames(df.2) <- df.2[, col.id]
                diff <- df.2[numerator.id, score.col] -
                  df.2[denominator.id, score.col]
                if(diff > 0) { return(1) }
                return(0)
              })))
    n.num / (n.tot - n.num)
  }

do.bootstrap.analysis <-
  function(res.input, bootstraps, 
           method.name.col, 
           subchallenge.col, measured.col, cell.type.col,
           dataset.name.col, sample.id.col, prediction.col,
           round.col, round = "latest") {
    
    submitter.tbl <- unique(res.input[, c(method.name.col, round.col, subchallenge.col), drop = FALSE])
    or <- order(submitter.tbl[, round.col])
    submitter.tbl <- submitter.tbl[or, ]
    submitter.tbl[, round.col] <- as.character(submitter.tbl[, round.col])
    flag <- submitter.tbl[, round.col] %in% c("latest", as.character(round))
    submitter.tbl <- submitter.tbl[flag, ]
    flag <- !duplicated(submitter.tbl[, c(method.name.col, subchallenge.col)], fromLast = FALSE)
    submitter.tbl <- submitter.tbl[flag, ]
    
    res.round <- merge(res.input, submitter.tbl, by = c(method.name.col, subchallenge.col, round.col))
    
    method.id.col <- method.name.col
    
    if(FALSE) {
      # Add ensemble here -- no need. We report bootstrapped scores
      # > colnames(results[["3"]][["res.round"]])
      # [1] "method.name"       "subchallenge"      "submission"       
      # [4] "dataset.name"      "sample.id"         "cell.type"        
      # [7] "objectId"          "comparator"        "submitterId"      
      # [10] "repo_name"         "tumor.type"        "distribution.type"
      # [13] "mixture.type"      "measured"          "prediction"       
      res.ensemble <-
        ddply(res.round[, !(colnames(res.round) %in% c("objectId", "comparator", "submittedId", "repo_name"))],
              .variables = c(round.col, subchallenge.col, dataset.name.col),
              .fun = function(df.in) {
                tres <- 
                  ddply(df.in, .variables = method.name.col,
                        .fun = function(df) {
                          df.ret <- data.frame(prediction = df[, prediction.col],
                                               pred.rank = rank(df[, prediction.col], ties.method="first"),
                                               sample.id = df[, sample.id.col])
                          colnames(df.ret) <- c(prediction.col, "pred.rank", sample.id.col)
                          df.ret
                        })
                # take the consensus (here, just mean) rank across methods
                cons <- 
                  ddply(tres, .variables = sample.id.col,
                        .fun = function(df) {
                          df.ret <- data.frame(cons.rank = mean(df$pred.rank))
                          colnames(df.ret)[1] <- prediction.col
                          df.ret
                        })
                tmp <- unique(df.in[, c(sample.id.col, measured.col, subchallenge.col, round.col, dataset.name.col, cell.type.col, "tumor.type", "distribution.type", "mixture.type")])
                cons <- merge(cons, tmp)
                cons$method.name <- "ensemble"
                cons$objectId <- "dummy"
                cons$comparator <- FALSE
                cons$submitterId <- "bwhite"
                cons$repo_name <- "dummy"
                ret <- cons[, c(method.name.col, subchallenge.col, round.col, dataset.name.col, sample.id.col, cell.type.col, 
                                "objectId", "comparator", "submitterId", "repo_name", "tumor.type", "distribution.type", "mixture.type",
                                measured.col, prediction.col)]
                ret
              })
      print(head(res.round))
      print(head(res.ensmble))
      res.round <- rbind(res.round, res.ensemble)
    } #end if(FALSE)
    res <- res.round
    
    tbls <-
      llply(sub.challenges,
            .fun = function(sub.challenge) {
              flag <- res[, subchallenge.col] == sub.challenge
              tmp <- res[flag, ]
              tmp
            })
    
    tbls.by.cell <-
      llply(sub.challenges,
            .fun = function(sub.challenge) {
              tmp <- tbls[[sub.challenge]]
              flag <- !is.na(tmp[, measured.col])
              tmp <- tmp[flag, ]
              cells <- unique(tmp[, cell.type.col])
              names(cells) <- cells
              llply(cells,
                    .fun = function(ct) {
                      flag <- tmp[, cell.type.col] == ct
                      ret <- tmp[flag, ]
                      ret$id <- paste0(ret[, dataset.name.col], "-", ret[, sample.id.col])
                      ret
                    })
            })
    
    ## Calculate both pearson and spearman correlation over bootstraps, within dataset and cell type
    cat(paste0("Calculating bootstrap correlations (ensemble)\n"))
    bootstrapped.cors.ensemble <-
      llply(sub.challenges, 
            .parallel = FALSE,
            .fun = function(sub.challenge) {
              n.bootstraps <- length(bootstraps[[sub.challenge]])
              indices <- 1:n.bootstraps
              # indices <- 1:10
              names(indices) <- indices
              
              ret.i <-
                ldply(indices,
                      .parallel = TRUE,
                      .fun = function(i) {
                        ret.cons <-
                          ddply(tbls[[sub.challenge]], .variables = c(cell.type.col, dataset.name.col),
                                .fun = function(df.in) {
                                  # method.name is row
                                  x <- acast(df.in[, c(prediction.col, method.name.col, sample.id.col)], 
                                             formula = paste0(method.name.col, " ~ ", sample.id.col), value.var = prediction.col)
                                  
                                  colnames(x) <- paste0(df.in[1, dataset.name.col], "-", colnames(x))
                                  sample.ids <- bootstraps[[sub.challenge]][[i]]$id
                                  sample.ids <- sample.ids[sample.ids %in% colnames(x)]
                                  x <- x[, sample.ids]
                                  m <- melt(x)
                                  colnames(m) <- c(method.name.col, sample.id.col, prediction.col)
                                  tres <- 
                                    ddply(m, .variables = method.name.col,
                                          .fun = function(df) {
                                            df.ret <- data.frame(prediction = df[, prediction.col],
                                                                 pred.rank = rank(df[, prediction.col], ties.method="first"),
                                                                 sample.id = df[, sample.id.col])
                                            colnames(df.ret) <- c(prediction.col, "pred.rank", sample.id.col)
                                            df.ret
                                          })
                                  # take the consensus (here, just mean) rank across methods
                                  cons <- 
                                    ddply(tres, .variables = sample.id.col,
                                          .fun = function(df) {
                                            df.ret <- data.frame(cons.rank = mean(df$pred.rank))
                                            df.ret
                                          })
                                  tmp <- unique(df.in[, c(sample.id.col, measured.col)])
                                  tmp[, sample.id.col] <- paste0(df.in[1, dataset.name.col], "-", tmp[, sample.id.col])
                                  cons <- merge(cons, tmp)
                                  ret <- data.frame(method.name = "ensemble", pearson = NA, rmse = NA, spearman = cor(cons$cons.rank, cons[, measured.col]), pearson.fc = NA)
                                  colnames(ret) <- c(method.name.col, "pearson", "rmse", "spearman", "pearson.fc")
                                  ret
                                })
                      })
              # "method.name"  "boot.i"       "dataset.name" "cell.type"    "pearson"   "spearman"     "rmse" "pearson.fc"
              colnames(ret.i)[1] <- "boot.i"
              ret.i <- ret.i[, c(method.name.col, "boot.i", dataset.name.col, cell.type.col, "pearson", "spearman", "rmse", "pearson.fc")]
              ret.i 
            })
    cat(paste0("Calculating bootstrap correlations\n"))
    bootstrapped.cors <-
      llply(sub.challenges, 
            .fun = function(sub.challenge) {
              tmp <- tbls[[sub.challenge]]                      
              flag <- !is.na(tmp[, measured.col])
              tmp <- tmp[flag, ]
              methods <- unique(tmp[, method.id.col])
              ##                      methods <- c("Aginome-XMU", "CIBERSORTx")
              names(methods) <- methods
              n.bootstraps <- length(bootstraps[[sub.challenge]])
              indices <- 1:n.bootstraps
              names(indices) <- indices
              ret.all <-
                ldply(methods, .parallel = TRUE,
                      .fun = function(method) {
                        print(method)
                        ret.i <-
                          ldply(indices,
                                .fun = function(i) {
                                  ret.method <-
                                    ldply(tbls.by.cell[[sub.challenge]],
                                          .fun = function(df.in) {
                                            flag <- df.in[, method.id.col] == method
                                            df <- df.in[flag, ]
                                            if(any(duplicated(df$id))) {
                                              print(head(df[my.dup(df$id),]))
                                              stop("stop")
                                            }
                                            rownames(df) <- df$id
                                            sample.ids <- bootstraps[[sub.challenge]][[i]]$id
                                            sample.ids <- sample.ids[sample.ids %in% df$id]
                                            if(!(all(sample.ids %in% rownames(df)))) {
                                              stop("Some sample ids not in df\n")
                                            }
                                            df <- df[sample.ids,]
                                            df
                                          })
                                  ret.method <- ret.method[, -1]
                                  score <-
                                    ddply(ret.method, .variables = c(dataset.name.col),
                                          .fun = function(df.ds) {
                                            tmp <- ddply(df.ds, .variables = c(cell.type.col),
                                                         .fun = function(df.ct) {
                                                           if(any(is.na(df.ct[,2]))) {
                                                             print(df.ct)
                                                           }
                                                           mts <- c("pearson", "spearman", "rmse", "pearson.fc")
                                                           names(mts) <- mts
                                                           vals <- llply(mts,
                                                                         .fun = function(comp.method) {
                                                                           pred <- as.numeric(df.ct[, prediction.col])
                                                                           measured <- as.numeric(df.ct[, measured.col])
                                                                           
                                                                           if(any(is.na(pred))) {
                                                                             if(!all(is.na(pred))) {
                                                                               stop("Some but not all NA\n")
                                                                             }
                                                                             return(NA)
                                                                           }
                                                                           if(var(pred) == 0) { return(0) }
                                                                           
                                                                           val <- NA
                                                                           if(comp.method %in% c("pearson", "spearman")) {
                                                                             val <- cor(pred, measured,
                                                                                        method = comp.method)
                                                                           } else if(comp.method == "pearson.fc") {
                                                                             # calculate the pearson correlation of log fold changes (across successive samples)
                                                                             # 1. Order the samples based on the their ground truth values g1 < g2 < … < gn (i.e., measured)
                                                                             # 2. Compute ground truth fold changes g2/g1, g3/g2, …, gn/gn-1
                                                                             # 3. Compute prediction fold changes p2/p1, p3/p2, …, pn/pn-1
                                                                             # 4. Compute the pearson correlation of these two vectors
                                                                             o <- order(measured, decreasing=FALSE)
                                                                             eps <- 10^-5
                                                                             gt.ordered <- measured[o] + eps
                                                                             # NB: the ordering is established by the _ground truth_ values and applied to those
                                                                             # values and the predicted values
                                                                             pred.ordered <- pred[o] + eps
                                                                             gt.fc <- unlist(lapply(2:length(gt.ordered), function(indx) gt.ordered[indx]/gt.ordered[indx-1]))
                                                                             pred.fc <- unlist(lapply(2:length(pred.ordered), function(indx) pred.ordered[indx]/pred.ordered[indx-1]))
                                                                             val <- cor(pred.fc, gt.fc, method = "pearson")
                                                                           } else if(comp.method == "rmse") {
                                                                             val <- sqrt(mean((pred-measured)^2))
                                                                           } else {
                                                                             stop(paste0("Unknown method ", comp.method, "\n"))
                                                                           }
                                                                           val
                                                                         })
                                                           as.data.frame(vals)
                                                         })
                                            colnames(tmp)[1] <- cell.type.col
                                            tmp
                                          })
                                  score
                                })
                        colnames(ret.i)[1] <- "boot.i"
                        ret.i
                      })
              colnames(ret.all)[1] <- method.id.col
              ret.all
            })
    
    for(nm in names(bootstrapped.cors)) {
      bootstrapped.cors[[nm]] <- rbind(bootstrapped.cors[[nm]], bootstrapped.cors.ensemble[[nm]])
    }
    
    print(colnames(bootstrapped.cors[[1]]))
    # "method.name"  "boot.i"       "dataset.name" "cell.type"    "pearson"   "spearman"     "rmse" "pearson.fc"
    print(head(bootstrapped.cors[[1]]))
    
    cat(paste0("Calculating bootstrapped scores\n"))
    bootstrapped.scores <-
      llply(sub.challenges, .parallel = TRUE,
            .fun = function(sub.challenge) {
              df <- bootstrapped.cors[[sub.challenge]]
              flag <- res.round[, subchallenge.col] == sub.challenge
              un <- unique(res.round[flag, unique(c(method.id.col, method.name.col))])
              df <- merge(df, un)
              ## Average over cell type (within method and dataset and bootstrap sample)
              df <- ddply(df,
                          .variables = c(method.id.col, dataset.name.col, "boot.i"),
                          .fun = function(tmp) {
                            data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman), rmse = mean(tmp$rmse), pearson.fc = mean(tmp$pearson.fc))
                          })
              ## Average over dataset (within method and bootstrap sample)
              df <- ddply(df,
                          .variables = c(method.id.col, "boot.i"),
                          .fun = function(tmp) {
                            data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman), rmse = mean(tmp$rmse), pearson.fc = mean(tmp$pearson.fc))
                          })
              df
            })
    
    cat(paste0("Calculating mean bootstrapped scores\n"))
    mean.bootstrapped.scores <-
      llply(bootstrapped.scores, .parallel = TRUE,
            .fun = function(df) {
              ## Average over bootstraps
              df <- ddply(df,
                          .variables = c(method.id.col),
                          .fun = function(tmp) {
                            data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman), rmse = mean(tmp$rmse), pearson.fc = mean(tmp$pearson.fc))
                          })
              o <- order(df$pearson, decreasing = TRUE)
              df <- df[o,]
              df
            })
    
    means.over.dataset <-
      llply(bootstrapped.cors,
            .fun = function(df) {
              methods <- c("pearson", "spearman", "rmse", "pearson.fc")
              na.rm <- FALSE
              names(methods) <- methods
              res <- llply(methods,
                           .fun = function(method) {
                             ## average over dataset
                             ret <- ddply(df, .variables = c(method.id.col, cell.type.col, "boot.i"),
                                          .fun = function(df) {
                                            data.frame(cor = mean(df[, method], na.rm=na.rm))
                                          })
                           })
            })
    
    means.over.bootstrap <-
      llply(bootstrapped.cors,
            .fun = function(df) {
              methods <- c("pearson", "spearman", "rmse", "pearson.fc")
              na.rm <- FALSE
              names(methods) <- methods
              res <- llply(methods,
                           .fun = function(method) {
                             ## average over dataset
                             ret <- ddply(df, .variables = c(method.id.col, cell.type.col, dataset.name.col),
                                          .fun = function(df) {
                                            data.frame(cor = mean(df[, method], na.rm=na.rm))
                                          })
                           })
            })
    
    cat(paste0("Calculating mean by cell type\n"))
    if(FALSE) {
      means.by.cell.type.method <-
        llply(bootstrapped.cors,
              .fun = function(df) {
                methods <- c("pearson", "spearman", "rmse", "pearson.fc")
                na.rm <- FALSE
                names(methods) <- methods
                res <- llply(methods,
                             .fun = function(method) {
                               ## first, average over bootstrap
                               ret <- ddply(df, .variables = c(method.id.col, cell.type.col, dataset.name.col),
                                            .fun = function(df) {
                                              data.frame(cor = mean(df[, method], na.rm=na.rm))
                                            })
                               ## now, average over dataset
                               ret2 <- ddply(ret, .variables = c(method.id.col, cell.type.col),
                                             .fun = function(df) {
                                               data.frame(cor = mean(df$cor, na.rm=na.rm))
                                             })
                             })
              })
    }
    
    means.by.cell.type.method <-
      llply(means.over.dataset,
            .fun = function(df) {
              methods <- c("pearson", "spearman", "rmse", "pearson.fc")
              na.rm <- FALSE
              names(methods) <- methods
              res <- llply(methods,
                           .fun = function(method) {
                             ## average over bootstrap (means.over.dataset has already been averaged over dataset)
                             ret <- ddply(df[[method]], .variables = c(method.id.col, cell.type.col),
                                          .fun = function(df) {
                                            data.frame(cor = mean(df$cor, na.rm=na.rm))
                                          })
                           })
            })
    
    cat(paste0("Calculating bayes factors\n"))
    top.performers <-
      llply(sub.challenges,
            .fun = function(sub.challenge) {
              mean.scores <- mean.bootstrapped.scores[[sub.challenge]]
              mean.scores <- na.omit(mean.scores)
              scores <- bootstrapped.scores[[sub.challenge]]
              numerator.indx <- 1
              numerator.id <- mean.scores[numerator.indx, method.id.col]
              numerator.id
            })
    
    bayes.factors <-
      llply(sub.challenges,
            .fun = function(sub.challenge) {
              mean.scores <- mean.bootstrapped.scores[[sub.challenge]]
              mean.scores <- na.omit(mean.scores)
              scores <- bootstrapped.scores[[sub.challenge]]
              numerator.indx <- 1
              numerator.id <- mean.scores[numerator.indx, method.id.col]
              indices <- 1:nrow(mean.scores)
              indices <- indices[!(indices == numerator.indx)]
              names(indices) <- mean.scores[indices, method.id.col]
              ret <-
                ldply(indices,
                      .fun = function(i) {
                        denominator.id <- mean.scores[i, method.id.col]
                        methods <- c("pearson", "spearman")
                        names(methods) <- methods
                        res <- llply(methods,
                                     .fun = function(method) {
                                       bf <- calculate.empirical.bayes(scores, col.id = method.id.col,
                                                                       numerator.id = numerator.id,
                                                                       denominator.id = denominator.id,
                                                                       sample.id.cols = "boot.i",
                                                                       score.col = method)
                                     })
                        ret.df <- as.data.frame(res)
                        ret.df
                      })
              ret
            })
    
    ret.list <- list("res.round" = res.round,
                     "bootstrapped.cors" = bootstrapped.cors,
                     "bootstrapped.scores" = bootstrapped.scores,
                     "mean.bootstrapped.scores" = mean.bootstrapped.scores,
                     "means.by.cell.type.method" = means.by.cell.type.method,
                     "means.over.dataset" = means.over.dataset,
                     "means.over.bootstrap" = means.over.bootstrap,                         
                     "top.performers" = top.performers,
                     "bayes.factors" = bayes.factors)
    
    return(ret.list)
    
  }

## Modified slightly from
## https://stackoverflow.com/questions/53170465/how-to-make-a-base-r-style-boxplot-using-ggplot2
geom_boxplotMod <- function(mapping = NULL, data = NULL, stat = "boxplot", 
                            position = "dodge2", ..., outlier.colour = NULL, outlier.color = NULL, 
                            outlier.fill = NULL, outlier.shape = 1, outlier.size = 1.5, 
                            outlier.stroke = 0.5, outlier.alpha = NULL, notch = FALSE, notchwidth = 0.5,
                            varwidth = FALSE, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,
                            linetype = "dashed") # to know how these come here use: args(geom_boxplot)
{
  list(geom_boxplot(
    mapping = mapping, data = data, stat = stat, position = position,
    outlier.colour = outlier.colour, outlier.color = outlier.color, 
    outlier.fill = outlier.fill, outlier.shape = outlier.shape, 
    outlier.size = outlier.size, outlier.stroke = outlier.stroke, 
    outlier.alpha = outlier.alpha, notch = notch, 
    notchwidth = notchwidth, varwidth = varwidth, na.rm = na.rm, 
    show.legend = show.legend, inherit.aes = inherit.aes, linetype = 
      linetype, ...),
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.25),
    #the width of the error-bar heads are decreased
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.25),
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), ...)
  )
}


make.round.text <- function(round, round.str = "Round") {
  round.text <- ""
  if(round == "latest") {
    round.text <- paste0("Latest ", round.str)
  } else if (round == "1") {
    round.text <- paste0(round.str, " 1")
  } else {
    round.text <- paste0("Latest ", round.str, " up to ", round.str, " ", round)
  }
  round.text
}


plot.bootstrap.analysis <-
  function(res, bootstrapped.scores, mean.boostrapped.scores, median.bootstrapped.scores,
           means.by.cell.type.method,
           means.over.dataset, method.anno.round,
           postfix, plot.spearman.distribution = FALSE) {
    
    top.performers <- c("Aginome-XMU", "DA_505", "mitten_TDC19", "Biogem")
    comparator.methods <- get.comparators.cap()
    print(comparator.methods) 
    priority.methods <- unique(c(top.performers, comparator.methods))
    
    for(sub.challenge in sub.challenges) {
      print(sub.challenge)
      print(head(method.anno.round))
      print(method.anno.round)
      print(subchallenge.col)
      flag <- is.na(method.anno.round[, subchallenge.col]) | ( method.anno.round[, subchallenge.col] == sub.challenge )
      print(flag)
      method.anno.round.sc <- method.anno.round[flag, c(method.name.col, "Output", "Method")]
      print(method.anno.round.sc)
      bootstrapped.scores[[sub.challenge]] <-
        merge(bootstrapped.scores[[sub.challenge]], method.anno.round.sc, all.x = TRUE)
      median.bootstrapped.scores[[sub.challenge]] <-
        merge(median.bootstrapped.scores[[sub.challenge]], method.anno.round.sc, all.x = TRUE)
      mean.bootstrapped.scores[[sub.challenge]] <-
        merge(mean.bootstrapped.scores[[sub.challenge]], method.anno.round.sc, all.x = TRUE)
    }
    
    barplots <- list()
    for(sub.challenge in sub.challenges) {
      scores <- bootstrapped.scores[[sub.challenge]]
      median.scores <- median.bootstrapped.scores[[sub.challenge]]
      
      median.scores[, method.name.col] <- as.character(median.scores[, method.name.col])
      flag <- median.scores[, method.name.col] == "ensemble"
      median.scores[flag, method.name.col] <- "consensus rank"
      
      scores[, method.name.col] <- as.character(scores[, method.name.col])
      flag <- scores[, method.name.col] == "ensemble"
      scores[flag, method.name.col] <- "consensus rank"
      
      
      o <- order(median.scores$pearson)
      median.scores <- median.scores[o, ]
      flag <- ( !is.na(scores$pearson) & !is.na(scores$spearman) ) | ( scores[,method.name.col] == "consensus rank")
      scores <- scores[flag, ] 
      scores[, method.name.col] <- factor(scores[, method.name.col], levels = median.scores[, method.name.col])
      
      lvls <- levels(scores[,method.name.col]) 
      lvls <- lvls[lvls %in% scores[, method.name.col]]
      bold.labels <- ifelse(lvls %in% comparator.methods, yes = "bold", no = "plain")
      cat("scores[, method.name.col]\n")
      print(lvls)
      print(bold.labels)
      print(table(bold.labels))
      #            scores <- na.omit(scores)
      flag <- ( !is.na(median.scores$pearson) & !is.na(median.scores$spearman) ) | ( median.scores[,method.name.col] == "consensus rank") | ( median.scores[,method.name.col] == "ensemble")
      median.scores[, method.name.col] <- factor(median.scores[, method.name.col], levels = median.scores[, method.name.col])
      #            median.scores <- na.omit(median.scores)
      median.scores <- median.scores[flag, ] 
      
      
      g1 <- ggplot(data = scores, aes_string(x = method.name.col, y = "pearson"))
      g1 <- g1 + geom_boxplotMod(fill = "#56B4E9")
      g1 <- g1 + coord_flip()
      g1 <- g1 + xlab("Method")
      ## g1 <- g1 + ylab("Pearson Correlation")
      g1 <- g1 + ylab("Pearson")
      g1 <- g1 + ylim(c(-0.25, 1))
      g1 <- g1 + theme(text = element_text(size=18), title = element_text(size = 20), axis.text.y = element_text(face = bold.labels))
      # g1 <- g1 + theme(text = element_text(size=18), title = element_text(size = 20))
      # g1 <- g1 + scale_y_discrete(labels = function(x) bold.highlight(x, comparator.methods))
      # g1 <- g1 + theme(axis.text.y=element_markdown())
      
      if(plot.spearman.distribution) {
        g2 <- ggplot(data = scores, aes_string(x = method.name.col, y = "spearman"))
        g2 <- g2 + geom_boxplotMod(fill = "#E69F00")
      } else{
        g2 <- ggplot(data = median.scores)
        g2 <- g2 + geom_col(aes_string(x = method.name.col, y = "spearman"), fill = "#E69F00")
      }
      g2 <- g2 + coord_flip()
      g2 <- g2 + xlab("Method")
      ## g2 <- g2 + ylab("Spearman Correlation")
      g2 <- g2 + ylab("Spearman")
      g2 <- g2 + ylim(c(-0.25, 1))            
      g2 <- g2 + theme(text = element_text(size=18))    
      g2 <- g2 + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                       axis.ticks.y = element_blank())
      
      g3 <- ggplot(data = scores, aes_string(x = method.name.col, y = "pearson.fc"))
      g3 <- g3 + geom_boxplotMod(fill = "#56B4E9")
      g3 <- g3 + coord_flip()
      g3 <- g3 + xlab("Method")
      ## g3 <- g3 + ylab("Pearson Correlation")
      g3 <- g3 + ylab("Pearson (Fold Change)")
      g3 <- g3 + ylim(c(-0.25, 1))
      g3 <- g3 + theme(text = element_text(size=18), title = element_text(size = 20), axis.text.y = element_text(face = bold.labels))
      # g3 <- g3 + theme(text = element_text(size=18), title = element_text(size = 20))
      # g3 <- g3 + scale_y_discrete(labels = function(x) bold.highlight(x, comparator.methods))
      # g3 <- g3 + theme(axis.text.y=element_markdown())
      
      
      tmp <- scores[, c(method.name.col, "Output", "Method")]            
      ret <- plot.anno.heatmap.with.multiple.legends(tmp, "method.name", c("Method", "Output"), c("Set3", "Set1"))
      
      full.plot <- ret[["full.plot"]]
      for.first.legend <- ret[["legends"]][["Method"]]
      for.second.legend <- ret[["legends"]][["Output"]]
      
      legs <- plot_grid(get_legend(for.first.legend), get_legend(for.second.legend), nrow = 2, align = "v", rel_heights = c(2,1))
      
      ## pg <- plot_grid(g1, g2, full.plot, legs, nrow=1, align="h", rel_widths = c(3,1,0.5,0.5))
      
      barplots[[paste0(sub.challenge,  "-pearson")]] <- g1
      barplots[[paste0(sub.challenge,  "-spearman")]] <- g2            
      barplots[[paste0(sub.challenge,  "-pearson.fc")]] <- g3            
      barplots[[paste0(sub.challenge,  "-anno")]] <- full.plot
      barplots[[paste0(sub.challenge,  "-legend")]] <- legs
      
      title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge")
      round.text <- make.round.text(round)            
      title <- paste0(title, " (", round.text, ")")
      ## png(paste0(figs.dir, "rerun-validation-score-box-and-barplots-", sub.challenge, postfix, ".png"), width = 2 * 480)
      ## grid.arrange(g1, g2, nrow=1, widths = c(3, 1), top = textGrob(title, gp = gpar(fontsize = 25)))
      ## d <- dev.off()
    }
    
    ##        ret.list <- list(
    ##            "barplots" = barplots)
    ##        return(ret.list)
    
    
    ## Spot check that the methods have the same scores for coarse- and fine-grained
    ## NB: some methods may differ between coarse and fine-grained; pick several
    ## baseline/comparator methods that we know to be the same
    check.methods <- c("CIBERSORT", "MCP-counter", "CIBERSORTx")
    check.methods <- sort(unique(c(as.character(means.by.cell.type.method[["coarse"]][["pearson"]][, method.name.col]),
                                   as.character(means.by.cell.type.method[["fine"]][["pearson"]][, method.name.col]))))
    for(meth in check.methods) {
      
      res.coarse <- means.by.cell.type.method[["coarse"]][["pearson"]]
      res.fine <- means.by.cell.type.method[["fine"]][["pearson"]]            
      
      meth.res.coarse <- res.coarse[res.coarse[, method.name.col] == meth, ]
      meth.res.fine <- res.fine[res.fine[, method.name.col] == meth, ]            
      m <- merge(meth.res.coarse, meth.res.fine, by = c(cell.type.col))
      if(nrow(m) == 0) { next }
      m$diff <- m$cor.x - m$cor.y
      m <- m[!is.na(m$cor.x),]
      eps <- 10^-4
      cat(paste0(meth, " ", postfix, " max diff between coarse and fine-grained is: ", max(abs(m$diff)), "\n"))
      flag <- abs(m$diff) > eps 
      if(any(flag)) {
        print(head(m[flag,,drop=F]))
        cat(paste0("Max diff exceeded for ", meth, " ", postfix, ": ", max(abs(m$diff)), "\n"))
      }
    }
    
    cat(paste0("Plotting heatmaps\n"))
    heatmaps <- list()
    spearman.heatmaps <- list()
    pearson.fc.heatmaps <- list()
    method.levels <- list()
    cell.type.levels <- list()
    for(sub.challenge in sub.challenges) {
      for(cor.type in c("pearson", "spearman", "pearson.fc")) {
        means <- means.by.cell.type.method[[sub.challenge]][[cor.type]]
        
        method.levels[[paste0(sub.challenge, "-", cor.type)]] <-
          calculate.method.levels(means, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
        cell.type.levels[[paste0(sub.challenge, "-", cor.type)]] <-
          calculate.cell.type.levels(means, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
        
        exclude.method <-
          ddply(means, .variables = c(method.name.col),
                .fun = function(df) {
                  exclude <- all(is.na(df$cor))
                  data.frame(exclude = exclude)
                })
        if(any(exclude.method$exclude)) {
          means <- means[!(means[, method.name.col] %in% subset(exclude.method, exclude == TRUE)[, method.name.col]),]
          method.levels[[paste0(sub.challenge, "-", cor.type)]] <- 
            method.levels[[paste0(sub.challenge, "-", cor.type)]][!(method.levels[[paste0(sub.challenge, "-", cor.type)]] %in% subset(exclude.method, exclude == TRUE)[, method.name.col])]
        }
        cat("means\n")
        print(unique(means[, method.name.col]))
        cat("levels\n")
        print(method.levels[[paste0(sub.challenge, "-", cor.type)]])
        cor.type.label <- "Pearson\nCorrelation"
        if(cor.type == "spearman") {
          cor.type.label <- "Spearman\nCorrelation"
        }
        if(cor.type == "pearson.fc") {
          cor.type.label <- "Pearson Correlation\n(Fold Change)"
        }
        g <- plot.cell.type.correlation.heatmap(means, show.corr.text = TRUE,
                                                id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor",
                                                second.col.summary.fun = "mean",
                                                cor.type.label = cor.type.label,
                                                method.levels = method.levels[[paste0(sub.challenge, "-", cor.type)]],
                                                cell.type.levels = cell.type.levels[[paste0(sub.challenge, "-", cor.type)]], ids.to.bold = comparator.methods)
        ##            g <- plot.cell.type.correlation.strip.plots(means, show.corr.text = TRUE, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
        if(cor.type == "pearson") { 
          heatmaps[[sub.challenge]] <- g
        } else if(cor.type == "spearman") {
          spearman.heatmaps[[sub.challenge]] <- g	    
        } else if(cor.type == "pearson.fc") {
          pearson.fc.heatmaps[[sub.challenge]] <- g
        } else {
          stop(paste0("Unknown cor.type ", cor.type, "\n"))
        }
        title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge")
        round.text <- make.round.text(round)
        title <- paste0(title, " (", round.text, ")")
        
        g <- g + ggtitle(title)
        ## png(paste0(figs.dir, "rerun-validation-bootstrap-cell-heatmap-", sub.challenge, postfix, ".png"), width = 2 * 480)
        ## print(g)
        ## d <- dev.off()
      } # for cor.type
    }
    
    for(cor.type in c("pearson", "spearman", "pearson.fc")) {
      coarse.means <- means.by.cell.type.method[["coarse"]][[cor.type]]
      fine.means <- means.by.cell.type.method[["fine"]][[cor.type]]        
      all.means <- rbind(coarse.means, fine.means)
      
      all.means <-
        ddply(all.means, .variables = c(method.name.col, cell.type.col),
              .fun = function(df) {
                data.frame(cor = summary.fun(df$cor))
              })
      
      method.levels[[paste0("merged", "-", cor.type)]] <-
        calculate.method.levels(all.means, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
      
      cell.type.levels[[paste0("merged", "-", cor.type)]] <-
        calculate.cell.type.levels(all.means, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
      
      exclude.method <-
        ddply(all.means, .variables = c(method.name.col),
              .fun = function(df) {
                exclude <- all(is.na(df$cor))
                data.frame(exclude = exclude)
              })
      cat(paste0("Exclude: ", cor.type, "\n"))
      print(exclude.method)
      if(any(exclude.method$exclude)) {
        all.means <- all.means[!(all.means[, method.name.col] %in% subset(exclude.method, exclude == TRUE)[, method.name.col]),]
        method.levels[[paste0("merged", "-", cor.type)]] <- 
          method.levels[[paste0("merged", "-", cor.type)]][!(method.levels[[paste0("merged", "-", cor.type)]] %in% subset(exclude.method, exclude == TRUE)[, method.name.col])]
      }
      
      cor.type.label <- "Pearson\nCorrelation"
      if(cor.type == "spearman") {
        cor.type.label <- "Spearman\nCorrelation"
      }
      if(cor.type == "pearson.fc") {
        cor.type.label <- "Pearson Correlation\n(Fold Change)"
      }
      
      g <- plot.cell.type.correlation.heatmap(all.means, show.corr.text = TRUE,
                                              id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor",
                                              cor.type.label = cor.type.label,						
                                              second.col.summary.fun = "mean",
                                              method.levels = method.levels[[paste0("merged", "-", cor.type)]],
                                              cell.type.levels = cell.type.levels[[paste0("merged", "-", cor.type)]], ids.to.bold = comparator.methods)
      merged.all.means <- all.means
      if(cor.type == "pearson") { 
        heatmaps[["merged"]] <- g
      } else if(cor.type == "spearman") {
        spearman.heatmaps[["merged"]] <- g	
      } else if(cor.type == "pearson.fc") {
        pearson.fc.heatmaps[["merged"]] <- g
      } else {
        stop(paste0("Unknown cor.type ", cor.type, "\n"))
      }
    } # for cor.type
    
    cat("Creating merged means.over.dataset\n")
    coarse.means <- means.over.dataset[["coarse"]][["pearson"]]
    fine.means <- means.over.dataset[["fine"]][["pearson"]]        
    all.means <- rbind(coarse.means, fine.means)
    
    all.means <-
      ddply(all.means, .variables = c(method.name.col, cell.type.col, "boot.i"),
            .fun = function(df) {
              data.frame(cor = summary.fun(df$cor))
            })
    
    cat(paste0("Plotting strip plots\n"))
    nms <- list("coarse" = "coarse", "fine" = "fine", "coarse-priority" = "coarse-priority", "fine-priority" = "fine-priority",
                "merged" = "merged", "merged-priority" = "merged-priority")
    strip.plots <-
      llply(nms,
            .parallel = FALSE,
            .fun = function(nm) {
              sub.challenge <- NA
              df <- NULL
              entry <- NULL
              lvl.entry <- NULL
              if(grepl(nm, pattern="coarse")) {
                entry <- "coarse"
                lvl.entry <- paste0(entry, "-pearson")
                df <- means.over.dataset[[entry]][["pearson"]]
              }
              if(grepl(nm, pattern="fine")) {
                entry <- "fine"
                lvl.entry <- paste0(entry, "-pearson")
                df <- means.over.dataset[[entry]][["pearson"]]                          
              }
              if(grepl(nm, pattern="merged")) {
                entry <- "merged"
                lvl.entry <- paste0(entry, "-pearson")
                df <- all.means
              }
              
              g <- NULL
              if(grepl(nm, pattern="priority")) {
                flag <- df[, method.name.col] %in% priority.methods
                ret <- plot.strip.plots(df[flag, ], id.var = method.name.col, cell.type.var = cell.type.col, var = "cor",
                                        method.levels = method.levels[[lvl.entry]],
                                        cell.type.levels = rev(cell.type.levels[[lvl.entry]]),
                                        label = "Pearson Correlation")
                g <- ret[["g"]]
                df <- ret[["df"]]
                lvls <- levels(df[,method.name.col])
                lvls <- lvls[lvls %in% df[,method.name.col]]
                y.bold.labels <- ifelse(lvls %in% comparator.methods, yes = "bold", no = "plain")
                print("here\n")
                print(comparator.methods)
                print(priority.methods)
                print(lvls)
                print(y.bold.labels)
              } else {
                # Exclude ensemble, which has NAs for pearson
                flag <- !(df[, method.name.col] %in% c("consensus rank", "ensemble"))
                ret <- plot.strip.plots(df[flag,], id.var = method.name.col, cell.type.var = cell.type.col, var = "cor",
                                        method.levels = method.levels[[lvl.entry]],
                                        cell.type.levels = rev(cell.type.levels[[lvl.entry]]),
                                        label = "Pearson Correlation")
                g <- ret[["g"]]
                df <- ret[["df"]]
                lvls <- levels(df[,method.name.col])
                lvls <- lvls[lvls %in% df[,method.name.col]]
                y.bold.labels <- ifelse(lvls %in% comparator.methods, yes = "bold", no = "plain")
                print("there\n")
                print(comparator.methods)
                print(lvls)
                print(y.bold.labels)
                
                
              }
              g <- g + theme(axis.text.y = element_text(face = y.bold.labels))
              # g <- g + scale_y_discrete(labels = function(x) bold.highlight(x, comparator.methods))
              # g <- g + theme(axis.text.y=element_markdown())
            })
    
    ## "boxplots" = boxplots,
    ret.list <- list("median.bootstrapped.scores" = median.bootstrapped.scores,
                     "mean.bootstrapped.scores" = mean.bootstrapped.scores,
                     "barplots" = barplots,
                     "strip.plots" = strip.plots,
                     "heatmaps" = heatmaps,
                     "spearman.heatmaps" = spearman.heatmaps,			 
                     "pearson.fc.heatmaps" = pearson.fc.heatmaps,			 
                     "merged.all.means" = merged.all.means)			 
    
    ret.list
    
  }

load.challenge.validation.data <- function(validation_path_base) { 
  # Location of Challenge validation admixtures on Synapse
  total_data_path <- paste0(validation_path_base,"/total_validation_data" )
  # List all files in the Challenge validation folder
  children <- list.files(total_data_path, full.names = F)
  l <- as.list(children)
  df <- do.call(rbind.data.frame, l)
  colnames(df) <-  c("name")
  # Load the "input.csv" file that lists the expression matrices in each
  # of the Challenge validation datasets
  input.tbl <- read.table(paste0(total_data_path, "/input.csv"), sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)
  # Subset to the expression data with Hugo symbols
  input.files <- input.tbl[, c("dataset.name", "hugo.expr.file")]
  df <- merge(df, input.files, by.x = "name", by.y = "hugo.expr.file")
  
  # Download each dataset from Synapse and read it in
  datasets <- 
    llply(df$name, 
          .fun = function(file_name) {
            obj <-read.csv(paste0(total_data_path, "/", file_name), row.names = 1)
            mat <- as.data.frame(obj)
          })
  names(datasets) <- df$dataset.name
  datasets
}


# Return a list with named entries:
# predictions: the predictions from the comparator and participant methods, with columns
#              dataset.name, subchallenge, sample.id, cell.type, prediction, method.name, and submission (1, 2, 3, or latest)
# dataset.anno: the description of each dataset, with columns
#              dataset.name, tumor.type, distribution.type, mixture.type
# ground.truth: the ground truth, with columns:
#              dataset.name, subchallenge, sample.id, cell.type, measured
# In all cases, limit results to target.submissions and target.subchallenges
load.challenge.validation.results <- function(validation_path_base, target.submissions = c("1", "2", "3", "latest"), 
                                              target.subchallenges = c("fine", "coarse"),
                                              method.name.col = "method.name") {
  data_path_to_read <-  paste0(validation_path_base, "/rerun-validation-sanitized-predictions.csv")  #"syn22320329"
  obj <- read.csv(data_path_to_read)
  res.all <- as.data.frame(obj)
  
  flag <- res.all[,method.name.col] %in% methods.to.exclude
  cat(paste0("Excluding methods: ", paste(unique(res.all[flag, method.name.col]), collapse = ", "), "\n"))
  res.all <- res.all[!flag,]
  
  res.all <- subset(res.all, submission %in% target.submissions)
  res.all <- subset(res.all, subchallenge %in% target.subchallenges)
  
  preds <- res.all[, c("dataset.name", "subchallenge", "sample.id", "cell.type", "prediction", "method.name", "submission")]
  anno <- unique(res.all[, c("dataset.name", "tumor.type", "distribution.type", "mixture.type")])
  rownames(anno) <- NULL
  gt <- unique(res.all[, c("dataset.name", "subchallenge", "sample.id", "cell.type", "measured")])
  rownames(gt) <- NULL
  
  ret <- list("predictions" = preds, "dataset.anno" = anno, "ground.truth" = gt)
  ret
}

## Ensure we have a prediction (even if it is NA) for all cell types in all datasets by all methods
make.missing.predictions.na <- function(res.all, cell.type.col = "cell.type", prediction.col = "prediction",
                                        measured.col = "measured", subchallenge.col = "subchallenge",
                                        sample.id.col = "sample.id", dataset.name.col = "dataset.name",
                                        round.col = "submission") {
  tmp <- unique(res.all[, !(colnames(res.all) %in% c(cell.type.col, prediction.col, measured.col))])
  cell.types.by.sub.challenge <- unique(res.all[, c(cell.type.col, subchallenge.col, sample.id.col, dataset.name.col)])
  tmp <- merge(tmp, cell.types.by.sub.challenge, all = TRUE)
  measured.tbl <- unique(res.all[, c(cell.type.col, subchallenge.col, sample.id.col, dataset.name.col, measured.col)])
  prediction.tbl <- res.all[, c(cell.type.col, subchallenge.col, sample.id.col, dataset.name.col, method.name.col, round.col, prediction.col)]
  tmp <- merge(tmp, measured.tbl, all = TRUE)
  tmp <- merge(tmp, prediction.tbl, all = TRUE)
  res.all <- tmp
  flag <- !is.na(res.all[, measured.col])
  res.all <- res.all[flag, ]
  
  res.all[, cell.type.col] <- as.character(res.all[, cell.type.col])
  res.all
}

# Aggregate predictions from raw cell types as output by a deconvolution method
# into those cell types expected of the Challenge.
# Aggregation is performed via summing.
# The translation between raw and Challenge subtypes is encoded in the translation table,
# where the raw.cell.type column gives the raw cell type name and cell.type gives the Challenge cell type.
# Results in res_df are described by columns method.name, dataset.name, sample.id, cell.type, and prediction
aggregate.cell.type.predictions <- function(res_df, translation_df) {
  nm <- unique(res_df$method.name)
  res_df <- res_df %>%
    dplyr::inner_join(translation_df) %>%
    dplyr::select(dataset.name, sample.id, cell.type, prediction) %>%
    dplyr::group_by(dataset.name, sample.id, cell.type) %>%
    dplyr::summarise(prediction = sum(prediction)) %>%
    dplyr::ungroup()
  res_df <- as.data.frame(res_df)
  res_df <- cbind(method.name = nm, res_df)
}

sub.challenges <- list("coarse" = "coarse", "fine" = "fine")

calculate.empirical.bayes <-
  function(df, col.id, numerator.id, denominator.id, sample.id.cols, score.col) {
    flag <- df[, col.id] %in% c(numerator.id, denominator.id)
    tmp <- df[flag, ]
    n.tot <- nrow(tmp) / 2
    n.num <-
      sum(unlist(
        dlply(tmp, .variables = sample.id.cols,
              .fun = function(df.2) {
                if(nrow(df.2) != 2) {
                  print(df.2)
                  print(c(numerator.id, denominator.id))
                  stop("Was expecting 2 rows\n")
                }
                if(any(is.na(df.2[, score.col]))) { stop("Got NA scores\n") }
                rownames(df.2) <- df.2[, col.id]
                diff <- df.2[numerator.id, score.col] -
                  df.2[denominator.id, score.col]
                if(diff > 0) { return(1) }
                return(0)
              })))
    n.num / (n.tot - n.num)
  }

do.bootstrap.analysis <-
  function(res.input, bootstraps, 
           method.name.col, 
           subchallenge.col, measured.col, cell.type.col,
           dataset.name.col, sample.id.col, prediction.col,
           round.col, round = "latest") {
    
    submitter.tbl <- unique(res.input[, c(method.name.col, round.col, subchallenge.col), drop = FALSE])
    or <- order(submitter.tbl[, round.col])
    submitter.tbl <- submitter.tbl[or, ]
    submitter.tbl[, round.col] <- as.character(submitter.tbl[, round.col])
    flag <- submitter.tbl[, round.col] %in% c("latest", as.character(round))
    submitter.tbl <- submitter.tbl[flag, ]
    flag <- !duplicated(submitter.tbl[, c(method.name.col, subchallenge.col)], fromLast = FALSE)
    submitter.tbl <- submitter.tbl[flag, ]
    
    res.round <- merge(res.input, submitter.tbl, by = c(method.name.col, subchallenge.col, round.col))
    
    method.id.col <- method.name.col
    
    if(FALSE) {
      # Add ensemble here -- no need. We report bootstrapped scores
      # > colnames(results[["3"]][["res.round"]])
      # [1] "method.name"       "subchallenge"      "submission"       
      # [4] "dataset.name"      "sample.id"         "cell.type"        
      # [7] "objectId"          "comparator"        "submitterId"      
      # [10] "repo_name"         "tumor.type"        "distribution.type"
      # [13] "mixture.type"      "measured"          "prediction"       
      res.ensemble <-
        ddply(res.round[, !(colnames(res.round) %in% c("objectId", "comparator", "submittedId", "repo_name"))],
              .variables = c(round.col, subchallenge.col, dataset.name.col),
              .fun = function(df.in) {
                tres <- 
                  ddply(df.in, .variables = method.name.col,
                        .fun = function(df) {
                          df.ret <- data.frame(prediction = df[, prediction.col],
                                               pred.rank = rank(df[, prediction.col], ties.method="first"),
                                               sample.id = df[, sample.id.col])
                          colnames(df.ret) <- c(prediction.col, "pred.rank", sample.id.col)
                          df.ret
                        })
                # take the consensus (here, just mean) rank across methods
                cons <- 
                  ddply(tres, .variables = sample.id.col,
                        .fun = function(df) {
                          df.ret <- data.frame(cons.rank = mean(df$pred.rank))
                          colnames(df.ret)[1] <- prediction.col
                          df.ret
                        })
                tmp <- unique(df.in[, c(sample.id.col, measured.col, subchallenge.col, round.col, dataset.name.col, cell.type.col, "tumor.type", "distribution.type", "mixture.type")])
                cons <- merge(cons, tmp)
                cons$method.name <- "ensemble"
                cons$objectId <- "dummy"
                cons$comparator <- FALSE
                cons$submitterId <- "bwhite"
                cons$repo_name <- "dummy"
                ret <- cons[, c(method.name.col, subchallenge.col, round.col, dataset.name.col, sample.id.col, cell.type.col, 
                                "objectId", "comparator", "submitterId", "repo_name", "tumor.type", "distribution.type", "mixture.type",
                                measured.col, prediction.col)]
                ret
              })
      print(head(res.round))
      print(head(res.ensmble))
      res.round <- rbind(res.round, res.ensemble)
    } #end if(FALSE)
    res <- res.round
    
    tbls <-
      llply(sub.challenges,
            .fun = function(sub.challenge) {
              flag <- res[, subchallenge.col] == sub.challenge
              tmp <- res[flag, ]
              tmp
            })
    
    tbls.by.cell <-
      llply(sub.challenges,
            .fun = function(sub.challenge) {
              tmp <- tbls[[sub.challenge]]
              flag <- !is.na(tmp[, measured.col])
              tmp <- tmp[flag, ]
              cells <- unique(tmp[, cell.type.col])
              names(cells) <- cells
              llply(cells,
                    .fun = function(ct) {
                      flag <- tmp[, cell.type.col] == ct
                      ret <- tmp[flag, ]
                      ret$id <- paste0(ret[, dataset.name.col], "-", ret[, sample.id.col])
                      ret
                    })
            })
    
    ## Calculate both pearson and spearman correlation over bootstraps, within dataset and cell type
    cat(paste0("Calculating bootstrap correlations (ensemble)\n"))
    bootstrapped.cors.ensemble <-
      llply(sub.challenges, 
            .parallel = FALSE,
            .fun = function(sub.challenge) {
              n.bootstraps <- length(bootstraps[[sub.challenge]])
              indices <- 1:n.bootstraps
              # indices <- 1:10
              names(indices) <- indices
              
              ret.i <-
                ldply(indices,
                      .parallel = TRUE,
                      .fun = function(i) {
                        ret.cons <-
                          ddply(tbls[[sub.challenge]], .variables = c(cell.type.col, dataset.name.col),
                                .fun = function(df.in) {
                                  # method.name is row
                                  x <- acast(df.in[, c(prediction.col, method.name.col, sample.id.col)], 
                                             formula = paste0(method.name.col, " ~ ", sample.id.col), value.var = prediction.col)
                                  
                                  colnames(x) <- paste0(df.in[1, dataset.name.col], "-", colnames(x))
                                  sample.ids <- bootstraps[[sub.challenge]][[i]]$id
                                  sample.ids <- sample.ids[sample.ids %in% colnames(x)]
                                  x <- x[, sample.ids]
                                  m <- melt(x)
                                  colnames(m) <- c(method.name.col, sample.id.col, prediction.col)
                                  tres <- 
                                    ddply(m, .variables = method.name.col,
                                          .fun = function(df) {
                                            df.ret <- data.frame(prediction = df[, prediction.col],
                                                                 pred.rank = rank(df[, prediction.col], ties.method="first"),
                                                                 sample.id = df[, sample.id.col])
                                            colnames(df.ret) <- c(prediction.col, "pred.rank", sample.id.col)
                                            df.ret
                                          })
                                  # take the consensus (here, just mean) rank across methods
                                  cons <- 
                                    ddply(tres, .variables = sample.id.col,
                                          .fun = function(df) {
                                            df.ret <- data.frame(cons.rank = mean(df$pred.rank))
                                            df.ret
                                          })
                                  tmp <- unique(df.in[, c(sample.id.col, measured.col)])
                                  tmp[, sample.id.col] <- paste0(df.in[1, dataset.name.col], "-", tmp[, sample.id.col])
                                  cons <- merge(cons, tmp)
                                  ret <- data.frame(method.name = "ensemble", pearson = NA, rmse = NA, spearman = cor(cons$cons.rank, cons[, measured.col]), pearson.fc = NA)
                                  colnames(ret) <- c(method.name.col, "pearson", "rmse", "spearman", "pearson.fc")
                                  ret
                                })
                      })
              # "method.name"  "boot.i"       "dataset.name" "cell.type"    "pearson"   "spearman"     "rmse" "pearson.fc"
              colnames(ret.i)[1] <- "boot.i"
              ret.i <- ret.i[, c(method.name.col, "boot.i", dataset.name.col, cell.type.col, "pearson", "spearman", "rmse", "pearson.fc")]
              ret.i 
            })
    cat(paste0("Calculating bootstrap correlations\n"))
    bootstrapped.cors <-
      llply(sub.challenges, 
            .fun = function(sub.challenge) {
              tmp <- tbls[[sub.challenge]]                      
              flag <- !is.na(tmp[, measured.col])
              tmp <- tmp[flag, ]
              methods <- unique(tmp[, method.id.col])
              ##                      methods <- c("Aginome-XMU", "CIBERSORTx")
              names(methods) <- methods
              n.bootstraps <- length(bootstraps[[sub.challenge]])
              indices <- 1:n.bootstraps
              names(indices) <- indices
              ret.all <-
                ldply(methods, .parallel = TRUE,
                      .fun = function(method) {
                        print(method)
                        ret.i <-
                          ldply(indices,
                                .fun = function(i) {
                                  ret.method <-
                                    ldply(tbls.by.cell[[sub.challenge]],
                                          .fun = function(df.in) {
                                            flag <- df.in[, method.id.col] == method
                                            df <- df.in[flag, ]
                                            if(any(duplicated(df$id))) {
                                              print(head(df[my.dup(df$id),]))
                                              stop("stop")
                                            }
                                            rownames(df) <- df$id
                                            sample.ids <- bootstraps[[sub.challenge]][[i]]$id
                                            sample.ids <- sample.ids[sample.ids %in% df$id]
                                            if(!(all(sample.ids %in% rownames(df)))) {
                                              stop("Some sample ids not in df\n")
                                            }
                                            df <- df[sample.ids,]
                                            df
                                          })
                                  ret.method <- ret.method[, -1]
                                  score <-
                                    ddply(ret.method, .variables = c(dataset.name.col),
                                          .fun = function(df.ds) {
                                            tmp <- ddply(df.ds, .variables = c(cell.type.col),
                                                         .fun = function(df.ct) {
                                                           if(any(is.na(df.ct[,2]))) {
                                                             print(df.ct)
                                                           }
                                                           mts <- c("pearson", "spearman", "rmse", "pearson.fc")
                                                           names(mts) <- mts
                                                           vals <- llply(mts,
                                                                         .fun = function(comp.method) {
                                                                           pred <- as.numeric(df.ct[, prediction.col])
                                                                           measured <- as.numeric(df.ct[, measured.col])
                                                                           
                                                                           if(any(is.na(pred))) {
                                                                             if(!all(is.na(pred))) {
                                                                               stop("Some but not all NA\n")
                                                                             }
                                                                             return(NA)
                                                                           }
                                                                           if(var(pred) == 0) { return(0) }
                                                                           
                                                                           val <- NA
                                                                           if(comp.method %in% c("pearson", "spearman")) {
                                                                             val <- cor(pred, measured,
                                                                                        method = comp.method)
                                                                           } else if(comp.method == "pearson.fc") {
                                                                             # calculate the pearson correlation of log fold changes (across successive samples)
                                                                             # 1. Order the samples based on the their ground truth values g1 < g2 < … < gn (i.e., measured)
                                                                             # 2. Compute ground truth fold changes g2/g1, g3/g2, …, gn/gn-1
                                                                             # 3. Compute prediction fold changes p2/p1, p3/p2, …, pn/pn-1
                                                                             # 4. Compute the pearson correlation of these two vectors
                                                                             o <- order(measured, decreasing=FALSE)
                                                                             eps <- 10^-5
                                                                             gt.ordered <- measured[o] + eps
                                                                             # NB: the ordering is established by the _ground truth_ values and applied to those
                                                                             # values and the predicted values
                                                                             pred.ordered <- pred[o] + eps
                                                                             gt.fc <- unlist(lapply(2:length(gt.ordered), function(indx) gt.ordered[indx]/gt.ordered[indx-1]))
                                                                             pred.fc <- unlist(lapply(2:length(pred.ordered), function(indx) pred.ordered[indx]/pred.ordered[indx-1]))
                                                                             val <- cor(pred.fc, gt.fc, method = "pearson")
                                                                           } else if(comp.method == "rmse") {
                                                                             val <- sqrt(mean((pred-measured)^2))
                                                                           } else {
                                                                             stop(paste0("Unknown method ", comp.method, "\n"))
                                                                           }
                                                                           val
                                                                         })
                                                           as.data.frame(vals)
                                                         })
                                            colnames(tmp)[1] <- cell.type.col
                                            tmp
                                          })
                                  score
                                })
                        colnames(ret.i)[1] <- "boot.i"
                        ret.i
                      })
              colnames(ret.all)[1] <- method.id.col
              ret.all
            })
    
    for(nm in names(bootstrapped.cors)) {
      bootstrapped.cors[[nm]] <- rbind(bootstrapped.cors[[nm]], bootstrapped.cors.ensemble[[nm]])
    }
    
    print(colnames(bootstrapped.cors[[1]]))
    # "method.name"  "boot.i"       "dataset.name" "cell.type"    "pearson"   "spearman"     "rmse" "pearson.fc"
    print(head(bootstrapped.cors[[1]]))
    
    cat(paste0("Calculating bootstrapped scores\n"))
    bootstrapped.scores <-
      llply(sub.challenges, .parallel = TRUE,
            .fun = function(sub.challenge) {
              df <- bootstrapped.cors[[sub.challenge]]
              flag <- res.round[, subchallenge.col] == sub.challenge
              un <- unique(res.round[flag, unique(c(method.id.col, method.name.col))])
              df <- merge(df, un)
              ## Average over cell type (within method and dataset and bootstrap sample)
              df <- ddply(df,
                          .variables = c(method.id.col, dataset.name.col, "boot.i"),
                          .fun = function(tmp) {
                            data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman), rmse = mean(tmp$rmse), pearson.fc = mean(tmp$pearson.fc))
                          })
              ## Average over dataset (within method and bootstrap sample)
              df <- ddply(df,
                          .variables = c(method.id.col, "boot.i"),
                          .fun = function(tmp) {
                            data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman), rmse = mean(tmp$rmse), pearson.fc = mean(tmp$pearson.fc))
                          })
              df
            })
    
    cat(paste0("Calculating mean bootstrapped scores\n"))
    mean.bootstrapped.scores <-
      llply(bootstrapped.scores, .parallel = TRUE,
            .fun = function(df) {
              ## Average over bootstraps
              df <- ddply(df,
                          .variables = c(method.id.col),
                          .fun = function(tmp) {
                            data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman), rmse = mean(tmp$rmse), pearson.fc = mean(tmp$pearson.fc))
                          })
              o <- order(df$pearson, decreasing = TRUE)
              df <- df[o,]
              df
            })
    
    means.over.dataset <-
      llply(bootstrapped.cors,
            .fun = function(df) {
              methods <- c("pearson", "spearman", "rmse", "pearson.fc")
              na.rm <- FALSE
              names(methods) <- methods
              res <- llply(methods,
                           .fun = function(method) {
                             ## average over dataset
                             ret <- ddply(df, .variables = c(method.id.col, cell.type.col, "boot.i"),
                                          .fun = function(df) {
                                            data.frame(cor = mean(df[, method], na.rm=na.rm))
                                          })
                           })
            })
    
    means.over.bootstrap <-
      llply(bootstrapped.cors,
            .fun = function(df) {
              methods <- c("pearson", "spearman", "rmse", "pearson.fc")
              na.rm <- FALSE
              names(methods) <- methods
              res <- llply(methods,
                           .fun = function(method) {
                             ## average over dataset
                             ret <- ddply(df, .variables = c(method.id.col, cell.type.col, dataset.name.col),
                                          .fun = function(df) {
                                            data.frame(cor = mean(df[, method], na.rm=na.rm))
                                          })
                           })
            })
    
    cat(paste0("Calculating mean by cell type\n"))
    if(FALSE) {
      means.by.cell.type.method <-
        llply(bootstrapped.cors,
              .fun = function(df) {
                methods <- c("pearson", "spearman", "rmse", "pearson.fc")
                na.rm <- FALSE
                names(methods) <- methods
                res <- llply(methods,
                             .fun = function(method) {
                               ## first, average over bootstrap
                               ret <- ddply(df, .variables = c(method.id.col, cell.type.col, dataset.name.col),
                                            .fun = function(df) {
                                              data.frame(cor = mean(df[, method], na.rm=na.rm))
                                            })
                               ## now, average over dataset
                               ret2 <- ddply(ret, .variables = c(method.id.col, cell.type.col),
                                             .fun = function(df) {
                                               data.frame(cor = mean(df$cor, na.rm=na.rm))
                                             })
                             })
              })
    }
    
    means.by.cell.type.method <-
      llply(means.over.dataset,
            .fun = function(df) {
              methods <- c("pearson", "spearman", "rmse", "pearson.fc")
              na.rm <- FALSE
              names(methods) <- methods
              res <- llply(methods,
                           .fun = function(method) {
                             ## average over bootstrap (means.over.dataset has already been averaged over dataset)
                             ret <- ddply(df[[method]], .variables = c(method.id.col, cell.type.col),
                                          .fun = function(df) {
                                            data.frame(cor = mean(df$cor, na.rm=na.rm))
                                          })
                           })
            })
    
    cat(paste0("Calculating bayes factors\n"))
    top.performers <-
      llply(sub.challenges,
            .fun = function(sub.challenge) {
              mean.scores <- mean.bootstrapped.scores[[sub.challenge]]
              mean.scores <- na.omit(mean.scores)
              scores <- bootstrapped.scores[[sub.challenge]]
              numerator.indx <- 1
              numerator.id <- mean.scores[numerator.indx, method.id.col]
              numerator.id
            })
    
    bayes.factors <-
      llply(sub.challenges,
            .fun = function(sub.challenge) {
              mean.scores <- mean.bootstrapped.scores[[sub.challenge]]
              mean.scores <- na.omit(mean.scores)
              scores <- bootstrapped.scores[[sub.challenge]]
              numerator.indx <- 1
              numerator.id <- mean.scores[numerator.indx, method.id.col]
              indices <- 1:nrow(mean.scores)
              indices <- indices[!(indices == numerator.indx)]
              names(indices) <- mean.scores[indices, method.id.col]
              ret <-
                ldply(indices,
                      .fun = function(i) {
                        denominator.id <- mean.scores[i, method.id.col]
                        methods <- c("pearson", "spearman")
                        names(methods) <- methods
                        res <- llply(methods,
                                     .fun = function(method) {
                                       bf <- calculate.empirical.bayes(scores, col.id = method.id.col,
                                                                       numerator.id = numerator.id,
                                                                       denominator.id = denominator.id,
                                                                       sample.id.cols = "boot.i",
                                                                       score.col = method)
                                     })
                        ret.df <- as.data.frame(res)
                        ret.df
                      })
              ret
            })
    
    ret.list <- list("res.round" = res.round,
                     "bootstrapped.cors" = bootstrapped.cors,
                     "bootstrapped.scores" = bootstrapped.scores,
                     "mean.bootstrapped.scores" = mean.bootstrapped.scores,
                     "means.by.cell.type.method" = means.by.cell.type.method,
                     "means.over.dataset" = means.over.dataset,
                     "means.over.bootstrap" = means.over.bootstrap,                         
                     "top.performers" = top.performers,
                     "bayes.factors" = bayes.factors)
    
    return(ret.list)
    
  }

## Modified slightly from
## https://stackoverflow.com/questions/53170465/how-to-make-a-base-r-style-boxplot-using-ggplot2
geom_boxplotMod <- function(mapping = NULL, data = NULL, stat = "boxplot", 
                            position = "dodge2", ..., outlier.colour = NULL, outlier.color = NULL, 
                            outlier.fill = NULL, outlier.shape = 1, outlier.size = 1.5, 
                            outlier.stroke = 0.5, outlier.alpha = NULL, notch = FALSE, notchwidth = 0.5,
                            varwidth = FALSE, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,
                            linetype = "dashed") # to know how these come here use: args(geom_boxplot)
{
  list(geom_boxplot(
    mapping = mapping, data = data, stat = stat, position = position,
    outlier.colour = outlier.colour, outlier.color = outlier.color, 
    outlier.fill = outlier.fill, outlier.shape = outlier.shape, 
    outlier.size = outlier.size, outlier.stroke = outlier.stroke, 
    outlier.alpha = outlier.alpha, notch = notch, 
    notchwidth = notchwidth, varwidth = varwidth, na.rm = na.rm, 
    show.legend = show.legend, inherit.aes = inherit.aes, linetype = 
      linetype, ...),
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.25),
    #the width of the error-bar heads are decreased
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.25),
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), ...)
  )
}


make.round.text <- function(round, round.str = "Round") {
  round.text <- ""
  if(round == "latest") {
    round.text <- paste0("Latest ", round.str)
  } else if (round == "1") {
    round.text <- paste0(round.str, " 1")
  } else {
    round.text <- paste0("Latest ", round.str, " up to ", round.str, " ", round)
  }
  round.text
}


plot.bootstrap.analysis <-
  function(res, bootstrapped.scores, mean.boostrapped.scores, median.bootstrapped.scores,
           means.by.cell.type.method,
           means.over.dataset, method.anno.round,
           postfix, plot.spearman.distribution = FALSE) {
    
    top.performers <- c("Aginome-XMU", "DA_505", "mitten_TDC19", "Biogem")
    comparator.methods <- get.comparators.cap()
    print(comparator.methods) 
    priority.methods <- unique(c(top.performers, comparator.methods))
    
    for(sub.challenge in sub.challenges) {
      print(sub.challenge)
      print(head(method.anno.round))
      print(method.anno.round)
      print(subchallenge.col)
      flag <- is.na(method.anno.round[, subchallenge.col]) | ( method.anno.round[, subchallenge.col] == sub.challenge )
      print(flag)
      method.anno.round.sc <- method.anno.round[flag, c(method.name.col, "Output", "Method")]
      print(method.anno.round.sc)
      bootstrapped.scores[[sub.challenge]] <-
        merge(bootstrapped.scores[[sub.challenge]], method.anno.round.sc, all.x = TRUE)
      median.bootstrapped.scores[[sub.challenge]] <-
        merge(median.bootstrapped.scores[[sub.challenge]], method.anno.round.sc, all.x = TRUE)
      mean.bootstrapped.scores[[sub.challenge]] <-
        merge(mean.bootstrapped.scores[[sub.challenge]], method.anno.round.sc, all.x = TRUE)
    }
    
    barplots <- list()
    for(sub.challenge in sub.challenges) {
      scores <- bootstrapped.scores[[sub.challenge]]
      median.scores <- median.bootstrapped.scores[[sub.challenge]]
      
      median.scores[, method.name.col] <- as.character(median.scores[, method.name.col])
      flag <- median.scores[, method.name.col] == "ensemble"
      median.scores[flag, method.name.col] <- "consensus rank"
      
      scores[, method.name.col] <- as.character(scores[, method.name.col])
      flag <- scores[, method.name.col] == "ensemble"
      scores[flag, method.name.col] <- "consensus rank"
      
      
      o <- order(median.scores$pearson)
      median.scores <- median.scores[o, ]
      flag <- ( !is.na(scores$pearson) & !is.na(scores$spearman) ) | ( scores[,method.name.col] == "consensus rank")
      scores <- scores[flag, ] 
      scores[, method.name.col] <- factor(scores[, method.name.col], levels = median.scores[, method.name.col])
      
      lvls <- levels(scores[,method.name.col]) 
      lvls <- lvls[lvls %in% scores[, method.name.col]]
      bold.labels <- ifelse(lvls %in% comparator.methods, yes = "bold", no = "plain")
      cat("scores[, method.name.col]\n")
      print(lvls)
      print(bold.labels)
      print(table(bold.labels))
      #            scores <- na.omit(scores)
      flag <- ( !is.na(median.scores$pearson) & !is.na(median.scores$spearman) ) | ( median.scores[,method.name.col] == "consensus rank") | ( median.scores[,method.name.col] == "ensemble")
      median.scores[, method.name.col] <- factor(median.scores[, method.name.col], levels = median.scores[, method.name.col])
      #            median.scores <- na.omit(median.scores)
      median.scores <- median.scores[flag, ] 
      
      
      g1 <- ggplot(data = scores, aes_string(x = method.name.col, y = "pearson"))
      g1 <- g1 + geom_boxplotMod(fill = "#56B4E9")
      g1 <- g1 + coord_flip()
      g1 <- g1 + xlab("Method")
      ## g1 <- g1 + ylab("Pearson Correlation")
      g1 <- g1 + ylab("Pearson")
      g1 <- g1 + ylim(c(-0.25, 1))
      g1 <- g1 + theme(text = element_text(size=18), title = element_text(size = 20), axis.text.y = element_text(face = bold.labels))
      # g1 <- g1 + theme(text = element_text(size=18), title = element_text(size = 20))
      # g1 <- g1 + scale_y_discrete(labels = function(x) bold.highlight(x, comparator.methods))
      # g1 <- g1 + theme(axis.text.y=element_markdown())
      
      if(plot.spearman.distribution) {
        g2 <- ggplot(data = scores, aes_string(x = method.name.col, y = "spearman"))
        g2 <- g2 + geom_boxplotMod(fill = "#E69F00")
      } else{
        g2 <- ggplot(data = median.scores)
        g2 <- g2 + geom_col(aes_string(x = method.name.col, y = "spearman"), fill = "#E69F00")
      }
      g2 <- g2 + coord_flip()
      g2 <- g2 + xlab("Method")
      ## g2 <- g2 + ylab("Spearman Correlation")
      g2 <- g2 + ylab("Spearman")
      g2 <- g2 + ylim(c(-0.25, 1))            
      g2 <- g2 + theme(text = element_text(size=18))    
      g2 <- g2 + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                       axis.ticks.y = element_blank())
      
      g3 <- ggplot(data = scores, aes_string(x = method.name.col, y = "pearson.fc"))
      g3 <- g3 + geom_boxplotMod(fill = "#56B4E9")
      g3 <- g3 + coord_flip()
      g3 <- g3 + xlab("Method")
      ## g3 <- g3 + ylab("Pearson Correlation")
      g3 <- g3 + ylab("Pearson (Fold Change)")
      g3 <- g3 + ylim(c(-0.25, 1))
      g3 <- g3 + theme(text = element_text(size=18), title = element_text(size = 20), axis.text.y = element_text(face = bold.labels))
      # g3 <- g3 + theme(text = element_text(size=18), title = element_text(size = 20))
      # g3 <- g3 + scale_y_discrete(labels = function(x) bold.highlight(x, comparator.methods))
      # g3 <- g3 + theme(axis.text.y=element_markdown())
      
      
      tmp <- scores[, c(method.name.col, "Output", "Method")]            
      ret <- plot.anno.heatmap.with.multiple.legends(tmp, "method.name", c("Method", "Output"), c("Set3", "Set1"))
      
      full.plot <- ret[["full.plot"]]
      for.first.legend <- ret[["legends"]][["Method"]]
      for.second.legend <- ret[["legends"]][["Output"]]
      
      legs <- plot_grid(get_legend(for.first.legend), get_legend(for.second.legend), nrow = 2, align = "v", rel_heights = c(2,1))
      
      ## pg <- plot_grid(g1, g2, full.plot, legs, nrow=1, align="h", rel_widths = c(3,1,0.5,0.5))
      
      barplots[[paste0(sub.challenge,  "-pearson")]] <- g1
      barplots[[paste0(sub.challenge,  "-spearman")]] <- g2            
      barplots[[paste0(sub.challenge,  "-pearson.fc")]] <- g3            
      barplots[[paste0(sub.challenge,  "-anno")]] <- full.plot
      barplots[[paste0(sub.challenge,  "-legend")]] <- legs
      
      title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge")
      round.text <- make.round.text(round)            
      title <- paste0(title, " (", round.text, ")")
      ## png(paste0(figs.dir, "rerun-validation-score-box-and-barplots-", sub.challenge, postfix, ".png"), width = 2 * 480)
      ## grid.arrange(g1, g2, nrow=1, widths = c(3, 1), top = textGrob(title, gp = gpar(fontsize = 25)))
      ## d <- dev.off()
    }
    
    ##        ret.list <- list(
    ##            "barplots" = barplots)
    ##        return(ret.list)
    
    
    ## Spot check that the methods have the same scores for coarse- and fine-grained
    ## NB: some methods may differ between coarse and fine-grained; pick several
    ## baseline/comparator methods that we know to be the same
    check.methods <- c("CIBERSORT", "MCP-counter", "CIBERSORTx")
    check.methods <- sort(unique(c(as.character(means.by.cell.type.method[["coarse"]][["pearson"]][, method.name.col]),
                                   as.character(means.by.cell.type.method[["fine"]][["pearson"]][, method.name.col]))))
    for(meth in check.methods) {
      
      res.coarse <- means.by.cell.type.method[["coarse"]][["pearson"]]
      res.fine <- means.by.cell.type.method[["fine"]][["pearson"]]            
      
      meth.res.coarse <- res.coarse[res.coarse[, method.name.col] == meth, ]
      meth.res.fine <- res.fine[res.fine[, method.name.col] == meth, ]            
      m <- merge(meth.res.coarse, meth.res.fine, by = c(cell.type.col))
      if(nrow(m) == 0) { next }
      m$diff <- m$cor.x - m$cor.y
      m <- m[!is.na(m$cor.x),]
      eps <- 10^-4
      cat(paste0(meth, " ", postfix, " max diff between coarse and fine-grained is: ", max(abs(m$diff)), "\n"))
      flag <- abs(m$diff) > eps 
      if(any(flag)) {
        print(head(m[flag,,drop=F]))
        cat(paste0("Max diff exceeded for ", meth, " ", postfix, ": ", max(abs(m$diff)), "\n"))
      }
    }
    
    cat(paste0("Plotting heatmaps\n"))
    heatmaps <- list()
    spearman.heatmaps <- list()
    pearson.fc.heatmaps <- list()
    method.levels <- list()
    cell.type.levels <- list()
    for(sub.challenge in sub.challenges) {
      for(cor.type in c("pearson", "spearman", "pearson.fc")) {
        means <- means.by.cell.type.method[[sub.challenge]][[cor.type]]
        
        method.levels[[paste0(sub.challenge, "-", cor.type)]] <-
          calculate.method.levels(means, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
        cell.type.levels[[paste0(sub.challenge, "-", cor.type)]] <-
          calculate.cell.type.levels(means, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
        
        exclude.method <-
          ddply(means, .variables = c(method.name.col),
                .fun = function(df) {
                  exclude <- all(is.na(df$cor))
                  data.frame(exclude = exclude)
                })
        if(any(exclude.method$exclude)) {
          means <- means[!(means[, method.name.col] %in% subset(exclude.method, exclude == TRUE)[, method.name.col]),]
          method.levels[[paste0(sub.challenge, "-", cor.type)]] <- 
            method.levels[[paste0(sub.challenge, "-", cor.type)]][!(method.levels[[paste0(sub.challenge, "-", cor.type)]] %in% subset(exclude.method, exclude == TRUE)[, method.name.col])]
        }
        cat("means\n")
        print(unique(means[, method.name.col]))
        cat("levels\n")
        print(method.levels[[paste0(sub.challenge, "-", cor.type)]])
        cor.type.label <- "Pearson\nCorrelation"
        if(cor.type == "spearman") {
          cor.type.label <- "Spearman\nCorrelation"
        }
        if(cor.type == "pearson.fc") {
          cor.type.label <- "Pearson Correlation\n(Fold Change)"
        }
        g <- plot.cell.type.correlation.heatmap(means, show.corr.text = TRUE,
                                                id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor",
                                                second.col.summary.fun = "mean",
                                                cor.type.label = cor.type.label,
                                                method.levels = method.levels[[paste0(sub.challenge, "-", cor.type)]],
                                                cell.type.levels = cell.type.levels[[paste0(sub.challenge, "-", cor.type)]], ids.to.bold = comparator.methods)
        ##            g <- plot.cell.type.correlation.strip.plots(means, show.corr.text = TRUE, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
        if(cor.type == "pearson") { 
          heatmaps[[sub.challenge]] <- g
        } else if(cor.type == "spearman") {
          spearman.heatmaps[[sub.challenge]] <- g	    
        } else if(cor.type == "pearson.fc") {
          pearson.fc.heatmaps[[sub.challenge]] <- g
        } else {
          stop(paste0("Unknown cor.type ", cor.type, "\n"))
        }
        title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge")
        round.text <- make.round.text(round)
        title <- paste0(title, " (", round.text, ")")
        
        g <- g + ggtitle(title)
        ## png(paste0(figs.dir, "rerun-validation-bootstrap-cell-heatmap-", sub.challenge, postfix, ".png"), width = 2 * 480)
        ## print(g)
        ## d <- dev.off()
      } # for cor.type
    }
    
    for(cor.type in c("pearson", "spearman", "pearson.fc")) {
      coarse.means <- means.by.cell.type.method[["coarse"]][[cor.type]]
      fine.means <- means.by.cell.type.method[["fine"]][[cor.type]]        
      all.means <- rbind(coarse.means, fine.means)
      
      all.means <-
        ddply(all.means, .variables = c(method.name.col, cell.type.col),
              .fun = function(df) {
                data.frame(cor = summary.fun(df$cor))
              })
      
      method.levels[[paste0("merged", "-", cor.type)]] <-
        calculate.method.levels(all.means, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
      
      cell.type.levels[[paste0("merged", "-", cor.type)]] <-
        calculate.cell.type.levels(all.means, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
      
      exclude.method <-
        ddply(all.means, .variables = c(method.name.col),
              .fun = function(df) {
                exclude <- all(is.na(df$cor))
                data.frame(exclude = exclude)
              })
      cat(paste0("Exclude: ", cor.type, "\n"))
      print(exclude.method)
      if(any(exclude.method$exclude)) {
        all.means <- all.means[!(all.means[, method.name.col] %in% subset(exclude.method, exclude == TRUE)[, method.name.col]),]
        method.levels[[paste0("merged", "-", cor.type)]] <- 
          method.levels[[paste0("merged", "-", cor.type)]][!(method.levels[[paste0("merged", "-", cor.type)]] %in% subset(exclude.method, exclude == TRUE)[, method.name.col])]
      }
      
      cor.type.label <- "Pearson\nCorrelation"
      if(cor.type == "spearman") {
        cor.type.label <- "Spearman\nCorrelation"
      }
      if(cor.type == "pearson.fc") {
        cor.type.label <- "Pearson Correlation\n(Fold Change)"
      }
      
      g <- plot.cell.type.correlation.heatmap(all.means, show.corr.text = TRUE,
                                              id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor",
                                              cor.type.label = cor.type.label,						
                                              second.col.summary.fun = "mean",
                                              method.levels = method.levels[[paste0("merged", "-", cor.type)]],
                                              cell.type.levels = cell.type.levels[[paste0("merged", "-", cor.type)]], ids.to.bold = comparator.methods)
      merged.all.means <- all.means
      if(cor.type == "pearson") { 
        heatmaps[["merged"]] <- g
      } else if(cor.type == "spearman") {
        spearman.heatmaps[["merged"]] <- g	
      } else if(cor.type == "pearson.fc") {
        pearson.fc.heatmaps[["merged"]] <- g
      } else {
        stop(paste0("Unknown cor.type ", cor.type, "\n"))
      }
    } # for cor.type
    
    cat("Creating merged means.over.dataset\n")
    coarse.means <- means.over.dataset[["coarse"]][["pearson"]]
    fine.means <- means.over.dataset[["fine"]][["pearson"]]        
    all.means <- rbind(coarse.means, fine.means)
    
    all.means <-
      ddply(all.means, .variables = c(method.name.col, cell.type.col, "boot.i"),
            .fun = function(df) {
              data.frame(cor = summary.fun(df$cor))
            })
    
    cat(paste0("Plotting strip plots\n"))
    nms <- list("coarse" = "coarse", "fine" = "fine", "coarse-priority" = "coarse-priority", "fine-priority" = "fine-priority",
                "merged" = "merged", "merged-priority" = "merged-priority")
    strip.plots <-
      llply(nms,
            .parallel = FALSE,
            .fun = function(nm) {
              sub.challenge <- NA
              df <- NULL
              entry <- NULL
              lvl.entry <- NULL
              if(grepl(nm, pattern="coarse")) {
                entry <- "coarse"
                lvl.entry <- paste0(entry, "-pearson")
                df <- means.over.dataset[[entry]][["pearson"]]
              }
              if(grepl(nm, pattern="fine")) {
                entry <- "fine"
                lvl.entry <- paste0(entry, "-pearson")
                df <- means.over.dataset[[entry]][["pearson"]]                          
              }
              if(grepl(nm, pattern="merged")) {
                entry <- "merged"
                lvl.entry <- paste0(entry, "-pearson")
                df <- all.means
              }
              
              g <- NULL
              if(grepl(nm, pattern="priority")) {
                flag <- df[, method.name.col] %in% priority.methods
                ret <- plot.strip.plots(df[flag, ], id.var = method.name.col, cell.type.var = cell.type.col, var = "cor",
                                        method.levels = method.levels[[lvl.entry]],
                                        cell.type.levels = rev(cell.type.levels[[lvl.entry]]),
                                        label = "Pearson Correlation")
                g <- ret[["g"]]
                df <- ret[["df"]]
                lvls <- levels(df[,method.name.col])
                lvls <- lvls[lvls %in% df[,method.name.col]]
                y.bold.labels <- ifelse(lvls %in% comparator.methods, yes = "bold", no = "plain")
                print("here\n")
                print(comparator.methods)
                print(priority.methods)
                print(lvls)
                print(y.bold.labels)
              } else {
                # Exclude ensemble, which has NAs for pearson
                flag <- !(df[, method.name.col] %in% c("consensus rank", "ensemble"))
                ret <- plot.strip.plots(df[flag,], id.var = method.name.col, cell.type.var = cell.type.col, var = "cor",
                                        method.levels = method.levels[[lvl.entry]],
                                        cell.type.levels = rev(cell.type.levels[[lvl.entry]]),
                                        label = "Pearson Correlation")
                g <- ret[["g"]]
                df <- ret[["df"]]
                lvls <- levels(df[,method.name.col])
                lvls <- lvls[lvls %in% df[,method.name.col]]
                y.bold.labels <- ifelse(lvls %in% comparator.methods, yes = "bold", no = "plain")
                print("there\n")
                print(comparator.methods)
                print(lvls)
                print(y.bold.labels)
                
                
              }
              g <- g + theme(axis.text.y = element_text(face = y.bold.labels))
              # g <- g + scale_y_discrete(labels = function(x) bold.highlight(x, comparator.methods))
              # g <- g + theme(axis.text.y=element_markdown())
            })
    
    ## "boxplots" = boxplots,
    ret.list <- list("median.bootstrapped.scores" = median.bootstrapped.scores,
                     "mean.bootstrapped.scores" = mean.bootstrapped.scores,
                     "barplots" = barplots,
                     "strip.plots" = strip.plots,
                     "heatmaps" = heatmaps,
                     "spearman.heatmaps" = spearman.heatmaps,			 
                     "pearson.fc.heatmaps" = pearson.fc.heatmaps,			 
                     "merged.all.means" = merged.all.means)			 
    
    ret.list
    
  }




