##########################
# Manhattan and QQ plots #
##########################

######################
# Configuration file #
######################
source("~/.saige_pipe.config")
library(harmonicmeanp)
library(ggplot2)
library(hash)
library(data.table)

#############
# Functions #
#############
# Plotting function
pl <- function(expr, filename = paste0("/Net/fs1/home/wilson/public_html/files/Rplots.png"), height = 6, width = 10, units = "in", res = 300, ...) {
    bitmap(filename, height = height, width = width, units = units, res = res, ...)
    tryCatch(expr, error = function(e) cat(e$message, "\n"))
    dev.off()
    invisible()
}
# Simes' test counterpart to HMP function
p_simes <- function(p, w = NULL, L = NULL, w.sum.tolerance = 1e-6, multilevel = TRUE) {
    if (is.null(L) & multilevel) {
        warning("L not specified: for multilevel testing set L to the total number of individual p-values")
        L <- length(p)
    }
    if (length(p) == 0) {
        return(NA)
    }
    if (length(p) > L) stop("The number of p-values cannot exceed L")
    if (is.null(w)) {
        w <- rep(1 / L, length(p))
    } else {
        if (any(w < 0)) stop("No weights can be negative")
        if (length(w) != length(p)) stop("When specified, length of w must equal length of p")
    }
    w.sum <- sum(w)
    if (w.sum > 1 + w.sum.tolerance) {
        stop("Weights cannot exceed 1")
    }
    pstar <- p / w
    return(c(p.Simes = min(sort(pstar) / (1:length(p)))))
}
# Genomic control inflation factor
calc_inflation_factor <- function(pval) {
    qchisq(median(pval), 1, low = FALSE) / qchisq(.5, 1)
}



##########################
# Command-line arguments #
##########################
help <- paste(
    "Usage: Rscript hgi-manhattan.R stem stratum (highlight)",
    sep = "\n"
)
# Argument
lg$args <- commandArgs(trailingOnly = TRUE)
print("Arguments:")
print(lg$args)
if (length(lg$args) != 2 && length(lg$args) != 3) {
    cat(help, sep = "\n")
    stop("\nIncorrect usage\n")
}
# outfile prefix
lg$stem <- as.character(lg$args[1])
if (is.na(lg$stem) || lg$stem == "") stop("Nonsensical stem")

# stratification group
lg$stratum <- as.character(lg$args[2])
if (is.na(lg$stratum) || lg$stratum == "") stop("Nonsensical stratum")

################
# Code version #
################
# Record GitLab version
print("User:")
print((lg$username <- Sys.info()["user"]))
print("Source directory:")
print((lg$srcdir <- config$srcdir))
print("Git repository version")
system(paste0("(cd ", lg$srcdir, " && git show --oneline -s)"))

##########################
# Input and output files #
##########################

lg$infile <- paste0(config$wrkdir, "/summary.", lg$stem, ".", lg$stratum, ".txt.gz")

# Plot titles
main <- paste0(lg$stem, ".", lg$stratum)
# Output files
lg$manhattan_30MAC_file <- paste0(config$wrkdir, "/Manhattan.30MAC.", lg$stem, ".", lg$stratum, "_regenie.png") # nolint: line_length_linter.
manhattan_filee <- paste0(config$wrkdir, "/Manhattan.", lg$stem, ".", lg$stratum, "_regenie.png")
lg$QQplot_file <- paste0(config$wrkdir, "/QQplot.", lg$stem, ".", lg$stratum, "_regenie.png")
lg$log_rds_outfile <- paste0(config$wrkdir, "/log.manhattan.", lg$stem, ".", lg$stratum, "_regenie.rds") # nolint: line_length_linter.
lg$highlight <- FALSE

# Enclose what follows in a tryCatch statement so the log is output
# even if R quits with an error

tryCatch(
    {
        ####################
        # Load the results #
        ####################
        if (lg$highlight) {
            hgi_top <- read.table(lg$highlight_file, sep = "\t", stringsAsFactors = F)
            colnames(hgi_top) <- c("analysis", "rsid", "chr", "pos")
        }
        res <- fread(
            lg$infile,
            header = TRUE,
            colClasses = c(
                CHROM = "character",
                GENPOS = "numeric",
                ID = "character",
                ALLELE0 = "character",
                ALLELE1 = "character",
                A1FREQ = "numeric",
                INFO = "numeric",
                N = "numeric",
                TEST = "character",
                BETA = "numeric",
                SE = "numeric",
                CHISQ = "numeric",
                LOG10P = "numeric",
                EXTRA = "numeric"
            ),
            sep = " ",
            nThread = 7
        )
        attach(res)
        print("Number of rows read:")
        print(nrow(res))

        # Prepare chromosome names and positions
        chr_ftr <- factor(CHROM, levels = c(1:23))
        chr_maxpos <- aggregate(GENPOS, by = list(chr_ftr), max)$x
        chr_cumpos <- cumsum(c(0, chr_maxpos))
        chr_midpos <- 0.5 * (chr_cumpos[-1] + chr_cumpos[-length(chr_cumpos)])
        chr_labels <- as.character(levels(chr_ftr))
        chr_num <- as.numeric(chr_ftr)
        cumpos <- GENPOS + chr_cumpos[chr_num]
        chrcol <- c("grey30", "blue3")[1 + chr_num %% 2]
        if (lg$highlight) {
            hgi_top$cumpos <- hgi_top$pos + chr_cumpos[hgi_top$chr]
        }

        # Pre-calculate 'significance' (-log10 p-value)
        signif <- LOG10P
        pval <- 10^(-1 * signif)
        mac <- pmin(A1FREQ * N * 2, (1 - A1FREQ) * N * 2)
        maf <- pmin(A1FREQ, 1 - A1FREQ)

        # Compute inflation factors
        lg$lambda <- calc_inflation_factor(pval)
        print("Inflation factor lambda (all variants):")
        print(lg$lambda)

        # Inflation factor by MAF
        lg$lambda.maf <- c(
            calc_inflation_factor(pval[0.1 <= maf]),
            calc_inflation_factor(pval[0.01 <= maf & maf < 0.1]),
            calc_inflation_factor(pval[0.001 <= maf & maf < 0.01]),
            calc_inflation_factor(pval[0.0001 <= maf & maf < 0.001]),
            calc_inflation_factor(pval[maf < 0.0001])
        )
        print("Inflation factor lambda (0.1    <= MAF         ):")
        print(lg$lambda.maf[1])
        print("Inflation factor lambda (0.01   <= MAF < 0.1   ):")
        print(lg$lambda.maf[2])
        print("Inflation factor lambda (0.001  <= MAF < 0.01  ):")
        print(lg$lambda.maf[3])
        print("Inflation factor lambda (0.0001 <= MAF < 0.001 ):")
        print(lg$lambda.maf[4])
        print("Inflation factor lambda (          MAF < 0.0001):")
        print(lg$lambda.maf[5])

        # Generally filter at 30 MAC
        lg$mac.threshold <- 30
        gd <- mac >= lg$mac.threshold
        print("Number of variants at 30 <= MAC:")
        print(sum(gd))

        # Minimum p-value
        print("Minimum p-value (30 <= MAC):")
        lg$min.p <- min(pval[gd])
        print(lg$min.p)

        # Harmonic mean p-value
        print("Raw harmonic mean p-value (30 <= MAC):")
        lg$raw.hmp <- 1 / mean(1 / pval[gd])
        print(lg$raw.hmp)
        print("Asymptotically exact harmonic mean p-value (30 <= MAC):")
        lg$p.hmp <- p.hmp(pval[gd], L = sum(gd))
        print(lg$p.hmp)
        print("Simes' test p-value (30 <= MAC):")
        lg$p_simes <- p_simes(pval[gd], L = sum(gd))
        print(lg$p_simes)

        # Manhattan plot: MAC 30
        gdgd <- gd
        data <- data.frame(
            cumpos = cumpos[gdgd],
            signif = signif[gdgd],
            col = chrcol[gdgd],
            rsid = ID[gdgd],
            stringsAsFactors = FALSE
        )

        data$is_highlight <- "none"
        if (lg$highlight) {
            rsid_hash <- hash(keys = data$rsid, values = 1:nrow(data))
            for (i in 1:nrow(hgi_top)) {
                data$is_highlight[rsid_hash[[hgi_top$rsid[i]]]] <- hgi_top$analysis[i]
            }
        }
        sub_data1 <- data[data$signif > 2 | data$is_highlight != "none", ]
        sub_data2 <- data[data$signif < 2 & data$is_highlight == "none", ]
        sub_data2 <- sub_data2[sample(1:nrow(sub_data2), 1000000), ]
        sub_data <- rbind(sub_data1, sub_data2)

        g <- ggplot(sub_data, aes(x = cumpos, y = signif)) +
            geom_point(aes(color = as.factor(col)), alpha = 0.8, size = 0.2) +
            scale_color_manual(values = rep(c("grey10", "blue"))) +
            geom_hline(yintercept = -log10(5e-8), , linetype = "dotted", size = 0.5) +
            scale_x_continuous(label = chr_labels[-24], breaks = chr_midpos[-24]) +
            # scale_x_continuous(label = chr_labels[-c(23, 24)], breaks = chr_midpos[-c(23, 24)] ) +
            theme_bw() +
            theme(
                aspect.ratio = 0.4,
                legend.position = "none",
                panel.border = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank()
            )
        if (lg$highlight) {
            analysis <- strsplit(strsplit(lg$stratum, "\\.")[[1]][1], "")[[1]][1]
            # highlight HGI top SNP
            g <- g + geom_point(
                data = subset(data, is_highlight == analysis),
                color = "orange", size = 0.3
            )
            # highlight relevant position
            g <- g + geom_vline(
                xintercept = hgi_top$cumpos[hgi_top$analysis == analysis],
                linetype = "dotted",
                alpha = 1, size = 0.15
            )
        }
        ggsave(lg$manhattan_30MAC_file, g, dev = "png", dpi = 400)


        # Manhattan plot: all frequencies
        gdgd <- TRUE | gd
        data <- data.frame(
            cumpos = cumpos[gdgd],
            signif = signif[gdgd],
            col = chrcol[gdgd],
            rsid = ID[gdgd],
            stringsAsFactors = FALSE
        )
        data$is_highlight <- "none"
        if (lg$highlight) {
            for (i in 1:nrow(hgi_top)) {
                data$is_highlight[which(data$rsid == hgi_top$rsid[i])] <- hgi_top$analysis[i]
            }
        }
        sub_data1 <- data[data$signif > 2 | data$is_highlight != "none", ]
        sub_data2 <- data[data$signif < 2 & data$is_highlight == "none", ]
        sub_data2 <- sub_data2[sample(1:nrow(sub_data2), 1000000), ]
        sub_data <- rbind(sub_data1, sub_data2)

        g <- ggplot(sub_data, aes(x = cumpos, y = signif)) +
            geom_point(aes(color = as.factor(col)), alpha = 0.8, size = 0.2) +
            scale_color_manual(values = rep(c("grey10", "blue"))) +
            geom_hline(yintercept = -log10(5e-8), linetype = "dotted", size = 0.5) +
            scale_x_continuous(label = chr_labels[-24], breaks = chr_midpos[-24]) +
            # scale_x_continuous(label = chr_labels[-c(23, 24)], breaks = chr_midpos[-c(23, 24)]) +
            theme_bw() +
            theme(
                aspect.ratio = 0.4,
                legend.position = "none",
                panel.border = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank()
            )
        if (lg$highlight) {
            analysis <- strsplit(strsplit(lg$stratum, "\\.")[[1]][1], "")[[1]][1]
            if (strsplit(lg$stratum, "\\.")[[1]][1] == "lenient") {
                analysis <- "A"
            }
            # highlight HGI top SNP
            g <- g + geom_point(
                data = subset(data, is_highlight == analysis),
                color = "orange",
                size = 0.3
            )
            # highlight relevant position
            g <- g + geom_vline(
                xintercept = hgi_top$cumpos[hgi_top$analysis == analysis],
                linetype = "dotted",
                alpha = 1, size = 0.15
            )
        }

        png(manhattan_filee)
        plot(g)
        dev.off()
        ggsave(manhattan_filee, g, dev = "png", dpi = 400)

        # QQplot
        qqfun <- list()
        qqfun_max <- list()
        maf_lwr <- c(0, 0.0001, 0.001, 0.01, 0.1)
        maf_upr <- c(0.0001, 0.001, 0.01, 0.1, 1)
        for (i in 1:length(maf_lwr)) {
            gdgd <- maf_lwr[i] <= maf & maf < maf_upr[i]
            if (sum(gdgd) > 0) {
                od <- order(signif[gdgd])
                expt <- -log10((length(od):1) / (length(od) + 1))
                qqfun[[i]] <- approxfun(expt, signif[gdgd][od], method = "linear")
                qqfun_max[[i]] <- max(c(max(expt), max(signif[gdgd])))
            }
        }
        pl(
            {
                xy_max <- max(unlist(qqfun_max), na.rm = T)
                col <- colorRampPalette(c("skyblue", "blue", "black"))(length(maf_lwr))
                plot(
                    c(0, xy_max),
                    c(0, xy_max),
                    type = "n",
                    xlab = "Expected -log10 p-values",
                    ylab = "Observed -log10 p-values",
                    main = main
                )
                abline(0, 1, col = 2)
                for (i in 1:length(maf_lwr)) {
                    if (!is.null(qqfun[[i]])) {
                        curve(qqfun[[i]](x), 0, xy_max, n = 1001, add = TRUE, col = col[i])
                        gdgd <- maf_lwr[i] <= maf & maf < maf_upr[i]
                        od <- order(signif[gdgd])
                        expt <- -log10((length(od):1) / (length(od) + 1))
                        gdgdgd <- expt > 5
                        if (sum(gdgdgd) > 0) {
                            points(expt[gdgdgd], signif[gdgd][od][gdgdgd], col = col[i], cex = 0.7)
                        }
                    }
                }
                legend("topleft", c(
                    "0.1       <= MAF         ",
                    "0.01     <= MAF < 0.1   ",
                    "0.001   <= MAF < 0.01  ",
                    "0.0001 <= MAF < 0.001 ",
                    "                MAF < 0.0001",
                    paste0("lambda = ", round(100 * lg$lambda.maf) / 100)
                ), ,
                lwd = rep(c(1, NA), each = 5), col = rep(rev(col), 2), bty = "n", ncol = 2
                )
            },
            filename = lg$QQplot_file
        )

        # Any significant variants (MAF >= 1%)
        print("Significant hits (0.01 <= MAF)")
        lg$hits <- res[pval < 5e-8 & gd, ]
        print(lg$hits)
    },
    finally = {
        # On error or clean exit
        saveRDS(lg, file = lg$log_rds_outfile)
    }
)
