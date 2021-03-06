### Environment ###############################################################
library("ggplot2");

significance <- 0.05;

### Functions #################################################################
# multiplot
# Description:
#   Multiple ggplot types on the same canvas
# Inputs:
#   ...:        ggplot objects can be passed in ...
#   plotlist:   list of ggplot objects (can be used instead of ...)
#   cols:       number of columns in the layout
#   layout:     matrix specifying the layout. If present, 'cols' is ignored
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
    require(grid);

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist);

    numPlots = length(plots);

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(
            seq(1, cols * ceiling(numPlots/cols)),
            ncol = cols,
            nrow = ceiling(numPlots/cols)
        );
    }

    if (numPlots == 1) {
        print(plots[[1]]);
    } else {
        # Set up the page
        grid.newpage();
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))));

        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE));

            print(
                plots[[i]],
                vp = viewport(
                    layout.pos.row = matchidx$row,
                    layout.pos.col = matchidx$col
                )
            );
        }
    }
}

# volcano_multi
# Description:
#   Volcano plot for multiple datasets
# Inputs:
#   input_data:     Data to be plotted. Should be in 2 ("X" and "Y") or 3 ("X", "Y", "Dataset") columns
#   filename:       Figure file name
#   title:          Plot title
#   xlab:           X label
#   ylab:           Y label
#   xlim:           X-axis boundary
#   transparent:    Transparent background or not
#   width:          Figure width
#   height:         Figure height
#   res:            Figure resolution
#   colours:        Array of colours
volcano_multi <- function(
    input_data,
    filename    = NULL,
    title       = NULL,
    xlab        = NULL,
    ylab        = NULL,
    xlim        = NULL,
    transparent = FALSE,
    width       = 3000,
    height      = 1800,
    res         = 300,
    colours     = NULL,
    ...) {

    p <- (
        ggplot(
            data = input_data,
            mapping = aes(x = X, y = Y, colour = Dataset)
        )
        + geom_point(alpha = 0.2)
        + geom_hline(yintercept = -log10(significance), linetype = 2)
    );
    if (!is.null(title)) {
        p <- p + ggtitle(title);
    }
    if (!is.null(xlab)) {
        p <- p + xlab(xlab);
    }
    if (!is.null(ylab)) {
        p <- p + ylab(ylab);
    }
    if (!is.null(xlim)) {
        p <- p + xlim(c(-xlim, xlim));
    }
    if (!is.null(input_data[,"Dataset"])) {
        p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1)));
    }
    if (transparent) {
        p <- p + theme(
            axis.text.x = element_text(angle = 90),
            panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
            panel.grid.minor = element_blank(), 
            panel.grid.major = element_blank(),
            plot.background = element_rect(fill = "transparent",colour = NA)
        );
        
    }

    if (is.null(filename)) {
        print(p);
    } else {
        png(filename, width = width, height = height, res = res);
        print(p);
        dev.off();
    }
}

# volcano_single
# Description:
#   Volcano plot for a single dataset
# Inputs:
#   input_data:     Data to be plotted. Should be in 2 ("X" and "Y") columns
#   filename:       Figure file name
#   title:          Plot title
#   xlab:           X label
#   ylab:           Y label
#   xlim:           X-axis boundary
#   transparent:    Transparent background or not
#   width:          Figure width
#   height:         Figure height
#   res:            Figure resolution
#   colour:         Plot colour
volcano_single <- function(
    input_data,
    GSE         = NULL,
    filename    = NULL,
    title       = NULL,
    xlab        = NULL,
    ylab        = NULL,
    xlim        = NULL,
    transparent = FALSE,
    width       = 3000,
    height      = 1800,
    res         = 300,
    colour      = "#000000",
    ...) {

    p <- (
        ggplot(
            data = subset(input_data, Dataset == GSE),
            mapping = aes(x = X, y = Y, colour = Dataset)
        )
        + geom_point(alpha = 0.2)
        + scale_colour_discrete(guide = FALSE, drop = TRUE, limits = levels(input_data$Dataset))
        + geom_hline(yintercept = -log10(significance), linetype = 2)
    );
    if (!is.null(title)) {
        p <- p + ggtitle(title);
    }
    if (!is.null(xlab)) {
        p <- p + xlab(xlab);
    }
    if (!is.null(ylab)) {
        p <- p + ylab(ylab);
    }
    if (!is.null(xlim)) {
        p <- p + xlim(c(-xlim, xlim));
    }
    if (transparent) {
        p <- p + theme(
            axis.text.x = element_text(angle = 90),
            panel.background = element_rect(fill = "transparent", colour = NA),
            panel.grid.minor = element_blank(), 
            panel.grid.major = element_blank(),
            plot.background = element_rect(fill = "transparent", colour = NA)
        );
        
    }

    if (is.null(filename)) {
        print(p);
    } else {
        png(filename, width = width, height = height, res = res);
        print(p);
        dev.off();
    }

    return(p);
}